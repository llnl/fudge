# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=0 data ('normal' reaction evaluations for various projectiles)
"""

import sys
import copy
import fractions

from pqu import PQU as PQUModule
from xData import enums as xDataEnumsModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs import alias as aliasPoPsModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.decays import misc as miscDecaysModule
from PoPs.decays import decayData as decayDataModule
from PoPs.decays import probability as probabilityModule
from PoPs.decays import product as productModule
from PoPs.chemicalElements import misc as chemicalElementMiscModule
from PoPs.families import nuclide as nuclideModule
from PoPs.families import gaugeBoson as gaugeBosonModule

from . import endfFileToGNDSMisc as endfFileToGNDSMiscModule
from .ENDF_ITYPE_0_Misc import BadResonances, getTotalOrPromptFission, getDelayedFission, readMF2, \
    readMF8, parseReaction, parseCovariances, promptToken, totalToken

from fudge import enums as enumsModule
from fudge import sums as sumsModule
from fudge import outputChannel as outputChannelModule
from fudge.reactions import reaction as reactionModule
from fudge.reactions import orphanProduct as orphanProductModule
from fudge.reactions import production as productionModule
from fudge.reactions import incompleteReaction as incompleteReactionModule
from fudge.reactions import fissionComponent as fissionComponentModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import unspecified as unspecifiedModule
from fudge.productData.distributions import branching3d as branching3dModule
from fudge.resonances import common as resonancesCommonModule

from .. import toGNDSMisc as toGNDSMiscModule
from .. import endf_endl as endf_endlModule


def deriveMT3MF3FromMT1_2(info, reactionSuite):
    totalCrossSection, elasticCrossSection = None, None
    for reaction in reactionSuite.sums.crossSectionSums:
        if reaction.ENDF_MT == 1:
            totalCrossSection = reaction.crossSection.evaluated
            break

    for reaction in reactionSuite.reactions:
        if reaction.ENDF_MT == 2:
            elasticCrossSection = reaction.crossSection.evaluated
            break
    if totalCrossSection is None: raise Exception('No total cross section for calculating non-elastic')
    if elasticCrossSection is None: raise Exception('No elastic cross section for calculating non-elastic')

    try:
        form = totalCrossSection - elasticCrossSection
    except Exception:
        totalCrossSection = totalCrossSection.toPointwise_withLinearXYs(accuracy=1e-3, lowerEps=1e-6)
        elasticCrossSection = elasticCrossSection.toPointwise_withLinearXYs(accuracy=1e-3, lowerEps=1e-6)
        form = totalCrossSection - elasticCrossSection

    form.label = info.style
    return form


def checkProductList(outputChannel, productList):
    if outputChannel is None:
        return
    for product in outputChannel:
        if product.pid in productList:
            try:
                multiplicity = int(product.multiplicity.getConstant())
            except:
                continue
            productList[product.pid] -= multiplicity


def ITYPE_0(MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip, parseCrossSectionOnly, doCovariances, verbose,
            reconstructResonances=True):
    warningList = []

    info.totalOrPromptFissionNeutrons = {}
    info.totalMF6_12_13Gammas = {}
    info.PoPsOverrides = {}

    if 452 in MTDatas:
        info.totalOrPromptFissionNeutrons[totalToken] = getTotalOrPromptFission(info, 452, MTDatas, totalToken,
                                                                                warningList)
    if 455 in MTDatas:
        info.delayedFissionDecayChannel = getDelayedFission(info, 455, MTDatas, warningList)
    if 456 in MTDatas:
        info.totalOrPromptFissionNeutrons[promptToken] = getTotalOrPromptFission(info, 456, MTDatas, promptToken,
                                                                                 warningList)
    if 458 in MTDatas:
        info.fissionEnergyReleaseData = MTDatas[458]
    sys.stdout.flush()
    if verbose > 0:
        for warning in warningList: info.logs.write("       WARNING: %s\n" % warning, stderrWriting=True)

    MTList = endfFileToGNDSMiscModule.niceSortOfMTs(list(MTDatas.keys()), verbose=0, logFile=info.logs)
    extraMTs = list(set(MTDatas.keys()).difference(set(MTList)))
    for MT in [151, 451, 452, 455, 456, 458]:
        if MT in extraMTs:
            extraMTs.pop(extraMTs.index(MT))

    haveTotalFission = (18 in MTList)
    fissionMTs = [mt for mt in MTList if mt in (19, 20, 21, 38)]

    continuumMTs = {91: range(50, 91), 649: range(600, 649), 699: range(650, 699), 749: range(700, 749),
                    799: range(750, 799), 849: range(800, 849), 999: range(900, 999)}
    summedReactions = {}
    summedReactionsInfo = {4: range(50, 92), 102: range(900, 1000), 103: range(600, 650), 104: range(650, 700),
                           105: range(700, 750), 106: range(750, 800), 107: range(800, 850)}

    keep102 = False
    for MT in summedReactionsInfo[102]:
        if MT in MTList:
            keep102 = True
    if not keep102:
        del summedReactionsInfo[102]

    for summedMT, partialReactions in summedReactionsInfo.items():
        if summedMT not in MTList: continue
        for MT in MTList:
            if MT in partialReactions:
                summedReactions[summedMT] = None
                break

    for summedMT in (1, 3):
        if summedMT in MTList: summedReactions[summedMT] = None

    try:
        inelasticMT = {1: 91, 1001: 649, 1002: 699, 1003: 749, 2003: 799, 2004: 849}[info.projectileZA]
    except:
        inelasticMT = -1

    MT5Reaction = None
    reactions = []
    fissionComponents = []
    productions = []
    incompleteReactions = []
    nonElastic = []
    delayInsertingSummedReaction = []
    linksToCheck = []  # links that may need to be updated after reading resonances
    SpecialLRProducts = []

    for MT in MTList:
        crossSection = None
        if MT in MTs2Skip: continue
        if (singleMTOnly is not None) and (MT != singleMTOnly): continue

        channelProcess = None
        if MT == inelasticMT:
            channelProcess = outputChannelModule.Processes.continuum
        elif MT in range(102, 108):
            channelProcess = outputChannelModule.Processes.inclusive
        elif MT in continuumMTs:
            for discreteMT in continuumMTs[MT]:
                if discreteMT in MTList: channelProcess = outputChannelModule.Processes.continuum

        warningList = []
        MTData = MTDatas[MT]

        fissionGenre = {18: enumsModule.FissionGenre.total,
                        19: enumsModule.FissionGenre.firstChance,
                        20: enumsModule.FissionGenre.secondChance,
                        21: enumsModule.FissionGenre.thirdChance,
                        38: enumsModule.FissionGenre.fourthChance}.get(MT, enumsModule.FissionGenre.none)

        # Sometimes excited states are identified in MF8.
        # Read this before reading distributions to make sure info is present.
        LMF, radioactiveDatas = readMF8(info, MT, MTData, warningList)

        doParseReaction = 3 in MTData
        if not doParseReaction:
            if MT == 3: doParseReaction = (12 in MTData) or (13 in MTData)
        if doParseReaction:  # normal reaction, with cross section and distributions
            try:
                crossSection, outputChannel, MFKeys, LRProducts = parseReaction(
                    info, info.target, info.projectileZA, info.targetZA, MT, MTData, warningList,
                    parseCrossSectionOnly=parseCrossSectionOnly, channelProcess=channelProcess)
            except KeyboardInterrupt:
                raise
            except:
                import traceback
                info.logs.write(traceback.format_exc(), stderrWriting=True)
                info.doRaise.append(traceback.format_exc())
                info.logs.write('\n')
                sys.stdout.flush()
                continue

            if LRProducts is not None: SpecialLRProducts.append([LRProducts, outputChannel])

            info.logs.write('\n')
            sys.stdout.flush()
            if len(MFKeys):
                warningList.append('For reaction MT = %d, the following MFs were not converted: %s\n' % (MT, MFKeys))
            if outputChannel is None: break

            if MT in summedReactions:
                summedReactions[MT] = [crossSection, outputChannel]
            else:
                if MT != 2: nonElastic.append(MT)
                reaction = reactionModule.Reaction(None, outputChannel.genre, ENDF_MT=MT, fissionGenre=fissionGenre)
                endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, reaction, outputChannel)
                if hasattr(info, 'dSigma_form'):
                    reaction.doubleDifferentialCrossSection.add(info.dSigma_form)
                    del info.dSigma_form
                    crossSection = crossSectionModule.CoulombPlusNuclearElastic(
                        link=reaction.doubleDifferentialCrossSection.evaluated,
                        label=info.style, relative=True)
                reaction.crossSection.add(crossSection)
                if MT == 5:
                    MT5Reaction = reaction
                elif MT in fissionMTs and haveTotalFission:
                    # This is 1st, 2nd, etc. chance fission but total is also present.
                    fissionComponents.append(reaction)
                else:
                    if MT in summedReactionsInfo:
                        delayInsertingSummedReaction.append(reaction)
                    else:
                        reactions.append([MT, reaction])
        else:
            if 12 in MTData:
                # MF=12-15 is occasionally provided without corresponding MF=3, e.g. ENDF-VIII.1 g-069_Tm_169
                from .ENDF_ITYPE_0_Misc import readMF12_13
                productList = []
                _dummyCrossSection = None
                # read branching-ratio data into info.MF12_LO2
                readMF12_13(info, MT, MTData, productList, warningList, crossSection, _dummyCrossSection,
                            gammaBRTolerance=1e-6)

            else:
                MFList = []
                for MF in [4, 5, 6]:
                    if MF in MTData: MFList.append('%d' % MF)
                if MFList: warningList.append(
                    'MT = %d has MF = %s data and no MF 3 data' % (MT, ', '.join(MFList)))

        for radioactiveData in radioactiveDatas:
            # Get radioactive production data (if any) from MF 8-10. Cross section form depends on value of LMF.
            if LMF in [3, 6, 9]:  # Cross-section is reference to MF3.
                productionCrossSection = reaction.crossSection.evaluated
                if MT in summedReactions:
                    productionCrossSection = crossSection
                productionCrossSection = crossSectionModule.Reference(link=productionCrossSection, label=info.style)
                linksToCheck.append(productionCrossSection)
            elif LMF == 10:  # MF10 data is cross-section. Product's multiplicity is 1.
                productionCrossSection = radioactiveData[4]
            else:
                raise Exception("Unknown LMF=%d encountered in MF=8 for MT=%d" % (LMF, MT))

            ZAP = radioactiveData[0]
            ELFS = radioactiveData[1]
            LFS = radioactiveData[2]

            if LMF in [9, 10]:
                Q = toGNDSMiscModule.returnConstantQ(info.style, radioactiveData[6], productionCrossSection)
            else:
                Q = outputChannel.Q[info.style]

            if MT == 18:
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                outputChannel.Q.add(Q)
                icr = incompleteReactionModule.IncompleteReaction('fission', outputChannel.genre, MT,
                                                                  fissionGenre=enumsModule.FissionGenre.total)
                endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, icr, outputChannel)

                icr.crossSection.add(productionCrossSection)
                incompleteReactions.append(icr)
                continue

            if LMF == 6:  # Product multiplicity is in MF6, so production channel multiplicity needs to refer to it:
                residual = toGNDSMiscModule.getTypeNameGamma(info, ZAP, level=ELFS, levelIndex=LFS)
                residualID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(residual.id,
                                                                                         info.specialNuclearParticleID)
                MF6prod = outputChannel.getProductsWithName(residualID)

                if len(MF6prod) != 1:  # problem appears in JEFF-3.2 Y90 and Y91
                    info.missingRadioactiveProduct.append(
                        'Unique MT%d radioactive product %s not found in product list!' % (MT, residualID))
                    continue

                multiplicity = multiplicityModule.Reference(label=info.style, link=MF6prod[0].multiplicity)
            else:
                multiplicity = radioactiveData[3]

            try:
                if crossSection is None:
                    crossSectionData = [[reactionSuite.styles[0].projectileEnergyDomain.min, 1.0],
                                        [reactionSuite.styles[0].projectileEnergyDomain.max, 1.0]]
                    crossSection = crossSectionModule.XYs1d(data=crossSectionData, axes=crossSectionModule.defaultAxes(
                        reactionSuite.styles[0].projectileEnergyDomain.unit))
                residual = toGNDSMiscModule.newGNDSParticle(
                    info, toGNDSMiscModule.getTypeNameGamma(info, ZAP, level=ELFS, levelIndex=LFS),
                    crossSection, multiplicity=multiplicity)
            except:
                info.logs.write('\nMT = %s\n' % MT)
                raise

            productionOutputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.production)
            productionOutputChannel.Q.add(Q)
            productionOutputChannel.products.add(productionOutputChannel.products.uniqueLabel(residual))
            productionOutputChannel.process = "%s%s" % (
                reactionSuite.target, endf_endlModule.endfMTtoC_ProductLists[MT].reactionLabel.replace(
                    'z,', reactionSuite.projectile + ','))

            production = productionModule.Production(None, productionOutputChannel.genre, ENDF_MT=MT)
            endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, production, productionOutputChannel)
            production.crossSection.add(productionCrossSection)
            productions.append(production)

        if verbose > 0:
            for warning in warningList: info.logs.write("       WARNING: %s\n" % warning, stderrWriting=True)

    for MT, reaction in reactions:
        if reaction.label in reactionSuite.reactions:
            warningList.append('Duplicate reaction label "%s".' % reaction.label)
            print('  WARNING: %s' % warningList[-1])
            reaction.outputChannel.process = 'ENDF MT %s' % MT
            reaction.updateLabel()
        if MT not in [5, 18, 19, 20, 21, 38]:
            productList = copy.copy(endf_endlModule.endfMTtoC_ProductLists[MT].productCounts)
            checkProductList(reaction.outputChannel, productList)
            missingProduct = []
            for product, multiplicity in productList.items():
                if product == IDsPoPsModule.photon:
                    continue
                if multiplicity > 0:
                    missingProduct.append(product)
            if len(missingProduct) > 0:
                warningList.append('For MT=%s, product list does not match MT. Missing products are %s.' % (
                    MT, ','.join(missingProduct)))
                print('  WARNING: %s' % warningList[-1])
        reactionSuite.reactions.add(reaction)

    for reaction in delayInsertingSummedReaction:
        reactionSuite.reactions.add(reaction)

    if MT5Reaction is not None:
        reactionSuite.reactions.add(MT5Reaction)

    # BRB, The channelIDs should be in a common area?????
    channelIDs = {1: 'total', 3: 'nonelastic', 4: '(z,n)', 103: '(z,p)', 104: '(z,d)', 105: '(z,t)', 106: '(z,He3)',
                  107: '(z,alpha)'}
    if 3 in summedReactions: summedReactionsInfo[3] = nonElastic
    if 1 in summedReactions:
        if 2 in MTList:
            summedReactionsInfo[1] = [2] + nonElastic
        else:
            summedReactionsInfo[1] = nonElastic
    summedReactionMTs = endfFileToGNDSMiscModule.niceSortOfMTs(list(summedReactions.keys()), verbose=0,
                                                               logFile=info.logs)
    for MT in (4, 3, 1):
        if MT in summedReactionMTs:
            summedReactionMTs.remove(MT)
            summedReactionMTs.insert(0, MT)
    for i1, MT in enumerate(summedReactionMTs):
        if summedReactions[MT] is None: continue
        crossSection, outputChannel = summedReactions[MT]
        omitWhenWritingENDF = False
        if (MT == 3) and (crossSection is None):
            crossSection = deriveMT3MF3FromMT1_2(info, reactionSuite)
            omitWhenWritingENDF = True
        summedCrossSection = sumsModule.CrossSectionSum(label=channelIDs[MT], ENDF_MT=MT)
        for reaction in reactionSuite.reactions:
            if reaction.ENDF_MT in summedReactionsInfo[MT]:
                summedCrossSection.summands.append(sumsModule.Add(link=reaction.crossSection))
        summedCrossSection.Q.add(outputChannel.Q[info.style])
        summedCrossSection.crossSection.add(crossSection)
        if omitWhenWritingENDF:
            info.ENDFconversionFlags.add(summedCrossSection, "omit")
        reactionSuite.sums.crossSectionSums.add(summedCrossSection)

        gammas = []
        for product in outputChannel:
            particle = reactionSuite.PoPs[product.pid]
            if isinstance(particle, gaugeBosonModule.Particle):
                gammas.append(product)
            else:
                if product.outputChannel is not None:
                    for product2 in product.outputChannel:
                        particle = reactionSuite.PoPs[product2.pid]
                        if isinstance(particle, gaugeBosonModule.Particle): gammas.append(product2)
        if len(gammas) > 0:
            productChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
            for QForm in outputChannel.Q: productChannel.Q.add(QForm)
            for gamma in gammas: productChannel.products.add(productChannel.products.uniqueLabel(gamma))
            productionReaction = orphanProductModule.OrphanProduct('Orphan product %s' % str(i1), productChannel.genre,
                                                                   ENDF_MT=MT)
            endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, productionReaction, productChannel)
            crossSectionLink = crossSectionModule.Reference(link=summedCrossSection.crossSection.evaluated,
                                                            label=info.style)
            linksToCheck.append(crossSectionLink)
            productionReaction.crossSection.add(crossSectionLink)
            reactionSuite.orphanProducts.add(productionReaction)

    for i1, reaction in enumerate(fissionComponents):
        # 1st-chance, 2nd-chance, etc. Convert them to fissionComponent instances:
        fissionComponent = fissionComponentModule.FissionComponent(
            None, reaction.outputChannel.genre, reaction.ENDF_MT, fissionGenre=reaction.fissionGenre)
        endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, fissionComponent, reaction.outputChannel)
        for crossSection in reaction.crossSection: fissionComponent.crossSection.add(crossSection)
        reactionSuite.fissionComponents.add(fissionComponent)

    unprocessedMTs = []
    for MT in extraMTs:
        if MT in range(201, 208):
            outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
            label = 'Total %s production' % \
                    {201: 'neutron', 202: 'gamma', 203: 'proton', 204: 'deuteron', 205: 'triton', 206: 'helion',
                     207: 'alpha'}[MT]
            productID = {201: IDsPoPsModule.neutron, 202: IDsPoPsModule.photon, 203: IDsPoPsModule.proton,
                         204: IDsPoPsModule.familiarDeuteron,
                         205: IDsPoPsModule.familiarTriton, 206: IDsPoPsModule.familiarHelion,
                         207: IDsPoPsModule.familiarAlpha}[MT]
            productionReaction = productionModule.Production(label, outputChannel.genre, ENDF_MT=MT)
            crossSection, outputChannel, MFKeys, LRProducts = parseReaction(
                info, info.target, info.projectileZA, info.targetZA, MT, MTDatas[MT], warningList,
                parseCrossSectionOnly=parseCrossSectionOnly, channelProcess=channelProcess)
            if len(outputChannel.products) == 0:
                productID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(
                    productID, specialNuclearParticleIDPoPsModule.Mode.nuclide)
                product = toGNDSMiscModule.newGNDSParticle(info, productID, crossSection)
                product.distribution.add(unspecifiedModule.Form(info.style, xDataEnumsModule.Frame.lab))
                outputChannel.products.add(product)
                info.ENDFconversionFlags.add(product, 'implicitProduct')
            if len(MFKeys) > 0:
                warningList.append('For reaction MT = %d, the following MFs were not converted: %s\n' % (MT, MFKeys))
            if LRProducts is not None:
                warningList.append('For reaction MT = %d, LRProducts is not None\n')
            productionReaction.crossSection.add(crossSection)
            endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, productionReaction, outputChannel)
            productions.append(productionReaction)
            info.ENDFconversionFlags.add(productionReaction, 'MF3production')
            if verbose > 0:
                print()  # put messages for each extraMT on a separate line
        elif MT in range(851, 871):
            continue  # dealt with while parsing covariances
        else:
            unprocessedMTs.append(MT)
    if len(unprocessedMTs) > 0:
        print('INFO: the following MTs were not converted: %s.' % unprocessedMTs)

    for production in productions:
        reactionSuite.productions.add(production)

    for incompleteReaction in incompleteReactions:
        reactionSuite.incompleteReactions.add(incompleteReaction)

    if hasattr(info, 'totalDelayedMultiplicity'):
        fissionOutputchannel = reactionSuite.getReaction('fission').outputChannel

        delayedNubar = sumsModule.MultiplicitySum(label="delayed fission neutron multiplicity", ENDF_MT=455)
        for neutron in fissionOutputchannel.fissionFragmentData.delayedNeutrons:
            delayedNubar.summands.append(sumsModule.Add(link=neutron.product.multiplicity))
        delayedNubar.multiplicity.add(info.totalDelayedMultiplicity)
        reactionSuite.sums.multiplicitySums.add(delayedNubar)

        totalNubar = sumsModule.MultiplicitySum(label="total fission neutron multiplicity", ENDF_MT=452)
        for product in fissionOutputchannel.products:
            if product.pid == IDsPoPsModule.neutron:
                totalNubar.summands.append(sumsModule.Add(product.multiplicity))
        totalNubar.summands.append(sumsModule.Add(delayedNubar.multiplicity))
        totalNubar.multiplicity.add(info.totalOrPromptFissionNeutrons[totalToken])
        reactionSuite.sums.multiplicitySums.add(totalNubar)

    warningList = []
    try:  # Parse resonance section.
        mf2 = None
        if 151 in MTDatas and not parseCrossSectionOnly:
            mf2 = MTDatas.get(151).get(2)  # Resonance data available.
        if mf2:
            info.logs.write('    Reading resonances (MF=2 MT=151)\n')
            resonances, resonanceMTs = readMF2(info, mf2, warningList)
            if info.LRP == 2:  # LRP was read in from first line of ENDF file
                if resonances.resolved: resonances.resolved.evaluated.useForSelfShieldingOnly = True
                if resonances.unresolved: resonances.unresolved.evaluated.useForSelfShieldingOnly = True
            reactionSuite.resonances = resonances

            if resonances.reconstructCrossSection:
                # modify cross sections for relevant channels to indicate resonance contribution is needed:
                resonanceLink = crossSectionModule.ResonanceLink(link=resonances)

                for MT in resonanceMTs:
                    MTChannels = [r1 for r1 in reactionSuite.reactions if (r1.ENDF_MT == MT)]
                    MTChannels += [r1 for r1 in reactionSuite.sums.crossSectionSums if (r1.ENDF_MT == MT)]
                    MTChannels += [r1 for r1 in reactionSuite.fissionComponents if (r1.ENDF_MT == MT)]
                    if len(MTChannels) == 0:
                        if MT in (3, 18, 19):
                            continue
                        else:
                            warningList.append('Unable to find channel corresponding to resonance data for MT%d' % MT)
                    elif len(MTChannels) == 1:
                        crossSectionComponent = MTChannels[0].crossSection

                        if isinstance(crossSectionComponent[info.style], crossSectionModule.CoulombPlusNuclearElastic):
                            continue  # Don't make any changes to Coulomb elastic scattering.

                        # break background up into resolved, unresolved and fast regions:
                        haveResolved = (resonances.resolved is not None and not
                        resonances.resolved.evaluated.useForSelfShieldingOnly)
                        haveUnresolved = (resonances.unresolved is not None and not
                        resonances.unresolved.evaluated.useForSelfShieldingOnly)

                        if MT in (18, 19) and not haveUnresolved:
                            # convert to resonancesWithBackground only if resolved region includes fission widths,
                            # not for threshold fissioners like Th232
                            resolved = resonances.resolved.evaluated
                            if isinstance(resolved, resonancesCommonModule.EnergyIntervals):
                                resolved = resolved[-1].evaluated
                            if hasattr(resolved, 'resonanceReactions'):
                                if not any([r.isFission() for r in resolved.resonanceReactions]):
                                    continue
                            else:  # Breit-Wigner, must resort to looking at column names
                                if not any(['fission' in c.name for c in resolved.resonanceParameters.table.columns]):
                                    continue

                        originalBackgroundForm = backgroundForm = crossSectionComponent[info.style]

                        if isinstance(backgroundForm, crossSectionModule.XYs1d):
                            tmp = crossSectionModule.Regions1d(axes=backgroundForm.axes)
                            tmp.append(backgroundForm)
                            backgroundForm = tmp

                        if MT == 3 and backgroundForm.domainMin > resonances.scatteringRadius.evaluated.domainMin:
                            # Address error appearing in several JEFF evaluations
                            warningList.append(
                                "MT3 background cross section does not span resonance region. Assuming background == 0")
                            backgroundForm.prepend(crossSectionModule.XYs1d(
                                axes=backgroundForm.axes,
                                data=[[resonances.scatteringRadius.evaluated.domainMin, 0],
                                      [backgroundForm.domainMin, 0]])
                            )

                        if haveUnresolved:
                            maxResonanceDomain = resonances.unresolved.domainMax
                        else:
                            maxResonanceDomain = resonances.resolved.domainMax
                        if backgroundForm.domainMin > maxResonanceDomain:
                            # assume resonances were provided for interference only, not for reconstructing cross section
                            continue

                        haveFast = backgroundForm.domainMax > maxResonanceDomain

                        backgroundBoundaries = [region.domainMax for region in backgroundForm[:-1]]

                        missingBoundaries = {}
                        if haveResolved and (
                                haveUnresolved or haveFast) and resonances.resolved.domainMax not in backgroundBoundaries:
                            missingBoundaries[resonances.resolved.domainMax] = ['mergeWithNextRegion']
                        if haveUnresolved and haveFast and resonances.unresolved.domainMax not in backgroundBoundaries:
                            missingBoundaries[resonances.unresolved.domainMax] = ['mergeWithNextRegion']

                        if missingBoundaries:
                            # ENDF manual doesn't insist on breaking the background into regions, but GNDS does
                            bkForm2 = crossSectionModule.Regions1d(axes=backgroundForm.axes)

                            domainMin = backgroundForm.domainMin
                            if isinstance(backgroundForm, crossSectionModule.XYs1d):
                                for boundary in sorted(missingBoundaries.keys()):
                                    if boundary not in backgroundForm.domainGrid:
                                        missingBoundaries[boundary].append('deletePoint')
                                    bkForm2.append(backgroundForm.domainSlice(domainMin, boundary))
                                    domainMin = boundary
                                bkForm2.append(backgroundForm.domainSlice(domainMin, backgroundForm.domainMax))

                            else:  # Regions1d
                                startIdx = 0
                                for boundary in sorted(missingBoundaries.keys()):
                                    for idx in range(startIdx, len(backgroundForm)):
                                        if backgroundForm[idx].domainMax > boundary: break
                                        bkForm2.append(backgroundForm[idx])
                                    # add the missing boundary:
                                    if boundary not in backgroundForm[idx].domainGrid:
                                        missingBoundaries[boundary].append('deletePoint')
                                    bkForm2.append(
                                        backgroundForm[idx].domainSlice(domainMin, boundary))
                                    domainMin = boundary
                                    startIdx = idx
                                bkForm2.append(
                                    backgroundForm[idx].domainSlice(domainMin, backgroundForm[idx].domainMax))

                                for region in backgroundForm[idx + 1:]:
                                    bkForm2.append(region)

                            backgroundForm = bkForm2

                        RRBack = URRBack = fastBack = None

                        idx = 0
                        if haveResolved:
                            while backgroundForm[idx].domainMax < resonances.resolved.domainMax:
                                idx += 1
                            if backgroundForm[idx].domainMax > resonances.resolved.domainMax:
                                warningList.append("Domain mismatch between resolved resonances and background")
                                info.doRaise.append(warningList[-1])
                            if idx == 0:
                                RRBack = backgroundForm[idx]
                                RRBack.index = None
                            else:
                                RRBack = crossSectionModule.Regions1d(axes=backgroundForm.axes)
                                for jdx in range(idx + 1):
                                    RRBack.append(backgroundForm[jdx])
                            RRBack = crossSectionModule.ResolvedRegion(RRBack)
                            if resonances.resolved.domainMax in missingBoundaries:
                                info.ENDFconversionFlags.add(RRBack,
                                                             ",".join(missingBoundaries[resonances.resolved.domainMax]))
                            idx += 1

                        startIdx = idx
                        if haveUnresolved:
                            while backgroundForm[idx].domainMax < resonances.unresolved.domainMax:
                                idx += 1
                            if backgroundForm[idx].domainMax != resonances.unresolved.domainMax:
                                warningList.append("Domain mismatch between unresolved resonances and background")
                                info.doRaise.append(warningList[-1])
                            if idx == startIdx:
                                URRBack = backgroundForm[idx]
                                URRBack.index = None
                            else:
                                URRBack = crossSectionModule.Regions1d(axes=backgroundForm.axes)
                                for jdx in range(startIdx, idx + 1):
                                    URRBack.append(backgroundForm[jdx])
                            URRBack = crossSectionModule.UnresolvedRegion(URRBack)
                            if resonances.unresolved.domainMax in missingBoundaries:
                                info.ENDFconversionFlags.add(URRBack, ",".join(
                                    missingBoundaries[resonances.unresolved.domainMax]))
                            idx += 1

                        startIdx = idx
                        if startIdx < len(backgroundForm):
                            if startIdx == len(backgroundForm) - 1:
                                fastBack = backgroundForm[startIdx]
                                fastBack.index = None
                            else:
                                fastBack = crossSectionModule.Regions1d(axes=backgroundForm.axes)
                                for jdx in range(startIdx, len(backgroundForm)):
                                    fastBack.append(backgroundForm[jdx])
                            fastBack = crossSectionModule.FastRegion(fastBack)

                        background_ = crossSectionModule.Background(RRBack, URRBack, fastBack)
                        crossSectionComponent.pop(info.style)
                        crossSectionComponent.add(crossSectionModule.ResonancesWithBackground(
                            info.style, resonanceLink, background_, backgroundForm.uncertainty))
                        for link in linksToCheck:
                            if link.link is originalBackgroundForm:
                                link.link = crossSectionComponent[info.style]
                    else:
                        raise NotImplementedError("Multiple reactions match resonance MT%d" % MT)

    except BadResonances as e:
        warningList.append('       ERROR: unable to parse resonances! Error message: %s' % e)
        info.doRaise.append(warningList[-1])

    MF12BaseMTsAndRange = [[50, 92], [600, 650], [650, 700], [700, 750], [750, 800], [800, 850], [900, 1000]]

    if singleMTOnly is None:
        for MTLO2, MF12_LO2 in sorted(info.MF12_LO2.items()):
            # The logic below assumes MTs are in ascending order per base MT.
            branchingBaseMT = None
            for MTBase, MTEnd in MF12BaseMTsAndRange:  # Determine base MT for this MTLO2
                if MTBase < MTLO2 < MTEnd:
                    branchingBaseMT = MTBase
                    break
            if branchingBaseMT is not None:
                residualZA = endf_endlModule.ENDF_MTZAEquation(info.projectileZA, info.targetZA, branchingBaseMT)[0][-1]
                residual = toGNDSMiscModule.getTypeNameENDF(info, residualZA, None)
                residualName = residual.id
                level = MTLO2 - branchingBaseMT
                levelEnergy = MF12_LO2[0]['ES']
                fullName = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex(residualName, level)

                decayProduct = None
                for reaction in reactionSuite.reactions:
                    if MTLO2 == reaction.ENDF_MT:
                        if reaction.outputChannel.genre == enumsModule.Genre.twoBody:
                            _decayProduct = reaction.outputChannel[1]
                            if _decayProduct.pid == fullName:
                                if _decayProduct.outputChannel is None:
                                    decayProduct = _decayProduct
                                    break
                if decayProduct is None:
                    pass
                else:
                    crossSection = reaction.crossSection[0]
                    decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                    decayChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, levelEnergy, crossSection))

                    decayResidualName = reactionSuite.PoPs[fullName].isotope.symbol
                    decayChannel.products.add(toGNDSMiscModule.newGNDSParticle(info, decayResidualName, crossSection))

                    multiplicity = multiplicityModule.Branching1d(info.style)
                    decayPhoton = toGNDSMiscModule.newGNDSParticle(
                        info, IDsPoPsModule.photon, crossSection, multiplicity=multiplicity)
                    decayPhoton.distribution.add(branching3dModule.Form(info.style, xDataEnumsModule.Frame.lab))
                    decayChannel.products.add(decayPhoton)

                    decayProduct.addOutputChannel(decayChannel)
                    reaction.updateLabel()

                if fullName not in reactionSuite.PoPs:
                    # add missing particle
                    finalParticle = toGNDSMiscModule.getPoPsParticle(info, residualZA, levelIndex=level,
                                                                     level=levelEnergy, levelUnit='eV')

                particleLevelEnergy_eV = reactionSuite.PoPs[fullName].energy[0].value
                # Compare this value to level energy from the particle list (from MF3 Q-value).
                if levelEnergy != particleLevelEnergy_eV:
                    if particleLevelEnergy_eV < 1e-12:
                        warningList.append("MF12 parent level energy (%s) set to zero?" % particleLevelEnergy_eV)
                        info.doRaise.append(warningList[-1])
                    elif abs(levelEnergy - particleLevelEnergy_eV) < 1e-4 * particleLevelEnergy_eV:
                        MFLabel = '3'
                        # Value with most precision wins.
                        str1 = PQUModule.floatToShortestString(
                            levelEnergy * 1e-20)  # 1e-20 to insure e-form is returned.
                        str2 = PQUModule.floatToShortestString(
                            particleLevelEnergy_eV * 1e-20)  # Want 1.23e-16 and not 12300 to differ
                        if len(str1) > len(str2):  # length from 1.2345e-16 and not 12345.
                            reactionSuite.PoPs[fullName].energy[0].value = levelEnergy
                            MFLabel = '12'
                        if str1 != str2:
                            warningList.append(
                                "MT%d MF12 level energy %s differs from MF3 value %s. Setting to MF%s value." %
                                (MTLO2, levelEnergy, particleLevelEnergy_eV, MFLabel))
                    else:
                        warningList.append("MT%d MF12 parent level energy (%s) doesn't match known level" %
                                           (MTLO2, particleLevelEnergy_eV))
                        info.doRaise.append(warningList[-1])
                for i1, MF12 in enumerate(MF12_LO2):
                    try:
                        finalLevelEnergy = MF12['ESk']
                        if finalLevelEnergy > 0.:  # Find particle in the particleList with energy = finalLevelEnergy
                            finalParticles = [lev for lev in reactionSuite.PoPs[residualName].ancestor
                                              if lev.energy.float('eV') == finalLevelEnergy]
                            if len(finalParticles) == 1:
                                finalParticle = finalParticles[0]
                            else:
                                # No exact match. First look for levels within .01% of the exact value:
                                idx = 0
                                while True:
                                    idx += 1
                                    finalParticleName = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex(
                                        residualName, idx)
                                    if not reactionSuite.hasParticle(finalParticleName):
                                        warningList.append(
                                            "MF12 final level energy (%s eV) doesn't match known level when decaying out of level %s " %
                                            (finalLevelEnergy, MTLO2))
                                        info.doRaise.append(warningList[-1])
                                    try:
                                        thisLevelEnergy = reactionSuite.PoPs[finalParticleName].energy.float('eV')
                                    except KeyError:
                                        raise Exception('Could not find nuclide of %s with desired energy level of %s.'
                                                        % (residualName, finalLevelEnergy))
                                    except:
                                        raise
                                    if abs(thisLevelEnergy - finalLevelEnergy) < 1e-4 * finalLevelEnergy:
                                        finalParticle = reactionSuite.PoPs[finalParticleName]
                                        break  # found it
                        else:
                            finalParticle = reactionSuite.PoPs[residualName]
                        gammaTransition = 1.
                        if len(MF12['branching']) > 2: gammaTransition = MF12['branching'][1]

                        if gammaTransition != 1: raise Exception('Fix me')
                        probability = probabilityModule.Double(info.PoPsLabel, MF12['branching'][0])

                        decayMode = decayDataModule.DecayMode(str(i1), miscDecaysModule.Mode.electroMagnetic)
                        decayMode.probability.add(probability)
                        _decay = decayDataModule.Decay(str(i1), decayDataModule.decayModesParticle)
                        _decay.products.add(productModule.Product(IDsPoPsModule.photon, IDsPoPsModule.photon))
                        _decay.products.add(productModule.Product(finalParticle.id, finalParticle.id))
                        decayMode.decayPath.add(_decay)
                        if MF12['LG'] == 2:  # internal conversion competes with gamma emission
                            from PoPs.decays import spectrum as spectrumModule
                            Pgamma = spectrumModule.Shell(MF12['branching'][1], label=spectrumModule.Shell.total,
                                                          unit='')
                            decayMode.photonEmissionProbabilities.add(Pgamma)
                        reactionSuite.PoPs[fullName].decayData.decayModes.add(decayMode)
                    except Exception as err:
                        warningList.append('raise somewhere in "for MF12 in MF12_LO2" loop: MT%d, %s' % (MT, err))
                        info.doRaise.append(warningList[-1])
            else:
                raise Exception("Could not determine base MT for MF=12's MT=%s" % MTLO2)

    if doCovariances:
        covarianceMFs = sorted(set([mf for mt in MTDatas.values() for mf in mt.keys() if mf > 30]))
        if covarianceMFs:
            info.logs.write('    Reading covariances (MFs %s)\n' % ','.join(map(str, covarianceMFs)))
        try:
            """ parse covariances. This also requires setting up links from data to covariances, so we
            must ensure the links are synchronized """

            MTdict = {}
            for reaction in (list(reactionSuite.reactions) + list(reactionSuite.sums.crossSectionSums) +
                             list(reactionSuite.productions) + list(reactionSuite.fissionComponents) +
                             list(reactionSuite.incompleteReactions)):
                MT = reaction.ENDF_MT
                if MT in MTdict:
                    MTdict[MT].append(reaction)
                else:
                    MTdict[MT] = [reaction]
            covarianceSuite, linkData = parseCovariances(
                info, MTDatas, MTdict, singleMTOnly=singleMTOnly, resonances=reactionSuite.resonances,
                verbose=verbose, formatVersion=info.formatVersion)
            if (len(covarianceSuite.covarianceSections) > 0) or (len(covarianceSuite.parameterCovariances) > 0):
                covarianceSuite.target = str(info.target)
                covarianceSuite.projectile = str(info.projectile)
                covarianceSuite.interaction = info.reactionSuite.interaction

                # Add same style as reactionSuite, but without documentation
                evaluated = info.evaluatedStyle.copy()
                evaluated.documentation.endfCompatible.body = None
                covarianceSuite.styles.add(evaluated)
                # covarianceSuite.removeExtraZeros() # disable for easier comparison to ENDF
            else:
                covarianceSuite = None
        except Exception as e:
            warningList.append("Couldn't parse covariances! Error message: %s" % e)
            info.doRaise.append(warningList[-1])
            covarianceSuite = None
            raise
    else:
        covarianceSuite = None

    for sgTable in info.extraFissionWidths:
        # remove 2nd fission width from resonance parameters. Must be done after processing MF32 covariances
        channels = sgTable.findAttributeInAncestry('channels')
        channels.pop(channels.labels()[-1])
        sgTable.removeColumn('fission width_2')

    def overridePoPsIfNecessary(AWRI, spinParity=None):
        """
        Resonance regions may use different particle properties than the rest of the evaluation.
        If so, create a new PoPs database overriding values as necessary.
        """
        pops = None
        massChanged = spinChanged = False

        if AWRI != info.massTracker.getMostCommonMassAWR(info.targetZA):
            massChanged = True
        if spinParity is not None:
            spinNow, parityNow = spinParity
            commonSpin, commonParity = info.particleSpins[info.target]
            if spinNow != commonSpin:
                spinChanged = True
            if parityNow is not None and parityNow != commonParity:
                spinChanged = True

        if massChanged or spinChanged:
            from PoPs import database as PoPsDatabaseModule
            from PoPs.families import nuclide as PoPsNuclideModule

            pops = PoPsDatabaseModule.Database("resolved resonances", version="1.0", formatVersion=info.formatVersion)
            target = info.PoPs[info.target]
            if isinstance(target, aliasPoPsModule.MetaStable):
                # metastable target. Mass goes with GS (adjusted for excitation energy), spin goes with actual level
                target = info.PoPs.final(target.pid)
                AWRI -= target.energy.float('amu*c**2') / info.massTracker.neutronMass
                groundState = target.isotope[0].copy()
                target = target.copy()
            else:
                target = target.copy()
                groundState = target

            if massChanged:
                groundState.mass.add(
                    massModule.Double(label=info.style, value=AWRI * info.massTracker.neutronMass, unit='amu')
                )
            else:
                groundState.mass.add(
                    massModule.Double(
                        label=info.style,
                        value=info.massTracker.getMostCommonMassAWR(info.targetZA) * info.massTracker.neutronMass,
                        unit='amu'))

            spin, parity = spinParity or info.particleSpins[info.target]
            target.nucleus.spin.add(spinModule.Fraction(label=info.style, value=spin, unit='hbar'))
            if parity:
                target.nucleus.parity.add(parityModule.Integer(label=info.style, value=parity, unit=''))

            pops.add(groundState)
            if target is not groundState:
                pops.add(target)

        return pops

    for key in info.PoPsOverrides:
        popsDB = overridePoPsIfNecessary(*info.PoPsOverrides[key])
        if popsDB is not None:
            key.PoPs = popsDB

    info.massTracker.useMostCommonAMUmasses()

    if info.level > 0:  # AWR is for isomer mass. Adjust info.ZAMasses to GS mass:
        groundStateMass = info.massTracker.getMassAMU(info.targetZA) - PQUModule.PQU(
            PQUModule.PQU_float.surmiseSignificantDigits(info.level), 'eV/c**2').getValueAs('amu')
        info.massTracker.addMassAMU(info.targetZA, groundStateMass)  # overwrite excited state mass

    for ZA in info.massTracker.amuMasses:
        if ZA in [1]: continue  # neutron already handled as special case
        mass = info.massTracker.amuMasses[ZA]
        elementSymbol = chemicalElementMiscModule.symbolFromZ[ZA // 1000]
        name = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA(elementSymbol, ZA % 1000)
        name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex(name, 0)
        mass = massModule.Double(info.PoPsLabel, mass, quantityModule.stringToPhysicalUnit('amu'))
        if name not in reactionSuite.PoPs: toGNDSMiscModule.getPoPsParticle(info, ZA, levelIndex=0)
        particle = reactionSuite.PoPs[name]
        particle.mass.add(mass)

    sys.stdout.flush()
    if verbose > 0:
        for warning in warningList: info.logs.write("       WARNING: %s\n" % warning, stderrWriting=True)

    targetID = reactionSuite.target
    if targetID in reactionSuite.PoPs.aliases: targetID = reactionSuite.PoPs[targetID].pid
    ignoreID = None
    for particleID, spinParity in info.particleSpins.items():
        if particleID == 'target': ignoreID = targetID
    for particleID, spinParity in info.particleSpins.items():
        if ignoreID == particleID: continue
        spin = spinModule.Fraction(info.PoPsLabel, fractions.Fraction(spinParity[0]), spinModule.baseUnit)
        if particleID == reactionSuite.target: particleID = targetID
        if particleID == 'target':
            particle = reactionSuite.PoPs[targetID]
        else:
            particle = reactionSuite.PoPs[particleID]

        if isinstance(particle, nuclideModule.Particle): particle = particle.nucleus

        if len(particle.spin) == 0: particle.spin.add(spin)
        if spinParity[1]:
            parity = spinParity[1].value
            particle.parity.add(parityModule.Integer(info.PoPsLabel, parity, parityModule.baseUnit))

    for reaction in reactionSuite.reactions:  # For two-body reactions, this sections decays an excited state
        if reaction.ENDF_MT == 2:
            continue
        if reaction.outputChannel.genre == enumsModule.Genre.twoBody:
            residual = reaction.outputChannel[1]
            if residual.outputChannel is None:
                particle = info.PoPs[residual.pid]
                if isinstance(particle, nuclideModule.Particle):
                    addGroundState = particle.nucleus.index != 0
                    if addGroundState:  # gamma data should be in orphanProducts.
                        groundState = particle.isotope.nuclides[0]
                        ZA = chemicalElementMiscModule.ZA(groundState)
                        multiplicity = residual.multiplicity[0].copy()
                        residual.addOutputChannel(outputChannelModule.OutputChannel(enumsModule.Genre.NBody))
                        residual.outputChannel.Q.add(
                            toGNDSMiscModule.returnConstantQ(
                                info.style, particle.nucleus.energy[0].value, multiplicity))
                        product = toGNDSMiscModule.newGNDSParticle(
                            info, toGNDSMiscModule.getPoPsParticle(info, ZA, groundState.id),
                            crossSection, multiplicity=multiplicity)
                        product.distribution.add(unspecifiedModule.Form(info.style, xDataEnumsModule.Frame.lab))
                        info.ENDFconversionFlags.add(product, 'implicitProduct')
                        residual.outputChannel.products.add(product)
                        reaction.updateLabel()

    if reconstructResonances and reactionSuite.resonances is not None and reactionSuite.resonances.reconstructCrossSection:
        info.logs.write('    Reconstructing resonances\n')
        try:
            reactionSuite.reconstructResonances(info.reconstructedStyle, info.reconstructedAccuracy, verbose=verbose,
                                                thin=True)
        except Exception as ex:
            warningList.append("Resonance reconstruction failed: %s" % ex)
            info.resonanceReconstructionFailed = str(ex)
            info.doRaise.append(warningList[-1])

    def adjustMF13Multiplicity2(multiplicity, crossSection):

        energyMultiplicity = []
        if multiplicity.domainMax > crossSection.domainMax:
            multiplicity = multiplicity.domainSlice(domainMax=crossSection.domainMax)
        for energyIn, multiplicityValue in multiplicity:
            crossSectionAtEnergy = crossSection.evaluate(energyIn)
            if crossSectionAtEnergy != 0: multiplicityValue /= crossSectionAtEnergy
            energyMultiplicity.append([energyIn, multiplicityValue])
        multiplicity.setData(energyMultiplicity)

    def adjustMF13Multiplicity(multiplicity, crossSection):

        if isinstance(multiplicity, multiplicityModule.XYs1d):
            adjustMF13Multiplicity2(multiplicity, crossSection)
        elif isinstance(multiplicity, multiplicityModule.Regions1d):
            for region in multiplicity: adjustMF13Multiplicity2(region, crossSection)
        else:
            raise Exception('Unsupported multiplicity type "%s"' % multiplicity.moniker)

    def adjustMF13Gammas(reaction):

        MT = reaction.ENDF_MT
        crossSection = None
        allproducts = list(reaction.outputChannel)
        for prod in reaction.outputChannel:
            if prod.outputChannel is not None:
                allproducts.extend(list(prod.outputChannel))

        for product in allproducts:
            multiplicity = product.multiplicity[info.style]
            if hasattr(multiplicity, '_temp_divideByCrossSection'):
                if crossSection is None:
                    try:
                        crossSection = reaction.crossSection.toPointwise_withLinearXYs(accuracy=1e-3, upperEps=1e-8)
                    except Exception as ex:
                        warningList.append("Could not get linear cross section to convert multiplicities: %s" % ex)
                        info.doRaise.append(warningList[-1])
                        continue
                adjustMF13Multiplicity(multiplicity, crossSection)
                del multiplicity._temp_divideByCrossSection

        if MT in info.totalMF6_12_13Gammas:
            MF, multiplicity = info.totalMF6_12_13Gammas[MT]
            if MF == 13: adjustMF13Multiplicity(multiplicity, crossSection)
            gammaProduction = [tmp for tmp in reactionSuite.reactions if tmp.ENDF_MT == MT]
            gammaProduction += [tmp for tmp in reactionSuite.orphanProducts if tmp.ENDF_MT == MT]
            if len(gammaProduction) != 1: raise ValueError("No unique match found.")
            gammaProduction = gammaProduction[0]
            summands = [sumsModule.Add(link=r.multiplicity) for r in
                        gammaProduction.outputChannel.getProductsWithName(IDsPoPsModule.photon)]
            if len(summands) == 0:
                for _product in gammaProduction.outputChannel:
                    if _product.outputChannel is not None:
                        summands += [sumsModule.Add(link=r.multiplicity) for r in
                                     _product.outputChannel.getProductsWithName(IDsPoPsModule.photon)]
            if MT in channelIDs:
                name = channelIDs[MT]
            else:
                name = gammaProduction.outputChannel.toString(MT=MT)

            multiplicitySum = sumsModule.MultiplicitySum(label=name + " total gamma multiplicity", ENDF_MT=MT)
            for summand in summands: multiplicitySum.summands.append(summand)
            multiplicitySum.multiplicity.add(multiplicity)
            reactionSuite.sums.multiplicitySums.add(multiplicitySum)

    for reaction in reactionSuite.reactions: adjustMF13Gammas(reaction)
    for reaction in reactionSuite.orphanProducts: adjustMF13Gammas(reaction)

    # Short-lived light isotopes sometimes listed implicitly in ENDF-6, without associated mass:
    for pid, ZA in (('He5', 2005), ('Li5', 3005), ('Be8', 4008)):
        if pid in reactionSuite.PoPs:
            particle = reactionSuite.PoPs[pid]
            try:
                particle.getMass('amu')
            except:
                from brownies.legacy.endl.structure import masses as AMEmasses
                particle.mass.add(massModule.Double(info.PoPsLabel, AMEmasses.getMassFromZA(ZA), 'amu'))

    for QI, outputChannel in SpecialLRProducts:
        Q = info.PoPs[info.projectile].getMass('eV/c**2') + info.PoPs[info.target].getMass('eV/c**2')
        Q -= info.PoPs[outputChannel[0].pid].getMass('eV/c**2')
        residual = info.PoPs[outputChannel[1].pid]
        residualGroundState = info.PoPs[info.PoPs[outputChannel[1].pid].isotope.symbol]
        Q -= residualGroundState.getMass('eV/c**2')
        residual.energy[0].value = Q - QI

    return covarianceSuite
