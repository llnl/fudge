# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from PoPs import database as databasePoPsModule
from PoPs import alias as PoPsAliasModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import mass as massModule
from PoPs.families import nuclide as nuclideModule
from PoPs.chemicalElements import chemicalElement as chemicalElementModule
from PoPs.chemicalElements import misc as chemicalElementMiscModule

from xData import enums as xDataEnumsModule
import fudge
from fudge import enums as enumsModule
from fudge import institution as institutionModule
from fudge import physicalQuantity as physicalQuantityModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule
from fudge.productData.distributions import unspecified as unspecifiedModule

from .ENDFToGNDS import endfFileToGNDSMisc, ENDF_ITYPE_0
from .ENDFToGNDS import ENDF_ITYPE_2, ENDF_ITYPE_3, ENDF_ITYPE_4, ENDF_ITYPE_6
from .ENDFToGNDS import ENDF_ITYPE_1 as ENDF_ITYPE_1Module
from .ENDFToGNDS import ENDF_ITYPE_5 as ENDF_ITYPE_5Module
from . import endf_endl as endf_endlModule
from . import toGNDSMisc

from brownies.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule


class LogFiles:

    def __init__(self, toStdOut=False, toStdErr=False, logFile=None, defaultIsStderrWriting=True):

        def addIfDifferent(newOutput):

            if newOutput.isatty():
                for output in self.outputs:
                    if output.isatty(): return False
            return True

        self.defaultIsStderrWriting = defaultIsStderrWriting
        self.outputs = []
        if toStdOut: self.outputs.append(sys.stdout)
        if toStdErr and addIfDifferent(sys.stderr): self.outputs.append(sys.stderr)
        if logFile:
            if isinstance(logFile, str):
                logFile = open(logFile, 'w')
            elif type(logFile) != type(sys.stdout):
                raise Exception('logFile type = %s is no supported' % (type(logFile)))
            if addIfDifferent(logFile): self.outputs.append(logFile)

    def write(self, msg, stderrWriting=None):

        if stderrWriting is None: stderrWriting = self.defaultIsStderrWriting
        for output in self.outputs:
            if (output == sys.stderr) and not stderrWriting: continue
            output.write(msg)


def readMF1MT451(_MAT, _MTDatas, formatVersion, specialNuclearParticleID, styleName='eval', logFile=None,
                 verboseWarnings=False, **kwargs):
    # We're going to save everything in the info instance
    info = toGNDSMisc.Infos(formatVersion, styleName)
    info.doRaise = []
    info.logs = logFile
    info.verboseWarnings = verboseWarnings
    info.ENDFconversionFlags = ENDFconversionFlagsModule.ENDFconversionFlags()
    info.specialNuclearParticleID = specialNuclearParticleID
    info.extraFissionWidths = []

    # Line #1
    targetZA, targetMass, LRP, LFI, NLIB, NMOD = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(_MTDatas[451][1][0], logFile=logFile)
    targetZA = int(targetZA)  # Target's ZA
    LRP = int(LRP)  # Resonance parameter data info
    LFI = int(LFI)  # Is fission present
    NLIB = int(NLIB)  # What library (e.g., 0 = ENDF/B
    NMOD = int(NMOD)  # Version modification flag
    isNaturalTarget = (targetZA % 1000) == 0

    # Line #2
    targetExcitationEnergy, STA, LIS, LISO, dummy, NFOR = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(_MTDatas[451][1][1], logFile=logFile)
    STA = int(STA)  # Is nucleus unstable
    LIS = int(LIS)  # Excitation number
    LISO = int(LISO)  # Isomeric state number
    if (targetExcitationEnergy != 0) or (LIS != 0) or (LISO != 0):
        if (targetExcitationEnergy == 0) or (LIS == 0) or (LISO == 0):
            print("    WARNING: target excitation energy = %s, LIS = %s and LISO = %s are not consistent." % (
                targetExcitationEnergy, LIS, LISO))
    NFOR = int(NFOR)  # Must be 6 for ENDF-6 format
    if NFOR != 6:
        print("    WARNING: endfFileToGNDS only supports ENDF-6 format. This file has unsupported NFOR=%d" % NFOR)

    info.ZA_massLineInfo.targetZA = targetZA
    info.ZA_massLineInfo.targetLIS = LIS
    info.ZA_massLineInfo.add(targetZA, targetMass, 451, 1, 0)
    # Line #3
    projectileMass, EMAX, LREL, dummy, NSUB, NVER = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(_MTDatas[451][1][2], logFile=logFile)
    NSUB = int(NSUB)  # 10 * ZA + iType for projectile
    NVER = int(NVER)  # Evaluation version number
    LREL = int(LREL)  # Evaluation sub-version number
    IPART, ITYPE = divmod(NSUB, 10)
    projectileZA = IPART
    if projectileZA == 11: projectileZA = 9
    info.ZA_massLineInfo.add(projectileZA, projectileMass, 451, 1, 2)

    # Line #4
    targetTemperature, dummy, LDRZ, dummy, NWD, NXC = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(_MTDatas[451][1][3], logFile=logFile)
    LDRZ = int(LDRZ)  # Primary or special evaluation of this material
    NWD = int(NWD)  #
    NXC = int(NXC)  #

    info.convertJENDL_stylePrimarygammas = kwargs.get('JENDL_stylePrimarygammas', NLIB == 6)
    # If True, treats MF=6 primary gamma energies as binding energy, otherwise as gamma energy.

    # Save the library name and version
    info.library = endf_endlModule.NLIBs.get(NLIB, 'Unknown')
    info.evaluation = info.library  # why do I need this?
    info.libraryVersion = "%d.%d.%d" % (NVER, LREL, NMOD)

    # Save the evaluation date
    try:
        info.Date = endfFileToGNDSMisc.getENDFDate(_MTDatas[451][1][4][22:33])
    except Exception as e:
        if not kwargs.get('ignoreBadDate'):
            info.doRaise.append(str(e))
        import datetime
        info.Date = datetime.datetime.today().strftime("%Y-%m-%d")

    # Save general information
    info.author = _MTDatas[451][1][4][33:66]
    info.documentation = '\n'.join(_MTDatas[451][1][4:4 + NWD])
    info.MAT = _MAT
    info.LISO = LISO
    info.LIS = LIS
    info.EMAX = EMAX
    info.levelIndex = LIS
    info.NSUB = NSUB
    info.LFI = LFI
    info.STA = STA
    info.targetTemperature = targetTemperature
    info.materialName = _MTDatas[451][1][4][:11].strip()
    info.targetMass = targetMass
    info.isNaturalTarget = isNaturalTarget
    info.LRP = LRP
    info.NLIB = NLIB
    info.NMOD = NMOD
    info.NVER = NVER
    info.LREL = LREL
    info.NXC = NXC
    info.NWD = NWD
    info.LDRZ = LDRZ
    info.ITYPE = ITYPE

    # All evaluations need a particle database and we fill in a lot of that information from MF1MT451
    info.PoPsLabel = styleName
    info.PoPs = databasePoPsModule.Database('protare_internal', '1.0', formatVersion=formatVersion)

    # Let's tell info about the target
    info.targetZA = targetZA
    info.targetLevel = LIS
    info.level = targetExcitationEnergy
    info.addMassAWR(info.targetZA, targetMass)

    # If there is a projectile, tell info about it too
    if ITYPE in [4, 5]:
        info.projectile = None
        info.projectileZA = None
    else:
        name = {0: 'g', 1: IDsPoPsModule.neutron, 11: 'e-',
                1001: 'H1', 1002: 'H2', 1003: 'H3', 2003: 'He3', 2004: 'He4'}[IPART]
        info.projectile = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name, specialNuclearParticleID)
        info.projectileZA = projectileZA
        info.addMassAWR(projectileZA, projectileMass, asTarget=False)

    return info


def printMassHistoryAndReturn(printMassHistory, info, stuffToReturn):
    if printMassHistory:
        info.massTracker.printHistory()
        info.ZA_massLineInfo.printInfo()

    return stuffToReturn


def endfFileToGNDS(fileName, useFilesQAlways=True, singleMTOnly=None, evaluation=None,
                   MTs2Skip=None, parseCrossSectionOnly=False, formatVersion=GNDS_formatVersionModule.default,
                   toStdOut=True, toStdErr=True, logFile=None, skipBadData=False, doCovariances=True,
                   verboseWarnings=False, verbose=1, reconstructResonances=True,
                   specialNuclearParticleID=specialNuclearParticleIDPoPsModule.Mode.nuclide, **kwargs):
    logs = LogFiles(toStdOut=toStdOut, toStdErr=toStdErr, logFile=logFile, defaultIsStderrWriting=False)
    header, MAT, MTDatas = endfFileToGNDSMisc.parseENDFByMT_MF(fileName, logFile=logs)

    styleName = 'eval'
    reconstructedStyleName = 'recon'

    if MTs2Skip is None: MTs2Skip = []

    # Parse the ENDF documentation section
    info = readMF1MT451(MAT, MTDatas, formatVersion, specialNuclearParticleID, styleName=styleName, logFile=logs,
                        verboseWarnings=verboseWarnings, **kwargs)

    info.missingRadioactiveProduct = []
    info.printBadNK14 = True
    info.continuumSpectraFix = False
    info.acceptBadMF10FissionZAP = False
    options = ['printBadNK14', 'continuumSpectraFix', 'ignoreBadDate', 'acceptBadMF10FissionZAP',
               'JENDL_stylePrimarygammas', 'printMassHistory']
    for option in kwargs:
        if option not in options: raise DeprecationWarning('invalid or deprecated option "%s"' % option)
        setattr(info, option, kwargs[option])

    printMassHistory = kwargs.get('printMassHistory', False)

    if info.ITYPE == 1:  # NFY
        return printMassHistoryAndReturn(printMassHistory, info,
                                         ENDF_ITYPE_1Module.ITYPE_1(MTDatas, info, verbose=verbose))
    elif info.ITYPE == 2:  # Thermal scattering law data
        return printMassHistoryAndReturn(printMassHistory, info,
                                         ENDF_ITYPE_2.ITYPE_2(fileName, MAT, MTDatas, info, evaluation,
                                                              verbose=verbose))
    elif info.ITYPE == 4:  # Decay data
        return printMassHistoryAndReturn(printMassHistory, info, ENDF_ITYPE_4.ITYPE_4(MTDatas, info, verbose=verbose))
    elif info.ITYPE == 5:  # SFY
        return printMassHistoryAndReturn(printMassHistory, info,
                                         ENDF_ITYPE_5Module.ITYPE_5(MTDatas, info, verbose=verbose))
    elif info.ITYPE == 6:  # Atomic relaxation data
        if info.NSUB not in [6]:
            raise ValueError('For ITYPE = %d, invalid NSUB = %s' % (info.ITYPE, info.NSUB))
        return printMassHistoryAndReturn(
            printMassHistory, info, ENDF_ITYPE_6.ITYPE_6(info.targetZA // 1000, MTDatas, info, verbose=verbose))

    # Second batch of ITYPE's are general transport data
    if info.LIS != 0: levelIndex = info.LIS
    if info.ITYPE in [0, 1, 9]:  # Nuclear transport data (g, n, p, d, t, h, a), including FPY
        target = toGNDSMisc.getTypeNameGamma(info, info.targetZA, level=info.level, levelIndex=info.levelIndex)
        targetID = target.id
        if info.STA != 0:
            halflife = halflifeModule.String(info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit)
            target = info.PoPs[targetID]
            if isinstance(target, nuclideModule.Particle): target = target.nucleus
            if len(target.halflife) == 0: target.halflife.add(halflife)
    elif info.ITYPE == 3:  # Atomic transport data (e, x)
        # target is a chemicalElement, but target mass is stored in an A=0 nuclide
        targetZ = info.targetZA // 1000
        targetID = chemicalElementMiscModule.symbolFromZ[targetZ]
        target = nuclideModule.Particle.buildFromClassAndRawData(targetID + "0",
                mass=(info.massTracker.getMassAMU(info.targetZA), 'amu'),
                label=info.PoPsLabel
                )
        info.PoPs.add(target)
    else:
        raise ValueError("Unsupported ITYPE = %s" % info.ITYPE)

    unrecognizedMAT = False
    try:
        ZA2, MAT2 = endf_endlModule.ZAAndMATFromParticleName(targetID)
        MAT2 += info.LISO
        if MAT2 != info.MAT:
            info.logs.write("       WARNING: ENDF MAT = %s not as expected (i.e., %s).\n" % (info.MAT, MAT2),
                            stderrWriting=True)
            unrecognizedMAT = True
    except KeyError:
        info.logs.write("        WARNING: could not compute MAT for particle %s.\n" % targetID)
        unrecognizedMAT = True
    info.targetID = targetID
    info.targetAliasID = targetID
    # Add alias if needed for metastable target
    if info.LISO != 0:
        aliasName = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(targetID, info.LISO)
        info.PoPs.add(PoPsAliasModule.MetaStable(aliasName, targetID, info.LISO))
        info.targetAliasID = aliasName

    # Compute the reconstructed and evaluated Styles
    projectileDomain = stylesModule.ProjectileEnergyDomain(
        1e-5, 2e+7, 'eV')  # will be overwritten after parsing all reactions
    evaluatedStyle = stylesModule.Evaluated(
        styleName, '', physicalQuantityModule.Temperature(
            PQUModule.PQU_float.surmiseSignificantDigits(info.targetTemperature), 'K'),
        projectileDomain, info.library, info.libraryVersion, date=info.Date)
    if evaluation is None: evaluation = "%s-%d.%d" % (info.library, info.NVER, info.LREL)
    info.reconstructedStyle = stylesModule.CrossSectionReconstructed(
        reconstructedStyleName, derivedFrom=evaluatedStyle.label, date="2016-11-06")

    # Stuff for handling transport data
    evaluatedStyle.documentation.endfCompatible.body = info.documentation
    endfFileToGNDSMisc.completeDocumentation(info, evaluatedStyle.documentation)
    info.evaluatedStyle = evaluatedStyle
    info.reconstructedAccuracy = 0.001
    info.MF12_LO2 = {}
    info.AWR_mode = None
    info.particleSpins = {}
    info.missingTwoBodyMasses = {}
    info.MF4ForNonNeutrons = []
    projectile = toGNDSMisc.getTypeNameGamma(info, info.projectileZA)
    projectileID = projectile.id
    if projectile.id != IDsPoPsModule.photon:
        projectileID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(
            projectile.id, specialNuclearParticleID)

    # Set up the reactionSuite
    interaction = enumsModule.Interaction.nuclear
    if info.ITYPE == 3: interaction = enumsModule.Interaction.atomic
    reactionSuite = reactionSuiteModule.ReactionSuite(projectileID, info.targetAliasID, evaluation,
                                                      interaction=interaction,
                                                      style=info.evaluatedStyle, MAT=MAT, PoPs=info.PoPs,
                                                      formatVersion=info.formatVersion)
    MTDatas[451][1] = MTDatas[451][1][:4 + info.NWD]
    info.reactionSuite = reactionSuite
    info.PoPs = reactionSuite.PoPs
    info.target = reactionSuite.target
    if unrecognizedMAT:
        info.ENDFconversionFlags.add(reactionSuite, "MAT=%d" % MAT)

    # Set up the covarianceSuite, if needed, and do other ITYPE-specific calculations
    covarianceSuite = None
    if info.ITYPE in [0, 9]:
        doRaise = info.NSUB not in {0: [0, 10, 10010, 10020, 10030, 20030, 20040], 9: [19]}[info.ITYPE]
        if doRaise: raise ValueError('For ITYPE = %d, invalid NSUB = %s' % (info.ITYPE, info.NSUB))
        covarianceSuite = ENDF_ITYPE_0.ITYPE_0(MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip,
                                               parseCrossSectionOnly, doCovariances, verbose,
                                               reconstructResonances=reconstructResonances)
    elif info.ITYPE in [3]:
        if info.NSUB not in [3, 11, 113]:
            raise ValueError('For ITYPE = %d, invalid NSUB = %s' % (info.ITYPE, info.NSUB))
        ENDF_ITYPE_3.ITYPE_3(MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose=verbose)

    # Fix incident energy domain limits:
    reactionList = reactionSuite.reactions
    if len(reactionList) == 0:
        reactionList = list(reactionSuite)

    if reactionList[0].ENDF_MT == 2:
        nonThresholdList = [reactionList[0]]
    else:
        nonThresholdList = []
        for reac in reactionList:
            if hasattr(reac, 'outputChannel'):
                if reac.ENDF_MT == 528:
                    # electron excitation reaction claims Q=0 but still has a threshold.
                    # Ignore it when determining projectileEnergyDomain
                    continue
                if reac.outputChannel.Q.getConstant() >= 0:
                    nonThresholdList.append(reac)
            elif hasattr(reac, 'Q'):
                print(type(reac.Q))
                print(len(reac.Q))
                print(type(reac.Q[0]))
                if reac.Q.getConstant() >= 0:
                    nonThresholdList.append(reac)
            else:
                pass  # print warning?

    if len(nonThresholdList) > 0:
        domainMin = max([reac.crossSection.domainMin for reac in nonThresholdList])
    else:
        domainMin = min([reac.crossSection.domainMin for reac in reactionList])
    domainMax = min([reac.crossSection.domainMax for reac in reactionList])

    info.evaluatedStyle.projectileEnergyDomain.min = domainMin
    info.evaluatedStyle.projectileEnergyDomain.max = domainMax

    if covarianceSuite is not None:
        # FIXME what if covariances span smaller domain than reactions?
        covarianceSuite.styles[info.style].projectileEnergyDomain.min = domainMin
        covarianceSuite.styles[info.style].projectileEnergyDomain.max = domainMax
        reactionSuite._loadedCovariances = [covarianceSuite]

    # Fill up the reactionSuite
    for reaction in reactionSuite.reactions: addUnspecifiedDistributions(info, reaction.outputChannel)
    for production in reactionSuite.productions: addUnspecifiedDistributions(info, production.outputChannel)

    if len(info.ENDFconversionFlags.flags) > 0:
        flagsToPurge = []
        for conversion in info.ENDFconversionFlags.flags:
            if conversion.link.rootAncestor is covarianceSuite:
                conversion.root = '$covariances'  # FIXME need better way of handling this
            elif conversion.link.rootAncestor is not reactionSuite:
                flagsToPurge.append(conversion)  # flag doesn't point to any GNDS node
        for flag in flagsToPurge:
            info.ENDFconversionFlags.flags.remove(flag)
        info.ENDFconversionFlags.flags.sort(key=lambda val: val.flags)
        LLNLdata = institutionModule.Institution("LLNL")
        LLNLdata.append(info.ENDFconversionFlags)
        reactionSuite.applicationData.add(LLNLdata)

    PoPs = reactionSuite.PoPs

    lightMasses = {  # atomic masses in amu, taken from AME2012 (http://www.nndc.bnl.gov/masses/mass.mas12)
        IDsPoPsModule.neutron: 1.00866491585,
        'H1': 1.00782503223,
        'H2': 2.01410177812,
        'H3': 3.01604927791,
        'He3': 3.01602932008,
        'He4': 4.00260325413,
    }
    for pid in lightMasses:
        if pid in info.PoPs:
            if len(info.PoPs[pid].mass) == 0:
                info.PoPs[pid].mass.add(
                    massModule.Double(info.PoPsLabel, lightMasses[pid], quantityModule.stringToPhysicalUnit('amu'))
                )

    for product1ID in info.missingTwoBodyMasses:
        if len(info.PoPs[product1ID].mass) != 0:  # mass was added to PoPs, likely in previous block
            product1Mass = info.PoPs[product1ID].getMass("amu")
        else:
            raise ValueError("Missing mass for 2-body reaction product %s" % product1ID)

        massPTmP1 = PoPs[reactionSuite.projectile].getMass('amu') + PoPs[targetID].getMass('amu') - \
                    product1Mass
        summedQM = 0
        for product2ID, QM in info.missingTwoBodyMasses[product1ID]: summedQM += QM
        mass = massPTmP1 - PQUModule.PQU(summedQM / len(info.missingTwoBodyMasses[product1ID]), 'eV/c**2').getValueAs(
            'amu')
        mass = massModule.Double(info.PoPsLabel, mass, quantityModule.stringToPhysicalUnit('amu'))
        residual = reactionSuite.PoPs[info.missingTwoBodyMasses[product1ID][0][0]]
        if len(residual.mass) == 0: residual.mass.add(mass)

    # Did we have trouble?
    if len(info.MF4ForNonNeutrons) > 0 and verbose > 0:
        print('    WARNING: MF=4 data for non-neutron product: MTs = %s' % sorted(info.MF4ForNonNeutrons))

    if len(info.doRaise) > 0 and not skipBadData:
        info.logs.write('\nRaising due to following errors:\n')
        for err in info.doRaise: info.logs.write('    ' + err + '\n', True)
        raise Exception('len( info.doRaise ) > 0')
    if info.resonanceReconstructionFailed != '':
        print('ERROR: Resonance reconstruction failed: "%s".' % info.resonanceReconstructionFailed)

    if len(info.missingRadioactiveProduct) > 0:
        print('Missing radioactive products:')
        for missingRadioactiveProduct in info.missingRadioactiveProduct:
            print('     ========', missingRadioactiveProduct)

    printMassHistory = kwargs.get('printMassHistory', False)
    if printMassHistory:
        info.massTracker.printHistory()
        info.ZA_massLineInfo.printInfo()

    return {'reactionSuite': reactionSuite, 'covarianceSuite': covarianceSuite, 'errors': info.doRaise, 'info': info}


def addUnspecifiedDistributions(info, outputChannel):
    if outputChannel is None: return
    for product in outputChannel:
        if len(product.distribution) == 0:
            product.distribution.add(unspecifiedModule.Form(info.style, productFrame=xDataEnumsModule.Frame.lab))
            info.ENDFconversionFlags.add(product, 'implicitProduct')
        addUnspecifiedDistributions(info, product.outputChannel)


if __name__ == '__main__':

    rce = endfFileToGNDS(fileName=sys.argv[1])
    x, c = rce['reactionSuite'], rce['covarianceSuite']
    f = open('test.xml', 'w')
    f.write('\n'.join(x.toXML_strList() + ['']))
    f.close()
    if c is not None:  # covariances
        f = open('test-covar.xml', 'w')
        f.write('\n'.join(c.toXML_strList() + ['']))
        f.close()
