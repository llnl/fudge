# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
ENDF ITYPE=2 files store thermal neutron scattering data.
This module is for translating ITYPE=2 from ENDF-6 into GNDS.

These methods are meant to be called from endfFileToGNDS.
"""

import os
import numpy

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import xDataArray as arrayModule
from xData import gridded as griddedModule
from xData import values as valuesModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as IDsPoPsModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.families import unorthodox as unorthodoxModule
from PoPs.chemicalElements import chemicalElement as chemicalElementModule

from fudge import enums as enumsModule, targetInfo
from fudge import physicalQuantity as physicalQuantityModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

from fudge.reactions import reaction as reactionModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import (
    coherentElastic as coherentElasticModule,
    incoherentElastic as incoherentElasticModule,
    incoherentInelastic as incoherentInelasticModule)
import fudge.productData.distributions.reference as referenceModule

from fudge import institution as institutionModule
from fudge import outputChannel as outputChannelModule
from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

from .. import endf_endl as endf_endlModule
from .. import toGNDSMisc as toGNDSMiscModule


def ITYPE_2(fileName, MAT, MTDatas, info, evaluation, verbose):
    """Parses ENDF Thermal neutron scattering law data."""

    errors = []
    EMin = 1e-5
    EMax = info.EMAX
    if EMax <= EMin:
        print('    WARNING: EMax = %s, setting to 5.0 eV' % EMax)
        EMax = 5.0
    if EMax > 5:
        print('    WARNING: EMax = %s is unphysical, setting to 5.0 eV' % EMax)
        EMax = 5.0

    projectile = toGNDSMiscModule.getTypeNameGamma(info, 1)
    neutronMass = projectile.getMass('amu')

    if evaluation is None: evaluation = "%s-%d.%d" % (info.library, info.NVER, info.LREL)

    projectileDomain = stylesModule.ProjectileEnergyDomain(EMin, EMax, 'eV')
    axes = crossSectionModule.defaultAxes('eV')
    domainXYs1d = XYs1dModule.XYs1d([[EMin, 1.0], [EMax, 1.0]], axes=axes)

    temperature = physicalQuantityModule.Temperature(
        PQUModule.PQU_float.surmiseSignificantDigits(info.targetTemperature), 'K')
    evaluatedStyle = stylesModule.Evaluated(info.PoPsLabel, '', temperature, projectileDomain,
                                            info.library, info.libraryVersion, date=info.Date)

    # FIXME next line assumes file starts with 'tsl-' or 'tsl_'
    fileBaseName = os.path.splitext(os.path.basename(fileName))[0][4:]
    targetID = fileBaseName
    scatterers = []

    # if multiple scattering atoms are present, MF=1 mass *should* correspond to the primary scatterer
    # however, some evaluations (e.g. benzene and CH4 in ENDF-VIII) enter the mass of the entire molecule instead!
    if fileBaseName[0].isdigit() and fileBaseName.count('_') == 2:
        Z, Symbol, A = fileBaseName.split('_')
        ZA = int(Z) * 1000 + int(A)
        targetID = 'tnsl-%s%d' % (Symbol, int(A))
    elif 'ice' in fileBaseName.lower():
        ZA = 1001
        scatterers = [('H', 2), ('O', 1)]
        if 'Oin' in fileBaseName:
            ZA = 8016
    elif 'Paraffin' in fileBaseName:
        ZA = 1001
        scatterers = [('H', 11), ('C', 15), ('C', 1), ('I', 1), ('O', 7)]
    elif 'para' in fileBaseName or 'ortho' in fileBaseName:
        ZA = 1001
        if fileBaseName[-1] == 'D': ZA = 1002
    elif 'graphite' in fileBaseName.lower():
        ZA = 6012
    elif 'Be-metal' in fileBaseName:
        ZA = 4009
    elif 'SiO2' in fileBaseName:
        ZA = 14028
        scatterers = [('Si', 1), ('O', 2)]
        if 'Oin' in fileBaseName:
            ZA = 8016
    elif 'Uin' in fileBaseName and 'HEU' in fileBaseName:
        ZA = 92235
    elif 'U-metal' in fileBaseName:
        ZA = 92238
        if 'HEU' in fileBaseName:
            ZA = 92235
    elif 'benzene' in fileBaseName.lower() or 'benzine' in fileBaseName.lower():  # ENDF-VII.1 had spelling error
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 6), ('C', 6)]
    elif 'ethanol' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 6), ('C', 2), ('O', 1)]
    elif 'triphenylmethane' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 16), ('C', 19)]
    elif 'CH4' in fileBaseName or 'methane' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 4), ('C', 1)]
    elif 'H2inCaH2' in fileBaseName:
        ZA = 1001  # misleading name: H2 refers to 2nd hydrogen atom, not deuterium
    elif 'toluene' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 8), ('C', 7)]
    elif 'xylene' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 10), ('C', 8)]
    elif 'mesitylene' in fileBaseName.lower():
        ZA = 1001
        if fileBaseName.startswith('Cin'):
            ZA = 6012
        scatterers = [('H', 12), ('C', 9)]
    elif 'in' in fileBaseName:
        targetSymbol, molecule = fileBaseName.split('in')
        ZA = {'H': 1001, 'H1': 1001, 'D': 1002, 'H2': 1002, 'Li': 3007, '7Li': 3007, 'Be': 4009, 'C': 6012, 'N': 7014,
              'O': 8016, 'O16': 8016, 'F': 9019, 'Mg': 12024, 'Al': 13027, 'Al27': 13027, 'Si': 14028, 'Ca': 20040,
              'Y': 39089, 'Zr': 40090, 'U': 92238, 'Pu': 94239}[targetSymbol]

        # try extract scattering atoms & number per molecule from molecule name:
        import re
        patt = re.compile("([A-Z][a-z]?)([1-9]?)")
        idx = 0
        while True:
            match = patt.match(molecule, idx)
            if not match: break
            idx = match.end()
            symbol, count = match.groups()
            if not count: count = 1
            scatterers.append((symbol, int(count)))
    elif 'BeO' in fileBaseName:  # ENDF-VII.0 and earlier didn't break up into Be and O components
        ZA = 4009
    elif 'Mg' in fileBaseName:
        targetID = "tnsl-Mg"
        ZA = 12024  # Mg24 is 79% of natural abundance
    elif 'Si' in fileBaseName:
        targetID = "tnsl-Si"
        ZA = 14028  # Si28 is 92% of natural abundance
    elif 'V-metal' in fileBaseName:
        targetID = "tnsl-V-metal"
        ZA = 23051  # V51, 99.8% of natural abundance
    elif 'Ni' in fileBaseName:
        targetID = "tnsl-Ni"
        ZA = 28058  # Ni58, 68% of natural abundance
    elif 'Cu' in fileBaseName:
        targetID = "tnsl-Cu"
        ZA = 29063  # Cu63 is 69% of natural abundance
    elif 'Mo' in fileBaseName:
        targetID = "tnsl-Mo"
        ZA = 42098  # Mo98, 24% of natural abundance
    elif 'W-metal' in fileBaseName:
        targetID = "tnsl-W-metal"
        ZA = 74184  # W184, 30% of natural abundance
    elif 'Pb' in fileBaseName:
        targetID = "tnsl-Pb"
        ZA = 82208  # Pb208, 52% of natural abundance
    else:
        # FIXME try using the MAT number or ZSYMAM to get targetID?
        # Unfortunately libraries use different conventions for TNSL MAT numbers, and ZSYMAM can only be 11 chars
        raise Exception("Cannot determine TNSL target ZA or id from filename %s" % fileName)

    targetInfo_ = None
    info.boundAtomCrossSections = []
    if 7 in MTDatas[451]:
        # ENDF-6 now (as of 2023) supports listing isotopic abundances within TNSL. Parse those into 'targetInfo'
        dat = MTDatas[451][7]
        targetInfo_, info.boundAtomCrossSections = readMF7_info(info, dat)

    target = toGNDSMiscModule.getTypeNameGamma(info, ZA)
    if targetID is None: targetID = "tnsl-" + target.id
    tnsl_target = unorthodoxModule.Particle(targetID)
    tnsl_target.mass.add(
        massModule.Double(info.PoPsLabel, info.targetMass * neutronMass,
                          quantityModule.stringToPhysicalUnit('amu'))
    )
    info.PoPs.add(tnsl_target)

    interaction = enumsModule.Interaction.TNSL
    info.principalScattererPid = target.id
    if len(info.boundAtomCrossSections) > 1:
        # switch principal scattering atom to chemical element:
        popsParticle = info.PoPs[target.id]
        info.principalScattererPid = popsParticle.chemicalElementSymbol
    info.scatterers = scatterers
    evaluatedStyle.documentation.endfCompatible.body = info.documentation
    endfFileToGNDSMisc.completeDocumentation(info, evaluatedStyle.documentation)

    reactions = []
    for MT in [2, 4]:
        if MT in MTDatas:  # elastic
            forms = readMF7(info, MT, MTDatas[MT][7])

            for form in forms:
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.twoBody, process=form.process)
                outputChannel.Q.add(toGNDSMiscModule.returnConstantQ(info.style, 0.0, domainXYs1d))

                product = toGNDSMiscModule.newGNDSParticle(info, projectile, domainXYs1d)
                outputChannel.products.add(outputChannel.products.uniqueLabel(product))

                reaction = reactionModule.Reaction(None, outputChannel.genre, ENDF_MT=MT)
                endf_endlModule.setReactionsOutputChannelFromOutputChannel(info, reaction, outputChannel)
                reaction.doubleDifferentialCrossSection.add(form)
                reactions.append(reaction)

                reaction.crossSection.add(
                    crossSectionModule.ThermalNeutronScatteringLaw1d(link=form, label=info.style, relative=True))
                product.distribution.add(
                    referenceModule.ThermalNeutronScatteringLaw(link=form, label=info.style, relative=True))

    reactionSuite = reactionSuiteModule.ReactionSuite(info.projectile, targetID, evaluation, style=evaluatedStyle,
                                                      interaction=interaction,
                                                      formatVersion=info.formatVersion, MAT=MAT, PoPs=info.PoPs)
    if targetInfo_ is not None:
        reactionSuite.styles.getEvaluatedStyle().targetInfo = targetInfo_

    for reaction in reactions:
        reactionSuite.reactions.add(reaction)

    # add applicationData section for help converting back to ENDF-6
    info.ENDFconversionFlags.add(reactionSuite, 'MAT=%d,ZA=%d' % (MAT, ZA))
    LLNLdata = institutionModule.Institution("LLNL")
    LLNLdata.append(info.ENDFconversionFlags)
    reactionSuite.applicationData.add(LLNLdata)

    if len(info.doRaise) > 0:
        info.logs.write('\nRaising due to following errors:\n')
        for err in info.doRaise: info.logs.write('    ' + err + '\n', True)
        raise Exception('len( info.doRaise ) > 0')

    return {'reactionSuite': reactionSuite, 'errors': errors, 'covarianceSuite': None}


def readSTable(line, MF7):
    line, dat = endfFileToGNDSMisc.getTAB1(line, MF7)
    temps, betas = [dat['C1']], [dat['C2']]

    # energy interpolation should be 'flat':
    e_interp = dat['interpolationInfo']
    if len(e_interp) > 1: raise ValueError("only one interpolation region allowed for S table")
    e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(e_interp[0][1])
    energies, data = map(list, zip(*dat['data']))

    # higher temperatures:
    ntemps = int(dat['L1'])
    t_interp = []
    for j in range(ntemps):
        line, dat = endfFileToGNDSMisc.getList(line, MF7)
        temps.append(dat['C1'])
        betas.append(dat['C2'])
        data.extend(dat['data'])
        if dat['L1'] == 0:
            # JENDL-5 issue, using unsupported temperature interpolation.
            dat['L1'] = 2
        t_interp.append(endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(dat['L1']))

    # sanity checks: beta and temperature interpolation should be identical for all T
    if len(set(betas)) != 1: raise ValueError("inconsistent beta values encountered!")
    beta = betas[0]
    if len(set(t_interp)) > 1: raise ValueError("inconsistent temperature interpolations encountered!")
    if t_interp:
        t_interp = t_interp[0]
    else:
        t_interp = xDataEnumsModule.Interpolation.linlin  # Only one temperature.

    return line, temps, energies, beta, data, e_interp, t_interp


def readT_effective(line, MF7):
    line, data = endfFileToGNDSMisc.getTAB1(line, MF7)
    interpolation = data['interpolationInfo']
    if len(interpolation) > 1: raise ValueError("Only one interpolation region allowed for T_eff data")
    interpolation = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(interpolation[0][1])

    t_axes = axesModule.Axes(2, labelsUnits={1: ('temperature', 'K'), 0: ('t_effective', 'K')})
    function1d = XYs1dModule.XYs1d(data=data['data'], axes=t_axes, interpolation=interpolation)
    return line, incoherentInelasticModule.T_effective(function1d)


def readMF7(info, MT, MF7):
    ZA, AWR, LTHR, LAT, LASYM, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF7[0], range(2, 6))
    info.ZA_massLineInfo.add(ZA, AWR, MT, 7, 0)
    if LTHR not in range(0, 4):
        raise NotImplementedError("MF7 with LTHR=%d" % LTHR)

    line = 1
    forms = []
    if LTHR in (1, 3):  # coherent elastic scattering

        line, temps, energies, dummy, data, e_interp, t_interp = readSTable(line, MF7)
        temps, energies = map(valuesModule.Values, (temps, energies))

        if e_interp != xDataEnumsModule.Interpolation.flat:
            raise ValueError("Unexpected energy interpolation encountered")

        axes = axesModule.Axes(3)
        axes[2] = axesModule.Grid('temperature', 2, 'K', xDataEnumsModule.GridStyle.points, values=temps,
                                  interpolation=t_interp)
        axes[1] = axesModule.Grid('energy_in', 1, 'eV', xDataEnumsModule.GridStyle.points, values=energies,
                                  interpolation=e_interp)
        axes[0] = axesModule.Axis('S_cumulative', 0, 'eV*b')
        array = arrayModule.Full(shape=(len(temps), len(energies)), data=data)
        Stab = griddedModule.Gridded2d(axes, array)

        forms.append(coherentElasticModule.Form(info.style, coherentElasticModule.S_table(Stab)))

    if LTHR in (2, 3):  # incoherent elastic
        line, dat = endfFileToGNDSMisc.getTAB1(line, MF7)
        SB = incoherentElasticModule.BoundAtomCrossSection(float(dat['C1']), 'b')

        e_interp = dat['interpolationInfo']
        if len(e_interp) > 1: raise ValueError("only one interpolation region allowed for Debye/Waller data")
        e_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(e_interp[0][1])

        t_axes = axesModule.Axes(2, labelsUnits={1: ('temperature', 'K'), 0: ('DebyeWallerIntegral', '1/eV')})
        if len(dat['data']) == 2 and dat['data'][0] == dat['data'][1]: dat['data'] = dat['data'][0:1]

        function1d = XYs1dModule.XYs1d(data=dat['data'], axes=t_axes, interpolation=e_interp)
        DbW = incoherentElasticModule.DebyeWallerIntegral(function1d)

        forms.append(incoherentElasticModule.Form(info.style, SB, DbW))

    if LTHR == 0:  # incoherent inelastic
        line, listy = endfFileToGNDSMisc.getList(line, MF7)
        LLN, NI, NS, b_n = listy['L1'], listy['NPL'], listy['N2'], listy['data']
        if LLN != 0:
            print("LLN=%g" % LLN)

        neutronMassAMU = info.PoPs[IDsPoPsModule.neutron].mass.float('amu')

        # b_n array contains information about principal / secondary scattering atoms
        pid = info.principalScattererPid
        massAMU = b_n[2] * neutronMassAMU
        numberPerMolecule = int(b_n[5])
        freeAtomCrossSection = b_n[0] / numberPerMolecule
        boundAtomCrossSection = incoherentInelasticModule.BoundAtomCrossSection(
            freeAtomCrossSection * ((massAMU + neutronMassAMU) / massAMU) ** 2, 'b')

        boundAtomCrossSectionByNuclide = None
        if hasattr(info, 'boundAtomCrossSections'):
            boundAtomCrossSectionByNuclide = incoherentInelasticModule.BoundAtomCrossSectionByNuclide()
            for bac in info.boundAtomCrossSections:
                boundAtomCrossSectionByNuclide.add(bac)

        principalAtom = {
            'pid': pid,
            'numberPerMolecule': numberPerMolecule,
            'mass': incoherentInelasticModule.Mass(massAMU, 'amu'),
            'e_max': incoherentInelasticModule.E_max(b_n[3], 'eV'),
            'boundAtomCrossSection': boundAtomCrossSection,
            'boundAtomCrossSectionByNuclide': boundAtomCrossSectionByNuclide,
            'selfScatteringKernel': None,  # principal scattering kernel is read after all scattering atoms
            'primaryScatterer': True,
            'e_critical': incoherentInelasticModule.E_critical(b_n[1], 'eV')
        }

        if pid in info.PoPs and len(info.PoPs[pid].mass) == 0:
            mass = massModule.Double(info.PoPsLabel, massAMU, quantityModule.stringToPhysicalUnit('amu'))
            info.PoPs[pid].mass.add(mass)

        secondaryAtoms = []
        for index in range(1, NS + 1):
            massAMU = b_n[6 * index + 2] * neutronMassAMU
            numberPerMolecule = int(b_n[6 * index + 5])
            freeAtomCrossSection = b_n[6 * index + 1] / numberPerMolecule
            boundAtomCrossSection = incoherentInelasticModule.BoundAtomCrossSection(
                freeAtomCrossSection * ((massAMU + neutronMassAMU) / massAMU) ** 2, 'b')
            if b_n[6 * index] not in (0., 1.):
                raise NotImplementedError("Unsupported kernel '%d' found for scattering atom %d" %
                                          (b_n[6 * index], index + 1))
            functionalForm = {
                0.0: incoherentInelasticModule.SCTApproximation,
                1.0: incoherentInelasticModule.FreeGasApproximation,
            }[b_n[6 * index]]()
            kernel = incoherentInelasticModule.SelfScatteringKernel(functionalForm)

            # try to figure out identity of secondary scatterer using #/molecule and mass
            scatteringElement = [c[0] for c in info.scatterers if c[1] == numberPerMolecule
                                 and not info.principalScattererPid.startswith(c[0])]
            if len(scatteringElement) == 1:
                pid = scatteringElement[0] + str(int(round(massAMU, 0)))
            else:
                print(
                    "    WARNING: could not determine identity of secondary scatterer! Will be identified as 'unknown'")
                pid = "unknown"

            secondaryAtoms.append(incoherentInelasticModule.ScatteringAtom(
                pid=pid,
                numberPerMolecule=numberPerMolecule,
                mass=incoherentInelasticModule.Mass(massAMU, 'amu'),
                e_max=incoherentInelasticModule.E_max(b_n[6 * index + 3], 'eV'),
                boundAtomCrossSection=boundAtomCrossSection,
                selfScatteringKernel=kernel)
            )

            # also add secondary scattering atom to PoPs
            if pid != "unknown":
                ZA, MAT = endf_endlModule.ZAAndMATFromParticleName(pid)
                particle = toGNDSMiscModule.getTypeNameGamma(info, ZA)
                mass = massModule.Double(info.PoPsLabel, massAMU, quantityModule.stringToPhysicalUnit('amu'))
                particle.mass.add(mass)

        line, t2header = endfFileToGNDSMisc.getTAB2Header(line, MF7)
        b_interp = t2header['interpolationInfo']
        if len(b_interp) > 1: raise ValueError("only one interpolation region allowed for S table")
        b_interp = endfFileToGNDSMisc.ENDFInterpolationToGNDS1d(b_interp[0][1])
        n_betas = int(t2header['NZ'])

        line, temps, alphas, beta, data, a_interp, t_interp = readSTable(line, MF7)  # Read first beta.
        betas = [beta]
        a_interps = [a_interp]
        t_interps = [t_interp]

        warnings = []
        for index in range(n_betas - 1):  # Read in remaining betas.
            nextLine, temps2, alphas2, beta, data2, a_interp, t_interp = readSTable(line, MF7)
            a_interps.append(a_interp)
            t_interps.append(t_interp)
            if temps != temps2:
                warnings.append(
                    "inconsistent temperatures for beta index %d starting on line %d of MF7 MT4" % (index, line))
            if alphas != alphas2:
                warnings.append(
                    "inconsistent alpha grids for beta index %d starting on line %d of MF7 MT4" % (index, line))
            betas.append(beta)
            data.extend(data2)
            line = nextLine

        if len(set(a_interps)) != 1:
            warnings.append("Inconsistent alpha interpolations in S_alpha_beta: %s" % set(a_interps))
        if len(set(t_interps)) != 1:
            warnings.append("Inconsistent temperature interpolations in S_alpha_beta: %s" % set(t_interps))
        if warnings:
            warnings = ["    WARNING: " + w for w in warnings]
            print("\n".join(warnings))
            info.doRaise.extend(warnings)

        temps, betas, alphas = map(valuesModule.Values, (temps, betas, alphas))
        for grid, array in (('temps', temps), ('betas', betas), ('alphas', alphas)):
            if list(array) != sorted(array):
                warnings.append(f"{grid} grid out of order in S_alpha_beta array!")
                info.doRaise.append(warnings[-1])

        axes = axesModule.Axes(4)
        axes[3] = axesModule.Grid('temperature', 3, 'K', xDataEnumsModule.GridStyle.points, values=temps,
                                  interpolation=t_interp)
        axes[2] = axesModule.Grid('beta', 2, '', xDataEnumsModule.GridStyle.points, values=betas,
                                  interpolation=b_interp)
        axes[1] = axesModule.Grid('alpha', 1, '', xDataEnumsModule.GridStyle.points, values=alphas,
                                  interpolation=a_interp)
        axes[0] = axesModule.Axis('S_alpha_beta', 0, '')

        arrayShape_orig = (len(betas), len(temps), len(alphas))
        data = numpy.array(data)
        if LLN == 1:  # ENDF stored ln(S) instead of S
            LLN_min = min(data)  # substitute for -inf when S=0
            data = numpy.exp(data)

            # FIXME assumes the original interpolation applied to ln(S) rather than to S. Assumption needs checking!
            for idx, interp, label in ((3, t_interp, 'T'), (2, b_interp, 'beta'), (1, a_interp, 'alpha')):
                if interp is None:
                    continue
                if not str(interp).startswith('lin-'):
                    info.doRaise.append(
                        "Expected linear interpolation for ln(S) along %s, got %s instead" % (label, interp))

                # modify interpolation since we're storing S instead of ln(S)
                axes[idx].interpolation = {
                    xDataEnumsModule.Interpolation.linlin: xDataEnumsModule.Interpolation.loglin,
                    xDataEnumsModule.Interpolation.linlog: xDataEnumsModule.Interpolation.loglog
                }.get(interp, interp)

        zeroFraction = sum(data == 0) / len(data)
        data = data.reshape(arrayShape_orig)
        # ENDF data are stored as 'beta,T,alpha'. Switch order to 'T,beta,alpha':
        data = numpy.transpose(data, axes=(1, 0, 2))
        arrayShape_new = (len(temps), len(betas), len(alphas))
        if zeroFraction > 0.1:
            array = arrayModule.Flattened.fromNumpyArray(data)
        else:
            array = arrayModule.Full(shape=arrayShape_new, data=data.flatten())
        Stab = griddedModule.Gridded3d(axes, array)

        principalAtom['selfScatteringKernel'] = incoherentInelasticModule.SelfScatteringKernel(Stab)
        principalAtom = incoherentInelasticModule.ScatteringAtom(**principalAtom)

        atoms = [principalAtom] + secondaryAtoms

        # last read T_eff tables. There must be at least one (for the principal scattering atom), plus more for
        # any other atoms that specify the short-collision-time approximation:
        line, t_eff = readT_effective(line, MF7)
        atoms[0].T_effective = t_eff
        for atom in atoms[1:]:
            if isinstance(atom.selfScatteringKernel.kernel, incoherentInelasticModule.SCTApproximation):
                line, t_eff = readT_effective(line, MF7)
                atom.T_effective = t_eff

        form = incoherentInelasticModule.Form(
            label=info.style,
            primaryScatterer=info.principalScattererPid,
            calculatedAtThermal=LAT != 0)
        for atom in atoms:
            form.scatteringAtoms.add(atom)

        if LLN == 1:
            info.ENDFconversionFlags.add(form, 'LLN=1,min=%f' % LLN_min)

        forms.append(form)

    if line != len(MF7): raise ValueError("Trailing data left in MT%i MF7!" % MT)

    return forms


def readMF7_info(info, dat):
    ZA, AWR, NA, dum, dum, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(dat[0], range(2, 6))
    neutronMassAMU = info.PoPs[IDsPoPsModule.neutron].mass.float('amu')

    targetInfo_ = targetInfo.TargetInfo()
    boundAtomCrossSections = []
    line = 1
    for element_idx in range(NA):
        dum, dum, NAS, dum, tmp, NI = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(dat[line], range(2, 6))
        assert tmp == NI * 6, f"Malformed line #{line} in MF=7 MT=451"
        line += 1

        element = None
        for nuclide_index in range(NI):
            ZAI, LISI, AFI, AWRI, SFI, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(dat[line])
            if LISI != 0:
                raise NotImplementedError("Isomeric state in MF=7 MT=451")
            pops_nuclide = toGNDSMiscModule.getTypeNameGamma(info, int(ZAI))
            info.PoPs.add(pops_nuclide)
            massAMU = AWRI * neutronMassAMU
            mass = massModule.Double(info.PoPsLabel, massAMU, quantityModule.stringToPhysicalUnit('amu'))
            info.PoPs[pops_nuclide.id].mass.add(mass)

            boundAtomCrossSections.append(incoherentInelasticModule.BoundAtomCrossSection(
                SFI * ((massAMU + neutronMassAMU) / massAMU) ** 2, 'b', label=pops_nuclide.id))
            nuclide = targetInfo.Nuclide(pops_nuclide.id, AFI)

            if element is None:
                symbol = pops_nuclide.findClassInAncestry(chemicalElementModule.ChemicalElement).symbol
                element = targetInfo.ChemicalElement(symbol)
            element.nuclides.add(nuclide)
            line += 1

        targetInfo_.isotopicAbundances.chemicalElements.add(element)

    return targetInfo_, boundAtomCrossSections
