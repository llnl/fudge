# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from fudge.resonances import resolved as resolvedModule

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from brownies.legacy.toENDF6 import gndsToENDF6 as gndsToENDF6Module


def toENDF6(self, flags, targetInfo, verbosityIndent=''):
    """
    Write BreitWigner sections back to ENDF6
    :param self:
    :param flags:
    :param targetInfo:
    :param verbosityIndent:
    :return:
    """

    endf = []
    AP = self.hardSphereRadius or self.getScatteringRadius()
    if AP.isEnergyDependent():
        scatRadius = AP.evaluated
        NR, NP = 1, len(scatRadius)
        endf.append(endfFormatsModule.endfHeadLine(0, 0, 0, 0, NR, NP))
        endf += endfFormatsModule.endfInterpolationList(
            (NP, gndsToENDF6Module.gndsToENDFInterpolationFlag(scatRadius.interpolation)))
        endf += endfFormatsModule.endfNdDataList(scatRadius.convertAxisToUnit(0, '10*fm'))
        AP = 0
    else:
        AP = AP.getValueAs('10*fm')
    L_list = self.resonanceParameters.table.getColumn('L')
    NLS = len(set(L_list))
    LAD = getattr(self, 'computeAngularDistribution') or 0
    endf += [endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, LAD, 0, NLS, 0)]

    table = [self.resonanceParameters.table.getColumn('energy', unit='eV'),
             self.resonanceParameters.table.getColumn('J')]
    NE = len(table[0])
    for attr in ('totalWidth', 'neutronWidth', 'captureWidth', 'fissionWidth'):
        column = self.resonanceParameters.table.getColumn(attr, unit='eV')
        if not column: column = [0] * NE
        table.append(column)
    CS = self.resonanceParameters.table.getColumn('channelSpin')
    if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
        targetSpin = targetInfo['spin']
        CS = [2 * (cs - targetSpin) for cs in CS]
        Js = [v[0] * v[1] for v in zip(table[1], CS)]
        table[1] = Js
    table = list(zip(*table))

    target = self.PoPs[targetInfo['targetID']]
    AWRI = target.getMass('amu') / targetInfo['massTracker'].neutronMass

    for L in set(L_list):
        APL = 0
        resonances = [table[i] for i in range(NE) if L_list[i] == L]
        NRS = len(resonances)
        endf.append(endfFormatsModule.endfHeadLine(AWRI, APL, L, 0, 6 * NRS, NRS))
        for row in resonances:
            endf.append(endfFormatsModule.endfDataLine(row))
    return endf

resolvedModule.BreitWigner.toENDF6 = toENDF6


#
# helper functions for RMatrix
#
def getENDFtuple(spin, parity):
    # ENDF combines spin & parity UNLESS spin==0. Then it wants (0,+/-1)
    if parity is None: raise ValueError("parity is None")
    if spin:
        return (spin * parity, 0)
    else:
        return (spin, parity)


def writeRMatrixParticlePairs(RMatrix, targetInfo):
    """
    Used for writing both MF=2 and MF=32 to ENDF6
    :param RMatrix:
    :param targetInfo:
    :return:
    """
    endf = []
    NPP = len(RMatrix.resonanceReactions)
    endf.append(endfFormatsModule.endfHeadLine(0, 0, NPP, 0, 12 * NPP, 2 * NPP))
    PoPs = RMatrix.PoPs

    def MZIP(particle, ignoreMissingJpi=False):  # helper: extract mass, z, spin and parity from particle list

        mass = particle.getMass('amu') / targetInfo['massTracker'].getMassAMU(1)

        if hasattr(particle, 'nucleus'):
            particle = particle.nucleus
        Z = chemicalElementMiscPoPsModule.ZAInfo(particle)[0]
        try:
            parity = particle.parity[0].value
            spin = particle.spin.float('hbar')
            I, P = getENDFtuple(spin, parity)
        except IndexError as err:
            if ignoreMissingJpi:
                I, P = (0, 0)  # ENDF omits spin/parity for compound nucleus from capture
            else:
                raise ValueError("When formatting particle %s, encountered %s" % (particle.id, err))
        return mass, Z, I, P

    for idx, pp in enumerate(RMatrix.resonanceReactions):
        reaction = pp.link.link
        MT = reaction.ENDF_MT
        if reaction.isFission():
            pA = PoPs[targetInfo['reactionSuite'].projectile]
            pB = PoPs[targetInfo['reactionSuite'].target]
        else:
            # get the PoPs instances for ejectile and residual:
            pA = PoPs[pp.ejectile]
            pB = PoPs[pp.residual]
        MA, ZA, IA, PA = MZIP(pA)
        MB, ZB, IB, PB = MZIP(pB, ignoreMissingJpi=(MT == 102))
        PNT = 1
        if MT in (18, 19, 102):
            PNT = 0  # special case

        if RMatrix.boundaryCondition == resolvedModule.BoundaryCondition.EliminateShiftFunction or pp.eliminated:
            SHF = 0
        elif RMatrix.boundaryCondition == resolvedModule.BoundaryCondition.Brune:
            SHF = 2
        else:   # 'Given' or 'NegativeOrbitalMomentum'
            SHF = 1
        if pp.Q is not None:
            Q = pp.Q.getConstantAs('eV')
        else:
            Q = reaction.getQ('eV')
            # getQ doesn't account for residual left in excited state:
            for particle in reaction.outputChannel:
                if hasattr(particle, 'getLevelAsFloat'):
                    Q -= particle.getLevelAsFloat('eV')
        if MT in (18, 19, 102):
            Q = 0
        endf.append(endfFormatsModule.endfDataLine([MA, MB, ZA, ZB, IA, IB]))
        endf.append(endfFormatsModule.endfDataLine([Q, PNT, SHF, MT, PA, PB]))
        pp.index = idx + 1  # 1-based index in ENDF

    return endf


def writeRMatrixSpinGroupHeader(RMatrix, spingrp, targetInfo):
    endf = []
    AJ, PJ = getENDFtuple(float(spingrp.spin), int(spingrp.parity))
    KBK = len([ch for ch in spingrp.channels if ch.externalRMatrix is not None
               or "LBK=0" in targetInfo["ENDFconversionFlags"].get(ch, "")])
    KPS = 0  # AFAIK this option is not used in any library
    NCH = len(spingrp.resonanceParameters.table.columns) - 1  # skip the 'energy' column
    try:
        endf.append(endfFormatsModule.endfHeadLine(AJ, PJ, KBK, KPS, 6 * NCH, NCH))
    except TypeError as err:
        raise TypeError("Got '%s' when formatting '%s'" % (err.message, str((AJ, PJ, KBK, KPS, 6 * NCH, NCH))))
    for chan in spingrp.channels:
        rreac = RMatrix.resonanceReactions[chan.resonanceReaction]
        PPI = rreac.index
        L = chan.L
        if type(L) not in [int, float]: L = L.value
        SCH = chan.channelSpin
        if chan.boundaryConditionValue is not None:
            BND = chan.boundaryConditionValue
        elif RMatrix.boundaryCondition == resolvedModule.BoundaryCondition.EliminateShiftFunction:
            BND = 0
        elif RMatrix.boundaryCondition == resolvedModule.BoundaryCondition.NegativeOrbitalMomentum:
            BND = -L
        else:
            raise NotImplementedError("Writing boundary condition '%s' to ENDF-6" % RMatrix.boundaryCondition)

        APT = chan.getScatteringRadius()
        APE = chan.hardSphereRadius or rreac.hardSphereRadius or APT
        APT = APT.getValueAs('10*fm')
        APE = APE.getValueAs('10*fm')
        if rreac.link.link.ENDF_MT == 102:
            APT, APE = 0, 0
        endf.append(endfFormatsModule.endfDataLine([PPI, L, SCH, BND, APE, APT]))
    return NCH, endf


def toENDF6(self, flags, targetInfo, verbosityIndent=''):
    """
    Write RMatrix section to ENDF6
    :param self:
    :param flags:
    :param targetInfo:
    :param verbosityIndent:
    :return:
    """

    if targetInfo.get("LRF7_as_LRF3"):
        return writeAsLRF3(self, flags, targetInfo, verbosityIndent=verbosityIndent)

    KRM = {resolvedModule.RMatrix.Approximation.ReichMoore: 3, resolvedModule.RMatrix.Approximation.RMatrix: 4}[self.approximation]
    try:
        endf = [endfFormatsModule.endfHeadLine(0, 0, self.reducedWidthAmplitudes, KRM,
                                               len(self.spinGroups), self.relativisticKinematics)]
    except TypeError as err:
        raise TypeError("Got '%s' when formatting '%s'" %
                        (err.message, str((0, 0, self.reducedWidthAmplitudes, KRM, len(self.spinGroups),
                                           self.relativisticKinematics))))

    endf.extend(writeRMatrixParticlePairs(self, targetInfo))

    for spingrp in self.spinGroups:
        NCH, spinGroupHeader = writeRMatrixSpinGroupHeader(self, spingrp, targetInfo)
        endf.extend(spinGroupHeader)

        # resonances:
        NRS = len(spingrp.resonanceParameters.table)
        NX = (NCH // 6 + 1) * NRS
        if NRS == 0:
            NX = 1  # special case
        endf.append(endfFormatsModule.endfHeadLine(0, 0, 0, NRS, 6 * NX, NX))
        for res in spingrp.resonanceParameters.table:
            for jidx in range(NCH // 6 + 1):
                endfLine = res[jidx * 6:jidx * 6 + 6]
                while len(endfLine) < 6:
                    endfLine.append(0)
                endf.append(endfFormatsModule.endfDataLine(endfLine))
        if NRS == 0:
            endf.append(endfFormatsModule.endfDataLine([0, 0, 0, 0, 0, 0]))

        for idx, channel in enumerate(spingrp.channels):
            if channel.externalRMatrix:
                from fudge.resonances import externalRMatrix as externalRMatrixModule
                terms = channel.externalRMatrix.terms

                def getTerm(key, unit):
                    if key not in terms:
                        return 0
                    return terms[key].float(unit)

                ED = getTerm('singularityEnergyBelow', 'eV')
                EU = getTerm('singularityEnergyAbove', 'eV')
                if isinstance(channel.externalRMatrix, externalRMatrixModule.SAMMY):
                    endf.append(endfFormatsModule.endfContLine(0, 0, idx + 1, 2, 0, 0))
                    endf.append(endfFormatsModule.endfContLine(ED, EU, 0, 0, 5, 0))
                    R0 = getTerm('constantExternalR', '')
                    R1 = getTerm('linearExternalR', '1/eV')
                    R2 = getTerm('quadraticExternalR', '1/eV**2')
                    S0 = getTerm('constantLogarithmicCoefficient', '')
                    S1 = getTerm('linearLogarithmicCoefficient', '1/eV')
                    endf.append(endfFormatsModule.endfDataLine([R0, R1, R2, S0, S1, 0]))
                elif isinstance(channel.externalRMatrix, externalRMatrixModule.Froehner):
                    endf.append(endfFormatsModule.endfContLine(0, 0, idx + 1, 3, 0, 0))
                    endf.append(endfFormatsModule.endfContLine(ED, EU, 0, 0, 3, 0))
                    R0 = getTerm('constantExternalR', '')
                    S0 = getTerm('poleStrength', '')
                    GA = getTerm('averageRadiationWidth', 'eV')
                    endf.append(endfFormatsModule.endfDataLine([R0, S0, GA, 0, 0, 0]))
            elif "LBK=0" in targetInfo["ENDFconversionFlags"].get(channel, ""):
                endf.append(endfFormatsModule.endfContLine(0, 0, idx + 1, 0, 0, 0))
    return endf

resolvedModule.RMatrix.toENDF6 = toENDF6


#
# LRF=3 is converted and stored in RMatrix. Logic below translates back to LRF=3:
#
def writeAsLRF3(RMatrix, flags, targetInfo, verbosityIndent=''):
    endf = []
    elastic, = [reac for reac in RMatrix.resonanceReactions if
                reac.link.link == targetInfo['reactionSuite'].getReaction('elastic')]
    capture, = [reac for reac in RMatrix.resonanceReactions if
                reac.link.link == targetInfo['reactionSuite'].getReaction('capture')]
    AP = elastic.getScatteringRadius()
    if AP.isEnergyDependent():
        scatRadius = AP.evaluated
        NR, NP = 1, len(scatRadius)
        endf.append(endfFormatsModule.endfHeadLine(0, 0, 0, 0, NR, NP))
        endf += endfFormatsModule.endfInterpolationList((NP,
                                                         gndsToENDF6Module.gndsToENDFInterpolationFlag(
                                                             scatRadius.interpolation)))
        endf += endfFormatsModule.endfNdDataList(scatRadius.convertAxisToUnit(0, '10*fm'))
        AP = 0
    else:
        AP = AP.getValueAs('10*fm')

    table = {'Ls': [], 'energies': [], 'Js': [], 'elastic': [], 'capture': [], 'fission width_1': [],
             'fission width_2': []}
    APLs = {}
    for sg in RMatrix.spinGroups:
        elasticChan, = [chan for chan in sg.channels if chan.resonanceReaction == elastic.label]
        if elasticChan.hardSphereRadius:
            APL = elasticChan.hardSphereRadius.getValueAs('10*fm')
            if APL != AP: APLs[elasticChan.L] = APL
        J = float(sg.spin)
        if 'ignoreChannelSpin' not in targetInfo["ENDFconversionFlags"].get(RMatrix, ""):
            if J != 0 and elasticChan.channelSpin < targetInfo['spin']: J *= -1

        energies = sg.resonanceParameters.table.getColumn('energy', unit='eV')
        NE = len(energies)
        table['energies'] += energies

        table['Ls'] += [elasticChan.L] * NE
        table['Js'] += [J] * NE
        table['elastic'] += sg.resonanceParameters.table.getColumn(elastic.label + ' width', unit='eV')
        table['capture'] += sg.resonanceParameters.table.getColumn(capture.label + ' width', unit='eV')
        for column in ('fission width_1', 'fission width_2'):
            vals = sg.resonanceParameters.table.getColumn(column, unit='eV')
            if not vals: vals = [0] * NE
            table[column] += vals

    L_list = sorted(set(table['Ls']))
    NLS = len(L_list)
    if APLs and 0 not in APLs:
        APLs[0] = AP
    LAD = int(RMatrix.supportsAngularReconstruction)
    NLSC = RMatrix.LValuesNeededForAngularConvergence or 0
    conversionFlags = targetInfo['ENDFconversionFlags'].get(RMatrix, "")
    if 'LvaluesNeededForConvergence' in conversionFlags:
        NLSC = int(conversionFlags.split('LvaluesNeededForConvergence=')[1].split(',')[0])
    APtmp = AP
    if 'AP=0' in conversionFlags:
        APtmp = 0
    endf += [endfFormatsModule.endfHeadLine(targetInfo['spin'], APtmp, LAD, 0, NLS, NLSC)]

    sortedTable = sorted(zip(*(
        table['Ls'], table['energies'], table['Js'], table['elastic'], table['capture'],
        table['fission width_1'], table['fission width_2']
    )))

    defaultAP = 0
    if 'explicitAPL' in conversionFlags: defaultAP = AP

    target = RMatrix.PoPs[targetInfo['targetID']]
    AWRI = target.getMass('amu') / targetInfo['massTracker'].neutronMass

    for L in L_list:
        APL = APLs.get(L, defaultAP)
        resonances = [sortedTable[i][1:] for i in range(len(sortedTable)) if sortedTable[i][0] == L]
        NRS = len(resonances)
        endf.append(endfFormatsModule.endfHeadLine(AWRI, APL, L, 0, 6 * NRS, NRS))
        for row in resonances:
            endf.append(endfFormatsModule.endfDataLine(row))

    targetInfo['LRF3conversion'] = {  # save useful information for covariance re-translation
        'table': table,
        'sortedTable': sortedTable,
        'AP': AP,
        'APLs': APLs
    }
    return endf
