# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from fudge.resonances import resonances as resonancesModule, common as commonResonancesModule, \
    resolved as resolvedModule

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule


def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """
    Write gnds.resonances.resonances to ENDF6

    :param self:
    :param endfMFList:
    :param flags:
    :param targetInfo:
    :param verbosityIndent:
    :return:
    """

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0
    ZAI = ZAM  # assuming only one isotope per file

    # get target spin from the particle list:
    reactionSuite = targetInfo['reactionSuite']
    targetID = reactionSuite.target
    if targetID in reactionSuite.PoPs.aliases:
        targetID = reactionSuite.PoPs[targetID].pid
    targetInfo['targetID'] = targetID
    target = reactionSuite.PoPs[targetID]
    if hasattr(target, 'nucleus'):
        targetInfo['spin'] = target.nucleus.spin[0].value
    else:
        targetInfo['spin'] = target.spin[0].value

    endf = [endfFormatsModule.endfHeadLine(ZAM, AWT, 0, 0, NIS, 0)]
    resolvedCount, unresolvedCount = 0, 0
    # resolved may have multiple energy regions:
    if self.resolved is not None:
        resolvedCount = 1
        if isinstance(self.resolved.evaluated, commonResonancesModule.EnergyIntervals):
            resolvedCount = len(self.resolved.evaluated)
    if self.unresolved is not None: unresolvedCount = 1

    # resonances may only contain a scattering radius:
    if not (resolvedCount + unresolvedCount) and self.scatteringRadius:
        scatRadius = self.scatteringRadius.evaluated
        lowerBound, upperBound = scatRadius.domainMin, scatRadius.domainMax
        endf.append(endfFormatsModule.endfHeadLine(ZAM, ABN, 0, 0, 1, 0))
        endf.append(endfFormatsModule.endfHeadLine(lowerBound, upperBound, 0, 0, 0, 0))
        AP = PQUModule.PQU(self.scatteringRadius.evaluated.value, self.scatteringRadius.evaluated.rangeUnit).getValueAs('10*fm')
        endf.append(endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, 0, 0, 0, 0))
        endf.append(endfFormatsModule.endfSENDLineNumber())
        endfMFList[2][151] = endf
        return

    # LFW only applies to unresolved, but must be written at the start of MF2
    LFW = int(unresolvedCount != 0 and reactionSuite.hasFission())
    if unresolvedCount != 0:
        LRFurr = 2  # default = all widths energy-dependent
        NAPSurr = 0 # default = calculate channel radius
        LRF_LFW = targetInfo['ENDFconversionFlags'].get(self.unresolved.evaluated)
        if LRF_LFW is not None:
            for flag in LRF_LFW.split(','):
                if 'LRF' in flag:
                    LRFurr = int(flag.replace('LRF', ''))
                elif 'LFW' in flag:
                    LFW = int(flag.replace('LFW', ''))
                elif 'NAPS' in flag:
                    NAPSurr = int(flag.replace('NAPS=', ''))
                else:
                    raise Exception("Unrecognized ENDF conversion flag %s for unresolved region" % flag)
    NER = resolvedCount + unresolvedCount
    endf.append(endfFormatsModule.endfHeadLine(ZAI, ABN, 0, LFW, NER, 0))
    for idx in range(resolvedCount):
        if resolvedCount == 1:
            region = self.resolved
        else:
            region = self.resolved.evaluated[idx]
        LRU = 1  # resolved
        if isinstance(region.evaluated, resolvedModule.BreitWigner):
            LRF = {
                resolvedModule.BreitWigner.Approximation.SingleLevel: 1,
                resolvedModule.BreitWigner.Approximation.MultiLevel: 2
            }[region.evaluated.approximation]
        elif isinstance(region.evaluated, resolvedModule.RMatrix):
            LRF = 7
            conversionFlags = targetInfo['ENDFconversionFlags'].get(region.evaluated)
            if conversionFlags is not None and 'LRF3' in conversionFlags:
                LRF = 3
                targetInfo['LRF7_as_LRF3'] = True
        EL, EH = region.domainMin, region.domainMax
        NRO = 0
        if LRF in (1, 2):
            # NRO=1 should mean energy-dependent hardSphereRadius, but GNDS-2.1 and earlier
            # incorrectly made the scatteringRadius energy-dependent instead
            for radius in (region.evaluated.hardSphereRadius, region.evaluated.scatteringRadius):
                if radius and radius.isEnergyDependent():
                    NRO = 1
        NAPS = not region.evaluated.calculateChannelRadius
        if LRF == 7:
            NAPS = 1    # never compute radii for LRF=7
        endf.append(endfFormatsModule.endfHeadLine(EL, EH, LRU, LRF, NRO, NAPS))
        endf += region.evaluated.toENDF6(flags, targetInfo, verbosityIndent)
    if unresolvedCount != 0:
        LRU = 2  # unresolved
        region = self.unresolved
        EL, EH = region.domainMin, region.domainMax
        NRO = 0
        if region.evaluated.getHardSphereRadius().isEnergyDependent(): NRO = 1
        NAPS = NAPSurr or not region.evaluated.calculateChannelRadius
        endf.append(endfFormatsModule.endfHeadLine(EL, EH, LRU, LRFurr, NRO, NAPS))
        # pass LFW/LRF so we don't have to calculate twice:
        targetInfo['unresolved_LFW'] = LFW
        targetInfo['unresolved_LRF'] = LRFurr
        targetInfo['LSSF'] = self.unresolved.evaluated.useForSelfShieldingOnly
        targetInfo['regionEnergyBounds'] = (region.domainMin, region.domainMax)
        endf += region.evaluated.toENDF6(flags, targetInfo, verbosityIndent)
    endf.append(endfFormatsModule.endfSENDLineNumber())
    endfMFList[2][151] = endf


resonancesModule.Resonances.toENDF6 = toENDF6
