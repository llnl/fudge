# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from xData import regions as regionsModule

from fudge.resonances import unresolved as unresolvedModule

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from brownies.legacy.toENDF6 import gndsToENDF6 as gndsToENDF6Module


def toENDF6(self, flags, targetInfo, verbosityIndent=''):
    """
    Write unresolved.TabulatedWidths back to ENDF6
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

    NLS = len(self.Ls)
    LFW = False
    reactionLabels = {'competitive': 'competitive', 'fission': 'fission'}
    for reaction in self.resonanceReactions:
        MT = reaction.link.link.ENDF_MT
        if MT == 2:
            reactionLabels['elastic'] = reaction.label
        elif MT == 102:
            reactionLabels['capture'] = reaction.label
        elif MT == 18:
            reactionLabels['fission'] = reaction.label
            LFW = True
        else:
            raise NotImplementedError("Unsupported reaction '%s' in unresolved resonanceReactions!" % reaction.label)

    LRF = targetInfo['unresolved_LRF']

    target = self.PoPs[targetInfo['targetID']]
    if hasattr(target, 'nucleus'):
        SPI = target.nucleus.spin[0].value
    else:
        SPI = target.spin[0].value

    target = self.PoPs[targetInfo['targetID']]
    AWRI = target.getMass('amu') / targetInfo['massTracker'].neutronMass

    if LFW == 0 and LRF == 1:  # 'Case A' from ENDF 2010 manual pg 70
        endf.append(endfFormatsModule.endfHeadLine(SPI, AP, targetInfo['LSSF'], 0, NLS, 0))
        for Lval in self.Ls:
            NJS = len(Lval.Js)
            endf.append(endfFormatsModule.endfHeadLine(AWRI, 0, Lval.L, 0, 6 * NJS, NJS))
            for Jval in Lval.Js:
                # here we only have one width per J:
                D = Jval.levelSpacing.data.evaluate(Jval.levelSpacing.data.domainMin)
                widths = {}
                DOFs = {}
                for label in reactionLabels:
                    for width in Jval.widths:
                        if reactionLabels[label] == width.resonanceReaction:
                            widths[label] = width.data.evaluate(width.data.domainMin)
                            DOFs[label] = width.degreesOfFreedom
                endf.append(endfFormatsModule.endfDataLine(
                    [D, Jval.J, DOFs['elastic'], widths['elastic'], widths['capture'], 0]))

    elif LFW == 1 and LRF == 1:  # 'Case B', all widths are constant except fission
        for width in self.Ls[0].Js[0].widths:
            if width.resonanceReaction == reactionLabels['fission']:
                energies = width.data.domainGrid
                break
        NE = len(energies)
        endf.append(endfFormatsModule.endfHeadLine(SPI, AP, targetInfo['LSSF'], 0, NE, NLS))
        nlines = int(math.ceil(NE / 6.0))
        for line in range(nlines):
            endf.append(endfFormatsModule.endfDataLine(energies[line * 6:line * 6 + 6]))
        for Lval in self.Ls:
            NJS = len(Lval.Js)
            endf.append(endfFormatsModule.endfHeadLine(AWRI, 0, Lval.L, 0, NJS, 0))
            for Jval in Lval.Js:
                D = Jval.levelSpacing.data.evaluate(Jval.levelSpacing.data.domainMin)
                elastic, = [width for width in Jval.widths if width.resonanceReaction == reactionLabels['elastic']]
                fission, = [width for width in Jval.widths if width.resonanceReaction == reactionLabels['fission']]
                capture, = [width for width in Jval.widths if width.resonanceReaction == reactionLabels['capture']]
                MUF = fission.degreesOfFreedom
                AMUN = elastic.degreesOfFreedom
                GNO = elastic.data.evaluate(elastic.data.domainMin)
                GG = capture.data.evaluate(capture.data.domainMin)

                endf.append(endfFormatsModule.endfHeadLine(0, 0, Lval.L, MUF, NE + 6, 0))
                endf.append(endfFormatsModule.endfDataLine([D, Jval.J, AMUN, GNO, GG, 0]))
                fissWidths = [xy[1] for xy in fission.data]
                for line in range(nlines):
                    endf.append(endfFormatsModule.endfDataLine(fissWidths[line * 6:line * 6 + 6]))

    elif LRF == 2:  # 'Case C', most common in ENDF
        endf.append(endfFormatsModule.endfHeadLine(SPI, AP, targetInfo['LSSF'], 0, NLS, 0))
        for Lval in self.Ls:
            NJS = len(Lval.Js)
            endf.append(endfFormatsModule.endfHeadLine(AWRI, 0, Lval.L, 0, NJS, 0))
            for Jval in Lval.Js:
                INT = gndsToENDF6Module.gndsToENDFInterpolationFlag(Jval.levelSpacing.data.interpolation)
                unionGrid = set(Jval.levelSpacing.data.domainGrid)
                for width in Jval.widths:
                    unionGrid.update(width.data.domainGrid)
                unionGrid = sorted(unionGrid)
                NE = len(unionGrid)
                endf.append(endfFormatsModule.endfHeadLine(Jval.J, 0, INT, 0, 6 * NE + 6, NE))

                levelSpacing = [Jval.levelSpacing.data.evaluate(x) for x in unionGrid]

                widths = {}
                DOF = {}
                for label in reactionLabels:
                    found = False
                    for width in Jval.widths:
                        if width.resonanceReaction == reactionLabels[label]:
                            DOF[label] = width.degreesOfFreedom
                            if isinstance(width.data, regionsModule.Regions1d):
                                widths[label] = []
                                for x in unionGrid:
                                    for region in width.data:
                                        if x <= region.domainMax:
                                            widths[label].append(region.evaluate(x))
                                            break
                            else:
                                widths[label] = [width.data.evaluate(x) for x in unionGrid]
                            found = True
                            break
                    if not found:
                        DOF[label] = 0
                        widths[label] = [0] * len(unionGrid)

                endf.append(endfFormatsModule.endfDataLine([0, 0,
                                                            DOF['competitive'], DOF['elastic'], DOF['capture'],
                                                            DOF['fission']]))
                table = [unionGrid, levelSpacing, widths['competitive'], widths['elastic'],
                         widths['capture'], widths['fission']]
                for row in zip(*table):
                    endf.append(endfFormatsModule.endfDataLine(row))
    return endf

unresolvedModule.TabulatedWidths.toENDF6 = toENDF6
