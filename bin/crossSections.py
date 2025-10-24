#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import pathlib
import shutil
import argparse

from pqu import PQU as PQUModule
from LUPY import argumentsForScripts as argumentsForScriptsModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as  PoPsIDsModule

from fudge import enums as enumsModule
from fudge import styles as stylesModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import coherent as coherentModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherentDoppler as incoherentDopplerModule

summaryDocString__FUDGE = '''Outputs the evaluated cross section for each reaction and total for a GNDS reactionSuite.  The processed data are written to files if the "-o" and "--processed" options are used.'''

description = '''
This script prints the total cross section data for the requested protare if input options "--outputDir" and "--plot" are not specified.
If the option "--outputDir" is used, cross section data for each reaction
and total are written to files in the directory specified by "--outputDir" and the protares input channel name.
The processed (i.e., heated and multi-group) cross section are also written to a file if the "--processed"
option is specified.
'''

energyUnitDefault = 'MeV'

parser = argparse.ArgumentParser(description=description)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)

parser.add_argument('--energyUnit', choices=('eV', 'MeV'), default=energyUnitDefault,   help='Energy unit to use. Default is "%s".' % energyUnitDefault)
parser.add_argument('-o', '--outputDir',                                                help='Write cross section for each reaction and total to the specified directory.')
parser.add_argument('--plot', action = 'store_true',                                    help='Draw plot of all reaction cross sections (requires PyQt5).')
parser.add_argument('--processed', action='store_true',                                 help='If prsent, cross sections for heated and multi-group styles are also written.')
parser.add_argument('-v', '--verbose', action='count', default = 0,                     help='Increase verbosity (-v, -vv, -vvv).')

args = parser.parse_args()

protare = singleProtareArguments.protare(args, verbosity=args.verbose, lazyParsing=True)
if args.energyUnit != protare.domainUnit: protare.convertUnits({protare.domainUnit: args.energyUnit})
protareLabel = '%s+%s' % (protare.projectile, protare.target)

if protare.interaction == enumsModule.Interaction.TNSL:
    if not args.processed:
        print('For TNSL need to include "--processed" option as only processed cross sections at specific temperatures are defined.')
        sys.exit()

labelsToWrite = []
if args.processed:
    for style in protare.styles:
        if isinstance(style, (stylesModule.Heated, stylesModule.HeatedMultiGroup, stylesModule.CoulombPlusNuclearElasticMuCutoff)):
            labelsToWrite.append(style.label)

outputDir = args.outputDir
if outputDir is not None:
    outputDir = pathlib.Path(outputDir) / protareLabel

class Gridded1dHistogram:
    """
    This class is used to create an instances of the multi-group cross section that adds
    the methods needed to write the multi-group cross section to a file with energy and
    cross section column data.
    """

    def __init__(self, crossSection):
        """
        :param crossSection:    A multi-group Gridded1d cross section instance.
        """

        self.gridded1d = crossSection

    def __len__(self):
        """
        This method returns the length of the multi-group cross section.

        :returns:       A python int.
        """

        return len(self.gridded1d)

    @property
    def axes(self):
        """
        This method returns the axes of the multi-group cross section.

        :returns:       An xData Axes instance.
        """

        return self.gridded1d.axes

    @property
    def label(self):
        """
        This method returns the label of the multi-group cross section.

        :returns:       A python str.
        """

        return self.gridded1d.label

    def toString(self, format):
        """
        This method returns a two column string representation of the multi-group cross section with
        the first column being the energy and the second column being the multi-group cross section.
        Points are added with a below each multi-group energy boundary, except for the first and
        last boundaries, so that lin-lin plotting gives the histogram shape of the multi-group cross section.

        :returns:       A lin-lin string representation of the histogram cross section.
        """

        lines = []
        y1 = None
        for x2, y2 in zip(*self.gridded1d.copyDataToXsAndYs()):
            if y1 is not None:
                lines.append('  %16.8e %16.8e' % (x2 * (1 - 1e-7), y1))
            lines.append('  %16.8e %16.8e' % (x2, y2))
            y1 = y2
        lines.pop(-2)

        return '\n'.join(lines)

def output(MT, reactionStr, prefix, crossSection, yLabel, subDir=None):
    """
    This function writes the *crossSection* to a file.
    """

    outputDir2 = outputDir
    if subDir is not None:
        outputDir2 = outputDir2 / subDir
        outputDir2.mkdir(parents=True, exist_ok=True)
    if len(crossSection) == 0: return

    if crossSection.label in labelsToWrite:
        subProcessDir = 'heated'
        if isinstance(crossSection, Gridded1dHistogram):
            subProcessDir = 'MultiGroup'
        outputDir2 = outputDir2 / subProcessDir/ crossSection.label

    fileName = prefix + '_%.3d.dat' % MT
    path = outputDir2 / fileName
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as fOut:
        fOut.write('# reaction = "%s"\n' % reactionStr)
        fOut.write('# x axes label = "Projectile energy [%s]"\n' % crossSection.axes[1].unit)
        fOut.write('# y axes label = "%s [%s]"\n' % (yLabel, crossSection.axes[0].unit))
        fOut.write(crossSection.toString(format='  %16.8e %16.8e'))

def outputLabel(crossSection, MT, reactionStr, reactionIndex):
    """
    This function writes the processed *crossSection* to a file.
    """

    if isinstance(crossSection, crossSectionModule.Gridded1d):
        crossSection = Gridded1dHistogram(crossSection)
    output(MT, reactionStr, reactionIndex, crossSection, 'Cross section')

if outputDir is not None:
    if outputDir.exists():
        shutil.rmtree(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    
    outputLog = open(outputDir / 'index', 'w')
    outputLog.write('index   MT  label                          \n')
    outputLog.write('-------------------------------------------\n')

crossSections = []
if args.verbose > 0: print(protare.sourcePath)
total = crossSectionModule.XYs1d(axes=crossSectionModule.defaultAxes(energyUnit=protare.domainUnit))
for reactionCounter, reaction in enumerate(protare.reactions):
    if args.verbose > 1: print('    %s' % reaction)
    if protare.interaction != enumsModule.Interaction.TNSL:
        try:
            crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
        except:
            print('WARNING: cross section lookup failed for reaction "%s".' % reaction)
        else:
            if not crossSection.areDomainsMutual(total):
                if args.verbose > 2: print('        Mutualifing domains.')
                total, crossSection = total.mutualify(1e-6, -1e-6, True, crossSection, 1e-6, -1e-6, True)
            crossSection.plotLabel = str(reaction)
            crossSections.append(crossSection)
            total = total + crossSection
            if outputDir is not None:
                outputLog.write('%5d  %3d  %s\n' % (reactionCounter, reaction.ENDF_MT, str(reaction)))
                output(reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter, crossSection, 'Cross section')

    if outputDir is not None:
        for labelToWrite in labelsToWrite:
            if labelToWrite in reaction.crossSection:
                outputLabel(reaction.crossSection[labelToWrite], reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter)

reactionCounterOffset = reactionCounter + 1

if len(protare.sums.crossSectionSums) > 0:
    if outputDir is not None:
        outputLog.write('\nData in crossSectionSums:\n')

for reactionCounter2, reaction in enumerate(protare.sums.crossSectionSums):
    reactionCounter = reactionCounterOffset + reactionCounter2
    if args.verbose > 1: print('    %s' % reaction)
    crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
    crossSection.plotLabel = str(reaction)
    crossSections.append(crossSection)
    if outputDir is not None:
        outputLog.write('%5d  %3d  %s\n' % (reactionCounter, reaction.ENDF_MT, str(reaction)))
        output(reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter, crossSection, 'Cross section', subDir='crossSectionSums')

if outputDir is not None:
    if len(protare.productions) > 0:
        outputLog.write('\nProduction reactions:\n')
    for reactionCounter, reaction in enumerate(protare.productions):
        if args.verbose > 1: print('    Production: %s' % reaction)
        crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
        crossSection.plotLabel = str(reaction)
        outputLog.write('%5d  %3d  %s\n' % (reactionCounter, reaction.ENDF_MT, str(reaction)))
        output(reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter, crossSection, 'Production cross section', subDir='Productions')

    if protare.interaction == enumsModule.Interaction.atomic and protare.projectile == PoPsIDsModule.photon:
        for reaction in protare.reactions:
            if reaction.ENDF_MT == 502:
                prefix = 'coherentFormFactor'
                formClass = coherentModule.Form
                dataName = 'formFactor'
            elif reaction.ENDF_MT == 504:
                prefix = 'incoherentScatteringFactor'
                formClass = incoherentModule.Form
                dataName = 'scatteringFactor'
            elif reaction.ENDF_MT >= 1534 and reaction.ENDF_MT <= 1572:
                prefix = 'incoherentDoppler'
                formClass = incoherentDopplerModule.Form
                dataName = 'ComptonProfile'
            else:
                continue
            for form in reaction.doubleDifferentialCrossSection:
                if isinstance(form, formClass):
                    data = getattr(form, dataName)
                    if hasattr(data, 'data'):
                        data = data.data
                    invLengthToEnergy = (PQUModule.PQU(1, data.axes[1].unit) * PQUModule.PQU('1 hplanck * c')).getValueAs(args.energyUnit)
                    data = data.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
                    data.setData(data.nf_pointwiseXY.scaleOffsetXAndY(xScale=invLengthToEnergy))
                    data.axes[1].unit = args.energyUnit
                    output(reaction.ENDF_MT, str(reaction), prefix, data, prefix)
                    break

total.plotLabel = 'total'
crossSections.insert(0, total)

if outputDir is None and not args.plot:
    print(total.toString(format='  %16.8e %16.8e'))
elif outputDir is not None:
    output(1, 'total', 'total', total, 'Cross section')
if args.plot: crossSections[0].multiPlot(crossSections, title=protareLabel, xLabel='Energy [%s]' % total.axes.axes[1].unit,
        yLabel='Cross section [%s]' % total.axes[0].unit)
