#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse
import pathlib
import shutil

from LUPY import argumentsForScripts as argumentsForScriptsModule

from PoPs import IDs as PoPsIDsModule
from PoPs import specialNuclearParticleID as PoPsSpecialNuclearParticleIDModule

from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from fudge import styles as stylesModule

from fudge import enums as enumsModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule

summaryDocString__FUDGE = '''For each reaction of a protare, writes available energy, each product's outgoing energy, energy balance, etc. to files.'''

description = '''
This script calculates, for each reaction, the available energy, the energy to each product and the energy balance. The results are
written into files.  Information (e.g., index, MT and label) about each reaction is written to the "index" file
in the output directory.  All files for a reaction are written into its own sub-directory under the output directory. The name 
of a reaction's sub-directory starts with its index in the index file.

In addition, within each reaction's directory, the directory CrossSectionWeighted is create and the data in the reaction's directory
are also written to the directory CrossSectionWeighted but with the data weighted by the reaction's cross section to the total
cross section.
'''

apdLabel = 'apd'
lowerUpperEpsilon = 1e-6
thinDefault = 1e-4

parser = argparse.ArgumentParser(description=description)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)

parser.add_argument('outputDir', type=pathlib.Path,                         help='Director where all output files are written.')
parser.add_argument('-o', '--output', type=pathlib.Path, default = '',      help='If present, file name to write the protare which contains the added energy data.')
parser.add_argument('--plot', action='store_true',                          help='If present, the energy balance by reaction curves are plotted.')
parser.add_argument('-w', '--crossSectionWeighted', action='store_true',    help='If present, the plot will be cross section weighted energy balance.')
parser.add_argument('-d', '--dump', action='store', default=None,           help='If present, the protare is written the specified file after average product data are calculated.')
parser.add_argument('-v', '--verbose', action='count', default=0,           help='If present, prints infomation.')
args = parser.parse_args()

sums = {}
reactionIndex = 0
energyBalanceCurves = []

inclusiveMT_products = {  4: PoPsIDsModule.neutron,         103: PoPsIDsModule.familiarProton, 104: PoPsIDsModule.familiarDeuteron, 
                        105: PoPsIDsModule.familiarTriton,  106: PoPsIDsModule.familiarHelion, 107: PoPsIDsModule.familiarAlpha}
for inclusiveMT_product in inclusiveMT_products:
    if inclusiveMT_product != 4:
        inclusiveMT_products[inclusiveMT_product] = PoPsSpecialNuclearParticleIDModule.nuclideID(inclusiveMT_products[inclusiveMT_product])

def doPrint(level, message):
    '''
    This function prints *message* if args.verbose >= *level*.
    '''

    if args.verbose >= level:
        print(message)

def addCurves1d(curve1, curve2):
    '''
    This function adds two XYs1d curves together. The to curves are lineaized and domains mutualized if needed.
    '''

    if not isinstance(curve1, XYs1dModule.XYs1d):
        curve1 = curve1.toPointwise_withLinearXYs(accuracy=1e-5, lowerEps=lowerUpperEpsilon, upperEps=lowerUpperEpsilon)
    if not isinstance(curve2, XYs1dModule.XYs1d):
        curve2 = curve2.toPointwise_withLinearXYs(accuracy=1e-5, lowerEps=lowerUpperEpsilon, upperEps=lowerUpperEpsilon)

    if len(curve1) > 0:
        if curve1.interpolation != curve2.interpolation:
            if curve1.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve1 = curve1.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, lowerUpperEpsilon, lowerUpperEpsilon)
            if curve2.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve2 = curve2.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, lowerUpperEpsilon, lowerUpperEpsilon)
    if not curve1.areDomainsMutual(curve2):
        curve1, curve2 = curve1.mutualify(lowerUpperEpsilon, lowerUpperEpsilon, True, curve2, lowerUpperEpsilon, lowerUpperEpsilon, True)

    return curve1 + curve2

def evaluateCurve(curve, domainValue):
    '''
    Returns the value of *curve* evaluated at *domainValue*. If the evaluated value is None or a raise
    is execute, the value is returned as 0.0.
    '''

    try:
        value = curve.evaluate(domainValue)
        if value is None:
            value = 0.0
    except:
        value = 0.0

    return value

def output(reactionOutputDir, fileName, curve, crossSection=None):
    '''
    This function writes XYs1d data to a file.
    '''

    outputName = reactionOutputDir / fileName
    with outputName.open('w') as fOut:
        print(curve.toString(), file=fOut, end='')

    if crossSection is not None:
        if isinstance(curve, regionsModule.Regions1d):
            curve = curve.toPointwise_withLinearXYs(lowerEps=1e-6)

        if fileName == 'availableEnergy.dat' and len(curve) == 2 and evaluateCurve(crossSection, curve.domainMin) == 0.0:
            priorEnergy = None                              # Add a point at the first non-zero cross section point.
            for energy, xSec in crossSection:
                if xSec != 0.0:
                    curve = curve.copy()
                    curve.setValue(energy, curve.evaluate(energy))
                    break
                priorEnergy = energy
            if priorEnergy is not None:                     # Add a point at the last zero cross section point.
                curve.setValue(priorEnergy, curve.evaluate(priorEnergy))
        elif fileName == 'Q.dat':
            curve = curve * crossSection.domainSlice(domainMin=curve.domainMin, domainMax=curve.domainMax)

        weighted = []
        for index, (xValue, yValue) in enumerate(curve):
            weighted.append([xValue, evaluateCurve(crossSection, xValue) * yValue])
        axes = curve.axes.copy()
        axes[0].unit = '%s * %s' % (axes[0].unit, crossSection.axes[0].unit)
        weighted = XYs1dModule.XYs1d(data=weighted, axes=axes)
        outputName = reactionOutputDir / 'CrossSectionWeighted' / fileName
        outputName.parent.mkdir(exist_ok=True)
        with outputName.open('w') as fOut:
            print(weighted.toString(), file=fOut, end='')

        if args.crossSectionWeighted:
            return weighted

        return curve

def outputProductData(reactionOutputDir, crossSection, outputChannel, productSums, level):
    '''
    This function loops over all products of *outputChannel*, outputting each products energy data and summing by product pid.
    '''

    for productIndex, product in enumerate(outputChannel):
        doPrint(2, (level + 1) * '    ' +  '    %s' % product)
        if product.outputChannel is None:
            if apdLabel in product.averageProductEnergy:
                averageProductEnergy = product.averageProductEnergy[apdLabel]
            else:
                averageProductEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
            fileName = 'averageProductEnergy_%.2d_%.2d_%s.dat' % (level, productIndex, product.pid)
            output(reactionOutputDir, fileName, averageProductEnergy, crossSection)
            if product.pid not in productSums:
                productSums[product.pid] = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
            productSums[product.pid] = addCurves1d(productSums[product.pid], averageProductEnergy)
        else:
            outputProductData(reactionOutputDir, crossSection, product.outputChannel, productSums, level + 1)

def checkTwoBody(reactionOutputDir, crossSection, outputChannel):
    """
    If the reaction is a two-body reaction, This function writes out energy balance data for the two two-body particles.
    """

    if outputChannel.genre == enumsModule.Genre.twoBody:
        product1, product2 = outputChannel.products[0], outputChannel.products[1]

        output(reactionOutputDir, 'twoBodyAverageProductEnergy_%s.dat' % product1.label, product1.averageProductEnergy[apdLabel], crossSection)
        if apdLabel not in product2.averageProductEnergy:
            print('WARNING: Two-body residual missing average product energy:' + product2.toXLink())
            return
        output(reactionOutputDir, 'twoBodyAverageProductEnergy_%s.dat' % product2.label, product2.averageProductEnergy[apdLabel], crossSection)

        averageProductEnergy = product1.averageProductEnergy[apdLabel] + product2.averageProductEnergy[apdLabel]
        output(reactionOutputDir, 'twoBodyAverageProductEnergy.dat', averageProductEnergy, crossSection)

        Q = outputChannel.Q[0].value

        domainMin = averageProductEnergy.domainMin
        domainMax = averageProductEnergy.domainMax
        availableEnergy = XYs1dModule.XYs1d(data=[[domainMin, domainMin + Q], [domainMax, domainMax + Q]], axes=averageProductEnergy.axes)
        energyBalance = addCurves1d(availableEnergy, -averageProductEnergy)
        output(reactionOutputDir, 'twoBodyEnergyBalance.dat', energyBalance, crossSection)

def process(reaction, isReaction):
    """
    This function processes energy data for one reaction.
    """

    global reactionIndex

    doPrint(1, '%5d: %s' % (reactionIndex, reaction))
    print(' %5d | %3d | %s' % (reactionIndex, reaction.ENDF_MT, reaction), file=fLog)

    MT = reaction.ENDF_MT

    reactionOutputDir = outputDir / ('%.4d_%.3d' % (reactionIndex, reaction.ENDF_MT))
    reactionOutputDir.mkdir()

    if reaction.isCoulombReaction():
        reactionIndex += 1
        print('        INFO: Skipping Coulomb scattering reaction.')
        return

    crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6)
    output(reactionOutputDir, 'crossSection.dat', crossSection)

    outputChannel = reaction.outputChannel
    if isReaction:
        Q = outputChannel.Q[0].toPointwise_withLinearXYs(accuracy=1e-5, lowerEps=1e-6)
        output(reactionOutputDir, 'Q.dat', Q, crossSection)
        availableEnergy = reaction.availableEnergy[apdLabel]
        output(reactionOutputDir, 'availableEnergy.dat', availableEnergy, crossSection)

        checkTwoBody(reactionOutputDir, crossSection, outputChannel)

    productSums = {}
    outputProductData(reactionOutputDir, crossSection, outputChannel, productSums, 0)

    totalProductEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
    for pid in productSums:
        totalProductEnergy = addCurves1d(totalProductEnergy, productSums[pid])
        output(reactionOutputDir, '_%s_averageProductEnergy.dat' % pid, productSums[pid], crossSection)
    plotCurve = output(reactionOutputDir, 'totalProductEnergy.dat', totalProductEnergy, crossSection)

    if isReaction:
        energyBalance = addCurves1d(availableEnergy, -totalProductEnergy)
        plotCurve = output(reactionOutputDir, 'energyBalance.dat', energyBalance, crossSection)

        for sumMT in sums:
            indices = sums[sumMT]['indices']
            if indices[0] <= MT <= indices[1]:
                sums[sumMT]['crossSection'] = addCurves1d(sums[sumMT]['crossSection'],  crossSection)
                sums[sumMT]['crossSection_availableEnergy'] = addCurves1d(sums[sumMT]['crossSection_availableEnergy'], crossSection * availableEnergy)
                for pid in productSums:
                    if pid not in sums[sumMT]['crossSection_energy']:
                        sums[sumMT]['crossSection_energy'][pid] = XYs1dModule.XYs1d(axes=crossSectionEnergyAxes)

                    crossSection2 = crossSection
                    productSum = productSums[pid]
                    if not crossSection2.areDomainsMutual(productSum):
                        crossSection2, productSum = crossSection2.mutualify(lowerUpperEpsilon, lowerUpperEpsilon, True, productSum, lowerUpperEpsilon, lowerUpperEpsilon, True)
                    sums[sumMT]['crossSection_energy'][pid] = addCurves1d(sums[sumMT]['crossSection_energy'][pid], crossSection2 * productSum)
    else:
        plotCurve = -totalProductEnergy
    plotCurve.plotLabel = str(reaction)
    energyBalanceCurves.append(plotCurve)

    reactionIndex += 1

if __name__=='__main__':
    doPrint(1, 'Reading protare:')
    protare = singleProtareArguments.protare(args)
    crossSectionAxes = crossSectionModule.defaultAxes(protare.domainUnit)
    averageProductEnergyAxes = averageProductEnergyModule.defaultAxes(protare.domainUnit)
    crossSectionEnergyAxes = averageProductEnergyModule.defaultAxes(protare.domainUnit)
    crossSectionEnergyAxes[0].unit = '%s * %s' % (crossSectionAxes[0].unit, averageProductEnergyAxes[0].unit)

    for MT, indices in [(4, (50, 91)), (103, (600, 649)), (104, (650, 699)), (105, (700, 749)), (106, (750, 799)), (107, (800, 849)), (102, (900, 999))]:
        sums[MT] = {'indices': indices, 'crossSection': XYs1dModule.XYs1d(axes=crossSectionAxes), 
                    'crossSection_energy': {},
                    'crossSection_availableEnergy': XYs1dModule.XYs1d(axes=crossSectionEnergyAxes)}

    rootStyle = protare.styles.preProcessingChainHead()

    AEPStyle = protare.styles.getStyleOfClass(stylesModule.AverageProductData)
    if AEPStyle is None:
        doPrint(1, 'Calculating average product data:')
        AEPStyle = stylesModule.AverageProductData(apdLabel, rootStyle.label)
        protare.styles.add(AEPStyle)
        protare.calculateAverageProductData(AEPStyle, indent='  ', skipCoulombPlusNuclearElasticError=True)
        if args.dump is not None:
            protare.saveToFile(args.dump)

    outputDir = args.outputDir / ('%s+%s' % (protare.projectile, protare.target))
    if outputDir.exists():
        shutil.rmtree(outputDir)
    outputDir.mkdir(parents=True)

    fLog = open(outputDir / 'index', 'w')
    print('# file name: %s' % protare.sourcePath, file=fLog)
    print('# absolute path: %s' % pathlib.Path(protare.sourcePath).resolve(), file=fLog)
    print('# energy unit: %s' % protare.domainUnit, file=fLog)
    print(file=fLog)

    print(' index |  MT | reaction label', file=fLog)
    print(50 * '=', file=fLog)

    doPrint(1, 'Processing:')
    for reaction in protare.reactions:
        process(reaction, True)
    for reaction in protare.orphanProducts:
        process(reaction, False)
    
    for sumMT in sorted(sums):
        indices = sums[sumMT]['indices']
        crossSection = sums[sumMT]['crossSection']
        if len(crossSection) > 0:
            sumOutputDir = outputDir / ('s%.3d_%.3d-%.3d' % (sumMT, indices[0], indices[1]))
            sumOutputDir.mkdir()

            output(sumOutputDir, 'crossSection.dat', crossSection)

            crossSection_availableEnergy = sums[sumMT]['crossSection_availableEnergy']
            if len(crossSection_availableEnergy) > 0:
                availableEnergy = crossSection_availableEnergy / crossSection
                output(sumOutputDir, 'availableEnergy.dat', availableEnergy.thin(thinDefault), crossSection)

                totalProductEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
                crossSection_energies = sums[sumMT]['crossSection_energy']
                for pid in crossSection_energies:
                    productEnergy = crossSection_energies[pid]
                    if len(productEnergy) == 0:
                        productEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
                    else:
                        try:
                            productEnergy = productEnergy / crossSection
                        except:
                            print('WARNING: for MT %s and product %s, division of product energy by cross section failed: continuing anyway.' % (sumMT, pid))
                            continue
                        totalProductEnergy = addCurves1d(totalProductEnergy, productEnergy)
                    output(sumOutputDir, '_%s_averageProductEnergy.dat' % pid, productEnergy.thin(thinDefault), crossSection)

                output(sumOutputDir, 'totalProductEnergy.dat', totalProductEnergy, crossSection)

                energyBalance = addCurves1d(availableEnergy, -totalProductEnergy)
                output(sumOutputDir, 'energyBalance.dat', energyBalance.thin(thinDefault), crossSection)

    if args.plot:
        title='Energy balance per reaction - %s' % str(protare)
        xLabel = 'Energy [%s]' % protare.domainUnit
        if args.crossSectionWeighted:
            yLabel = 'Cross section * (energy available - sum of product energy) [%s]' % energyBalanceCurves[0].axes[0].unit
            total = XYs1dModule.XYs1d()
            for curve in energyBalanceCurves:
                total = addCurves1d(total, curve)
            total.plotLabel = 'total'
            energyBalanceCurves.append(total)
        else:
            yLabel = 'Energy available - sum of product energy) [%s]' % energyBalanceCurves[0].axes[0].unit
        XYs1dModule.XYs1d.multiPlot(energyBalanceCurves, title=title, xLabel=xLabel, yLabel=yLabel)

    if str(args.output) != '.':
        protare.saveToFile(args.output)

    fLog.close()
