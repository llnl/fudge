#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import shutil
import argparse
import pathlib

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as PoPsID_Module
from PoPs.decays import misc as PoPsDecayMiscModule

from LUPY import times as timesModule
from LUPY import argumentsForScripts as argumentsForScriptsModule
from fudge import styles as stylesModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import energy as energyModule

energyUnitDefault = 'MeV'
temperatureUnitDefault = 'MeV/k'
discreteGammaResolutionDefault = 1e-2
twoBodyCOMResolutionDefault = 1e-2
twoBodyCOM_energyResolutionDefault = 0.01
reactionCounter = -1
outputLog = None
accuracy = 1e-3
lowerUpperEpsilon = -1e-7

summaryDocString__FUDGE = '''For the specified projectile energy and product, outputs energy spectra by reaction and also summed spectra.'''

description = '''
    Outputs the energy spectrum for the specified outgoing particle at the specified incident energy from 
    a GNDS protare file.  If option "--outputDir" is present, spectra for total and each reaction are written to the
    directory specified by the next command parameter. Otherwise, output for the total spectrum is printed.

    With the "--outputDir" option, three curves for each spectrum are written to the outputDir. They are the 
    spectrum (cross section times multiplicity times pdf) with suffix ".spec", the pdf with suffix "._pdf" 
    and the cdf with suffix "._cdf". For each reaction the file name prefix is the index of the reaction 
    in the protare and the ENDF MT value as "III_MMM" where "III" is the reaction's index and "MMM" is the
    reaction's ENDF MT number.  For example, "004_054._pdf" may specify data for the 4rd (n,n') level. 
    Information about each reaction is given in an "index" file in the output directory.

    If discrete branching gammas are present in the PoPs node, their delta functions are treated as a
    triangle with energy width twice the photons energy times the parameter after the
    "--discreteGammaResolution" option. For example, if "--discreteGammaResolution" is 0.01 and a discrete 
    photon has energy 3.2 MeV then the base of the triangle will be 2 * 0.01 * 3.2 MeV = 0.064 MeV wide.

    By default, the outputted energy spectra are in the lab frame. However, if the "--com" option is present, 
    only outgoing particles with a distribution in the center-of-mass frame are processed and the outputted 
    energy spectra are in the center-of-mass frame.

    For thermal neutron scattering law (TNSL) data, a temperature must be specified if there is more than one 
    temperature in the TSNL protare.
'''

parser = argparse.ArgumentParser( description = description )

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments( parser )

parser.add_argument( 'product', type = str,                                                     help = 'Product PoPs id (e.g., n, photon).' )
parser.add_argument( 'energy', type = float,                                                    help = 'Incident energy in units of "--energyUnit".' )
parser.add_argument( '--energyUnit', choices = ( 'eV', 'MeV' ), default = energyUnitDefault,    help = 'Energy unit to use. Default is "%s".' % energyUnitDefault )
parser.add_argument('--muMin', action='store', default=-1, type=float,                          help = 'The outgoing distribution is integrated from muMin to muMax. The default is -1.')
parser.add_argument('--muMax', action='store', default=1, type=float,                           help = 'The outgoing distribution is integrated from muMin to muMax. The default is 1.')
parser.add_argument( '--com', action = 'store_true',                                            help = 'Produce the energy spectra in the center-of-mass frame. The default is the lab frame.' )
parser.add_argument( '--temperature', type = float,                                             help = 'Specifies the temperature of the target material.' )
parser.add_argument( '--temperatureUnit', type = str, default = temperatureUnitDefault,         help = 'Temperature unit to convert to. Default is "%s".' % temperatureUnitDefault )
parser.add_argument( '--outputDir', action = 'store', type=pathlib.Path,                        help = 'If present, output energy spectrum for each reaction and total are written to the directory specified by this option.' )
parser.add_argument( '--twoBodyCOMResolution', action = 'store', type = float, default = twoBodyCOMResolutionDefault,
                                                                                                help = 'This options is deprecated, please use --twoBodyCOMEnergyResolution. When option "--com" is present, two-body reaction product a delta function for the outgoing particles. Their pdfs are treated as a triangle with energy width equal to twoBodyCOMResolution * its energy. Default value is %s' % twoBodyCOMResolutionDefault )
parser.add_argument( '--twoBodyCOM_energyResolution', action = 'store', type = float, default = None,
                                                                                                help = 'When option "--com" is present, two-body reaction product a delta function for the outgoing particles. Their pdfs are treated as a triangle with energy width equal to twoBodyCOM_energyResolution. Default value is %s MeV.' % twoBodyCOM_energyResolutionDefault )
parser.add_argument( '--discreteGammaResolution', action = 'store', type = float, default = discreteGammaResolutionDefault,
                                                                                                help = 'All discrete gammas pdfs are treated as a triangle with energy width equal to discreteGammaResolution * its energy. Default value is %s' % discreteGammaResolutionDefault )
parser.add_argument( '--plot', action = 'store_true',                                           help = 'If present, the total spectrum is plotted.' )
parser.add_argument( '--skipDelayedNeutrons', action = 'store_true',                            help = 'If present, delayed neutrons are not included in the energy spectra.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,                          help = 'Determines the amount of vervosity - the more the merrier.' )
parser.add_argument( '--selfCheck', action = 'store_true',                                      help = 'If true, no spectrum is printed to the screen.' )

args = parser.parse_args( )
outputDir = args.outputDir

frame = xDataEnumsModule.Frame.lab
if args.com:
    frame = xDataEnumsModule.Frame.centerOfMass

muMin = min(1.0, max(-1, args.muMin))
muMax = min(1.0, max(-1, args.muMax))
if muMin >= muMax:
    raise Exception('muMin = %s must be less than muMax = %s.' % (muMin, muMax))

minimumLogPoint = PQUModule.PQU(1.e-11, 'MeV').getValueAs(args.energyUnit)
if args.twoBodyCOM_energyResolution is None:
    args.twoBodyCOM_energyResolution = PQUModule.PQU(twoBodyCOM_energyResolutionDefault, 'MeV').getValueAs(args.energyUnit)

productID = { PoPsID_Module.neutron : PoPsID_Module.neutron,
            'H1' : 'H1',            'p' : 'H1',
            'H2' : 'H2',            'd' : 'H2',
            'H3' : 'H3',            't' : 'H3',
            'He3' : 'He3',          'h' : 'He3',
            'He4' : 'He4',          'a' : 'He4',
            PoPsID_Module.photon : PoPsID_Module.photon,    'gamma' : PoPsID_Module.photon,     'g' : PoPsID_Module.photon }[args.product]

readTime = timesModule.Times( )
protare = singleProtareArguments.protare(args, verbosity = args.verbose, lazyParsing = True)
readTime = readTime.toString( current = False )

if outputDir is not None:
    outputDir = args.outputDir / ('%s+%s' % (protare.projectile, protare.target))
    if outputDir.exists():
        shutil.rmtree(outputDir)
    outputDir.mkdir(parents=True)

unitConversionTime = None
if( args.energyUnit != protare.domainUnit ) :
    unitConversionTime = timesModule.Times( )
    protare.convertUnits( { protare.domainUnit : args.energyUnit } )
    unitConversionTime = unitConversionTime.toString( current = False )

styleLabel = protare.styles.preProcessingChainHead( ).label
crossSectionUnit = 'b'
for reaction in protare.reactions:
    if hasattr(reaction.crossSection[styleLabel], 'axes'):
        crossSectionUnit = reaction.crossSection[styleLabel].axes[0].unit
        break

energyAxes = energyModule.defaultAxes( args.energyUnit )
spectrumAxes = energyAxes.copy( )
spectrumAxes.axes[0].unit = '%s/%s' % ( crossSectionUnit, args.energyUnit )

sums = {}
for MT, indices in [(4, (50, 91)), (103, (600, 649)), (104, (650, 699)), (105, (700, 749)), (106, (750, 799)), (107, (800, 849)), (102, (900, 999))]:
    sums[indices] = [ 0.0, XYs1dModule.XYs1d( axes = spectrumAxes ), MT ]

reactionTiming = []

temperatures = {}
for style in protare.styles :
    if( isinstance( style, stylesModule.Heated ) ) :
        temperatures[style.temperature.getValueAs( args.temperatureUnit )] = style.label 

if( protare.isThermalNeutronScatteringLaw( ) ) :
    if( len( temperatures ) == 0 ) : raise Exception( 'For a TNSL protare, procssing to a temperature is required; this protare has not been processed.' )
    if( args.temperature is None ) :
        if( len( temperatures ) > 1 ) : raise Exception( 'For this TNSL protare the "--temperature" options is required as this protare has been processed to multiple temperatures.' )
        args.temperature = list( temperatures.keys( ) )[0]
    closestTemperature = -1.0
    for temperature in temperatures :
        if( abs( temperature - args.temperature ) < abs( closestTemperature - args.temperature ) ) :
            closestTemperature = temperature
            
    styleLabel = temperatures[closestTemperature]

def discreteGammaSpectrumToPDF( discreteGammaData ) :

    discreteGammaSpectrum = XYs1dModule.XYs1d( axes = spectrumAxes )
    if( args.product == PoPsID_Module.photon ) :
        for energy in sorted( discreteGammaData ) :
            energy1 = energy * ( 1.0 - args.discreteGammaResolution )
            energy2 = energy * ( 1.0 + args.discreteGammaResolution )
            height = 2.0 * discreteGammaData[energy] / ( energy2 - energy1 )        # discreteGammaData[energy] has cross section included.
            gamma = XYs1dModule.XYs1d( data = [ [ energy1, 0.0 ], [ energy, height ], [ energy2, 0.0 ] ], axes = spectrumAxes )
            discreteGammaSpectrum += gamma
    return( discreteGammaSpectrum )

def addCurves(curve1, curve2):

    if len(curve1) > 0:
        if curve1.interpolation != curve2.interpolation:
            if curve1.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve1 = curve1.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, abs(lowerUpperEpsilon), abs(lowerUpperEpsilon))
            if curve2.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve2 = curve2.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, abs(lowerUpperEpsilon), abs(lowerUpperEpsilon))
    if not curve1.areDomainsMutual(curve2):
        curve1, curve2 = curve1.mutualify(lowerUpperEpsilon, lowerUpperEpsilon, True, curve2, lowerUpperEpsilon, lowerUpperEpsilon, True)
    return curve1 + curve2

def output( MT, reactionStr, prefix, spectrum, crossSection ) :

    def write( curve, suffix, yLabel ) :

        fOut = open( os.path.join( outputDir, prefix + '_%.3d.' % MT + suffix ), 'w' )
        fOut.write( '# reaction = "%s"\n' % reactionStr )
        fOut.write( '# cross section = %.5g\n' % crossSection )
        fOut.write( '# x axes label = "Outgoing %s energy [%s]"\n' % ( args.product, args.energyUnit ) )
        fOut.write( '# y axes label = "%s"\n' % yLabel )
        if len(curve) > 0:
            fOut.write( curve.toString(format=" %17.9e %16.8e"))
        fOut.close( )

    def addPointForLogPlotting( spectrum ) :

        x1, y1 = spectrum[0]
        x2, y2 = spectrum[1]                                    # Does not handle the case where y2 is 0.0.
        if( x2 - x1 < 0.2 * x2 ) : return( spectrum )           # Mainly happends for discrete gammas.
        if( ( x1 == 0.0 ) or ( y1 == 0.0 ) ) :
            xMid = max( ( 1.0 + 1e-7 ) * x1, minimumLogPoint )
            if( xMid < ( 1 - 1e-7 ) * x2 ) :
                yMid = spectrum.evaluate( xMid )
                spectrum = XYs1dModule.XYs1d( data = spectrum, interpolation = spectrum.interpolation )
                spectrum.setValue( xMid, yMid )

        return( spectrum )

    if( len( spectrum ) == 0 ) : return

    spectrum = addPointForLogPlotting( spectrum )

    write( spectrum, 'spec', 'spectrum [%s/%s]' % ( crossSectionUnit, args.energyUnit ) )

    pdf = spectrum
    try :
        pdf = spectrum.normalize( )
    except :
        return
    write( pdf, '_pdf', 'pdf [1/%s]' % args.energyUnit )

    cdf = XYs1dModule.XYs1d( data = [ pdf.domainGrid, pdf.runningIntegral( ) ], dataForm = 'xsandys' )
    write( cdf, '_cdf', 'cdf' )

def productSpectrum(self, pid, energy, parentMultiplicity, spectrum, discreteGammaData):

    def branchingGammas(initialState, photonBranchingData, probability, discreteGammaData):
        '''This needs to be replaced by the "completePhotons" data now returned by the call to PoPsDecayMiscModule.photonBranchingData.'''

        if initialState in photonBranchingData:
            gammas = photonBranchingData[initialState]['photons']
            for branchingRatio, gammaEnergy, finalState, photonEmissionProbability in gammas:
                gammaEnergy = float(gammaEnergy)
                if gammaEnergy not in discreteGammaData:
                    discreteGammaData[gammaEnergy] = 0.0
                discreteGammaData[gammaEnergy] += branchingRatio * probability * photonEmissionProbability
                branchingGammas(finalState, photonBranchingData, branchingRatio * probability, discreteGammaData)

    if( isinstance( self.multiplicity[0], multiplicityModule.Branching1d ) ) :
        parentProduct = self.parentProduct
        photonBranchingData = PoPsDecayMiscModule.photonBranchingData( protare.PoPs, parentProduct.pid )
        branchingGammas( parentProduct.pid, photonBranchingData, 1.0, discreteGammaData )
        return( 0.0 )

    multiplicityAtEnergy = self.multiplicity.evaluate( energy )
    if specialNuclearParticleIDPoPsModule.sameSpecialNuclearParticle(pid, self.pid) and multiplicityAtEnergy != 0.0:
        spectrum2 = self.distribution.energySpectrumAtEnergy(energy, frame, twoBodyCOMResolution=args.twoBodyCOMResolution, 
                twoBodyCOM_energyResolution=args.twoBodyCOM_energyResolution, styleLabel=styleLabel, muMin=muMin, muMax=muMax)
        spectrum2 *= multiplicityAtEnergy
        spectrum2 = addCurves( spectrum, spectrum2 )
        spectrum.setData( spectrum2 )
        spectrum.interpolation = spectrum2.interpolation
    else :
        multiplicityAtEnergy = 0.0

    if( self.outputChannel is not None ) :
        multiplicityAtEnergy += outputChannelSpectrum( self.outputChannel, pid, energy, multiplicityAtEnergy * parentMultiplicity, spectrum, discreteGammaData )

    return( multiplicityAtEnergy )

def outputChannelSpectrum( self, pid, energy, parentMultiplicity, spectrum, discreteGammaData ) :

    multiplicityAtEnergy = 0.0
    for product in self.products :
        multiplicityAtEnergy += productSpectrum( product, pid, energy, parentMultiplicity, spectrum, discreteGammaData )

    if not args.skipDelayedNeutrons:
        for delayedNeutron in self.fissionFragmentData.delayedNeutrons:
            multiplicityAtEnergy += productSpectrum( delayedNeutron.product, pid, energy, parentMultiplicity, spectrum, discreteGammaData )

    return( multiplicityAtEnergy )

def reactionSpectrum(self, pid, energy, totalSpectrum, totolDiscreteGammaData, totalCrossSection, reactionSuffix):

    global reactionCounter, outputLog

    timing = timesModule.Times( )

    MT = self.ENDF_MT
    crossSection = 0.0
    multiplicityAtEnergy = 0.0
    if( ( energy >= self.domainMin ) or ( energy <= self.domainMax ) ) : crossSection = self.crossSection.evaluate( energy ) 
    reactionString = str(self) + reactionSuffix
    if args.verbose == 1: print('    %s' % reactionString)
    if args.verbose  > 1: print('    %-40s' % reactionString, crossSection)
    if( crossSection > 0.0 ) :
        spectrum = XYs1dModule.XYs1d( axes = energyAxes )
        discreteGammaData = {}
        multiplicityAtEnergy = outputChannelSpectrum( self.outputChannel, pid, energy, 1.0, spectrum, discreteGammaData )
        spectrum *= crossSection
        spectrum.axes = spectrumAxes
        totalSpectrum = addCurves( totalSpectrum, spectrum )
        if( pid == PoPsID_Module.photon ) :
            for energy in discreteGammaData :
                multiplicityAtEnergy += discreteGammaData[energy]
                discreteGammaData[energy] *= crossSection
                if( energy not in totolDiscreteGammaData ) : totolDiscreteGammaData[energy] = 0.0
                totolDiscreteGammaData[energy] += discreteGammaData[energy]

    if( ( outputDir is not None ) or args.selfCheck ) :
        if( crossSection > 0.0 ) :
            spectrum2 = addCurves( spectrum, discreteGammaSpectrumToPDF( discreteGammaData ) )
            if( outputDir is not None ) : output( MT, str( self ), '%3.3d' % reactionCounter, spectrum2, crossSection )
            integral = spectrum2.integrate() / crossSection
            for indicies in sums :
                if( indicies[0] <= MT <= indicies[1] ) :
                    sums[indicies][0] += crossSection
                    sums[indicies][1] = addCurves( sums[indicies][1], spectrum2 )
        else :
            integral = 0.0
        total = multiplicityAtEnergy
        error = 0.0
        if( ( total + integral ) != 0.0 ) : error = ( total - integral ) / max( total, integral )
        if outputDir is not None: outputLog.write('%5d  %3d  %-40s  %12.5e %12.5e %12.5e %8.1e\n' % ( reactionCounter, MT, reactionString, crossSection, total, integral, error ))
        if( abs( error ) > 2e-5 ) :
            if( ( args.verbose > 0 ) or args.selfCheck ) :
                if( args.verbose == 0 ) : print( '    %-32s' % self, crossSection )
                print( '        relative error = %8.1e' % error )

    reactionTiming.append( timing.toString( current = False, names = False ) )

    return( totalSpectrum, crossSection )

def getSpectrum(reactions, reactionSuffix):

    global reactionCounter

    spectrum = XYs1dModule.XYs1d( axes = spectrumAxes )
    discreteGammaData = {}
    totalCrossSection = 0.0
    for reactionIndex, reaction in enumerate( reactions ) :
        reactionCounter += 1
        if isinstance(reaction.crossSection.evaluated, crossSectionModule.CoulombPlusNuclearElastic):
            continue    # FIXME ignoring CP elastic scattering
        spectrum, crossSection = reactionSpectrum(reaction, productID, args.energy, spectrum, discreteGammaData, totalCrossSection, reactionSuffix)
        totalCrossSection += crossSection
    return( spectrum, discreteGammaData, totalCrossSection )

if outputDir is not None:
    outputLog = open(os.path.join(outputDir, 'index'), 'w')
    outputLog.write( 'index   MT  label                                     crossSection  multiplicity    integral    error\n' )
    outputLog.write( '-----------------------------------------------------------------------------------------------------\n' )

totalSpectraTime = timesModule.Times( )
reactionSpectrumNonDiscrete, discreteGammaData, totalCrossSection = getSpectrum(protare.reactions, '')
orphanSpectrum, dummy, dummy = getSpectrum(protare.orphanProducts, ' (orphan product)')
totalSpectraTime = totalSpectraTime.toString( current = False )

discreteGammaSpectrum = discreteGammaSpectrumToPDF( discreteGammaData )

spectrum = addCurves( addCurves( reactionSpectrumNonDiscrete, orphanSpectrum ), discreteGammaSpectrum )
if( args.plot ) :
    curves = []
    reactionTotal = reactionSpectrumNonDiscrete + discreteGammaSpectrum

    reactionSpectrumNonDiscrete = reactionSpectrumNonDiscrete.thicken( fDomainMax = 2.0, sectionSubdivideMax = 10 )
    reactionSpectrumNonDiscrete.plotLabel = 'reactions - non discrete'
    if( len( reactionSpectrumNonDiscrete ) > 0 ) : curves.append( reactionSpectrumNonDiscrete )

    discreteGammaSpectrum.plotLabel = 'reactions - discrete'
    if( len( discreteGammaSpectrum ) > 0 ) : curves.append( discreteGammaSpectrum )

    reactionTotal = reactionTotal.thicken( fDomainMax = 2.0, sectionSubdivideMax = 10 )
    reactionTotal.plotLabel = 'reactions - total '
    if( len( reactionTotal ) > 0 ) : curves.append( reactionTotal )

    orphanSpectrum = orphanSpectrum.thicken( fDomainMax = 2.0, sectionSubdivideMax = 10 )
    orphanSpectrum.plotLabel = 'orphanProducts'
    if( len( orphanSpectrum ) > 0 ) : curves.append( orphanSpectrum )

    spectrumPlot = spectrum.thicken( fDomainMax = 2.0, sectionSubdivideMax = 10 )
    spectrumPlot.plotLabel = 'total'
    if( len( spectrum ) > 0 ) : curves.append( spectrumPlot )

    if( len( curves ) > 0 ) :
        xLabel = 'Outgoing %s energy [%s]' % ( args.product, args.energyUnit )
        yLabel = '%s spectrum [%s/%s]' % ( args.product, crossSectionUnit, args.energyUnit )
        title = str( protare )
        curves[0].multiPlot( curves, xLabel = xLabel, yLabel = yLabel, title = title, domainMin = minimumLogPoint )

if outputDir is None:
    if not args.selfCheck:
        print(spectrum.toString())
else :
    with open(os.path.join(outputDir, 'info.txt'), 'w') as file:
        file.write( 'Working directory: %s\n' % os.path.realpath( os.curdir ) )
        file.write( 'Input arguments: %s\n' % ' '.join( sys.argv ) )
        file.write( 'Protare file path: %s\n' % os.path.realpath( protare.sourcePath ) )

        file.write( 'Timing:\n' )
        file.write( '  Read time: %s\n' % readTime )
        if( unitConversionTime is not None ) : file.write( '  Unit coversion time: %s\n' % unitConversionTime )
        file.write( '  Total spectra time: %s\n' % totalSpectraTime )
        file.write( '  Reaction spectra time:\n' )
        file.write( '    index | Spectrum time per reaction.\n' )
        file.write( '          |   cpu    wall\n' )
        file.write( '    %s\n' % ( 60 * '-' ) )
        for index, timing in enumerate( reactionTiming ) :
            file.write( '     %4d | %s\n' % ( index, timing ) )

    output( 1, 'total', 'total', spectrum, totalCrossSection )
    for indicies in sums :
        crossSection, spectrum, MT = sums[indicies]
        reactionStr = '%s to %s' % ( indicies[0], indicies[1] )
        prefix = '%3.3d-%3.3d' % ( indicies[0], indicies[1] )
        output( MT, reactionStr, prefix, spectrum, crossSection )
