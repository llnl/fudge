# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Issues:
    -) When delayed neutrons are present, their multiplicities are added to prompt, and prompt and total nu_bar
        are written. However, their decay constants, etc are currently not written.
    -) Only neutrons distribution data are outputted.
    -) The average heating numbers (KERMA) in the ESZ block are not calculated and are written as zeros. There is
        currently no equivalent of NJOY's HEATR module in FUDGE. A warning is emitted when the zeros are written.
"""

import time
import warnings

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import Ys1d as Ys1dModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as miscPoPsModule

from fudge import styles as stylesModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import energy as energyModule
from brownies.legacy.converting import endf_endl as endf_endlModule

from . import specialMF6 as specialMF6Module
from . import URR_probabilityTables as URR_probabilityTablesModule

floatFormatOthers = '%20.12E'
floatFormatEnergyGrid = '%21.13E'
floatFormatEnergyGrid = floatFormatOthers

def toACE(self, styleLabel, cdf_style, fileName, evaluationId, productData, delayedNeutronRateAndDatas, addAnnotation, skipURR=False, skipILF_logic=False, verbose=0):

    style = self.styles[styleLabel]
    if not isinstance( style, stylesModule.GriddedCrossSection ):
        raise TypeError( 'style labelled "%s" not a "GriddedCrossSection" style' % styleLabel )

    massUnit = 'eV/c**2'
    neutron = self.PoPs[IDsPoPsModule.neutron]
    target = self.PoPs[self.target]
    metastable = 0
    if hasattr(target, 'metaStableIndex'):
        metastable = target.metaStableIndex
    if target.id in self.PoPs.aliases: target = self.PoPs[target.pid]
    targetZ, targetA, targetZA, level = miscPoPsModule.ZAInfo( target )
    if metastable:
        targetZA += 100 * (3 + metastable)

    neutronMass = neutron.getMass( massUnit )
    targetMass = target.getMass( massUnit )
    target2NeutronMass = targetMass / neutronMass

    temperature = style.temperature
    temperature_K = temperature.getValueAs( 'K' )
    temperature_MeV = temperature.getValueAs( 'MeV/k' )

    processingTime = time.localtime( )
    strRecords = [ "%6d.%.2dc%12.7f %11.5E %.2d/%.2d/%.2d" % ( targetZA, evaluationId, target2NeutronMass,
            temperature_MeV, processingTime[1], processingTime[2], processingTime[0] ) ]

    nonFissionEnergyDependentNeutronMultiplicities = {}
    nonFissionEnergyDependentNeutronMultiplicitiesIndex = 100

    MAT = endf_endlModule.ZAAndMATFromParticleName( target.id )[1]
    MAT_ID = '  mat %d' % MAT
    evaluated = self.styles.getEvaluatedStyle( )
    HK = '%2d-%s at %.1fK from %s-%s Fudge.toACE' % \
            ( targetZ, target.id, temperature_K, evaluated.library, evaluated.version )
    HK = "%-70s" % HK[:70]
    strRecords.append( HK + MAT_ID )
    for i in range( 4 ) : strRecords.append( '      0   0.000000      0   0.000000      0   0.000000      0   0.000000' )

    NXS = 16 * [ 0 ]
    JXS = 32 * [ 0 ]
    MTR = []
    NU_prompt = None
    NU_total = None
    LQR = []
    TYR = []
    SigData = {}
    neutronAngular = []
    neutronEnergies = []

    NXS[2-1] = targetZA

    energyGrid = style.grid.values

    prior = -1
    for i1, value in enumerate(energyGrid.values):
        if value == prior:
            print('WARNING: cross section energy grid values the same at index', i1, value)
        elif abs(value - prior) < 1e-12 * value:
            print('WARNING: cross section energy grid values too close at index', i1, "%.16e" % prior, "%.16e" % value, value - prior, 1.0 - prior / value)
        prior = value

    totalXSec = Ys1dModule.Ys1d(Ys=valuesModule.Values([]))
    absorptionXSec = Ys1dModule.Ys1d(Ys=valuesModule.Values([]))

    sortedMTs = list(sorted([[MT_MTData[0], i1] for i1, MT_MTData in enumerate(productData)]))  # Sort MTs like NJOY.
    extraMTs = []
    extraData = {}

    if skipURR:
        ILF, IOA, URRPT = 0, 0, []
    else:
        ILF, IOA, URRPT = URR_probabilityTablesModule.URR_probabilityTable(self, styleLabel, fromPDF=False, verbose=verbose)
        if ILF == 4:
            crossSectionSum4 = None
            for crossSectionSum in self.sums.crossSectionSums:
                if crossSectionSum.ENDF_MT == 4:
                    crossSectionSum4 = crossSectionSum
            if crossSectionSum4 is None:
                print('WARNING: need missing MT 4 cross section to URRPT tables.')
            else:
                extraMTs.append(4)
                extraData[4] = {'ESZ': crossSectionSum4.crossSection[styleLabel], 'Q_initial': 0, IDsPoPsModule.neutron: []}

    for MT, i1 in sortedMTs :
        printMessage(verbose > 0, 'Processing MT %s' % MT)
        _MT, MTData = productData[i1]
        XSec = MTData['ESZ']
        neutronDatas = MTData[IDsPoPsModule.neutron]
        multiplicity = 0

        totalXSec += XSec

        if MT == 2:
            elasticXSec = XSec
            if len( neutronDatas ) > 1: raise Exception( 'Only one type of neutron data are supported for the elastic reaction.' )
            if len( neutronDatas ) > 0: neutronAngular.append( ( MT, neutronDatas[0]['angularData'] ) )
        else :
            MTR.append( MT )
            LQR.append( MTData['Q_initial'] )
            SigData[MT] = XSec
            if len( neutronDatas ) == 0:
                absorptionXSec += XSec
                neutronMultiplicity = 0
            else :
                NXS[5-1] += 1
                product = neutronDatas[0]['product']
                frame = neutronDatas[0]['frame']
                angularData = neutronDatas[0]['angularData']
                energyData = neutronDatas[0]['energyData']
                if MTData['isFission']:
                    neutronMultiplicity = 19
                    if len( neutronDatas ) != 1: raise Exception( 'Fission must have one neutron product.' )
                    NU_prompt = neutronDatas[0]['multiplicity']
                else :
                    if len( neutronDatas ) == 1:
                        neutronMultiplicity = neutronDatas[0]['multiplicity']
                    else :
                        neutronMultiplicity, angularData, energyData = specialMF6Module.neutrons( neutronDatas )
                    if isinstance( neutronMultiplicity, float ):
                        if neutronMultiplicity != int( neutronMultiplicity ):
                            raise Exception( 'Bad neutronMultiplicity = %e' % neutronMultiplicity )
                        neutronMultiplicity = int( neutronMultiplicity )
                    if not isinstance( neutronMultiplicity, int ):
                        nonFissionEnergyDependentNeutronMultiplicitiesIndex += 1
                        nonFissionEnergyDependentNeutronMultiplicities[MT] = [len(TYR), nonFissionEnergyDependentNeutronMultiplicitiesIndex, neutronMultiplicity, 0]
                        neutronMultiplicity = nonFissionEnergyDependentNeutronMultiplicitiesIndex
 
                    if frame == xDataEnumsModule.Frame.centerOfMass:
                        neutronMultiplicity *= -1  # Negative for COM frame.
                    if energyData is None:
                        if 50 <= MT <= 91:
                            energyData = N_nPrimeEnergyData( MT, target2NeutronMass, MTData['Q_initial'] )
                        else :
                            raise Exception('raise hell, but first implement energy for MT = %d' % MT)
                    else :
                        if isinstance( energyData, energyModule.NBodyPhaseSpace ): angularData = None
                neutronAngular.append( ( MT, angularData ) )
                neutronEnergies.append( [ MT, energyGrid[XSec.Ys.start], energyGrid[-1], energyData ] )
            TYR.append(neutronMultiplicity)

    for MT in extraMTs:
        MTData = extraData[MT]
        MTR.append(MT)
        LQR.append( MTData['Q_initial'] )
        TYR.append(0)
        SigData[MT] = MTData['ESZ']

    printMessage(verbose > 0, 'Processing delayed neutrons.')
    delayedNeutronDatas = []
    delayedNeutronRates = []
    for delayedNeutronRate, delayedNeutronData in delayedNeutronRateAndDatas :
        delayedNeutronRates.append( delayedNeutronRate )
        delayedNeutronDatas.append( delayedNeutronData )

    doDelayedNeutrons = len( delayedNeutronDatas ) > 0
    for neutronData in delayedNeutronDatas :
        if isinstance( neutronData['multiplicity'], multiplicityModule.Unspecified ): doDelayedNeutrons = False
    if not doDelayedNeutrons: delayedNeutronDatas = []

    if len( delayedNeutronDatas ) > 0:
        NU_total = NU_prompt
        for neutronData in delayedNeutronDatas :
            NU_delayed = neutronData['multiplicity']
            if NU_delayed.domainMax < NU_total.domainMax: NU_delayed = NU_delayed.dullEdges( upperEps = 1e-8 )
            NU_total = NU_total + NU_delayed

    annotates = []
    XSS = []
    
# 1) Add the ESZ block.
    NXS[3-1] = len( energyGrid )
    updateXSSInfo( 'energyGrid', annotates, XSS, energyGrid.values )
    updateXSSInfo( 'totalXSec', annotates, XSS, totalXSec )
    if len( absorptionXSec ) == 0:
        updateXSSInfo( 'absorption cross section', annotates, XSS, len( energyGrid ) * [ 0. ] )
    else:
        updateXSSInfo( 'absorption cross section', annotates, XSS, mapCrossSectionToGrid(absorptionXSec) )
    updateXSSInfo( 'elastic cross section', annotates, XSS, elasticXSec )
    warnings.warn('Average heating numbers (KERMA) are not calculated and are written to the ESZ block as zeros.')
    averageHeating = len( energyGrid ) * [ 0. ]
    updateXSSInfo( 'average heating', annotates, XSS, averageHeating )

# 2) Add the NU block.
    if NU_prompt is not None:
        JXS[2-1] = len( XSS ) + 1
        NU_prompt = NU_prompt.toACE( )
        if( NU_total is not None ) : NU_prompt.insert( 0, -len( NU_prompt ) )
        updateXSSInfo( 'NU', annotates, XSS, NU_prompt )
        if NU_total is not None:
            NU_total = NU_total.toACE( )
            updateXSSInfo( 'NU(total)', annotates, XSS, NU_total )

# 3) Add the MTR block.
    NXS[4-1] = len( MTR )
    JXS[3-1] = len( XSS ) + 1
    updateXSSInfo( 'MTR', annotates, XSS, MTR )

# 4) Add the LQR block.
    JXS[4-1] = len( XSS ) + 1
    updateXSSInfo( 'LQR', annotates, XSS, LQR )

# 5) Add the TYR block.
    JXS[5-1] = len( XSS ) + 1
    updateXSSInfo( 'TYR', annotates, XSS, TYR )

# 6 and 7) Add the LSIG and SIG blocks.
    SIG = []
    for MT, i1 in sortedMTs:
        if MT not in [2]:
            addSigData(MT, SIG, SigData[MT])
    for MT in extraMTs:
        addSigData(MT, SIG, SigData[MT])
    JXS[6-1] = len( XSS ) + 1
    LSIG = [ 1 ]
    for MT, firstNonZero, reactionSIG in SIG[:-1] : LSIG.append( LSIG[-1] + len( reactionSIG ) + 2 )
    updateXSSInfo( 'LSIG', annotates, XSS, LSIG )
    JXS[7-1] = len( XSS ) + 1
    for MT, firstNonZero, reactionSIG in SIG :
        if( MT == 18 ) : JXS[21-1] = len( XSS ) + 1
        _reactionSIG = [ firstNonZero + 1, len( reactionSIG ) ] + reactionSIG
        updateXSSInfo( 'SIG(MT=%s)' % MT, annotates, XSS, _reactionSIG )

# 8 and 9) Add the LAND and AND blocks.
    LAND, MT_AND, length = [], [], 0
    for MT, angular in neutronAngular :
        if angular is None:
            LAND.append( -1 )
        elif angular.isIsotropic( ):
            LAND.append( 0 )
        elif isinstance( angular, ( angularModule.XYs2d, angularModule.Regions2d ) ):
            if isinstance( angular, angularModule.Regions2d ):
                regions = angular
                angular = []
                for region in regions :
                    for xs_pdf_cdf in region : angular.append( xs_pdf_cdf )
            LAND.append( length + 1 )
            length += 2 * len( angular ) + 1
            energies_in = [ len( angular ) ]
            LCs = []
            Ps = []
            for xs_pdf_cdf in angular :
                energies_in.append( xs_pdf_cdf.outerDomainValue )
                LCs.append( -length - 1 )
                Ps += [ 2, len( xs_pdf_cdf ) ] + xs_pdf_cdf.xs.values.values + xs_pdf_cdf.pdf.values.values + xs_pdf_cdf.cdf.values.values
                length += 3 * len( xs_pdf_cdf ) + 2
            MT_AND.append( ( MT, energies_in + LCs + Ps ) )
        else :
            raise Exception( 'Unsupport neutron angular distribution type = %s' % type( angular ) )
    JXS[8-1] = len( XSS ) + 1
    updateXSSInfo( 'LAND', annotates, XSS, LAND )
    JXS[9-1] = len( XSS ) + 1
    for MT, AND in MT_AND : updateXSSInfo( 'AND(MT=%s)' % MT, annotates, XSS, AND )

# 10 and 11) Add the LDLW and DLW blocks.
    LDLW, MT_DLW = processEnergyData( massUnit, neutronMass, neutronEnergies, JXS, XSS, nonFissionEnergyDependentNeutronMultiplicities )
    JXS[10-1] = len( XSS ) + 1
    updateXSSInfo( 'LDLW', annotates, XSS, LDLW )
    JXS[11-1] = len( XSS ) + 1
    for MT, DLW in MT_DLW : updateXSSInfo( 'DLW(MT=%s)' % MT, annotates, XSS, DLW )
    for MT in nonFissionEnergyDependentNeutronMultiplicities :
        TYR_index, index, multiplicity, offset = nonFissionEnergyDependentNeutronMultiplicities[MT]
        offset += 101
        JXS5 = JXS[5-1]
        if XSS[JXS5+TYR_index-1] < 0: offset *= -1
        XSS[JXS5+TYR_index-1] = offset

# 23) Add LUNR block - URR probability tables.

    if len( URRPT ) > 0:
        JXS[23-1] = len( XSS ) + 1
        updateXSSInfo( 'UNR', annotates, XSS, URRPT )

# 24, 25, 26 and 27) Added the DNU, BDD, DNEDL and DNED blocks

    if len( delayedNeutronDatas ) > 0:
        skipDelayed = False
        for delayedNeutronData in delayedNeutronDatas :
            if delayedNeutronData['energyData'] is None: skipDelayed = True
        if skipDelayed: delayedNeutronDatas = []

    if len( delayedNeutronDatas ) > 0:
        NXS[8-1] = len( delayedNeutronDatas )
        delayedFissionNeutrons = []
        axes = axesModule.Axes( labelsUnits = { 1 : ( '', 'MeV' ), 0 : ( '', '' ) } )
        totolDelayedMultiplicity = XYs1dModule.XYs1d( data = [], axes = axes )
        for i1, delayedNeutronData in enumerate( delayedNeutronDatas ) :
            multiplicity = delayedNeutronData['multiplicity']
            if not isinstance( multiplicity, multiplicityModule.XYs1d ): raise Exception( 'fix me' )
            if multiplicity.interpolation != xDataEnumsModule.Interpolation.linlin:
                raise Exception('fix me')
            delayedFissionNeutrons.append( [ delayedNeutronRates[i1], multiplicity ] )
            totolDelayedMultiplicity += multiplicity
            
        energies, multiplicities = [], []
        for energy, multiplicity in totolDelayedMultiplicity :
            energies.append( energy )
            multiplicities.append( multiplicity )
        DNU = [ 2, 0, len( totolDelayedMultiplicity ) ] + energies + multiplicities
        JXS[24-1] = len( XSS ) + 1
        updateXSSInfo( 'DNU', annotates, XSS, DNU )

        BDD = []
        for decayRate, multiplicity in delayedFissionNeutrons :
            probability = multiplicity / totolDelayedMultiplicity
            probability = probability.thin( 1e-3 )
            energies, probabilities = [], []
            for energy, _probability in probability :
                energies.append( energy )
                probabilities.append( _probability )

            BDD += [ decayRate, 0, len( probability ) ] + energies + probabilities
        JXS[25-1] = len( XSS ) + 1
        updateXSSInfo( 'BDD', annotates, XSS, BDD )

        delayedNeutronEnergies = []
        for delayedNeutronData in delayedNeutronDatas :
            multiplicity = delayedNeutronData['multiplicity']
            delayedNeutronEnergies.append( [ 18, multiplicity.domainMin, multiplicity.domainMax, delayedNeutronData['energyData'] ] )
        LDLW, MT_DLW = processEnergyData( massUnit, neutronMass, delayedNeutronEnergies, JXS, XSS, {} )
        JXS[26-1] = len( XSS ) + 1
        updateXSSInfo( 'DNEDL', annotates, XSS, LDLW )
        JXS[27-1] = len( XSS ) + 1
        for MT, DLW in MT_DLW : updateXSSInfo( 'DNED(MT=%s)' % MT, annotates, XSS, DLW )

# 12) GPD block. Not needed.

# 13) MTRP block (gammas).

# 14) LSIGP block (gammas).

# 15) SIGP block (gammas).

# 16) LANDP block (gammas).

# 17) ANDP block (gammas).

# 18) LDLWP block (gammas).

# 19) DLWP block (gammas).

# 20) YP block (gammas).

# 21) FIS block. Not needed.

# 22) UNR block.

# Time to wrap it all up.
    NXS[1-1] = len( XSS )
    JXS[1-1] = 1
    JXS[22-1] = len( XSS ) + 1      # Location of the last word of XSS (i.e., its length).

    strRecords += intArrayToRecords( NXS )
    strRecords += intArrayToRecords( JXS )
    strRecords += XSSToStrings( annotates, XSS, addAnnotation, len( energyGrid ) )

    strRecords.append( '' )
    printMessage(verbose > 0, 'Writing ACE file.')
    with open( fileName, 'w' ) as fOut:
        fOut.write( '\n'.join( strRecords ) )

def updateXSSInfo( label, annotates, XSS, data ) :

    annotates.append( ( len( XSS ), label ) )
    XSS += data

def intArrayToRecords( _array ) :

    s = [ "%9d" % d for d in _array ]
    records = []
    while( len( s ) ) :
        record = []
        for i1 in range( 8 ) : record.append( s.pop( 0 ) )
        records.append( ''.join( record ) )
    return records

def XSSToStrings( annotates, XSS, addAnnotation, numberOfEnergiesInGrid ) :

    strData, record, i3 = [], [], 0
    annotates.append( ( len( XSS ), '' ) )
    i2, label = annotates[i3]
    for i1, datum in enumerate( XSS ) :
        if type( datum ) == int:
            record.append( "%20d" % datum )
        else :
            try :
                floatFormat = floatFormatOthers
                if( i1 < numberOfEnergiesInGrid ) : floatFormat = floatFormatEnergyGrid
                record.append( floatFormat % datum )
            except :
                print( i1, len(XSS), datum, annotates[i3] )
                raise
        if len( record ) == 4:
            annotate = ''
            if( i2 <= i1 ) :
                annotate = ' !'
                sep = ' '
                while( i2 <= i1 ) :
                    annotate += sep + label + ' (%s)' % ( i2 - i1 + 3 )
                    sep = ', '
                    i3 += 1
                    i2, label = annotates[i3]
            if( not addAnnotation ) :   annotate = ''
            strData.append( ''.join( record ) + annotate )
            record = []
    if( record ) : strData.append( ''.join( record ) )
    return strData

def mapCrossSectionToGrid(XSec):

    Ys = XSec.Ys
    return Ys.start * [0] + Ys.values

def addSigData( MT, SIG, SigData ) :

    SIG.append( ( MT, SigData.Ys.start, SigData.Ys.values ) )

class N_nPrimeEnergyData :

    def __init__( self, MT, A, Q ) :

        self.LAW = 3
        self.MT = MT
        self.A = float( A )
        self.Q = float( Q )

    def toACE( self, label, offset, weight, **kwargs ) :

        r1 = self.A / ( 1. + self.A )
        return [ 0, self.LAW, offset + len( weight ) + 4 ] + weight + [ abs( self.Q ) / r1, r1 * r1 ]

def processEnergyData( massUnit, neutronMass, energyDatas, JXS, XSS, nonFissionEnergyDependentNeutronMultiplicities ) :

    LDLW, MT_DLW, length = [], [], 0
    for MT, EMin, EMax, energyData in energyDatas :
        n_multiplicity = []
        if MT in nonFissionEnergyDependentNeutronMultiplicities:        # Fixup the TYR data whose abs( value ) is greater than 100.
            TYR_index, index, multiplicity, dummy = nonFissionEnergyDependentNeutronMultiplicities[MT]
            JXS5 = JXS[5-1]
            if abs( XSS[JXS5+TYR_index-1] ) != index: raise Exception( 'Neutron multiplicity for index %s not found in TYR table' % index )
            n_multiplicity = multiplicity.toACE( )
            LNU, n_multiplicity = n_multiplicity[0], n_multiplicity[1:]
            if LNU != 2: raise Exception( 'Only tabular neutron multiplicity is allowed, not LNU = %d' % LNU )
            nonFissionEnergyDependentNeutronMultiplicities[MT][3] = length
            length += len( n_multiplicity )
        LDLW.append( length + 1 )
        defaultWeight = [ 0, 2, EMin, EMax, 1.0, 1.0 ]
        label = 'MT = %s for product = "neutron"' % MT
        _DLW = energyData.toACE( label, length, defaultWeight, massUnit = massUnit, neutronMass = neutronMass )
        MT_DLW.append( ( MT, n_multiplicity + _DLW ) )
        length += len( _DLW )
    return (LDLW, MT_DLW)


def printMessage(doPrint, message):
    '''Prints *message* if *doPrint* is **True**.'''

    if doPrint:
        print(message)
