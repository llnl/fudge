# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the method toGNDS for the class endlZA. This module must only be import by the
toGNDS method of the endlZA class.

Nuclide particle name is composed of SA[_L]
        S is the element's symbol (e.g., 'He' for helium).
        A is the nucleon number (e.g., 'He4' for He-4). For natural (whatever that means) A = '0' (e.g., 'C0').
        L level designation for level greater than 0 (e.g., 'O16_e3' for the third excition state of O16).
            The level name is 'e#' where '#' is an integer (e.g., 'e2' for the second excited state).

gamma: name is simply given as 'photon'.
"""

import os
import numpy
from xml.etree import ElementTree

from pqu import PQU as PQUModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import enums as xDataEnumsModule
from xData import xDataArray as arrayModule
from xData import axes as axesModule
from xData import gridded as griddedModule
from xData import link as linkModule
from xData import values as valuesModule
from xData import uncertainties as uncertaintiesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule

from fudge import enums as enumsModule
import fudge.covariances.covarianceMatrix as covarianceMatrixModule
import fudge.covariances.covarianceSection as covarianceSectionModule
import fudge.covariances.covarianceSuite as covarianceSuiteModule
import fudge.covariances.mixed as covarianceMixedModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs import database as databasePoPsModule
from PoPs import alias as PoPsAliasModule
from PoPs.quantities import quantity as quantityPoPsModule
from PoPs.quantities import mass as massPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.chemicalElements import chemicalElement as chemicalElementPoPsModule
from PoPs.families import nuclide as nuclideModule
from PoPs.families import unorthodox as unorthodoxPoPsModule
from PoPs.fissionFragmentData import rate as rateModule

from brownies.legacy.endl import endlmisc, endlIClasses, endl2, endl_C
from . import endf_endl as endf_endlModule
from . import toGNDSMisc as toGNDSMiscModule

from fudge import physicalQuantity as physicalQuantityModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import outputChannel as outputChannelModule
from fudge import sums as sumsModule
from fudge.reactions import reaction as reactionModule
from fudge.reactions import orphanProduct as orphanProductModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection import base as baseModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import coherent as coherentModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import \
    nuclearPlusInterference as nuclearPlusInterferenceModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import RutherfordScattering as RutherfordScatteringModule

from fudge.outputChannelData import Q as QModule
from fudge.outputChannelData.fissionFragmentData import delayedNeutron as delayedNeutronModule

from fudge.productData import averageProductEnergy as averageProductEnergyModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule

from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import energy as energyModule
from fudge.productData.distributions import Legendre as LegendreModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import LLNL_angularEnergy as LLNL_angularEnergyModule
from fudge.productData.distributions import reference as referenceModule
from fudge.productData.distributions import unspecified as unspecifiedModule
from fudge.productData.distributions import photonScattering as photonScatteringModule

from brownies.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule

from .ENDFToGNDS import ENDF_ITYPE_3_6_Misc as ENDF_ITYPE_3_6_MiscModule

FUDGE_EPS = 1e-8
energyUnit = "MeV"

crossSectionAxes = crossSectionModule.defaultAxes( energyUnit )
multiplicityAxes = multiplicityModule.defaultAxes( energyUnit )
QAxes = QModule.defaultAxes( energyUnit )
averageProductEnergyAxes = averageProductEnergyModule.defaultAxes( energyUnit )
averageProductMomentumAxes = averageProductMomentumModule.defaultAxes( energyUnit, energyUnit + '/c' )
angularAxes = angularModule.defaultAxes( energyUnit )
energyAxes = energyModule.defaultAxes( energyUnit )
LLNLPointwiseAxes = LegendreModule.defaultAxes( energyUnit )
angularEnergAxes = LLNL_angularEnergyModule.defaultAxes( "MeV", 'P(energy_out|energy_in,mu)' )

metaTargets = { 21046 : ( 1, 2 ), 95242 : ( 1, 2 ) }

def uncorrelatedIsotropicGivenEnergySpectrum( info, product, ENDL_energySpectrum ) :

    angularSubform = angularModule.Isotropic2d( )
    energySubform = energyModule.XYs2d(axes=energyAxes, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
    for energy, probability_list in ENDL_energySpectrum.data :
        energySubform.append( energyModule.XYs1d( data = probability_list, axes = energyAxes, outerDomainValue = energy ) )
    form = uncorrelated(info.style, xDataEnumsModule.Frame.lab, angularSubform, energySubform)

    product.distribution.add( form )

def covarianceDict(rootdir):
    """
    Creates a dictionary of covariance files in rootdir, keyed by the pair (C,S).
    """
    covar_dictionary = {}
    for folder, dirs, files in os.walk(rootdir):
        for covFile in files:
            if covFile.endswith('_cov.xml'):
                C = int(covFile.split('c')[1].split('i')[0])
                S = int(covFile.split('s')[1].split('_')[0])
                if (C,S) not in covar_dictionary:
                    covar_dictionary[(C,S)] = []
                covar_dictionary[(C,S)].append(covFile)

    return ( covar_dictionary )

def parse_endl_covariance(covFile): 
    """
    Returns a list of (energy bins,  covariance matrix, covariance_type, [enminmax]) tuples.
    The list may be empty (some ENDL cov.xml files are empty).
    """
    xdoc = ElementTree.parse(covFile)
    root = xdoc.getroot()
    ebins, covariances, covariance_types, enminmax = [],[],[],[]
    # should I just set the I# here and return it?
    # it should be based on the covariance type (cross section vs. energy dist. vs. ???)
    if root.tag in ('cross_section_covariance', 'nubar_covariance'):
        for child in root:
            if child.tag == 'covariance':
                xmltype = child.get('type')
                if xmltype == 'relative':
                    covariance_types.append( xmltype )
                else:
                    covariance_types.append( 'absolute' )
                matrixNode = child.find('matrix')
                rows = int(matrixNode.get('dim1'))
                columns = int(matrixNode.get('dim2'))
                covariances.append( numpy.array(list(map(float, matrixNode.text.split()))).reshape(rows,columns) )
            elif child.tag == 'histogram':
                ebins.append( numpy.array( list(map(float, child.text.split())) ) )
        return list(zip(ebins, covariances, covariance_types))
    
    elif root.tag in ('energy_dist_covariance'): # Prompt fission neutron spectrum covariance
        for child in root:
            if child.tag == 'ein_grid':
                matrixNode = child.find('histogram')
                enminmax.append( numpy.array( list(map(float, matrixNode.text.split())) ) )
            elif child.tag =='outgoing_energy_region':
                covarSec = child.find('covariance')
                xmltype = covarSec.get('type')
                if xmltype == 'relative':
                    covariance_types.append( xmltype )
                else:
                    covariance_types.append( 'absolute' )
                matrixNode = covarSec.find('matrix')
                rows = int(matrixNode.get('dim1'))
                columns = int(matrixNode.get('dim2'))
                covariances.append( numpy.array(list(map(float, matrixNode.text.split()))).reshape(rows,columns) )
                eoutGrid = child.find('eout_grid')
                ebinNode = eoutGrid.find('histogram')
                ebins.append( numpy.array( list(map(float, ebinNode.text.split())) ) )
        return list(zip(ebins, covariances, covariance_types, enminmax))

    else:
        raise NotImplementedError("Unsupported covariance type %s" % root.tag)
# nubar: rowdata needs to point to outgoing neutron multiplicity 

def getXYInterpolation( data ) :

    if data.columns != 2:
        raise Exception('Only 2 column data supported, # of columns = %s for %s' % (data.columns, str(data)))
    if data.interpolation == 0:
        return xDataEnumsModule.Interpolation.linlin
    if data.interpolation == 1:
        return xDataEnumsModule.Interpolation.linlog
    if data.interpolation == 2:
        return xDataEnumsModule.Interpolation.loglin
    if data.interpolation == 3:
        return xDataEnumsModule.Interpolation.loglog
    raise Exception('Unsupported interpolation = "%s" for %s' % (data.interpolation, str(data)))

def returnConstantQ( label, Q_MeV, crossSection ) :

    return( toGNDSMiscModule.returnConstantQ( label, Q_MeV, crossSection ) )

def uncorrelated( style, frame, angularSubform, energySubform ) :

    _angularSubform = uncorrelatedModule.AngularSubform( angularSubform )
    _energySubform = uncorrelatedModule.EnergySubform( energySubform )
    return( uncorrelatedModule.Form( style, frame, _angularSubform, _energySubform ) )

def toGNDS(self, evaluationLibrary, evaluationVersion, formatVersion=GNDS_formatVersionModule.default, excludeAverageProductData=True,
        TNSL_include_all_reactions=True, verbose=0, specialNuclearParticleID=specialNuclearParticleIDPoPsModule.Mode.nuclide):
    """
    Returns an reactionSuite for self where self is an endlZA class.
    @param self: brownies.legacy.endl.endlZAClass instance to be translated to GNDS
    @param evaluationLibrary: library name (str)
    @param evaluationVersion: library version (str)
    @param excludeAverageProductData: if false, include C=10 / C=11 data in translation
    @param TNSL_include_all_reactions: if false and endlZA is a TNSL target, only translate thermal scattering reactions
    @param verbose: verbosity level (int)

    @return:    reactionSuite and covarianceSuite
    """
#
#   A discrete N-body reaction is one where all the yo's are known and they have energy independent multiplicities.
#

    def getDistributionIs( yo, yosDatas ) :
        """Gets the list of I-values matching the list [ 1, 3, 4, 7, 9 ] for yo."""

        Is = []
        for data in yosDatas[yo] :
            if( data.I in ( 1, 3, 4, 7, 9 ) ) : Is.append( data.I )
        return( Is )

    def getIDistributionData( yo, yosDatas, I, S = None ) :
        """Return the data for the requested yo/I/S or None if not data present."""

        for data in yosDatas[yo] :
            if( S is None ) : S = data.S
            if( ( data.I == I ) and ( data.S == S ) ) : return( data )
        return( None )

    def getAllDistributionIs( yoDatas, IMax = 9, IExtras = () ) :
        """For a yo get all distribution I-values (i.e., I in [1,9] and in IExtras)."""

        Is = []
        for yoData in yoDatas :
            I = yoData.I
            if( ( ( I in IExtras ) or ( I <= IMax ) ) and ( I not in Is ) ) : Is.append( I )
        Is.sort( )
        return( Is )

    def isYoAndResidaulDataTheSame( yo, yosDatas ) :
        """Test if yo and the residual have the same distribution data (i.e., I = 1, 3, 4, 7 and 9)."""

        yoData = yosDatas[yo]
        resData = yosDatas[yo+10]
        if( ( len( yoData ) == 0 ) or ( len( resData ) == 0 ) ) : return( False )
        yoIs = getDistributionIs( yo, yosDatas )
        resIs = getDistributionIs( yo + 10, yosDatas )
        if( yoIs != resIs ) : return( False )
        for I in yoIs :
            yoIData = getIDistributionData( yo, yosDatas, I )
            resIData = getIDistributionData( yo + 10, yosDatas, I )
            if( yoIData.data != resIData.data ) : return( False )
        return( True )

    def getMultiplicityYos( self, yos, residuals, yosDatas ) :
        """Gets the (multiplicty, yo) pairs and adjust residuals if needed."""

        if( len( residuals ) > 1 ) :
            raise Exception( 'len( residuals ) > 1 not supported: residual = %s' % repr(residuals) )
        yoOld = None
        mYos = []
        n = 0
        for yo in yos :
            if( yo == yoOld ) :
                n += 1
            else :
                if( yoOld is not None ) : mYos.append( [ n, yoOld ] )
                yoOld = yo
                n = 1
        if( len( yos ) > 0 ) :
            if( residuals[0] == yo ) :
                rYo = endl2.ZAToYo(residuals[0])
                if( isYoAndResidaulDataTheSame( rYo, yosDatas ) ) :
                    residuals = residuals[1:]           # residuals should now be empty.
                    n += 1
                    yosDatas[rYo + 10] = []
                    for yosData in yosDatas[rYo] :
                        if( yosData.I == 10 ) : yosData.set( yosData * ( n / float( n - 1 ) ) )

            mYos.append( [ n, yo ] )
        if( len( residuals ) > 0 ) : mYos.append( [ 1, residuals[0] ] )
        return( mYos )

    def addDistributionDataAndRemove( particle, yo, yosDatas, crossSection, promptNeutronParticle = None ) :

        def angular( data ) :

            subform = angularModule.XYs2d( axes = angularAxes )
            for i1, EMuP in enumerate( data.data ) :
                E1, muP = EMuP
                subform[i1] = angularModule.XYs1d( data = muP, outerDomainValue = E1, axes = angularAxes )
            return( subform )

        def energy( data ) :

            subform = energyModule.XYs2d(axes=energyAxes, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
            for E1, EpP in data.data[0][1] : subform.append( energyModule.XYs1d( data = EpP, outerDomainValue = E1, axes = energyAxes ) )
            return( subform )

        def LLNLLegendrePointwise( data ) :

            axes = LLNLPointwiseAxes
            subform = LegendreModule.LLNLPointwise( axes )
            for l, EEpPs in data.data :
                w_xys = multiD_XYsModule.XYs2d(axes=axes, outerDomainValue=l, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
                for index, E_EpPs in enumerate( EEpPs ) :
                    xPrior = -1
                    E1, EpPs = E_EpPs
                    for xy in EpPs :
                        if( xy[0] == xPrior ) : xy[0] *= ( 1 + FUDGE_EPS )
                        xPrior = xy[0]
                    w_xys.append( XYs1dModule.XYs1d( data = EpPs, axes = axes, outerDomainValue = E1 ) )
                subform.append( w_xys )
            return( subform )

        def LLNLLegendrePointwise_L0_only( info, data ) :
            """Has only L=0 form, can be converted to uncorrelated angular (isotropic) and energy distributions."""

            angularSubform = angularModule.Isotropic2d( )
            axes = energyAxes
            energySubform = energyModule.XYs2d(axes=axes, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
            for index, ( energy, probability_list ) in enumerate( data.data[0][1] ) :
                energySubform.append( energyModule.XYs1d( data = probability_list, axes = axes, outerDomainValue = energy ) )
            return uncorrelated(info.style, xDataEnumsModule.Frame.lab, angularSubform, energySubform)

        def LLNLAngularEnergy( data ) :

            axes = angularEnergAxes
            subform = LLNL_angularEnergyModule.XYs3d(axes=axes, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
            for i, E_MuEpPs in enumerate( data.data ) :
                w_xys = LLNL_angularEnergyModule.XYs2d( outerDomainValue = E_MuEpPs[0], axes = axes )
                for j, Mu_EpPs in enumerate( E_MuEpPs[1] ) :
                    w_xys.append( LLNL_angularEnergyModule.XYs1d( data = Mu_EpPs[1], axes = axes, outerDomainValue = Mu_EpPs[0] ) )
                subform.append( w_xys )
            return( subform )

        def getSubform( I, datas, func ) :

            for i1, data in enumerate( datas ) :
                if( data.I == I ) :
                    del datas[i1]
                    return( func( data ) )
            raise Exception( 'I = %d not in datas' % I )

        if( yo not in yosDatas ) : return
        if( excludeAverageProductData ) :
            _datas = []
            for data in yosDatas[yo] :
                if( data.I not in [ 10, 11, 13 ] ) : _datas.append( data )
            yosDatas[yo] = _datas
        datas = yosDatas[yo]
        Is = getAllDistributionIs( datas, IMax = 4, IExtras = [ 941, 942 ] )
        form = None
        doubleDifferentialCrossSectionForm = None
        if( Is == [] ) :
            if( len( datas ) ) :
                if( ( len( datas ) == 1 ) and ( datas[0].S == 7 ) and ( datas[0].I == 7 ) ) :   # Special case for endl I, S = 7, 7 with no distribution data.
                    form = referenceModule.Form( link = promptNeutronParticle.distribution[info.style], label = info.style )
                else :
                    print()
                    for data in datas : print(data)
                    raise Exception( "Unsupported data for particle %s: I's = %s." % (particle.pid, repr(Is)))
        elif( Is == [ 941 ] ) :
            for i1 in range( len( datas ) ) :
                if( datas[i1].I == 941 ) : break
            axes = baseModule.defaultAxes( 'coherent form factor', energyUnit = '1/cm' )
            formFactor = coherentModule.Regions1d( axes = axes )
            formFactor.append( coherentModule.XYs1d( data = datas[i1].data[:2], axes = axes ) )
            formFactor.append(coherentModule.XYs1d(data=datas[i1].data[1:], axes=axes, interpolation=xDataEnumsModule.Interpolation.loglog))
            formFactor = coherentModule.FormFactor( formFactor )
            doubleDifferentialCrossSectionForm = coherentModule.Form( particle.pid, info.style, xDataEnumsModule.Frame.lab, formFactor, None, None )
            form = photonScatteringModule.CoherentPhotonScattering.Form( label = info.style, link = doubleDifferentialCrossSectionForm )
            del datas[i1]
        elif( 942 in Is ) :
            for i1 in range( len( datas ) ) :
                if( datas[i1].I == 942 ) : break
            axes = baseModule.defaultAxes( 'incoherent scattering function', energyUnit = '1/cm' )
            scatteringFactor = incoherentModule.Regions1d( axes = axes )
            scatteringFactor.append( incoherentModule.XYs1d( data = datas[i1].data[:2], axes = axes ) )
            scatteringFactor.append(incoherentModule.XYs1d(data=datas[i1].data[1:], axes=axes, interpolation=xDataEnumsModule.Interpolation.loglog))
            scatteringFactor = incoherentModule.ScatteringFactor(scatteringFactor)
            doubleDifferentialCrossSectionForm = incoherentModule.Form( particle.pid, info.style, xDataEnumsModule.Frame.lab, scatteringFactor )
            form = photonScatteringModule.IncoherentPhotonScattering.Form( label = info.style, link = doubleDifferentialCrossSectionForm )
            del datas[i1]
        elif( Is == [ 1 ] ) : 
            if( datas[0].S == 3 ) :
                energySubform = energyModule.DiscreteGamma( datas[0].getX1( ), crossSection.domainMin,
                        crossSection.domainMax, axes = energyAxes )
                angularSubform = getSubform( 1, datas, angular )
                form = uncorrelated( info.style, xDataEnumsModule.Frame.lab, angularSubform, energySubform )
            elif( 1.99e-12 < datas[0].X1 < 2.01e-12 ) :
                # Thermal neutron scattering law coherent elastic, needs to be translated as uncorrelated
                from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import \
                    base as baseTNSLModule
                productFrame = xDataEnumsModule.Frame.lab
                angularSubform = getSubform( 1, datas, angular )
                energySubform = baseTNSLModule.energyDelta2d(angularSubform.domainMin, angularSubform.domainMax, FUDGE_EPS, energyUnit )
                form = uncorrelated( info.style, xDataEnumsModule.Frame.lab, angularSubform, energySubform )
            else :
                for data in datas :
                    if( data.I == 1 ) : productFrame = xDataEnumsModule.Frame.centerOfMass
                subform = getSubform( 1, datas, angular )
                form = angularModule.TwoBody( info.style, productFrame, subform )
        elif( Is == [ 4 ] ) :
            for i1 in range( len( datas ) ) :
                if( datas[i1].I == 4 ) : break
            if( ( len( datas[i1].data ) == 1 ) and ( datas[i1].data[0][0] == 0 ) ) :
                form = LLNLLegendrePointwise_L0_only( info, datas[i1] )
                del datas[i1]
            else:
                subform = getSubform( 4, datas, LLNLLegendrePointwise )
                form = LegendreModule.Form( info.style, xDataEnumsModule.Frame.lab, subform )
        elif( Is == [ 1, 3 ] ) :
            angularSubform = getSubform( 1, datas, angular )
            if any([xys.domainMin != -1 for xys in angularSubform]):
                angularSubform.interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase
            angularSubform = LLNL_angularEnergyModule.LLNLAngularOfAngularEnergySubform( angularSubform )
            angularEnergySubform = getSubform( 3, datas, LLNLAngularEnergy )
            angularEnergySubform = LLNL_angularEnergyModule.LLNLAngularEnergyOfAngularEnergySubform( angularEnergySubform )
            form = LLNL_angularEnergyModule.LLNLAngularEnergyForm( info.style, xDataEnumsModule.Frame.lab,
                    angularSubform, angularEnergySubform )
        elif( Is == [ 1, 4 ] ) :
            angularSubform = getSubform( 1, datas, angular )
            energySubform = getSubform( 4, datas, energy )
            form = uncorrelated( info.style, xDataEnumsModule.Frame.lab, angularSubform, energySubform )
        else :
            raise Exception( "I's = %s not supported for particle %s." % (repr(Is), particle.pid))

        for idx in range( len( datas ) - 1, -1, -1 ) :
            data = datas[idx]
            if( data.I in [ 7, 9, 10, 13, 941, 942 ] ) :
                if( excludeAverageProductData and ( data.I in [ 10, 13 ] ) ) : continue
                interpolation = getXYInterpolation( data )
                if( data.I in [ 7, 9 ] ) : 
                    axes = multiplicityAxes
                    particle.multiplicity.remove( info.style )
                    particle.multiplicity.add( multiplicityModule.XYs1d( data = data.data, label = info.style, 
                            axes = axes, interpolation = interpolation ) )
                else :          # for I = 10 and 13
                    if( data.I == 10 ) :
                        axes = averageProductEnergyAxes
                        dataForm = averageProductEnergyModule.XYs1d( label = info.style, axes = axes, 
                                data = data.data, interpolation = interpolation )
                        particle.averageProductEnergy.add( dataForm )
                    elif( data.I == 13 ) :
                        axes = averageProductMomentumAxes
                        dataForm = averageProductMomentumModule.XYs1d( label = info.style, axes = axes, 
                                data = data.data, interpolation = interpolation )
                        particle.averageProductMomentum.add( dataForm )
                del datas[idx]

        if( not form is None ) : particle.distribution.add( form )

        return( doubleDifferentialCrossSectionForm )

    def addRecoilDistributionAndRemove( particle, recoilPartner, yo, yosDatas ) :
        """
        For 2-body reactions involving two light products, ENDL explicitly stores distributions for both
        rather than computing one as the recoil of the other. When translating, discard the 2nd distribution
        and store it as 'recoil' instead.
        """
        if yo in yosDatas:
            assert yosDatas[yo][0].I == 1
            yosDatas[yo] = []
        productFrame = xDataEnumsModule.Frame.centerOfMass
        subform = angularModule.Recoil( link = recoilPartner.distribution[info.style], relative=True )
        form = angularModule.TwoBody( info.style, productFrame, subform )
        particle.distribution.add( form )

    def makeNBodyChannelFrom_mYos( channelClass, genre, mYos, ZAsYos, yosDatas, Q_MeV, crossSection,
            specialCase = None, levelIndex = None, process = None ) :
        """
        This routine does not support yo's (including the residual) in an excited state. This is consistent with endl99.
        """

        metaZA = 0
        metaName = None
        if( ( len( yosDatas[7] ) > 0 ) and ( hasattr( yosDatas[7][0], 'metaInfo' ) ) ) : metaZA, metaName = yosDatas[7][0].metaInfo
        s = ' '
        channel = channelClass( genre, process = process )
        channel.Q.add( returnConstantQ( info.style, Q_MeV, crossSection ) )
        i1 = 0
        n1 = len( mYos ) - 1

        try :
            index = mYos.index( [ 1, 4008 ] )
        except :
            index = -1
        if( index > -1 ) :
            if( [ 1, 2004 ] in mYos ) :
                mYos[index] = [ 2, 2004 ]
            else :
                if( ( 6 in yosDatas ) and ( 16 in yosDatas ) ) :
                    mYos[index] = [ 1, 2004 ]
                    mYos.insert( index, [ 1, 2004 ] )
                else :
                    mYos[index] = [ 2, 2004 ]

        He4_counter = 0
        for m, ZA in mYos :
            if( ZA == 2004 ) : He4_counter += 1
            decayChannel = None
            if( m > 1 ) :
                particle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, ZA ), crossSection, multiplicity = m, outputChannel = decayChannel )
            else :
                if( ( metaZA == ZA ) and( metaName is not None ) ) :
                    particle = metaName
                else :
                    particle = toGNDSMiscModule.getTypeName( info, ZA )
                
                particle = toGNDSMiscModule.newGNDSParticle( info, particle, crossSection, outputChannel = decayChannel )

            s = ' + '
            if( ZA in ZAsYos ) :
                if( ( specialCase == 'S8' ) and ( ZA == 1 ) ) :
                    yo = 11
                else :
                    yo = ZAsYos[ZA]
                    ICounter = 0
                    for data in yosDatas[yo] :
                        if( data.I < 10 ) : ICounter += 1
                    if( ICounter == 0 ) :
                        if( yo == 6 ) and False :
                            if( He4_counter == 2 ) : yo += 10
                        else :
                            yo += 10
                addDistributionDataAndRemove( particle, yo, yosDatas, crossSection )
            channel.products.add( channel.products.uniqueLabel( particle ) )
            i1 += 1

        return( channel )

    isThermalNeutronScatteringLawTarget = self.ZA in [ 1701, 1801, 1901, 1902, 4809, 4909, 6912, 8916, 40900 ]
    targetZA = self.ZA
    targetName = None
    if( isThermalNeutronScatteringLawTarget ) :
        targetZA -= 800
        if( self.ZA in [ 1901, 1902, 4909, 6912, 8916 ] ) : targetZA -= 100
        if( self.ZA in [ 1701 ] ) : targetZA += 100
        targetName = { 1701 : "HinZrH", 1801 : "HinH2O", 1901 : "HinCH2",
                       1902 : "DinD2O", 4809 : "Be-metal", 4909 : "BeinBeO",
                       6912 : "crystalline-graphite", 8916 : "OinBeO", 40900 : "ZrinZrH" }[self.ZA]

        I0 = self.findData( C = 10, I = 0 )
        for energy, xSec in I0.data :
            if( xSec > 0 ) : break
            TNSL_EMax = energy
    targetZ = targetZA // 1000

    styleName = 'eval'
    yosDatas = {  0 : [],  1 : [],  2 : [],  3 : [],  4 : [],  5 : [],  6 : [],  7 : [],  8 : [],  9 : [],
                     10 : [], 11 : [], 12 : [], 13 : [], 14 : [], 15 : [], 16 : [], 17 : [], 18 : [], 19 : [] }
    yosZAs = {}
    for yo in range( 1, 17 ) : yosZAs[yo] = abs(endl2.yoToZA(yo))
    ZAsYos = {}
    for yo in range( 1, 7 ) : ZAsYos[endl2.yoToZA(yo)] = yo
    ZAsYos[7] = 7
    ZAsYos[9] = 9
    if self.projectileName == IDsPoPsModule.neutron:
        projectileName = IDsPoPsModule.neutron
    elif self.projectileName == 'gamma':
        projectileName = IDsPoPsModule.photon
    else:
        projectileName = specialNuclearParticleIDPoPsModule.nuclideID(self.projectileName)
    ENDLCS_To_ENDFMT = endf_endlModule.ENDLCS_To_ENDFMT(projectileName)
    projectileID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(projectileName, specialNuclearParticleID)

    level = None
    I0 = self.findDatas( I = 0 )[0]
    if( I0.getELevel( ) > 0 ) : level = I0.getELevel( )

    info = toGNDSMiscModule.Infos( formatVersion, styleName )
    info.PoPsLabel = styleName
    info.PoPs = databasePoPsModule.Database( 'protare_internal', '1.0', formatVersion = info.formatVersion )
    info.ENDFconversionFlags = ENDFconversionFlagsModule.ENDFconversionFlags()
    info.specialNuclearParticleID = specialNuclearParticleID

    I0s = self.findDatas( C = 46, S = 0, I = 0 )
    if( len( I0s ) == 2 ) :
        I0m1, I0g = I0s
        if( I0m1.getQ( ) > I0g.getQ( ) ) : I0m1, I0g = I0g, I0m1
        metaLevel = I0g.getQ( ) - I0m1.getQ( )
        dataM1s = self.findDatas( C = 46, S = 0, Q = I0m1.getQ( ) )
        metaZA = self.ZA + 1
        if( metaZA not in metaTargets ) : raise Exception( 'Need meta level index for ZA = %s' % metaZA )

        excitedLevel = toGNDSMiscModule.getPoPsParticle( info, metaZA, name = None, levelIndex = metaTargets[metaZA][1], level = metaLevel, levelUnit = 'MeV' )
        metaName = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(excitedLevel.id, metaTargets[metaZA][0])
        info.PoPs.add(PoPsAliasModule.MetaStable(metaName, excitedLevel.id, metaTargets[metaZA][0]))
        for dataM1 in dataM1s : dataM1.metaInfo = [ metaZA, metaName ]

    projectile = toGNDSMiscModule.getTypeName(info, abs(endl2.yoToZA(self.yi)))

    inelasticC = { 1 : 11, 2 : 40, 3 : 41, 4 : 42, 5 : 44, 6 : 45, 7 : -1, 9 : -1 }[self.yi]
    CsAndSs, residualExcitationIndexLevels = {}, {}
    if( level is not None ) : residualExcitationIndexLevels[targetZA] = [ [ None, level ] ]
    reactionsDatas = self.findReactionsDatas( )
    if( 1 in self.CList( ) ) : reactionsDatas += [ self.findDatas( C = 1 ) ]

    compoundZA = endl2.yoToZA(self.yi) + targetZA

    for reactionDatas in reactionsDatas :               # Determine a unique set of levels for excited residuals.
        C, S = reactionDatas[0].C, reactionDatas[0].S
        if( C not in CsAndSs ) : CsAndSs[C] = []
        if( S not in CsAndSs[C] ) : CsAndSs[C].append( S )
        if( S == 1 ) :
            yos = endl_C.endl_C_yoInfo(C)
            residualZA = compoundZA - yos[0]
            if( residualZA not in residualExcitationIndexLevels ) : residualExcitationIndexLevels[residualZA] = []
            residualExcitationIndexLevels[residualZA].append( [ None, reactionDatas[0].X1 ] )
    for residualZA in residualExcitationIndexLevels :
        residualExcitationIndexLevels[residualZA].sort( )
        indexLevels, indexOffset = {}, 1
        if( residualExcitationIndexLevels[residualZA][0][1] == 0 ) : indexOffset = 0       # Should only happend for meta-stable targets.
        for index, indexAndLevel in enumerate( residualExcitationIndexLevels[residualZA] ) : 
            indexLevel = index + indexOffset
            indexLevels[indexAndLevel[1]] = indexLevel
        residualExcitationIndexLevels[residualZA] = indexLevels

    levelIndex = None
    if( level is not None ) : levelIndex = residualExcitationIndexLevels[targetZA][level]
    if( self.yi == 7 and targetZA not in (99120, 99125) ) :
        # target is a chemicalElement, but target mass is stored in an A=0 nuclide
        targetID = chemicalElementMiscPoPsModule.symbolFromZ[targetZ]
        target = nuclideModule.Particle.buildFromClassAndRawData(targetID + "0",
                mass=(info.massTracker.getMassAMU(targetZA), 'amu'),
                label=info.PoPsLabel)
        info.PoPs.add( target )
        info.targetID = targetID
    elif( self.yi == 9 ) :
        info.targetID = chemicalElementMiscPoPsModule.symbolFromZ[targetZ]
        targetName = chemicalElementMiscPoPsModule.nameFromZ[targetZ]
        info.PoPs.add( chemicalElementPoPsModule.ChemicalElement( info.targetID, targetZ, targetName ) )
    elif isThermalNeutronScatteringLawTarget:
        info.targetID = targetName

        # TNSL target (and mass) needs to be in PoPs for processProtare and other codes to work
        TNSL_target = unorthodoxPoPsModule.Particle( targetName )
        mass = self.findDatas()[0].getMass()
        TNSL_target.mass.add(
            massPoPsModule.Double( info.PoPsLabel, mass, quantityPoPsModule.stringToPhysicalUnit( 'amu' ) ) )
        info.PoPs.add(TNSL_target)
    else :
        target = toGNDSMiscModule.getTypeName( info, targetZA, name = targetName, level = level, levelIndex = levelIndex )
        info.targetID = target.id

    if( levelIndex is not None and levelIndex > 0 ):
        metastableIndex = 1 # FIXME can this be read from ENDL instead?
        aliasName = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(target.id, metastableIndex)
        info.PoPs.add(PoPsAliasModule.MetaStable(aliasName, target.id, metastableIndex))
        info.targetID = aliasName

    maxDate = 19590101
    allData = self.findDatas( )
    for oneDatum in allData :
        if( 0 <= oneDatum.yo <= 16 ) :
            if( oneDatum.I in [ 80, 81, 84, 89, 90, 91, 92 ] ) : continue
            if( oneDatum.I not in [ 10, 11, 13, 20 ] ) :
                date = oneDatum.getDate( ) + 19 * 1000000
                if( date < 600000 ) : date += 1000000   # Must treat dates after 2000 different then those before (600000 is arbitrary).
                maxDate = max( date, maxDate )
    maxDate = str( maxDate )
    maxDate = '%s-%s-%s' % ( maxDate[:4], maxDate[4:6], maxDate[6:] )
    projectileDomain = stylesModule.ProjectileEnergyDomain( 1e-11, 20, 'MeV' )  # will be overwritten after parsing all reactions
    evaluatedStyle = stylesModule.Evaluated( info.style, '',
            physicalQuantityModule.Temperature( PQUModule.PQU_float.surmiseSignificantDigits( I0.getTemperature( ) ), 'MeV/k' ),
            projectileDomain,
            evaluationLibrary, evaluationVersion, date = maxDate )

    evaluation = "%s-%s" % ( evaluationLibrary, evaluationVersion ) if evaluationVersion != '' else evaluationLibrary
    interaction = enumsModule.Interaction.nuclear
    if( projectile.id == IDsPoPsModule.photon ) : interaction = enumsModule.Interaction.atomic
    if( ( self.A / 100 ) > 3 ) : interaction = enumsModule.Interaction.LLNL_TNSL
    reactionSuite = reactionSuiteModule.ReactionSuite( projectileID, info.targetID, evaluation, style = evaluatedStyle, PoPs = info.PoPs,
            interaction = interaction, formatVersion = formatVersion )
    covarianceSuite = None  # initialized later if covariances found
    info.PoPs = reactionSuite.PoPs
    info.reactionSuite = reactionSuite

    reactionsByCS = {}
    I20s = []                                           # Some I20 data has bad Q values, so let's ignore all I20 data.
    for i, reactionDatas in enumerate( reactionsDatas ) :
        for j, reactionData in enumerate( reactionDatas ) :
            if( reactionData.I == 20 ) : I20s.append( [ i, j ] )
    I20s.reverse( )
    for i, j in I20s :
        del reactionsDatas[i][j]
        if( len( reactionsDatas[i] ) == 0 ) : del reactionsDatas[i]

    channelList, delayedNeutrons, c55s = [], {}, []
    iChannel = 0

    if( self.yi == 7 ) :            # Make reaction order the same as for ENDF.
        C73s = []
        C_Others = []
        for reactionDatas in reactionsDatas :
            if( reactionDatas[0].C == 73 ) :
                C73s.append( reactionDatas )
            else :
                C_Others.append( reactionDatas )
        reactionsDatas = C_Others + C73s

    if( self.yi == 9 ) :            # Make reaction order the same as for ENDF.
        C81s = []
        C_Others = []
        for reactionDatas in reactionsDatas :
            if( reactionDatas[0].C == 81 ) :
                C81s.append( reactionDatas )
            else :
                C_Others.append( reactionDatas )
        reactionsDatas = C_Others + C81s

    crossSectionC8 = None
    for reactionDatas in reactionsDatas :
        doubleDifferentialCrossSectionForm = None
        isThermalNeutronScatteringLawReaction = False
        C, S = reactionDatas[0].C, reactionDatas[0].S   # Phase 1

        if( isThermalNeutronScatteringLawTarget ) :
            I0 = None
            for reactionData in reactionDatas :
                if( reactionData.I == 0 ) : I0 = reactionData
            if( ( C, S ) == ( 11, 1 ) and ( 0.9e-12 < I0.getX1( ) < 2.1e-12 ) ) :
                isThermalNeutronScatteringLawReaction = True
            else:
                if( not TNSL_include_all_reactions ) :
                    continue
        elif C == 73 and S == 91:
            # fluorescence / atomic relaxation data, skip for now
            continue
        elif( C == 81 ) :
            if( verbose > 0 ) : print( '   ', C, S )
            X1_MTs = {}
            for index, X1 in enumerate( [ 1, 3, 5, 6, 8, 10, 11, 13, 14, 16, 18, 19, 21, 22, 24, 25, 27, 29, 30, 32, 33, 35, 36, 41, 43, 44, 46, 47, 58 ] ) : X1_MTs[X1] = 534 + index
            X1 = int( reactionDatas[0].X1 )
            try :
                MT = X1_MTs[X1]
            except :
                print( MT, X1 )
                continue
                raise
            targetID = info.targetID + '{%s}' % ENDF_ITYPE_3_6_MiscModule.MT_AtomicConfigurations[MT]

            for data in reactionDatas :
                if( data.I == 0 ) :
                    I0 = data
                elif( data.I == 21 ) :
                    I21 = data

            reaction = reactionModule.Reaction(str(MT), enumsModule.Genre.NBody, ENDF_MT=MT)
            crossSection = crossSectionModule.XYs1d( data = I0.data, interpolation = getXYInterpolation( I0 ), axes = crossSectionAxes, label = info.style )
            reaction.crossSection.add( crossSection )

            outputChannel = reaction.outputChannel
            outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )

            electron1 = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 9 ), crossSection )
            uncorrelatedIsotropicGivenEnergySpectrum( info, electron1, I21 )
            outputChannel.products.add( outputChannel.products.uniqueLabel( electron1 ) )

            electron2 = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 9 ), crossSection )
            outputChannel.products.add( outputChannel.products.uniqueLabel( electron2 ) )

            residual = toGNDSMiscModule.newGNDSParticle( info, targetID, crossSection )
            outputChannel.products.add( outputChannel.products.uniqueLabel( residual ) )

            reaction.updateLabel()
            reactionSuite.reactions.add( reaction )
            continue

        channelProcess = None
        if( S == 0 ) :
            if( 1 in CsAndSs[C] ) :
                channelProcess = outputChannelModule.Processes.continuum
            else :
                if( C == inelasticC ) : channelProcess = outputChannelModule.Processes.continuum
        I0, I12, I20 = None, None, None
        specialCase = None
        for yo in yosDatas : yosDatas[yo] = []

        I942Present = False
        for data in reactionDatas :
            if( data.I == 942 ) : I942Present = True
        if( I942Present ) :
            for i1 in range( len( reactionDatas ) ) :
                if( reactionDatas[i1].I == 4 ) :
                    del reactionDatas[i1]
                    break

        for data in reactionDatas :                         # yosDatas should only contain I = 1, 3, 4, 7, 9, 10, 13, 20, 21, 22, 941 and 942 type data.
            if( excludeAverageProductData and ( data.I in [ 10, 13 ] ) ) :
                if( not( ( C == 82 ) and ( data.yo == 9 ) ) ) : continue
            data.frame = (endlmisc.getNumberOfColumns_(data.I, '') * ('%s,' % xDataEnumsModule.Frame.lab))[:-1]
            if( data.S == 7 ) :                             # Special treatment for delayed neutron data.
                decayRate = data.getX1( )
                if( decayRate not in delayedNeutrons ) : delayedNeutrons[decayRate] = []
                delayedNeutrons[decayRate].append( data )
                continue
            if( 0 <= data.yo <= 19 ) :
                if( data.I in [ 80, 81, 84, 89, 90, 91, 92 ] ) : continue
                if( data.I == 0 ) :
                    I0 = data
                elif( data.I == 11 ) :
                    continue
                elif( data.I == 12 ) :
                    I12 = data
                elif( data.I == 20 ) :
                    I20 = data
                elif( ( data.yo == 0 ) and ( data.I == 13 ) ) :
                    continue
                else :
                    if( data.I not in [ 1, 3, 4, 7, 9, 10, 13, 20, 21, 22, 941, 942 ] ) :
                        print(data)
                        raise Exception( 'I = %d data not currently supported' % data.I )
                    if( ( data.S == 8 ) and ( data.yo == 1 ) and ( data.I == 4 ) and ( data.C == 12 ) ) :   # This should only be for endl n+d->2n+p
                        if( data.ZA != 1002 ) :
                            print(data)
                            raise Exception( 'S = 8, yo = 1, I = 4 only supported for ZA = 1002' )
                        yosDatas[data.yo + 10].append( data )
                        specialCase = 'S8'
                    else :
                        yosDatas[data.yo].append( data )
            else :
                print(data)
                raise Exception( 'Unsupported yo = %d' % data.yo )

        if( specialCase == 'S8' ) :                         # May need to add I = 10 data to yo = 11 and readjust yo = 1's I = 10 data.
            if( len( yosDatas[11] ) ) :
                I10Present = False
                for data in yosDatas[11] :
                    if( data.I == 10 ) : I10Present = True
                if( not I10Present ) :
                    yo1I10 = None
                    for data in yosDatas[1] :
                        if( data.I == 10 ) : yo1I10 = data
                    yo12I10 = None
                    for data in yosDatas[12] :
                        if( data.I == 10 ) : yo12I10 = data
                    if( ( yo1I10 is not None ) and ( yo12I10 is not None ) ) :
                        yosDatas[11].append( yo12I10 )
                        yo1I10.set( yo1I10 - yo12I10 )
        for yo in yosDatas :
            if( ( len( yosDatas[yo] ) == 1 ) and ( yosDatas[yo][0].I == 10 ) ) : yosDatas[yo] = []
        if( I0 is None ) : raise Exception( 'Missing cross section data for reaction C = %d S = %d.' % ( C, S ) )   # End of Phase 1

        axes = crossSectionAxes
        xgrid = I0.xArray()
        hasRegions = len(xgrid) > len(set(xgrid))
        if I0.C in (26, 38, 43, 47, 48) and len( self.findDatas( C=I0.C, I=0 ) ) > 1:
            yos = [dat.yo for dat in self.findDatas( C=I0.C, X1=I0.X1 )]
            if max(yos) == 16:  # alpha as residual: assume this is a breakup, and recast as a different reaction
                I0.C = {26:11, 38:44, 43:42, 47:41, 48:40}[ I0.C ]
                print( "    WARNING: recasting C=%d to C=%d" % ( C, I0.C ) )
        if( hasRegions ):
            regions1d = crossSectionModule.Regions1d( label = info.style, axes = axes )
            data = []
            x1 = -1
            for x2, y2 in I0.data :
                if( x1 == x2 ) :
                    regions1d.append( crossSectionModule.XYs1d( data = data, interpolation = getXYInterpolation( I0 ), axes = axes ) )
                    data = []
                data.append( [ x2, y2 ] )
                x1 = x2
            crossSection = crossSectionModule.XYs1d( data = data, interpolation = getXYInterpolation( I0 ), axes = axes )
            if( len( regions1d ) > 0 ) :
                regions1d.append( crossSection )
                crossSection = regions1d
            crossSection.label = info.style
        elif( I0.C == 74 ) :                # First value is 0 so make it first two points lin-lin like mcfgen does.
            crossSection = crossSectionModule.Regions1d( label = info.style, axes = axes )
            crossSection.append(crossSectionModule.XYs1d(data=I0.data[:2], interpolation=xDataEnumsModule.Interpolation.linlin, axes=axes))
            crossSection.append( crossSectionModule.XYs1d( data = I0.data[1:], interpolation = getXYInterpolation( I0 ), axes = axes ) )
        else :
            crossSection = crossSectionModule.XYs1d( data = I0.data, label = info.style, 
                    interpolation = getXYInterpolation( I0 ), axes = axes )
                
        if( C == 1 ) :
            yos = ( 0, )
        elif( isThermalNeutronScatteringLawReaction ) :
            residualZA, yos = targetZA, ( 1, targetZA )

            points = I0.data[:]
            while( points[-1][1] == 0 ) : del points[-1]
            priorEnergy = -1
            for energyXSec in points :
                if( energyXSec[0] == priorEnergy ) : energyXSec[0] += 1e-6 * energyXSec[0]
                priorEnergy = energyXSec[0]
            I0 = endlIClasses.endlI0(None, I0.yo, I0.C, I0.I, I0.S, I0.h, points, I0.bdflsFile)
            if( TNSL_EMax != I0.xMax( ) ) :
                print( '    WARNING: TNSL_EMax = %s != I0.xMax( ) = %s: %s' % (TNSL_EMax, I0.xMax( ), repr(I0)))
        elif( ( C == 12 ) and ( targetZA == 4009 ) and ( self.yi == 1 ) ) :  # Special case for Be_9(n,2n) which should really be Be_9(n,2n a)
            specialCase = 'Be_9(n, 2n)'
            yos = ( 1, 1, 2004 )
        elif( ( C == 13 ) and ( targetZA == 4010 ) and ( self.yi == 1 ) ) :  # Special case for Be10(n,3n)
            specialCase = 'Be_10(n, 3n)'
            yos = ( 1, 1, 1, 2004 )
        elif( ( C == 14 ) and ( targetZA == 4011 ) and ( self.yi == 1 ) ) :  # Special case for Be_11(n,4n)
            specialCase = 'Be_11(n, 4n)'
            yos = ( 1, 1, 1, 1, 2004 )
        else :
            yos = endl_C.endl_C_yoInfo(C)                 # Phase 2. Determine products (i.e., yos and residuals).
        for yo in yos :
            if( yo in [ 8, 9, 10 ] ) :
                if( C in [ 73, 74, 82 ] ) : continue
                raise Exception( 'Electron, positron or electron capture (i.e., yo = 8, 9 or 10) not supported: C = %d' % C )
        residuals = []
        if( len( yos ) == 0 ) :
            if( C != 1 ) : print('C = %d with no yos is currently not handled by toGNDS' % C)
            continue
        elif( yos[0] == -1 ) :                              # Unknown yos.
            if( C != 1 ) : print('C = %d has unknown products which is currently not handled by toGNDS' % C)
            continue
        elif( yos[0] == -2 ) :                              # yo = yi
            if( len( yos ) != 1 ) :
                print('product list = ', yos)
                raise Exception( 'For C = %d, endl_C.endl_C_yoInfo returned -2 plus other products, only the -2 is supported.' % C )
            yos = ( yosZAs[self.yi], )
            residuals = [ targetZA ]
        elif( yos[0] == -3 ) :                                          # fission
            residuals = [ 99120, 99120 ]
        elif( yos[0] == -4 ) :                                          # nX?
            pass                                                        # Leave residuals empty as residual is undefined.
        else :
            yiZA = endl2.yoToZA(self.yi)
            if( self.A == 0 ) : yiZA = 1000 * ( yiZA // 1000 )
            productsZA = -yiZA
            for yo in yos :
                if( yo != 7 ) : productsZA += yo
            if( self.A == 0 ) : productsZA = 1000 * ( productsZA // 1000 )
            residuals = [ targetZA - productsZA ]        # End of phase 2

        #FIXME do we still need the next 4 lines? 4008 isn't in residuals for Be9(n,2n), Be10(n,3n) or Be11(n,4n) due to specialCase logic above
        if( 4008 in residuals ) :                                       # Special case for Be_8 breakup into two alphas
            if( ( 2004 not in yos ) and ( yosDatas[6] != [] ) and ( yosDatas[16] == [] ) ) :
                yosDatas[16] = yosDatas[6]
                yosDatas[6] = []

        if( 99000 < targetZA < 99200 ) :
            residuals = [ targetZA ]
        elif(( C == 30 ) and (endl2.yoToZA(self.yi) in [1002, 1003]) and (targetZA in [1002, 1003])) :   # Special treatment for C = 30 as
            if( yosDatas[6] != [] ) : raise Exception( 'C = 30 as yo = 6 data' )                    # yos should be [ g, n ] and not [ g, n, a ].
            yos = ( 7, 1 )
            residuals = ( 2004, )

        if( verbose > 0 ) : print('   ', C, S, end='')
        if( isThermalNeutronScatteringLawReaction ) :
            if( 0.9e-12 < reactionDatas[0].getX1( ) < 1.1e-12 ) :
                TNSL_MT = 4
                from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentInelastic
                channelProcess = incoherentInelastic.Form.process
            elif( 1.9e-12 < reactionDatas[0].getX1( ) < 2.1e-12 ) :
                TNSL_MT = 2
                from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import coherentElastic
                channelProcess = coherentElastic.Form.process
                if( self.ZA == 1901 ) : channelProcess = "Inc" + channelProcess[1:]
            else :
                raise Exception("Encountered TNSL data with unexpected X1 = %f!" % reactionDatas[0].getX1( ))
            outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process=channelProcess)
            outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )

            firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 1 ), crossSection )
            addDistributionDataAndRemove( firstParticle, 1, yosDatas, crossSection )
            outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )

            for yo in yosDatas : yosDatas[yo] = []
        elif( ( yos[0] > 0 ) and ( C not in [ 71, 72, 73, 74, 82 ] ) ) :  # Reactions that are two-body, two-body + break-up or discrete N-body except C = 71 - 74.
            outputChannel = None
            IExtras = []
            if( ( C == 8 ) and ( self.yi == 9 ) ) : IExtras = [ 22 ]
            Is = getAllDistributionIs( yosDatas[ZAsYos[yos[0]]], IExtras = IExtras )
            if( ( Is == [ 1 ] ) or ( S == 1 ) ) :     # Two-body initial product state.
                if( verbose > 0 ) : print( ' 2Body', end = '' )
                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, yos[0] ),
                        crossSection )    # Assume first particle is in the ground state.
                addDistributionDataAndRemove( firstParticle, ZAsYos[yos[0]], yosDatas, crossSection )
                if crossSectionC8 is not None and C == 9:
                    Q = returnConstantQ( info.style, I0.getQ( ), crossSectionC8 )
                    crossSectionC8 = None
                else:
                    Q = returnConstantQ( info.style, I0.getQ( ), crossSection )
                if( len( yos ) == 1 ) :                                 # Simple two-body, yos[0] + Residual [ + gamma].
                    resYo = residuals[0]
                    if( resYo in ZAsYos ) : resYo = ZAsYos[resYo] + 10
                    if( S == 1 ) :
                        Q = returnConstantQ( info.style, I0.getQ( ) - I0.getX1( ), crossSection )
                        resParticle = toGNDSMiscModule.newGNDSParticle(info, toGNDSMiscModule.getTypeName(info, residuals[0]), crossSection)

                        addResidual = False
                        decayChannel = None
                        if( 7 in yos ) : print( ' 7 in yos' )

                        if( yosDatas[7] != [] ) :
                            decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                            decayChannel.Q.add( returnConstantQ( info.style, I0.getX1( ), crossSection ) )  # Assume residual returns to ground state, hence Q = getX1( )
                            gammaParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 7 ),
                                    crossSection )
                            Is = sorted( [ data.I for data in yosDatas[7] ] )
                            if( Is == [ 10, 13 ] ) :
                                print('  gamma has I = 10 and 13 but no others', end='')
                                yosDatas[7] = []
                            else :
                                addDistributionDataAndRemove( gammaParticle, 7, yosDatas, crossSection )
                            decayChannel.products.add( decayChannel.products.uniqueLabel( gammaParticle ) )
                            addResidual = True

                        if( ( residuals == [ 4008 ] ) and ( yosDatas[16] != [] ) ) :
                            if( decayChannel is None ) :
                                decayChannel =  outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                                decayChannel.Q.add( returnConstantQ( info.style, I0.getX1( ), crossSection ) )  # Assume residual returns to ground.
                            He4Particle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 2004 ),
                                    crossSection, multiplicity = 2 )
                            addDistributionDataAndRemove( He4Particle, 16, yosDatas, crossSection )
                            decayChannel.products.add( decayChannel.products.uniqueLabel( He4Particle ) )
                            addResidual = False

                        if( addResidual ) : decayChannel.products.add( decayChannel.products.uniqueLabel( resParticle ) )
                        level = I0.getX1( )
                        levelIndex = residualExcitationIndexLevels[residuals[0]][level]
                        if( levelIndex is None ) : level = None
                        secondParticle = toGNDSMiscModule.newGNDSParticle( info,
                                toGNDSMiscModule.getTypeName( info, residuals[0], level = level, levelIndex = levelIndex ),
                                crossSection, outputChannel = decayChannel )
                        if( Is == [ 1 ] ) : addRecoilDistributionAndRemove( secondParticle, firstParticle, resYo, yosDatas )

                        if levelIndex > 0 and not secondParticle.outputChannel:
                            particle = info.PoPs[ secondParticle.pid ]
                            groundState = particle.isotope.nuclides[0]
                            ZA = chemicalElementMiscPoPsModule.ZA( groundState )
                            multiplicity = secondParticle.multiplicity[0].copy()
                            secondParticle.addOutputChannel(outputChannelModule.OutputChannel(enumsModule.Genre.NBody))
                            secondParticle.outputChannel.Q.add( toGNDSMiscModule.returnConstantQ(
                                info.style, particle.nucleus.energy[0].value, multiplicity ) )
                            product = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getPoPsParticle( info, ZA, groundState.id ),
                                    crossSection, multiplicity = multiplicity )
                            product.distribution.add( unspecifiedModule.Form(info.style, xDataEnumsModule.Frame.lab) )
                            secondParticle.outputChannel.products.add( product )
                    else :
                        decayChannel, level, levelIndex = None, None, None
                        if( S == 0 ) :
                            level = I0.getELevel( )
                            levelIndex = None
                            if( C not in [ 8, 9, 10 ] ) :
                                if( ( yosDatas[7] != [] ) and ( residuals[0] != 7 ) ) :
                                    decayChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                                    decayChannel.Q.add( returnConstantQ( info.style, I0.getX1( ), crossSection ) )
                                    resParticle = toGNDSMiscModule.newGNDSParticle( info,
                                            toGNDSMiscModule.getTypeName( info, residuals[0] ), crossSection )       # Residual in ground state.
                                    addDistributionDataAndRemove( resParticle, resYo, yosDatas, crossSection )
                                    decayChannel.products.add( decayChannel.products.uniqueLabel( resParticle ) )
                                    gammaParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 7 ),
                                            crossSection )
                                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas, crossSection )
                                    decayChannel.products.add( decayChannel.products.uniqueLabel( gammaParticle ) )
                            if( ( level == 0.0 ) and ( levelIndex is None ) ) : level = None
                        if( residuals[0] <= 2004 ) : level, levelIndex = None, None
                        if( C==10 and I0.getELevel( ) != 0 ) :        # Elastic scattering off an isomer, use alias for reaction label
                            secondParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                        else:
                            secondParticle = toGNDSMiscModule.newGNDSParticle( info,
                                    toGNDSMiscModule.getTypeName( info, residuals[0], level = level, levelIndex = levelIndex ),
                                    crossSection, outputChannel = decayChannel )
                        addRecoilDistributionAndRemove( secondParticle, firstParticle, resYo, yosDatas )
                else :                                          # Two-body with residual breaking up.
                    if( verbose > 0 ) : print( ' w/breakup', end = '' )
                    if( self.A == 0 ) : raise Exception( 'For C = %d, breakup of product not supported for two-body reaction for natural target.' )
                    mYos = getMultiplicityYos( self, yos[1:], residuals, yosDatas )
                    residualZA = targetZA + endl2.yoToZA(self.yi) - yos[0]
                    level = 0.0
                    levelIndex = None
                    if( I0.S == 1 ) :
                        level = I0.getX1( )
                        levelIndex = residualExcitationIndexLevels[residualZA][level]
                    elif( I0.S == 8 ) : 
                        levelIndex, level = 1, I0.getX1( )
                    Q_MeV = endl2.reactionQByZAs([endl2.yoToZA(self.yi), targetZA], [yos[0], residualZA]) - level
                    residualName = toGNDSMiscModule.getTypeName( info, residualZA, level = level, levelIndex = levelIndex )
                    decayChannel = makeNBodyChannelFrom_mYos(outputChannelModule.OutputChannel, enumsModule.Genre.NBody, mYos, ZAsYos, yosDatas, 
                            I0.getQ( ) - Q_MeV, crossSection, specialCase=specialCase)
                    secondParticle = toGNDSMiscModule.newGNDSParticle( info, residualName, crossSection,
                            outputChannel = decayChannel )
                    Q = returnConstantQ( info.style, Q_MeV, crossSection )
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.twoBody)
                outputChannel.Q.add( Q )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            elif( Is == [ 22 ] ) :
                if( verbose > 0 ) : print( ' 2Body', end = '' )
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.twoBody)
                Q = returnConstantQ( info.style, 0.0, crossSection )
                outputChannel.Q.add( Q )

                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, yos[0] ), crossSection )
                I22 = yosDatas[9].pop( )
                points = []
                for energy, PofMu in I22.data : points.append( [ energy, list( reversed( [ [ 1 - mu, P ] for mu, P in PofMu ] ) ) ] )
                I22.setI( 1 )
                I1 = endlIClasses.endlI1(None, I22.yo, I22.C, 1, I22.S, I22.h, points, I22.bdflsFile)
                addDistributionDataAndRemove( firstParticle, 9, { 9 : [ I1 ] }, crossSection )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )

                secondParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            else :                                              # Mutli-particle breakup (i.e., not two body).
                if( verbose > 0 ) : print( ' NBody', end = '' )
                mYos = getMultiplicityYos( self, yos, residuals, yosDatas )
                gammaPresent = False
                for m, yo in mYos : gammaPresent = gammaPresent or ( yo == 7 )
                if( ( not gammaPresent ) and ( yosDatas[7] != [] ) ) : mYos.append( [ 1, 7 ] )

                if( specialCase == 'Be_9(n, 2n)' ) :                                            # Special case for endl Be_9(n, 2n)
                    if( yosDatas[6] == [] ) : mYos = [[2, 1], [2, 2004] ]
                elif( specialCase == 'Be_10(n, 3n)' ) :
                    if( yosDatas[6] == [] ) : mYos = [[3, 1], [2, 2004] ]
                elif( specialCase == 'Be_11(n, 4n)' ) :
                    if( yosDatas[6] == [] ) : mYos = [[4, 1], [2, 2004] ]

                if( C == 83 ) :
                    channelProcess = 'excitation'
                    mYos = mYos[:-1]
                outputChannel = makeNBodyChannelFrom_mYos(outputChannelModule.OutputChannel, enumsModule.Genre.NBody, mYos, ZAsYos, yosDatas, I0.getQ(), 
                        crossSection, specialCase=specialCase, process=channelProcess)
                if( C == 83 ) :
                    electron = outputChannel.products[0]
                    for data in reactionDatas :
                        if( ( data.yo == 9 ) and ( data.I == 10 ) ) :
                            averageProductEnergy = averageProductEnergyModule.XYs1d( label = info.style, axes = averageProductEnergyAxes, data = data.data )
                            electron.averageProductEnergy.add(averageProductEnergy)
                    particle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                    outputChannel.products.add( particle )

        else :  # Reactions that are not two-body, two-body + break-up or discrete N-body except C = 71 to 74.
            if( I0.C == 1 ) :                                 # (n,total)
                outputChannel = None
            elif( I0.C == 15 ) :                              # Fission
                if( verbose > 0 ) : print( ' Fission ', end = '' )

                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )
                particle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 1 ), crossSection )
                addDistributionDataAndRemove( particle, 1, yosDatas, crossSection )
                outputChannel.products.add( outputChannel.products.uniqueLabel( particle ) )

                decayedLinks = []
                totalDelayedNuBar = multiplicityModule.XYs1d( axes = multiplicityAxes, label = info.style )

                decayRates = sorted( delayedNeutrons.keys( ) )
                for decayRate in decayRates :
                    for delayedDatum in delayedNeutrons[decayRate] :
                        if( delayedDatum.I == 7 ) : totalDelayedNuBar += multiplicityModule.XYs1d( data = delayedDatum.data, axes = multiplicityAxes )
                decayRates.reverse( )

                for i1, decayRate in enumerate( decayRates ) :
                    gammasRemoved = []                          # Remove any delayed gamma data.
                    for index, delayedDatum in enumerate( delayedNeutrons[decayRate] ) :
                        if( delayedDatum.yo == 7 ) : gammasRemoved.append( index )
                    if( len( gammasRemoved ) > 0 ) : print( )
                    while( len( gammasRemoved ) > 0 ) :
                        gammaRemoved = delayedNeutrons[decayRate].pop( gammasRemoved[-1] )
                        print( '        Removing delayed gamma data', repr( gammaRemoved ) )
                        gammasRemoved.pop( )

                    delayedData = { 1 : delayedNeutrons[decayRate] }        # 1 is for neutron.
                    delayedParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 1 ), crossSection )

                    Is = sorted( [ data.I for data in delayedData[1] ] )
                    if( ( 4 not in Is ) and ( 3 not in Is ) ) : # If distribution data incomplete, use prompt neutron distribution data.
                        delayedNeutronData = delayedData[1]
                        delayedData[1] = []
                        for data in delayedNeutronData :
                            if( data.I == 7 ) : delayedData[1].append( data )
                        Is = sorted( [ data.I for data in delayedData[1] ] )

                    if( Is == [ 7, 10, 13 ] ) :
                        if( verbose ) : print('  delayed neutron has I = 7, 10 and 13 but no others', end = '')
                    else :
                        addDistributionDataAndRemove( delayedParticle, 1, delayedData, crossSection, promptNeutronParticle = particle )

                    product = delayedNeutronModule.Product( IDsPoPsModule.neutron, label = IDsPoPsModule.neutron )
                    product.multiplicity.add( delayedParticle.multiplicity[info.style] )
                    product.distribution.add( delayedParticle.distribution[info.style] )
                    delayedNeutron = delayedNeutronModule.DelayedNeutron( str( i1 ), product )
                    delayedNeutron.rate.add( rateModule.Double( info.style, decayRate, '1/s' ) )
                    outputChannel.fissionFragmentData.delayedNeutrons.add( delayedNeutron )
                    decayedLinks.append( sumsModule.Add( link = delayedNeutron.product.multiplicity ) )

                if( len( decayRates ) > 0 ) :
                    delayedNubar = sumsModule.MultiplicitySum( label = "delayed fission neutron multiplicity", ENDF_MT = 455 )
                    for summand in decayedLinks : delayedNubar.summands.append( summand )
                    delayedNubar.multiplicity.add( totalDelayedNuBar )
                    reactionSuite.sums.multiplicitySums.add( delayedNubar )

                    total = particle.multiplicity[info.style] + totalDelayedNuBar
                    total.label = info.style
                    totalNubar = sumsModule.MultiplicitySum( label = "total fission neutron multiplicity", ENDF_MT = 452 )
                    totalNubar.summands.append( sumsModule.Add( link = particle.multiplicity ) )
                    totalNubar.summands.append( sumsModule.Add( link = delayedNubar.multiplicity ) )
                    totalNubar.multiplicity.add( total )
                    reactionSuite.sums.multiplicitySums.add( totalNubar )

                if( yosDatas[7] != [] ) :
                    gammaParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 7 ), crossSection )
                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas, crossSection )
                    outputChannel.products.add( outputChannel.products.uniqueLabel( gammaParticle ) )

            elif( I0.C in [ 71, 72 ] ) :
                if( verbose > 0 ) : print( " 7[12] ", end = '' )
                for idx in range( len( yosDatas[7] ) - 1, -1, -1 ) :
                    if( ( yosDatas[7][idx].I == 4 ) or ( ( I0.C == 74 ) and ( yosDatas[7][idx].I == 10 ) ) ) : del yosDatas[7][idx]
                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 7 ), crossSection )
                doubleDifferentialCrossSectionForm = addDistributionDataAndRemove( firstParticle, 7, yosDatas, crossSection )
                secondParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                process = { 71 : 'coherent', 72 : 'incoherent' }[I0.C]
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.twoBody, process=process)
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            elif( I0.C == 73 ) :
                if( verbose > 0 ) : print( " 73", end = '' )
                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 9 ), crossSection )
                secondParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process="photo-electric")
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            elif( I0.C == 74 ) :
                if( verbose > 0 ) : print( " 74 ", end = '' )
                yosDatas[7] = []
                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 9 ), crossSection )
                secondParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 8 ), crossSection )
                thirdParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process='pair production: electron field')
                outputChannel.Q.add( returnConstantQ( info.style, -crossSection.domainMin, crossSection ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( thirdParticle ) )
            elif( I0.C == 82 ) :
                if( verbose > 0 ) : print( " 82", end = '' )
                outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody, process="bremsstrahlung")
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )

                product = toGNDSMiscModule.getTypeNameGamma( info, 0 )
                thirdParticle = toGNDSMiscModule.newGNDSParticle( info, product, crossSection )
                energyData = yosDatas[7].pop( )
                uncorrelatedIsotropicGivenEnergySpectrum( info, thirdParticle, energyData )
                outputChannel.products.add( outputChannel.products.uniqueLabel( thirdParticle ) )

                firstParticle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, 9 ), crossSection )
                for data in reactionDatas :
                    if( ( data.yo == 9 ) and ( data.I == 10 ) ) :
                        averageProductEnergy = averageProductEnergyModule.XYs1d( label = info.style, axes = averageProductEnergyAxes, data = data.data )
                        firstParticle.averageProductEnergy.add(averageProductEnergy)
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )

                secondParticle = toGNDSMiscModule.newGNDSParticle( info, info.targetID, crossSection )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            else :
                if( ( I0.C >= 50 ) and ( I0.C < 57 ) ) : 
                    if( verbose > 0 ) : print( " production", end = '' )
                    if( yos[0] != -4 ) : raise Exception(
                        'Does not appear to be production channel. Internal error. C = %d, yos = %s' %(C, repr(yos)) )
                    multiplicity = multiplicityModule.PartialProduction( info.style )
                    particle = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeName( info, yos[1] ),
                            crossSection, multiplicity = multiplicity )
                    addDistributionDataAndRemove( particle, ZAsYos[yos[1]], yosDatas, crossSection )
                    outputChannel = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
                    outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )
                    outputChannel.products.add( outputChannel.products.uniqueLabel( particle ) )
                else :
                    if( verbose > 0 ) : print()
                    print( " NONO:", I0 )
                    continue
        if( verbose > 0 ) : print()
        tmp = str( outputChannel )
        if( S == 2 ) :
            tmp += ' S2'
        elif( C == 9 ) :
            tmp += ' n+i'
        elif( C == 55 ) :
            tmp += ' discreteEnergy = %s' % I0.getX1( )
        if( tmp in channelList ) :
            endlmisc.printWarning('    WARNING: %s already in list for %d: C=%s, S=%s, X1=%s\n' % (tmp, self.ZA, C, S, I0.getX1()))
        else :
            channelList.append( tmp )

        if I12 is not None:
            Q = QModule.XYs1d( label = info.style, data = I12.data, axes = QAxes, interpolation = getXYInterpolation( I12 ) )
            outputChannel.Q.remove( info.style )    # replace constant approximation with energy-dependent version
            outputChannel.Q.add( Q )

        CCounts = len( self.findDatas( C = I0.C, I = 0 ) )
        if( isThermalNeutronScatteringLawReaction ) :
            MT = TNSL_MT
        elif( I0.C in ( 8, 9 ) ) :
            MT = 2
            if( self.yi == 9 ) :
                MT = 526
                outputChannel.process = 'large angle Coulomb scattering'
        elif( I0.C == 20 and len(self.findDatas(C=21))==0 ) :
            MT = 28
        elif( I0.C == 21 and len(self.findDatas(C=20))==0 ) :
            MT = 28
        elif( I0.C == 82 ) : 
            MT = 527
            outputChannel.process = 'bremsstrahlung'
        elif( I0.C == 83 ) : 
            MT = 528
        else :
            I0S = I0.S
            if outputChannel is not None and outputChannel.genre == enumsModule.Genre.twoBody:
                I0S = 1
            MT = ENDLCS_To_ENDFMT.getMTFromCS( I0.C, I0S, CCounts = CCounts )
        if( I20 is not None ) : endlmisc.printWarning('    WARNING: I20 URR data not currently supported.')

        if( C == 1 ) :
            reaction = sumsModule.CrossSectionSum( label = "total", ENDF_MT = MT )
            for _reaction in reactionSuite.reactions : reaction.summands.append( sumsModule.Add( link = _reaction.crossSection ) )
            iChannel += 1
            reaction.Q.add( returnConstantQ( info.style, I0.getQ( ), crossSection ) )
        else :
            fissionGenre = enumsModule.FissionGenre.none
            if C == 15:
                fissionGenre = enumsModule.FissionGenre.total
            reaction = reactionModule.Reaction( str(MT), outputChannel.genre, ENDF_MT=MT, fissionGenre=fissionGenre )
            endf_endlModule.setReactionsOutputChannelFromOutputChannel( info, reaction, outputChannel )
            if( S == 2 ) : reaction.outputChannel.process = 'ENDL:S2'

        reaction.crossSection.add( crossSection )
        reactionsByCS[(C, S)] = reaction
        if( doubleDifferentialCrossSectionForm is not None ) : reaction.doubleDifferentialCrossSection.add( doubleDifferentialCrossSectionForm )

        # keep track of TNSL reactions in order to exclude them from non-elastic sum if C=55 data present
        reaction._isThermalScatteringReaction = isThermalNeutronScatteringLawReaction

        if( ( self.yi != 9 ) and ( C in ( 8, 9 ) ) ) :
            if C == 8 and self.findDatas(C=9):
                crossSectionC8 = crossSection
                continue                            # If both C=8 and C=9 are present, follow C==9 logic.

            # convert cross section / distribution into nuclearPlusInterference instance:

            identicalParticles = target is projectile
            firstProduct = reaction.outputChannel[0]
            if( C==9 ):
                cutoffs = [tmp.trim().domainMax for tmp in firstProduct.distribution.evaluated.angularSubform]
                muCutoff = max( cutoffs )
                effXsc = reaction.crossSection.evaluated.copy()
                effDist = firstProduct.distribution.evaluated.angularSubform
                for idx in range(len(effDist)):
                    if identicalParticles:  # ENDL covers mu=-cutoff .. cutoff, but only need to store 0 .. cutoff
                        tmpDist = effDist[idx].domainSlice(0, muCutoff)
                        tmpDist = tmpDist.normalize()
                    else:
                        tmpDist = effDist[idx].domainSlice(-1, muCutoff)
                    if ( tmpDist[-1][1] == 0 and
                            tmpDist[-2][0] / muCutoff > 0.9999 ):   # Get rid of 'blurred edge' at cutoff angle
                        rawDat = tmpDist.copyDataToXYs()
                        rawDat[-1][1] = rawDat[-2][1]
                        rawDat.pop(-2)
                        tmpDist.setData( rawDat )
                    effDist[idx] = tmpDist

                effXsc.label = None # label unnecessary since it will reside inside CoulombPlusNuclearElastic
                nuclearPlusInterference = nuclearPlusInterferenceModule.NuclearPlusInterference( muCutoff=muCutoff,
                        crossSection=nuclearPlusInterferenceModule.CrossSection( effXsc),
                        distribution=nuclearPlusInterferenceModule.Distribution( effDist)
                        )
                Rutherford = RutherfordScatteringModule.RutherfordScattering()
            else:   # Only have C=8 data. Only give the Rutherford scattering term
                nuclearPlusInterference = None
                Rutherford = RutherfordScatteringModule.RutherfordScattering(
                        reaction.crossSection.domainMin, reaction.crossSection.domainMax,
                        reaction.crossSection.domainUnit )

            CPElastic = CoulombPlusNuclearElasticModule.Form( projectileID, info.style, identicalParticles = identicalParticles,
                    RutherfordScattering = Rutherford, nuclearPlusInterference = nuclearPlusInterference )
            reaction.doubleDifferentialCrossSection.add( CPElastic )
            reaction.crossSection.remove( info.style )
            reaction.crossSection.add( crossSectionModule.CoulombPlusNuclearElastic(
                    link = reaction.doubleDifferentialCrossSection.evaluated, label = info.style, relative = True ) )
            firstProduct.distribution.remove( info.style )
            firstProduct.distribution.add( referenceModule.CoulombPlusNuclearElastic(
                    link = reaction.doubleDifferentialCrossSection.evaluated, label = info.style, relative = True ) )

            secondProduct = reaction.outputChannel[1]
            secondProduct.distribution.evaluated.angularSubform.link = firstProduct.distribution.evaluated

        if( C == 1 ) :
            reactionSuite.sums.crossSectionSums.add( reaction )
        elif( C == 55 ) :
            # wait until all other reactions are done to handle C55
            c55s.append( reaction )
        else :
            reaction.updateLabel()
            reactionSuite.reactions.add( reaction )
        if( excludeAverageProductData ) :
            for yo in yosDatas :
                datas = yosDatas[yo]
                IsToDelete = []
                for i1, data in enumerate( datas ) :
                    if( data.I in [ 10, 13 ] ) : IsToDelete.insert( 0, i1 )
                for i1 in IsToDelete : del datas[i1]
        for yo in yosDatas :
            if( len( yosDatas[yo] ) ) :
                message = [repr(yosData) for yosData in yosDatas[yo]]
                endlmisc.printWarning('    WARNING: unused data for ZA = %d, yo = %d\n        %s' % (self.ZA, yo, '\n        '.join(message)))

    if c55s:
        # C=55 stores (total xsc)*multiplicity, but we want multiplicity for nonelastic only.
        # sigma_total * mult_total = sigma_nonelastic * M, so M = (C55 data) / sigma_nonelastic.
        # Also need to store non-elastic as a cross section sum

        nonElastic = sumsModule.CrossSectionSum( label = "nonelastic", ENDF_MT = 3 )
        summands = [ sumsModule.Add( link = r.crossSection ) for r in reactionSuite.reactions
                if r is not reactionSuite.getReaction('elastic') and not r._isThermalScatteringReaction ]
        for summand in summands : nonElastic.summands.append( summand )

        # calculate nonElastic cross section by addition  (total - elastic sometimes has numeric issues)
        nonElasticXSec = summands[0].link.toPointwise_withLinearXYs(upperEps=1e-8)
        for summand in summands[1:]:
            addend = summand.link.evaluated
            if not nonElasticXSec.areDomainsMutual( addend ):
                nonElasticXSec, addend = nonElasticXSec.mutualify( -1e-7,0,0, addend, -1e-7,0,0 )
            nonElasticXSec += addend

        nonElasticXSec.label = nonElasticXSec.style = info.style
        nonElastic.crossSection.add( nonElasticXSec )
        nonElastic.Q.add( returnConstantQ( info.style, 0, crossSection ) )
        reactionSuite.sums.crossSectionSums.add( nonElastic )

        gammaProduction = orphanProductModule.OrphanProduct('Orphan product 0', enumsModule.Genre.NBody, ENDF_MT = 3)
        gammaProduction.crossSection.add( crossSectionModule.Reference( link = nonElastic.crossSection, label = info.style ) )
        gammaProduction.outputChannel.Q.add( c55s[0].outputChannel.Q.evaluated )
        for C55 in c55s:
            orphanGamma = C55.outputChannel.products[0]
            orphanGamma.multiplicity.remove( info.style )
            if( False ) :           # Old way
                if not nonElasticXSec.areDomainsMutual( C55.crossSection.evaluated ):
                    num, den = C55.crossSection.evaluated.mutualify(-1e-8,0,0, nonElasticXSec, 0,0,0)
                    newMultiplicity = num / den
                else:
                    newMultiplicity = C55.crossSection.evaluated / nonElasticXSec
                newMultiplicity = multiplicityModule.XYs1d( data = newMultiplicity.thin( 0.001 ), label = info.style, axes = multiplicityAxes )
            else :                  # New way.
                newMultiplicity = []
                for energy, multiplicity in C55.crossSection.evaluated :
                    crossSection = nonElasticXSec.evaluate( energy )
                    if( crossSection != 0.0 ) : multiplicity /= crossSection
                    newMultiplicity.append( [ energy, multiplicity ] )
                newMultiplicity = multiplicityModule.XYs1d( data = newMultiplicity, label = info.style, axes = multiplicityAxes )
            orphanGamma.multiplicity.add( newMultiplicity )
            gammaProduction.outputChannel.products.add( gammaProduction.outputChannel.products.uniqueLabel(orphanGamma) )

        reactionSuite.orphanProducts.add( gammaProduction )

    # add a few other summed cross sections if applicable:
    for (MT, label, MT_range) in (
        (4, '(z,n)', range(51,92)),
        (103, '(z,p)', range(600,650)),
        (104, '(z,d)', range(650,700)),
        (105, '(z,t)', range(700,750)),
        (106, '(z,h)', range(750,800)),
        (107, '(z,a)', range(800,850)),
    ):
        reacs = [reac for reac in reactionSuite.reactions if reac.ENDF_MT in MT_range]
        if len(reacs) > 1:

            newSum = sumsModule.CrossSectionSum(label = label, ENDF_MT = MT)
            summands = [sumsModule.Add(link = r.crossSection) for r in reacs]
            for summand in summands : newSum.summands.append( summand )

            newSumXsc = summands[0].link.toPointwise_withLinearXYs(upperEps=1e-8)
            for summand in summands[1:]:
                addend = summand.link.evaluated
                if not newSumXsc.areDomainsMutual( addend ):
                    newSumXsc, addend = newSumXsc.mutualify( -1e-7,0,0, addend, -1e-7,0,0 )
                newSumXsc += addend

            newSumXsc.label = newSumXsc.style = info.style
            newSum.crossSection.add( newSumXsc )
            newQ = max([reac.getQ('MeV', final=True) for reac in reacs])
            newSum.Q.add( returnConstantQ( info.style, newQ, newSumXsc ) )
            reactionSuite.sums.crossSectionSums.add( newSum )

    covarianceFileCS = covarianceDict(self.source)   #  Dictionary covariance filenames from the user-input directory, key: (C,S)

    if len(covarianceFileCS) > 0:
        covarianceLabelsUnits = {0 : ( 'matrix_elements', '' ),
                                 1 : ( 'column_energy_bounds', 'MeV' ),
                                 2 : ( 'row_energy_bounds', 'MeV' ) }

        sectionList = []
        for (C,S) in covarianceFileCS.keys():
            if (C,S) in reactionsByCS:
                reaction = reactionsByCS[(C,S)]
            else:
                # likely means summed covariance was given for S=0 while reactions are listed as S=1.
                # FIXME skipping for now
                print("  WARNING: no matching reaction found for covariance with C=%d and S=%d" % (C,S))
                continue

            for thisCovarianceEntry in covarianceFileCS[C, S]:

                thisCovariancePath = os.path.join(self.source, thisCovarianceEntry)
                covarianceDatas = parse_endl_covariance(thisCovariancePath)
                if not covarianceDatas:
                    continue    # ignore empty covariance files

                covar_label = reaction.label

                I = int(thisCovarianceEntry.split('i')[1].split('s')[0])
                if( I == 0 ):
                    # cross section covariance
                    rowdat = covarianceSectionModule.RowData(root = '$reactions', link = reaction.crossSection.evaluated,
                            ENDF_MFMT=f'33,{reaction.ENDF_MT}') 
                elif( I == 4 and C == 15 ):
                    # Prompt fission neutron spectrum covariance
                    covar_label += ' [spectrum] energy range '
                    rowdat = covarianceSectionModule.RowData(root = '$reactions',
                            link = reaction.outputChannel.products[0].distribution.evaluated.energySubform,
                            ENDF_MFMT=f'35,{reaction.ENDF_MT}')
                elif( I==7 ):
                    # nubar covariance
                    covar_label += ' neutron multiplicity [nubar]'
                    rowdat = covarianceSectionModule.RowData(root = '$reactions',
                            link = reaction.outputChannel.products[0].multiplicity.evaluated,
                            ENDF_MFMT=f'31,{reaction.ENDF_MT}')
                else:
                    raise NotImplementedError('Covariances for I=%d, C=%d, S=%d' % (I,C,S))
                
                # trying to link to the covariance inside the mean-value protare
                centralValue = rowdat.link 
                

                def createCovarianceMatrix(covarianceEntry, label):
                    # helper method
                    axes = axesModule.Axes(3, labelsUnits = covarianceLabelsUnits)
                    energyBounds = covarianceEntry[0]
                    rawMatrix = numpy.asarray(covarianceEntry[1])
                    Type=covarianceEntry[2]

                    if numpy.count_nonzero(rawMatrix - numpy.diag(numpy.diagonal(rawMatrix))):
                        lowerDiagonal = rawMatrix[numpy.tril_indices(len(rawMatrix))]
                        matrix = arrayModule.Full( rawMatrix.shape, lowerDiagonal, symmetry=arrayModule.Symmetry.lower)
                    else:   # matrix is diagonal
                        matrix = arrayModule.Diagonal( rawMatrix.shape, rawMatrix.diagonal() )
                    axes[2] = axesModule.Grid( axes[2].label, axes[2].index, axes[2].unit,
                            style = xDataEnumsModule.GridStyle.boundaries, values = valuesModule.Values( energyBounds ) )
                    axes[1] = axesModule.Grid( axes[1].label, axes[1].index, axes[1].unit,
                            style = xDataEnumsModule.GridStyle.boundaries, values = linkModule.Link( link = axes[2].values, relative = True ) )
                    covmatrix = covarianceMatrixModule.CovarianceMatrix(label = label, type = Type,
                            matrix = griddedModule.Gridded2d( axes, matrix ))
                    return covmatrix

                if (I == 4):
                    for index, covarianceEntry in enumerate(covarianceDatas):
                        label = covar_label + str(index)
                        covarianceSection = covarianceSectionModule.CovarianceSection(label = label, rowData = rowdat.copy())
                        domain = covarianceEntry[3]
                        covarianceSection.rowData.slices.add(
                            covarianceSectionModule.Slice(2, domainUnit="MeV", domainMin=domain[0], domainMax=domain[1])
                        )
                        covmatrix = createCovarianceMatrix(covarianceEntry, info.style)
                        covarianceSection.add(covmatrix)
                        sectionList.append(covarianceSection)
                else:
                    covarianceSection=covarianceSectionModule.CovarianceSection(label = covar_label, rowData = rowdat)
                    if (len(covarianceDatas) > 1):
                        covarsThisSection = []
                        for index, covarianceEntry in enumerate(covarianceDatas):
                            covmatrix = createCovarianceMatrix(covarianceEntry, str(index))
                            covarsThisSection.append(covmatrix)
                        covarianceForm = covarianceMixedModule.MixedForm(label = info.style, components=covarsThisSection)
                        covarianceSection.add(covarianceForm)

                    elif (len(covarianceDatas) == 1):
                        covmatrix = createCovarianceMatrix(covarianceDatas[0], info.style)
                        covarianceSection.add(covmatrix)
                    
                    sectionList.append(covarianceSection)

                covarianceLink = uncertaintiesModule.Covariance(link=covarianceSection[info.style], root="$covariances")  # Change cov[] to whatever is the equivalent in here
                centralValue.uncertainty = uncertaintiesModule.Uncertainty( functional=covarianceLink )

        if len(sectionList) > 0:
            covarianceSuite = covarianceSuiteModule.CovarianceSuite( projectileID, info.targetID, evaluation, interaction, formatVersion = formatVersion )
            
            for cov in sectionList:
                covarianceSuite.covarianceSections.add(cov)

    evaluatedStyle.projectileEnergyDomain.min = min( [ reaction.crossSection.domainMin for reaction in reactionSuite.reactions ] )
    if( projectile.id == IDsPoPsModule.neutron ) :          # Ignore case where the last point's y-value is 0.
        maxs = []
        for reaction in reactionSuite.reactions :
            if( reaction.crossSection[0][-1][1] != 0 ) : maxs.append( reaction.crossSection.domainMax )
        evaluatedStyle.projectileEnergyDomain.max = min( maxs )
    else :
        evaluatedStyle.projectileEnergyDomain.max= min( [ reaction.crossSection.domainMax for reaction in reactionSuite.reactions ] )

    if covarianceSuite is not None:
        covarianceSuite.styles.add(evaluatedStyle.copy())   # don't include documentation in covariances
        reactionSuite._loadedCovariances = [covarianceSuite]

    documentation = self.getDocumentation( )
    if( documentation is None ) : documentation = "From ENDL: no documentation specified."
    evaluatedStyle.documentation.body.body = documentation

    for reaction in reactionSuite.reactions:
        frame = xDataEnumsModule.Frame.lab
        if reaction.outputChannel.genre == enumsModule.Genre.twoBody and self.yi != 7:
            frame = xDataEnumsModule.Frame.centerOfMass
        toGNDSMiscModule.addUnspecifiedDistributions(info, reaction.outputChannel, frame)
    for production in reactionSuite.productions:
        toGNDSMiscModule.addUnspecifiedDistributions( info, production.outputChannel, xDataEnumsModule.Frame.lab)

    for chemicalElement in reactionSuite.PoPs.chemicalElements :
        for isotope in chemicalElement :
            nuclide = isotope[0]
            if( len( nuclide.mass ) > 0 ) : continue
            ZA = chemicalElementMiscPoPsModule.ZA( nuclide )
            mass = self.bdflsFile.mass( ZA )
            mass = massPoPsModule.Double( info.PoPsLabel, mass, quantityPoPsModule.stringToPhysicalUnit( 'amu' ) )
            nuclide.mass.add( mass )
    
    return (reactionSuite, covarianceSuite)
