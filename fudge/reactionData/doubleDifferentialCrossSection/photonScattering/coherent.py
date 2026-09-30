# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon coherent double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule

from fudge.core.math import fudgemath as fudgemathModule
from ... import crossSection as crossSectionModule
from .. import base as baseModule


class XYs1d( baseModule.XYs1d ) :

    pass

class Regions1d( baseModule.Regions1d ) :

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class CoherentFunctionBase( ancestryModule.AncestryIO ) :

    def __init__( self, data ) :

        ancestryModule.AncestryIO.__init__( self )
        if not isinstance(data, (XYs1d, Regions1d)):
            raise TypeError( "Invalid data type: got %s." % type(data))

        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ):
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def check( self, info ) :

        return( [] )

    def toXML_strList( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )

        data = element[0]
        if( data.tag == baseModule.XYs1d.moniker ) :
            data = XYs1d.parseNodeUsingClass(data, xPath, linkData, **kwargs)
        elif( data.tag == baseModule.Regions1d.moniker ) :
            data = Regions1d.parseNodeUsingClass(data, xPath, linkData, **kwargs)
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( data.tag, cls.tag ) )

        instance = cls(data)

        xPath.pop( )

        return instance

class FormFactor( CoherentFunctionBase ) :

    moniker = 'formFactor'
    ENDFMT = 502

class ImaginaryAnomalousFactor( CoherentFunctionBase ) :

    moniker = 'imaginaryAnomalousFactor'
    ENDFMT = 505

class RealAnomalousFactor( CoherentFunctionBase ) :

    moniker = 'realAnomalousFactor'
    ENDFMT = 506

class Form( baseModule.Form ) :

    moniker = 'coherentPhotonScattering'
    keyName = 'label'

    subformAttributes = ( 'formFactor', 'anomalousScatteringFactor_realPart', 'anomalousScatteringFactor_imaginaryPart' )

    def __init__( self, pid, label, productFrame, _formFactor, realPart, imaginaryPart ) :

        if not isinstance( _formFactor, FormFactor ):
            raise Exception( 'Instance is class "%s" and not formFactor' % FormFactor.__class__ )

        if realPart is not None:
            if not isinstance( realPart, RealAnomalousFactor ):
                raise Exception( 'Instance is class "%s" and not RealAnomalousFactor' % realPart.__class__ )

        if imaginaryPart is not None:
            if not isinstance( imaginaryPart, ImaginaryAnomalousFactor ):
                raise Exception( 'Instance is class "%s" and not ImaginaryAnomalousFactor' % imaginaryPart.__class__ )

        baseModule.Form.__init__( self, pid, label, productFrame, ( _formFactor, realPart, imaginaryPart ) )

    def check( self, info ) :

        return []

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        multiplicity = kwargs['multiplicity']
        aveEnergy = [ [ multiplicity.domainMin, multiplicity.domainMin ], [ multiplicity.domainMax, multiplicity.domainMax ] ]
        aveMomentum = [ [ multiplicity.domainMin, 0 ], [ multiplicity.domainMax, 0 ] ]

        return( [ aveEnergy ], [ aveMomentum ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule
        from fudge.processing import group as groupModule

        verbosity = tempInfo['verbosity']
        if verbosity > 2: print('%sGrouping %s' % (indent, self.moniker))

        anomalousScatteringFactor = XYs1dModule.XYs1d( [ [ 0.0, 0.0 ], [ 21., 0.0 ] ] )
        anomalousScatteringFactor_realPart = anomalousScatteringFactor
        if self.anomalousScatteringFactor_realPart is not None:
            anomalousScatteringFactor_realPart = self.anomalousScatteringFactor_realPart.data
        anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor
        if self.anomalousScatteringFactor_imaginaryPart is not None:
            anomalousScatteringFactor_imaginaryPart = self.anomalousScatteringFactor_imaginaryPart.data

        TM_1, TM_E = transferMatricesModule.wholeAtomScattering( style, tempInfo, self.productFrame, self.formFactor.data,
                realAnomalousFactor = anomalousScatteringFactor_realPart, imaginaryAnomalousFactor = anomalousScatteringFactor_imaginaryPart,
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
        return groupModule.TMs2Form( style, tempInfo, TM_1, TM_E )

    def crossSection(self, tolerance=1e-3, interpolation=xDataEnumsModule.Interpolation.linlin):
        """
        Integrates the double differential cross section to get the cross section :math:`\\sigma(E)` where :math:`E` is the projectile's energy.
        """

        class Parameters :

            pass

        class Tester :

            def __init__( self, parameters, relativeTolerance, absoluteTolerance, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) :

                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart
                self.anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart

            def evaluateAtX( self, energy ) :

                return( integratedCrossSection( energy, self.parameters, self.anomalousScatteringFactor_realPart, self.anomalousScatteringFactor_imaginaryPart ) )

        def integrateSub( _n, _a, logX, x1, _y1, x2, _y2 ) :

            epsilon = _a + _n + 1
            if abs( epsilon ) < 1e-3:
                epsilon_logX = epsilon * logX
                integral = _y1 * math.pow( x1, _n + 1 ) * logX * ( 1 + 0.5 * epsilon_logX * ( 1 + epsilon_logX / 3.0 * ( 1 + 0.25 * epsilon_logX ) ) )
            else :
                integral = ( _y2 * math.pow( x2, _n + 1 ) - _y1 * math.pow( x1, _n + 1 ) ) / epsilon
            return integral
            

        def integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) :

            parameters.anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.evaluate( energy )
            if parameters.anomalousScatteringFactor_realPart == None: parameters.anomalousScatteringFactor_realPart = 0

            parameters.anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.evaluate( energy )
            if parameters.anomalousScatteringFactor_imaginaryPart == None: parameters.anomalousScatteringFactor_imaginaryPart = 0.0

            formFactorLower = parameters.formFactor[0]

            integralPrime = 0.0
            integralUpper = 0.0
            if energy <= formFactorLower.domainMax:
                muMin = -1.0
                integralLower = ( formFactorLower[0][1] + parameters.anomalousScatteringFactor_realPart )**2 + parameters.anomalousScatteringFactor_imaginaryPart**2
                integralLower *= 8.0 / 3.0
            else:
                ratio2 = 2.0 * ( formFactorLower.domainMax / energy )**2
                muMin = 1.0 - ratio2
                integralPrime = 8.0 / 3.0 * ( parameters.anomalousScatteringFactor_realPart**2 + parameters.anomalousScatteringFactor_imaginaryPart**2 )
                integralLower = ratio2 * ( 4 + muMin * ( 1 + muMin ) ) / 3.0 * \
                        ( formFactorLower[0][1] * ( formFactorLower[0][1] + 2 * parameters.anomalousScatteringFactor_realPart ) )

            i1 = 0
            if muMin > -1.0:
                formFactorLower = parameters.formFactor[1]
                for i1, xy in enumerate( formFactorLower ) :
                    _x2, y2 = xy
                    x2 = _x2
                    if( i1 > 0 ) :
                        xRatio = x2 / x1
                        logX = math.log( xRatio )
                        logY = math.log( y2 / y1 )
                        _a = logY / logX

                        if( x2 > energy ) :
                            x2 = energy
                            y2 = y1 * math.pow( energy / x1, _a )
                        x1_e = x1 / energy
                        x2_e = x2 / energy

                        integralSub  = 0.5 * integrateSub( 1, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub -=       integrateSub( 3, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub +=       integrateSub( 5, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub *= 2.0 * parameters.anomalousScatteringFactor_realPart

                        y1_2 = y1 * y1
                        y2_2 = y2 * y2
                        _2_a = 2 * _a
                        integralSub2  = 0.5 * integrateSub( 1, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )
                        integralSub2 -=       integrateSub( 3, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )
                        integralSub2 +=       integrateSub( 5, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )

                        integralUpper += integralSub + integralSub2
                    x1 = x2
                    y1 = y2
                    if _x2 >= energy: break
                integralUpper *= 16

            return (integralPrime + integralLower + integralUpper) * parameters.factor_crossSection

        parameters = Parameters( )

        product = self.ancestor.ancestor

        tolerance = max( 1e-6, min( 0.1, tolerance ) )
        parameters.tolerance = 1e-4 * tolerance

        factor_crossSection = PQUModule.PQU( 1.0, 'e**2 / ( 4. * pi * eps0 * me * c**2 )' )**2      # r_e^2 where r_e is the classical electron radius.
        factor_crossSection = math.pi * factor_crossSection.getValueAs( 'b' )
        parameters.factor_crossSection = factor_crossSection

        formFactor = self.formFactor.data
        if isinstance(formFactor, Regions1d):
            _formFactor = formFactor.copy()
        elif isinstance(formFactor, XYs1d):
            assert formFactor[0][1] == formFactor[1][1], "Non-constant form factor at low wave number!"
            for idx in range(len(formFactor)-1):
                if formFactor[idx+1][1] != formFactor[0][1]:
                    break
            assert idx < len(formFactor) // 10, "Form factor should decrease for longer wave number!"
            cutoff = formFactor[idx][0]

            _formFactor = Regions1d(axes=formFactor.axes)
            _formFactor.append(formFactor.domainSlice(None, cutoff))
            _formFactor.append(formFactor.domainSlice(cutoff, None))
        else :
            raise Exception( 'Form factor must be a Regions1d or XYs1d instance' )
        _formFactor.axes[1].unit = product.domainUnit
        factor_E2x = 1.0 / PQUModule.PQU( 1.0, '%s / hplanck / c' % product.domainUnit ).getValueAs( formFactor.axes[1].unit )
        for region in _formFactor : region.scaleOffsetXAndY( 0.0, factor_E2x, 0.0, 1.0, True )
        parameters.formFactor = _formFactor

        domainMin = product.domainMin
        domainMax = product.domainMax
        anomalousScatteringFactor_realPart = self.anomalousScatteringFactor_realPart
        if anomalousScatteringFactor_realPart is None:
            anomalousScatteringFactor_realPart = XYs1d( [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ], axes = _formFactor.axes.copy( ) )
        else :
            anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.data

        anomalousScatteringFactor_imaginaryPart = self.anomalousScatteringFactor_imaginaryPart
        if anomalousScatteringFactor_imaginaryPart is None:
            anomalousScatteringFactor_imaginaryPart = XYs1d( [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ], axes = _formFactor.axes.copy( ) )
        else :
            anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.data

        transitionEnergy = PQUModule.PQU( 1e3, 'eV' ).getValueAs( product.domainUnit )
        transitionEnergy2 = PQUModule.PQU( 1e4, 'eV' ).getValueAs( product.domainUnit )
        energy = domainMin
        initialEnergies = []
        # Fine grid below 1 keV (1% steps)
        while energy < transitionEnergy:
            initialEnergies.append( energy )
            energy *= 1.01
        # Medium grid from 1 keV to 10 keV (5% steps for better sampling of absorption edges)
        energy = transitionEnergy
        while energy < transitionEnergy2:
            initialEnergies.append( energy )
            energy *= 1.05
        # Coarse grid above 10 keV (decade steps)
        energy = transitionEnergy2
        while energy < domainMax:
            initialEnergies.append( energy )
            energy *= 10
        initialEnergies.append( domainMax )

        energy_crossSection = []
        for energy in initialEnergies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) ] )

        _tester = Tester( parameters, tolerance, tolerance * energy_crossSection[-1][1], 
                anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 10, interpolation = interpolation )
        return crossSectionModule.XYs1d( data = energy_crossSection,
                                         axes = crossSectionModule.defaultAxes( product.domainUnit ),
                                         interpolation = interpolation )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )

        subForms = { FormFactor.moniker : None, RealAnomalousFactor.moniker : None, ImaginaryAnomalousFactor.moniker : None }
        for child in element :
            for _class in ( FormFactor, RealAnomalousFactor, ImaginaryAnomalousFactor ) :
                if( child.tag == _class.moniker ) : break

            if( child.tag != _class.moniker ) : raise TypeError( "Invalid element '%s' encountered '%s'" % ( child.tag, _class.moniker ) )
            subForms[_class.moniker] = _class.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        _form = cls( element.get( 'pid' ), element.get( 'label' ), element.get( 'productFrame' ),
                        subForms[FormFactor.moniker], subForms[RealAnomalousFactor.moniker],
                        subForms[ImaginaryAnomalousFactor.moniker] )

        xPath.pop( )

        return( _form )
