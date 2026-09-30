# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon incoherent double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule
from numericalFunctions import integration as nf_integrationModule

from LUPY import ancestry as ancestryModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import enums as xDataEnumsModule

from fudge.core.math import fudgemath as fudgemathModule
from ... import crossSection as crossSectionModule
from .. import base as baseModule

class XYs1d(baseModule.XYs1d):

    pass

class Regions1d(baseModule.Regions1d):

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class ScatteringFactor(ancestryModule.AncestryIO):

    moniker = 'scatteringFactor'
    ENDFMT = 504

    def __init__(self, data):

        ancestryModule.AncestryIO.__init__(self)

        if not isinstance(data, (XYs1d, Regions1d)):
            raise TypeError( "Invalid data type: got %s." % type(data))

        self.__data = data
        self.__data.setAncestor(self)

    @property
    def data(self):

        return self.__data

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.__data.convertUnits(unitMap)

    def check(self, info):

        return []

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)

        if formatVersion == GNDS_formatVersionModule.version_1_10:
            return self.__data.toXML_strList(indent = indent, **kwargs)

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        xmlStringList += self.__data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param node:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *node*.
        """

        xPath.append( node.tag )

        data = node[0]
        if( data.tag == baseModule.XYs1d.moniker ) :
            data = XYs1d.parseNodeUsingClass(data, xPath, linkData, **kwargs)
        elif( data.tag == baseModule.Regions1d.moniker ) :
            data = Regions1d.parseNodeUsingClass(data, xPath, linkData, **kwargs)
        else :
            raise TypeError('Invalid data "%s" for "%s"' % (data.tag, cls.tag))

        scatteringFactor = cls(data)

        xPath.pop()

        return scatteringFactor

class Form( baseModule.Form ) :

    moniker = 'incoherentPhotonScattering'
    keyName = 'label'

    subformAttributes = ( 'scatteringFactor', )

    def __init__( self, pid, label, productFrame, scatteringFactor ) :

        if not isinstance(scatteringFactor, (ScatteringFactor)):
            raise Exception('Instance is class "%s" and not scatteringFactor' % scatteringFactor.__class__)

        baseModule.Form.__init__( self, pid, label, productFrame, ( scatteringFactor, ) )

    def check(self, info):
        return []

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        class Parameters :

            pass

        def integrand( mu, parameters ) :

            _x = min( parameters.energy_hc * math.sqrt( 0.5 * ( 1.0 - mu ) ), parameters.scatteringFactorDomainMax )
            kPrime_k = 1.0 / ( 1 + parameters.k * ( 1.0 - mu ) )
            integral = kPrime_k * kPrime_k * ( 1 + mu * mu + kPrime_k * parameters.k**2 * ( 1 - mu )**2 )
            integral *= parameters.scatteringFactor.evaluate( _x )
            if( parameters.energyWeight ) : integral *= kPrime_k
            if( parameters.calcualateMomentum ) : integral *= mu
            return( integral )

        def averageEnergy( parameters ) :

            parameters.energyWeight = False
            norm = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0]
            parameters.energyWeight = True
            energyPrime = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0]

            return( energyPrime / norm )

        fraction = 1.2
        energyMin, energyMax = kwargs['multiplicity'].domain( )
        N = int( math.log( energyMax / energyMin ) / math.log( fraction ) )
        energy = energyMin
        energies = []
        for i in range( N ) :
            energies.append( energy )
            energy *= fraction
        if( 1.2 * energy < energyMax ) : energies.append( energy )
        energies.append( energyMax )

        scatteringFactor = self.scatteringFactor.data

        factor_electronMass = 1.0 / PQUModule.PQU( 1.0, 'me * c**2' ).getValueAs( kwargs['incidentEnergyUnit'] )
        factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % kwargs['incidentEnergyUnit'] ).getValueAs( scatteringFactor.axes[1].unit )

        parameters = Parameters( )
        parameters.tolerance = 1e-2 * kwargs['energyAccuracy']
        parameters.scatteringFactor = scatteringFactor
        parameters.scatteringFactorDomainMax = parameters.scatteringFactor.domainMax

        aveEnergy = []
        aveMomenta = []
        for energy in energies :
            parameters.k = energy * factor_electronMass
            parameters.energy_hc = energy * factor_E2x

            parameters.calcualateMomentum = False
            aveEnergy.append( [ energy, energy * averageEnergy( parameters ) ] )

            parameters.calcualateMomentum = True
            aveMomenta.append( [ energy, energy * averageEnergy( parameters ) ] )

        return( [ aveEnergy ], [ aveMomenta ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule
        from fudge.processing import group as groupModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

        TM_1, TM_E = transferMatricesModule.comptonScattering( style, tempInfo, self.productFrame, self.scatteringFactor.data,
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def crossSection( self, domainMin, domainMax, domainUnit, tolerance=1e-3, interpolation=xDataEnumsModule.Interpolation.linlin):
        """
        Integrates the double differential cross section to get the cross section :math:`\\sigma(E)` where :math:`E` is the projectile's energy.
        """

        class Parameters :

            pass

        class Tester :

            def __init__( self, parameters, relativeTolerance, absoluteTolerance ) :

                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                return( integratedCrossSection( energy, self.parameters ) )

        def integrand( mu, parameters ) :

            oneMinusMu = 1.0 - mu
            _x = parameters.energy_hc * math.sqrt( 0.5 * oneMinusMu )
            kPrime_k = 1.0 / ( 1.0 + parameters.k * oneMinusMu )

            integral = kPrime_k * kPrime_k * ( 1.0 + mu * mu + kPrime_k * parameters.k**2 * oneMinusMu**2 )
            integral *= parameters.scatteringFactor.evaluate( _x )

            return( integral )

        def integratedCrossSection( energy, parameters ) :

            parameters.energy_hc = energy * parameters.factor_E2x
            parameters.k = energy * factor_electronMass
            return( parameters.constant * nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0] )

        scatteringFactor = self.scatteringFactor.data

        factor_electronMass = 1.0 / PQUModule.PQU( 1.0, 'me * c**2' ).getValueAs( domainUnit )
        factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % domainUnit ).getValueAs( scatteringFactor.axes[1].unit )

        parameters = Parameters( )
        tolerance = min( 0.1, max( tolerance, 1e-6 ) )
        parameters.tolerance = 1e-2 * tolerance
        parameters.factor_E2x = factor_E2x
        parameters.scatteringFactor = scatteringFactor
        electronRadius = PQUModule.PQU( 1.0, 'e**2 / ( 4 * pi * eps0 * me * c**2 )' ).getValueAs( '1e-12 * cm' )
        parameters.constant = math.pi * electronRadius**2

        domainMax = min( domainMax, scatteringFactor.domainMax / factor_E2x )
        energies = [ domainMin, math.sqrt( domainMin * domainMax ), domainMax ]
        energy_crossSection = []
        for energy in energies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters ) ] )
        _tester = Tester( parameters, tolerance, tolerance * energy_crossSection[-1][1] )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 20, interpolation = interpolation )

        return( crossSectionModule.XYs1d( data = energy_crossSection, axes = crossSectionModule.defaultAxes( domainUnit ),
                interpolation = interpolation ) )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param node:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        child = node[0]
        data = None
        if child.tag == ScatteringFactor.moniker:
            scatteringFactor = ScatteringFactor.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        elif child.tag == XYs1d.moniker:        # Pre GNDS 2.0.
            data = XYs1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        elif child.tag == Regions1d.moniker:    # Pre GNDS 2.0.
            data = Regions1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( child.tag, form.moniker ) )
        if data is not None:                    # Pre GNDS 2.0.
            scatteringFactor = ScatteringFactor(data)

        form = cls(node.get('pid'), node.get('label'), node.get('productFrame'), scatteringFactor)

        xPath.pop()

        return form
