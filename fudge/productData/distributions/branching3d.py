# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module contains the 'branching3d' form for distribution."""

"""
This module contains class and functions needed for storing and processing a GNDS *branching3d* distribution. It
contains the following classes:

    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+
    | Form          | Class representing GNDS unspecified distribution.                     |
    +---------------+-----------------------------------------------------------------------+
"""

from PoPs import IDs as IDsPoPsModule

from xData import enums as xDataEnumsModule
from fudge.processing import group as groupModule

from . import base as baseModule

class Form( baseModule.Form ) :
    """
    Class representing a GNDS unspecified distribution form.

    :param label:                   The label for this form.
    :param productFrame:            The frame the product data are specified in.
    """

    moniker = 'branching3d'
    subformAttributes = [ ]

    def __init__(self, label, productFrame=xDataEnumsModule.Frame.lab):

        baseModule.Form.__init__( self, label, productFrame, [] )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        pass

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :
        """         
        This method calculates the average energy and momentum of the outgoing particle as a function of projectile energy.
        Average means over all outgoing angle and energy.
        
        :param style:   The style instance which the calculated data will belong to.
        :param indent:  If this method does any printing, this is the amount of indentation of the printed line.
        :param kwargs:  A dictionary that contains data not in *self* that is needed to calculate the average energy and momentum.
        """

        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = kwargs['momentumUnit']
        domainMin = kwargs['EMin']
        domainMax = kwargs['EMax']

        Q = kwargs['outputChannel'].Q[0]
        QValue = Q.evaluate( domainMin )

        energyData = [ [ [ domainMin, QValue ], [ domainMax, QValue ] ] ]       # Gives all the energy to the gammas and non to the residual.
        momentumData = [ [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ] ]           # Assumes isotropic scattering.

        return( [ energyData, momentumData ] )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        The method does nothing.
        """

        return 0

    def check(self, info):
        """Check that the excited state supplying this photon cascade has decay data."""

        from fudge import warning as warningModule

        parentProduct = self.product.parentProduct
        if parentProduct is None:
            return []

        particle = info['reactionSuite'].PoPs[parentProduct.pid]
        if getattr(particle, 'index', 0) > 0 and len(particle.decayData.decayModes) == 0:
            return [warningModule.MissingDecayData(parentProduct.pid, self)]

        return []

    def processMC_cdf( self, style, tempInfo, indent ) :
        """
        This methods returns does nothing and returns None.
    
        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.
    
        :return:                None
        """     

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Returns a multi-group representation of *self*.

        :param style:           The style for the returned data.
        :param tempInfo:        Dictionary of data needed to calculate the data.
        :param indent:          The indentation for any verbosity.

        :return:                A multi-group representation of *self*.
        """

        from .. import multiplicity as multiplicityModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        def processDecayModesSetup( product, probability, gammas ) :

            initialEnergy = product.nucleus.energy[0].value
            for decayMode in product.decayData.decayModes :
                branchingRatio = decayMode.probability[0].value * probability

                decayPath = decayMode.decayPath[0]
                residual = decayPath.products[0].pid
                if( residual == IDsPoPsModule.photon ) : residual = decayPath.products[1].pid
                residual = PoPs[residual]

                gammaEnergy = initialEnergy - residual.nucleus.energy[0].value
                if( gammaEnergy not in gammas ) : gammas[gammaEnergy] = 0.0
                gammas[gammaEnergy] += branchingRatio

                processDecayModesSetup( residual, branchingRatio, gammas )

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        productLabel = tempInfo['productLabel']

        if( verbosity > 2 ) : print( '%sGrouping %s' % ( indent, self.moniker ) )

        crossSection = tempInfo['crossSection']
        domainMin = crossSection.domainMin
        domainMax = crossSection.domainMax

        PoPs = tempInfo['reactionSuite'].PoPs
        multiplicity = self.ancestor.ancestor.multiplicity[self.label]
        nuclide = PoPs[self.product.parentProduct.pid]

        gammas = {}
        processDecayModesSetup( nuclide, 1.0, gammas )

        TM_1s = []
        TM_Es = []
        multipliticyAxes = multiplicityModule.defaultAxes( tempInfo['incidentEnergyUnit'] )
        if len(gammas) == 0:                # Kludge to handle case where there no gamma is emitted.
            gammas[0.0] = 0.0
        for gammaEnergy in gammas :
            branchingRatio = gammas[gammaEnergy]
            multiplicity = multiplicityModule.XYs1d( data = [ [ domainMin, branchingRatio ], [ domainMax, branchingRatio ] ], axes = multipliticyAxes )

            TM_1, TM_E = transferMatricesModule.discreteGammaAngularData( style, tempInfo, gammaEnergy, crossSection, None, 
                    multiplicity = multiplicity, comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            TM_1s.append( TM_1 )
            TM_Es.append( TM_E )
        TM_1 = transferMatricesModule.addTMs( TM_1s )
        TM_E = transferMatricesModule.addTMs( TM_Es )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of self.
    
        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame

        xmlString = [ '%s<%s%s/>' % ( indent, self.moniker, attributeStr ) ]

        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into :py:class:`Form` instance.

        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )

        _form = cls( element.get( 'label' ), productFrame = element.get( 'productFrame' ) )

        xPath.pop( )
        return( _form )
