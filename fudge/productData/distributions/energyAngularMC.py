# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains Energy/angular double differential distribution classes used for storing
distribution data that have been pre-processed for use in Monte Carlo transport codes.
The distribution is a product of an energy probability P(E'|E) and an angular probability P(mu|E,E') where
E is the projectile's energy, and mu and E' are the product's angular and energy variables. The inter most
1d functions P(mu) for P(mu|E) and P(E') for P(E'|E,mu) are stored with an xs, pdf, and cdf
(i.e., :py:class:`Xs_pdf_cdf1d`) function.
"""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import energy as energyModule


class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :
    """
    Class for storing the inter most of P(mu|E,E') for each E'.
    """

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    Class for storing a list of :py:class:`Xs_pdf_cdf1d` instances for each projectile's energy.
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs2d` instance.
        """

        return( ( Xs_pdf_cdf1d, ) )

class XYs3d( multiD_XYsModule.XYs3d ) :
    """
    Class for storing a list of :py:class:`XYs2d` instances (i.e., the P(mu|E,E')).
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 1-d function) of an :py:class:`XYs3d` instance.
        """

        return( ( XYs2d, ) )

class Subform( baseModule.Subform ) :
    """
    Abstract base class for energy and angular sub-nodes of an EnergyAngularMC instance.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the data (i.e., function).             |
    +-----------+-----------------------------------------------------------+

    :param data:    Data (i.e, function) for the angular and energy sub-nodes of an EnergyAngularMC instance.
    """

    moniker = 'dummy'                   # This is not used but added to stop lylint from reporting.

    def __init__( self, data ) :

        baseModule.Subform.__init__( self )
        if( not( isinstance( data, self.allowedSubElements ) ) ) : raise TypeError( 'Invalid instance: %s' % type( data ) )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.data.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.data.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

class Energy( Subform ) :
    """ 
    Class for storing the energy data (i.e., P(E'|E)).
    """

    moniker = 'energy'
    allowedSubElements = ( energyModule.XYs2d, )
    ancestryMembers = ( 'energy', )

    @property
    def energy(self):

        return self.data

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {energyModule.XYs2d.moniker: energyModule.XYs2d}.get(subformElement.tag)
        if( subformClass is None ) : raise Exception( 'unknown energy subform "%s"' % subformElement.tag )
        energySubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _energy = cls(energySubform)

        xPath.pop( )
        return( _energy )

class EnergyAngular( Subform ) :
    """ 
    Class for storing the energy data (i.e., P(mu|E,E')).
    """

    moniker = 'energyAngular'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'energyAngular', )

    @property
    def energyAngular(self):

        return self.data

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        subformElement = node[0]
        subformClass = {XYs3d.moniker: XYs3d}.get(subformElement.tag)
        if( subformClass is None ) : raise Exception( 'unknown energyAngular subform "%s"' % subformElement.tag )
        energyAngularSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _energyAngular = cls(energyAngularSubform)

        xPath.pop( )
        return( _energyAngular )

class Form( baseModule.Form ) :
    """
    Class for storing the distribution that is a product of an energy probability P(E'|E) and an angular
    probability P(mu|E,E') where E is the projectile's energy, and mu and E' are the product's angular and energy variables.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | label         | Unique label for the form within the distribution instance.   |
    +---------------+---------------------------------------------------------------+
    | productFrame  | The frame the product data are in.                            |
    +---------------+---------------------------------------------------------------+
    | _energy       | The energy child node containing P(E'|E).                     |
    +---------------+---------------------------------------------------------------+
    | _energyAngular| The angular child node containing P(mu|E,E').                 |
    +---------------+---------------------------------------------------------------+

    :param label:           Unique label for the form within the distribution instance.
    :param productFrame:    The frame the product data are in.
    :param _energy:         The energy child node containing P(E'|E).
    :param _energyAngular:  The angular child node containing P(mu|E,E').
    """

    moniker = 'energyAngularMC'
    subformAttributes = ( 'energy', 'energyAngular' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _energy, _energyAngular ) :

        if( not( isinstance( _energy, Energy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energy ) )
        if( not( isinstance( _energyAngular, EnergyAngular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _energyAngular ) )
        baseModule.Form.__init__( self, label, productFrame, ( _energy, _energyAngular ) )

    @property
    def domainMin( self ) :
        """Returns the minimum projectile energy for the energy child node."""

        return( self.energy.data.domainMin )

    @property
    def domainMax( self ) :
        """Returns the maximum projectile energy for the energy child node."""

        return( self.energy.data.domainMax )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        _energy = None
        _energyAngular = None
        for child in node:
            if( child.tag == Energy.moniker ) :
                _energy = Energy.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == EnergyAngular.moniker ) :
                _energyAngular = EnergyAngular.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        energyAngularMC = cls(node.get("label"), node.get("productFrame"), _energy, _energyAngular)

        xPath.pop( )
        return( energyAngularMC )
