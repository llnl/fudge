# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains Angular/Energy double differential distribution classes used for storing
distribution data that have been pre-processed for use in Monte Carlo transport codes.
The distribution is a product of an angular probability P(mu|E) and an energy  probability P(E'|E,mu) where
E is the projectile's energy, and mu and E' are the product's angular and energy variables. The inter most
1d functions P(mu) for P(mu|E) and P(E') for P(E'|E,mu) are stored with an xs, pdf, and cdf
(i.e., :py:class:`Xs_pdf_cdf1d`) function.
"""

from xData import xs_pdf_cdf as xs_pdf_cdfModule
from xData import multiD_XYs as multiD_XYsModule

from . import base as baseModule
from . import angular as angularModule


class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :
    """
    Class for storing the inter most of P(E'|E,mu) for each mu.
    """

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    Class for storing a list of :py:class:`Xs_pdf_cdf1d` instances for each projectile energy E.
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., a 1d function) of an :py:class:`XYs2d` instance.
        """

        return( ( Xs_pdf_cdf1d, ) )

class XYs3d( multiD_XYsModule.XYs3d ) :
    """
    Class for storing a list of :py:class:`XYs2d` instances (i.e., the P(E'|E,mu)).
    """

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns the list of classes that can be sub-nodes (i.e., 2d function) of an :py:class:`XYs3d` instance.
        """

        return( ( XYs2d, ) )

class Subform( baseModule.Subform ) :
    """
    Abstract base class for angular and energy sub-nodes of an AngularEnergyMC instance.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | data      | This member stores the data (i.e., function).             |
    +-----------+-----------------------------------------------------------+

    :param data:    Data (i.e, function) for the angular and energy sub-nodes of an AngularEnergyMC instance.
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

class Angular( Subform ) :
    """
    Class for storing the angular data (i.e., P(mu|E)).
    """

    moniker = 'angular'
    allowedSubElements = ( angularModule.XYs2d, )
    ancestryMembers = ( 'angular', )

    @property
    def angular(self):

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
        subformClass = {angularModule.XYs2d.moniker: angularModule.XYs2d}.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angular subform "%s"' % subformElement.tag )
        angularSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _angular= cls(angularSubform)

        xPath.pop( )
        return( _angular )

class AngularEnergy( Subform ) :
    """
    Class for storing the energy data (i.e., P(E'|E,mu)).
    """

    moniker = 'angularEnergy'
    allowedSubElements = ( XYs3d, )
    ancestryMembers = ( 'angularEnergy', )

    @property
    def angularEnergy(self):

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
        subformClass = {XYs3d.moniker: XYs3d}.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( 'unknown angularEnergy subform "%s"' % subformElement.tag )
        angularEnergySubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _angularEnergy = cls(angularEnergySubform)

        xPath.pop( )
        return( _angularEnergy )

class Form( baseModule.Form ) :
    """
    Class for storing the distribution that is a product of an angular probability P(mu|E) and an energy 
    probability P(E'|E,mu) where E is the projectile's energy, and mu and E' are the product's angular and energy variables.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | label         | Unique label for the form within the distribution instance.   |
    +---------------+---------------------------------------------------------------+
    | productFrame  | The frame the product data are in.                            |
    +---------------+---------------------------------------------------------------+
    | _angular      | The angular child node containing P(mu|E).                    |
    +---------------+---------------------------------------------------------------+
    | _angularEnergy| The energy child node containing P(E'|E,mu).                  |
    +---------------+---------------------------------------------------------------+

    :param label:           Unique label for the form within the distribution instance.
    :param productFrame:    The frame the product data are in.
    :param _angular:        The angular child node containing P(mu|E).
    :param _angularEnergy:  The energy child node containing P(E'|E,mu).
    """

    moniker = 'angularEnergyMC'
    subformAttributes = ( 'angular', 'angularEnergy' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _angular, _angularEnergy ) :

        if( not( isinstance( _angular, Angular ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angular ) )
        if( not( isinstance( _angularEnergy, AngularEnergy ) ) ) : raise TypeError( 'Invalid instance: %s' % type( _angularEnergy ) )
        baseModule.Form.__init__( self, label, productFrame, ( _angular, _angularEnergy ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.angular.convertUnits( unitMap )
        self.angularEnergy.convertUnits( unitMap )

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

        _angular = None
        _angularEnergy = None
        for child in node:
            if( child.tag == Angular.moniker ) :
                _angular = Angular.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == AngularEnergy.moniker ) :
                _angularEnergy = AngularEnergy.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise TypeError( "Encountered unknown yAngular subform: %s" % child.tag )

        angularEnergyMC = cls(node.get("label"), node.get('productFrame'), _angular, _angularEnergy)

        xPath.pop( )
        return angularEnergyMC
