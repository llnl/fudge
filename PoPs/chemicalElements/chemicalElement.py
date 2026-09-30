# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class for the suite chemicalElement.
"""

from .. import misc as miscModule

from ..atomic import atomic as atomicDataModule

from . import misc as chemicalElementMiscModule
from . import isotope as isotopeModule

class ChemicalElement( miscModule.ClassWithSymbolKey ) :
    """
    This class represents a chemical element (i.e., Oyxgen, Iron, Uranium). A chemical element
    stores all isotopes that contain the same number of protons (i.e., O16, O17, and O18). 
    The members of a chemicalElement are:

    .. table:: Chemical element members
       :widths: auto

       =========== =========== =========================================================================
       member      type        description
       =========== =========== =========================================================================
       symbol      string      This is the symbol for the chemical element (e.g., 'O' for oxygen).
       Z           int         Number of protons for this chemical element (e.g., 8 for 'O').
       name        string      The name of the chemical element (i.e., 'Oxygen')
       atomicData  atomicData  Container for atomic data, including ionized states and atomic relaxation
       isotopes    suite       The list of isotopes
       =========== =========== =========================================================================

    """

    moniker = 'chemicalElement'

    def __init__( self, symbol, Z, name ) :
        """
        :param symbol: Chemical symbol for the element, e.g. 'Fe' (string)
        :param Z: Atomic number, e.g. 26 (int)
        :param name: Full name of the element, e.g. 'Iron'
        """

        miscModule.ClassWithSymbolKey.__init__( self, symbol )

        base, anti = miscModule.baseAntiFromID( symbol )

        self.__anti = anti == miscModule.antiSuffix

        chemicalElementMiscModule.checkZ( Z )
        self.__Z = Z
        if chemicalElementMiscModule.symbolFromZ[Z] != base:
            raise ValueError('Z = %s and chemical element symbol = "%s" are not consistent' %
                ( Z, miscModule.toLimitedString( symbol ) ) )

        if not isinstance(name, str): TypeError( 'name must be a string' )
        if name not in chemicalElementMiscModule.namesFromZ[Z]:
            names = chemicalElementMiscModule.namesFromZ[Z]
            message = "%s" % names[0]
            if len(names) > 1:
                message = f"one of {', '.join(names)}."
            raise ValueError( 'Z = "%s" and chemical element name = "%s" are not consistent, should be %s.' %
                    (Z, miscModule.toLimitedString(name), message))
        self.__name = name

        self.__atomicData = None

        self.__isotopes = isotopeModule.Isotopes( )
        self.__isotopes.setAncestor( self ) 

    def __contains__( self, key ) : 

        return( key in self.__isotopes )

    def __getitem__( self, key ) :

        return( self.__isotopes[key] )

    def __iter__( self ) :
        """ Iterates over isotopes of this element (but not over the excited states of those isotopes) """

        for isotope in self.__isotopes : yield isotope

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def isAnti( self ) :

        return( self.__anti )

    @property
    def name( self ) :

        return( self.__name )

    @property
    def atomicData( self ) :

        return( self.__atomicData )

    @atomicData.setter
    def atomicData( self, value ) :

        if( not( isinstance( value, atomicDataModule.Atomic ) ) ) : raise TypeError( 'value must be an instance of atomicData.' )
        self.__atomicData = value
        self.__atomicData.setAncestor( self )

    @property
    def isotopes( self ) :

        return( self.__isotopes )

    def add( self, item ) :
        """
        :param item: isotope instance, will be added to the list of isotopes
        """

        self.__isotopes.add( item )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        if( self.__atomicData is not None ) : self.__atomicData.convertUnits( unitMap )
        self.__isotopes.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _chemicalElement = ChemicalElement( self.symbol, self.Z, self.name )
        if( self.__atomicData is not None ) : _chemicalElement.atomicData = self.__atomicData.copy( )
        for item in self.__isotopes : _chemicalElement.add( item.copy( ) )
        return( _chemicalElement )

    def sortCompare( self, other ) :

        if( not( isinstance( other, ChemicalElement ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.Z - other.Z )

    def check(self, info):
        """ See check documentation in PoPs.database """

        from .. import warning as warningModule
        warnings = []

        for isotope in self:
            isotopeWarnings = isotope.check(info)
            if isotopeWarnings:
                warnings.append(warningModule.Context('Isotope %s' % isotope.symbol, isotopeWarnings))

        return warnings

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s symbol="%s" Z="%s" name="%s">' % ( indent, self.moniker, self.symbol, self.Z, self.name ) ]
        if( self.__atomicData is not None ) : XMLStringList += self.__atomicData.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__isotopes.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        for child in element :
            if( child.tag == isotopeModule.Isotopes.moniker ) :
                self.__isotopes.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == atomicDataModule.Atomic.moniker ) :
                self.atomicData = atomicDataModule.Atomic.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise ValueError( 'Invalid child = "%s" for %s' % ( child.tag, self.moniker ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        self = cls( element.get('symbol'), int( element.get('Z') ), element.get('name') )
        xPath.pop( )

        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )
