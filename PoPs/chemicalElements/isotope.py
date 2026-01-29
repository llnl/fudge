# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the Isotope and Isotopes classes.
"""

from .. import suite as suiteModule
from .. import misc as miscModule

from ..chemicalElements import misc as chemicalElementMiscModule

from ..families import nuclide as nuclideModule
from ..families import nucleus as nucleusModule

class Isotope( miscModule.ClassWithSymbolKey ) :
    """
    An isotope that contains a list of nuclides. The symbol is of the form  *SA[_anti]*
    where *S* is the istopes's chemical element symbol and *A* is the isotope's nucleon number.
    """

    moniker = 'isotope'

    def __init__( self, symbol, A ) :
        """
        :param symbol: string consisting of the chemical symbol + total number of nucleons, e.g. 'Mn55'
        :param A: integer equal to the total number of nucleons, e.g. 55
        """

        miscModule.ClassWithSymbolKey.__init__( self, symbol )

        baseID, chemicalElementSymbol, _A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( symbol )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        if( symbol != isotopeSymbol ) : ValueError( 'invalid isotope symbol = "%s"' % symbol )

        if( A != _A ) : raise ValueError( 'symbol = "%s" does not agree with A = "%s"' % 
                ( miscModule.toLimitedString( symbol ), miscModule.toLimitedString( A ) ) )

        self.__chemicalElementSymbol = chemicalElementSymbol
        self.__A = A

        self.__Z = chemicalElementMiscModule.ZFromSymbol[chemicalElementSymbol]

        self.__atomicData = None

        self.__nuclides = nuclideModule.Suite( )
        self.__nuclides.setAncestor( self )

    def __len__( self ) :

        return( len( self.__nuclides ) )

    def __contains__( self, key ) :

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A, levelID ] ) : return( False )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        nuclideID = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeSymbol, levelID )
        if( nuclideID in self.__nuclides ) :
            if( not( isNucleus ) ) : return( True )
            nuclide = self.__nuclides[nuclideID]
            return( nuclide.nucleus is not None )
        return( False )

    def __getitem__( self, key ) :
        """
        Get child particle by id or by integer index.
        :param key: may be a string (corresponding to a nuclide or nucleus inside the isotope) or an integer.
        :return: nuclide or nucleus instance, or returns KeyError if key not found
        """

        if( isinstance( key, int ) ) : return( self.__nuclides[0] ) # FIXME shouldn't that be self.__nuclides[key]?
        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A, levelID ] ) : raise KeyError( 'key "%s" not found' % key )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        nuclideID = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeSymbol, levelID ) + anti
        if( nuclideID in self.__nuclides ) :
            nuclide = self.__nuclides[nuclideID]
            if( not( isNucleus ) ) : return( nuclide )
            if( nuclide.nucleus is not None ) : return( nuclide.nucleus )

        raise KeyError( 'key "%s" not found' % key )

    def __iter__( self ) :
        """ Iterates over nuclides, but not over the nuclei contained inside each nucleus """

        for nuclide in self.__nuclides : yield nuclide

    @property
    def A( self ) :

        return( self.__A )

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__chemicalElementSymbol )

    @property
    def nuclides( self ) :

        return( self.__nuclides )

    def add( self, particle ) :

        self.__nuclides.add( particle )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        emax = -1
        for nuclide in self:
            nuclideWarnings = nuclide.check(info)
            if nuclideWarnings:
                warnings.append(warningModule.Context('nuclide %s' % nuclide.id, nuclideWarnings))

            if len(nuclide.nucleus.energy) == 0:
                warnings.append(warningModule.UnknownEnergy(nuclide.nucleus))
                continue

            enow = nuclide.nucleus.energy.float('eV')
            if enow <= emax:
                warnings.append(warningModule.DiscreteLevelsOutOfOrder(nuclide.nucleus.index, nuclide))
            else:
                emax = enow

        return warnings

    def convertUnits( self, unitMap ) :

        if( self.__atomicData is not None ) : self.__atomicData.convertUnits( unitMap )
        self.__nuclides.convertUnits( unitMap )

    def copy( self ) :

        _isotope = Isotope( self.symbol, self.A )
        for nuclide in self.__nuclides : _isotope.add( nuclide.copy( ) )
        return( _isotope )

    def sortCompare( self, other ) :

        if( not( isinstance( other, Isotope ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.A - other.A )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s symbol="%s" A="%s">' % ( indent, self.moniker, self.symbol, self.A ) ]
        if( self.__atomicData is not None ) : XMLStringList += self.__atomicData.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.__nuclides.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )
    
        for child in element :
            if( child.tag == nuclideModule.Suite.moniker ) :
                self.nuclides.parseNode(child,  xPath, linkData, **kwargs)
            else :
                raise Exception( 'Fix me: tag = %s' % child.tag )
    
        xPath.pop( )
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        self = cls( element.get( 'symbol' ), int( element.get( 'A' ) ) )
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )

class Isotopes( suiteModule.SortedSuite ) :

    moniker = 'isotopes'

    def __init__( self, replace = True ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( Isotope, ), replace = replace )

    def __contains__( self, key ) :

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A ] ) : return( False )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        if( suiteModule.SortedSuite.__contains__( self, isotopeSymbol ) ) :
            return( key in suiteModule.SortedSuite.__getitem__( self, isotopeSymbol ) )
        return( False )

    def __getitem__( self, key ) :

        if( isinstance( key, int ) ) : return( suiteModule.SortedSuite.__getitem__( self, key ) )
        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A ] ) : raise KeyError( 'key "%s" not found' % key )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        return( suiteModule.SortedSuite.__getitem__( self, isotopeSymbol )[key] )

    def add( self, item ) :

        if( isinstance( item, Isotope ) ) :
            suiteModule.SortedSuite.add( self, item )
        elif( isinstance( item, ( nuclideModule.Particle, nucleusModule.Particle ) ) ) :
            isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( item.chemicalElementSymbol, item.A )
            if( not( suiteModule.SortedSuite.__contains__( self, isotopeSymbol ) ) ) : suiteModule.SortedSuite.add( self, Isotope( isotopeSymbol, item.A ) )
            suiteModule.SortedSuite.__getitem__( self, isotopeSymbol ).add( item )
        else :
            raise TypeError( "Object not an Isotope, Nuclide or Nucleus particle." )

    def convertUnits( self, unitMap ) :

        for isotope in self : isotope.convertUnits( unitMap )

    def getSymbol( self, key ) :

        return( suiteModule.SortedSuite.__getitem__( self, key ) )

    def hasSymbol( self, key ) :

        return( suiteModule.SortedSuite.__contains__( self, key ) )
