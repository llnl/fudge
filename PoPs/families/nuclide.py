# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the nuclear level classes.
"""

from ..chemicalElements import misc as chemicalElementMiscModule
from ..fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

from . import particle as particleModule
from . import nucleus as nucleusModule

class Particle( particleModule.Particle ) :

    moniker = 'nuclide'
    familyOrder = 4

    def __init__( self, id ) :

        particleModule.Particle.__init__( self, id )

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( id )
        self.__nucleus = nucleusModule.Particle( chemicalElementMiscModule.nucleusIDFromNuclideID( id ), levelID )
        self.__nucleus.setAncestor(self)

        self.__fissionFragmentData = fissionFragmentDataModule.FissionFragmentData( )
        self.__fissionFragmentData.setAncestor( self )

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )
        if( self.familyOrder != other.familyOrder ) : return( False )
        if( self.Z < other.Z ) : return( True )
        if( self.Z != other.Z ) : return( False )
        if self.A < other.A: return True
        if self.A != other.A: return False
        return self.index < other.index

    @property
    def fissionFragmentData( self ) :

        return( self.__fissionFragmentData )

    @property
    def nucleus( self ) :

        return( self.__nucleus )

    @property
    def A( self ) :

        return( self.__nucleus.__A )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__nucleus.chemicalElementSymbol )

    @property
    def index( self ) :

        return( self.__nucleus.index )

    @property
    def isotope( self ) :

        return( self.ancestor.ancestor )

    @property
    def energy( self ) :

        return( self.__nucleus.energy )

    @property
    def Z( self ) :

        return( self.__nucleus.Z )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        subWarnings = self.__nucleus.check(info)
        if subWarnings:
            warnings.append( warningModule.Context('nucleus', subWarnings) )
        # FIXME other checks to perform on the nuclide? Will decay info ever live in the nuclide?

        return warnings

    def convertUnits( self, unitMap ) :

        particleModule.Particle.convertUnits( self, unitMap )
        self.__nucleus.convertUnits( unitMap )
        self.__fissionFragmentData.convertUnits( unitMap )

    def copy( self ) :

        _particle = Particle( self.id )
        self.__copyStandardQuantities( _particle )
        _particle.__nucleus.replicate( self.__nucleus )
        _particle.__fissionFragmentData.replicate( self.__fissionFragmentData )

        return( _particle )

    def extraXMLElements( self, indent, **kwargs ) :

        XMLStringList = self.__nucleus.toXML_strList( indent, **kwargs )
        XMLStringList += self.__fissionFragmentData.toXML_strList( indent, **kwargs )

        return( XMLStringList )

    def getMass(self, unit):

        if len(self.mass) > 0: return self.mass[0].float(unit)
        if self.index == 0:
            raise Exception('Recursion detected as ground-state does not have a mass: ID = %s.' % self.id)
        return self.ancestor[0].mass[0].float(unit) + self.energy[0].float(unit + ' * c**2')

    def intid(self, intidDB={}):
        '''
        Converts the particle id into a unique integer dubbed an INTeger ID (INTID).
        '''

        sign = -1 if self.isAnti else 1

        return sign * (1000 * (1000 * self.index + self.Z) + self.A)

    def parseExtraXMLElement(self, element, xPath, linkData, **kwargs):

        if( element.tag == nucleusModule.Particle.moniker ) :
            nucleus = nucleusModule.Particle.parseNodeUsingClass(element, xPath, linkData, **kwargs)
            self.__nucleus.replicate( nucleus )
            return( True )
        elif( element.tag == fissionFragmentDataModule.FissionFragmentData.moniker ) :
            fissionFragmentData = fissionFragmentDataModule.FissionFragmentData.parseNodeUsingClass(element, xPath, linkData, **kwargs)
            self.__fissionFragmentData.replicate( fissionFragmentData )
            return( True )

        return( False )

    def sortCompare(self, other):

        if not isinstance(other, Particle):
            if isinstance(other, particleModule.Particle):
                raise TypeError('Invalid other of instance %s is not an instance of %s.' % (other.moniker, self.moniker))
            raise TypeError('Invalid other, not a particle.')

        if self.index - other.index == 0:
            if self.isAnti:
                return 0 if other.isAnti else 1
            return -1 if other.isAnti else 0

        return self.index - other.index

    def particleCompare(self, other):
        """Compares *self* to particle *other* which can be from an particle family."""

        if self.familyOrder != other.familyOrder:
            return self.familyOrder - other.familyOrder
        if self.Z != other.Z:
            return self.Z - other.Z
        if self.A != other.A:
            return self.A - other.A
        return self.index - other.index

class Suite( particleModule.Suite ) :

    moniker = 'nuclides'
    particle = Particle
