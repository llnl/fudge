# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the nucleus class definitions.
The nucleus consists only of baryons, so quantities like spin, charge, etc. do not include atomic electrons.
For the compound particle formed of a nucleus + electrons, see the `nuclide` module.
"""

from .. import misc as miscModule

from .. import IDs as IDsModule


from ..chemicalElements import misc as chemicalElementMiscModule

from ..quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from . import particle as particleModule

from pqu import PQU as PQUModule

class Particle( particleModule.Particle ) :

    moniker = 'nucleus'
    familyOrder = 3

    def __init__( self, id, index ) :
        """
        :param id: Unique id for this nucleus (string). Per naming convention, should be similar to 'he4'
        :param index: Nuclear excited level index (int). 0 = ground state, 1 = 1st excited, etc.
        """

        particleModule.Particle.__init__( self, id )

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( id )
        if( not( isNucleus ) ) : raise ValueError( 'Invalid nucleus id = "%s"' % id )

        index = int( index )
        if( index != levelID ) : raise ValueError( 'id = "%s" does not agree with index = "%s"' % 
                ( miscModule.toLimitedString( id ), miscModule.toLimitedString( index ) ) )

        self.__chemicalElementSymbol = chemicalElementSymbol
        self.__Z = chemicalElementMiscModule.ZFromSymbol[chemicalElementSymbol]
        self.__A = chemicalElementMiscModule.checkA( A )
        self.__index = chemicalElementMiscModule.checkIndex( index )
        self.__energy = self.addSuite( nuclearEnergyLevelModule )

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        if( self.Z < other.Z ) : return( True )
        if( self.Z != other.Z ) : return( False )
        if self.A < other.A: return True
        if self.A != other.A: return False
        return self.index < other.index

    @property
    def A( self ) :

        return( self.__A )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__chemicalElementSymbol )

    @property
    def index( self ) :

        return( self.__index )

    @property
    def energy( self ) :

        return( self.__energy )

    @property
    def nuclide( self ) :

        return( self.ancestor )

    @property
    def Z( self ) :

        return( self.__Z )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        particleModule.Particle.convertUnits( self, unitMap )
        self.__energy.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _particle = Particle( self.id, self.index )
        self.__copyStandardQuantities( _particle )
        for item in self.__energy : _particle.energy.add( item.copy( ) )
        return( _particle )

    def getMass( self, unit ) :
        """
        Evaluate the nucleus mass in the desired unit. Mass is adjusted for excitation energy,
        but not for electron masses / binding energy.
        :param unit: desired unit (string)
        :return: mass (float)
        """

# Correction for electron mass and atomic binding energy according to Eq. (2) of AME2003 presented in Nuclear Physics A 729 (2003) 337–676.
        if len(self.mass) > 0:
            return self.mass[0].float(unit)
        else:
            nuclideMass = PQUModule.PQU(self.nuclide.getMass(unit), unit)
            Zme = PQUModule.PQU(self.Z, 'me')
            AtomicBindingEnergy = PQUModule.PQU( 14.4381 * pow(self.Z, 2.39) + 1.55468e-6 * pow(self.Z, 5.35), 'eV/c**2')
            nuclearMass = nuclideMass - Zme + AtomicBindingEnergy
            return nuclearMass.getValueAs(unit)

    def intid(self, intidDB={}):
        '''
        Converts the particle id into a unique integer dubbed an INTeger ID (INTID).
        '''

        sign = -1 if self.isAnti else 1

        return sign * (1000 * (1000 * (self.index + 500) + self.Z) + self.A)

    def replicate( self, other ) :
        """
        Copy data from other into self
        :param other: another nucleus.Particle instance
        """

        particleModule.Particle.replicate( self, other )
        self.__energy.replicate( other.energy )

    def extraXMLAttributes( self ) :

        return( ' index="%s"' % self.index )

    def extraXMLElements( self, indent, **kwargs ) :

        return( self.energy.toXML_strList( indent, **kwargs ) )

    def parseExtraXMLElement(self, element, xPath, linkData, **kwarg):

        if( element.tag == nuclearEnergyLevelModule.Suite.moniker ) :
            nuclearEnergyLevelModule.Suite.parseNode(self.energy, element, xPath, linkData, **kwarg)
            return( True )

        return( False )

    def sortCompare( self, other ) :

        if( not( isinstance( other, Particle ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.index - other.index )

    def particleCompare(self, other):
        """Compares *self* to particle *other* which can be from an particle family."""

        if self.familyOrder != other.familyOrder:
            return self.familyOrder - other.familyOrder
        if self.Z != other.Z:
            return self.Z - other.Z
        if self.A != other.A:
            return self.A - other.A
        return self.index - other.index
