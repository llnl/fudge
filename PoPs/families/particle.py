# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the abstract base class for all types of particle appearing in PoPs.
"""

import abc

from .. import misc as miscModule
from .. import suite as suiteModule
from .. import warning as warningModule

from ..quantities import quantity as quantityModule
from ..quantities import mass as massModule
from ..quantities import spin as spinModule
from ..quantities import parity as parityModule
from ..quantities import charge as chargeModule
from ..quantities import halflife as halflifeModule

from ..decays import decayData as decayDataModule

class Particle(miscModule.ClassWithIDKey, abc.ABC):
    """
    This is the abstract base class for all particles.
    """

    def __init__( self, id ) :
        """
        :param id: The particle id. Every particle in a PoPs database must have a unique id.
        """

        miscModule.ClassWithIDKey.__init__( self, id )

        base, anti, qualifier = miscModule.baseAntiQualifierFromID( id )

        self.__anti = anti == miscModule.antiSuffix

        self.__mass = self.addSuite( massModule )
        self.__spin = self.addSuite( spinModule )
        self.__parity = self.addSuite( parityModule )
        self.__charge = self.addSuite( chargeModule )
        self.__halflife = self.addSuite( halflifeModule )
        self.__decayData = decayDataModule.DecayData( )
        self.__decayData.setAncestor( self )

    def __eq__( self, other ) :

        return( self.id == other.id )

    @property
    @abc.abstractmethod
    def familyOrder( self ) :

        pass

    @property
    def isAnti( self ) :
        """
        :return: whether this is an anti-particle
        """

        return( self.__anti )

    @property
    def family( self ) :
        """
        :return: family name of the particle
        """

        return( self.moniker )

    @property
    def mass( self ) :
        """
        :return: particle mass
        """

        return( self.__mass )

    @property
    def spin( self ) :
        """
        :return: particle spin
        """

        return( self.__spin )

    @property
    def parity( self ) :
        """
        :return: particle parity
        """

        return( self.__parity )

    @property
    def charge( self ) :
        """
        :return: particle charge
        """

        return( self.__charge )

    @property
    def halflife( self ) :
        """
        :return: particle halflife
        """

        return( self.__halflife )

    @property
    def decayData( self ) :
        """
        :return: particle decay information
        """

        return( self.__decayData )

    def addSuite( self, module ) :      # FIXME rename this __addSuite? Not intended for public use

        suite = module.Suite( )
        suite.setAncestor( self )
        return( suite )

    def check( self, info ):
        """
        Basic physics checking for particles.
        """
        warnings = []

        decayWarnings = self.__decayData.check(info)
        if decayWarnings:
            warnings.append(warningModule.Context('decayData', decayWarnings))

        return warnings

    def buildFromRawData( self, mass = None, spin = None, parity = None, charge = None, halflife = None, label = 'default' ) :
        """
        Helper method to conviniently add the common properties to a particle.

        :param mass:  tuple(float value, string unit)
        :param spin:  tuple(int or fraction value, string unit)
        :param parity: tuple(int value, string unit)
        :param charge: tuple(int value, string unit)
        :param halflife: tuple(float or string value, string unit)
        :param label: style label (string) to apply to each quantity
        """

        def getUnit( unit ) :

            if( isinstance( unit, str ) ) : unit = quantityModule.stringToPhysicalUnit( unit )
            return( unit )

        if( mass is not None ) : self.mass.add( massModule.Double( label, mass[0], getUnit( mass[1] ) ) )
        if( spin is not None ) :
            self.spin.add( spinModule.Fraction( label, spin[0], spin[1] ) )
        if( parity is not None ) : self.parity.add( parityModule.Integer( label, parity[0], getUnit( parity[1] ) ) )
        if( charge is not None ) : self.charge.add( chargeModule.Integer( label, charge[0], getUnit( charge[1] ) ) )
        if( halflife is not None ) :
            if( isinstance( halflife[0], float ) ) :
                _halflife = halflifeModule.Double( label, halflife[0], getUnit( halflife[1] ) )
            elif( isinstance( halflife[0], str ) ) :
                _halflife = halflifeModule.String( label, halflife[0], getUnit( halflife[1] ) )
            else :
                TypeError( 'halflife neither double or str instance' )
            self.halflife.add( _halflife )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        self.__mass.convertUnits( unitMap )
        self.__spin.convertUnits( unitMap )
        self.__parity.convertUnits( unitMap )
        self.__charge.convertUnits( unitMap )
        self.__halflife.convertUnits( unitMap )
        self.__decayData.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _particle = self.__class__( self.id )
        self.__copyStandardQuantities( _particle )
        return( _particle )

    def __copyStandardQuantities( self, other ) :

        if( other.moniker != self.moniker ) : TypeError( '''other's family = "%s" not same as self's family = "%s".''' %
            ( other.moniker, self.moniker ) )
        for item in self.__mass : other.mass.add( item.copy( ) )
        for item in self.__spin : other.spin.add( item.copy( ) )
        for item in self.__parity : other.parity.add( item.copy( ) )
        for item in self.__charge : other.charge.add( item.copy( ) )
        for item in self.__halflife : other.halflife.add( item.copy( ) )
        self.__decayData.copyItems( other.decayData )

    def familyOrderLessThan( self, other ) :

        if( not( isinstance( other, Particle ) ) ) : raise TypeError( 'Other must be a Particle: %s.' % other.__class__ )
        if( self.familyOrder < other.familyOrder ) : return( True )
        return( False )

    def replicate( self, other ) :
        """
        Copy data from other into self

        :param other: another Particle instance
        """

        self.__mass.replicate( other.mass )
        self.__spin.replicate( other.spin )
        self.__parity.replicate( other.parity )
        self.__charge.replicate( other.charge )
        self.__halflife.replicate( other.halflife )
        self.__decayData = other.decayData.copy( )

    def getMass( self, unit ) :
        """
        Evaluate the particle mass in the desired unit.

        :param unit: desired unit (string)
        :return: mass (float)
        """

        return( self.mass[0].float( unit ) )

    def extraXMLAttributes( self ) :

        return( '' )

    def extraXMLElements( self, indent = '', **kwargs ) :

        return( [] )

    def isAlias( self ) :

        return( False )

    def sortCompare( self, other ) :
        """This should be defined in each sub-class."""

        if( not( isinstance( other, self.__class__ ) ) ) : raise TypeError( 'Invalid other.' )
        if( self.id > other.id ) : return( 1 )
        if( self.id == other.id ) : return( 0 )
        return( -1 )

    def particleCompare(self, other):
        """Compares *self* to particle *other* which can be from an particle family."""

        if self.familyOrder != other.familyOrder:
            return self.familyOrder - other.familyOrder
        if self.id > other.id: return 1
        if self.id == other.id: return 0
        return -1 

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s id="%s"%s>' % ( indent, self.moniker, self.id, self.extraXMLAttributes( ) ) ]

        XMLStringList += self.mass.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.spin.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.parity.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.charge.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.halflife.toXML_strList( indent = indent2, **kwargs )

        XMLStringList += self.decayData.toXML_strList( indent = indent2, **kwargs )

        XMLStringList += self.extraXMLElements( indent2, **kwargs )

        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseExtraXMLElement(self, element, xPath, linkData, **kwarg):

        return( False )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@id="%s"]' % ( element.tag, element.get( 'id' ) ) )

        children = { 'mass' : massModule,       'spin' : spinModule,            'parity' : parityModule, 
                     'charge' : chargeModule,   'halflife' : halflifeModule }
        for child in element :
            if( child.tag in children ) :
                children[child.tag].Suite.parseNode(getattr(self, child.tag), child, xPath, linkData, **kwargs)
            elif( child.tag == decayDataModule.DecayData.moniker ) :  # not a suite so can't be treated like other children
                self.decayData.parseNode(child, xPath, linkData, **kwargs)
            else :
                if not self.parseExtraXMLElement(child, xPath, linkData, **kwargs):
                    raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def buildFromClassAndRawData(cls, id, mass = None, spin = None, parity = None, charge = None, halflife = None, label = 'default'):
        """
        Helper method to create a particle from class *cls* and to add common properties to the particle.

        :param id:  the string id of the partical
        :param mass:  tuple(float value, string unit)
        :param spin:  tuple(int or fraction value, string unit)
        :param parity: tuple(int value, string unit)
        :param charge: tuple(int value, string unit)
        :param halflife: tuple(float or string value, string unit)
        :param label: style label (string) to apply to each quantity
        """

        particle = cls(id)
        particle.buildFromRawData(mass=mass, spin=spin, parity=parity, charge=charge, halflife=halflife, label=label)

        return particle

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@id="%s"]' % ( element.tag, element.get( 'id' ) ) )

        attrs = {v[0]:v[1] for v in element.items()}
        del attrs['id']
        self = cls( element.get('id'), **attrs )
        xPath.pop()

        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )

class Suite( suiteModule.SortedSuite ) :

    def __init__( self, replace = True ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( self.particle, ), replace = replace )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            self.add(self.particle.parseNodeUsingClass(child,  xPath, linkData, **kwargs))

        xPath.pop( )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        self = cls( )
        for child in element :
            self.add(self.particle.parseNodeUsingClass(child,  xPath, linkData, **kwargs))

        xPath.pop( )
        return( self )

class ParticleSorter:
    """
    """

    def __init__(self, overwrite=False):
        """Constructor."""

        self.__overwrite = overwrite
        self.__particles = []

    def __len__(self):
        """Returns the number of particles in *self*."""

        return len(self.__particles)

    def __getitem__(self, index):
        """Returns the particle(s) specified by *index*."""

        return self.__particles[index]

    def __iter__(self):
        """Iterates over the particles in *self*."""

        n1 = len(self)
        for i1 in range(n1): yield self.__particles[i1]

    @property
    def overwrite(self):
        """Returns the value of self.__overwrite."""

        return self.__overwrite

    def add(self, particle):
        """Determines the location of *particle* in *self* and inserts it at that location."""

        if not isinstance(particle, Particle):
            raise TypeError('Invalid particle of type "%s".' % type(particle))

        cmp = 0
        index = 0
        for index, item in enumerate(self):
            cmp = particle.particleCompare(item)
            if cmp == 0:
                if not self.__overwrite:
                    raise Exception('ParticleSorter instance does not allow overwriting (i.e., replacing) a particle (id = "%s"). %s %s' % (particle.id, item.id, cmp))
                del self.__particles[index]
                break
            elif cmp < 0:
                break
        if cmp > 0: index += 1

        self.__particles.insert(index, particle)
