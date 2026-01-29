# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
Define some helper methods and some base classes. These are mostly meant for internal use by the PoPs module.
"""

from LUPY import ancestry as ancestryModule

maxLength = 32
antiSuffix = '_anti'

def toLimitedString( object, maxLength = maxLength ) :
    """
    Returns a str instance of object that is truncated to no more than 32 characters.
    """

    string = str( object )
    if( len( string ) <= maxLength ) : return( string )
    return( string[:maxLength-4] + ' ...' )

def parseParticleList(particleList):
    '''
    Breaks a particle list where each particle id can have a leading integer multiplicity into an equivalent list of 
    particle ids.  That is, each particle in *particleList* to be of the form "id" or "#id" where "#" is an integer.
    Examples of the form "#id" include "2n" and "3n". For example, the paticle list ["2n", "Pu238"] is returned as ["n", "n", "Pu238"].
    '''

    particles = []
    for outgoingParticle in particleList:
        outgoingParticleId = outgoingParticle.lstrip('0123456789')
        multiplicity = outgoingParticle[:len(outgoingParticle)-len(outgoingParticleId)]
        if multiplicity == '':
            multiplicity = 1;
        multiplicity = int(multiplicity)
        particles += multiplicity * [outgoingParticleId]

    return particles


def baseAntiQualifierFromID( id, qualifierAllowed = False ) :

    if( not( isinstance( id, str ) ) ) : raise TypeError( 'id is not a str: %s' % type( id ) )
    anti, qualifier = '', ''
    base = id.split( '{' )
    if( len( base ) > 1 ) :
        if( len( base ) > 2 ) : raise ValueError( 'Invalid id = "%s"' % toLimitedString( id, 64 ) )
        base, qualifier = base
        if( qualifier[-1] != '}' ) : raise ValueError( 'Invalid qualifier "%s": id = "%s"' % \
                ( toLimitedString( qualifier, 64 ), toLimitedString( id, 64 ) ) )
        qualifier = qualifier[:-1]
    else :
        base = base[0]
    if( base[-5:] == antiSuffix ) : base, anti = base[:-5], base[-5:]

    if( not( qualifierAllowed ) and ( qualifier != '' ) ) : 
        raise ValueError( 'Particle id cannot have a qualifier: id = "%s"' % toLimitedString( id ) )

    return( base, anti, qualifier )

def baseAntiFromID( id ) :
    '''
    This function returns the base and anti parts of a PoPs GNDS id for *pid*.

    :param pid:     A PoPs GNDS id.

    :return:        A python tuple of two str instances.
    '''

    base, anti, qualifier = baseAntiQualifierFromID(id, qualifierAllowed=True)

    return( base, anti )

def idWithQualifierRemoved(pid):
    '''
    This function returns the PoPs GNDS id for *pid* with the qualifer part removed.

    :param pid:     A PoPs GNDS id.

    :return:        A python str instance.
    '''

    base, anti, qualifier = baseAntiQualifierFromID(pid, qualifierAllowed=True)

    return base + anti

def buildParticleFromRawData( cls, ID, mass = None, spin = None, parity = None, charge = None, halflife = None,
        nucleus = None, index = None, energy = None, generation = None, label = 'default' ) :

    from .quantities import quantity as quantityModule
    from .quantities import charge as chargeModule
    from .quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
    from .families import particle as particleModule
    from .families import lepton as leptonModule
    from .families import nucleus as nucleusModule
    from .families import nuclide as nuclideModule

    def getUnit( unit ) :

        if( isinstance( unit, str ) ) : unit = quantityModule.stringToPhysicalUnit( unit )
        return( unit )

    if( issubclass( cls, leptonModule.Particle ) ) :
        particle = cls( ID, generation = generation )
    elif( issubclass( cls, nucleusModule.Particle ) ) :
        ID = ID[0].lower( ) + ID[1:]
        if( index is None ) : raise ValueError( 'index must be defined for nuclide to be built' )
        particle = cls( ID, index )
        if( energy is not None ) : particle.energy.add( nuclearEnergyLevelModule.Double( label, energy[0], getUnit( energy[1] ) ) )
    elif( issubclass( cls, nuclideModule.Particle ) ) :
        particle = cls( ID )
        if( nucleus is not None ) : particle.nucleus.replicate( nucleus )
        if( charge is None ) : charge = ( 0, chargeModule.baseUnit )
        if( len( particle.nucleus.charge ) == 0 ) :
            particle.nucleus.charge.add( chargeModule.Integer( label, particle.Z, chargeModule.baseUnit ) )
        if( energy is not None ) : particle.nucleus.energy.add( nuclearEnergyLevelModule.Double( label, energy[0], getUnit( energy[1] ) ) )
    elif( issubclass( cls, particleModule.Particle ) ) :
        particle = cls( ID )
    else :
        raise TypeError( 'Invalid class.' )

    particle.buildFromRawData( mass = mass, spin = spin, parity = parity, charge = charge, halflife = halflife, label = label )

    return( particle )

def returnAntiParticleIDFromId(particleID):
    """
    This function returns the anti particle id of **particleID**. Note, if **particleID** is the id of an anti particle, 
    then the particle id is returned. For example, if **particleID** is 'O16', 'O16_anti' is returned, but if 
    **particleID** is 'O16_anti then 'O16' is returned'.
    """

    n1 = len(antiSuffix)
    if particleID[-n1:] == antiSuffix:
        return particleID[-n1:]

    return particleID + antiSuffix

def returnAntiParticleID( particle ) :

    return( returnAntiParticleIDFromId( particle.ID ) )

class ClassWithIDKey(ancestryModule.AncestryIO):   # FIXME should classes below all be declared abstract?

    __keyName = 'id'

    def __init__( self, id ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( id, str ) ) ) : raise TypeError( 'id not str' )
        self.__id = id

    @property
    def id( self ) :

        return( self.__id )

    @property
    def key( self ) :

        return( self.__id )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'id must be a string instance.' )
        self.__id = value

    @property
    def keyName( self ) :

        return( self.__keyName )

    @property
    def keyValue(self):
        '''Returns self's keyValue. If keyName is *None*, then *None* is returned.'''

        if self.keyName is None:
            return None

        return getattr(self, self.keyName)

class ClassWithSymbolKey(ancestryModule.AncestryIO):

    __keyName = 'symbol'

    def __init__( self, symbol ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( symbol, str ) ) ) : raise TypeError( 'symbol not str' )
        self.__symbol = symbol

    @property
    def symbol( self ) :

        return( self.__symbol )

    @property
    def key( self ) :

        return( self.__symbol )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'symbol must be a string instance.' )
        self.__symbol = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class ClassWithLabelKey(ancestryModule.AncestryIO):

    __keyName = 'label'

    def __init__( self, label ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label not str' )
        self.__label = label

    @property
    def label( self ) :
    
        return( self.__label )

    @label.setter
    def label( self, value ) :

        self.key = value

    @property
    def key( self ) :
    
        return( self.__label )
        
    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.__label = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class ClassWithIndexKey(ancestryModule.AncestryIO):

    __keyName = 'index'

    def __init__( self, index ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( index, str ) ) ) : raise TypeError( 'index not str' )
        self.key = index

    @property
    def index( self ) :

        return( self.__index )

    @property
    def key( self ) :

        return( self.__index )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'indexmust be a string instance.' )
        self.__index = value

    @property
    def keyName( self ) :

        return( self.__keyName )

class ClassWithSubshellKey( ancestryModule.AncestryIO):

    __keyName = 'subshell'

    def __init__( self, subshell ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( subshell, str ) ) ) : raise TypeError( 'subshell not str' )
        self.__subshell = subshell

    @property
    def subshell( self ) :

        return( self.__subshell )

    @property
    def key( self ) :

        return( self.__subshell )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'subshell must be a string instance.' )
        self.__subshell = value

    @property
    def keyName( self ) :

        return( self.__keyName )
