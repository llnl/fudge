# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the PoPs database class, which is the top-level container in PoPs.
A database is a collection of particles, each with a unique id.
It may also store documentation and 'style' information.

Within the database, particles are organized by family (i.e. gauge bosons, leptons, baryons or chemical elements).
For nuclear data users, the largest part of the database will likely be inside the chemical elements section
where isotopes, nuclides and excited nuclear states are defined.

The database supports looking up particles by their id. For example, if 'pops' is a PoPs.database instance that
contains a particle with id 'Al27',

>>> 'Al27' in pops
True
>>> type( pops['Al27'] )
<class PoPs.families.nuclide.particle>
>>> pops.keys()   # list all particle ids in the database
['photon', 'n', ... 'Al28', 'Al26_m1', 'Na24_m1']

The particle database may also define aliases. These can be used to refer to a particle by another name.
For example, the id for the electron is normally "e-", but could also be aliased to 'e' or 'electron'.
In the example above, 'Al26_m1' and 'Na24_m1' are aliases pointing to excited nuclides.

>>> type( pops['Na24_m1'] )
<class 'PoPs.alias.MetaStable'>

Each alias has a 'pid' field with the id of the actual particle it points to:

>>> pops['Na24_m1'].pid
'Na24_e1'

In this case, the first metastable state of Na24 is also the first excited state (as indicated by the '_e1' suffix).

A PoPs database has methods to write itself out to an XML-formatted file or string, to read back in from XML,
to convert units, etc.

One PoPs database may point to another parent database. In that case, if a particle is not found in the current
database it will also check in the parent database before throwing a KeyError.
This option allows an evaluator to override part of the parent database (i.e. to define a different mass
for one or more particles).
"""

import os

from LUPY import ancestry as ancestryModule

from xData.Documentation import documentation as documentationModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from .chemicalElements import chemicalElements as chemicalElementsModule
from .chemicalElements import chemicalElement as chemicalElementModule
from .chemicalElements import isotope as isotopeModule
from .chemicalElements import misc as chemicalElementsMiscModule
from .families import particle as particleModule
from .families import gaugeBoson as gaugeBosonModule
from .families import lepton as leptonModule
from .families import baryon as baryonModule
from .families import nuclide as nuclideModule
from .families import nucleus as nucleusModule
from .families import unorthodox as unorthodoxModule

from . import IDs as IDsModule
from . import misc as miscModule
from . import styles as stylesModule
from . import alias as aliasModule

class Database(ancestryModule.AncestryIO):

    moniker = 'PoPs'
    ancestryMembers = ( 'styles', 'documentation', 'aliases', 'gaugeBosons', 'leptons', 'baryons', 'chemicalElements', 'unorthodoxes' )

    def __init__( self, name, version, formatVersion = GNDS_formatVersionModule.default ) :
        """
        :param name: string naming the database
        :param version: database version
        :param formatVersion: the GNDS format the data are represented in.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__format = formatVersion

        self.name = name

        self.version = version

        self.__styles = stylesModule.Styles( )
        self.__styles.setAncestor( self )

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

        self.__aliases = aliasModule.Suite( )
        self.__aliases.setAncestor( self )

        self.__gaugeBosons = gaugeBosonModule.Suite( )
        self.__gaugeBosons.setAncestor( self )

        self.__leptons = leptonModule.Suite( )
        self.__leptons.setAncestor( self )

        self.__baryons = baryonModule.Suite( )
        self.__baryons.setAncestor( self )

        self.__chemicalElements = chemicalElementsModule.ChemicalElements( )
        self.__chemicalElements.setAncestor( self )

        self.__unorthodoxes = unorthodoxModule.Suite( )
        self.__unorthodoxes.setAncestor( self )

    def __contains__( self, key ) :

        if( key in self.__aliases ) : return( True )
        if( key in self.__gaugeBosons ) : return( True )
        if( key in self.__leptons ) : return( True )
        if( key in self.__baryons ) : return( True )
        if( key in self.__chemicalElements ) : return( True )
        if( key in self.__unorthodoxes ) : return( True )
        parentDB = self.parentDatabase()
        if parentDB is not None:
            return key in parentDB
        return( False )

    def __getitem__( self, key ) :
        """
        :param key: id or symbol of a particle or particle group

        :return: particle with given key. KeyError is raised if the key is not found.
        """

        if( '{' in key ) : key = key.split( '{' )[0]
        return( self.__getitem( key, 0 ) )

    def __getitem(self, key, level):

        if not isinstance(key, str):
            raise TypeError('key must be a str')
        if key in self.__aliases:
            return self.__aliases[key]
        for suite in [self.__gaugeBosons, self.__leptons, self.__baryons, self.__chemicalElements, self.__unorthodoxes]:
            try:
                return suite[key]
            except KeyError:
                continue
            except:
                raise
        parentDB = self.parentDatabase()
        if parentDB is not None:
            return parentDB[key]
        raise KeyError('id = "%s" not found' % key)

    def __iter__( self ) :
        """
        Iterates over all particles in the database. Particle groups like the
        chemicalElement and isotope classes are not included, although the particles they contain are.
        """

        for gaugeBoson in self.__gaugeBosons : yield gaugeBoson
        for lepton in self.__leptons : yield lepton
        for baryon in self.__baryons : yield baryon
        for chemicalElement in self.__chemicalElements :
            for isotope in chemicalElement :
                for nuclide in isotope :
                    yield nuclide
                    yield nuclide.nucleus
        for item in self.__unorthodoxes : yield item
        for item in self.__aliases : yield item
        parentDB = self.parentDatabase()
        if parentDB is not None:
            for item in parentDB: yield item    # FIXME skip particles that are overridden in child database?

    def keys( self ) :
        """Return a list of ids for all particles in self."""

        return( [ particle.id for particle in self ] )

    @property
    def format( self ) :

        return( self.__format )

    @property
    def name( self ) :

        return( self.__name )

    @name.setter
    def name( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'name must be a string' )
        self.__name = value

    @property
    def version( self ) :

        return( self.__version )

    @version.setter
    def version( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = value

    @property
    def styles( self ) :

        return( self.__styles )

    @property
    def gaugeBosons( self ) :

        return( self.__gaugeBosons )

    @property
    def leptons( self ) :

        return( self.__leptons )

    @property
    def baryons( self ) :

        return( self.__baryons )

    @property
    def chemicalElements( self ) :

        return( self.__chemicalElements )

    @property
    def unorthodoxes( self ) :

        return( self.__unorthodoxes )

    @property
    def aliases( self ) :

        return( self.__aliases )

    def parentDatabase(self):
        """
        Search up ancestry chain for another PoPs database.
        Return PoPs.database if found, otherwise return None.
        """
        parentDB = None
        try:
            grandparent = self.ancestor.ancestor
            parentDB = grandparent.findAttributeInAncestry(Database.moniker)
        except Exception as ex:
            # no parent database found
            pass

        return parentDB

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    def add( self, particle ) :
        """
        Add a particle to the database, inside the appropriate particle family.

        :param particle: instance of families.particle.Particle or derived class
        """

        if( isinstance( particle, gaugeBosonModule.Particle ) ) :
            self.__gaugeBosons.add( particle )
        elif( isinstance( particle, leptonModule.Particle ) ) :
            self.__leptons.add( particle )
        elif( isinstance( particle, baryonModule.Particle ) ) :
            self.__baryons.add( particle )
        elif( isinstance( particle, ( chemicalElementModule.ChemicalElement, isotopeModule.Isotope, nuclideModule.Particle ) ) ) :
            self.__chemicalElements.add( particle )
        elif( isinstance( particle, unorthodoxModule.Particle ) ) :
            self.__unorthodoxes.add( particle )
        elif( isinstance(particle, aliasModule.BaseAlias)) :
            if( particle.id not in self.__aliases ) :      # Do not allow an alias id to overwrite a particle id.
                if( particle.id in self ) : raise ValueError( 'particle with id = "%s" already in database' % particle.id )
            self.__aliases.add( particle )
        else :
            raise TypeError( 'Unsupported particle family = "%s"' % particle.moniker )

    def addFile( self, fileName, replace = False ) :
        """
        Opens the file *fileName* and adds its particles to self. If replace is False, particles already in self are
        not replaced.
        """

        pops = self.read(fileName)
        for particle in pops :
            if( not( replace ) and particle.id in self ) : continue
            if( isinstance( particle, nucleusModule.Particle ) ) : continue
            self.add( particle )

    def check( self, **kwargs ) :
        """
        Check for physics problems in the database

        :param kwargs: options controlling i.e. how strict checks should be. Currently unused.
        :return: list of warnings along with context about where they appear in the database
        """

        from . import warning as warningModule

        warnings = []
        info = kwargs
        info['PoPs'] = self
        for suite in ( self.__gaugeBosons, self.__leptons, self.__baryons, self.__chemicalElements, self.__unorthodoxes, self.__aliases ) :
            suiteWarnings = suite.check( info )
            if suiteWarnings:
                warnings.append( warningModule.Context("%s" % suite.moniker, suiteWarnings) )

        result = warningModule.Context('PoPs', warnings)
        result.info = info
        return result

    def convertUnits( self, unitMap ) :
        """
        Recursively searches for units within all sections, converting units that appear in the unitMap.

        :param unitMap: dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        self.__styles.convertUnits( unitMap )
        self.__gaugeBosons.convertUnits( unitMap )
        self.__leptons.convertUnits( unitMap )
        self.__baryons.convertUnits( unitMap )
        self.__chemicalElements.convertUnits( unitMap )
        self.__unorthodoxes.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of the entire database.
        """

        _database = Database( self.name, self.version, formatVersion = self.format )
        for item in self.__styles : _database.styles.add( item.copy( ) )
            # FIXME, need to add documentation.
        for item in self.__aliases : _database.add( item.copy( ) )
        for item in self.__gaugeBosons : _database.add( item.copy( ) )
        for item in self.__leptons : _database.add( item.copy( ) )
        for item in self.__baryons : _database.add( item.copy( ) )
        for item in self.__chemicalElements : _database.add( item.copy( ) )
        for item in self.__unorthodoxes : _database.add( item.copy( ) )

        return( _database )

    def estimateMass(self, pid, unit, maxASpan=4):
        '''
        Returns the mass of the particle specified by *pid* if it exists in *self* and has a mass; otherwise,
        the mass is estimate if possible. Currenlty, estimates are only performed if *pid* is a nuclide. In 
        which case, its mass is estimated by trying to find a near by nuclide with the same Z. If this fails, 
        None is returned.

        If the particle exists and has a mass (ergo, the mass is not estimated), the returned value is positive. 
        If the mass is estimated, a negative value is returned whose absolute value is the estimated mass. 
        Otherwise, None is returned.

        If the mass is estimated for a nuclide, the returned mass is for the ground state of the nuclide. Ergo, the *pid*
        specifies a nuclear level, the level information is ignored.

        :param pid:             The PoPs id of the particle whose mass is to be returned.
        :param unit:            The unit of the returned mass.
        :maxASpan:              If estimating a nuclide mass, only consider A's within this value of the A for *pid*.

        :return:                A positive float if a mass for *pid* exist in *self*, a negative float if the mass is estimated and None otherwise.
        '''

        mass = None

        try:
            particle = self[pid]
            try:
                return particle.getMass(unit)
            except:
                pass
        except:
            pass

        dummy, symbol, A, level, dummy, dummy, dummy = chemicalElementsMiscModule.chemicalElementALevelIDsAndAnti(pid)
        if symbol is None:
            return None
        try:
            chemicalElement = self.chemicalElements[symbol]
        except:
            return None

        lowerIsotopes = [999, None]
        upperIsotopes = [999, None]
        neutronMass = self[IDsModule.neutron].getMass(unit)
        for index, isotope in enumerate(chemicalElement):
            A_diff = A - isotope.A
            if A_diff > 0:
                if A_diff < lowerIsotopes[0]:
                    lowerIsotopes = [A_diff, isotope]
            else:
                A_diff = abs(A_diff)
                if A_diff < upperIsotopes[0]:
                    upperIsotopes = [A_diff, isotope]

        if lowerIsotopes[0] < maxASpan:
            isotope = lowerIsotopes[1]
            mass = isotope[isotope.symbol].getMass(unit) + neutronMass * lowerIsotopes[0]
        if upperIsotopes[0] < maxASpan:
            numberOfMass = 2
            if mass is None:
                mass = 0.0
                numberOfMass = 1
            isotope = upperIsotopes[1]
            mass = (isotope[isotope.symbol].getMass(unit) - neutronMass * upperIsotopes[0] + mass) / numberOfMass
        if mass is not None:
            mass = -mass

        return mass

    def final(self, id, returnAtMetaStableAlias=False):
        """
        Returns the particle matching id. If id is an alias, will follow the alias pid until a non-alias particle is found and that particle is returned.
        """

        particle = self[id]
        if isinstance(particle, aliasModule.BaseAlias):
            if not isinstance(particle, aliasModule.MetaStable) or not returnAtMetaStableAlias:
                particle = self.final(particle.pid, returnAtMetaStableAlias)

        return particle

    def hasAlias( self, id ) :

        if id in self.__aliases:
            return True
        parentDB = self.parentDatabase()
        if parentDB is not None:
            return parentDB.hasAlias(id)
        return False

    def QValue(self, ingoingParticles, outgoingParticles, unit='MeV'):
        '''
        Calculates the Q value from the list of ingoing and outgoing particle ids.

        :param ingoingParticles:    List of ingoing particle ids.
        :param outgoingParticles:   List of outgoing particle ids.
        :param unit:                Unit of returned Q-value.

        :return:                    Q-value.
        '''

        Q = 0.0

        ingoingParticles = list(ingoingParticles)
        ingoingParticles = miscModule.parseParticleList(ingoingParticles)
        outgoingParticles = miscModule.parseParticleList(outgoingParticles)

        removeDuplicates = True
        while(removeDuplicates):
            removeDuplicates = False
            for particle in ingoingParticles:
                if particle in outgoingParticles:
                    removeDuplicates = True
                    break
            if removeDuplicates:
                ingoingParticles.remove(particle)
                outgoingParticles.remove(particle)

        for particle in ingoingParticles:
            Q += self[particle].getMass(unit + '/c**2')

        for particle in outgoingParticles:
            Q -= self[particle].getMass(unit + '/c**2')

        return Q

    def saveToFile( self, fileName, **kwargs ) :

        dirname = os.path.dirname( fileName )
        if( dirname != '' ) :
            if( not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        fOut = open( fileName, 'w' )
        fOut.write( ancestryModule.XML_declaration )
        fOut.write( self.toXML( **kwargs ) )
        fOut.write( '\n' )
        fOut.close( )

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        if formatVersion in (GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4):
            print('INFO: converting GNDS format from "%s" to "%s".' % (formatVersion, GNDS_formatVersionModule.version_2_0))
            formatVersion = GNDS_formatVersionModule.version_2_0
        if formatVersion not in GNDS_formatVersionModule.allowed:
            raise Exception("Unsupported GNDS structure '%s'!" % str(formatVersion))
        kwargs['formatVersion'] = formatVersion

        xmlString = ['%s<%s name="%s" version="%s" format="%s">' % (indent, self.moniker, self.name, self.version, formatVersion)]
        # FIXME for v1.10 also print link to parent database?
        xmlString += self.__styles.toXML_strList(indent2, **kwargs)

        if formatVersion == '0.1' and len(self.styles) != 0:
            xmlString += self.styles.getEvaluatedStyle().documentation.toXML_strList(indent2, **kwargs)
        # FIXME: currently both PoPs and PoPs/styles/evaluated have documentation nodes according to GNDS-2.0 specs
        xmlString += self.__documentation.toXML_strList(indent2, **kwargs)

        xmlString += self.__aliases.toXML_strList(indent2, **kwargs)
        xmlString += self.__gaugeBosons.toXML_strList(indent2, **kwargs)
        xmlString += self.__leptons.toXML_strList(indent2, **kwargs)
        xmlString += self.__baryons.toXML_strList(indent2, **kwargs)
        xmlString += self.__chemicalElements.toXML_strList(indent2, **kwargs)
        xmlString += self.__unorthodoxes.toXML_strList(indent2, **kwargs)

        if len(xmlString) == 1:
            xmlString = []
        else:
            xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        formatVersion = element.get( 'format' )
        if( formatVersion == '0.1' ) :                             # PoPs 0.1 and 2.0 are the same so allow for now. There was no 1.10 for PoPs.
            formatVersion = GNDS_formatVersionModule.version_1_10
        else :
            if( formatVersion not in GNDS_formatVersionModule.allowedPlus ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )

        self.__format = formatVersion
        kwargs['formatVersion'] = formatVersion

        children = {}
        for child in ( gaugeBosonModule, leptonModule, baryonModule, unorthodoxModule, aliasModule ) :
            children[child.Suite.moniker] = child.Suite

        children[chemicalElementsModule.ChemicalElements.moniker] = chemicalElementsModule.ChemicalElements

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseNode(getattr( self, child.tag ), child, xPath, linkData, **kwargs)
            elif( child.tag == self.__styles.moniker ) :
                self.__styles.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == 'parentDatabase' and formatVersion in (GNDS_formatVersionModule.version_1_10,
                        GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4) ) :
                pass  # ignore
            elif( child.tag == self.documentation.moniker ) :
                if formatVersion == GNDS_formatVersionModule.version_1_10 :
                    self.styles.getEvaluatedStyle().documentation.parseNode(child, xPath, linkData, **kwargs)
                else:
                    self.documentation.parseNode(child, xPath, linkData, **kwargs)
            else :
                raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )

        return self

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        formatVersion = element.get( 'format' )
        if formatVersion in [GNDS_formatVersionModule.version_2_0_LLNL_3, GNDS_formatVersionModule.version_2_0_LLNL_4]:
            formatVersion = GNDS_formatVersionModule.version_2_0
        if( formatVersion == '0.1' ) :                             # PoPs 0.1 and 2.0 are the same so allow for now. There was no 1.10 for PoPs.
            formatVersion = GNDS_formatVersionModule.version_1_10
        else :
            if( formatVersion not in GNDS_formatVersionModule.allowedPlus ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )

        self = cls( element.get( 'name' ), element.get( 'version' ), formatVersion = formatVersion )
        self.parseNode(element, xPath, linkData, **kwargs)
        return( self )

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **database** instance.
        """

        return Database.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    """
    Reads in the file name *fileName* and returns a **database** instance.
    """

    return Database.read(fileName)
