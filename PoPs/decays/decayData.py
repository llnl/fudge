# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines classes for storing particle decay information.
Decay data are organized as follows::

   - particle
      - decayData
         - averageEnergies
         - decayModes
            - decayMode
               - probability
               - Q
               - spectra
               - decayPath
                  - decay

All decay info is stored inside the <decayData> section.

<averageEnergies> stores a summary of average outgoing energy for various types of decay radiation,
summed over all possible decay modes.

<decayModes> contains a list of <decayMode>, each with a corresponding probability.
A <decayMode> may also contain a <Q>  (i.e. the decay Q-value), a <spectra> (describing the outgoing energy spectrum
for each types of decay product), and a <decayPath> (listing each step of the decay along with outgoing products).

Ideally the outgoing spectra would be associated with specific decay products under the <decayPath>, but
often the spectra are summed over multiple products.  For example, in B-nn decay (beta-delayed 2-neutron emission),
the spectra of the 2 outgoing neutrons cannot be resolved, so the evaluation usually stores only the summed spectrum
for the two neutrons.
"""

import abc  # FIXME remove unused import?

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule


from .. import suite as suiteModule
from .. import misc as miscModule

from ..decays import probability as probabilityModule
from ..decays import Q as QModule
from ..decays import product as productModule
from ..decays import spectrum as spectrumModule
from ..decays import averageEnergy as averageEnergyModule

decayModesIT = "isomeric transition"
decayModesSF = "spontaneous fission"
decayModesParticle = ""

decayModeTypes = [decayModesIT, decayModesSF, decayModesParticle]

class Decay( miscModule.ClassWithIndexKey ) :
    """
    Describes one step in a multi-step decay. For example, in beta-delayed 2-n emission
    the first <decay> is for the beta-decay, and the 2nd and 3rd <decay> elements describe the two neutrons
    being emitted sequentially.
    """

    moniker = 'decay'

    def __init__( self, index, mode, complete = True ) :
        """
        :param index: int, identifies the position of this decay within the decayPath. 0 = first, 1 = second, etc.
        :param mode: string identifying a decay mode, such as "isomeric transition"
        :param complete: boolean, whether all outgoing products are listed
        """

        miscModule.ClassWithIndexKey.__init__( self, index )

        if( not( isinstance( mode, str ) ) ) : raise TypeError( 'decay mode type not an instance of str' )
        if( mode not in decayModeTypes ) : raise ValueError( 'decay mode type invalid' )
        self.__mode = mode

        if( not( isinstance( complete, bool ) ) ) : raise TypeError( 'complete not an instance of bool' )
        self.__complete = complete

        self.__products = productModule.Suite( )
        self.__products.setAncestor( self )

    @property
    def complete( self ) :

        return( self.__complete )

    @property
    def products( self ) :

        return( self.__products )

    @property
    def mode( self ) :

        return( self.__mode )

    def convertUnits( self, unitMap ) :
        """See documentation in PoPs.database.convertUnits"""

        self.__products.convertUnits( unitMap )

    def copy( self ) :

        _decay = self.__class__( self.index, self.mode, self.complete )
        for item in self.__products : _decay.products.add( item.copy( ) )
        return( _decay )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        if self.mode != '':
            attrName = 'mode'
            if kwargs.get('formatVersion') == GNDS_formatVersionModule.version_1_10:
                attrName = 'type'
            attributes = ' %s="%s"' % (attrName, self.mode)
        if not self.complete: attributes = ' complete="false"'
        XMLStringList = [ '%s<%s index="%s"%s>' % ( indent, self.moniker, self.index, attributes ) ]
        XMLStringList += self.products.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        self.products.parseNode(element.find(self.products.moniker), xPath, linkData, **kwargs)

        xPath.pop()
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        complete = True
        if element.get('complete') == 'false': complete = False
        mode = element.get('mode','')
        if kwargs['formatVersion'] == GNDS_formatVersionModule.version_1_10:
            mode = element.get('type','')
        self = cls( element.get('index'), mode, complete )
        xPath.pop( )
        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )

class DecayPath( suiteModule.Suite ) :
    """
    Stores all <decay> steps involved in a decayMode
    """

    moniker = 'decayPath'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( Decay, ) )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            self.add(Decay.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop( )
        return( self )

class DecayMode( miscModule.ClassWithLabelKey ) :
    """
    Describes one way a particle can decay. The description includes the probability, Q-value, outgoing products
    and spectra. For electro-magnetic decays, may also include a set of internal conversion coefficients
    and/or photon emission probabilities.
    """

    moniker = 'decayMode'

    def __init__( self, label, mode ) :
        """
        :param label: string, must be unique within this <decayModes> section
        :param mode: string, describes the decay mode. Examples: "beta+ or e.c.", "beta-,n", "alpha", etc.
        """

        from ..decays import misc as miscDecaysModule

        miscModule.ClassWithLabelKey.__init__( self, label )

        if not isinstance(mode, str):
            if not isinstance(mode, miscDecaysModule.Mode):
                raise TypeError('mode is not a str or an enum of %s' % miscDecaysModule.Mode)
        self.__mode = mode

        self.__probability = probabilityModule.Suite( )
        self.__probability.setAncestor( self )

        self.__internalConversionCoefficients = spectrumModule.InternalConversionCoefficients()
        self.__internalConversionCoefficients.setAncestor( self )

        self.__photonEmissionProbabilities = spectrumModule.PhotonEmissionProbabilities()
        self.__photonEmissionProbabilities.setAncestor( self )

        self.__Q = QModule.Suite( )
        self.__Q.setAncestor( self )

        self.__decayPath = DecayPath( )
        self.__decayPath.setAncestor( self )

        self.__spectra = spectrumModule.Spectra( )
        self.__spectra.setAncestor( self )

    @property
    def mode( self ) :

        return( self.__mode )

    @property
    def probability( self ) :

        return( self.__probability )

    @property
    def internalConversionCoefficients( self ):

        return( self.__internalConversionCoefficients )

    @property
    def photonEmissionProbabilities( self ):

        return( self.__photonEmissionProbabilities )

    @property
    def decayPath( self ) :

        return( self.__decayPath )

    @property
    def Q( self ) :

        return( self.__Q )

    @property
    def spectra( self ) :

        return( self.__spectra )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        self.__probability.convertUnits( unitMap )
        self.__Q.convertUnits( unitMap )
        self.__decayPath.convertUnits( unitMap )
        self.__spectra.convertUnits( unitMap )

    def copy(self):
        """
        :return: deep copy of self
        """

        other = self.__class__(self.label, self.mode)
        for item in self.__probability:
            other.probability.add(item.copy())
        for item in self.__photonEmissionProbabilities:
            other.photonEmissionProbabilities.add(item.copy())
        for item in self.__Q:
            other.Q.add(item.copy())
        for item in self.__decayPath:
            other.decayPath.add(item.copy())
        for item in self.__spectra:
            other.spectra.add(item.copy())
        return other

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s" mode="%s">' % ( indent, self.moniker, self.label, self.mode ) ]
        XMLStringList += self.probability.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.internalConversionCoefficients.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.photonEmissionProbabilities.toXML_strList(indent=indent2, **kwargs)
        XMLStringList += self.Q.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.decayPath.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.spectra.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( "%s[@label='%s']" % (element.tag, element.get('label')) )

        for child in element :
            if( child.tag == probabilityModule.Suite.moniker ) :
                self.probability.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == spectrumModule.InternalConversionCoefficients.moniker ) :
                self.internalConversionCoefficients.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == spectrumModule.PhotonEmissionProbabilities.moniker ) :
                self.photonEmissionProbabilities.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == QModule.Suite.moniker ) :
                self.Q.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == DecayPath.moniker ) :
                self.decayPath.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == spectrumModule.Spectra.moniker ) :
                self.spectra.parseNode(child, xPath, linkData, **kwargs)
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        from ..decays import misc as miscDecaysModule

        xPath.append( node.tag )

        mode = node.get('mode')
        try:
            mode = miscDecaysModule.Mode.fromString(mode)
        except:
            pass

        self = cls(node.get('label'), mode)
        self.parseNode(node, xPath, linkData, **kwargs)

        xPath.pop( )
        return( self )

class DecayModes( suiteModule.Suite ) :
    """
    Contains a list of decayMode instances
    """

    moniker = 'decayModes'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( DecayMode, ) )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            self.add(DecayMode.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop( )
        return( self )

class DecayData(ancestryModule.AncestryIO_bare):
    """
    Contains all decay information for a particle, including average energies for decay products
    and a list of decay modes.
    """

    moniker = 'decayData'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__decayModes = DecayModes( )
        self.__decayModes.setAncestor( self )

        self.__averageEnergies = averageEnergyModule.AverageEnergies( )
        self.__averageEnergies.setAncestor( self )

    @property
    def decayModes( self ) :

        return( self.__decayModes )

    @property
    def averageEnergies( self ) :

        return( self.__averageEnergies )

    def convertUnits( self, unitMap ) :
        """See convertUnits documentation in PoPs.database"""

        self.__decayModes.convertUnits( unitMap )
        self.__averageEnergies.convertUnits( unitMap )

    def check(self, info):
        """Check that decay-mode branching ratios sum to one."""

        from .. import warning as warningModule

        if len(self.__decayModes) == 0:
            return []

        branchingSumTolerance = info.get('branchingRatioSumTolerance', 1e-5)
        branchingSum = sum([decayMode.probability[0].value for decayMode in self.__decayModes])
        if abs(branchingSum - 1) > branchingSumTolerance:
            return [warningModule.UnnormalizedDecayProbabilities(branchingSum, self)]

        return []

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _decayData = DecayData( )
        self.copyItems( _decayData )
        return( _decayData )

    def copyItems( self, other ) :
        """
        Copy all items in self to other
        :param other: decayData instance where contents of self will be copied
        """

        for item in self.__decayModes : other.decayModes.add( item.copy( ) )
        for item in self.__averageEnergies : other.averageEnergies.add( item.copy( ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList_subs  = self.__decayModes.toXML_strList( indent = indent2, **kwargs )
        XMLStringList_subs += self.averageEnergies.toXML_strList( indent = indent2, **kwargs )
        if( len( XMLStringList_subs ) == 0 ) : return( [] )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ] + XMLStringList_subs
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        children = { 'decayModes' : self.decayModes, 'averageEnergies' : self.averageEnergies }

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseNode(child, xPath, linkData, **kwargs)
            else:
                if not self.parseExtraXMLElement(child, xPath, linkData, **kwargs): raise ValueError('Sub-element = "%s" not allowed' % child.tag)

        xPath.pop()

        return self
