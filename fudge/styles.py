# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import datetime
import textwrap
import shutil

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import axes as axesModule
from xData.Documentation import documentation as documentationModule

from fudge import physicalQuantity as physicalQuantityModule
from fudge import targetInfo as targetInfoModule

from fudge.processing import flux as fluxModule
from fudge.processing import transportables as transportablesModule
from fudge.processing import inverseSpeed as inverseSpeedModule


def getClassDict():
    """Determine supported styles through module introspection"""

    classDict = {}
    from . import styles as tmpStyles
    for thing in [ x for x in dir( tmpStyles ) if not x.startswith('__')  ] :
        try:
            thingClass = tmpStyles.__dict__[thing]
            if issubclass(thingClass, Style) and thingClass != Style:
                classDict[thingClass.moniker]=thingClass
        except TypeError: pass # catches imported modules that are irrelevant
    return classDict

class TemperatureInfo:
    """
    For a given temperature, stores the labels for each type of processed style. The list of processed styles represented 
    are Heated, GriddedCrossSection, URR_probabilityTables, HeatedMultiGroup and SnElasticUpScatter. If a style is not 
    present, its label will be an empty string.
    """

    def __init__(self, temperature, heated, griddedCrossSection, URR_probabilityTables, heatedMultiGroup, SnElasticUpScatter):
        """
        :temperature:               The temperature of the data represented by the labels.
        :heated:                    The label for the Heated.
        :griddedCrossSection:       The label for the GriddedCrossSection.
        :URR_probabilityTables:     The label for the URR_probabilityTables.
        :heatedMultiGroup:          The label for the HeatedMultiGroup.
        :SnElasticUpScatter:        The label for the SnElasticUpScatter.
        """

        self.temperature = temperature
        self.heated = heated
        self.griddedCrossSection = griddedCrossSection
        self.URR_probabilityTables = URR_probabilityTables
        self.heatedMultiGroup = heatedMultiGroup
        self.SnElasticUpScatter = SnElasticUpScatter

    def info(self):
        '''
        This method prints each member of *self* and information about its.
        '''

        width = shutil.get_terminal_size().columns
        subsequent_indent = len('    URR_probabilityTables: ') * ' '

        def wrap(message):
            '''
            For internal use only. This function wraps *message* to fit onto the screen.
            '''

            for line in textwrap.wrap(message, width=width, subsequent_indent=subsequent_indent):
                print(line)

        print()
        wrap('The TemperatureInfo class members:')
        wrap('    temperature:           The temperature associated with the following labels.')
        wrap('    heated:                The label for the cross sections heated to the specified temperature.')
        wrap('    griddedCrossSection:   The label for the cross sections mapped onto a unified grid for all reactions. These data are for Monte Carlo transport.')
        wrap('    URR_probabilityTables: The label for the cross section URR probability tables.')
        wrap('    heatedMultiGroup:      The label for the multi-group data.')
        wrap('    SnElasticUpScatter:    The label for the elastic multi-group average-product-energy and transfer-matrices with upscatter data.')

class Styles(ancestryModule.AncestryIO_base):
    """
    Stores the list of nuclear data styles that appear inside a file.

    The list generally includes one 'evaluated' style, plus 0 or more derived styles (Heated, HeatedMultiGroup, etc.).
    """

    moniker = 'styles'
    """Label for use as the tag or identifier for this style"""

    def __init__( self ) :

        self.__styles = []
        ancestryModule.AncestryIO_base.__init__(self)

    def __contains__( self, label ) :

        for item in self :
            if( item.label == label ) : return( True )
        return( False )

    def __len__( self ) :

        return( len( self.__styles ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__styles[i1]

    def __getitem__( self, label ) :

        if( isinstance( label, int ) ) : return( self.__styles[label] )
        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        for _style in self.__styles :
            if( _style.label == label ) : return( _style )
        raise IndexError( 'No style labelled == "%s"' % label )

    def add( self, _style ) :
        """
        Append a style to the list of styles.

        @param _style: style instance
        """

        if( not( isinstance( _style, Style ) ) ) : raise TypeError( 'invalid style instance' )
        for __style in self :
            if( __style.label == _style.label ) : raise ValueError( 'style labeled "%s" already exists' % _style.label )

        self.__styles.append( _style )
        _style.setAncestor(self)

    def keys(self):
        """ Returns a list of all style labels. """
        return [style.label for style in self]

    def chains( self, ends = False, _styles = None ) :
        """For self, determine the derived from chains for each style. If ends is True, trains which are a part of a longer train are removed from the list."""

        if( _styles is None ) :
            _chains = [ _style.chain( ) for _style in self ]
        else :
            _chains = [ _style.chain( ) for _style in _styles ]

        if( ends ) :
            removeIndices = set( )
            for i1 in range( len( _chains ) ) :
                style1 = _chains[i1][0]
                for i2 in range( len( _chains ) ) :
                    if( style1 in _chains[i2][1:] ) : removeIndices.add( i1 )
            for index in reversed( sorted( removeIndices ) ) : _chains.pop( index )

        return( _chains )

    def chainWithLabel( self, label ) :

        chains = self.chains( ends = True )
        for chain in chains :
            for style in chain :
                if( style.label == label ) : return( chain )
        raise ValueError( 'No style with label = "%s" found.' % label )

    def preProcessingHeadInChainWithLabel(self, label, include_reconstructed=True):
        """
        Returns the pre-processing head of the chain containing the label *label*.
        """

        chain = self.chainWithLabel( label )
        preProcessingStyles = self.preProcessingStyles(include_reconstructed)
        for style in chain :
            if( style.__class__ in preProcessingStyles ) : return( style )
        raise ValueError( 'No pre-processing style in chain with label = "%s" found.' % label )

    def preProcessingOnly(self, include_reconstructed=True) :
        """
        Returns True if self only contains pre-processed styles and False otherwise.

        @param include_reconstructed: whether to count reconstructed resonances as pre-processed styles
        @type include_reconstructed: bool
        """

        preProcessingStyles = self.preProcessingStyles(include_reconstructed)
        clean = True                            # True means no processed data.
        for style in self :
            if( not( isinstance( style, preProcessingStyles ) ) ) :
                clean = False
                break

        return( clean )

    def convertUnits( self, unitMap ) :

        for style in self:
            style.convertUnits(unitMap)

    def preProcessingChainHead(self, label=None):
        """
        Returns the pre-processing style that heads the chain containing the pre-processing label *label*. If label is None, there must be
        only one pre-processing chain and that head is returned.
        """

        chains = self.preProcessingChains(ends=True)
        if len(chains) == 0:
            raise Exception('No pre-processing style exist.')
        if label is None:
            if len(chains) > 1:
                labels = []
                for chain in chains:
                    labels.append([style.label for style in reversed(chain)])
                raise Exception('Must have only one pre-processing chain, have %s: %s.' % (len(chains), labels))
            return chains[0][0]
        else:
            for chain in chains:
                for style in chain:
                    if style.label == label:
                        return chain[0]
            raise Exception('No pre-processing chain containing label "%s" found.' % label)

    def preProcessingChains(self, ends=True):
        """
        Returns chains for each style that can be processed (i.e., Evaluated, CrossSectionReconstructed, 
        AngularDistributionReconstructed and Realization).
        """

        preProcessingStyles = self.preProcessingStyles(include_reconstructed=True)
        preProcessingStyleInstances = []
        for style in self :
            if isinstance(style, preProcessingStyles): preProcessingStyleInstances.append(style)

        return self.chains(ends=ends, _styles=preProcessingStyleInstances)

    def projectileEnergyDomain(self):
        """
        Returns the projectileEnergyDomain for the evaluated style.
        """

        return self.preProcessingChainHead().projectileEnergyDomain

    def findEntity(self, entityName, attribute = None, value = None):

        for style in self:
            if style.label == value: return style

        raise ValueError('Could not find style "%s" with label "%s".' % (entityName, value))

    def findDocumentationEntriesWithKey(self, keyValue):
        '''
        Calls **findDocumentationEntriesWithKey** on each style in *self* and returns the results.
        '''

        entries = []
        for style in self:
            entries += style.findDocumentationEntriesWithKey(keyValue)

        return entries

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for item in self :
            if isinstance(item, cls): foundInstances.append(item)
            foundInstances += item.findInstancesOfClassInChildren( cls, level )

        return( foundInstances )

    def fixDomains(self, energyMax):
        """
        For all **evaluated* styles, calls their **fixDomains** method.
        """

        numberOfFixes = 0
        labels = []
        for entry in self:
            if hasattr(entry, 'fixDomains'):
                labels.append(entry.label)
                numberOfFixes += entry.fixDomains(energyMax)

        return(numberOfFixes, labels)

    def remove( self, _style ):
        """
        Remove the specified style. Accepts either a string or a style instance
        """

        index = None
        for idx,tmp in enumerate(self):
            if( ( tmp is _style ) or ( tmp.label == _style ) ) : index = idx
        if index is None:
            raise KeyError("style '%s' not found in styles" % _style)
        self.__styles.pop(index)

    def getStylesOfClass( self, cls ) :
        """
        Returns a list of all styles of class cls.
        """

        styleList = []
        for _style in self :
            if( isinstance( _style, cls ) ) : styleList.append( _style )
        return( styleList )

    def getStyleOfClass( self, cls ) :
        """
        Returns the style of class cls. If no style is found, None is returned. If more than 1 style is found, an exception is raised.
        """

        styleList = self.getStylesOfClass( cls )
        if( len( styleList ) == 0 ) : return( None )
        if( len( styleList )  > 1 ) : raise KeyError( '''multiple (%d) style's of type "%s" found''' % ( len( styleList ), cls.moniker ) )
        return( styleList[0] )

    def getTempStyleNameOfClass( self, cls ):

        styleNames = [ _style.label for _style in self ]
        index = 0
        while( True ) :
            tmpLabel = 'tmp%d_%s' % ( index, cls.moniker )
            if( tmpLabel not in styleNames ) : return( tmpLabel )
            index += 1

    def getEvaluatedStyle( self ) :

        _evaluateds = self.getStylesOfClass( Evaluated )
        _styles = self.chains( ends = True, _styles = _evaluateds )
        if( len( _styles ) == 0 ) : raise Exception( '''No evaluated style found.''' )
        if( len( _styles ) > 1 ) : raise Exception( '''Multiple (%d) evaluated styles.''' % ( len( _styles ) ) )

        return( self[_styles[0][0].label] )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        for label in styleLabels :
            for index, style in enumerate( self.__styles ) :
                if( style.label == label ) :
                    del self.__styles[index]
                    break

    def temperatures( self, unit = 'MeV/k' ) :
        """
        Returns a listed of [ temperature, dictionary ] instances sorted by temperature. The directory keys are the
        monikers for the heated process styles. The associated values are the labels for that heated style. If the
        label is an empty string (i.e., ''), that moniker has not data.
        """

        heatedProcessedStyles = self.heatedProcessedStyles( )
        temperatureInfos = {}
        for instance in self :
            for style in heatedProcessedStyles :
                if( instance.moniker == style.moniker ) :
                    temperature = instance.temperature.getValueAs( unit )
                    if( temperature not in temperatureInfos ) :
                        temperatureInfos[temperature] = {}
                        for style2 in heatedProcessedStyles : temperatureInfos[temperature][style2.moniker] = ''
                    temperatureInfos[temperature][style.moniker] = instance.label
                    break

        temperatureInfos2 = []
        for temperature in sorted(temperatureInfos):
            temperatureInfo = temperatureInfos[temperature]
            temperatureInfos2.append(TemperatureInfo(temperature, temperatureInfo['heated'], temperatureInfo['griddedCrossSection'], 
                temperatureInfo['URR_probabilityTables'], temperatureInfo['heatedMultiGroup'], temperatureInfo['SnElasticUpScatter']))

        return temperatureInfos2

    def toXML_strList(self, indent = '', **kwargs):
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        do_multiGroup = kwargs['formatVersion'] == GNDS_formatVersionModule.version_1_10

        parameters = ''
        if do_multiGroup:
            heatedMultiGroups = [_style for _style in self if isinstance(_style, HeatedMultiGroup)]
            do_multiGroup = len(heatedMultiGroups) > 0
            if do_multiGroup:
                parameters = 'MultiGroup'
                multiGroup1 = MultiGroup(parameters, 3)
                for transportable in heatedMultiGroups[0].transportables:
                    multiGroup1.transportables.add(transportable)

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for _style in self:
            if do_multiGroup and isinstance(_style, Heated):
                xmlStringList += multiGroup1.toXML_strList(indent2, **kwargs)
                do_multiGroup = False
            xmlStringList += _style.toXML_strList(indent2, **kwargs, parameters = parameters)

        if len(xmlStringList) == 1:
            if kwargs.get('showEmpty', False) or kwargs.get('showEmptySuite', False):
                xmlStringList = [xmlStringList[0][:-1] + '/>']
            else:
                xmlStringList = []
        else:
            xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    def parseNode(self, stylesElement, xPath, linkData, **kwargs):
        """
        Generate instances of allowed classes from the available XML nodes with the corresponding styles.
        """

        xPath.append( stylesElement.tag )

        classDict = {}
        for _style in ( Evaluated, CrossSectionReconstructed, AngularDistributionReconstructed, Realization,
                        AverageProductData, MonteCarlo_cdf, Heated, HeatedMultiGroup,
                        GriddedCrossSection, SnElasticUpScatter, CoulombPlusNuclearElasticMuCutoff, URR_probabilityTables ) :
            classDict[_style.moniker] = _style

        multiGroupChildren = []
        for child in stylesElement :
            if( child.tag == MultiGroup.moniker ) :       # Only in GNDS 1.10.
                multiGroupChildren.append( child )
            else :
                _class = classDict.get( child.tag, None )
                if( _class is None ) : raise TypeError( 'Encountered unknown style "%s".' % child.tag )
                self.add(_class.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        for child in multiGroupChildren :
            multiGroup1 = MultiGroup.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            for style1 in self :
                if( isinstance( style1, HeatedMultiGroup ) ) :
                    if( multiGroup1.label == style1.parameters ) :
                        for transportable in multiGroup1.transportables : style1.transportables.add( transportable )

        xPath.pop()

    @staticmethod
    def preProcessingStyles(include_reconstructed):
        """
        Returns the list of style classes that can be processed. These are known as pre-processed styles.

        @param include_reconstructed: whether to include CrossSectionReconstructed and AngularDistributionReconstructed
        @type include_reconstructed: bool
        """

        preprocessed_styles = [Evaluated, Realization]
        if include_reconstructed:
            preprocessed_styles += [CrossSectionReconstructed, AngularDistributionReconstructed]
        return tuple(preprocessed_styles)

    @staticmethod
    def heatedProcessedStyles():

        return (Heated, GriddedCrossSection, URR_probabilityTables, HeatedMultiGroup, SnElasticUpScatter)

class Style(ancestryModule.AncestryIO, abc.ABC):
    """
    Abstract base class for the classes that represent actual GNDS style nodes.
    """

    keyName = 'label'

    def __init__( self, label, derivedFrom, date = None ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        if( date is None ) : date = str( datetime.date.today( ) )
        self.__date = date

        if( not( isinstance( derivedFrom, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__derivedFrom = derivedFrom

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

    @property
    def date( self ) :
        """
        Property for the style release date.
        """

        return( self.__date )

    @property
    def derivedFrom( self ) :
        """
        Property indicating the parent style label from which the data in style is derived.
        """

        return( self.__derivedFrom )

    @property
    def derivedFromStyle( self ) :
        """
        Property for the style (not just the label) from which the data in this style is derived.
        """

        for _style in self.ancestor :
            if( _style.label == self.__derivedFrom ) : return( _style )
        return( None )

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    @property
    def label( self ) :
        """
        Property for the string identifier of this style.
        """

        return( self.__label )

    @property
    def projectileEnergyDomain( self ) :
        """
        Property for the projectileEnergyDomain at which the data in this style is evaluated.
        """

        return( self.derivedFromStyle.projectileEnergyDomain )

    @property
    def temperature( self ) :
        """
        Property for the temperature at which the data in this style is evaluated.
        """

        return( self.derivedFromStyle.temperature )

    @property
    @abc.abstractmethod
    def sortOrderIndex( self ) :

        pass

    def chain( self ) :
        """Returns the list of derivedFrom style instances for self, ordered from self to root."""

        _chain = [ self ]
        while( True ) :
            _derivedFrom = _chain[-1].derivedFromStyle
            if( _derivedFrom is None ) : break
            _chain.append( _derivedFrom )

        return( _chain )

    def findDocumentationEntriesWithKey(self, keyValue):
        '''
        Calls **findEntriesWithKey** the *self*'s documentation member and returns a list of [self, entry] for each entry in *self*'s documentation
        that **findEntriesWithKey** returns.
        '''

        entries = self.__documentation.findEntriesWithKey(keyValue)
        entries = [[self, entry] for entry in entries]

        return entries

    def findFormMatchingDerivedStyle( self, component, styleFilter = None ) :
        """
        This method searches the link of derivedFroms, starting with self's derivedFrom,
        to find a form in component matching one of the derivedFroms. If a form is found
        matching one of the derivedFroms, that form is returned. If no match is found, None is returned.
        """

        def alwaysTrue( a_style ) : return( True )

        if( styleFilter is None ) : styleFilter = alwaysTrue

        parent = self.ancestor
        derivedFrom = self.derivedFrom
        while( derivedFrom != '' ) :
            for form in component :
                if( ( derivedFrom == form.label ) and styleFilter( parent[derivedFrom] ) ) : return( form )
            derivedFrom = parent[derivedFrom].derivedFrom
        return( None )
        
    def findDerivedFromStyle( self, cls ) :
        """
        Method to find the derived style that corresponds to the input argument style.
        """

        _style = self
        while( _style is not None ) :
            _style = _style.derivedFromStyle
            if( ( _style is not None ) and ( isinstance( _style, cls ) ) ) : return( _style )
        return( None )

    def hasTemperature( self ) :
        """
        Base method to indicate whether style has a temperature property.
        """

        return( False )

    def sibling( self, label ) :
        """Returns the sibling of the current style which has a label that correspond to the input argument."""

        return( self.ancestor[label] )

    def convertUnits( self, unitMap ) :
        """
        Base method to convert the style's units
        """

        pass

    def XMLCommonAttributes( self ) :
        """
        Return a string with the style label, and if available, the date and derivedFrom properties.
        """

        XMLCommon = 'label="%s"' % self.label
        if( self.derivedFrom != "" ) : XMLCommon += ' derivedFrom="%s"' % self.derivedFrom
        if( self.date is not None ) : XMLCommon += ' date="%s"' % self.date
        return( XMLCommon )

    def toXML_strListCommon(self, indent = '', **kwargs):

        xmlStringList = self.documentation.toXML_strList( indent, **kwargs )

        return xmlStringList

    def parseCommonChildNodes(self, node, xPath, linkData, **kwargs):

        documentation = node.find(documentationModule.Documentation.moniker)
        if documentation is not None:
            self.documentation.parseNode(documentation, xPath, linkData, **kwargs)

    @staticmethod
    def parseNodeBase(node, xPath, linkData, **kwargs):
        """
        Update the XPath with the current node label and return the attributes label, derivedFrom and date.
        """

        label = node.get('label')
        xPath.append('%s[@label="%s"]' % (node.tag, label))

        derivedFrom = node.get('derivedFrom', '')
        date = node.get('date', None)

        return label, derivedFrom, date

class StyleWithTemperature( Style ) :
    """
    Base class for the styles with a temperature attribute.
    """

    def __init__( self, label, derivedFrom, temperature, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( temperature, ( physicalQuantityModule.Temperature, None.__class__ ) ) ) ) : raise TypeError( 'invalid temperature object' )
        self.__temperature = temperature

    @property
    def temperature( self ) :
        """
        Property for the style temperature.
        """

        if( self.__temperature is None ) : return( self.derivedFromStyle.temperature )
        return( self.__temperature )

    def convertUnits(self, unitMap):
        """
        Method to convert the units of the style temperature.
        """

        Style.convertUnits(self, unitMap)
        if self.__temperature is not None: self.__temperature.convertUnits(unitMap)

    def hasTemperature( self ) :
        """
        Method indicating if the style has a temperature attribute.
        """

        return( self.__temperature is not None )

    def toXML_strListCommon(self, indent = '', **kwargs):

        xmlStringList = Style.toXML_strListCommon(self, indent, **kwargs)
        if self.hasTemperature(): xmlStringList += self.temperature.toXML_strList(indent, **kwargs)

        return xmlStringList

    @staticmethod
    def parseNodeBase(node, xPath, linkData, **kwargs):
        """
        This methods reads the label, derivedFrom and date attributes, and the temperature child from *node* and
        return them as a tuple. The temperature is tranformed to an xData.physicalQuantity.PhysicalQantity
        object before it is returned in the tuple.
        """

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)
        temperature = node.find( 'temperature' )
        if temperature is not None: temperature = physicalQuantityModule.Temperature.parseNodeUsingClass(temperature, xPath, linkData, **kwargs)

        return label, derivedFrom, date, temperature

class Evaluated( StyleWithTemperature ) :
    """
    Style for data entered by an evaluator.
    """

    moniker = 'evaluated'
    sortOrderIndex = 0

    def __init__( self, label, derivedFrom, temperature, projectileEnergyDomain, library, version,
                  date = None, targetInfo = None ) :

        StyleWithTemperature.__init__( self, label, derivedFrom, temperature, date = date )

        if( not( isinstance( library, str ) ) ) : raise TypeError( 'library must be a string' )
        self.__library = library

        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = version

        self.projectileEnergyDomain = projectileEnergyDomain

        self.targetInfo = targetInfo

    @property
    def library( self ) :
        """
        Library property, i.e. the name of the library for this evaluation.
        """

        return( self.__library )

    @library.setter
    def library( self, value ) :
        """
        Setter for the library property, i.e. the name of the library for this evaluation.
        """

        self.__library = value

    @property
    def projectileEnergyDomain( self ) :
        """
        Projectile energy domain property.
        """

        if( self.__projectileEnergyDomain is None ) : return( self.derivedFromStyle.projectileEnergyDomain )
        return self.__projectileEnergyDomain

    @projectileEnergyDomain.setter
    def projectileEnergyDomain( self, value ) :
        """
        Setter for the projectile energy domain property.
        """

        if( not( isinstance( value, ( ProjectileEnergyDomain, None.__class__ ) ) ) ) :
            raise TypeError("Expected ProjectileEnergyDomain instance, got '%s' instead" % type(value))
        if( value is not None ) : value.setAncestor( self )
        self.__projectileEnergyDomain = value

    @property
    def targetInfo(self):
        """
        Target info property (stores isotopic abundances for TNSL evaluations).
        """

        if self.__targetInfo is None and self.derivedFromStyle is not None:
            return self.derivedFromStyle.targetInfo
        return self.__targetInfo

    @targetInfo.setter
    def targetInfo(self, value):
        """
        Setter for target info property
        """

        if not isinstance(value, (targetInfoModule.TargetInfo, type(None))):
            raise TypeError("Expected TargetInfo instance, got '%s' instead" % type(value))
        if value is not None: value.setAncestor(self)
        self.__targetInfo = value

    @property
    def version( self ) :
        """
        Version property, i.e. the version of this library's release.
        """

        return( self.__version )

    @version.setter
    def version( self, value ) :
        """
        Setter for the version property, i.e. the version of this library's release.
        """

        self.__version = value

    def convertUnits(self, unitMap):
        """
        Convert units for the style temperature and projectile energy domain.
        """

        StyleWithTemperature.convertUnits(self, unitMap)
        if self.__projectileEnergyDomain is not None: self.__projectileEnergyDomain.convertUnits(unitMap)

    def fixDomains(self, energyMax):
        """
        Calls the *projectileEnergyDomain's **fixDomains** method.
        """

        if self.__projectileEnergyDomain is not None:
            return self.__projectileEnergyDomain.fixDomains(energyMax)

        return 0

    def toXML_strList(self, indent='', **kwargs):
        """Return styles XML node and its child nodes as a list"""

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlStringList = ['%s<%s %s library="%s" version="%s">' %
                         (indent, self.moniker, self.XMLCommonAttributes(), self.library, self.version)]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        if self.__projectileEnergyDomain is not None:
            xmlStringList += self.projectileEnergyDomain.toXML_strList(indent2, **kwargs)
        if self.targetInfo is not None:
            # targetInfo was introduced in GNDS-2.1
            formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
            versionTuple = GNDS_formatVersionModule.versionComponents(formatVersion)
            if versionTuple[:2] >= (2, 1):
                xmlStringList += self.targetInfo.toXML_strList(indent2, **kwargs)

        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date, temperature = StyleWithTemperature.parseNodeBase(node, xPath, linkData, **kwargs)

        library = node.get('library')
        version = node.get('version')
        _projectileEnergyDomain = node.find(ProjectileEnergyDomain.moniker)
        if _projectileEnergyDomain is not None:
            _projectileEnergyDomain = ProjectileEnergyDomain.parseNodeUsingClass(
                _projectileEnergyDomain, xPath, linkData, **kwargs)

        evaluated = cls(label, derivedFrom, temperature, _projectileEnergyDomain, library, version, date=date)
        evaluated.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        _targetInfo = node.find(targetInfoModule.TargetInfo.moniker)
        if _targetInfo is not None:
            _targetInfo = targetInfoModule.TargetInfo.parseNodeUsingClass(
                _targetInfo, xPath, linkData, **kwargs)
            evaluated.targetInfo = _targetInfo

        xPath.pop()

        return evaluated


class ProjectileEnergyDomain( ancestryModule.AncestryIO ) :
    """
    Style for the storage of the energy domain and the unit of the projectile in the laboratory frame.
    """

    moniker = 'projectileEnergyDomain'

    def __init__(self, min, max, unit):

        ancestryModule.AncestryIO.__init__(self)
        self.min = min
        self.max = max
        self.unit = unit

    def convertUnits( self, unitMap ) :
        """
        Method for the conversion of the energy units.
        """

        if( self.unit in unitMap ) :
            factor = PQUModule.PQU( 1, self.unit ).getValueAs( unitMap[self.unit] )
            self.min *= factor
            self.max *= factor
            self.unit = unitMap[self.unit]

    def fixDomains(self, energyMax):
        """
        Sets *max* to the lesser of its current value and *energyMax*.
        """

        if self.max > energyMax:
            self.max = energyMax
            return 1

        return 0

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        return ['%s<%s min="%s" max="%s" unit="%s"/>'
                % (indent, self.moniker, PQUModule.floatToShortestString( self.min ), PQUModule.floatToShortestString( self.max ), self.unit)]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        xPath.append(node.tag)

        PED = cls( float(node.get('min')), float(node.get('max')), node.get('unit') )

        xPath.pop()

        return PED


class CrossSectionReconstructed(StyleWithTemperature):
    """
    Style for cross section data reconstructed from resonance parameters.
    """

    moniker = 'crossSectionReconstructed'
    sortOrderIndex = 2

    def __init__(self, label, derivedFrom, temperature=None, date=None):

        StyleWithTemperature.__init__(self, label, derivedFrom, temperature, date=date)

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date, temperature = StyleWithTemperature.parseNodeBase(node, xPath, linkData, **kwargs)

        crossSectionReconstructed = cls(label, derivedFrom, temperature=temperature, date=date)
        crossSectionReconstructed.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return crossSectionReconstructed

class AngularDistributionReconstructed(StyleWithTemperature):
    """
    Style for angular distribution data reconstructed from resonance parameters. 
    """

    moniker = 'angularDistributionReconstructed'
    sortOrderIndex = 3

    def __init__(self, label, derivedFrom, temperature=None, date=None ) :

        StyleWithTemperature.__init__(self, label, derivedFrom, temperature, date=date)

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date, temperature = StyleWithTemperature.parseNodeBase(node, xPath, linkData, **kwargs)

        angularDistributionReconstructed = cls(label, derivedFrom, temperature=temperature, date=date)
        angularDistributionReconstructed.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return angularDistributionReconstructed

class CoulombPlusNuclearElasticMuCutoff( Style ) :
    """
    This is a derived data style for Coulomb interactions which have a singularity at 
    :math:`\\mu = 1.0`. The cross sections may consequently represented by the integral of 
    :math:`\\frac{d\\sigma \\left(E, \\mu\\right)}{d\\mu}` from :math:`\\mu = -1.0` to
    :math:`\\mu_{\\textrm{cutoff}}`, i.e. the cutoff scattering cosine which satisfies the condition 
    :math:`\\mu_{\\textrm{cutoff}} < 1.0`. 
    """

    moniker = 'CoulombPlusNuclearElasticMuCutoff'
    sortOrderIndex = 10

    def __init__( self, label, derivedFrom, muCutoff, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( muCutoff, float ) ) ) : raise TypeError( 'invalid muCutoff object' )
        self.__muCutoff = muCutoff

    @property
    def muCutoff( self ) :
        """Return the value for :math:`\\mu_{\\textrm{cutoff}}`"""

        return( self.__muCutoff )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s muCutoff="%s" %s>' % ( indent, self.moniker, self.__muCutoff, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        muCutoff = float( node.get( 'muCutoff' ) )
        coulombPlusNuclearElasticMuCutoff = cls( label, derivedFrom, muCutoff, date = date )
        coulombPlusNuclearElasticMuCutoff.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return coulombPlusNuclearElasticMuCutoff

class Heated( StyleWithTemperature ) :
    """
    Style for data that has been heated to a specific temperature, e.g. Doppler broaded cross sections.
    """

    moniker = 'heated'
    sortOrderIndex = 100

    def __init__( self, label, derivedFrom, temperature, date = None ) :

        StyleWithTemperature.__init__( self, label, derivedFrom, temperature, date = date )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date, temperature = StyleWithTemperature.parseNodeBase(node, xPath, linkData, **kwargs)

        heated = cls( label, derivedFrom, temperature=temperature, date = date )
        heated.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return heated

class AverageProductData(StyleWithTemperature):
    """
    Style for data calculated for the averageProductEnergy and averageMomentum components.
    """

    moniker = 'averageProductData'
    sortOrderIndex = 11

    def __init__(self, label, derivedFrom, temperature=None, date=None):

        StyleWithTemperature.__init__(self, label, derivedFrom, temperature, date=date)

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date, temperature = StyleWithTemperature.parseNodeBase(node, xPath, linkData, **kwargs)

        averageProductData = cls(label, derivedFrom, temperature=temperature, date=date)
        averageProductData.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return averageProductData

class MultiGroup( Style ) :
    """
    Style for the multigroup data. This style is deprecated (i.e., it does not exists it GNDS 2.0 or higher).
    """

    moniker = 'multiGroup'
    sortOrderIndex = -1

    def __init__( self, label, lMax, date = None ) :

        Style.__init__( self, label, "", date = date )
        self.lMax = int( lMax )

        self.transportables = transportablesModule.Transportables( )
        self.transportables.setAncestor( self )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s lMax="%s">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.lMax ) ]
        xmlStringList += self.transportables.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        lMax = node.get( 'lMax', 3 )

        _multiGroup = cls( label, lMax, date = date)
        for child in node:
            if child.tag == transportablesModule.Transportables.moniker:
                _multiGroup.transportables.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Encountered unexpected node '%s'" % child.tag)

        xPath.pop( )
        return( _multiGroup )

class HeatedMultiGroup(Style):
    """
    Style for data that has been heated and multi-grouped.
    """

    moniker = 'heatedMultiGroup'
    sortOrderIndex = 120

    def __init__(self, label, derivedFrom, flux, date=None, parameters=''):

        Style.__init__(self, label, derivedFrom, date=date)

        self.__transportables = transportablesModule.Transportables( )
        self.__transportables.setAncestor( self )

        if( not( isinstance( flux, fluxModule.Flux ) ) ) : raise TypeError( 'Invalid flux instance' )
        self.__flux = flux
        self.__flux.setAncestor( self )

        self.__inverseSpeed = None

        if( not( isinstance( parameters, str ) ) ) : raise TypeError( 'parameters must be a str instance' )
        self.parameters = parameters              # This is only used if GNDS 1.10.

        self.__multiGroupFlux = None                # Only for internal use, it is not written to GNDS file.

    @property
    def transportables( self ) :
        """
        Property for the transportables, i.e. the list of particles for which the multi-group processing was performed.
        """

        return( self.__transportables )

    @property
    def flux( self ) :
        """
        Property for the flux (used to collapse to the multi-group energy structure) as a fudge.processing.flux.flux object.
        """

        return( self.__flux )

    @property
    def inverseSpeed( self ) :
        """
        Property for the inverse speed.

        The group multi-group inverse speeds (:math:`\\left(\\upsilon_g\\right)`) are calculated from 
        :math:`\\upsilon_g = \\int_g dE \\frac{\\phi\\left(E\\right)}{\\upsilon\\left(E\\right)}`,
        where :math:`E` is the projectile energy, :math:`\\phi\\left(E\\right)` is the scalar flux,
        :math:`\\upsilon\\left(E\\right)` is the projectiles velocity, and the integral is over the
        energy range in group :math:`g`.
        """

        return( self.__inverseSpeed )

    @inverseSpeed.setter
    def inverseSpeed( self, value ) :
        """
        Setter for the inverseSpeed property.
        """

        self.__inverseSpeed = value

    @property
    def multiGroupFlux( self ) :
        """
        Property for the flux (used to collapse to the multi-group energy structure) as a fudge.processing.flux.Gridded2d object.
        """

        return( self.__multiGroupFlux )

    def convertUnits(self, unitMap):
        '''
        Converts unit per *unitMap*.
        '''

        Style.convertUnits(self, unitMap)
        self.__transportables.convertUnits(unitMap)
        self.__flux.convertUnits(unitMap)
        self.__inverseSpeed.convertUnits(unitMap)

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Method to set the multiGroupFlux and inverseSpeed properties.
        """

        self.__multiGroupFlux = self.flux.processMultiGroup( self, tempInfo, indent )
        self.inverseSpeed = inverseSpeedModule.InverseSpeed( inverseSpeedModule.multiGroupInverseSpeed( self, tempInfo ) )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        formatVersion = kwargs['formatVersion']
        parameters = kwargs.get( 'parameters', '' )

        attributes = self.XMLCommonAttributes( )
        if( parameters != '' ) : attributes += ' parameters="%s"' % parameters

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, attributes ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        if( formatVersion != GNDS_formatVersionModule.version_1_10 ) : xmlStringList += self.transportables.toXML_strList( indent2, **kwargs )
        xmlStringList += self.flux.toXML_strList( indent2, **kwargs )
        if( self.inverseSpeed is not None ) : xmlStringList += self.inverseSpeed.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)
        parameters = node.get( 'parameters', '' )

        transportables = None
        flux = None
        inverseSpeed = None
        for child in node:
            if( child.tag == transportablesModule.Transportables.moniker ) :
                transportables = child
            elif( child.tag == fluxModule.Flux.moniker ) :
                flux = fluxModule.Flux.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif( child.tag == inverseSpeedModule.InverseSpeed.moniker ) :
                inverseSpeed = inverseSpeedModule.InverseSpeed.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        heatedMultiGroup = cls(label, derivedFrom, flux, date = date, parameters = parameters)
        heatedMultiGroup.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        if transportables is not None: heatedMultiGroup.transportables.parseNode(transportables, xPath, linkData, **kwargs)
        heatedMultiGroup.inverseSpeed = inverseSpeed

        xPath.pop()

        return heatedMultiGroup

class MonteCarlo_cdf( Style ) :
    """
    This is a processed data style with distributions processed for use in Monte Carlo transport codes which
    requires:

        - Conversion of :math:`P\\left(\\mu, E ^{\\prime} \\mid E \\right)` to
          :math:`P\\left(\\mu \\mid E\\right) \\times P\\left(E ^{\\prime} \\mid E, \\mu\\right)`; and
        - Conversion of :math:`P\\left(E ^{\\prime}, \\mu \\mid E \\right)` to
          :math:`P\\left(E ^{\\prime} \\mid E \\right) \\times P\\left(\\mu \\mid E, E ^{\\prime} \\right)`

    """

    moniker = 'MonteCarlo_cdf'
    sortOrderIndex = 12

    def __init__( self, label, derivedFrom, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        monteCarlo_cdf = cls(label, derivedFrom, date=date )
        monteCarlo_cdf.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return monteCarlo_cdf

class SnElasticUpScatter( Style ) :
    """
    This is a thermal elastic upscatter style for data that has been heated and multi-grouped.
    This is similar to the heatedMultigroup style but with an upscatter correction included for
    elastic scattering.
    """

    moniker = 'SnElasticUpScatter'
    sortOrderIndex = 121

    def __init__( self, label, derivedFrom, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )
        
        self.__upperCalculatedGroup = 0 
        

    @property
    def upperCalculatedGroup( self ) :
        """Return or set  the index of the highest-energy group (counting from 0) impacted by the upscatter correction."""

        return( self.__upperCalculatedGroup )
    
    @upperCalculatedGroup.setter
    def upperCalculatedGroup( self, value ) :
        """Return or set  the index of the highest-energy group (counting from 0) impacted by the upscatter correction."""

        #if( not( isinstance( upperCalculatedGroup, int ) ) ) : raise TypeError( 'upperCalculatedGroup must be an integer' )
        self.__upperCalculatedGroup = value

    @property
    def flux( self ) :
        """Property for the flux used to collapse to group structure."""

        return( self.derivedFromStyle.flux )
    
    @property
    def transportables( self ) :
        """Property for the transportables, i.e. the list of particles for which the multi-group processing was performed."""

        return( self.derivedFromStyle.transportables )

    @property
    def inverseSpeed( self ) :
        """
        Property for the inverse speed.

        The group multi-group inverse speeds (:math:`\\left(\\upsilon_g\\right)`) are calculated from 
        :math:`\\upsilon_g = \\int_g dE \\frac{\\phi\\left(E\\right)}{\\upsilon\\left(E\\right)}`,
        where :math:`E` is the projectile energy, :math:`\\phi\\left(E\\right)` is the scalar flux,
        :math:`\\upsilon\\left(E\\right)` is the projectiles velocity, and the integral is over the
        energy range in group :math:`g`.
        """

        return( self.derivedFromStyle.inverseSpeed )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s upperCalculatedGroup="%d">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.upperCalculatedGroup,  ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        snElasticUpScatter = cls( label, derivedFrom, date = date )
        snElasticUpScatter.parseCommonChildNodes(node, xPath, linkData, **kwargs)
        snElasticUpScatter.upperCalculatedGroup = int(node.get('upperCalculatedGroup'))
        
        xPath.pop()

        return snElasticUpScatter

class GriddedCrossSection( Style ) :
    """
    Style for cross section data that has been converted to a common energy grid as used in Monte Carlo Transport.
    """

    moniker = 'griddedCrossSection'
    sortOrderIndex = 110

    def __init__( self, label, derivedFrom, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

        self.__grid = None

    @property
    def grid( self ) :
        """
        Common energy grid (as used in Monte Carlo transport) property.
        """

        return( self.__grid )

    @grid.setter
    def grid( self, grid ) :
        """
        Setter for the common energy grid (as used in Monte Carlo transport) property.
        """

        if( not( isinstance( grid, axesModule.Grid ) ) ) : raise TypeError( 'grid not instance of axesModule.Grid' )
        self.__grid = grid
        self.__grid.setAncestor( self )

    def convertUnits(self, unitMap):
        """
        Convert unit per *unitMap*.
        """

        Style.convertUnits(self, unitMap)
        self.__grid.convertUnits(unitMap)

    def findEntity(self, entityName, attribute = None, value = None):

        if entityName == axesModule.Grid.moniker:
            return self.grid

        raise ValueError('Could not find grid with label "%s".' % (entityName, value))

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList += self.__grid.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        grid = None
        for child in node:
            if( child.tag == axesModule.Grid.moniker ) :
                grid = axesModule.Grid.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        griddedCrossSection = cls( label, derivedFrom, date = date )
        griddedCrossSection.grid = grid
        griddedCrossSection.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return griddedCrossSection

class URR_probabilityTables( Style ) :
    """
    This is a precessed style for cross section probability tables :math:`P\\left(\\sigma \\mid E\\right)` 
    to better capture the possible rapid cross section fluctuations in the unresolved resonance region.
    """

    moniker = 'URR_probabilityTables'
    sortOrderIndex = 111

    def __init__( self, label, derivedFrom, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = ['%s<%s %s>' % (indent, self.moniker, self.XMLCommonAttributes())]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        URR_probability_tables = cls( label, derivedFrom, date = date )
        URR_probability_tables.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return URR_probability_tables

class Realization( Style ) :
    """
    This style represents a change to the mean value of some of the data.
    """

    moniker = 'realization'
    sortOrderIndex = 1

    def __init__( self, label, derivedFrom, date = None ) :

        Style.__init__( self, label, derivedFrom, date = date )

    def toXML_strList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = ['%s<%s %s>' % (indent, self.moniker, self.XMLCommonAttributes())]
        xmlStringList += self.toXML_strListCommon(indent2, **kwargs)
        xmlStringList[-1] += '</%s>' % self.moniker

        return xmlStringList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Generate an instance of this class from *node*."""

        label, derivedFrom, date = Style.parseNodeBase(node, xPath, linkData, **kwargs)

        realization = cls( label, derivedFrom, date = date )
        realization.parseCommonChildNodes(node, xPath, linkData, **kwargs)

        xPath.pop()

        return realization

def findEvaluated( forms ) :
    """
    Find the "evaluated" style class from a list of style classes. It calls a method that ensures that 
    the list of style classes contains exactly one occurence of the "evaluated" style.
    """

    return( findStyle( Evaluated, forms ) )

def findAllOfStyle( cls, forms ) :
    """
    Find all the occurences of a guven style class from a list of style classes.
    """

    formsFound = []
    for form in forms :
        if( isinstance( form, cls ) ) : formsFound.append( form )
    return( formsFound )

def findStyle( cls, forms ) :
    """
    Find all occurences of a given style class and ensure that there is only a single occurence.
    """

    formsFound = findAllOfStyle( cls, forms )
    if( len( formsFound )  > 1 ) : raise Exception( 'multiple styles found' )
    if( len( formsFound ) == 0 ) : raise Exception( 'no style found' )
    return( formsFound[0] )
