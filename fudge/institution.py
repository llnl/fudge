# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
This module defines the 'institution' class, used inside GNDS to store institution-specific data
that does not necessarily conform to GNDS specifications.

Within a reactionSuite (or possibly covarianceSuite), institution data is stored inside an <applicationData> section.
Each institution has a unique label and a list of data types. The only restriction on those data types is that
they be well-formed xml. Codes reading in GNDS files are free to ignore any unrecognized institution or data type.
"""

from xml.etree import ElementTree
from xml.sax import saxutils

from LUPY import ancestry as ancestryModule

from fudge.reactions import reaction as reactionModule
from fudge import outputChannel as outputChannelModule
from fudge import product as productModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule
from fudge.processing import nuclearPlusCoulombInterference as nuclearPlusCoulombInterferenceModule
from fudge.processing.deterministic import tokens as deterministicTokensModule
from fudge.resonances import probabilityTables as probabilityTablesModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherentDoppler as photoAtomicIncoherentDopplerModule

photoAtomicIncoherentDoppler = 'LLNL::photoAtomicIncoherentDoppler'

class Institution(ancestryModule.AncestryIO):

    moniker = "institution"
    keyName = 'label'

    def __init__(self, label):

        ancestryModule.AncestryIO.__init__(self)
        self.__label = label
        self.__data = []

    def __getitem__(self, item):

        return self.__data[item]

    @property
    def label(self): return self.__label

    @property
    def data(self): return self.__data

    def append( self, item ) :

        self.__data.append( item )
        item.setAncestor( self )
        setattr(self, item.moniker, item)

    def convertUnits(self, unitMap):
        '''
        This method calls convertUnits on the item of *self* if it has a *convertUnits* method.
        '''

        for item in self.__data:
            if hasattr(item, 'convertUnits'):
                item.convertUnits(unitMap)

    def findLinks( self, links ) :

        for item in self :
            if( hasattr( item, 'findLinks' ) ) : getattr( item, 'findLinks' )( links )

    def parseAndReplace(self, index, cls):
        '''
        Replaces the unknown item at *index* with a parsed version of its self using class *cls*.
        The parsed instance is also returned.  This currently assumes that the data are from XML.

        :param index:       The index of *self.__data* of the data to parse and replace.
        :param cls:         The class used to parse the data at *index*.

        :return:            The parsed instance.
        '''

        self.__data[index] = cls.parseXMLString('\n'.join(self.__data[index].node.stringList))

        return self.__data[index]

    def toXML_strList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s label="%s">' % (indent,self.moniker,self.label)]
        for data in self:
            xml += data.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        label = node.get('label')
        xPath.append('%s[@label="%s"]' % (node.tag, label))

        if label == 'LLNL':
            institution = cls(label)
            for child in node:
                if child.tag == nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference.moniker:
                    institution.append(nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference.parseNodeUsingClass(child, xPath, linkData))
                else:
                    institution.append(UnknownLLNL_child(child.tag, UnknownLLNL_XML_child(child)))
        elif label == deterministicTokensModule.multiGroupReactions:
            institution = cls(label)
            for child in node:
                if child.tag == reactionModule.Reaction.moniker:
                    institution.append(reactionModule.Reaction.parseNodeUsingClass(child, xPath, linkData))
                else:                       # This should not happen.
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        elif label == deterministicTokensModule.multiGroupIncompleteProducts:
            institution = cls(label)
            for child in node:
                if child.tag == productModule.Products.moniker:
                    products = productModule.Products()
                    products.parseNode(child, xPath, linkData)
                    institution.append(products)
                else:                       # This should not happen.
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        elif label == deterministicTokensModule.multiGroupDelayedNeutrons:
            institution = cls(label)
            for child in node:
                if child.tag == outputChannelModule.OutputChannel.moniker:
                    institution.append(outputChannelModule.OutputChannel.parseNodeUsingClass(child, xPath, linkData))
                else:                       # This should not happen.
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        elif label == deterministicTokensModule.pointwiseAverageProductEnergies:
            institution = cls(label)
            for child in node:
                if child.tag == averageProductEnergyModule.Component.moniker:
                    institution.append(averageProductEnergyModule.Component.parseNodeUsingClass(child, xPath, linkData))
                else:                       # This should not happen.
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        elif label == photoAtomicIncoherentDoppler:
            institution = cls(label)
            for child in node:
                if child.tag == photoAtomicIncoherentDopplerModule.Reactions.moniker:
                    institution.append(photoAtomicIncoherentDopplerModule.Reactions.parseNodeUsingClass(child, xPath, linkData))
                else:                       # This should not happen.
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        elif label == probabilityTablesModule.LLNLProbabilityTablesToken:
            institution = cls(label)
            for child in node:
                if child.tag == probabilityTablesModule.ProbabilityTables.moniker:
                    institution.append(probabilityTablesModule.ProbabilityTables.parseNodeUsingClass(child, xPath, linkData))
                else:
                    raise Exception('Unsupported child "%s" for institution "%s".' % (child.tag, label))
        else:
            institution = UnknownInstitution(UnknownInstitutionXML_node(node))

        xPath.pop()

        return institution

class UnknownInstitution(ancestryModule.AncestryIO):

    moniker = 'institution2'
    keyName = 'label'

    def __init__(self, node):

        ancestryModule.AncestryIO.__init__(self)

        if not (hasattr(node, 'label') and hasattr(node, 'toXML_strList')): raise TypeError('Node is missing a label and/or toXML_strList method.')
        self.__node = node

    @property
    def label(self):
        """This method returns the label of self."""

        return self.__node.label

    @property
    def node(self):
        """This method returns the node instance of self."""

        return self.__node

    def convertUnits(self, unitMap):
        """ Prints a warning that convertUnits is not implemented for this institution. """

        print('WARNING: convertUnits not implemented for institution "%s".' % self.label)

    def toXML_strList(self, indent = '', **kwargs):
        """This method returng the XML text that was not parsed."""

        return self.__node.toXML_strList(indent = indent, **kwargs)

    def parseNodeWithInstitutionClass(self, cls, **kwargs):

        return self.__node.parseNodeWithInstitutionClass(cls, **kwargs)

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        return cls(UnknownInstitutionXML_node(node))

class UnknownInstitutionXML_node:

    def __init__(self, XML_node):

        self.__label = XML_node.get( 'label' )
        self.__stringList = ElementTree.tostring(XML_node.data, encoding='unicode').rstrip().split('\n')

    @property
    def label( self ) :
        """This method returng the XML text that was not parsed."""

        return( self.__label )

    @property
    def stringList( self ) :
        """This method returng the XML text that was not parsed."""

        return( self.__stringList )

    def toXML_strList( self, indent = '', **kwargs ) :
        """This method returng the XML text that was not parsed."""

        strings = self.stringList
        if( len( strings ) > 0 ) : strings[0] = indent + strings[0]

        return( self.stringList )

    def parseNodeWithInstitutionClass(self, cls, **kwargs):

        string = '\n'.join(self.__stringList)
        string = string[string.find('>')+1:]
        string = string[:string.find('</institution>')]

        return cls.parseXMLString(string, **kwargs)

class UnknownLLNL_child(ancestryModule.Ancestry):
    """
    This class stores a child node of 'applicationData/institution['[@label="LLNL"]' which is not supported by a class under
    the fudge module (i.e., directory). Currently, only the ENDFconversionFlags class is not supported by the fudge module.
    """

    def __init__(self, moniker, node):

        self.__moniker = moniker
        ancestryModule.Ancestry.__init__(self)

        if not isinstance(node, UnknownLLNL_XML_child):
            raise Exception('Invalid UnknownLLNL_child node "%s".' % type(self))

        self.__node = node

    @property
    def moniker(self):
        """This method returns the moniker of self."""

        return self.__moniker

    @property
    def node(self):
        """This method returns the node instance of self."""

        return self.__node

    def toXML_strList(self, indent = '', **kwargs):
        """This method returng the XML text that was not parsed."""

        return self.__node.toXML_strList(indent=indent, **kwargs)

class UnknownLLNL_XML_child:
    """
    An XML sub-node for use with an UnknownLLNL_child instance.
    """

    def __init__(self, XML_node):

        self.__stringList = saxutils.unescape(ElementTree.tostring(XML_node.data, encoding='unicode')).rstrip().split('\n')

    @property
    def stringList(self):
        """This method returng the XML text that was not parsed."""

        return self.__stringList

    def toXML_strList(self, indent='', **kwargs):
        """This method returng the XML text that was not parsed."""

        stringList = self.__stringList
        if len(stringList) > 0:
            if stringList[0][0] != ' ': stringList[0] = indent + stringList[0]

        return stringList
