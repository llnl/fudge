# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation class.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Documentation             | This is the class for the GNDS documentation node.                                |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from .. import date as dateModule
from . import author as authorModule
from . import contributor as contributorModule
from . import collaboration as collaborationsModule
from . import dates as datesModule
from . import copyright as copyrightModule
from . import acknowledgement as acknowledgementModule
from . import computerCode as computerCodeModule
from . import keyword as keywordModule
from . import relatedItem as relatedItemModule
from . import bibitem as bibitemModule
from . import exforDataSet as exforDataSetModule
from . import texts as textsModule
from .. import text as textParentModule


class Documentation(ancestryModule.AncestryIO):
    """
    The class for the GNDS documentation node.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | doi                   | This is the doi if the data/publication/code were published.              |
    +-----------------------+---------------------------------------------------------------------------+
    | version               | This is the version of the data/publication/code for this documentation.  |
    +-----------------------+---------------------------------------------------------------------------+
    | publicationDate       | This is the date the data/publication/code was published.                 |
    +-----------------------+---------------------------------------------------------------------------+
    | authors               | This is the list of authors for the data/publication/code.                |
    +-----------------------+---------------------------------------------------------------------------+
    | contributors          | This is the list of contributors for the data/publication/code.           |
    +-----------------------+---------------------------------------------------------------------------+
    | collaborations        | This is the list of collaborations for the data/publication/code.         |
    +-----------------------+---------------------------------------------------------------------------+
    | dates                 | This is the list of dates for the data/publication/code.                  |
    +-----------------------+---------------------------------------------------------------------------+
    | copyright             | This is the copy right for the data/publication/code.                     |
    +-----------------------+---------------------------------------------------------------------------+
    | acknowledgements      | This is the list of acknowledgements for the data/publication/code.       |
    +-----------------------+---------------------------------------------------------------------------+
    | keywords              | This is the list of keywords for the data/publication/code.               |
    +-----------------------+---------------------------------------------------------------------------+
    | relatedItems          | This is the list of relatedItems for the data/publication/code.           |
    +-----------------------+---------------------------------------------------------------------------+
    | title                 | This is the title for the data/publication/code.                          |
    +-----------------------+---------------------------------------------------------------------------+
    | abstract              | This is the abstract for the data/publication/code.                       |
    +-----------------------+---------------------------------------------------------------------------+
    | body                  | This is the body for the data/publication/code.                           |
    +-----------------------+---------------------------------------------------------------------------+
    | computerCodes         | This is the list of computer codes used for the data/publication.         |
    +-----------------------+---------------------------------------------------------------------------+
    | experimentalDataSets  | This is the list of experimental datasets used for the data/publication.  |
    +-----------------------+---------------------------------------------------------------------------+
    | bibliography          | This is the list of bibliographies for the data/publication/code.         |
    +-----------------------+---------------------------------------------------------------------------+
    | endfCompatible        | This is text that is limited to 66 characters per line that can be        |
    |                       | written to an ENDF-6 formatted file.                                      |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = 'documentation'
    ancestryMembers = (
        'authors', 'contributors', 'collaborations', 'dates', 'copyright', 'acknowledgements', 'keywords',
        'relatedItems', 'title', 'abstract', 'body', 'computerCodes', 'experimentalDataSets', 'bibliography',
        'endfCompatible')

    def __init__(self, doi='', version='', publicationDate=None):
        """
        :param doi:
        :param version:
        :param publicationDate:
        """

        ancestryModule.AncestryIO.__init__(self)

        self.doi = doi
        self.version = version
        self.publicationDate = publicationDate

        self.__authors = authorModule.Authors()
        self.__authors.setAncestor(self)

        self.__contributors = contributorModule.Contributors()
        self.__contributors.setAncestor(self)

        self.__collaborations = collaborationsModule.Collaborations()
        self.__collaborations.setAncestor(self)

        self.__dates = datesModule.Dates()
        self.__dates.setAncestor(self)

        self.__copyright = copyrightModule.Copyright()
        self.__copyright.setAncestor(self)

        self.__acknowledgements = acknowledgementModule.Acknowledgements()
        self.__acknowledgements.setAncestor(self)

        self.__keywords = keywordModule.Keywords()
        self.__keywords.setAncestor(self)

        self.__relatedItems = relatedItemModule.RelatedItems()
        self.__relatedItems.setAncestor(self)

        self.__title = textsModule.Title()
        self.__title.setAncestor(self)

        self.__abstract = textsModule.Abstract()
        self.__abstract.setAncestor(self)

        self.__body = textsModule.Body()
        self.__body.setAncestor(self)

        self.__computerCodes = computerCodeModule.ComputerCodes()
        self.__computerCodes.setAncestor(self)

        self.__experimentalDataSets = exforDataSetModule.ExperimentalDataSets()
        self.__experimentalDataSets.setAncestor(self)

        self.__bibliography = bibitemModule.Bibliography()
        self.__bibliography.setAncestor(self)

        self.__endfCompatible = textsModule.EndfCompatible()
        self.__endfCompatible.setAncestor(self)

    @property
    def doi(self):
        """
        This method returns the *doi* member of *self*.

        :returns:       An instance of :py:class:`textParentModule.Text`.
        """

        return self.__doi

    @doi.setter
    def doi(self, doi):
        """
        This method sets the *doi* member of *self* to *doi*.

        :param doi:     An instance of :py:class:`textParentModule.Text`.
        """

        self.__doi = textParentModule.raiseIfNotString(doi, 'doi')

    @property
    def version(self):
        """
        This method returns the *version* member of *self*.

        :returns:           An instance of :py:class:`textParentModule.Text`.
        """

        return self.__version

    @version.setter
    def version(self, version):
        """
        This method sets the *version* member of *self* to *version*.

        :param version:     An instance of :py:class:`textParentModule.Text`.
        """

        self.__version = textParentModule.raiseIfNotString(version, 'version')

    @property
    def publicationDate(self):
        """
        This method returns the *publicationDate* member of *self*.

        :returns:           An instance of :py:class:`dateModule.Date`.
        """

        return self.__publicationDate

    @publicationDate.setter
    def publicationDate(self, date):
        """
        This method sets the *date* member of *self* to *date*.

        :param date:        An instance of :py:class:`dateModule.Date`.
        """

        if date is None:
            self.__publicationDate = dateModule.Date(None, dateModule.Resolution.undefined)
        else:
            if not isinstance(date, dateModule.Date): raise TypeError('Invalid date.')
            self.__publicationDate = date

    @property
    def authors(self):
        """
        This method returns the *authors* member of *self*.

        :returns:           An instance of :py:class:`authorModule.Authors`.
        """

        return self.__authors

    @property
    def contributors(self):
        """
        This method returns the *contributors* member of *self*.

        :returns:           An instance of :py:class:`contributorModule.Contributor`.
        """

        return self.__contributors

    @property
    def collaborations(self):
        """
        This method returns the *collaborations* member of *self*.

        :returns:           An instance of :py:class:`collaborationsModule.Collaborations`.
        """

        return self.__collaborations

    @property
    def dates(self):
        """
        This method returns the *dates* member of *self*.

        :returns:           An instance of :py:class:`datesModule.Dates`.
        """

        return self.__dates

    @property
    def copyright(self):
        """
        This method returns the *copyright* member of *self*.

        :returns:           An instance of :py:class:`copyrightModule.Copyright`.
        """

        return self.__copyright

    @property
    def acknowledgements(self):
        """
        This method returns the *acknowledgements* member of *self*.

        :returns:           An instance of :py:class:`acknowledgementModule.Acknowledgements`.
        """

        return self.__acknowledgements

    @property
    def keywords(self):
        """
        This method returns the *keywords* member of *self*.

        :returns:           An instance of :py:class:`keywordModule.Keywords`.
        """

        return self.__keywords

    @property
    def relatedItems(self):
        """
        This method returns the *relatedItems* member of *self*.

        :returns:           An instance of :py:class:`relatedItemModule.RelatedItems`.
        """

        return self.__relatedItems

    @property
    def title(self):
        """
        This method returns the *title* member of *self*.

        :returns:           An instance of :py:class:`textsModule.Title`.
        """

        return self.__title

    @property
    def abstract(self):
        """
        This method returns the *abstract* member of *self*.

        :returns:           An instance of :py:class:`textsModule.Abstract`.
        """

        return self.__abstract

    @property
    def body(self):
        """
        This method returns the *body* member of *self*.

        :returns:           An instance of :py:class:`textsModule.Body`.
        """

        return self.__body

    @property
    def computerCodes(self):
        """
        This method returns the *computerCodes* member of *self*.

        :returns:           An instance of :py:class:`computerCodeModule.ComputerCodes`.
        """

        return self.__computerCodes

    @property
    def experimentalDataSets(self):
        """
        This method returns the *experimentalDataSets* member of *self*.

        :returns:           An instance of :py:class:`exforDataSetModule.ExperimentalDataSets`.
        """

        return self.__experimentalDataSets

    @property
    def bibliography(self):
        """
        This method returns the *bibliography* member of *self*.

        :returns:           An instance of :py:class:`bibitemModule.Bibliography`.
        """

        return self.__bibliography

    @property
    def endfCompatible(self):
        """
        This method returns the *endfCompatible* member of *self*.

        :returns:           An instance of :py:class:`textsModule.EndfCompatible`.
        """

        return self.__endfCompatible

    def convertUnits(self, unitMap):
        """This method does nothing as there are no units to convert in documentation."""

        return

    def findEntriesWithKey(self, keyValue):
        """
        This method returns the list of each entry in *self*'s authors, contributors, collaborations, acknowledgements,
        relatedItems and computerCodes which have the key *keyValue*.

        :param keyValue:        The value of the key.
        """

        entries = []
        if keyValue in self.__authors:
            entries.append(self.__authors[keyValue])
        if keyValue in self.__contributors:
            entries.append(self.__contributors[keyValue])
        if keyValue in self.__collaborations:
            entries.append(self.__collaborations[keyValue])
        if keyValue in self.__acknowledgements:
            entries.append(self.__acknowledgements[keyValue])
        if keyValue in self.__relatedItems:
            entries.append(self.__relatedItems[keyValue])
        if keyValue in self.__computerCodes:
            entries.append(self.__computerCodes[keyValue])

        return entries

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:       The minimum amount of indentation.
        :param kwargs:       dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:             List of str instances representing the XML lines of self.
        """

        if kwargs.get('skipDocumentation', False):
            return []

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XMLList = self.__authors.toXML_strList(indent2, **kwargs)
        XMLList += self.__contributors.toXML_strList(indent2, **kwargs)
        XMLList += self.__collaborations.toXML_strList(indent2, **kwargs)
        XMLList += self.__dates.toXML_strList(indent2, **kwargs)
        XMLList += self.__copyright.toXML_strList(indent2, **kwargs)
        XMLList += self.__acknowledgements.toXML_strList(indent2, **kwargs)
        XMLList += self.__keywords.toXML_strList(indent2, **kwargs)
        XMLList += self.__relatedItems.toXML_strList(indent2, **kwargs)
        XMLList += self.__title.toXML_strList(indent2, **kwargs)
        XMLList += self.__abstract.toXML_strList(indent2, **kwargs)
        XMLList += self.__body.toXML_strList(indent2, **kwargs)
        XMLList += self.__computerCodes.toXML_strList(indent2, **kwargs)
        XMLList += self.__experimentalDataSets.toXML_strList(indent2, **kwargs)
        XMLList += self.__bibliography.toXML_strList(indent2, **kwargs)
        XMLList += self.__endfCompatible.toXML_strList(indent2, **kwargs)

        attributes = '' if self.doi == '' else ' doi="%s"' % self.doi
        attributes += '' if self.version == '' else ' version="%s"' % self.version
        attributes += self.__publicationDate.asXML_attribute(name='publicationDate')

        if len(XMLList) == 0 and len(attributes) == 0:
            return []

        XMLList.insert(0, '%s<%s%s>' % (indent, self.moniker, attributes))
        XMLList[-1] += '</%s>' % self.moniker

        return XMLList

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method fills *self* by parsing the data in *node*.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      dictionary of extra arguments controlling how *self* is converted to a list of XML strings.
        """

        xPath.append(node.tag)

        self.parseAncestryMembers(node, xPath, linkData, **kwargs)

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        doi = node.get('doi', '')
        version = node.get('version', '')
        publicationDate = dateModule.Date.parse(node.get('publicationDate', ''))

        return cls(doi, version, publicationDate)
