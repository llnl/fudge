# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.

This module contains the following classes:

    +---------------------------+-------------------------------------------------------------------------------+
    | Class                     | Description                                                                   |
    +===========================+===============================================================================+
    | CovarianceScript          | This class represents a GNDS documentation/covarianceScript node.             |
    +---------------------------+-------------------------------------------------------------------------------+
    | CorrectionScript          | This class represents a GNDS documentation/correctionScript node.             |
    +---------------------------+-------------------------------------------------------------------------------+
    | Note                      | This is the suite class for the GNDS documentation/note node.                 |
    +---------------------------+-------------------------------------------------------------------------------+
    | ExforDataSet              | This class represents a GNDS                                                  |
    |                           |     documentation/experimentalDataSets/exforDataSets/exforDataSet node.       |
    +---------------------------+-------------------------------------------------------------------------------+
    | ExforDataSets             | This is the suite class for the GNDS                                          |
    |                           |     documentation/experimentalDataSets/exforDataSets node.                    |
    +---------------------------+-------------------------------------------------------------------------------+
    | ExperimentalDataSets      | This is the suite class for the GNDS documentation/experimentalDataSets node. |
    +---------------------------+-------------------------------------------------------------------------------+
"""

import datetime

from LUPY import ancestry as ancestryModule

from .. import text as textModule
from .. import suite as suiteModule
from .. import date as dateModule


class CovarianceScript(textModule.Text):
    """
    This class represents a GNDS documentation/covarianceScript node. 
    """

    moniker = 'covarianceScript'


class CorrectionScript(textModule.Text):
    """
    This class represents a GNDS documentation/correctionScript node. 
    """

    moniker = 'correctionScript'


class Note(textModule.Text):
    """
    This is the suite class for the GNDS documentation/note node.
    """

    moniker = 'note'


class ExforDataSet(ancestryModule.AncestryIO):
    """
    This class represents a GNDS documentation/experimentalDataSets/exforDataSets/exforDataSet node.

    The following table list the primary members of this class:

    +-------------------+-------------------------------------------------------------------+
    | Member            | Description                                                       |
    +===================+===================================================================+
    | subentry          | This is the EXFOR subentry for the dataset.                       |
    +-------------------+-------------------------------------------------------------------+
    | retrievalDate     | This is the date the dataset was retrieved.                       |
    +-------------------+-------------------------------------------------------------------+
    | covarianceScript  | Script used to computing covariance matrices based off the        |
    |                   | statistical and systematic errors presented in an EXFOR entry.    |
    +-------------------+-------------------------------------------------------------------+
    | correctionScript  | Script used for data correction of an EXFOR entry.                |
    +-------------------+-------------------------------------------------------------------+
    | note              | A :py:class:`textModule.Text` for storing notes.                  |
    +-------------------+-------------------------------------------------------------------+
    """

    moniker = 'exforDataSet'
    keyName = 'subentry'
    ancestryMembers = ('covarianceScript', 'correctionScript', 'note')

    def __init__(self, _subentry, _retrievalDate):
        """
        :param _subentry:           This is the EXFOR subentry for the dataset.
        :param _retrievalDate:      This is the date the dataset was retrieved.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__subentry = textModule.raiseIfNotString(_subentry, 'subentry')
        self.__retrievalDate = dateModule.raiseIfNotDate(_retrievalDate)
        self.__covarianceScript = CovarianceScript()
        self.__correctionScript = CorrectionScript()
        self.__note = Note()

    @property
    def subentry(self):
        """
        The method returns the subentry.

        :returns:       An instance of :py:class:`textModule.Text`.
        """

        return self.__subentry

    @property
    def covarianceScript(self):
        """
        This method returns a reference to the covarianceScript member.

        :returns:       An instance of :py:class:`CovarianceScript`.
        """

        return self.__covarianceScript

    @property
    def correctionScript(self):
        """
        This method returns a reference to the correctionScript member.

        :returns:       An instance of :py:class:`CorrectionScript`.
        """

        return self.__correctionScript

    @property
    def note(self):
        """
        This method returns a reference to the note member.

        :returns:       An instance of :py:class:`Note`.
        """

        return self.__note

    @property
    def retrievalDate(self):
        """
        This method returns the retrieval date.

        :returns:       An instance of :py:class:`dateModule.Date`.
        """

        return self.__retrievalDate.__str__

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:       The minimum amount of indentation.
        :param kwargs:       dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:             List of str instances representing the XML lines of self.
        """

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XMLList = [
            '%s<%s subentry="%s" retrievalDate="%s">' % (indent, self.moniker, self.__subentry, self.retrievalDate())]
        XMLList += self.__covarianceScript.toXML_strList(indent2, **kwargs)
        XMLList += self.__correctionScript.toXML_strList(indent2, **kwargs)
        XMLList += self.__note.toXML_strList(indent2, **kwargs)
        XMLList[-1] += '</%s>' % self.moniker

        return XMLList

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

        _subentry = node.get('subentry')
        _retrievalDate = dateModule.Date(datetime.datetime.strptime(node.get('retrievalDate'), '%Y-%m-%d'))

        exforDataSet = cls(_subentry, _retrievalDate)
        exforDataSet.parseAncestryMembers(node, xPath, linkData, **kwargs)

        return exforDataSet

class ExforDataSets(suiteModule.Suite):
    """
    This is the suite class for the GNDS documentation/experimentalDataSets/exforDataSets node.
    """

    moniker = 'exforDataSets'
    suiteName = 'subentry'

    def __init__(self):
        suiteModule.Suite.__init__(self, [ExforDataSet])


class ExperimentalDataSets(ancestryModule.AncestryIO_base):
    """
    This is the suite class for the GNDS documentation/experimentalDataSets node. 

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | exforDataSets     | A suite of EXFOR datasets.                                    |
    +-------------------+---------------------------------------------------------------+
    """

    moniker = 'experimentalDataSets'
    ancestryMembers = ('exforDataSets',)

    def __init__(self):
        ancestryModule.AncestryIO_base.__init__(self)

        self.__exforDataSets = ExforDataSets()
        self.__exforDataSets.setAncestor(self)

    @property
    def exforDataSets(self):
        """
        This method returns a reference to the exforDataSets member.

        :returns:       instance of :py:class`ExforDataSets`.
        """

        return self.__exforDataSets

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:       The minimum amount of indentation.
        :param kwargs:       dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:             List of str instances representing the XML lines of self.
        """

        if not len(self.__exforDataSets): return []

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XMLList = ['%s<%s>' % (indent, self.moniker)]
        XMLList += self.__exforDataSets.toXML_strList(indent2, **kwargs)
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

        for child in node:
            if child.tag == ExforDataSets.moniker:
                self.__exforDataSets.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise Exception('Unsupported sub-child "%s" for node "%s".' % (child.tag, node.tag))
