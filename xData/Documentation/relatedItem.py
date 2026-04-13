# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation TBD.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | RelationType              | This enum class represents the allowed relation types.                            |
    +---------------------------+-----------------------------------------------------------------------------------+
    | RelatedItem               | This is the suite class for the GNDS TBD.                                         |
    +---------------------------+-----------------------------------------------------------------------------------+
    | RelatedItems              | This is the suite class for the GNDS TBD.                                         |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule


class RelationType(enumsModule.Enum):
    """
    This enum class represents the allowed relation types. 
    """

    none = enumsModule.auto()
    isCitedBy = 'IsCitedBy'
    cites = 'Cites'
    isSupplementTo = 'IsSupplementTo'
    isSupplementedBy = 'IsSupplementedBy'
    isContinuedBy = 'IsContinuedBy'
    continues = 'Continues'
    describes = 'Describes'
    isDescribedBy = 'IsDescribedBy'
    hasMetadata = 'HasMetadata'
    isMetadataFor = 'IsMetadataFor'
    hasVersion = 'HasVersion'
    isVersionOf = 'IsVersionOf'
    isNewVersionOf = 'IsNewVersionOf'
    isPreviousVersionOf = 'IsPreviousVersionOf'
    isPartOf = 'IsPartOf'
    hasPart = 'HasPart'
    isPublishedIn = 'IsPublishedIn'
    isReferencedBy = 'IsReferencedBy'
    references = 'References'
    isDocumentedBy = 'IsDocumentedBy'
    documents = 'Documents'
    isCompiledBy = 'IsCompiledBy'
    complies = 'Complies'
    isVariantFormOf = 'IsVariantFormOf'
    isOriginalFormOf = 'IsOriginalFormOf'
    isIdenticalTo = 'IsIdenticalTo'
    isReviewedBy = 'IsReviewedBy'
    reviews = 'Reviews'
    isDerivedFrom = 'IsDerivedFrom'
    isSourceOf = 'IsSourceOf'
    isRequiredBy = 'IsRequiredBy'
    requires = 'Requires'
    obsoletes = 'Obsoletes'
    isObsoletedBy = 'IsObsoletedBy'


class RelatedItem(ancestryModule.AncestryIO):
    """
    TBD.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | name          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    | href          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    | relationType  | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'relatedItem'
    keyName = 'name'

    def __init__(self, name, href, relationType):
        """
        :param name:            TBD.
        :param href:            TBD.
        :param relationType:    TBD.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__name = textModule.raiseIfNotString(name, 'name')
        self.__href = textModule.raiseIfNotString(href, 'href')
        self.__relationType = RelationType.checkEnumOrString(relationType)

    @property
    def name(self):
        """
        This method returns the *name* member.

        :returns:       instance of :py:class:`textModule.Text`.
        """

        return self.__name

    @property
    def href(self):
        """
        This method returns the *href* member.

        :returns:       instance of :py:class:`textModule.Text`.
        """

        return self.__href

    @property
    def relationType(self):
        """
        This method returns the *relationType* member.

        :returns:       instance of :py:class:`RelationType`.
        """

        return self.__relationType

    def toXML_strList(self, **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param kwargs:       dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:             List of str instances representing the XML lines of self.
        """

        indent = kwargs.get('indent', '')

        attributes = ' name="%s"' % self.__name
        if len(self.__href) > 0: attributes += ' href="%s"' % self.__href
        if len(self.__relationType) > 0: attributes += ' relationType="%s"' % self.__relationType

        return ['%s<%s%s/>' % (indent, self.moniker, attributes)]

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

        name = node.get('name')
        href = node.get('href', '')
        relationType = node.get('relationType', '')

        return cls(name, href, relationType)


class RelatedItems(suiteModule.Suite):
    """
    This is the suite class for the GNDS TBD.
    """

    moniker = 'relatedItems'
    suiteName = 'name'

    def __init__(self):
        suiteModule.Suite.__init__(self, [RelatedItem])
