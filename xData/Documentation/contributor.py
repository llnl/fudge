# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation contributors node and its child nodes.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | ContributorType           | This enum class represents the allowed contributor types.                         |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Contributor               | This is the class for the GNDS documentation/contributors/contributor node.       |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Contributors              | This is the suite class for the GNDS documentation/contributors node.             |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import enums as enumsModule

from .. import suite as suiteModule
from . import abstractClasses as abstractClassesModule


class ContributorType(enumsModule.Enum):
    """
    This enum class represents the allowed contributor types.
    """

    contactPerson = 'ContactPerson'
    dataCollector = 'DataCollector'
    dataCurator = 'DataCurator'
    dataManager = 'DataManager'
    distributor = 'Distributor'
    editor = 'Editor'
    hostingInstitution = 'HostingInstitution'
    producer = 'Producer'
    projectLeader = 'ProjectLeader'
    projectManager = 'ProjectManager'
    projectMember = 'ProjectMember'
    registrationAgency = 'RegistrationAgency'
    registrationAuthority = 'RegistrationAuthority'
    relatedPerson = 'RelatedPerson'
    researcher = 'Researcher'
    researchGroup = 'ResearchGroup'
    rightsHolder = 'RightsHolder'
    sponsor = 'Sponsor'
    supervisor = 'Supervisor'
    workPackageLeader = 'WorkPackageLeader'
    other = 'Other'


class Contributor(abstractClassesModule.AuthorAbstract):
    """
    This is the class for the GNDS documentation/contributors/contributor node.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | name              | The author's name.                                            |
    +-------------------+---------------------------------------------------------------+
    | orcid             | The orcid of the author.                                      |
    +-------------------+---------------------------------------------------------------+
    | email             | The email address of the author.                              |
    +-------------------+---------------------------------------------------------------+
    | contributorType   | The contribution type.                                        |
    +-------------------+---------------------------------------------------------------+
    """

    moniker = 'contributor'
    keyName = 'name'

    def __init__(self, name, orcid, email, contributorType):
        """
        :param name:                The author's name.
        :param orcid:               The orcid of the author.
        :param email:               The email address of the author.
        :param contributorType:     The contribution type.
        """

        abstractClassesModule.AuthorAbstract.__init__(self, name, orcid, email)

        self.__contributorType = ContributorType.checkEnumOrString(contributorType)

    @property
    def contributorType(self):
        """
        This method returns the contributorType.

        :returns:       An instance of :py:class:`ContributorType`.
        """

        return self.__contributorType

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        attributes = abstractClassesModule.AuthorAbstract.XML_extraAttributes(self, **kwargs)
        attributes += ' contributorType="%s"' % self.__contributorType
        return attributes

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
        orcid = node.get('orcid', '')
        email = node.get('email', '')
        contributorType = node.get('contributorType', '')

        return cls(name, orcid, email, contributorType)


class Contributors(suiteModule.Suite):
    """
    This is the suite class for the GNDS documentation/contributors node.
    """

    moniker = 'contributors'
    suiteName = 'name'

    def __init__(self):
        suiteModule.Suite.__init__(self, [Contributor])
