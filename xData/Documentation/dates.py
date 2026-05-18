# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes date and dates.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | DateType                  | This enum class represents the allowed data types.                                |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Date                      | This is the suite class for the GNDS documentation/dates/date node.               |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Dates                     | This is the suite class for the GNDS documentation/dates node.                    |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import date as dateModule


class DateType(enumsModule.Enum):
    """
    This enum class represents the allowed data types.
    """

    accepted = enumsModule.auto()
    available = enumsModule.auto()
    copyrighted = enumsModule.auto()
    collected = enumsModule.auto()
    created = enumsModule.auto()
    issued = enumsModule.auto()
    submitted = enumsModule.auto()
    updated = enumsModule.auto()
    valid = enumsModule.auto()
    withdraw = enumsModule.auto()


class Date(ancestryModule.AncestryIO):
    """
    This is the suite class for the GNDS documentation/dates/date node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | start         | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    | dateType      | The type of date.                                             |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'date'
    keyName = 'dateType'

    def __init__(self, start, dateType):
        """
        :param start:           TBD.
        :param dateType:        An instance of :py:class:`DateType`.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__start = dateModule.raiseIfNotDate(start)
        self.__dateType = DateType.checkEnumOrString(dateType)

    @property
    def start(self):
        """
        FIXME: What is this?
        """

        return self.__start

    value = start

    @property
    def dateType(self):
        """
        This method returns the data type.

        :returns:       instance of :py:class:`DateType`.
        """

        return self.__dateType

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        return ' %s dateType="%s"' % (self.start.asXML_attribute(name='value'), self.dateType)

    def toXML_strList(self, **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param kwargs:       dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:             List of str instances representing the XML lines of self.
        """

        indent = kwargs.get('indent', '')

        return [
            '%s<%s%s dateType="%s"/>' % (indent, self.moniker, self.start.asXML_attribute(name='value'), self.dateType)]

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

        start = dateModule.Date.parse(node.get('value'))  # FIXME2, what is this? How is Date to have a date?
        dateType = node.get('dateType')

        return cls(start, dateType)


class Dates(suiteModule.Suite):
    """
    This is the suite class for the GNDS documentation/dates node.
    """

    moniker = 'dates'
    suiteName = 'dateType'

    def __init__(self):
        suiteModule.Suite.__init__(self, [Date])
