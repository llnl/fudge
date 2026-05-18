# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes and functions for storing a date/time per GNDS specifications. A data/time can be in one of the 
following three forms::

    -) date only (e.g., "2021-01-21"),
    -) date and time (e.g., "2021-01-21T08:42:04") or
    -) date and time with UTC offset (e.g., "2021-01-21T08:42:04-08:00").

This module contains the following classes:

    +-----------------------------------+------------------------------------------------------------------------+
    | Class                             | Description                                                            |
    +===================================+========================================================================+
    | Resolution                        | This class is an enum of the supported resolutions for a data/time.    |
    +-----------------------------------+------------------------------------------------------------------------+
    | Date                              | This class represents a GNDS date/time.                                |
    +-----------------------------------+------------------------------------------------------------------------+

This module contains the following functions:

    +-----------------------------------+------------------------------------------------------------------------+
    | Function                          | Description                                                            |
    +===================================+========================================================================+
    | raiseIfNotDate                    | This function checks the type of *date* and executes a raise if it     |
    |                                   | is not a :py:class:`Date` instance.                                    |
    +-----------------------------------+------------------------------------------------------------------------+
    | raiseIfNotPythonDateTime          | This function checks the type of *_datetime* and executes a raise      |
    |                                   | if it is not a :py:class:`datetime.datetime` instance.  If *_datetime* |
    |                                   | is an instance of :py:class:`datetime.datetime`, a :py:class:`Date`    |
    |                                   | representation of it is returned.                                      |
    +-----------------------------------+------------------------------------------------------------------------+
"""

import sys
import datetime

from LUPY import enums as enumsModule

class Resolution(enumsModule.Enum):
    """
    This class is an enum of the supported resolutions for a data/time.
    """

    date = enumsModule.auto()
    time = enumsModule.auto()
    undefined = enumsModule.auto()

class Date :
    """
    The class supports a date/time in one of following three3 forms::

        -) date only (e.g., "2021-01-21"),
        -) date and time (e.g., "2021-01-21T08:42:04") or
        -) date and time with UTC offset (e.g., "2021-01-21T08:42:04-08:00").

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | datetime      | This is an instance of :py:class:`datetime`.                  |
    +---------------+---------------------------------------------------------------+
    | resolution    | This is an instance of :py:class:`Resolution`.                |
    +---------------+---------------------------------------------------------------+

    If *resolution* is **undefined** then *datetime* is None. If *resolution* is **date**, then
    the time part of *datetime* is set to "00:00:00".
    """

    def __init__( self, _datetime = None, resolution = Resolution.date, UTC = False ) :
        """
        If *resolution* is "undefied", *datetime* is set to None, independent of *_datetime*.
        If *_datetime* is None, the current date/time is used.
        If *resolution* is **date**, then the time part of *datetime* is set to "00:00:00".
        The microsecond part of *datetime* is always set to 0.

        :param _datetime:   User supplied date/time.
        :param resolution:  An instance of :py:class:`Resolution`.
        :param UTC:         If True, then the date/time has a UTC offset.
        """

        self.__resolution = Resolution.checkEnumOrString(resolution)

        self.__datetime = None
        if( resolution == Resolution.undefined ) : return

        if( _datetime is None ) :
            _datetime = datetime.datetime.now( ).replace( microsecond = 0 )
            if( resolution == Resolution.date ) :
                _datetime = _datetime.replace( hour = 0, minute = 0, second = 0, microsecond = 0 )
            elif( ( sys.version_info.major == 3 ) and UTC ) :
                _datetime = _datetime.astimezone( )

        if( not( isinstance( _datetime, datetime.datetime ) ) ) : raise TypeError( 'Invalid _datetime instance.' )

        self.__datetime = _datetime

    def __str__( self ) :
        """
        The method converts the interal data/time into a python str to the specified *resolution*.
        If *resolution* is "undefined", and empty string is returned.

        :returns:       A python str.
        """

        if( self.__resolution == Resolution.undefined ) : return( '' )

        stringDate = self.__datetime.isoformat( sep = 'T', timespec = 'seconds' )
        if( self.__resolution == Resolution.date ) : stringDate = stringDate.split( 'T' )[0]

        return( stringDate )

    @property
    def date( self ) :
        """
        This method returns a reference to *self*'s *datatime* memeber.

        :returns:       A :py:class:`datetime.datetime` instance or None.
        """

        return( self.__datetime )

    @property
    def resolution( self ) :
        """
        This method returns a reference to *self*'s *resolution* memeber.

        :returns:       A :py:class:`Resolution` instance.
        """

        return( self.__resolution )

    def asXML_attribute( self, name = 'date' ) :

        if( self.resolution == Resolution.undefined ) : return( '' )
        return( ' %s="%s"' % ( name, self ) )

    @staticmethod
    def parse( stringDate ) :
        """
        This method converts a python str into a :py:class:`Date` instance.

        :param stringDate:      A data/time str.

        :returns:               An instance of :py:class:`Date`.
        """

        resolution = Resolution.date
        datetime1 = None

        if( len( stringDate ) == 0 ) :
            resolution = Resolution.undefined
        elif( len( stringDate ) == 10 ) :
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%d" )
        elif( len( stringDate ) == 19 ) :
            resolution = Resolution.time
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%dT%H:%M:%S" )
        else :
            resolution = Resolution.time
            datetime1 = datetime.datetime.strptime( stringDate, "%Y-%m-%dT%H:%M:%S%z" )

        return( Date( datetime1, resolution = resolution ) )

    @staticmethod
    def fromYearMonthDay( year, month, day ) :
        """
        This method converts *year*, *month* and *day* into a :py:class:`Date` instance.

        :param year:        The year for the date.
        :param month:       The month for the date.
        :param day:         The day for the date.

        :returns:       An instance of :py:class:`Date`.
        """

        return( Date( datetime.datetime( year, month, day ) ) )

def raiseIfNotDate(date):
    """
    This function checks the type of *date* and executes a raise if it is not a :py:class:`Date` instance.

    :param date:            The instance to check.

    :returns:               The instance *date*.
    """

    if not isinstance(date, Date): raise TypeError('Invalid date instance.')

    return date

def raiseIfNotPythonDateTime(_datetime):
    """
    This function checks the type of *_datetime* and executes a raise if it is not a :py:class:`datetime.datetime` instance.
    If *_datetime* is an instance of :py:class:`datetime.datetime`, a :py:class:`Date` representation of it is returned.

    :param _datetime:       The instance to check.

    :returns:               An instance of :py:class:`Date`.
    """

    if not isinstance(_datetime,datetime.datetime): raise TypeError('Invalid python datetime instance.')

    return( Date(_datetime) )
