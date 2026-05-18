# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes for handling a GNDS values container.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Values                            | This class represents a GNDS values container.                        |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import copy

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from numericalFunctions import listOfDoubles_C as listOfDoubles_CModule

from . import enums as enumsModule
from . import base as baseModule

class Values(baseModule.XDataCoreMembers):
    """
    This class is used to store a list of numbers. The type of the numbers (e.g., floats, ints) stored is specified by 
    the *valueType* member.  The list of all numbers the instance represents may be greater than the number of values 
    in the *values* member as leading and trailing zero values do not have to be stored in the *values* member. 
    The *start* member specifies the index of the first value in the *values* member with all values before *start* being zero. 
    The *length* member specifies the total number of values the instance represents, including leading and trailing zero values.

    The following table list the primary members of this class:

    +---------------+-----------------------------------------------------------------------------------+
    | Member        | Description                                                                       |
    +===============+===================================================================================+
    | start         | The index of the data represented by the first data value in the values member.   |
    +---------------+-----------------------------------------------------------------------------------+
    | length        | The total number of data values the instance represents, including leading and    |
    |               | trailing zero values that are not stored in the *values* member.                  |
    +---------------+-----------------------------------------------------------------------------------+
    | values        | A list of the values not including the first *start* values which are 0.0 or      |
    |               | trailing 0.0 values.                                                              |
    +---------------+-----------------------------------------------------------------------------------+
    | valueType     | The type of data stored in the values member.                                     |
    +---------------+-----------------------------------------------------------------------------------+
    """

    moniker = 'values'

    def __init__(self, _values, valueType=enumsModule.ValueType.float64, start=0, length=None, index=None, label=None):
        """
        :param _values:         The list of the values to store.
        :param valueType:       The type of the data stored in the values member.
        :param start:           The index of the data represented by the first data value in the values member.
        :param length:          The total number of data values the instance represents, including leading and trailing zero values that are not stored in the *values* member.
        :param index:           This is not used.
        :param label:           The label of the values instance.
        """

        baseModule.XDataCoreMembers.__init__(self, index=index, label=label)

        self.__valueType = enumsModule.ValueType.checkEnumOrString(valueType)

        self.__start = int(start)
        numberOfValues = len(_values)
        if length is not None:
            length = int(length)
            if length < self.__start + numberOfValues: raise ValueError('length = %d < start + numberOfValues = %d' % ( length, self.__start + numberOfValues ))
            self.__length = length
        else:
            self.__length = self.__start + numberOfValues

        self.__values = _values

    def __len__(self):
        """
        This method returns the number of stored values in *self*. The actual number of values that *self* represents is returned by the 
        *length* method which includes leading and trailing zeros not stored by *self*.

        :returns:           A python int.
        """

        return len(self.__values)

    def __getitem__(self, index):
        """
        This returns the stored value a index *index* - 1. Note, this will not be the value represented by *self* if *start* is not 0.

        :param index:   An index into the stored numbers.

        :returns:       A number of type *valueType*.
        """

        return self.__values[index]

    def copy(self, untrimZeros = False):
        """
        Returns a :py:class:`Values` instance that is a copy of *self*.

        :returns:       An instance of :py:class:`Values`.
        """

        start = self.start
        length = self.length
        _values = copy.copy(self.__values)
        if untrimZeros and ( self.start != 0 or self.end != self.length ):
            start = 0
            length = None
            _values = self.start * [ 0 ] + self.__values + ( self.length - self.end ) * [ 0 ]

        return Values(_values, self.valueType, start = start, length = length, label = self.label)

    __copy__ = copy

    @property
    def end(self):
        """
        This method returns the number of values not including the trailing zeros that are not stored. Ergo, returns *start* + the number of stored values.

        :returns:       A python int.
        """

        return self.__start + len(self)

    @property
    def length(self):
        """
        This method returns the actual number of values represented by *self* which includes leading and trailing zeros not stored by *self*.

        :returns:       A python int.
        """

        return self.__length

    @property
    def start(self):
        """
        Returns the index represent by the first value of *self*. All values before *start* are zeros and are not store by *self*.

        :returns:       A python int.
        """

        return self.__start

    @property
    def values(self):
        """
        This method returns a reference to *self*'s values member.

        :returns:       A python list.
        """

        return self.__values

    @values.setter
    def values(self, _values):
        """
        This member sets the *values* member of *self* to *_values*

        :param _values:     The new list of values for *self*.
        """

        if self.valueType == enumsModule.ValueType.float64:
            checker = float
        elif self.valueType == enumsModule.ValueType.float32:
            checker = float
        else:
            checker = int

        self.__length = len(_values)
        self.__values = [ checker(value) for value in _values ]

    @property
    def valueType(self):
        """
        This method returns the type of numbers stored in *self*.
        """

        return self.__valueType


    def offsetScaleValues(self, offset, scale):
        """
        This method multiplies each value of *self* by *scale* and then adds *offset* to the value.

        :param offset:      The offset for each value.
        :param scale:       The multiplying factor for each value.
        """

        for i1, value in enumerate(self.__values): self.__values[i1] = value * scale + offset

    def toString(self):
        """
        This method returns a simple string representation of the stored values of *self*.

        :returns:       A python str.
        """

        strList = [ "%s" % value for value in self.__values ]
        return ' '.join(strList)

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        def intValueFormatter(value, significantDigits = 0):     # Ignores significantDigits.

            return str(value)

        formatVersion = kwargs.get('formatVersion')
        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        valuesPerLine = kwargs.get('valuesPerLine', 100)
        valueFormatter = kwargs.get('valueFormatter', floatToShortestString)
        significantDigits = kwargs.get('significantDigits', 15)
        outline = kwargs.get('outline', False)
        if len(self) < 14: outline = False

        strList, line = [], None
        attributeStr = ''

        if self.start != 0: attributeStr += ' start="%d"' % self.start
        if self.end != self.length: attributeStr += ' length="%d"' % self.length

        if self.valueType != enumsModule.ValueType.float64: attributeStr += ' valueType="%s"' % self.valueType
        if self.valueType == enumsModule.ValueType.integer32: valueFormatter = intValueFormatter
        attributeStr += baseModule.XDataCoreMembers.attributesToXMLAttributeStr(self)

        HDF_opts = kwargs.get("HDF_opts")
        if HDF_opts is not None and len(self) >= HDF_opts['minLength']:

            startIndexName = 'startIndex'
            if formatVersion == GNDS_formatVersionModule.version_1_10: startIndexName = 'offset'
            if HDF_opts['flatten']:
                if self.valueType == enumsModule.ValueType.integer32:
                    datasetName = 'iData'
                else:
                    datasetName = 'dData'
                flatData = HDF_opts[datasetName]
                startLength = len(flatData)
                flatData.extend(list(self))
                XML_strList = ['%s<%s href="HDF#/%s" %s="%d" count="%d"%s/>' %
                        (indent, self.moniker, datasetName, startIndexName, startLength, len(self), attributeStr)]
            else:
                datasetName = "data%d" % HDF_opts['index']
                HDF_opts['index'] += 1
                HDF_opts['h5file'].create_dataset(datasetName, data=list(self), **HDF_opts['compression'])
                XML_strList = ['%s<%s href="HDF#/%s"%s/>' % (indent, self.moniker, datasetName, attributeStr)]
        else:
            XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
            if outline:                     # Logic above guarantees more than 14 elements in self.
                line = []
                for i1 in range(6): line.append(valueFormatter(self[i1], significantDigits = significantDigits))
                line.append(' ... ')
                for i1 in range(-6, 0): line.append(valueFormatter(self[i1], significantDigits = significantDigits))
                strList.append(' '.join(line))
            else:
                dataToString = kwargs.get('dataToString', None)
                if dataToString is None:
                    for index, value in enumerate(self.__values):
                        if ( index % valuesPerLine ) == 0:
                            if line is not None: strList.append(' '.join(line))
                            line = []
                        line.append(valueFormatter(value, significantDigits = significantDigits))
                    if line is not None and len(line) > 0: strList.append(' '.join(line))
                else:
                    kwargs['valueFormatter'] = valueFormatter
                    XML_strList += kwargs['dataToString'](self, kwargs['dataToStringParent'], indent = indent2, **kwargs)
            if len(strList) > 0: XML_strList[-1] += strList.pop(0)

            for line in strList: XML_strList.append("%s%s" % ( indent2, line ))
            XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        attrs = {      'valueType': enumsModule.ValueType.float64,      'start': 0,   'length': None, 'label': None }
        attributes = { 'valueType': str,                                'start': int, 'length': int,  'label': str }
        extraAttributes = {}
        for key, item in list(node.items()):
            if key in attributes: 
                attrs[key] = attributes[key](item)
            else:
                extraAttributes[key] = item

        if 'href' in extraAttributes:
            if 'offset' in extraAttributes:
                startIndex = int(extraAttributes.get('offset'))
            else:
                startIndex = int(extraAttributes.get('startIndex', 0))
            if 'startIndex' in extraAttributes: extraAttributes.pop('startIndex')

            count = -1
            if 'count' in extraAttributes: count = int(extraAttributes.pop('count'))

            HDF = linkData['HDF']
            datasetName = extraAttributes.pop('href').split('#')[-1]
            if not datasetName in HDF: HDF[datasetName] = HDF['h5File'].get(datasetName)
            data = HDF[datasetName]
            if count == -1: count = len(data)

            return Values(data[startIndex:startIndex+count].tolist(), **attrs)

        else:
            if len(extraAttributes) > 0: raise TypeError('Invalid attributes:"%s"' % ' '.join([ str(key) for key in extraAttributes ]))

        if node.text is None: node.text = ''
        if attrs['valueType'] == enumsModule.ValueType.integer32.value:
            values1 = map(int, node.text.split())
        else:
            values1, extraCharacters = listOfDoubles_CModule.createFromString(node.text, sep=" ")
            if len(extraCharacters) > 0: raise ValueError('Extra character at the of string of number: "%s".' % extraCharacters[:40])

        return Values(list(values1), **attrs)

    @staticmethod
    def sortAndThin(_values, rel_tol = 0.0, abs_tol = 0.0):
        """
        This function sorts the list of float in *_values* into an ascending order and then thins the list. 
        If there are 2 or fewer points in _values, nothing is done.

        :param rel_tol:     Relative tolerance used to determine if two consecutive points are too close.
        :param abs_tol:     Abolute tolerance used to determine if two consecutive points are too close.

        :returns:       A pythoh list of the sorted, thinned values.
        """

        __values = [ float(value) for value in _values ]
        __values.sort()

        length = len(__values)
        ___values = __values[:1]

        if length < 2: return ___values

        x1 = ___values[0]
        for i1 in range(1, length):
            x2 = __values[i1]
            deltaMax = max(abs_tol, rel_tol * min(abs(x1), abs(x2)))
            if x2 - x1 <= deltaMax: continue
            ___values.append(x2)
            x1 = x2

        ___values[-1] = __values[-1]                # If last point of __values did not get added, replace last point with it, otherwise this does nothing.
            
        return ___values

    @staticmethod
    def valuesWithTrimmedZeros(values1, valueType=enumsModule.ValueType.float64):
        """
        This method determines the *start* and *end* of the non-zero data in *values1* and returns a :py:class`Values` instance with
        *start* and *length* set accordingly.

        :returns:       A :py:class`Values` instance.
        """

        length = len(values1)
        start, end = 0, length - 1
        while end >= start:
            if values1[end] != 0: break
            end -= 1
        end += 1
        if end != 0:
            for value in values1:
                if value != 0: break
                start += 1
        if start != 0 or end != length: values1 = values1[start:end]
        return Values(values1, valueType = valueType, start = start, length = length)
