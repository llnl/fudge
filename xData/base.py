# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the xData class that must be the base class for all xData classes.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | XDataBase                         | This class is only inherited by XData so can be deleted.              |
    +-----------------------------------+-----------------------------------------------------------------------+
    | XData                             | This class is only inherited by XDataCoreMembers so can be deleted.   |
    +-----------------------------------+-----------------------------------------------------------------------+
    | XDataCoreMembers                  | This is the base class inherited by many xData class.                 |
    +-----------------------------------+-----------------------------------------------------------------------+
    | XDataFunctional                   |                                                                       |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Function                          | Description                                                           |
    +===================================+=======================================================================+
    | getArguments                      | This function replaces the values in one dictionary with              |
    |                                   | corresponding values in another dictionary.                           |
    +-----------------------------------+-----------------------------------------------------------------------+
    | isXDataFunctional                 | This function returns True if object is an instance of XDataBase and  |
    |                                   | False otherwise.                                                      |
    +-----------------------------------+-----------------------------------------------------------------------+
    | getDomainValue                    | This function returns a number or a default.                          |
    +-----------------------------------+-----------------------------------------------------------------------+
    | getDomainLimits                   | This function returns the limits of a domain.                         |
    +-----------------------------------+-----------------------------------------------------------------------+
    | processUnits                      | This function returns a multiplication or diviision of two units.     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | getDomainValue2                   | This function returns either a float or a :py:class:`PQUModule.PQU`   |
    |                                   | instance depending if its argument has a unit or not.                 |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import abc

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from . import enums as enumsModule

def getArguments(kwargs, arguments):
    """
    This function replaces every value for a key in *arguments* with the value of the key in *kwargs*. If a key in *arguments* 
    is not in *kwargs*, *arguments* value for that key is not changed.

    :param kwargs:          A python dictionry containing some replacement values for the key/value entries of *arguments*.
    :param arguments:       A python dictionry whose values are replace by the values in *kwargs* when their keys appear in *kwargs*.

    :returns:               The modified *arguments*.
    """

    for argument in arguments: arguments[argument] = kwargs.get(argument, arguments[argument])
    return arguments

class XDataBase:
    """This class is only inherited by XData so can be deleted."""

    pass

class XData(XDataBase, ancestryModule.AncestryIO):
    """This class is only inherited by XDataCoreMembers so can be deleted."""


    def __init__(self):

        ancestryModule.AncestryIO.__init__(self)

class XDataCoreMembers(XData):
    """
    This is the base class inherited by many xData class and defined the members index and label.
    
    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | index     | This is the index member use by some xData classes.       |
    +-----------+-----------------------------------------------------------+
    | label     | This is the label member use by some xData classes.       |
    +-----------+-----------------------------------------------------------+
    """

    def __init__(self, index=None, label=None):
        """
        :param index:       This argument defines index member.
        :param label:       This argument defines the label member.
        """

        XData.__init__(self)

        if index is not None: index = int(index)
        self.index = index
        self.label = label

    def attributesToXMLAttributeStr(self):
        """
        This functions returns the XML attributes string for the index and label members.
        """

        attributeStr = ''
        if self.index is not None: attributeStr += ' index="%s"' % self.index
        if self.label is not None: attributeStr += ' label="%s"' % self.label
        return attributeStr

    @property
    def index(self):
        """The method returns the index member."""

        return self.__index

    @index.setter
    def index(self, value):
        """The method sets the index member to *value*."""

        if value is not None: value = int(value)
        self.__index = value

    @property
    def label(self):
        """The method returns the label member."""

        return self.__label

    @label.setter
    def label(self, value):
        """The method sets the label member to *value*."""

        if value is not None:
            if not isinstance(value, str): raise TypeError('label must be a string, got %s' % type(value))
        self.__label = value

    @staticmethod
    def getArguments(kwargs, arguments):
        """
        This function calls :py:func:`getArguments` and returns its results.

        :param kwargs:          A python dictionry containing some replacement values for the key/value entries of *arguments*.
        :param arguments:       A python dictionry whose values are replace by the values in *kwargs* when their keys appear in *kwargs*.

        :returns:               The modified *arguments*.
        """

        return getArguments(kwargs, arguments)

class XDataFunctional(XDataCoreMembers):
    """
    This is the base class inherited by the xData classes that can have an axes, outerDomainValue, index and label members.
    
    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | index             | This is the index member use by some xData classes.           |
    +-------------------+---------------------------------------------------------------+
    | label             | This is the label member use by some xData classes.           |
    +-------------------+---------------------------------------------------------------+
    | axes              | This is the axes member.                                      |
    +-------------------+---------------------------------------------------------------+
    | outerDomainValue  | This is the domain value for the next higher dimension for    |
    |                   | a function that is embedded in a high dimensional function.   |
    +-------------------+---------------------------------------------------------------+
    """

    ancestryMembers = ( 'axes', 'uncertainty' )

    def __init__(self, axes, index=None, valueType=enumsModule.ValueType.float64, outerDomainValue=None, label=None):
        """
        :param axes:                This is the axes member.
        :param valueType:           This describes the type of data (i.e., float, int) returned by the function.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is 
                                    embedded in a high dimensional function.
        :param index:               See documentation for constructor for :py:class:`XDataCoreMembers`.
        :param label:               See documentation for constructor for :py:class:`XDataCoreMembers`.
        """

        from . import axes as axesModule
        from . import uncertainties as uncertaintiesModule

        XDataCoreMembers.__init__(self, index=index, label=label)

        self.valueType = valueType

        if outerDomainValue is not None: outerDomainValue = float(outerDomainValue)
        self.__outerDomainValue = outerDomainValue

        if axes is None: axes = axesModule.Axes(0)
        self.axes = axes
        self.uncertainty = uncertaintiesModule.Uncertainty()                # FIXME2, uncertainty should be a suite like object.
        self.uncertainty.setAncestor(self)

    @property
    @abc.abstractmethod
    def dimension(self):
        """This method must be over-written by the sub-class with the correct dimension (i.e., the number of independent axes) of the function."""

        pass

    def attributesToXMLAttributeStr(self):
        """
        This functions returns the XML attributes string for the index, label, and *outerDomainValue*  members.
        """

        attributeStr = XDataCoreMembers.attributesToXMLAttributeStr(self)
        if self.outerDomainValue is not None: attributeStr += ' outerDomainValue="%s"' % PQUModule.floatToShortestString(self.outerDomainValue, 12)
        return attributeStr

    @property
    def axes(self):
        """This method returns self's axes."""

        return self.__axes

    @axes.setter
    def axes(self, axes1):
        """
        This method sets self's axes to *axes1*.

        :param axes1:       The new axes for *self*.
        """

        from . import axes as axesModule

        if not isinstance(axes1, axesModule.Axes): raise TypeError('axes in not an axes instance')
        if len(axes1) > 0:
            if len(axes1) <= self.dimension:
                raise Exception('len( axes ) = %d != ( self.dimension + 1 ) = %d' % (len(axes1), self.dimension + 1))
        if axes1.ancestor is not None:
            axes1 = axes1.copy()
        axes1.setAncestor(self)
        self.__axes = axes1

    @property
    def outerDomainValue(self):
        """This method returns self's *outerDomainValue* value."""

        return self.__outerDomainValue

    @outerDomainValue.setter
    def outerDomainValue(self, outerDomainValue):
        """
        This method sets self's *outerDomainValue*.

        :param outerDomainValue:        The new *outerDomainValue* value for self.
        """

        if outerDomainValue is not None: outerDomainValue = float(outerDomainValue)
        self.__outerDomainValue = outerDomainValue

    @property
    def valueType(self):
        """This method returns self's *valueType*."""

        return self.__valueType

    @valueType.setter
    def valueType(self, valueType):
        """
        This method sets self's *valueType*.

        :param valueType:       The new *valueType* for self.
        """

        valueType = enumsModule.ValueType.checkEnumOrString(valueType)
        if valueType != enumsModule.ValueType.float64:
            raise TypeError('valueType = "%s" not supported' % valueType)
        self.__valueType = valueType

    def fixValuePerUnitChange(self, factors):
        """
        This methods changes self's *outerDomainValue* per a unit change if self's *outerDomainValue* is not None.
        This method is called from a sub-class's **convertUnits** method.

        :param factors:         A list of at least dimenstion + 2 floats that are the unit scaling factors..
        """

        if self.__outerDomainValue is not None: self.__outerDomainValue *= factors[self.dimension+1]

    def getAxisIndexByIndexOrName(self, indexOrName):
        """
        This method returns the index of self's axis that is associated with *indexOrName*. If *indexOrName* is an int, it is returned.
        Otherwise, *indexOrName* must be a string that matches the label of one of self's axis whose index will be returned.

        :param indexOrName:     The index or label of the axis whose index is requested.

        :returns:               A python int for the index matching *indexOrName*.
        """

        if isinstance(indexOrName, int): return indexOrName
        if isinstance(indexOrName, str):
            for index, axis in enumerate(self.axes):
                if axis.label == indexOrName: return index
        raise TypeError('argument must be an integer or a string')

    def getAxisUnitSafely(self, index):
        """
        This method returns the unit for the axis at *index*. If no axis' are present, an empty string is returned.

        :param index:       The index of the axis whose unit is being requested.

        :returns:           A python str for the unit of the requested axis.
        """

        if len(self.axes) == 0: return ''
        return self.axes[index].unit

    def getPrimaryXData(self):
        """
        Thie method returns a reference to the top level :py:class:`XDataFunctional` function. That is, if *self* is embedded in a 
        higher dimensional :py:class:`XDataFunctional` function, a reference to the outer most function is returned, otherwise, *self* is returned.

        :returns:           A reference to the top level :py:class:`XDataFunctional` function.
        """

        ancestor = self
        while not ancestor.isPrimaryXData(): ancestor = ancestor.ancestor
        return ancestor

    def isPrimaryXData(self):
        """This method returns False if self is contained in a higher dimension :py:class:`XDataFunctional` and False otherwise."""

        ancestry = self.ancestor
        if ancestry is None: return True
        return not isinstance(ancestry, XDataFunctional)

    def buildXML_strList(self, indent, startTag, extraData, **kwargs):
        """
        This methods builds the list of XML strings by adding (if present) the axes, uncertainty data, the end tag name and
        the additional data sub-class adds.  This method is for internal use.

        :param indent:          The minimum amount of indentation.
        :param startTag:        This is the first XML line containing that must contain the start tag and the attributes.
        :param extraData:       This is the extra data as list of XML strings that the sub-class adds.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:               The list of XML strings.
        """

        XML_strList = startTag
        if self.isPrimaryXData():
            if len(self.axes) > 0 or kwargs.get('showEmpty', False):
                XML_strList += self.axes.toXML_strList(indent=indent, **kwargs)
        XML_strList += extraData
        XML_strList += self.uncertainty.toXML_strList(indent=indent, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    def parseNodeStandardChildren(self, node, xPath, linkData, **kwargs):
        """
        This method parses the common child nodes (i.e., 'axes' and 'uncertainty') of *self*, and returns the other child nodes, 
        which were not parsed, as a dictionary with keys of the child.tag values.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Python dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is parsed.

        :returns:           A dictionary of the child nodes not parsed by this method.
        """

        extraNodes = []
        for child in node:
            if child.tag == 'axes':
                self.axes.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == 'uncertainty':
                self.uncertainty.parseNode(child, xPath, linkData, **kwargs)
            else:
                extraNodes.append(child)

        if len(self.axes) == 0 and 'axes' in kwargs: self.axes = kwargs['axes']

        return extraNodes

    @staticmethod
    def parseBareNodeCommonAttributes(node, xPath, allowInterapolation=False):
        """
        This method returns two dictionaries. One dictionary contains all the common attributions of *self*
        and the other dictionary constains the other attributions.

        :param node:                    Node to parse.
        :param xPath:                   List containing xPath to current node, useful mostly for debugging.
        :param allowInterapolation:     If True, an interpolation attribution is added to the list of common attributions.

        :returns:                       Two python dictionaries.
        """

        keyExist = False
        for attrName in ( 'label', 'outerDomainValue' ):
            if node.get(attrName) is not None:
                keyExist = True
                xPath.append('%s[@%s="%s"]' % (node.tag, attrName, node.get(attrName) ))
                break
        if not keyExist: xPath.append(node.tag)

        extraAttributes = {}
        attributes =     { 'label': None, 'index': None, 'outerDomainValue': None  }
        attributeTypes = { 'label': str,  'index': int,  'outerDomainValue': float }
        if allowInterapolation:
            attributes['interpolation'] = enumsModule.Interpolation.linlin
            attributeTypes['interpolation'] = str

        for key, item in list(node.items()):
            if key in attributes:
                attributes[key] = attributeTypes[key](item)
            else:
                extraAttributes[key] = item

        return attributes, extraAttributes

def isXDataFunctional(object):
    """
    This function returns True if object is an instance of XDataBase and False otherwise.

    :param object:      The instance to check.

    :returns:           Python boolean.
    """

    return isinstance(object, XDataBase)

def getDomainValue(value, unit, default):
    """
    This function returns a number which will be *default* if *value* is None. Otherwise, the number is *value* (in unit of *unit* if 
    *value* is a :py:class:`PQUModule.PQU` instance).

    :param value:       This number to return if not None.
    :param unit:        The unit of the returned number.
    :param default:     The default number to return if *value* is None.

    :returns:           A number.
    """

    if value is None: return default
    if isinstance(value, ( str, PQUModule.PQU )): return PQUModule.PQU(value).getValueAs(unit)
    return value

def getDomainLimits(self, domainMin, domainMax, unit):
    """
    This function returns the limits of a domain in unit of *unit*.

    :param domainMin:   The lower limit of the domain.
    :param domainMax:   The upper limit of the domain.
    :param unit:        The unit of the returned domain limits.

    :returns:           Two numbers.
    """

    defaultMin, defaultMax = self.domainMin, self.domainMax
    return getDomainValue(domainMin, unit, defaultMin), getDomainValue(domainMax, unit, defaultMax)

def processUnits(unit1, unit2, operator):
    """
    This function returns a unit this is *unit1* operating on *unit2*. The *operator* can only be "*" or "/".

    :param unit1:       The unit to the left of *operator*.
    :param unit2:       The unit to the right of *operator*.
    :param operator:    The operator.

    :returns:           A python str instance..
    """

    if operator not in [ '*', '/' ]: raise ArithmeticError('unsupported unit operation "%s"' % operator)
    result = eval('PQUModule.PQU(1, unit1) %s PQUModule.PQU(1, unit2)' % operator)
    return result.getUnitSymbol()

def getDomainValue2(domainValue):
    """
    This function returns *domainValue* as either a float if *domainValue* is a number or a :py:class:`PQUModule.PQU` if *domainValue*
    is a number with a unit.

    :domainValue:       A number, :py:class:`PQUModule.PQU` or str representing the number with optional unit to return.

    :returns:           A python float or :py:class:`PQUModule.PQU` instance.
    """

    if isinstance(domainValue, PQUModule.PQU): return domainValue
    if isinstance(domainValue, str):
        try:
            return float(domainValue)
        except:
            return PQUModule.PQU(domainValue)
    return float(domainValue)
