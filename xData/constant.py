# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the xData classes that represent functions that are constant over their independent axes.

This module contains the following classes:

    +------------------------+--------------------------------------------------------------------------+
    | Class                  | Description                                                              |
    +========================+==========================================================================+
    | Constant               | This class is the base class for all constant classes.                   |
    +------------------------+--------------------------------------------------------------------------+
    | Constant1d             | This class represents a 1-d function that is constant over its domain.   |
    +------------------------+--------------------------------------------------------------------------+
"""

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from pqu import PQU as PQUModule

from . import enums as enumsModule
from . import base as baseModule
from . import XYs1d as XYs1dModule

class Constant( baseModule.XDataFunctional )  :
    """
    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | value             | This is the value of the function at its points.              |
    +-------------------+---------------------------------------------------------------+
    | domainMin         | This is the lower domain value for which the function defined.|
    +-------------------+---------------------------------------------------------------+
    | domainMax         | This is the upper domain value for which the function defined.|
    +-------------------+---------------------------------------------------------------+
    | axes              | This is the axes member.                                      |
    +-------------------+---------------------------------------------------------------+
    | outerDomainValue  | This is the domain value for the next higher dimension for    |
    |                   | a function that is embedded in a high dimensional functions.  |
    +-------------------+---------------------------------------------------------------+
    | index             | This is the index member use by some xData classes.           |
    +-------------------+---------------------------------------------------------------+
    | label             | This is the label member use by some xData classes.           |
    +-------------------+---------------------------------------------------------------+
    """

    def __init__(self, _value, domainMin, domainMax, axes=None, label=None, index=None, outerDomainValue=None):
        """
        :param _value:              This is the value of the function within its defined domain.
        :param domainMin:           This is the lower domain value for which the function defined.
        :param domainMax:           This is the upper domain value for which the function defined.
        :param axes:                This is the axes member.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is 
                                    embedded in a high dimensional function.
        :param index:               This is the index member.
        :param label:               This is the label member.
        """

        baseModule.XDataFunctional.__init__(self, label=label, axes=axes)

        self.value = _value

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : TypeError( 'domainMin not a float instance' )
        self.__domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : TypeError( 'domainMax not a float instance' )
        self.__domainMax = domainMax

    def copy( self ) :
        """
        This method returns a new instance that is a copy of *self*.

        :returns:           An instance the same type as *self*.
        """

        axes = self.axes
        if( axes is not None ) : axes = self.axes.copy( )
        return( self.__class__( self.value, self.domainMin, self.domainMax, axes = axes, label = self.label ) )

    __copy__ = copy

    @property
    def value( self ) :
        """
        This method returns *self*'s value.

        :returns:           An instance of type float.
        """

        return( self.__value )

    @value.setter
    def value( self, _value ) :
        """
        This method sets *self*'s value to *_value*.

        :param _value:  The new value for *self*.
        """

        if( isinstance( _value, int ) ) : _value = float( _value )
        if( not( isinstance( _value, float ) ) ) : TypeError( 'value not a float instance' )
        self.__value = _value

    @property
    def domainMin( self ) :
        """
        This method returns the minimum domain value for *self*.

        :returns:       A float.
        """

        return( self.__domainMin )

    @domainMin.setter
    def domainMin(self, domainMin):
        """
        This method sets self's domainMin to *domainMin*. *domainMin* must be less than self.domainMax or a **raise** is executed.

        :param domainMin:   The new lower limit for the domain.
        """

        if domainMin >= self.domainMax:
            raise ValueError('domainMin = %.17 >= self.domainMax = %.17e' % (domainMin, self.domainMax))
        self.__domainMin = domainMin

    @property
    def domainMax( self ) :
        """
        This method returns the maximum domain value for *self*.

        :returns:       A float.
        """

        return( self.__domainMax )

    @domainMax.setter
    def domainMax(self, domainMax):
        """
        Thie method sets *self*'s domainMax to *domainMax*. *domainMax* must be greater than self.domainMin or a **raise** is executed.

        :param domainMin:   The new lower limit for the domain.
        """

        if domainMax <= self.domainMin:
            raise ValueError('domainMax = %.17 <= self.domainMin = %.17e' % (domainMax, self.domainMin))
        self.__domainMax = domainMax

    @property
    def domainGrid(self):
        """
        This method returns *self*'s *domainMin* and *domainMax* values.

        :returns:           A python tuple.
        """

        return self.__domainMin, self.__domainMax

    @property
    def domainUnit( self ) :
        """
        This method returns the domain unit for *self*.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :
        """
        This method always returns *self*'s value.

        :returns:       A python float.
        """

        return( self.__value )

    rangeMax = rangeMin

    @property
    def rangeUnit( self ) :
        """
        This method returns the unit for the dependent variable.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( 0 ) )

    def fixDomainPerUnitChange( self, factors ) :
        """
        This method multiplies *domainMin* and *domainMax* by factors[self.dimension].

        :param factors:     This is a list of scaling factors and must contain at least two floats.
        """

        self.__domainMin *= factors[self.dimension]
        self.__domainMax *= factors[self.dimension]

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        This method sets *domainMin* and *domainMax* per the arguments.

        :param domainMin:       The lower limit of the domain.
        :param domainMax:       The upper limit of the domain.
        :param fixToDomain:     An instance of :py:class:`enumsModule.FixDomain` that specifies which limits are to be fixed.

        :returns:               This method returns 0 if no domain limit was moved and 1 if at least one was moved.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        domainMin = max(domainMin, self.domainMin)
        domainMax = min(domainMax, self.domainMax)
        if fixToDomain == enumsModule.FixDomain.lower:
            self.__domainMin = domainMin
        elif fixToDomain == enumsModule.FixDomain.upper:
            self.__domainMax = domainMax
        else:
            self.__domainMin = domainMin
            self.__domainMax = domainMax

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

class Constant1d( Constant ) :
    """
    This class represents a 1d-function that has a constant value over its domain.
    """

    moniker = 'constant1d'
    dimension = 1

    def __truediv__(self, other):
        """
        This method returns a :py:class:`Constant1d` that is *self* divided by *other* where *other* must be
        a python int of float.

        :param other:       A python int or float.

        :returns:           An instance of :py:class:`Constant1d`.
        """

        if isinstance(other, (int, float)):
            return Constant1d(self.value / other, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __rtruediv__(self, other):
        """
        This method returns a :py:class:`Constant1d` that is *other* divided by *self* where *other* must be
        a python int of float.

        :param other:       A python int or float.

        :returns:           An instance of :py:class:`Constant1d`>
        """

        if isinstance(other, (int, float)):
            return Constant1d(other / self.value, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __itruediv__(self, other):
        """
        This method divides *self* by *other* where *other* must be a python int of float.

        :param other:       A python int or float.

        :returns:           An instance of :py:class:`Constant1d` which is *self*.
        """

        if isinstance(other, (int, float)):
            self.value /= other
            return self
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __mul__(self, other):
        """
        This method returns a :py:class:`Constant1d` that is *self* multiplied by *other* where *other* must be
        a python int of float.

        :param other:       A python int or float.

        :returns:           An instance of :py:class:`Constant1d`.
        """

        if isinstance(other, (int, float)):
            return Constant1d(self.value * other, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Multiplying Constant1d by {type(other)}")

    __rmul__ = __mul__

    def __imul__(self, other):
        """
        This method multiplies *self* by *other* where *other* must be a python int of float.

        :param other:       A python int or float.

        :returns:           An instance of :py:class:`Constant1d` which is *self*.
        """

        if isinstance(other, (int, float)):
            self.value *= other
            return self
        else:
            raise NotImplementedError(f"Multiplying Constant1d by {type(other)}")

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance. 

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps:    Used to dull the lower point for "flat" interpolation.
        :param upperEps:    Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self.toPointwise_withLinearXYs(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        return xys1d

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if len(self.axes) == 0: return

        factors = self.axes.convertUnits( unitMap )
        self.value *= factors[0]
        self.fixDomainPerUnitChange( factors )
        self.fixValuePerUnitChange( factors )

    def convertAxisToUnit( self, indexOrName, newUnit ) :
        """
        This method converts the axis at *indexOrName* to unit *newUnit* and returns a :py:class:`Constant1d` instance with 
        that axis data scaled to the new unit.

        :param indexOrName:     The index or name of the axis to convert to the new unit.
        :param newUnit:         New unit for axis *indexOrName*.

        :returns:               An instance of :py:class:`Constant1d`.
        """

        if len(self.axes) == 0:
            n = self
        else:
            n = self.copy()
            index = n.getAxisIndexByIndexOrName( indexOrName )
            axis = n.axes[index]
            factor = PQUModule.PQU( '1 ' + axis.unit ).getValueAs( newUnit )
            if index == 0:
                n.fixDomainPerUnitChange({0: factor})
            elif index == 1:
                n.value *= factor
            axis.unit = newUnit
        return n

    def evaluate( self, x ) :
        """
        This method returns the value of *self* at the point *x*.

        :param x:           The domain point where *self* is evaluated.

        :returns:           A python float.
        """

        return( self.value )

    def toPointwise_withLinearXYs(self, **kwargs):
        """
        This method returns an :py:class:`XYs1dModule.XYs1d` representation of *self*.

        :param kwargs:      Not used but present to be compatible with other similar methods.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        xys1d = XYs1dModule.XYs1d(data=[[self.domainMin, self.value], [self.domainMax, self.value]],
                      axes=self.axes, outerDomainValue=self.outerDomainValue, label=self.label)
        return xys1d

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        valueFormatter = kwargs.get('valueFormatter', floatToShortestString)
        significantDigits = kwargs.get('significantDigits', 15)

        attributeStr = baseModule.XDataCoreMembers.attributesToXMLAttributeStr(self)
        startTag = [ '%s<%s%s value="%s" domainMin="%s" domainMax="%s">' % ( indent, self.moniker, attributeStr, 
                valueFormatter(self.value, significantDigits = significantDigits),
                valueFormatter(self.domainMin, significantDigits = significantDigits),
                valueFormatter(self.domainMax, significantDigits = significantDigits) ) ]

        return self.buildXML_strList(indent2, startTag, [], **kwargs)

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:            An instance of *cls* representing *node*.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.

        value = float(extraAttributes.pop('value'))
        domainMin = float(extraAttributes.pop('domainMin'))
        domainMax = float(extraAttributes.pop('domainMax'))

        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        instance = cls(value, domainMin, domainMax, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(instance, node, xPath, linkData, **kwargs)
        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()                                     # Per comment above, parseBareNodeCommonAttributes adds to xPath.

        return instance
