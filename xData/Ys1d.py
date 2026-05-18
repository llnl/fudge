# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the :py:class:`Ys1d` class, used to store union-grid functions such as gridded cross sections.
When multiple functions can be mapped onto the same independent axis grid, the :py:class:`Ys1d` class can saves 
space by linking to the common independent grid (i.e., x-values) rather than repeating them.
Otherwise, Ys1d functions should behave like XYs1d.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Ys1d                              | This class represents a 1d function as a list of y-value.             |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

from pqu import PQU as PQUModule

from . import enums as enumsModule
from . import base as baseModule
from . import axes as axesModule
from . import values as valuesModule
from . import XYs1d as XYs1dModule


class Ys1d(baseModule.XDataFunctional):
    r"""
    This class presents a 1d function :math:`y(x)` as a list of y-values with the x-values specified by a :py:class:`axesModule.Grid` axis.
    Also an interpolation rule is used to determine :math:`y(x)` between to consecutive x-values. This class has methods
    for adding, subtracting, multipling and dividing a :py:class:`XYs1d` by a number or another :py:class:`XYs1d` instance,
    as well as other operations common to a 1d function.

    In addition to the members in the inherited class :py:class:`baseModule.XDataFunctional`, 
    the following table list the primary members of this class:

    +---------------+-----------------------------------------------------------------------------------+
    | Member        | Description                                                                       |
    +===============+===================================================================================+
    | Ys            | This is the list of :math:`y_i` points representing the y-value of the function.  |
    +---------------+-----------------------------------------------------------------------------------+
    | interpolation | The interpolation rule used to determine :math:`y(x)` between to consecutive      |
    |               | x-values.                                                                         |
    +---------------+-----------------------------------------------------------------------------------+
    """

    moniker = 'Ys1d'
    dimension = 1

    ancestryMembers = ('Ys',)

    def __init__(self, Ys=None, interpolation=enumsModule.Interpolation.linlin, axes=None,
                 index=None, valueType=enumsModule.ValueType.float64, outerDomainValue=None, label=None):
        """
        :param Ys:                  An instance of :py:class:`valuesModule.Values` representing the y-values of *self.*
        :param interpolation:       The interpolation rule used to determine :math:`y(x)` between to consecutive x-values.
        """

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType,
                                            outerDomainValue=outerDomainValue, label=label)

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if Ys is None:
            Ys = valuesModule.Values([])
        if not isinstance(Ys, valuesModule.Values): raise TypeError('Ys must be an instance of values.values.')
        self.__Ys = Ys
        self.__Ys.setAncestor(self)

    def __len__(self):
        """
        This method returns the number of y-values started in *self*.

        :returns:       A python int.
        """
        # FIXME should this account for self.__Ys.start?

        return len(self.__Ys)

    def __getitem__(self, index):
        """
        This method returned the y-value at index *index*.

        :param index:       The index of the requested y-value.

        :returns:           A python float.
        """
        # FIXME should this account for self.__Ys.start?

        return self.__Ys[index]

    def __add__(self, other):
        """
        This method adds two Ys1d instances, checking to ensure they are on the same x grid.

        :param other:   An instance of :py:class:`Ys1d` to added to *self*.

        :return:        A new :py:class:`Ys1d` instance
        """

        if len(self.__Ys) == 0: return other.copy()

        assert isinstance(other, Ys1d), "Adding Ys1d to %s not supported" % type(other)
        if self.axes is not None and other.axes is not None:
            if self.grid != other.grid:
                raise NotImplementedError("Adding Ys1d with different grids")
        if self.Ys.length != other.Ys.length:
            raise Exception('self.Ys.length = %d != other.Ys.length = %d' % (self.Ys.length, other.Ys.length))

        if self.__Ys.start <= other.Ys.start:
            ys1d_1 = self.copy()
            ys1d_2 = other
        else:
            ys1d_1 = other.copy()
            ys1d_2 = self 

        offset = ys1d_2.Ys.start - ys1d_1.Ys.start
        values = [y for y in ys1d_1.Ys.values]
        for i1, y2 in enumerate(ys1d_2.Ys): values[i1+offset] += y2
        ys1d_1.Ys.values = values

        return ys1d_1

    @property
    def grid(self):
        """
        This method returns the X-grid (following link if necessary). Raises AttributeError if no grid found.

        @return:    An instance of :py:class:`axesModule.Grid`.
        """

        if isinstance(self.axes[1], axesModule.Grid):
            return self.axes[1]
        elif hasattr(self.axes[1], 'link') and isinstance(self.axes[1].link, axesModule.Grid):
            return self.axes[1].link
        else:
            raise AttributeError("No Grid found for Ys1d!")

    @property
    def Ys(self):
        """
        This method returns a reference to the *Ys* member of *self*.

        :returns:       An instance of :py:class:`valuesModule.Values`.
        """

        return self.__Ys

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance.

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps:    Used to dull the lower point for "flat" interpolation.
        :param upperEps:    Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self.toPointwise_withLinearXYs(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps, cls=self.toLinearXYsClass())

        return xys1d

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        factors = self.axes.convertUnits(unitMap)
        yFactor = factors[0]
        if yFactor != 1:
            self.__Ys.values = [yFactor * value for value in self.__Ys]
        self.fixValuePerUnitChange(factors)

    def copy(self):
        """
        This method returns a copy of *self*.

        :returns:       An instance of class of *self*.
        """

        axes = self.axes
        if axes is not None: axes = axes.copy()
        Ys = self.__class__(Ys=self.__Ys.copy(), interpolation=self.interpolation, axes=axes,
                            index=self.index, outerDomainValue=self.outerDomainValue, label=self.label)

        return Ys

    __copy__ = copy

    def evaluate(self, domainValue, extrapolation=enumsModule.Extrapolation.none, epsilon=0):
        """
        This method returns the y-value of *self* evaluated at *domainValue*. This method is currently not implemented and instead
        execute a 'raise NotImplementedError'..

        :param domainValue:     python float.
        :param extrapolation:   TBD.
        :param epsilon:         TBD.

        :returns:               A python float.
        """


        raise NotImplementedError("Still TBD. As a temporary work-around, call toPointwise_withLinearXYs to get XYs1d")

    @property
    def domainMin(self):
        """
        This method returns the minimum domain value for *self*.

        :returns:               A python float.
        """

        return self.grid.domainMin

    @property
    def domainMax(self):
        """
        This method returns the maximum domain value for *self*.

        :returns:               A python float.
        """

        return self.grid.domainMax

    @property
    def domainUnit(self):
        """
        This method returns the domain unit for *self*.

        :returns:       A python str.
        """

        return self.grid.domainUnit

    def domainUnitConversionFactor(self, unitTo):
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A python float.
        """

        return self.grid.domainUnitConversionFactor(unitTo)

    @property
    def domainGrid(self):
        """
        This method returns all domain values for *self* as a python list.

        :returns:           A python list.
        """

        return self.grid.domainGrid

    @property
    def rangeMin(self):
        """
        This method returns the minimum y-value of *self*.

        :returns:           A python float.
        """

        return min(self.__Ys.values)

    @property
    def rangeMax(self):
        """
        This method returns the maximum y-value of *self*.

        :returns:           A python float.
        """

        return max(self.__Ys.values)

    @property
    def rangeUnit(self):
        """
        This method returns the unit for the y-values of *self*. 
 
        :returns:       A python str.
        """


        return self.getAxisUnitSafely(0)

    def rangeUnitConversionFactor(self, unitTo):
        """
         This method returns the factor needed to convert self' y-values to unit *unitTo*.
        
        :param unitTo:      The unit for converting self's y-values.

        :returns:           A python float.
        """

        if unitTo is None: return 1.
        return PQUModule.PQU('1 ' + self.rangeUnit).getValueAs(unitTo)

    def toPointwise_withLinearXYs(self, **kwargs):
        """
        This method returns an :py:class:`XYs1dModule.XYs1d` representation of *self*'.
        The returned class can be changed by the 'cls' key in *kwargs*.

        :param kwargs:      A dictionary.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        cls = kwargs.pop('cls', XYs1dModule.XYs1d)

        xys = [self.domainGrid[self.__Ys.start:self.__Ys.end], self.__Ys]
        return cls(data=xys, dataForm='xsandys', interpolation=self.interpolation, axes=self.axes, index=self.index,
                   valueType=self.valueType, outerDomainValue=self.outerDomainValue, label=self.label)

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        outline = kwargs.get('outline', False)
        if len(self) < 6: outline = False

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if self.interpolation != enumsModule.Interpolation.linlin:
            attributeStr += ' interpolation="%s"' % self.interpolation

        XML_strList = ['%s<%s%s>' % (indent, self.moniker, attributeStr)]
        if self.isPrimaryXData() and self.axes is not None: XML_strList += self.axes.toXML_strList(indent2)
        XML_strList += self.__Ys.toXML_strList(indent2, **kwargs)
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

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath, True) # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % (', '.join(list(extraAttributes.keys()))))

        values = None
        for child in node:
            if child.tag == 'values': values = valuesModule.Values.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        if values is None: raise Exception('No "values" child node found.')

        ys1d = cls(Ys=values, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(ys1d, node, xPath, linkData, **kwargs)
        if len(extraNodes) == 1: values = extraNodes.pop()

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()                             # Per comment above, parseBareNodeCommonAttributes adds to xPath.

        return ys1d

    @staticmethod
    def defaultAxes(labelsUnits=None):
        """
        This static method returns an instance of :py:class:`axesModule.Axes` with two :py:class:`axesModule.Axis` instances.
        The argument *labelsUnits* is of the form::

            { 0: ( 'dependent label',   'dependent unit' ),
              1: ( 'independent label', 'independent unit' ) }

        :param labelsUnits:     A python dict of the x and y labels.

        :return:                An instance of :py:class:`axesModule.Axes`.
        """

        return axesModule.Axes(2, labelsUnits=labelsUnits)
