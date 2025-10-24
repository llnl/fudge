# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
This module contains the class ``XYs1d`` which treats a list of :math:`(x_i, y_i)` pairs as if they were 
the function :math:`y(x)`.  That is, it is a numerical representation of :math:`y(x)`. The function :math:`y(x)` 
is also called a 1-dimensional or univariant function. As an example, let 

    .. math::

        y_1(x) = a + b * x

and

    .. math::

        y_2(x) = c * \exp( - d * x )
            
then, :math:`y_1` and :math:`y_2` can be added, subtracted, multiplied and divided, for example. Similiarly, 
two XYs1d instances can be added, subtracted, multiplied and divided. The only restriction on the :math:`(x_i, y_i)` 
pairs is that :math:`x_i < x_{i+1}`. This is, it most represent a math function; ergo, single valued.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | XYS1d                             | This class represents a GNDS axis node.                               |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions:

    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | Function                                      | Description                                                                                           |
    +===============================================+=======================================================================================================+
    | return_pointwiseXY_AsXYs1d                    | This function returns an instance that is the same class as *self*` from the input arguments.         |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | otherToSelfsUnits                             | This function scales the data and changes the units of an :py:class:`XYs1d` instance                  |
    |                                               | if they are not the same as the units of another :py:class:`XYs1d` instance                           |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | getValueAsUnit                                | This function returns the value of a :py:class:`PQUModule.PQU` instance in unit a requested unit.     |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | getValueAsAxis                                | This function returns the value of a :py:class:`PQUModule.PQU` instance in unit of an axis instance.  |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | getOtherAndUnit                               | This function is used by the multiply and divide methods of the class :py:class:`XYs1d` to ensure     |
    |                                               | that is unit of another instance is the same as *self*.                                               |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
    | allow_XYsWithSameUnits_orNumberWithSameUnit   | This function returns a pointwiseXY instance or a number whose units match that of a                  |
    |                                               | :py:class:`XYs1d` instance                                                                            |
    +-----------------------------------------------+-------------------------------------------------------------------------------------------------------+
"""

"""
Notes:
    1) plot method still using fudge stuff.
    2) need list of values with unit. Used, for example, for method groupOneFunction.
"""

import math

from pqu import PQU as PQUModule

from numericalFunctions import pointwiseXY_C, pointwiseXY

from . import enums as enumsModule
from . import base as baseModule
from . import link as linkModule
from . import axes as axesModule
from . import values as valuesModule

defaultAccuracy = 1e-3
xAxisIndex = 1
yAxisIndex = 0
domainEpsilon = 1e-15

def return_pointwiseXY_AsXYs1d( self, data, units = None, index = None, outerDomainValue = None, axes = None ) :
    """
    This function returns an instance that is the same class as *self*` from the input arguments. For arugments that are None, the member of *self* will be used. 
    For internal use only.

    :param self:                A derived class of :py:class:`XYs1d` whose class is used to construct the returned instance. 
                                Also, its members are used to fill in for arguemnts that are None.
    :param data:                Data points for the returned instance.
    :param units:               Default index for the returned instance.
    :param index:               Default index for the returned instance.
    :param outerDomainValue:    Default outerDomainValue for the returned instance.
    :param axes:                Default axes for the returned instance.
    """

    if( not( isinstance( self, XYs1d ) ) ) : raise TypeError( 'Self must be an instance of XYs1d.' )

    if( index is None ) : index = self.index
    if( axes is None ) : axes = self.axes
    if( outerDomainValue is None ) : outerDomainValue = self.outerDomainValue
    xys1d = self.__class__( data = data, axes = axes, infill = True, safeDivide = False, index = index, outerDomainValue = outerDomainValue )
    if len(xys1d.axes) > 0 and units is not None:
        for k in units : xys1d.axes[k].unit = units[k]

    return( xys1d )

def otherToSelfsUnits(self, other, checkXOnly=False):
    """
    This function scales *other*'s data and changes its units if they are not the same units as in *self*.
    For internal use only.

    :param self:        A :py:class:`XYs1d` whose axes is used to scale, if needed, *other*'s data.
    :param other:       A :py:class:`XYs1d` whose data may be scaled and whose axes units may be changed.
    :param checkXOnly:  If True, only the X-axis data of *other* are scaled to that of *self*.
    """

    if not isinstance(other, XYs1d):
        raise TypeError('other instance not XYs1d instance; is type "%s".' % type(other))

    if len(self.axes) == 0 or len(other.axes) == 0:
        return other.nf_pointwiseXY

    yScale = 1
    xUnitSelf = PQUModule._getPhysicalUnit(self.axes[xAxisIndex].unit)
    xScale = xUnitSelf.conversionFactorTo(other.axes[xAxisIndex].unit)
    if not checkXOnly:
        yUnitSelf = PQUModule._getPhysicalUnit(self.axes[yAxisIndex].unit)
        yScale = yUnitSelf.conversionFactorTo(other.axes[yAxisIndex].unit)

    other = other.nf_pointwiseXY
    if xScale != 1 or yScale != 1:
        other = other.scaleOffsetXAndY(xScale=xScale, yScale=yScale)

    return other

def getValueAsUnit( unit, quantity ) :
    """
    This function returns the value of *quantity* in unit of *unit*.  For internal use only.

    :param unit:        The unit of the returned value.
    :param quantity:    An instance of :py:class:`PQUModule.PQU`.
    """

    if( not( isinstance( quantity, PQUModule.PQU ) ) ) : raise TypeError( 'Quantity is not an instance of PQUModule.PQU' )
    return( quantity.getValueAs( unit ) )

def getValueAsAxis( axis, quantity ) :
    """
    This function returns the value of *quantity* in unit of *axis*.  For internal use only.

    :param axis:        The axis whose unit is used for the returned value.
    :param quantity:    An instance of :py:class:`PQUModule.PQU`.
    """

    return( getValueAsUnit( axis.unit, quantity ) )

def getOtherAndUnit( self, other ) :
    """
    This function is used by the multiply and divide methods of the class :py:class:`XYs1d`. Other must be a :py:class:`XYs1d` instance, 
    a pointwiseXY instance or an object convertible to a :py:class:`PQUModule.PQU` instance. This function returns a tuple of two items.
    The first item returned ia a pointwiseXY instance or number. The second item is the y-axis unit of *other*. For internal use only.

    :param self:    A :py:class:`XYs1d` instance used to set *other*'s units and data if needed.
    :param other:   A :py:class:`XYs1d`, :py:class:`PQUModule.PQU` or pointwiseXY instance whose data (and unit if present) may be alters to agree with the unit of *self*.
    """

    if( isinstance( other, XYs1d ) ) :
        yUnit = other.getAxisUnitSafely( yAxisIndex )
        other = otherToSelfsUnits( self, other, checkXOnly = True )                     # Returned other is a pointwiseXY instance
    elif( isinstance( other, pointwiseXY ) ) :
        yUnit = ''
    else :                  # Must be an object convertible to a PQU instance.
        if( not( isinstance( other, PQUModule.PQU ) ) ) : other = PQUModule.PQU( other )
        yUnit = other.getUnitSymbol( )
        other = other.getValue( )

    return( other, yUnit )

def allow_XYsWithSameUnits_orNumberWithSameUnit( self, other ) :
    """
    This function returns a pointwiseXY instance or a number whose units match that of *self*.
    For internal use only.

    :param self:        A :py:class:`XYs1d` instance 
    :param other:       A :py:class:`XYs1d` or :py:class:`PQUModule.PQU` instance whose data are used in the returned instance.

    :returns:           A :py:class:`PQUModule.PQU` or pointwiseXY instance.
    """

    if( isinstance( other, XYs1d ) ) :
         other = otherToSelfsUnits( self, other )                                       # Returned other is a pointwiseXY instance
    else :
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        if isinstance(other, PQUModule.PQU):
            other = PQUModule.PQU( other, checkOrder = False ).getValueAs( yUnit )
    return( other )

class XYs1d(baseModule.XDataFunctional):
    r"""
    This class presents a 1d function :math:`y(x)` as a list of points :math:`(x_i, y_i)` with :math:`x_i < x_{i+1}`
    and with an interpolation rule for determining :math:`y(x)` between to consecutive x-values. This class has methods
    for adding, subtracting, multipling and dividing a :py:class:`XYs1d` by a number or another :py:class:`XYs1d` instance,
    as well as other operations common to a 1d function.

    In addition to the members in the inherited class :py:class:`baseModule.XDataFunctional`, the following table list the primary members of this class:

    +---------------+-----------------------------------------------------------------------------------+
    | Member        | Description                                                                       |
    +===============+===================================================================================+
    | data          | This is the list of :math:`(x_i, y_i)` points representing the function.          |
    |               | This member is not directly available to the user.                                |
    +---------------+-----------------------------------------------------------------------------------+
    | interpolation | The interpolation rule used to determine :math:`y(x)` between to consecutive      |
    |               | x-values.                                                                         |
    +---------------+-----------------------------------------------------------------------------------+
    """

    moniker = 'XYs1d'
    dimension = 1
    mutableYUnit = True     # For __imul__ and __idiv__.

    def __init__( self, data=None, dataForm="xys", interpolation=enumsModule.Interpolation.linlin, axes=None,
            index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None, 
            initialSize = 10, overflowSize = 10, infill = True, safeDivide = False ) :
        r"""
        :param data:
        :param dataForm:            This specifies how the data are layed out in *data*. Valid values are 'xys', 'xsandys' or 'list'.
        :param interpolation:       The interpolation rule used to determine :math:`y(x)` between to consecutive x-values.
        :param axes:                See documentation for constructor for :py:class:`baseModule.XDataFunctiona`.
        :param index:               See documentation for constructor for :py:class:`baseModule.XDataFunctiona`.
        :param valueType:           See documentation for constructor for :py:class:`baseModule.XDataFunctiona`.
        :param outerDomainValue:    See documentation for constructor for :py:class:`baseModule.XDataFunctiona`.
        :param label:               See documentation for constructor for :py:class:`baseModule.XDataFunctiona`.
        :param initialSize:         This argument sets the initial size of memory for *data* of the member *nf_pointwiseXY*.
        :param overflowSize:        This arguments sets the size of the overflow memory for *data* of the member *nf_pointwiseXY*.
                                    In general, the default value will suffice.
        :param infill:              This arguement/parameter is deprecated.
        :param safeDivide:          This arguement/parameter is deprecated.
        """

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

        interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if data is None: data = []
        if( isinstance( data, XYs1d ) ) : data = data.nf_pointwiseXY

        initialSize = max( initialSize, len( data ) )
        self.__nf_pointwiseXY = pointwiseXY(data = data, dataForm=dataForm, initialSize=initialSize, overflowSize=overflowSize,
            accuracy=defaultAccuracy, biSectionMax=16, interpolation=str(interpolation), infill=infill, safeDivide=safeDivide)

    def __getstate__( self ) :
        """
        This method returns a pickled version ov *self*.

        :returns:       A python dict.
        """

        state = self.__dict__.copy( )

        for key in state :
            item = getattr( self, key )
            if( isinstance( item, pointwiseXY ) ) : break

        nf_pointwiseXY_pickled = state.pop( key )
        state['nf_pointwiseXY_pickled'] = { 'name' : key, 'interpolation' : str(self.interpolation), 'data' : [ xy for xy in self.__nf_pointwiseXY ] }

        return( state )

    def __setstate__( self, state ) :
        """
        This method unpickles *state* into *self*.

        :param state:   The data to unpickle.
        """

        nf_pointwiseXY_pickled = state.pop( 'nf_pointwiseXY_pickled' )
        self.__dict__ = state

        name = nf_pointwiseXY_pickled['name']
        interpolation = nf_pointwiseXY_pickled['interpolation']
        data = nf_pointwiseXY_pickled['data']

        nf_pointwiseXY_pickled = pointwiseXY( data = data, initialSize = len( data ), accuracy = defaultAccuracy,
                biSectionMax = 16, interpolation = interpolation, infill = True, safeDivide = False )

        setattr( self, name, nf_pointwiseXY_pickled )

    def __len__( self ) :
        """
        This method returns the number of points (i.e., x, y pairs) stored within *self*.

        :returns:       A python int.
        """

        return( len( self.nf_pointwiseXY ) )

    def __abs__( self ) :
        """
        This method returns a copy of *self* but with all y-values the absolute of all y-values of *self*.

        :returns:       A new instance that is the same class as *self*.
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.__abs__( ) ) )

    def __neg__( self ) :
        """
        This method returns a copy of *self* but with all y-values the negative of all y-values of *self*.

        :returns:       A new instance that is the same class as *self*.
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.__neg__( ) ) )

    def __add__( self, other ) :
        """
        This method adds *self* and *other*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A new instance that is the same class as *self*.
        """

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        return( self.returnAsClass( self, self.nf_pointwiseXY.__add__( other ) ) )

    __radd__ = __add__

    def __iadd__( self, other ) :
        """
        This method adds *other* to *self*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A reference to *self*.
        """

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        self.nf_pointwiseXY.__iadd__( other )
        return( self )

    def __sub__( self, other ) :
        """
        This method subtracts *other* from *self*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A new instance that is the same class as *self*.
        """

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        return( self.returnAsClass( self, self.nf_pointwiseXY.__sub__( other ) ) )

    def __rsub__( self, other ) :
        """
        This method subtracts *self* from *other*.

        :param other:   A number.

        :returns:       A new instance that is the same class as *self*.
        """

        sub = self.__sub__( other )
        return( sub.__neg__( ) )

    def __isub__( self, other ) :
        """
        This method subtracts *other* from *self*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A reference to *self*.
        """

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        self.nf_pointwiseXY.__isub__( other )
        return( self )

    def __mul__( self, other ) :
        """
        This method multiplies *self* and *other*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '*' )
        points = self.nf_pointwiseXY.__mul__( other )
        return( return_pointwiseXY_AsXYs1d( self, points, units = { yAxisIndex : unit } ) )

    __rmul__ = __mul__

    def __imul__( self, other ) :
        """
        This method multiplies *self* and *other*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A reference to *self*.
        """

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__imul__( other )
        if len(self.axes) > 0 and otherUnit1 != '': self.axes[yAxisIndex].unit = baseModule.processUnits(self.axes[yAxisIndex].unit, otherUnit1, '*')
        return( self )

    # division operators for Python 2.7:
    def __div__( self, other ) :
        """
        This method divides *self* by *other*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '/' )
        points = self.nf_pointwiseXY.__div__( other )
        return( return_pointwiseXY_AsXYs1d( self, points, units = { yAxisIndex : unit } ) )

    def __rdiv__( self, other ) :
        """
        This method divides *other* by *self*.

        :param other:   A number.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        unit2 = self.getAxisUnitSafely( yAxisIndex )
        other, unit1 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        points = self.nf_pointwiseXY.__rdiv__( other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs1d( self, points, units = { yAxisIndex : unit } ) )

    def __idiv__( self, other ) :
        """
        This method divides *self* by *other*.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A reference to *self*.
        """

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__idiv__( other )
        if len(self.axes) > 0 and otherUnit1 != '': self.axes[yAxisIndex].unit = baseModule.processUnits(self.axes[yAxisIndex].unit, otherUnit1, '/')
        return( self )

    # division operators for Python 3.x:
    def __truediv__( self, other ) :
        """
        This method divides *self* by *other*. There is no difference between division and true division for this class.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '/' )
        points = self.nf_pointwiseXY.__truediv__( other )
        return( return_pointwiseXY_AsXYs1d( self, points, units = { yAxisIndex : unit } ) )

    def __rtruediv__( self, other ) :
        """
        This method divides *other* by *self*. There is no difference between division and true division for this class.

        :param other:   A number.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        unit2 = self.getAxisUnitSafely( yAxisIndex )
        other, unit1 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        points = self.nf_pointwiseXY.__rtruediv__( other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs1d( self, points, units = { yAxisIndex : unit } ) )

    def __itruediv__( self, other ) :
        """
        This method divides *self* by *other*. There is no difference between division and true division for this class.

        :param other:   A number or another :py:class:`XYs1d` instance.

        :returns:       A reference to *self*.
        """

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__itruediv__( other )
        if len(self.axes) > 0 and otherUnit1 != '': self.axes[yAxisIndex].unit = baseModule.processUnits(self.axes[yAxisIndex].unit, otherUnit1, '/')
        return( self )

    def __getitem__( self, indexOrSlice ) :
        """
        This method returns a point or slice of the points of *self*.

        :param indexOrSlice:    The index or slice of *self* to return.

        :returns:               Either the pair [x, y] if *indexOrSlice* is an index or a pointwiseXY instance if *indexOrSlice* is a slice.
        """

        if( isinstance( indexOrSlice, slice ) ) :
            start, stop, step = indexOrSlice.indices( len( self ) )
            if( step != 1 ) :  raise ValueError( "For slicing, only a step of 1 (or None) is allowed. Entered step = %d." % step )
            return( self.__getslice__( start, stop ) )
        else :
            return( self.nf_pointwiseXY.__getitem__( indexOrSlice ) )

    def __setitem__( self, indexOrSlice, xy ) :
        """
        This method sets a point or a slice of points of *self* to *xy*.

        :param indexOrSlice:    The index or slice of *self* to set.
        :param xy:              A point to list of points.
        """

        if( isinstance( indexOrSlice, slice ) ) :
            start, stop, step = indexOrSlice.indices( len( self ) )
            if( step != 1 ) :  raise ValueError( "For slicing, only a step of 1 (or None) is allowed. Entered step = %d." % step )
            return( self.__setslice__( start, stop, xy ) )
        else :
            if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
            self.nf_pointwiseXY.__setitem__( indexOrSlice, xy )

    def __getslice__( self, index1, index2 ) :
        """
        This method returns the slice of the points of *self* between *index1* to *index2*.

        :param index1:      The lower index of the slice.
        :param index2:      The upper index of the slice.

        :returns:           A pointwiseXY instance.
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.getslice( index1, index2 ) ) )

    def __setslice__( self, index1, index2, slice1 ) :
        """
        This method sets a slice of points of *self* to *slice1*.

        :param index1:      The lower index where to add the points.
        :param index2:      The upper index where to add the points.
        :param slice1:      A list of points to insert into *self*.
        """

        slice1 = otherToSelfsUnits( self, slice1 )                                  # Returned other is a pointwiseXY instance
        self.nf_pointwiseXY.__setslice__( index1, index2, slice1 )

    @property
    def nf_pointwiseXY( self ) :
        """
        Returns a reference to self's hidden nf_pointwiseXY member.

        :returns:           A pointwiseXY instance.
        """

        return( self.__nf_pointwiseXY )

    @property
    def interpolation( self ) :
        """
        This method returns the interpolation string of *self*.

        :returns:   A python str.
        """

        return enumsModule.Interpolation.checkEnumOrString(self.nf_pointwiseXY.getInterpolation())

    @interpolation.setter
    def interpolation( self, interpolation ) :
        """
        This method sets the interpolation of *self* to *interpolation* without changing the data..

        :param interpolation:   The new interpolation.
        """

        interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)
        self.nf_pointwiseXY.setInterpolation(str(interpolation))

    def applyFunction( self, f, parameters, accuracy = -1, biSectionMax = -1, checkForRoots = False ) :
        """
        This method maps the function *f* onto each y-value of *self*, returning a new :py:class:`XYs1d` instance.
        Additional points may be added to preserve the desired accuracy as determined by *accuracy* and *biSectionMax*.

        For example, the following will transform all negative y-values in an A :py:class:`XYs1d` object to zero:
        >>> newXYs = XYs1d.applyFunction( lambda y,tmp : 0 if y < 0 else y, None )

        Extra parameters needed by *f* are be passed via the *parameters* argument. The function *f* must take a float and *parameters* as arguments,
        and returns a float.

        :param f:               The function which returns a new value for each y-vlaue.
        :param parameters:      A list of parameters pasted to *f*.
        :param accuracy:        The desired accuracy of the return function.
        :param biSectionMax:    This limits the number of bisections.
        :param checkForRoots:   If True, a point is added between two consecutive points if the curve passes through zero between the two points.

        :returns:               A :py:class:`XYs1d` instance.
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.applyFunction( f, parameters, accuracy = accuracy, biSectionMax = biSectionMax, checkForRoots = checkForRoots ) ) )

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance. 
        Beware, this method returns *self* if data have lin-lin interpolation.

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps     Used to dull the lower point for "flat" interpolation.
        :param upperEps     Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self
        if asLinLin and xys1d.interpolation != enumsModule.Interpolation.linlin:
            xys1d = xys1d.changeInterpolation(enumsModule.Interpolation.linlin, accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        return xys1d

    def changeInterpolation( self, interpolation, accuracy, lowerEps = 0., upperEps = 0., cls = None ) :
        """
        This methods returns :py:class:`XYs1d` instance that is *self* converted to a new interpolation. Points are added as needed to
        maintain the accuracy of the returned instance to *self*. Generally, *interpolation* can only be "lin-lin".
        If the interpolation of *self* is "flat", a point of *self* is converted into two points around the point in *self* to ensure that
        the returned instande is not multi-values (which is not supported by the :py:class:`XYs1d` class. The relative distance of the
        new points from the point in *self* is determined by the arguments *lowerEps* and *upperEps*. At least of these arguments must be
        non-zero or an exception is raised.

        :param interpolation:   The interpolation of the returned function.
        :param accuracy:        The desired accuracy of the returned function.
        :param lowerEps:        For flat interpolation, this is the fractional espilon for the position of a point below a point in *self*.
        :param upperEps:        For flat interpolation, this is the fractional espilon for the position of a point above a point in *self*.
        :param cls:             The class of the returned object.

        :returns:               An instance of *cls* or *self* if *cls* is None.
        """

        interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)
        if interpolation != enumsModule.Interpolation.linlin:
            raise ValueError('Only "%s" interpolation currently supported: not %s' % (enumsModule.Interpolation.linlin, interpolation))
        c1 = self.nf_pointwiseXY.changeInterpolation(interpolation=str(interpolation), accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        axes = self.axes
        c1 = return_pointwiseXY_AsXYs1d( self, c1, axes = axes, outerDomainValue = self.outerDomainValue )
        if( cls is None  ) : cls = self

        return( cls.returnAsClass( self, c1, axes = axes, interpolation = interpolation ) )

    def changeInterpolationIfNeeded( self, allowedInterpolations, accuracy, lowerEps = 0, upperEps = 0, cls = None ) :
        """
        This method returns *self* if the interpolation of *self* is the same as one of the interpolations in *allowedInterpolations*.
        Otherwise, it returns the results of :py:func:`changeInterpolation` with the first interpolation in *allowedInterpolations*.

        :param allowedInterpolations:   The interpolation of the returned function.
        :param accuracy:                The desired accuracy of the returned function.
        :param lowerEps:                For flat interpolation, this is the fractional espilon for the position of a point below a point in *self*.
        :param upperEps:                For flat interpolation, this is the fractional espilon for the position of a point above a point in *self*.
        :param cls:                     The class of the returned object.

        :returns:                       Either *self* or the instance returned by :py:func:`changeInterpolation`.
        """

        for interpolation in allowedInterpolations :
            if( interpolation == self.interpolation ) : return( self )
        return( self.changeInterpolation( allowedInterpolations[0], accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = cls ) )

    def cloneToInterpolation( self, interpolation ) :
        """
        This method returns a clone of *self* but with the interpolation of the returned instance set to *interpolation*.
        The data of the returned instance are the same as *self*.

        :returns:       An new instance of class *self*.
        """

        interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)
        c1 = self.nf_pointwiseXY.cloneToInterpolation(str(interpolation))
        return( self.__class__.returnAsClass( self, c1, axes=self.axes, interpolation=interpolation ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if len(self.axes) == 0: return
        factors = self.axes.convertUnits( unitMap )
        if( factors[:2] != [ 1., 1. ] ) : self.nf_pointwiseXY.scaleOffsetXAndY( xScale = factors[1], yScale = factors[0], insitu = True )
        self.fixValuePerUnitChange( factors )

    def clip( self, rangeMin = None, rangeMax = None ) :
        """
        This method returns a instance that is the same class as *self*, with the same data except that the range of the returned instance
        is limited to the range [rangeMin: rangeMax]. Points are added as needed to keep the returned function the same as *self* between
        the range [rangeMin: rangeMax]. For example, rangeMin = 1.5 and rangeMax = 5, and if one point in *self* is [10, 2] and the next two
        points [12, 1] and [14, 3], then the points [11], 1.5] and [12.5, 1.5] are added and the point at x value = 12 is set to [12, 1.5].

        :param rangeMin:    The minimum allowed y-value of the returned instance.
        :param rangeMax:    The maximum allowed y-value of the returned instance.

        :returns:           An new instance of class *self*.
        """

        if( rangeMin is None ) : rangeMin = self.rangeMin
        if( rangeMax is None ) : rangeMax = self.rangeMax
        return( self.returnAsClass( self, self.nf_pointwiseXY.clip( rangeMin, rangeMax ) ) )

    def commonDomainGrid( self, others ) :
        """
        This method returns copies of *self* and *others* that are mapped to the same X-grid. That is, a union is made
        of all the x-values for *self* and *others*, and are mapped to it. *Others* must be
        a list of instances of :py:class:`XYs1d`. All :py:class:`XYs1d` instances must have the same domain.

        :param others:      A python list of :py:class:`XYs1d` instances.

        :returns:           A python list of :py:class:`XYs1d` instances of length len(*others*) + 1.
        """

# BRB
        domainMin = self.domainMin                  # Need to check units. Currently union is raising if not the same.
        grid = self
        for i1, other in enumerate( others ) :
            domainMinO = other.domainMin
            if( domainMin != domainMinO ) :
                raise Exception( "domainMin = %e != other's domainMin = %e for other index = %d" % ( domainMin, domainMinO, i1 ) )
            grid = grid.union( other, fillWithSelf = False )
        return( [ o1.union( grid, fillWithSelf = True ) for o1 in [ self ] + others ] )

    def convertAxisToUnit( self, indexOrName, newUnit ) :
        """
        This method converts the axis at *indexOrName* to unit *newUnit* and returns a :py:class:`XYs1d` instance with
        the data for that axis scaled to the new unit.

        :param indexOrName:     The index or name of the axis to convert to the new unit.
        :param newUnit:         New unit for axis *indexOrName*.

        :returns:               An instance of :py:class:`XYs1d`.
        """

        if len(self.axes) == 0:
            n = self
        else:
            index = self.getAxisIndexByIndexOrName( indexOrName )
            axis = self.axes[index]
            factor = PQUModule.PQU( '1 ' + axis.unit ).getValueAs( newUnit )
            data = []
            xyIndex = 1 - index
            for xy in self :
                xy[xyIndex] *= factor
                data.append( xy )
            n = return_pointwiseXY_AsXYs1d( self, data, units = { index : newUnit } )
        return( self.returnAsClass( self, n, axes = n.axes ) )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       An instance of class of *self*.
        """

        xys = self.nf_pointwiseXY.copy( )
        axes = self.axes.copy()
        theCopy = self.returnAsClass(self, xys, index = self.index, outerDomainValue = self.outerDomainValue, axes = axes)
        theCopy.label = self.label
        theCopy.setAncestor(self.ancestor)

        return theCopy

    __copy__ = copy

    def __deepcopy__( self, memodict = {} ) :
        """
        This method returns a deep copy of *self*.

        :returns:       An instance of class of *self*.
        """

        copy_ = self.copy( )
        memodict[ id(self.axes) ] = copy_.axes
        return copy_

    def copyDataToNestedLists(self):
        """
        This is the same as the method :py:func:`copyDataToXYs`.

        :returns:       A python is list of [float, float].
        """

        return self.copyDataToXYs()

    def areDomainsMutual( self, other ) :
        """
        This method returns True if the domains of *self* and *other* are mutual, and False otherwise.

        :param other:       A :py:class:`XYs1d` instance.

        :returns:           A python boolean.
        """

        return( self.nf_pointwiseXY.areDomainsMutual( other.nf_pointwiseXY ) )

    def copyDataToXYs( self ) :
        """
        This method returns a python list of the points of *self* as [x, y] pairs.

        :returns:       A python is list of [float, float].
        """

        return( self.__nf_pointwiseXY.copyDataToXYs( ) )

    def copyDataToXsAndYs( self ) :
        """
        This method returns two python lists; one for the x-value and one for the y-values.

        :returns:       Two  python is list.
        """

        return( self.__nf_pointwiseXY.copyDataToXsAndYs( ) )

    def domain( self ) :
        """This method is deprecated."""

        return( self.nf_pointwiseXY.domain( ) )

    def dullEdges( self, lowerEps = 0., upperEps = 0., positiveXOnly = 0 ) :
        """
        This method returns a copy of *self* but ensures that the first and last y-values are 0.0, by adding points if needed.

        A point (or two) is added near the first point if; its y-value is not 0.0 and *lowerEps* is not 0.0. If a point is added then::

            -) if *lowerEps* is positive, the first point is set to 0.0 and a point is added at x[0] * ( 1 + lowerEps ) if room available.
            -) if *lowerEps* is negative, a point with 0.0 y-value is added below the first point at x[0] * ( 1 + lowerEps ). Also,
               a point is added at x[0] * ( 1 + abs(lowerEps) ) if room available.

        "If room available" means a point is only added if it fits between the first two original points.
        Similar logic is used for *upperEps* at the last point.

        :param lowerEps:        The fractional espilon for the first point of *self*.
        :param upperEps:        The fractional espilon for the last point of *self*.
        :param positiveXOnly1:  If True, for negative *lowerEps* a point is not added below the first point it would result in a negative x-value.

        :returns:               An instance of class of *self*.
        """

        print('positiveXOnly = ', positiveXOnly)
        dulled = self.nf_pointwiseXY.dullEdges( lowerEps = lowerEps, upperEps = upperEps, positiveXOnly = positiveXOnly )
        return( self.returnAsClass( self, dulled ) )

    def cumulativeBins(self, other):
        """
        This method returns two lists of bins that are determined by two sets of evaluations to perform cumulative points interpolation.

        :param other:   Another :py:class:`XYs1d` instance.

        :returns:       Two python lists. Each list is a list of x-values for *self* and *other*.
        """

        if (self.rangeMin < 0.0 or other.rangeMin < 0.0): raise Exception('Cannot calculate cumulative bins for functions with negative values.')

        selfTrimmed = self.trim()
        otherTrimmed = other.trim()
        if (len(selfTrimmed) == 0 or len(otherTrimmed) == 0): raise Exception('Function does not have any area - 1.')

        selfRunningIntegrals = selfTrimmed.runningIntegral()
        otherRunningIntegrals = otherTrimmed.runningIntegral()
        selfRunningIntegralMax = selfRunningIntegrals[-1]
        otherRunningIntegralMax = otherRunningIntegrals[-1]

        if (selfRunningIntegrals == 0.0 or otherRunningIntegrals == 0.0): raise Exception('Function does not have any area - 2.')

        normalizedSelfRunningIntegrals = [x / selfRunningIntegralMax for x in selfRunningIntegrals]
        normalizedOtherRunningIntegrals = [x / otherRunningIntegralMax for x in otherRunningIntegrals]

        jointedCdfList = sorted(normalizedSelfRunningIntegrals+normalizedOtherRunningIntegrals)
        cleanedCdfList = []
        [cleanedCdfList.append(x) for x in jointedCdfList if x not in cleanedCdfList]

        selfUnnormalizedCdfList = [x * selfRunningIntegralMax for x in cleanedCdfList]
        otherUnnormalizedCdfList = [x * otherRunningIntegralMax for x in cleanedCdfList]

        cumuBins1 = [ selfTrimmed[0][0] ]
        indexOfBin = 1
        numberOfBins = len(selfUnnormalizedCdfList)-1
        integral = selfUnnormalizedCdfList[1]
        for index, runningIntegral in enumerate(selfTrimmed.runningIntegral()):
            if indexOfBin >= numberOfBins: break
            x2, y2 = selfTrimmed[index]
            if index > 0:
                while integral <= runningIntegral:
                    deltaArea = integral - priorRunningIntegral
                    if deltaArea == 0.0:
                        nextX = x2
                    else:
                        if y1 == y2:
                            if y1 == 0.0:
                                nextX = None
                            else:
                                nextX = x1 + deltaArea / y1
                        else:
                            if self.interpolation == enumsModule.Interpolation.flat:
                                nextX = x1 + deltaArea / y1
                            elif self.interpolation == enumsModule.Interpolation.linlin:
                                slope = ( y2 - y1 ) / ( x2 - x1 )
                                sqrtArgument = y1 * y1 + 2.0 * slope * deltaArea
                                if sqrtArgument <= 0:
                                    nextX = x2
                                else:
                                    nextX = x1 + 2.0 * deltaArea / ( y1 + math.sqrt( sqrtArgument ) )
                            else: raise NotImplementedError( 'cumulativeBins is not implemented for '+ str(self.interpolation) )
                    if nextX is not None: cumuBins1.append(nextX)
                    indexOfBin += 1
                    if indexOfBin >= numberOfBins: break
                    integral = selfUnnormalizedCdfList[indexOfBin]
            x1 = x2
            y1 = y2
            priorRunningIntegral = runningIntegral
        cumuBins1.append(selfTrimmed[-1][0])

        cumuBins2 = [ otherTrimmed[0][0] ]
        indexOfBin = 1
        integral = otherUnnormalizedCdfList[1]
        for index, runningIntegral in enumerate(otherTrimmed.runningIntegral()):
            if indexOfBin >= numberOfBins: break
            x2, y2 = otherTrimmed[index]
            if index > 0:
                while integral <= runningIntegral:
                    deltaArea = integral - priorRunningIntegral
                    if deltaArea == 0.0:
                        nextX = x2
                    else:
                        if y1 == y2:
                            if y1 == 0.0:
                                nextX = None
                            else:
                                nextX = x1 + deltaArea / y1
                        else:
                            if other.interpolation == enumsModule.Interpolation.flat:
                                nextX = x1 + deltaArea / y1
                            elif other.interpolation == enumsModule.Interpolation.linlin:
                                slope = ( y2 - y1 ) / ( x2 - x1 )
                                sqrtArgument = y1 * y1 + 2.0 * slope * deltaArea
                                if sqrtArgument <= 0:
                                    nextX = x2
                                else:
                                    nextX = x1 + 2.0 * deltaArea / ( y1 + math.sqrt( sqrtArgument ) )
                            else: raise NotImplementedError( 'cumulativeBins is not implemented for '+ str(other.interpolation) )
                    if nextX is not None: cumuBins2.append(nextX)
                    indexOfBin += 1
                    if indexOfBin >= numberOfBins: break
                    integral = otherUnnormalizedCdfList[indexOfBin]
            x1 = x2
            y1 = y2
            priorRunningIntegral = runningIntegral  
        cumuBins2.append(otherTrimmed[-1][0])

        return cumuBins1, cumuBins2

    def equalProbableBins(self, numberOfBins):
        """
        This method returns the list of x-values that are equal probable bins of self. Currently, only supports lin-lin and histogram interpolation.

        :param numberOfBins:    The number of equal probable bins to generate.

        :returns:       A python list of x-values.
        """

        if self.rangeMin < 0.0: raise Exception('Cannot calculate equal probable bins for a function with negative values.')

        trimmed = self.trim()
        if len(trimmed) == 0: raise Exception('Self does not have any area - 1.')
        runningIntegrals = trimmed.runningIntegral()
        runningIntegralMax = runningIntegrals[-1]
        if runningIntegralMax == 0.0: raise Exception('Self does not have any area - 2.')
        epbs = [ trimmed[0][0] ]
        indexOfBin = 1
        integral = runningIntegralMax / numberOfBins
        absDomainMax = max(abs(trimmed[0][0]), abs(trimmed[-1][0])) # This variable is never used. Should be removed?
        for index, runningIntegral in enumerate(runningIntegrals):
            if indexOfBin >= numberOfBins: break
            x2, y2 = trimmed[index]
            if index > 0:
                while integral <= runningIntegral:
                    deltaArea = integral - priorRunningIntegral
                    if deltaArea == 0.0:
                        nextX = x2
                    else:
                        if y1 == y2:
                            if y1 == 0.0:
                                nextX = None
                            else:
                                nextX = x1 + deltaArea / y1
                        else:
                            if self.interpolation == enumsModule.Interpolation.flat:
                                nextX = x1 + deltaArea / y1
                            elif self.interpolation == enumsModule.Interpolation.linlin:
                                slope = ( y2 - y1 ) / ( x2 - x1 )
                                sqrtArgument = y1 * y1 + 2.0 * slope * deltaArea
                                if sqrtArgument <= 0:
                                    nextX = x2
                                else:
                                    nextX = x1 + 2.0 * deltaArea / ( y1 + math.sqrt( sqrtArgument ) )
                            else: raise NotImplementedError( 'equalProbableBins is not implemented for '+ str(self.interpolation) )
                    if nextX is not None: epbs.append(nextX)
                    indexOfBin += 1
                    if indexOfBin >= numberOfBins: break
                    integral = runningIntegralMax * indexOfBin / numberOfBins
            x1 = x2
            y1 = y2
            priorRunningIntegral = runningIntegral

        epbs.append(trimmed[-1][0])

# FIXME: Should look for case when one of the epbs is finite and very close to 0.0 (i.e., probably should be 0.0). Logic would look something like:
# if epbs[index] * epbs[index-1] < 0.0 # going from a negative value to a positive value.
# Now compare epbs[index] to 1e-12 * epbs[index-1] and to 1e-12 * epbs[index+1]. If smaller, set epbs[index] 0.0.
# Note, in the above epbs[index-1] may be the smaller value (in magnitude) and should be the one testing.

        return epbs

    def evaluate(self, x, **kwargs):
        """
        This method returns the y-value of *self* evaluated at *x*.

        :param x:           A python float.
        :param kwargs:      This argument is not used but kept to be compatible with other evaluate methods.

        :returns:           A python float.
        """

        return( self.__nf_pointwiseXY.evaluate( x ) )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Sets the domain minimum and maximum of *self* per the arguments. If the domain of *self* is within [*domainMin*,
        *domainMax*], *self* is unchanged.

        :param domainMin:       The lower domain value.
        :param domainMax:       The upper domain value.
        :param fixToDomain:     A value of the enum :py:class:`enumsModule.FixDomain`.

        :returns:       0 if domains are not adjusted and 1 otherwise.
        """

        oldDomainMin = self.domainMin
        oldDomainMax = self.domainMax
        oldLengthFlag = min(1, len(self))

        if fixToDomain == enumsModule.FixDomain.lower:
            self.setData(self.domainSlice(domainMin = domainMin, fill = True))
        elif fixToDomain == enumsModule.FixDomain.upper:
            self.setData(self.domainSlice(domainMax = domainMax, fill = True))
        else:
            self.setData(self.domainSlice(domainMin = domainMin, domainMax = domainMax, fill = True))

        if len(self) == 0:
            if oldLengthFlag > 0:
                print('WARNING: %s.fixDomains has removed all data from %s.' % (self.moniker, self.toXLink()))
            return oldLengthFlag
        if oldDomainMin == self.domainMin and oldDomainMax == self.domainMax: return 0
        return 1

    def lowerIndexBoundingX( self, x ) :
        """
        Thie method returns the index of *self* for which *x* greater than or equal to self[index][0]. If *x* is outside the domain of *self*,
        a -1 is returned.

        :param x:   The x-value.

        :returns:   A python int.
        """

        return( self.__nf_pointwiseXY.lowerIndexBoundingX( x ) )

    def multiply(self, other, accuracy, biSectionMax, lowerEps, upperEps, asXYs1d, allowFlat, cls=None):
        """
        The method returns the multiplication of *self* by another xData 1d-function. Currently, this method only works for
        *other* of class :py:class:`XYs1d`. Also, *allowFlat* is currently not implemented.

        :param other:               An xData 1d-function. Currently, must be an instance of :py:class:`XYs1d`.
        :param accuracy:            The accuracy of the multiplication as allowed by *biSectionMax*.
        :param biSectionMax:        The maximum amount of besecting between two initial union points.
        :param lowerEps:            This parameter is passed onto the method :py:func:`changeInterpolationIfNeeded`.
        :param upperEps:            This parameter is passed onto the method :py:func:`changeInterpolationIfNeeded`.
        :param asXYs1d:             If True, the returned instance will be an :py:class:`XYs1d`. Otherwise, it may be an instance of the class of *other*.
        :param allowFlat:           If True, will try to preserve flat interapolation if possible.
        :param cls:                 The class of the returned instance. Default is :py:class:`XYs1d`.
        """

        if cls is None:
            cls = XYs1d

        self2  = self.changeInterpolationIfNeeded([enumsModule.Interpolation.linlin], accuracy, lowerEps=lowerEps, upperEps=upperEps, cls=cls)
        self2.nf_pointwiseXY.setAccuracy(accuracy)
        self2.nf_pointwiseXY.setBiSectionMax(biSectionMax)

        other2 = other.changeInterpolationIfNeeded([enumsModule.Interpolation.linlin], accuracy, lowerEps=lowerEps, upperEps=upperEps, cls=cls)
        other2.nf_pointwiseXY.setAccuracy(accuracy)
        other2.nf_pointwiseXY.setBiSectionMax(biSectionMax)

        product = self2 * other2
        product = cls(product, axes=product.axes)

        return product

    def mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 ) :
        """
        Thie method returns copies of *self* and *other*s which have mutual domains.

        :param lowerEps1:       The fractional espilon for the first point of *self*.
        :param upperEps1:       The fractional espilon for the last point of *self*.
        :param positiveXOnly1:  If True, for negative *lowerEps1* a point is not added below the first point of *self* if it would result in a negative x-value.
        :param other:           An instance of :py:class:XYs1d`.
        :param lowerEps2:       The fractional espilon for the first point of *other*.
        :param upperEps2:       The fractional espilon for the last point of *other*.
        :param positiveXOnly1:  If True, for negative *lowerEps2* a point is not added below the first point of *other* if it would result in a negative x-value.

        :returns:               A copy of *self* and *other*s which have mutualified domains.

        .. note:: Need to check that x units are the same.
        """

        m1, m2 = self.nf_pointwiseXY.mutualify( lowerEps1, upperEps1, positiveXOnly1, other.nf_pointwiseXY, lowerEps2, upperEps2, positiveXOnly2 )
        return( self.returnAsClass( self, m1 ), other.returnAsClass( other, m2 ) )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        Thie method returns a "copy" of *self* whose integral is 1.
        The dimension argument is ignored. Only here to be compatable with calling from XYsnd.normalize.

        :param insitu:          If True, *self* is normalized and returned. Otherwise, a copy of *self* is normalized and returned.
        :param dimension:       This argument is ignored.

        :returns:               A instance that is the same class as *self*. Returns *self* if *insitu* is True.
        """

        xys = self.nf_pointwiseXY.normalize( insitu = insitu )
        if( insitu ) : return( self )
        return( self.returnAsClass( self, xys ) )

    def setData( self, points ) :
        """
        This method replaces the points in *self* with *points*.

        :param points:  A python list of [x, y] pairs.
        """

        self.nf_pointwiseXY.setData( points )

    def setDataFromList(self, pointsAsList):
        """
        This method replaces the points in *self* the list of values in *pointsAsList*. *pointsAsList* must have an even number of points
        of [x1, y1, x2, y2, ... xn, yn] values.

        :param pointsAsList:    A python list of floats.
        """

        self.nf_pointwiseXY.setDataFromList(pointsAsList)

    def setDataFromXsAndYs(self, xs, ys):
        """
        Replaces the points in *self* the list of *xs* and *ys*.
        """

        self.nf_pointwiseXY.setDataFromXsAndYs(xs, ys)

    def setValue( self, xValue, yValue ) :
        """
        """

        self.nf_pointwiseXY.setValue( xValue, yValue )

    def thicken( self, sectionSubdivideMax = 1, dDomainMax = 0., fDomainMax = 1. ) :
        """
        .. note:: Need unit for dDomainMax.
        """

        t = self.nf_pointwiseXY.thicken( sectionSubdivideMax = sectionSubdivideMax, dDomainMax = dDomainMax, fDomainMax = fDomainMax )
        return( self.returnAsClass( self, t ) )

    def thin( self, accuracy ) :
        """
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.thin( accuracy ) ) )

    def thinToNumberOfPoints( self, maximumNumber, numberFraction = 0.95 ) :
        """
        Returns a tuple containing an accuracy and a new instance of self that hopefully has no more than 
        maximumNumber points.  The returned accuracy is the one passed to the thin method that obtained the 
        returned XYs1d instance.  No thinning is performed if the number of points in self is less than 6.
        A scan of accuracies is performed until the number of points is between minimumNumber and maximumNumber.
        minimumNumber is calculated as int( maximumNumber * numberFraction ) with the restriction that
        minimumNumber must be at least 2 less than maximumNumber. numberFraction is limited to be
        between 0.75 and 0.98. If accuracy is negative, self length is not greater than maximumNumber
        and a copy of self is returned as the XYs1d instance.

        Note, there is no guarantee that the returned instance will have less than maximumNumber points.
        The caller should check the returned instance's length.  For example, if self contains y-values 
        that oscillate between 1 and -1, no thinning is possible.
        """

        if( ( len( self ) > maximumNumber ) and ( len( self ) > 5 ) ) :
            numberFraction = min( .98, max( 0.75, numberFraction ) )
            minimumNumber = int( numberFraction * maximumNumber )
            if( minimumNumber > ( maximumNumber - 2 ) ) : minimumNumber = maximumNumber - 2

            accuracyMin = 1e-12
            accuracyMax = 1

            minSelf = self.thin( accuracyMin )
            maxSelf = self.thin( accuracyMax )
            for i1 in range( 10 ) :
                accuracyMid = math.sqrt( accuracyMin * accuracyMax )
                midSelf = self.thin( accuracyMid )
                if( len( midSelf ) < maximumNumber ) :
                    if( len( midSelf ) > minimumNumber ) : break
                    accuracyMax = accuracyMid
                else :
                    accuracyMin = accuracyMid
        else :
            accuracyMid = -1
            midSelf = self.copy( )

        return( accuracyMid, self.returnAsClass( self, midSelf ) )

    def trim( self ) :
        """
        """

        trimmed = self.returnAsClass( self, self.nf_pointwiseXY.trim( ) )
        if (self.interpolation == enumsModule.Interpolation.flat):
            while trimmed[0][1] == 0:
                trimmed = trimmed[1:]
        return trimmed

    def union( self, other, fillWithSelf = 1, trim = 0 ) :
        """
        """

        other = otherToSelfsUnits( self, other, checkXOnly = True )                         # Returned other is a pointwiseXY instance
        t = self.nf_pointwiseXY.union( other, fillWithSelf = fillWithSelf, trim = trim  )
        return( self.returnAsClass( self, t ) )

    @property
    def domainMin( self ) :
        """
        This method returns the minimum domain value for *self*.
        """

        return( self.nf_pointwiseXY.domainMin( ) )

    @property
    def domainMax( self ) :
        """
        This method returns the maximum domain value for *self*.
        """

        return( self.nf_pointwiseXY.domainMax( ) )

    @property
    def domainUnit( self ) :
        """
        This method returns the domain unit for *self*.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( xAxisIndex ) )

    def domainUnitConversionFactor( self, unitTo ) :
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A python float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainGrid( self ) :
        """
        This method returns all domain values for *self* as a python list.

        :returns:           A python list.
        """

        return( self.nf_pointwiseXY.domainGrid( 1.0 ) )

    @property
    def rangeMin( self ) :
        """
        This method returns the minimum y-value of *self*.

        :returns:           A python float.
        """

        return( self.nf_pointwiseXY.rangeMin( ) )

    @property
    def rangeMax( self ) :
        """
        This method returns the maximum y-value of *self*.

        :returns:           A python float.
        """

        return( self.nf_pointwiseXY.rangeMax( ) )

    @property
    def rangeUnit( self ) :
        """
        This method returns the unit for the y-values of *self*.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( yAxisIndex ) )

    def rangeUnitConversionFactor( self, unitTo ) :
        """
        This method returns the factor needed to convert self' y-values to unit *unitTo*.

        :param unitTo:      The unit for converting self's y-values.

        :returns:           A python float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def range(self):
        """Returns the list (rangeMin, rangeMax)."""

        return self.nf_pointwiseXY.range()

    def domainSlice(self, domainMin=None, domainMax=None, fill=1, dullEps=0.0, **kwargs):
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   The lower x-value of the slice, default is domain minimum of self.
        :param domainMax:   The upper x-value of the slice, default is domain maximum of self.
        :param fill:        If True, points are added at domainMin and domainMax if they are not in self, 
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     (Currently not implemented) the lower and upper points are dulled, default is 0.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        if( domainMax is None ) : domainMax = self.domainMax
        s = self.nf_pointwiseXY.domainSlice( domainMin = domainMin, domainMax = domainMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

    def __mod__( self, other ) : raise NotImplementedError( 'Currently, mod is not implemented' )

    def __pow__( self, other ) : raise NotImplementedError( 'Currently, pow is not implemented' )

    def __exp__( self, other ) : raise NotImplementedError( 'Currently, __exp__ is not implemented' )

    def convolute( self, func ) : 
        """
        Method method convolutes *self* with *func*. This function uses pointwiseXY's convolute function to do the grunt work.  

        :param func:    A :py:class:`XYs1d` instance.

        :returns:       An instance that is the same class as *self*.
        """

        if( isinstance( func, XYs1d ) ) : func = func.nf_pointwiseXY
        if( not( isinstance( func, pointwiseXY ) ) ) : raise TypeError( 'func argument of convolute() must be an instance of pointwiseXY' )
        return( self.returnAsClass( self, self.nf_pointwiseXY.convolute( func, 0 ) ) )

    def group( self, xs, f2 = None, f3 = None, norm = None, asXYs = False ) :
        r"""
        This method multi-group *self* with possibly *f2* and/or *f3.
        The argument ``xs`` is a list of x-values with ``xs[i] < xs[i+1]``. This function calculates the integrals

        .. math::

            \int_{xs[i]}^{xs[i+1]} dx \; { f(x) \over n[i] }

        where :math:`f(x)` is defined as

            +------------------------------------+--------------------+
            | :math:`f(x)`                       |  conditional       |
            +====================================+====================+
            | :math:`f_1(x)`                     | for f2 = f3 = None |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_2(x)`           | for f3 = None      |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_3(x)`           | for f2 = None      |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_2(x) \, f_3(x)` | otherwise          |
            +------------------------------------+--------------------+

        :math:`f_1(x)` is self evaluated at `x` and `n[i]` is a normalization determined from the `norm` arugment as

            +---------------+---------------+
            | norm          | n[i]          |
            +===============+===============+
            | None          | 1             |
            +---------------+---------------+
            | 'dx'          | x[i+1] - x[i] |
            +---------------+---------------+
            | a python list | norm[i]       |
            +---------------+---------------+

        The arguments f2 (and f3) must be None or an :py:class:`XYs1d` instance.

        If ``asXYs`` is ``False``, then ``len( xs ) - 1`` integrals are returned.
        If ``asXYs`` is ``True``, the last integral's values is appended to the end to make a list of length ``len( xs )``, and 
        an instance of class ``XYs1d`` is returned with the x-values from ``xs``, the y-values from the integrals and the 
        interpolation is 'flat'.

        Historical note: the word group comes from deterministic neutron transport (e.g., transport used to simulate nuclear reactors).

        :param xs:      The list of multi-group boundaries.
        :param f2:      An instance of :py:class:`XYs1d` or None.
        :param f2:      An instance of :py:class:`XYs1d` or None.
        :param norm:    Can be None, the python str 'dx', or a python list of floats of lenght the number of groups.
        :param asXYs:   A python boolean.
        """

        yUnit = PQUModule.PQU( 1, self.getAxisUnitSafely( yAxisIndex ) )
        if( f2 is None ) :
            if( f3 is None ) : 
                grouped = self.groupOneFunction( xs, norm = norm )
            else :
                grouped = self.groupTwoFunctions( xs, f3, norm = norm )
                yUnit = yUnit * PQUModule.PQU( 1,  f3.getAxisUnitSafely( yAxisIndex ) )
        else :
            yUnit = yUnit * PQUModule.PQU( 1, f2.getAxisUnitSafely( yAxisIndex ) )
            if( f3 is None ) :
                grouped = self.groupTwoFunctions( xs, f2, norm = norm )
            else :
                grouped = self.groupThreeFunctions( xs, f2, f3, norm = norm )
                yUnit = yUnit * PQUModule.PQU( 1, f3.getAxisUnitSafely( yAxisIndex ) )
        if( norm is None ) :
            yUnit = PQUModule.PQU( 1, self.getAxisUnitSafely( xAxisIndex ) ) * yUnit
        elif( norm != 'dx' ) :
            pass                    # Need to add units to norm. That is, norm, grouped and xs should be an instance of Ys.
        if( asXYs ) :
            grouped.append( grouped[-1] )
            axes = axesModule.Axes(labelsUnits = {1: [self.axes[xAxisIndex].label, self.axes[xAxisIndex].unit], 0: ["", yUnit.getUnitSymbol()]})
            grouped = XYs1d([xs, grouped], dataForm='xsandys', interpolation=enumsModule.Interpolation.flat, axes=axes)
        return( grouped )

    def groupOneFunction( self, xs, norm = None ) :
        """
        This function multi-groups *self* using the boundaries *xs*. Since only one function is involved, each group is just the value
        of the integral over its domain divided by the *norm* for that group.
        The argument *norm* specifies how to normalize each group. If *norm* is None, no normalization happens. If
        *norm* is 'dx', each group integral is divided by its domain width. If *norm* is a list of float, there must be a one-to-one
        mapping of the number of float to the number of groups, and the value in the *norm* at the same index as the group is used as
        the norm for that group.  The number of groups is the lenght of *xs* - 1.

        .. note:: Need unit of xs.

        :param xs:      The list of multi-group boundaries.
        :param norm:    Can be None, the python str 'dx', or a python list of floats of lenght the number of groups.

        :returns:       A python list of floats.
        """

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupOneFunction( boundaries, norm = norm ) )

    def groupTwoFunctions( self, xs, f2, norm = None ) :
        """
        This function multi-groups *self*, *f2* and *f3* using the boundaries *xs*. Each group is just the value
        of the integral of the product of *self*, *f2* and *f3* over the group's domain divided by the *norm* for that group.
        The argument *norm* specifies how to normalize each group. If *norm* is None, no normalization happens. If
        *norm* is 'dx', each group integral is divided by its domain width. If *norm* is a list of float, there must be a one-to-one
        mapping of the number of float to the number of groups, and the value in the *norm* at the same index as the group is used as
        the norm for that group.  The number of groups is the lenght of *xs* - 1.

        .. note:: Need unit of xs.

        :param xs:      The list of multi-group boundaries.
        :param f2:      An instance of :py:class:`XYs1d`.
        :param norm:    Can be None, the python str 'dx', or a python list of floats of lenght the number of groups.

        :returns:       A python list of floats.
        """

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( f2, XYs1d ) ) : f2 = f2.nf_pointwiseXY
        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return self.nf_pointwiseXY.groupTwoFunctions( boundaries, f2, norm = norm)

    def groupThreeFunctions( self, xs, f2, f3, norm = None ) :
        """
        This function multi-groups *self*, *f2* and *f3* using the boundaries *xs*. Each group is just the value
        of the integral of the product of *self*, *f2* and *f3* over the group's domain divided by the *norm* for that group.
        The argument *norm* specifies how to normalize each group. If *norm* is None, no normalization happens. If
        *norm* is 'dx', each group integral is divided by its domain width. If *norm* is a list of float, there must be a one-to-one
        mapping of the number of float to the number of groups, and the value in the *norm* at the same index as the group is used as
        the norm for that group.  The number of groups is the lenght of *xs* - 1.

        .. note:: Need unit of xs.

        :param xs:      The list of multi-group boundaries.
        :param f2:      An instance of :py:class:`XYs1d`.
        :param f3:      An instance of :py:class:`XYs1d`.
        :param norm:    Can be None, the python str 'dx', or a python list of floats of lenght the number of groups.

        :returns:       A python list of floats.
        """

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( f2, XYs1d ) ) : f2 = f2.nf_pointwiseXY
        if( isinstance( f3, XYs1d ) ) : f3 = f3.nf_pointwiseXY
        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupThreeFunctions(  boundaries, f2, f3, norm = norm ) )

    def groupFourFunctions( self, xs, f2, f3, f4, norm = None ) :
        """
        This function multi-groups *self*, *f2*, *f3* and *f4* using the boundaries *xs*. Each group is just the value
        of the integral of the product of *self*, *f2*, *f3* and *f4* over the group's domain divided by the *norm* for that group.
        The argument *norm* specifies how to normalize each group. If *norm* is None, no normalization happens. If
        *norm* is 'dx', each group integral is divided by its domain width. If *norm* is a list of float, there must be a one-to-one
        mapping of the number of float to the number of groups, and the value in the *norm* at the same index as the group is used as
        the norm for that group.  The number of groups is the lenght of *xs* - 1.

        .. note:: Need unit of xs.

        :param xs:      The list of multi-group boundaries.
        :param f2:      An instance of :py:class:`XYs1d`.
        :param f3:      An instance of :py:class:`XYs1d`.
        :param f4:      An instance of :py:class:`XYs1d`.
        :param norm:    Can be None, the python str 'dx', or a python list of floats of lenght the number of groups.

        :returns:       A python list of floats.
        """

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( f2, XYs1d ) ) : f2 = f2.nf_pointwiseXY
        if( isinstance( f3, XYs1d ) ) : f3 = f3.nf_pointwiseXY
        if( isinstance( f4, XYs1d ) ) : f4 = f4.nf_pointwiseXY
        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupFourFunctions(  boundaries, f2, f3, f4, norm = norm ) )

    def hasData(self):
        """
        This method returns True if self's len is greater than 1 and False otherwise.

        :returns:   A python boolean.
        """

        return len(self) > 1

    def lowerIndexContainingDomainValue(self, domainValue):
        """
        The method returns the index within *self* for which self[index] <= *domainValue* < self[index+1].
        If *domainValue* is not in the domain of *self*, then -1 is returned.

        :param domainValue: The x-value.

        :returns:           A python3 int.
        """

        return self.nf_pointwiseXY.lowerIndexBoundingX(domainValue)

    def integrate( self, domainMin = None, domainMax = None ) :
        r"""
        This method returns the integral of *self* from *domainMin* to *domainMax*.
        
        .. math::

            \int_{domainMin}^{domainMax} dx \; f(x)

        If *domainMin* (*domainMax*) is None, it is taken from the domain minimum (maximum) of *self*'s domain.

        :returns:       A python float.
        """

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return self.nf_pointwiseXY.integrate( domainMin = domainMin, domainMax = domainMax )

    def integrateWithWeight_x( self, domainMin = None, domainMax = None ) :
        r"""
        This method returns the integral of *self* weighted by :math:`x` from *domainMin* to *domainMax*.
        
        .. math::

            \int_{domainMin}^{domainMax} dx \; x \, f(x)

        If *domainMin* (*domainMax*) is None, it is taken from the domain minimum (maximum) of *self*'s domain.

        :returns:       A python float.
        """

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( self.__nf_pointwiseXY.integrateWithWeight_x( domainMin, domainMax ) )

    def integrateWithWeight_sqrt_x( self, domainMin = None, domainMax = None ) :
        r"""
        This method returns the integral of *self* weighted by :math:`\sqrt{x}` from *domainMin* to *domainMax*.

        .. math::

            \int_{domainMin}^{domainMax} dx \; \sqrt{x} \, f(x)

        If *domainMin* (*domainMax*) is None, it is taken from the domain minimum (maximum) of *self*'s domain.

        :returns:       A python float.
        """

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( self.__nf_pointwiseXY.integrateWithWeight_sqrt_x( domainMin, domainMax ) )

    def indefiniteIntegral( self, domainMin = None, domainMax = None ) :
        r"""
        This method returns an :py:class:`XYs1d` instance with the same x-values as *self* and whose y-values are the running integral of *self*
        to the corresponding x-value.

        The new :py:class:`XYs1d` instance is defined on the range of the old one and the units are wrong.

        :param domainMin:       This limits the lower limit of the integral.
        :param domainMax:       This limits the upper limit of the integral.

        :returns:               An instance of :py:class:`XYs1d`.
        """

# BRB: I think this is just a running sum and may be implemented already. Needs units.
        myAxes = self.axes
        myData = [ [ self.domainMin, 0.0 ] ]
        for i in range( len( self ) - 1 ):
            domainMin = self[i][0]
            domainMax = self[i+1][0]
            myData.append( [ domainMax, myData[-1][1] + self.nf_pointwiseXY.integrate( domainMin = domainMin, domainMax = domainMax ) ] )
        return XYs1d( myData, axes = myAxes )

    def integrateTwoFunctions( self, f2, domainMin = None, domainMax = None ) :
        r"""
        This method returns the integral of *self* times *f2* from *domainMin* to *domainMax*.

        .. math::

            \int_{domainMin}^{domainMax} dx \; f(x) \; f2(x)

        :param f2:          An instance of :py:class:`XYs1d`.
        :param domainMin:   The lower limit of the integral, default is domain minimum of self.
        :param domainMax:   The upper limit of the integral, default is domain maximum of self.

        :returns:           A python float.
        """

        if( not isinstance( f2, XYs1d ) ) : raise TypeError( "f2 must be an instance of an XYs1d" )

        unit = ''
        if len(self.axes) > 0 and len(f2.axes) > 0:
            unit = self.axes[xAxisIndex].unit
            if f2.axes[xAxisIndex].unit != unit:
                f2 = f2.copy( )
                f2.convertAxisToUnit(xAxisIndex, unit)

        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin, f2.domainMin )
        domainMax = min( domainMax, self.domainMax, f2.domainMax )

        return self.groupTwoFunctions([domainMin, domainMax], f2.nf_pointwiseXY)[0]

    def integrateThreeFunctions( self, f2, f3, domainMin = None, domainMax = None ) :
        r"""
        This method returns the integral of *self* times *f2* times *f3* from *domainMin* to *domainMax*.

        .. math::

            \int_{domainMin}^{domainMax} dx \; f(x) \; f2(x) \, f3(x)

        :param f2:          An instance of :py:class:`XYs1d`.
        :param f3:          An instance of :py:class:`XYs1d`.
        :param domainMin:   The lower limit of the integral, default is domain minimum of self.
        :param domainMax:   The upper limit of the integral, default is domain maximum of self.

        :returns:           A python float.
        """

        if( not isinstance( f2, XYs1d ) ) : raise TypeError( "f2 must be an instance of an XYs1d" )
        if( not isinstance( f3, XYs1d ) ) : raise TypeError( "f3 must be an instance of an XYs1d" )

        unit = ''
        if len(self.axes) > 0 and len(f2.axes) > 0 and len(f3.axes) > 0:
            unit = self.axes[xAxisIndex].unit
            if( f2.axes[xAxisIndex].unit != unit ) :
                f2 = f2.copy( )
                f2.convertAxisToUnit( xAxisIndex, unit )
            if( f3.axes[xAxisIndex].unit != unit ) :
                f3 = f3.copy( )
                f3.convertAxisToUnit( xAxisIndex, unit )

        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin, f2.domainMin, f3.domainMin )
        domainMax = min( domainMax, self.domainMax, f2.domainMax, f3.domainMax )

        return self.groupThreeFunctions([domainMin, domainMax], f2, f3 )[0]

    def inverse( self ) :
        """
        This method returns an instance of :py:class:`XYs1d` with pairs (x, y) from *self* as (y, x) in the returned instance.
        This fails if the y-values are not in ascending or descending value.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        return( self.returnAsClass( self, self.nf_pointwiseXY.inverse( ) ) )

    def pdfOfY( self, epsilon, domainMin = None, domainMax = None ) :
        """     
        This method calculates the pdf of the y-values of self (i.e., pdf(y)).

        :param epsilon:
        :param domainMin:   This is used to limit the domain of *self*.
        :param domainMax:   This is used to limit the domain of *self*.

        :returns:       An instance of :py:class:`XYs1d`.
        """

        def addData( pdf, data, epsilon ) :
            """
            """

            data = XYs1d(data=data, axes=axesModule.Axes(2))
            if pdf is None: return data
            pdf, data = pdf.mutualify( epsilon, epsilon, False, data, epsilon, epsilon, False )
            pdf += data
            return pdf
                
        region = self
        
        if( domainMin is None ) : domainMin = region.domainMin
        if( domainMax is None ) : domainMax = region.domainMax
            
        if( len( region ) < 2 ) : raise ValueError( 'self must have at least 2 points' )
        if( domainMin < region.domainMin ) :
            if( region[0][1] != 0.0 ) : raise ValueError( "Specified domainMin is less than self's domainMin." )
            region[0] = [ domainMin, 0.0 ]
        if( domainMax > region.domainMax ) : 
            if( region[-1][1] != 0.0 ) : raise ValueError( "Specified domainMax is greater than self's domainMax." )
            region[len(self)] = [ domainMax, 0.0 ]

        region = region.domainSlice( domainMin, domainMax, fill = True )

        domainWidth = domainMax - domainMin
        deltas = {} 
        pdf = None

        if region.interpolation == enumsModule.Interpolation.flat:
            x1 = None
            for x2, y2 in region :
                if( x1 is not None ) :
                    if( y1 not in deltas ) : deltas[y1] = 0.0
                    deltas[y1] += ( x2 - x1 ) / domainWidth
                x1 = x2
                y1 = y2
        elif region.interpolation == enumsModule.Interpolation.linlin:
            slope = 1
            data = []
            for i1, ( x2, y2 ) in enumerate( region ) :
                if( i1 > 0 ) :
                    if( y1 == y2 ) :
                        if( y1 not in deltas ) : deltas[y1] = 0.0
                        deltas[y1] += ( x2 - x1 ) / domainWidth
                    else :
                        dy = y2 - y1
                        dXdY = ( x2 - x1 ) / ( abs( dy ) * domainWidth )
                        if( i1 == 1 ) :
                            data = [ [ y1, dXdY ] ]
                            slope = dy
                        if( dy * slope < 0 ) :
                            pdf = addData( pdf, data, epsilon )
                            data = [ [ y1, dXdY ] ]
                        if( dy < 0 ) :
                            data.insert( 0, [ y2, dXdY ] )
                        else :
                            data.append( [ y2, dXdY ] )
                        slope = dy
                x1 = x2
                y1 = y2
            pdf = addData( pdf, data, epsilon )
        else :      
            raise Exception( 'Unsupported interpolation = "%s"' % self.interpolation )

        for yValue in sorted( deltas ) :
            if( yValue == 0.0 ) :
                data = [ [ 0.0, 1.0 ], [ epsilon, 0.0 ] ]
            else :
                delta = [ [ yValue * ( 1.0 - epsilon ), 0.0 ], [ yValue, 1.0 ], [ yValue * ( 1.0 + epsilon ), 0.0 ] ]
            delta = XYs1d( data ).normalize( )
            delta *= deltas[yValue]
            pdf = addData( pdf, delta, epsilon )

        pdf.normalize( insitu = True )

        return( pdf )

    def runningIntegral( self ) :
        """
        This method returns a python list of the running integral for each point in *self*.
        The first value is 0.0 as it represents the integral from the first point to itself.

        :returns:       A python list of floats.
        """

        return( self.__nf_pointwiseXY.runningIntegral( ) )

    def scaleDependent( self, value, insitu = False ) :
        """
        This method scales all y-values of self by *value*.

        :param value:           A python float to scale all y-values by.
        :param insitu:          If True, *self* is normalized and returned. Otherwise, a copy of *self* is normalized and returned.

        :returns:               A :py:class:`XYs1d` instance.
        """

        xys = self
        if( not( insitu ) ) : xys = self.copy( )
        xys.nf_pointwiseXY.scaleOffsetXAndY( yScale = value, insitu = True )
        return( xys )

    def scaleOffsetXAndY(self, xOffset=0.0, xScale=1.0, yOffset=0.0, yScale=1.0, insitu=False):
        """
        This method scales and offsets the x- and y-values of *self*. Each value is scaled and then the offset added.

        :param xOffset:         The offset for the x-values.
        :param xScale:          The scale for the x-values.
        :param yOffset:         The offset for the y-values.
        :param yScale:          The scale for the y-values.
        :param insitu:          If True, *self* is normalized and returned. Otherwise, a copy of *self* is normalized and returned.

        :returns:               A :py:class:`XYs1d` instance.
        """

        xys = self
        if not insitu: xys = self.copy()
        xys.nf_pointwiseXY.scaleOffsetXAndY(xOffset=xOffset, xScale=xScale, yOffset=yOffset, yScale=yScale, insitu=True)

        return xys

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) : 
        """
        This method splits *self* at *domainValue* into two XYs1d instances. This method returns either None if domainValue not in
        domain or two :py:class:`XYs1d` instances.

        :param domainValue:     The domain value where to split *self*.
        :param epsilon:         This allows for a relative fuzz at the domain limits.

        :returns:               Two :py:class:`XYs1d` instances.
        """

        domainMin, domainMax = self.domainMin, self.domainMax
        if( domainValue <= ( 1 + epsilon ) * domainMin ) : return( None )
        if( domainValue >= ( 1 - epsilon ) * domainMax ) : return( None )
        return( self.domainSlice( domainMax = domainValue ), self.domainSlice( domainMin = domainValue ) )

    def toLinearXYsClass( self ) :
        """
        This method always returns a reference to the :py:class:`XYs1d` class.

        :returns:           A reference to the :py:class:`XYs1d` class.
        """

        return( XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        """
        This method returns a new :py:class:`XYs1d` with the data of *self* converted to lin-lin interpolation and with 
        points added to maintain desired accuracy.  For more details see :py:func:`toPointwise_withLinearXYs`.

        :param kwargs:      A dictionary of arguments that controls how *self* is converted to a lin-lin instance.

        :returns:           A :py:class:`XYs1d` instance.
        """

        return( self.toPointwise_withLinearXYs( **kwargs ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method returns a new instance, converted to lin-lin interpolation with points added to maintain desired accuracy.

        Optional keys in *kwargs*:

        :param accuracy:    This is the desired accuracy of the returned instance when converting to lin-lin interpolation.
        :param lowerEps:    This is used to dull lower edges for flat interpolation.
        :param upperEps:    This is used to dull upper edges for flat interpolation.
        :param cls:         This is the class of the returned instance.

        :returns:           A :py:class:`XYs1d` instance.
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : defaultAccuracy, 'lowerEps' : 0, 'upperEps' : 0, 'cls' : None } )
        accuracy = arguments['accuracy']
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']
        cls = arguments['cls']
        return self.changeInterpolation(enumsModule.Interpolation.linlin, accuracy, lowerEps=lowerEps, upperEps=upperEps, cls=cls)

    def toString( self, pairsPerLine = 1, format = " %16.8e %16.8e", pairSeparator = "" ) :
        """
        This method returns a string representation *self*.

        :param pairsPerLine:        This is a python int that represents the number of (x, y) pairs per line.
        :param format:              This is the format used to convert an (x, y) pair to a string.
        :param pairSeparator:       This is the separator added after each (x, y) pair except for the last pair.

        :returns:                   A python str.
        """

        return( self.nf_pointwiseXY.toString( pairsPerLine = pairsPerLine, format = format, pairSeparator = pairSeparator ) )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        oneLine = kwargs.get('oneLine', False)
        indent2 = indent + incrementalIndent
        if oneLine: indent2 = ''

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if self.interpolation != enumsModule.Interpolation.linlin: attributeStr += ' interpolation="%s"' % self.interpolation
        startTag = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 

        XML_strList = []
        xys = []
        for xy in self:
            xys += xy
        XML_strList += valuesModule.Values(xys, valueType = self.valueType).toXML_strList(indent2, **kwargs)

        XML_strList = self.buildXML_strList(indent2, startTag, XML_strList, **kwargs)

        if oneLine: XML_strList = [ ''.join(XML_strList) ]

        return XML_strList

    def tweakDomain( self, domainMin = None, domainMax = None, epsilon = domainEpsilon ) :
        """
        This method tweaks the domain limits to match *domainMin* and *domainMax* it they are close.
        This is, this method tweaks the lower (upper) domain of *self* to match *domainMin* (*domainMax*) if the current
        lower (upper) domain value of *self* is within *epsilon*, in a relative sense, of *domainMin* (*domainMax*).
        This method is designed to help treat cases when the domains of two :py:class:`XYs1d` as very close but not the same.

        :param domainMin:       A lower domain value.
        :param domainMax:       An uppoer domain value.
        :param epsilon:         This argument specifies how close the input limits must be to the current domain limits for the 
                                domain limits to be tweak.
        """

        if( len( self ) == 0 ) : return
        if( domainMin is not None ) : 
            x, y = self[0]
            if( abs( domainMin - x ) < epsilon * ( max( abs( domainMin ), abs( x ) ) ) ) : self[0] = domainMin, y
        if( domainMax is not None ) : 
            x, y = self[-1]
            if( abs( domainMax - x ) < epsilon * ( max( abs( domainMax ), abs( x ) ) ) ) : self[-1] = domainMax, y

    def plot( self, logs = 0, domainMin = None, domainMax = None, rangeMin = None , rangeMax = None, title = '',
              multiPlot = False, dataKey = None ) :
        """
        This method plots *self*.
        The x-axis is linear unless *logs* is odd.
        The y-axis is linear unless *logs//2* is odd.

        :param logs:        This argument sets the linearly or log scale of th x and/or y axis.
        :param domainMin:   This set the lower limit of the plotted x-axis.
        :param domainMax:   This set the upper limit of the plotted x-axis.
        :param rangeMin:    This set the lower limit of the plotted y-axis.
        :param rangeMax:    This set the upper limit of the plotted y-axis.
        :param title:       This is the title for the plot.
        :param multiPlot:   If True, not plotting is done, instead plotOptions and plotData are returned.
        :param dataKey:     This is the key for the plot.
        """

        from .interactivePlot import plotbase as plotbaseModule
        from .interactivePlot import multiplot as multiplotModule

        def getUnitlessNumber( value, unit, default ) :
            """
            This method allows the user domain and range limits to have units.
            This function is for internal use.
            """

            if value is None: return default
            if isinstance(value, float): return value
            if not isinstance(value, PQUModule.PQU): value = PQUModule.PQU(value)
            return value.getValueAs(unit)

        if len(self.axes) == 0:
            xUnit, yUnit = '', ''
            xLabel, yLabel = '', ''
        else :
            xUnit = self.getAxisUnitSafely( xAxisIndex )
            yUnit = self.getAxisUnitSafely( yAxisIndex )
            xLabel = self.axes[xAxisIndex].plotLabel( )
            yLabel = self.axes[yAxisIndex].plotLabel( )

        domainMin = getUnitlessNumber( domainMin, xUnit, self.domainMin )
        domainMax = getUnitlessNumber( domainMax, xUnit, self.domainMax )
        rangeMin = getUnitlessNumber( rangeMin, yUnit, self.rangeMin )
        rangeMax = getUnitlessNumber( rangeMax, yUnit, self.rangeMax )

        plotOptions = plotbaseModule.parsePlotOptions( domainMin, domainMax, rangeMin, rangeMax, xLabel, yLabel, title )
        if dataKey is None:
            if self.label:
                dataKey = self.label
            elif not multiPlot:
                dataKey = 'noLabel' if title == '' else title
        plotOptions['logs'] = logs

        plotData = {dataKey: self.copyDataToXsAndYs()}

        if multiPlot:
            return plotOptions, plotData

        else:
            multiplotModule.MultiPlotWithPyQt5(plotOptions, plotData)

    def ysMappedToXs( self, cls, grid, label = None, extendToEnd = False ) :
        """
        This method creates a list of the y-values of *self* evaluated at the x-values specified by *grid*.
        If the values of *grid* extend beyond the domain of *self*, padding of 0.0 is only done if *extendToEnd* is True.

        :param cls:             A derived class of Ys1d.
        :param grid:            An instance of :py:class:`axesModule.Grid`.
        :param label:           The label for the returned instance.
        :param extendToEnd:     If True padding of 0.0's is done if needed.
        """

        offset, Ys = self.nf_pointwiseXY.ysMappedToXs( grid.values )
        if( extendToEnd ) : Ys += ( len( grid.values ) - ( offset + len( Ys ) ) ) * [ 0.0 ]
        ys1d = cls(Ys=valuesModule.Values(Ys, start=offset), axes=self.axes, label=label)
        ys1d.axes[1] = linkModule.Link2(grid.moniker, instance=grid)
        return ys1d

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

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath, True)         # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % (', '.join(list(extraAttributes.keys()))))

        xys1d = cls(**attributes)

        extraNodes = xys1d.parseNodeStandardChildren(node, xPath, linkData, **kwargs)

        if len(extraNodes) == 1:
            values = extraNodes.pop()
            values = valuesModule.Values.parseNodeUsingClass(values, xPath, linkData, **kwargs)
        else:
            values = []
        xys1d.__nf_pointwiseXY.setDataFromList(values)

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()

        return xys1d

    @classmethod
    def returnAsClass( cls, self, other, index = None, outerDomainValue = None, axes = None, interpolation = None ) :
        """
        This class method returns xy data of other as a class of cls. Other must be a sub-class pointwiseXY. cls must be a sub-class 
        of XYs1d.  If index and axes are not specified, they are taken from self. The main use of this classmethod is for 
        methods like __add__ where the addends may be a class derived from XYs1d. For example, the crossSection.XYs1d
        class is derived from The XYs1d class.  If two instances xSec1 and xSec2 of the crossSection.pointwise class 
        are added together (i.e., xSec1 + xSec2) then, since the __add__ method used this classmethod, the returned 
        instance will also be an instance of crossSection.pointwise class.

        :param cls:                 The class of the returned instance.
        :param self:                The source of default meta-data if not specified with arguments.
        :param other:               The xy data.
        :param index:               The index of the returned instance.
        :param outerDomainValue:    The outerDomainValue of the returned instance.
        :param axes:                The axes of the returned instance.
        :param interpolation:       The interpolation for the data.
        """

        if( not( isinstance( self, XYs1d ) ) ) : raise TypeError( 'Self is not a XYs1d instance.' )
        if( isinstance( other, XYs1d ) ) : other = other.nf_pointwiseXY

        if( index is None ) : index = self.index
        if( outerDomainValue is None ) : outerDomainValue = self.outerDomainValue
        if( axes is None ) : axes = self.axes
        if( interpolation is None ) :
            interpolation = self.interpolation
            if( len( self ) == 0 ) : interpolation = other.getInterpolation( )

        xys = cls( data = other, interpolation = interpolation, axes = axes, overflowSize = 10, infill = other.getInfill( ), 
                safeDivide = other.getSafeDivide( ), index = index, outerDomainValue = outerDomainValue )

        return( xys )

    @classmethod
    def createFromFunction( cls, axes, Xs, func, parameters, accuracy, biSectionMax, checkForRoots = False, infill = 1, safeDivide = 1 ) :
        """
        Given an ascending list of x-values (Xs) and a function (func), thie class method creates a linear pointwise 
        representation of the function over the domain of Xs. There must be at least 2 values in Xs.
        See pointwiseXY_C.createFromFunction for all other arguments. The function *func* is called as::

            func(x, parameters)

        :param cls:                 The class of the returned instance.
        :param axes:                The axes of the returned instance.
        :param Xs:                  A list of initial x-values to evalaute *func* at.
        :param func:                The function to convert to an :py:class:`XYs1d` instance.
        :param parameters:          The second argument passed to *func*.
        :param accuracy:            The desired accuracy of the returned function.
        :param biSectionMax:        The maximum number of bisections to preform between each pair in *Xs*.
        :param checkForRoots:       If True, a point is added at zero crossings of *func*.
        :param infill:              TBD.
        :param safeDivide:          TBD.
        """

        xys = pointwiseXY_C.createFromFunction( Xs, func, parameters, accuracy, biSectionMax, checkForRoots = checkForRoots, infill = infill, safeDivide = safeDivide )
        return( cls( data = xys, axes = axes, infill = infill, safeDivide = safeDivide ) )

    @staticmethod
    def defaultAxes( labelsUnits = None ) :
        """
        This static method returns an instance of :py:class:`axesModule.Axes` with two :py:class:`axesModule.Axis` instances.

        :param labelsUnits:     A python dict of the x and y labels.
        """

        return( axesModule.Axes(2, labelsUnits = labelsUnits or {} ) )

    @staticmethod
    def multiPlot( curve1ds, **kwargs ) :
        """
        Thie static method plots a list of 1d curves on the same plot. Uses each curve's 'plotLabel' as its legend key. Each curve 
        must have a copyDataToXsAndYs method.

        :param curve1ds:    A python list of 1d curves to plot.
        :param kwargs:      A python dict of additions parameters (e.g., 'title', 'logs', 'xLabel') used to initialize the plot.
        """

        def argumentValue( name, default ) :
            """
            For internal use only.

            :param name:        The the key for the data in *kwargs*.
            :param default:     The default value if not in *kwargs*.
            """

            kwargsValue = kwargs.get( name, None )
            if( kwargsValue is not None ) : return( kwargsValue )
            return( default )

        from .interactivePlot import plotbase as plotbaseModule
        from .interactivePlot import multiplot as multiplotModule

        if( len( curve1ds ) == 0 ) : return

        plotData = {}
        domainMin = []
        domainMax = []
        rangeMin = []
        rangeMax = []
        for index, curve1d in enumerate( curve1ds ) :
            if( hasattr( curve1d, 'plotLabel' ) ) :
                plotLabel = curve1d.plotLabel
            else :
                plotLabel = 'curve %d' % index
            if not hasattr(curve1d, 'copyDataToXsAndYs'):
                raise TypeError('Curve does not have a "copyDataToXsAndYs" method. Curve type is %s.' % type(curve1d))
            plotData[plotLabel] = curve1d.copyDataToXsAndYs( )
            domainMin.append( curve1d.domainMin )
            domainMax.append( curve1d.domainMax )
            rangeMin.append( curve1d.rangeMin )
            rangeMax.append( curve1d.rangeMax )
        domainMin = argumentValue( 'domainMin', min( domainMin ) )
        domainMax = argumentValue( 'domainMax', max( domainMax ) )
        rangeMin =  argumentValue( 'rangeMin',  min( rangeMin ) )
        rangeMax =  argumentValue( 'rangeMax',  max( rangeMax ) )

        title = kwargs.get( 'title', 'No title' )
        xLabel = kwargs.get( 'xLabel', 'x' )
        yLabel = kwargs.get( 'yLabel', 'y' )

        plotOptions = plotbaseModule.parsePlotOptions( domainMin, domainMax, rangeMin, rangeMax, xLabel, yLabel, title )
        plotOptions['logs'] = kwargs.get( 'logs', 0 )

        multiplotModule.MultiPlotWithPyQt5( plotOptions, plotData )

    def checkEqualProbableBinsResult(self, epbs, epsilon=1e-12, printResults=False):
        """This method checks the results of XYs1d.equalProbableBins and executes a raise if it finds an issue.

        :param self:            The XYs1d instance that equalProbableBins was call for.
        :param epbs:            The results returned by equalProbableBins.
        :param epsilon:         The relative accuracy to compare the values in epbs to.
        :param printResults:    If true, prints information on each value in epbs.

        :return:                The number of errors found.
        """

        errCounter = 0

        trimmed = self.trim()
        runningIntegrals = trimmed.runningIntegral()
        try:
            inverse = XYs1d( data = [ [ x, runningIntegrals[index] ] for index, (x, y) in enumerate(trimmed) ] ).inverse()
        except:
            print( runningIntegrals )
            print( self.toString())
            raise

        numberOfBins = len(epbs) - 1
        for index in range(numberOfBins+1):
            runningIntegral = index * runningIntegrals[-1] / numberOfBins
            integral = self.integrate(domainMax=epbs[index])
            delta = runningIntegral - integral
            ratio = 0.0
            if( delta != 0 ):
                ratio = delta / runningIntegral 
                if abs(ratio) > epsilon:
                    errCounter += 1
            if printResults:
                s = ''
                if abs(ratio) > epsilon: s = ' *********'
                print( '%4d %20.13e %10.10e %18.10e %10.2e, %10.2e' % (index, epbs[index], runningIntegral, integral, delta, ratio ), s )

        return errCounter
