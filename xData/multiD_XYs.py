# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes for representing a 

    -) 2d function :math:`f(x,y)` as a list of :math:`x_i` for :math:`x_i < x_{i+1}` with each :math:`x_i`
        having an associated 1d function :math:`f_i(x)`.
    -) 3d function :math:`g(x,y,z)` as a list of :math:`x_i` for :math:`x_i < x_{i+1}` with each :math:`x_i`
        having an associated 2d function :math:`g_i(y,z)`.

This module contains the following classes:

    +-----------+-------------------------------------------------------------------------------------------------------+
    | Class     | Description                                                                                           |
    +===========+=======================================================================================================+
    | XYsnd     | Base class for classes :py:class:`XYs2d` and :py:class:`XYs3d`.                                       |
    +-----------+-------------------------------------------------------------------------------------------------------+
    | XYs2d     | This class stores a 2d function :math:`f(x,y)` as a list of :math:`x_i` for :math:`x_i < x_{i+1}`     |
    |           | with each :math:`x_i` having an associated 1d function :math:`f_i(x)`.                                |
    +-----------+-------------------------------------------------------------------------------------------------------+
    | XYs3d     | This class stores a 3d function :math:`g(x,y,z)` as a list of :math:`x_i` for :math:`x_i < x_{i+1}`   |
    |           | with each :math:`x_i` having an associated 2d function :math:`g_i(y,z)`.                              |
    +-----------+-------------------------------------------------------------------------------------------------------+


This module contains the following functions:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | flatInterpolationToLinearPoint    | This function returns a domain for converting the flat interpolation  |
    |                                   | function into a smeared pointwise lin-lin interpolation function.     |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

"""
Missing methods

    copyDataToGridWsAndXsAndYs
    def getValue( self, *values ) :
    setFromList
    setFromW_XYs
    thin
    toPointwise_withLinearXYs
    toString
    plot
"""

import abc

from pqu import PQU as PQUModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from . import enums as enumsModule
from . import base as baseModule
from . import axes as axesModule
from . import XYs1d as XYs1dModule
from . import regions as regionsModule
from . import series1d as series1dModule
from . import xs_pdf_cdf as xs_pdf_cdfMoudle

def flatInterpolationToLinearPoint( lowerDomain, upperDomain, domain, epsilon ) :
    """
    This function returns a domain for converting the flat interpolation function into a smeared pointwise lin-lin interpolation function.
    The new point will be in the range (*lowerDomain*, *upperDomain*).  For internal use.

    :param lowerDomain:     The point at or below *domain*.
    :param upperDomain:     The point at of above *domain*.
    :param domain:          The point that will be moved.
    :param epsilon:         The fractional amount to move the new domain value from *domain*.

    :returns:       A python float.
    """

    if( domain < 0 ) :
        domain *= ( 1.0 - epsilon )
    elif( domain > 0 ) :
        domain *= ( 1.0 + epsilon )
    else :
        domain = epsilon

    if( epsilon != 0.0 ) :
        if( ( domain <= lowerDomain ) or ( domain >= upperDomain ) ) : domain = 0.5 * ( lowerDomain + upperDomain )

    return( domain )
            
class XYsnd( baseModule.XDataFunctional ) :
    """
    This class is the base class for classes :py:class:`XYs2d` and :py:class:`XYs3d`. These classes store a function
    as a list of (domain-value, lower dimentional function). For example, let :math:`f_d` be a function of dimension
    :math:`d` and let :math:`x` be the outer most independent axes of :math:`f_d` (i.e., :math:`f_d(x, y, ...)`, then :math:`f_d`
    is effectively stored as a list of (math:`x, f_{d-1}(y, ...)`) where :math:`f_{d-1}(y, ...)` is a function
    of one less dimension then :math:`f_d` and is not a function of :math:`x`. For each math:`f_{d-1}(y, ...)`
    its associated :math:`x` is its outerDomainValue member. The function math:`f_{d-1}(y, ...)` are called
    sub-function of *self*.

    The following table list the primary members of this class:

    +---------------------------+---------------------------------------------------------------+
    | Member                    | Description                                                   |
    +===========================+===============================================================+
    | label                     | This is the label member use by some xData classes.           |
    +---------------------------+---------------------------------------------------------------+
    | axes                      | This is the axes member.                                      |
    +---------------------------+---------------------------------------------------------------+
    | interpolation             | This describes how to interpolation between two consecutive   |
    |                           | domain points.                                                |
    +---------------------------+---------------------------------------------------------------+
    | interpolationQualifier    | This is an additional interpolation rule mainly needed for    |
    |                           | distributions in ensure that interpolation produces physical  |
    |                           | results.                                                      |
    +---------------------------+---------------------------------------------------------------+
    | outerDomainValue          | This is the domain value for the next higher dimension for    |
    |                           | a function that is embedded in a high dimensional functions.  |
    +---------------------------+---------------------------------------------------------------+
    | index                     | This is the index member use by some xData classes.           |
    +---------------------------+---------------------------------------------------------------+
    """

    def __init__(self, interpolation=enumsModule.Interpolation.linlin, axes=None,
            index=None, valueType=enumsModule.ValueType.float64, outerDomainValue=None, label=None, 
            interpolationQualifier=enumsModule.InterpolationQualifier.none):
        """
        :param interpolation:               This is the interpolation rule.
        :param interpolationQualifier:      An additional interpolation rule mainly needed for distributions.
        :param axes:                        This is the axes for the function.
        :param index:                       This is the index member use by some xData classes.
        :param valueType:                   This is currently not used.
        :param outerDomainValue:            This is the domain value for the next higher dimension for a function that is embedded in a high dimensional functions.
        :param label:                       This is the label member.
        """

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

        self.ancestryMembers = baseModule.XDataFunctional.ancestryMembers     # + (self.functionNdsName, ) Need to add this when __functionals is a suite with name __functionNdsName.

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if interpolationQualifier is None:
            interpolationQualifier = enumsModule.InterpolationQualifier.none
        self.interpolationQualifier = enumsModule.InterpolationQualifier.checkEnumOrString(interpolationQualifier)

        self.__functionals = []

    def __len__( self ) :
        """
        This method returns the number of sub-functions *self* contains.

        :returns:       A python int.
        """

        return( len( self.__functionals ) )

    def __getitem__( self, index ) :
        """
        This method returns the number of sub-functions (i.e., function of one less dimension) that *self* contains.

        :param index:   The index of the sub-function to return.

        :returns:       A function of one less dimension than *self*.
        """

        return( self.__functionals[index] )

    def __setitem__( self, index, functional ) :
        """
        This method adds *functional* to *self* at the specified index.

        :param index:   The index of the sub-function to return.
        """

        index1, functional1 = self._set_insertCommon( index, functional.outerDomainValue, functional )
        if( index1 is not None ) :
            if( ( index1 > 0 ) and ( functional1.outerDomainValue <= self.__functionals[index1-1].outerDomainValue ) ) :
                raise ValueError( 'functional.outerDomainValue = %s is <= prior functional.outerDomainValue = %s' % ( functional1.outerDomainValue , self.__functionals[index1-1].outerDomainValue ) )
            if( ( index1 < ( len( self ) - 1 ) ) and ( functional1.outerDomainValue >= self.__functionals[index1+1].outerDomainValue ) ) :
                raise ValueError( 'functional.outerDomainValue = %s is >= next functional.outerDomainValue = %s' % ( functional1.outerDomainValue, self.__functionals[index1+1].outerDomainValue ) )
            self.__functionals[index1] = functional1

    def  __mul__( self, other ) :
        """
        This method multiplies *self* by *other* and returns a new instance with the results.

        :param other:   A python float that *self* is multiplied by.

        :returns:       An instance that is the same class as *self*.
        """

        try :
            value = float( other )
        except :
            raise ValueError( "Other must be a number." )

        functionNd = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, outerDomainValue = self.outerDomainValue,
                label = self.label, interpolationQualifier = self.interpolationQualifier )
        for function in self.__functionals : functionNd.append( function * value )

        return( functionNd )

    __rmul__ = __mul__

    @property
    def domainMin( self ) :
        """
        This method returns the minimum domain value for *self*.

        :returns:       A python float.
        """

        return( self.__functionals[0].outerDomainValue )

    @property
    def domainMax( self ) :
        """
        This method returns the maximum domain value for *self*.

        :returns:       A python float.
        """

        return( self.__functionals[-1].outerDomainValue )

    @property
    def domainUnit( self ) :
        """
        This method returns the domain unit for *self*.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def domainGrid( self ) :
        """
        This method returns all domain values for *self* as a python list.

        :returns:           A python list.
        """

        return( [ functional.outerDomainValue for functional in self ] )

    @property
    def rangeMin( self ) :
        """
        This method returns the minimum range value for *self*.

        :returns:       A python float.
        """

        return( min( [ func.rangeMin for func in self ] ) )

    @property
    def rangeMax( self ) :
        """
        This method returns the maximum range value for *self*.

        :returns:       A python float.
        """

        return( max( [ func.rangeMax for func in self ] ) )

    @property
    def rangeUnit( self ) :
        """
        This method returns the range unit for *self*.

        :returns:       A python str.
        """

        return( self.getAxisUnitSafely( 0 ) )

    @property
    def functionals( self ) :
        """
        This method returns a reference to self's __functionals which is currently a python list.

        :returns:           A python list.
        """

        return( self.__functionals )

    @property
    def functionNdsName( self ) :
        """
        Returns the node name for the child "function#ds" which depends on the dimension of *self*.

        :returns:       A python str.
        """

        return( "function%dds" % ( self.dimension - 1 ) )

    def append( self, functional ) :
        """
        Thie method appends *functional* to the list of self's functions by calling the :py:func:`insert` method.

        :param functional:  The function of one less dimension than *self* to append to *self*.
        """

        self.insert( len( self ), functional )

    def insert( self, index, functional, outerDomainValue = None ) :
        """
        Inserts *functional* at *index*. If *outerDomainValue* is None, outerDomainValue is take from the outerDomainValue of *functional*.

        :param index:               The index where *functional* is placed in *self*'s list of functions.
        :param functional:          The function of one less dimension than *self* to append to *self*.
        :param outerDomainValue:    The outer domain value for *functional*.
        """

        if( outerDomainValue is None ) : outerDomainValue = functional.outerDomainValue
        index1, functional1 = self._set_insertCommon( index, outerDomainValue, functional )
        if( index1 is not None ) :
            if( ( index1 > 0 ) and ( outerDomainValue <= self.__functionals[index1-1].outerDomainValue ) ) :
                raise Exception( 'outerDomainValue = %s is <= prior functionals.outerDomainValue = %s' % ( outerDomainValue, self.__functionals[index1-1].outerDomainValue ) )
            if( outerDomainValue >= self.__functionals[index1].outerDomainValue ) :
                raise Exception( 'outerDomainValue = %s is >= next functionals.outerDomainValue = %s. index = %d' % ( outerDomainValue, self.__functionals[index1].outerDomainValue, index1 ) )
            self.__functionals.insert( index1, functional1 )

    def insertAtValue( self, functional, outerDomainValue = None ) :
        """
        This method inserts *functional* at the appropriate index for *outerDomainValue*. If *outerDomainValue* is not None,
        the inserted functional instance will have outerDomainValue set to *outerDomainValue*, even if functional has an outerDomainValue.

        :param functional:          The function of one less dimension than *self* to append to *self*.
        :param outerDomainValue:    The outerDomainValue for *functional*.
        """

        if( outerDomainValue is None ) : outerDomainValue = functional.outerDomainValue
        outerDomainValue = float( outerDomainValue )
        index = -1               # Set in case self is empty and next line does not set index or functional.
        for index, functional1 in enumerate( self ) :
            if( functional1.outerDomainValue >= outerDomainValue ) : break
        if( index == -1 ) :
            index = 0
        else :
            if( functional1.outerDomainValue == outerDomainValue ) :
                del self.__functionals[index]
            elif( functional1.outerDomainValue < outerDomainValue ) :     # Happens when outerDomainValue is greater than last items outerDomainValue.
                index += 1
        self.insert( index, functional, outerDomainValue = outerDomainValue )

    def pop( self, index ):
        """
        This method removes (pops) the function at *index*.

        :param index:   The index of the function to remove from *self*.

        :returns:       The popped function.
        """

        return self.__functionals.pop( index )

    def _set_insertCommon( self, index, outerDomainValue, functional ) :
        """
        This method checks if it is valid to add *functional* to *self*.
        If (None, None) is returned, then this method added *functional* to self and no further action is needed;
        otherwise, the caller is responsible for add *functional* at the proper location in *self*.
        If the users needs to add *functional*, they must use the index returned by *self*.
        For internal use only.

        :param index:               The index if *self*'s sub-function list with *functional* will reside.
        :param outerDomainValue:    The domain value of *self* where *functional* will reside.
        :param functional:          The sub-function to be added.

        :returns:                   The tuple (None, None) or (index, functional).
        """

        if( not( isinstance( functional, self.allowedSubElements( ) ) ) ) :
            raise TypeError( 'Invalid class "%s" for insertion into "%s".' % ( functional.__class__, self.__class__ ) )
        outerDomainValue = float( outerDomainValue )
        if( not( isinstance( functional, baseModule.XDataFunctional ) ) ) :
            raise TypeError( 'right-hand-side must be instance of XDataFunctional' )
        if( functional.dimension != ( self.dimension - 1 ) ) :
            raise Exception( 'functional dimension = %d not one less than self diemension = %d'
                             % ( functional.dimension, self.dimension ) )
        n1 = len( self )
        if( n1 < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n1 ) )
        index1 = index
        if( index1 < 0 ) : index1 += n1
        if( index1 < 0 ) : raise IndexError( 'index = %s' % index )
        functional.setAncestor( self )
        if( n1 == 0 ) :
            self.__functionals.append( functional )
            return( None, None )
        elif( n1 == index1 ) :
            if( outerDomainValue <= self.__functionals[-1].outerDomainValue ) :
                raise Exception( 'outerDomainValue = %s is <= prior functional.outerDomainValue = %s' % ( outerDomainValue, self.__functionals[-1].outerDomainValue ) )
            self.__functionals.append( functional )
            return( None, None )
        return( ( index1, functional ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        for functional in self : functional.convertUnits( unitMap )

        if len(self.axes) == 0: return
        factors = self.axes.convertUnits( unitMap )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       A new instance of *self*.
        """

        axes = self.axes.copy()
        multid_xys = self.__class__( interpolation = self.interpolation, index = self.index,
                outerDomainValue = self.outerDomainValue, axes = axes, interpolationQualifier = self.interpolationQualifier )
        for i1, functional in enumerate( self ) : multid_xys[i1] = functional.copy( )
        return( multid_xys )

    __copy__ = copy

    def copyDataToNestedLists( self ) :
        """
        This method returns the data in *self* as a nested list of python lists. The depth of the nesting is *self*'s dimension + 1.

        :returns:   A nested list of python lists with the nesting being *self*'s dimension + 1.
        """

        return( [ [ subData.outerDomainValue, subData.copyDataToNestedLists( ) ] for subData in self ] )

    def evaluate(self, domainValue, extrapolation=enumsModule.Extrapolation.none, interpolationQualifier=None, **kwargs):
        """
        Thid method evaluates the *self* at the domain point *domainValue*.
        *Self*'s Interpolation is used if domainValue is between two sub-functions. However, if one of the
        sub-functions is within *domainValue* * epsilon of *domainValue* then that sub-function is returned where
        *epsilon* is a value in *kwargs* with key "epsilon". The default value if 0.0 is used if *kwargs* does not
        contain the key "epsilon".  If two sub-functions are within *domainValue* * epsilon of domainValue, the closest is returned.

        :param domainValue:             The domain value to evaulate *self* at.
        :param extrapolation:           If *domainValue* is outside the domain of *self*, this argument defines how to extrapolate.
        :param interpolationQualifier:  If not None, this interpolation qualifier is used instead of that of *self*.

        :returns:                       A function of one less dimension than *self*.
        """
        epsilon = kwargs.get('epsilon', 0.0)
        if interpolationQualifier is None:
            interpolationQualifier = self.interpolationQualifier
        else:
            interpolationQualifier = enumsModule.InterpolationQualifier.checkEnumOrString(interpolationQualifier)

        outerDomainValue = baseModule.getDomainValue2( domainValue )
        if extrapolation not in enumsModule.Extrapolation:
            raise ValueError( 'Invalid extrapolation outerDomainValue = "%s"' % extrapolation )
        position, function1, function2, frac, interpolation, interpolationQualifier2 = self.getBoundingSubFunctions( domainValue )
        if( position is None ) : raise Exception( "No data to interpolate" )

        fracRel = frac
        if function2 is not None: fracRel = abs(function1.outerDomainValue / function2.outerDomainValue - 1)
        fracRel = min(frac, fracRel)

        if fracRel <= epsilon:                # If close to first point pick it.
            function = function1.copy()
        elif abs(1 - fracRel) <= epsilon:     # If close to second point pick it.
            function = function2.copy()
        else :
            if position in ('=', '<', '>'):
                if position != '=':
                    if extrapolation != enumsModule.Extrapolation.flat:
                        index = { '<' : 0, '>' : -1 }[position]
                        raise Exception( "evaluation point = %s %s than %s" % 
                                (outerDomainValue, { '<' : 'less', '>' : 'greater' }[position], self[index].outerDomainValue))
                function = function1.copy()
            else:
                if isinstance(function1, series1dModule.Series1d):
                    function = ( 1.0 - frac ) * function1 + frac * function2
                    function.outerDomainValue = outerDomainValue
                    return function
                else:
                    if not isinstance(function1, XYs1dModule.XYs1d):      # FIXME, accuracy, lowerEps and upperEps should not be hardwired.
                        if hasattr(function1, 'toPointwiseLinear'):
                            function1 = function1.toPointwiseLinear(accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6)
                        else:
                            function1 = function1.toPointwise_withLinearXYs(accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6)
                    if not isinstance(function2, XYs1dModule.XYs1d):
                        if hasattr(function1, 'toPointwiseLinear'):
                            function2 = function2.toPointwiseLinear(accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6)
                        else:
                            function2 = function2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                    if function1.interpolation != function2.interpolation:
                        raise NotImplementedError("Interpolation between two functions with different interpolation rules")
                    if interpolationQualifier == enumsModule.InterpolationQualifier.unitBase: 
                        if( function1.dimension == 1 ) :
                            if (self.interpolation == enumsModule.Interpolation.linlin):
                                xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1.nf_pointwiseXY,
                                        function2.outerDomainValue, function2.nf_pointwiseXY, 1 )
                            elif (function1.interpolation == enumsModule.Interpolation.flat):
                                xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1.trim().nf_pointwiseXY,
                                        function2.outerDomainValue, function2.trim().nf_pointwiseXY, 1 )
                                function = function1.returnAsClass(function1, xy, outerDomainValue = outerDomainValue, interpolation = function1.interpolation)
                                xylist = regionsModule.Regions1d()
                                xylist.append(function)
                                firstbin = XYs1dModule.XYs1d(data=[(0.,0.),(xylist.copyDataToXsAndYs()[0][0],xylist.copyDataToXsAndYs()[1][0])],dataForm="xys")
                                xylist.prepend(firstbin)
                                xy = xylist.toPointwise_withLinearXYs(lowerEps=1e-15)
                            else :
                                raise Exception( str(self.interpolation)+" is not supported in unit-base interpolation at the moment." )
                        elif( function1.dimension == 2 ) :
                            frac1 = 1.0 - frac
                            frac2 = frac
                            function1 = function1.copy( )
                            function2 = function2.copy( )
                            EPrime1_1 = function1.domainMin
                            EPrime2_1 = function1.domainMax
                            EPrime1_2 = function2.domainMin
                            EPrime2_2 = function2.domainMax
                            EPrime1 = frac1 * EPrime1_1 + frac2 * EPrime1_2
                            EPrime2 = frac1 * EPrime2_1 + frac2 * EPrime2_2

                            energyPrimes = set( )

                            for function1d in function1 :
                                frac = ( EPrime2_1 - function1d.outerDomainValue ) / ( EPrime2_1 - EPrime1_1 )
                                function1d.outerDomainValue = frac * EPrime1 + ( 1.0 - frac ) * EPrime2
                                energyPrimes.add( function1d.outerDomainValue )

                            for function1d in function2 :
                                frac = ( EPrime2_2 - function1d.outerDomainValue ) / ( EPrime2_2 - EPrime1_2 )
                                function1d.outerDomainValue = frac * EPrime1 + ( 1.0 - frac ) * EPrime2
                                energyPrimes.add( function1d.outerDomainValue )

                            energyPrimes = sorted( list( energyPrimes ) )
                            function = function1.__class__( outerDomainValue = outerDomainValue )
                            scale1 = ( EPrime2_1 - EPrime1_1 ) / ( EPrime2 - EPrime1 )
                            scale2 = ( EPrime2_2 - EPrime1_2 ) / ( EPrime2 - EPrime1 )
                            for energyPrime in energyPrimes :
                                function1d1 = function1.evaluate( energyPrime )
                                function1d2 = function2.evaluate( energyPrime )
                                function1d = scale1 * frac1 * function1d1 + scale2 * frac2 * function1d2
                                function.append( function1d )
                            return( function )
                        else :
                            raise ValueError( "Unitbase interpolate of %d dimensional function not supported." % function1.dimension )
                    elif interpolationQualifier == enumsModule.InterpolationQualifier.unitBaseUnscaled:
                        xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1.nf_pointwiseXY,
                                                                                 function2.outerDomainValue, function2.nf_pointwiseXY, 0 )
                    elif interpolationQualifier == enumsModule.InterpolationQualifier.correspondingPoints:
                        if( function1.dimension == 1 ) :
                            nBins = kwargs.get('nBins', 32)
                            xylist = regionsModule.Regions1d()
                            f1eqb = function1.equalProbableBins(nBins)
                            f2eqb = function2.equalProbableBins(nBins)
                            for i in range(nBins):
                                function1partial = function1.copy()
                                function1partial.fixDomains(f1eqb[i], f1eqb[i+1],enumsModule.FixDomain.both)
                                function2partial = function2.copy()
                                function2partial.fixDomains(f2eqb[i], f2eqb[i+1],enumsModule.FixDomain.both)
                                xypartial = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1partial.nf_pointwiseXY,
                                    function2.outerDomainValue, function2partial.nf_pointwiseXY, 1 )
                                function = function1partial.returnAsClass( function1partial, xypartial, outerDomainValue = outerDomainValue, interpolation = function1.interpolation )
                                xylist.append(function)
                            if function1.interpolation == enumsModule.Interpolation.flat:
                                xs, ys = xylist.copyDataToXsAndYs()
                                if xs[0] != 0.0:
                                    firstbin = XYs1dModule.XYs1d(data=[(0.,0.), (xs[0], ys[0])], dataForm="xys")
                                    xylist.prepend(firstbin)
                            xy = xylist.toPointwise_withLinearXYs(lowerEps=1e-15)
                        else :
                            raise ValueError( "CorrespondingPoints interpolate of %d dimensional function not supported." % function1.dimension )
                    elif interpolationQualifier == enumsModule.InterpolationQualifier.cumulativePoints:
                        if( function1.dimension == 1 ) :
                            xylist = regionsModule.Regions1d()
                            cumuBins1, cumuBins2 = function1.cumulativeBins(function2)
                            for i in range(len(cumuBins1)-1):
                                function1partial = function1.copy()
                                function1partial.fixDomains(cumuBins1[i], cumuBins1[i+1], enumsModule.FixDomain.both)
                                function2partial = function2.copy()
                                function2partial.fixDomains(cumuBins2[i], cumuBins2[i+1], enumsModule.FixDomain.both)
                                xypartial = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1partial.nf_pointwiseXY,
                                    function2.outerDomainValue, function2partial.nf_pointwiseXY, 1 )
                                function = function1partial.returnAsClass( function1partial, xypartial, outerDomainValue = outerDomainValue, interpolation = function1.interpolation )
                                xylist.append(function)
                            if (function1.interpolation == enumsModule.Interpolation.flat):
                                firstbin = XYs1dModule.XYs1d(data=[(0.,0.),(xylist.copyDataToXsAndYs()[0][0],xylist.copyDataToXsAndYs()[1][0])],dataForm="xys")
                                xylist.prepend(firstbin)
                            xy = xylist.toPointwise_withLinearXYs(lowerEps=1e-15)
                        else :
                            raise ValueError( "CumulativePoints interpolate of %d dimensional function not supported." % function1.dimension )
                    else :
                        xy = ( 1.0 - frac ) * function1 + frac * function2
                    try :
                        interpolation = xy.interpolation
                    except :
                        interpolation = xy.getInterpolation( )
                    function = function1.returnAsClass( function1, xy, outerDomainValue = outerDomainValue, interpolation = function1.interpolation )

        function.outerDomainValue = outerDomainValue

        return( function )

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :
        """
        This method finds all instances of class *cls* in self's children, grand-children, etc.

        :param cls:     The class to find.
        :param level:   The maximum depth to search in *self*.

        :returns:       A python list of the instances found.
        """

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for functional in self :
            if( isinstance( functional, cls ) ) : foundInstances.append( functional )
            foundInstances += functional.findInstancesOfClassInChildren( cls, level = level )

        return( foundInstances )

    def integrate( self, **limits ):
        """
        This method returns the integral of *self* over the domain defined in *limits*. Supports limits for each axis.

        Example:
        >XYsnd.integrate( energy_in = ('1e-5 eV', '10 eV'), energy_out = ('1 keV', '10 keV') )

        :param limits:  A dictionary containing limits for each independent axis (keyed by axis label or index).
                        If an independent axis is missing from the dictionary, integrate over the entire domain of that axis.

        :return:        A python float.
        """

        domainMin, domainMax = None, None
        if( len( limits ) > 0 ) :
            if self.axes[-1].label in limits :
                domainMin, domainMax = limits.pop( self.axes[-1].label )
            elif self.axes[-1].index in limits :
                domainMin, domainMax = limits.pop( self.axes[-1].index )

        xys_ = []
        for functional in self :
            if isinstance( functional, ( XYs1dModule.XYs1d, series1dModule.Series1d ) ) :
                xys_.append( [ functional.outerDomainValue, functional.integrate( domainMin = domainMin, domainMax = domainMax ) ] )
            elif isinstance( functional, ( XYsnd, regionsModule.Regions ) ) :
                xys_.append( [ functional.outerDomainValue, functional.integrate( **limits ) ] )
            else :
                raise TypeError( "Unsupported class for integration: %s" % type( functional ) )
        xys = [ [ x, float( y ) ] for x, y in xys_ ]

        unit = self.getAxisUnitSafely( self.dimension )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        value = XYs1dModule.XYs1d(xys, interpolation = self.interpolation).integrate(domainMin, domainMax)

        return value

    def interpolateAtValue(self, value, unitBase=False, extrapolation=enumsModule.Extrapolation.none):
        """
        This method returns a functional with dimension one less than self that is the interpolation of self at value. 
        If value is outside the domain of *self* and *extrapolation* is none, a raise is executed. Otherwise,
        a flat interpolated functional is returned.  If *unitBase* is True, then unit base interpolation is performed on 
        the lowest dimension and the dependent data.  This method is deprecated (see evaluate).

        :param value:           The domain value to evaluate *self* at.
        :param unitBase:        A python boolean.
        :param extrapolation:   If *domainValue* is outside the domain of *self*, this argument defines how to extrapolate.

        :returns:                       A function of one less dimension than *self*.
        """

        if extrapolation not in enumsModule.Extrapolation: raise ValueError( 'Invalid extrapolation value = "%s"' % extrapolation )
        if( len( self ) == 0 ) : raise Exception( "No data to interpolate" )
        if( value < self[0].outerDomainValue ) :
            if extrapolation == enumsModule.Extrapolation.flat:
                function = self[0].copy( )
                function.outerDomainValue = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s less than %s" % ( value, self[0].outerDomainValue ) )
        if( value > self[-1].outerDomainValue ) :
            if extrapolation == enumsModule.Extrapolation.flat:
                function = self[-1].copy( )
                function.outerDomainValue = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s greater than %s" % ( value, self[-1].outerDomainValue ) )
        for index, functional2 in enumerate( self ) :
            if( functional2.outerDomainValue >= value ) : break
        if( value == functional2.outerDomainValue ) :
            function = functional2.copy( )
            function.outerDomainValue = value
            return( function )
        functional1 = self[index-1]
# FIXME: following logic only works if functional1 and functional2 are both XYs1d:
        if( unitBase ) :
            xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( value, functional1.outerDomainValue, functional1.nf_pointwiseXY, 
                    functional2.outerDomainValue, functional2.nf_pointwiseXY, 1 )
        else :
            f = ( functional2.outerDomainValue - value ) / ( functional2.outerDomainValue - functional1.outerDomainValue )
            xy = f * functional1 + ( 1. - f ) * functional2
        xyp = functional1.returnAsClass( functional1, xy, outerDomainValue = value )
        return( xyp )

    def getBoundingSubFunctions( self, value ) :
        """
        This method returns the two sub-functions of *self* whose outerDomainValues bound *value*.
        Additonal information is also returned.
        The returned inforamtion is the tuple (flag, functional1, functional2, frac, interpolation, interpolationQualifier).

        Flag is one of
            +-------+---------------------------------------------------------------------------+
            | None  | no data,                                                                  |
            +-------+---------------------------------------------------------------------------+
            | '<'   | value below domainMin,                                                    |
            +-------+---------------------------------------------------------------------------+
            | '>'   | value above domainMax,                                                    |
            +-------+---------------------------------------------------------------------------+
            | '='   | value at functional1 or                                                   |
            +-------+---------------------------------------------------------------------------+
            | ''    | functional1.outerDomainValue <= value < functional2.outerDomainValue.     |
            +-------+---------------------------------------------------------------------------+

        If flag is None then functional1, functional2 and frac are also None.  If flag is not '' then functional2 is None.
        The data frac, interpolation and interpolationQualifier give the fraction position of *value* between the 
        two sub-functions, and how to interpolate between the two returned sub-functions. this method is mainly for internal use.

        :param value:   The domain value which is between the two returned sub-functions outerDomainValues.

        :returns:       A tuple containing information about the bounding sub-functions.
        """

        interpolation = self.interpolation
        interpolationQualifier = self.interpolationQualifier
        if( len( self ) == 0 ) :
            return( None, None, None, None, interpolation, interpolationQualifier )
        elif( len( self ) == 1 ) :
            if( value == self[0].outerDomainValue ) :
                return( '=', self[0], None, 0.0, interpolation, interpolationQualifier )
            symbol = '<'
            if( value > self[0].outerDomainValue ) : symbol = '>'
            return( symbol, self[0], None, 0.0, interpolation, interpolationQualifier )
        elif( value < self[0].outerDomainValue ) :
            frac = ( self[0].outerDomainValue - value ) / max( abs( value ), abs( self[0].outerDomainValue ) )
            return( '<', self[0], None, frac, interpolation, interpolationQualifier )
        elif( value > self[-1].outerDomainValue ) :
            frac = ( value - self[-1].outerDomainValue ) / max( abs( value ), abs( self[-1].outerDomainValue ) )
            return( '>', self[-1], None, frac, interpolation, interpolationQualifier )

        for index, functional2 in enumerate( self ) :
            if( functional2.outerDomainValue >= value ) : break
            functional1 = functional2
        if( value == functional2.outerDomainValue ) : return( '=', functional2, None, 0, interpolation, interpolationQualifier )
        frac = ( value - functional1.outerDomainValue ) / ( functional2.outerDomainValue - functional1.outerDomainValue )
        return( '', functional1, functional2, frac, interpolation, interpolationQualifier )

    def normalize( self, insitu = True, dimension = None ) :
        """
        This method normalizes the data in *self* at dimension. If *dimension* is greater than or equal to *self*'s dimension,
        then *self* is normalized. Otherwise, normalize is called on each sub-fuction of *self*.

        :param insitu:          If True, *self* is normalized and returned.
        :param dimension:       The argument determines at what dimension the data are normalized.

        :returns:               A *self* like instance.
        """

        selfsDimension = self.dimension
        if( dimension is None ) : dimension = selfsDimension
        if( dimension < 0 ) : dimension += selfsDimension
        if( dimension < 1 ) : raise Exception( 'Dimension %d out of range, must be greater than 1' % dimension )
        multid_xys = self
        if( not( insitu ) ) : multid_xys = self.copy( )
        if( dimension == 0 ) : return( multid_xys )
        if( dimension >= selfsDimension ) :
            multid_xys.scaleDependent( 1. / multid_xys.integrate( ), insitu = True )
        else :
            for functional in multid_xys.__functionals : functional.normalize( insitu = True, dimension = dimension )
        return( multid_xys )

    def domainUnitConversionFactor( self, unitTo ) :
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def domainSlice(self, domainMin=None, domainMax=None, fill=1, dullEps = 0., **kwargs):
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   The lower domain value of the slice, default is domain minimum of self,
        :param domainMax:   The upper domain value of the slice, default is domain maximum of self,
        :param fill:        If True, points are added at domainMin and domainMax if they are not in self,
                            otherwise only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     (Currently not implemented) the lower and upper points are dulled, default is 0.

        :returns:           A *self* like instance.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        domainMin = max( domainMin, self.domainMin )
        if( domainMax is None ) : domainMax = self.domainMax
        domainMax = min( domainMax, self.domainMax )

        newMultiD = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, valueType = self.valueType, 
                outerDomainValue = self.outerDomainValue, label = self.label, interpolationQualifier = self.interpolationQualifier )

        domainGrid = [ tmp.outerDomainValue for tmp in self ]
        for idx1, val in enumerate( domainGrid ) :
            if( val >= domainMin ) : break
        for idx2, val in enumerate( domainGrid ) :
            if( val >= domainMax ) : break

        if( domainGrid[idx1] == domainMin ) :
            newMultiD.append( self[idx1].copy( ) )
            idx1 += 1
        else :
            if self.dimension == 3:
                newMultiD.append(self.evaluate(domainMin), **kwargs)
            else :
                newMultiD.append(self.evaluate(domainMin))

        for idx in range( idx1, idx2 ) : newMultiD.append( self[idx].copy( ) )

        if( domainGrid[idx2] == domainMax ) :
            newMultiD.append( self[idx2].copy( ) )
        else: 
            if self.dimension == 3:
                newMultiD.append(self.evaluate(domainMax, **kwargs))
            else:
                newMultiD.append(self.evaluate(domainMax))

        return( newMultiD )

    def fixDomains(self, domainMin, domainMax, fixToDomain, tweakLower=False, **kwargs):
        """
        This method calls :py:func:`domainSlice` to trim the domain of *self* per the arguments.

        :param domainMin:       The lower limit of the domain.
        :param domainMax:       The upper limit of the domain.
        :param fixToDomain:     An instance of :py:class:`enumsModule.FixDomain` that specifies which limits are to be fixed.
        :param tweakLower:      If True and the domain minimum only needs to be moved a little, then the first sub-function's outerDomainValue is set to *domainMin*.
        :param kwargs:          A dictionary of extra arguments that controls how the domain is fixed.

        :returns:               1 if *self*'s domain changed and 0 otherwise.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        if fixToDomain != enumsModule.FixDomain.upper:
            if tweakLower and domainMin < self.__functionals[1].outerDomainValue:
                if abs(domainMin - self.__functionals[0].outerDomainValue) < 0.1 * ( self.__functionals[1].outerDomainValue - self.__functionals[0].outerDomainValue ):
                    self.__functionals[0].outerDomainValue = domainMin
            if domainMin < self.__functionals[0].outerDomainValue:
                minFunctional = self.__functionals[0].copy()
                minFunctional.outerDomainValue = domainMin
                self.insert(0, minFunctional)

        if fixToDomain == enumsModule.FixDomain.lower:
            functionals = self.domainSlice(domainMin=domainMin, fill=True, **kwargs)
        elif fixToDomain == enumsModule.FixDomain.upper:
            functionals = self.domainSlice(domainMax=domainMax, fill=True, **kwargs)
        else:
            functionals = self.domainSlice(domainMin=domainMin, domainMax=domainMax, fill=True, **kwargs)

        self.__functionals = functionals

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

    def rangeUnitConversionFactor( self, unitTo ) :
        """
        This method returns the scale factor needed to convert the *self*'s range to unit *unitTo*.

        :returns:   A python float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def scaleDependent( self, value, insitu = False ) :
        """
        This method scales the range of *self* by *value*. If *insitu* is True, *self* is scaled and returned; otherwise, a 
        copy of *self* is scaled and returned.

        :param value:       The scale fractor.
        :param insitu:      If True, *self* is scales and returned.

        :returns:           A *self* like instance.
        """

        multid_xys = self
        if( not( insitu ) ) : multid_xys = self.copy( )
        for functional in multid_xys : functional.scaleDependent( value, insitu = True )

    def toLinearXYsClass(self):
        """
        This method returns self.__class__ since that is the linear class for this dimensional data.

        :returns:   Always returns the class for *self*.
        """

        return self.__class__

    def toPointwiseLinear( self, **kwargs ) :
        """
        This method returns an instance of *self*'s that is lin-lin interpolable. This method only works in *self*'s
        interpolation is "lin-lin" or "flat" and executes a raise if not..

        :param kwargs:      A dictionary of extra arguments that controls how linearize is done.

        :returns:           A *self* like instance.
        """

        if self.interpolation not in [enumsModule.Interpolation.flat, enumsModule.Interpolation.linlin]:
            raise TypeError( 'Unsupported interpolation = "%s".' % self.interpolation )

        flatInterpolation = self.interpolation == enumsModule.Interpolation.flat
        if( flatInterpolation ) :
            lowerEps = kwargs.get( 'lowerEps', 0.0 )
            upperEps = kwargs.get( 'upperEps', 0.0 )
            lowerEps = abs( lowerEps )
            upperEps = abs( upperEps )
            if( ( lowerEps <= 0.0 ) and ( upperEps <= 0.0 ) ) :
                raise ValueError( 'For "%s" interpolation, one or both of lowerEps and/or upperEps must be greater than zero.' % self.interpolation )
        lowerFlatInterpolation = flatInterpolation and ( lowerEps > 0 )
        upperFlatInterpolation = flatInterpolation and ( upperEps > 0 )

        pointwiseLinear = self.toLinearXYsClass( )
        pointwiseLinear = pointwiseLinear(interpolation=enumsModule.Interpolation.linlin, axes=self.axes, index=self.index, 
                valueType = self.valueType, outerDomainValue = self.outerDomainValue, interpolationQualifier = self.interpolationQualifier )

        lastIndex = len( self.__functionals ) - 1
        for index, subFunction in enumerate( self ) :
            if( hasattr( subFunction, 'toPointwiseLinear' ) ) :
                subFunctionPointwiseLinear = subFunction.toPointwiseLinear( **kwargs )
            else :
                subFunctionPointwiseLinear = subFunction.toPointwiseLinear( **kwargs )

            if( flatInterpolation and ( index > 0 ) ) :
                lowerDomainValue = priorSumPointwiseLinear.outerDomainValue
                upperDomainValue = subFunction.outerDomainValue
                if( index == lastIndex ) : lowerEps = 0.0

                priorSumPointwiseLinear.outerDomainValue = flatInterpolationToLinearPoint( lowerDomainValue, upperDomainValue, upperDomainValue, -lowerEps )
                pointwiseLinear.append( priorSumPointwiseLinear )

                if( index < lastIndex ) :
                    nextDomainValue = self[index+1].outerDomainValue
                    subFunctionPointwiseLinear.outerDomainValue = flatInterpolationToLinearPoint( upperDomainValue, nextDomainValue, upperDomainValue, upperEps )
                    pointwiseLinear.append( subFunctionPointwiseLinear )
                    priorSumPointwiseLinear = subFunctionPointwiseLinear.copy( )
            else :
                pointwiseLinear.append( subFunctionPointwiseLinear )
                if( flatInterpolation ) : priorSumPointwiseLinear = subFunctionPointwiseLinear.copy( )

        return( pointwiseLinear )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This function returns a instance of *self*'s class with the lowest independent axis converted to an XYs1d instance.

        :param kwargs:      A dictionary of extra arguments that controls how linearize is done.

        :returns:           A *self* like instance.
        """

        arguments = self.getArguments( kwargs, { 'cls' : None } )
        cls = arguments['cls']
        kwargs.pop( 'cls', None )

        if( cls is None ) : cls = self.__class__
        newMultiD = cls( interpolation = self.interpolation, axes = self.axes, index = self.index, valueType = self.valueType, 
                outerDomainValue = self.outerDomainValue, label = self.label, interpolationQualifier = self.interpolationQualifier )
        for subsec in self:
            newPW = subsec.toPointwise_withLinearXYs( cls = subsec.toLinearXYsClass( ), **kwargs )
            newPW.outerDomainValue = subsec.outerDomainValue
            newMultiD.append( newPW )

        return newMultiD

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)

        indent2 = indent  + kwargs.get('incrementalIndent', '  ')
        indent3 = indent2 + kwargs.get('incrementalIndent', '  ')
        if formatVersion == GNDS_formatVersionModule.version_1_10: indent3 = indent2

        outline = kwargs.get('outline', False)
        if len(self) < 6: outline = False

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if self.interpolation != enumsModule.Interpolation.linlin:
            attributeStr += ' interpolation="%s"' % self.interpolation
        if self.interpolationQualifier != enumsModule.InterpolationQualifier.none:
            attributeStr += ' interpolationQualifier="%s"' % self.interpolationQualifier

        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 
        if self.isPrimaryXData() and self.axes is not None: XML_strList += self.axes.toXML_strList(indent2)
        if 'oneLine' not in kwargs:
            if self.dimension == 2: kwargs['oneLine'] = True

        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList.append('%s<%s>' % (indent2, self.functionNdsName))
        if outline:
            XML_strList += self.__functionals[0].toXML_strList(indent3, **kwargs)
            XML_strList += self.__functionals[1].toXML_strList(indent3, **kwargs)
            XML_strList += [ '%s    ... ' % indent3 ]
            XML_strList += self.__functionals[-2].toXML_strList(indent3, **kwargs)
            XML_strList += self.__functionals[-1].toXML_strList(indent3, **kwargs)
        else:
            for functional in self.__functionals: XML_strList += functional.toXML_strList(indent3, **kwargs)
        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList[-1] += "</%s>" % (self.functionNdsName)
        if self.uncertainty: XML_strList += self.uncertainty.toXML_strList(indent2, **kwargs)
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
        attributes['interpolationQualifier'] = extraAttributes.pop('interpolationQualifier', enumsModule.InterpolationQualifier.none)
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        multid_xys = cls(**attributes)
        
        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(multid_xys, node, xPath, linkData, **kwargs)

        if len(extraNodes) > 0:                                                     # The next few line support GNDS 1.10 and 2.0.
            for index, extraNode in enumerate(extraNodes):
                if extraNode.tag == multid_xys.functionNdsName: break
            if extraNode.tag == multid_xys.functionNdsName: extraNodes = extraNodes.pop(index)

        if 'axes' not in kwargs: kwargs['axes'] = multid_xys.axes
        allowedSubElements = cls.allowedSubElements()
        for child in extraNodes:
            childClass = None
            for allowedChildClass in allowedSubElements:
                if allowedChildClass.moniker == child.tag:
                    childClass = allowedChildClass
                    break
            if childClass is None: raise TypeError('unknown sub-element "%s" in element "%s"' % ( child.tag, cls.moniker ))
            xdata = childClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            multid_xys.append(xdata)

        xPath.pop()                             # Per comment above, parseBareNodeCommonAttributes adds to xPath.

        return( multid_xys )

    @classmethod
    def defaultAxes( cls, labelsUnits ) :
        """
        This class method returns an :py:class:`axesModule.Axes` suitable for class *cls* whose labels and units can be
        set via *labelsUnits*. See :py:class:`axesModule.Axes` for how to use *labelsUnits*.


        :param labelsUnits:     A dictionary that can specify labels and units for some of the axis.

        :return: new axes instance
        """

        return( axesModule.Axes(cls.dimension + 1, labelsUnits = labelsUnits ) )

class XYs2d( XYsnd ) :
    """
    This classes represents a 2d function :math:`f(x,y)` effective as a list of :math:`x_i, f_i(y)` objects where :math:`x_i < x_{i+1}` and
    :math:`f_i(y)` is a 1d function for :math:`f(x,y)` at :math:`x_i`. The :math:`x_i` is stored as the outerDomainValue of :math:`f_i(y)`.
    The interpolation and interpolationQualifier define has to obtain f:math:`f(x,y)` between :math:`x_i` and :math:`x_{i+1}`.
    """

    moniker = 'XYs2d'
    dimension = 2

    @staticmethod
    def allowedSubElements():

        return ( XYs1dModule.XYs1d, series1dModule.LegendreSeries, series1dModule.Polynomial1d, series1dModule.LinearSpline1d, 
                regionsModule.Regions1d, xs_pdf_cdfMoudle.Xs_pdf_cdf1d )

class XYs3d( XYsnd ) :
    """
    This classes represents a 3function :math:`f(x,y,z)` effective as a list of :math:`x_i, f_i(y,z)` objects where :math:`x_i < x_{i+1}` and
    :math:`f_i(y,z)` is a 2d function for :math:`f(x,y,z)` at :math:`x_i`. The :math:`x_i` is stored as the outerDomainValue of :math:`f_i(y,z)`.
    The interpolation and interpolationQualifier define has to obtain f:math:`f(x,y,z)` between :math:`x_i` and :math:`x_{i+1}`.
    """

    moniker = 'XYs3d'
    dimension = 3

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, regionsModule.Regions2d ) )
