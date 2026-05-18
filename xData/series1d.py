# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
This module contains 1d functions that store a function as a polynomial expansion.

This module contains the following classes:

    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | Series1d          | This is the base class for the other polynomial classes.                                              |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | LegendreSeries    | This class stores the function as a Legendre polynomials.                                             |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Polynomial1d      | This class stores the function as a simple polynomial expansion.                                      |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | LinearSpline1d    | This class is a simple linear spline.                                                                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+

This module contains the following functions:
    
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | Legendre          | This function evaluates :math:`P_l(\mu)`.                                                             |
    +-------------------+-------------------------------------------------------------------------------------------------------+
"""

import abc
import math

from pqu import PQU as PQUModule

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from numericalFunctions import Legendre as Legendre_C
maxLegendreOrder = Legendre_C.maxMaxOrder( )

from . import enums as enumsModule
from . import base as baseModule
from . import values as valuesModule
from . import axes as axesModule
from . import XYs1d as XYs1dModule

def Legendre( n, mu, checkXRange = True ) :
    r"""
    Thid function returns the value of the Legendre polynomial of order *n* at *mu* using the recursive relationship.

    :param n:               The order of the Legendre polynomial.
    :param mu:              The :math:`\mu` value.
    :param checkXRange:     If True, :math:`\mu` is check to ensure that it is in the domain [-1, 1].

    :returns:               A python float.
    """

    if( n < 0 ) : raise ValueError( "\nError, n = %d < 0" % n )
    if( checkXRange and ( abs( mu ) > 1 ) ) : raise ValueError( "Legendre: |mu| > 1; mu = %g" % mu )
    Pn = 0.
    Pnp1 = 1.
    n_ = 0
    twoNp1 = 1
    while( n_ < n ) :
        Pnm1 = Pn
        Pn = Pnp1
        n_p1 = n_ + 1
        Pnp1 = ( twoNp1 * mu * Pn - n_ * Pnm1 ) / n_p1
        twoNp1 += 2
        n_ = n_p1
    return( Pnp1 )


class Series1d( baseModule.XDataFunctional ) :
    r"""
    This class is the base class for storing a 1d function as a polynomial series. The function store the 
    polynomial coefficients and has methods for manipulations of the coefficients that are generic to
    all polynomial series (e.g., simple polynomial, Legendre).

    Mathematically and from the table below, the function can be written as

    :math:`\sum_{i=l}^{l+n} \, {\rm Cs[i-l]} \, x^i`.

    where Cs are the coefficients, l is *lowerIndex*, n is the number of coefficients and :math:`x` is the domain value
    for :math:`{\rm domainMin} \le x \le {\rm domainMax}`.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | coefficients      | This is a python list of the regions that have been added.    |
    +-------------------+---------------------------------------------------------------+
    | domainMin         | This is a python list of the regions that have been added.    |
    +-------------------+---------------------------------------------------------------+
    | domainMax         | This is a python list of the regions that have been added.    |
    +-------------------+---------------------------------------------------------------+
    | lowerIndex        | This is the power for the first coefficient.                  |
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
    | valueType         | This describes the type of data in array (i.e., float, int).  |
    +-------------------+---------------------------------------------------------------+
    """

    dimension = 1
    ancestryMembers = baseModule.XDataFunctional.ancestryMembers # + ( 'coefficients', )

    def __init__(self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):
        """
        :param coefficients:        This is the list of coefficients for the 1d polymonial function.
        :param domainMin:           This is the minimum domain value that the function is valid for.
        :param domainMax:           This is the maximum domain value that the function is valid for.
        :param lowerIndex:          This is the power for the first coefficient.
        :param axes:                This is the axes for the function.
        :param index:               This is the index member use by some xData classes.
        :param valueType:           This is currently not used.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is embedded in a high dimensional function.
        :param label:               This is the label member.
        """

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

        self.__domainMin = float( domainMin )
        self.__domainMax = float( domainMax )
        self.__lowerIndex = int( lowerIndex )

        self.coefficients = coefficients

    def __len__( self ) :
        """
        This method returns the number of coefficients in the instance (e.g., for Legendre series it is lMax + 1).

        :returns:       A python int.
        """

        return( len( self.__coefficients ) )

    def __getitem__( self, l ) :
        """
        This method returns the (l+1)^th coefficient.

        :param l:       The index for the coefficient to return.

        :returns:       A python float.
        """

        return( self.__coefficients[l] )

    def __setitem__( self, l, c_l ) :
        """
        This method sets the (l+1)^th coefficient to c_l. l must be less than or equal to lMax.

        :param l:       The index for the coefficient to set.
        :param c_l:     The new value for the coefficient.
        """

        if( l == len( self ) ) :
            self.__coefficients.append( float( c_l ) )
        else :
            self.__coefficients[l] = float( c_l ) 

    def __add__( self, other ) :
        """
        This method returns a series that is the sum of *self* and *other*. *Other* must be a number or the same type as *self*.

        :param other:   A number or an instace that is the same class as *self*.

        :returns:       A new instance that is the same class as *self*.
        """

        try :
            other = float( other )
            c_ls = self.copy( )
            for l, c_l in enumerate( c_ls ) : c_ls.coefficients[l] += other
        except :
            self.checkSameSeriesType( other, 'add' )
            c_l1, c_l2 = self, other
            if( len( self ) < len( other ) ) : c_l1, c_l2 = other, self
            c_ls = c_l1.copy( )
            for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __sub__( self, other ) :
        """
        This method returns a series that is the difference of *self* and *other*. *Other* must be a number or the same type as *self*.

        :param other:   A number or an instace that is the same class as *self*.

        :returns:       A new instance that is the same class as *self*.
        """

        try :
            other = float( other )
            c_ls = self.copy( )
            for l, c_l in enumerate( c_ls ) : c_ls.coefficients[l] -= other
        except :
            self.checkSameSeriesType( other, 'subtract' )
            c_l1, c_l2 = self.__coefficients, other.coefficients
            if( len( self ) < len( other ) ) : c_l1, c_l2 = c_l2, c_l1
            c_ls = c_l1.copy( )
            for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __mul__( self, value ) :
        """
        This method returns a new series that is each coefficient of self multiplied by *value*. *Value* must be convertible to a float.

        :param other:   A number to multipliy each coefficient of *self* by.

        :returns:       A new instance that is the same class as *self*.
        """

        value_ = float( value )
        c_ls = self.copy( )
        for l, c_l in enumerate( self ) : c_ls[l] *= value
        return( c_ls )

    __rmul__ = __mul__

    def __div__(self, value):
        """
        This method returns a new series that is each coefficient of *self* divided by *value*. *Value* must be convertible to a float.

        :param other:   A number to divide each coefficient of *self* by.

        :returns:       A new instance that is the same class as *self*.
        """

        value_ = float(value)
        c_ls = self.copy()
        for l, c_l in enumerate(self):
            c_ls[l] /= value
        return c_ls 

    def __eq__(self, other):
        """
        This method returns True if *self* and *other* have the same coefficients, domain and lowerIndex, and False otherwise.

        :param other:   An instance of the same class as *self*.

        :returns:       A python boolean.
        """

        return self.areDomainsMutual(other) and self.lowerIndex == other.lowerIndex and self.coefficients == other.coefficients

    def __str__(self):
        """
        This method returns a string representation of the coefficients of *self*.

        :returns:       A python str,
        """

        return ' '.join(["%g" % c_l for c_l in self.__coefficients])

    @property
    def coefficients(self):
        """
        This method returns a reference to *self*'s coefficients.

        :returns:       A python list of floats.
        """

        return self.__coefficients

    @coefficients.setter
    def coefficients(self, coefficients):
        """
        This method sets *self*'s coefficients to *coefficients*.

        :param coefficients:    A list of floats.
        """

        self.__coefficients = list(map(float, coefficients))

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance. 

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps:    Used to dull the lower point for "flat" interpolation.
        :param upperEps:    Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self.toPointwise_withLinearXYs(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps, biSectionMax=biSectionMax)

        return xys1d

    def areDomainsMutual(self, other):
        """
        This method returns True if the domain minimum for *self* and *other* are the same, and if the domain maximum for *self* and *other*
        are the same. It returns False otherwise.

        :param other:       Another 1d function.

        :returns:           A python boolean.
        """

        return self.domainMin == other.domainMin and self.domainMax == other.domainMax

    def evaluate(self, mu):
        """
        This method returns the results of the series evaluated at the mu.

        :param mu:              The domain value.

        :returns:               A python float
        """

        total = 0.0
        for i,c in enumerate(self.__coefficients):
            total += c * self.evaluateBasisFunction(mu, i)
        return total

    def evaluateBasisFunction(self, x, i):
        """
        This is a dummy method that needs to be overwritten by the derived class.

        :param x:       The domain value.
        :param i:       The order of the function to evaluate.
        """

        raise NotImplementedError("Implement in derived classes")

    def checkSameSeriesType( self, other, operator ) :
        """
        This method executes a raise if *self* and *other* not the same type of series.

        :param other:       The instance to check against *self*.
        :param operator:    A python str that represents the operation of the calling method.
        """

        if( not( isinstance( other, Series1d ) ) ) : raise TypeError( 'other of type "%s"' % type( other ) )
        if( self.moniker != other.moniker ) : raise TypeError( 'Cannot %s series %s to series %s' % ( operator, self.moniker, other.moniker ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if len(self.axes) == 0: return
        factors = self.axes.convertUnits( unitMap )
        if( factors[:2] !=[ 1., 1. ] ) :
            self.__domainMin *= factors[1]
            self.__domainMax *= factors[1]
            for l, c_l in enumerate( self ) : self[l] *= factors[0]
            factor = math.pow( factors[1], -self.__lowerIndex )
            for l, c_l in enumerate( self ) :
                if( ( l + self.__lowerIndex ) == 0 ) : factor = 1
                self[l] *= factor
                factor *= 1 / factors[1]
        self.fixValuePerUnitChange( factors )

        if( self.uncertainty ) : self.uncertainty.convertUnits( unitMap )

    def copy( self ) :
        """
        Thie method creates a new series that is a copy of *self*.

        :returns:       An instance that is the same class as *self*.
        """

        return( self.returnAsClass( self, self.__coefficients, index = self.index, outerDomainValue = self.outerDomainValue ) )

    __copy__ = copy

    @property
    def domainMin( self ) :
        """
        This method returns the minimum domain value for *self*.

        :returns:       A python float.
        """

        return( self.__domainMin )

    @property
    def domainMax( self ) :
        """
        This method returns the maximum domain value for *self*.

        :returns:       A python float.
        """

        return( self.__domainMax )

    @property
    def lowerIndex( self ) :
        """
        This method returns the lowerIndex for *self*.

        :returns:       A python int.
        """

        return( self.__lowerIndex )

    @property
    def upperIndex( self ) :
        """
        This method returns the upper index for *self* which is lowerIndex plus the number of coefficients of *self*.

        :returns:       A python int.
        """


        return( self.__lowerIndex + len( self ) )

    def domainUnitConversionFactor( self, unitTo ) :
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def getCoefficientSafely( self, l ) :
        """
        This method returns the (l+1)^th coefficient or 0 if *l* is greater than the number of coefficients of *self*. This is like
        __getitem__ but allows for *l* to be greater than the number of coefficients of *self*.
        """

        if( l >= len( self ) ) : return( 0. )
        return( self.__coefficients[l] )

    def hasData(self):
        """
        This method returns True if the length of self's coefficients is greater than 0.

        :returns:   A python int.
        """

        return len(self.__coefficients) > 0

    def invert( self ) :
        """
        This method returns a series instance that is the mirror of self about x1 = 0 (i.e., negates all odd coefficients).
        For example, for a Legendre series return equilavent of pdf(-mu).

        :returns:       An instance that is the same class as *self*.
        """

        series_ = self.copy( )
        for l in range( 1, len( series_ ), 2 ) : series_.coefficients[l] *= -1
        return( series_ )

    def setData( self, data ) :
        """
        This method sets *self* coefficients to *data*.

        :param data:        A list of floats.
        """

        self.__coefficients = data

    @property
    def rangeMin( self ) :
        """
        This method should return the miminum range of *self*; however, since this depends on the type of polynomial series, this
        method just executes a raise. Ergo, this method needs to be overwritten by the derived class.

        :raises NotImplementedError:    This is always raised.
        """

        raise NotImplementedError( )

    @property
    def rangeMax( self ) :
        """
        This method should return the maxinum range of *self*; however, since this depends on the type of polynomial series, this
        method just executes a raise. Ergo, this method needs to be overwritten by the derived class.

        :raises NotImplementedError:    This is always raised.
        """

        raise NotImplementedError( )

    def toString( self ) :
        """
        This method returns a simple string representation of *self*.

        :returns:       A python str.
        """

        str = [ "%d %16.8g" % ( l, coefficient ) for l, coefficient in enumerate( self.__coefficients ) ]
        return( '\n'.join( str ) )

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

        attributesStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if None in self.fixedDomain(): attributesStr += ' domainMin="%s" domainMax="%s"' % (
                valueFormatter(self.domainMin, significantDigits = significantDigits),
                valueFormatter(self.domainMax, significantDigits = significantDigits) )
        if self.lowerIndex != 0: attributesStr += ' lowerIndex="%s"' % self.lowerIndex

        # FIXME2: converting self.__coefficients to values for printing. Should it be stored as values in the first place?
        coefs = valuesModule.Values(self.__coefficients, valueType = self.valueType)

        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr) ]
        if self.isPrimaryXData() and len(self.axes) > 0:
            XML_strList += self.axes.toXML_strList(indent2, **kwargs)
            XML_strList += coefs.toXML_strList(indent2, **kwargs)
            if self.uncertainty: XML_strList += self.uncertainty.toXML_strList(indent2, **kwargs)
            XML_strList[-1] += '</%s>' % self.moniker
        else:
            XML_strList += coefs.toXML_strList('', **kwargs)
            XML_strList[-1] += '</%s>' % self.moniker
            XML_strList = [ ''.join(XML_strList) ]

        return XML_strList

    @classmethod
    def returnAsClass( cls, self, coefficients, index = None, outerDomainValue = None ) :
        """
        This method returns an instance of *cls* that gets *domainMin*, *domainMax*, *lowerIndex*, *axes*, *label and *valueType* from *self*.
        This is for internal use.

        :param cls:                 The class of the instance to return.
        :param self:                An instance of class *cls* where meta-data not given by the arguments are taken from.
        :param coefficients:        This is the list of coefficients for the 1d polymonial function.
        :param index:               This is the index member use by some xData classes.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is embedded in a high dimensional function.

        :returns:                   An instance of class *cls*.
        """

        return cls(coefficients, self.domainMin, self.domainMax, lowerIndex = self.lowerIndex, axes = self.axes, index = index, 
                valueType = self.valueType, outerDomainValue = outerDomainValue, label = self.label)

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        This method parses *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.

        domainMin, domainMax = cls.fixedDomain()
        if domainMin is None: domainMin = float(extraAttributes['domainMin'])
        extraAttributes.pop('domainMin', None)
        if domainMax is None: domainMax = float(extraAttributes['domainMax'])
        extraAttributes.pop('domainMax', None)

        series = cls([], domainMin=domainMin, domainMax=domainMax, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(series, node, xPath, linkData, **kwargs)

        if len(extraNodes) == 1:
            values = extraNodes.pop()
            values = valuesModule.Values.parseNodeUsingClass(values, xPath, linkData, **kwargs)
            series.coefficients = values                                    # FIXME store self.__coefficients as values instance instead of list?

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()

        return series

    @staticmethod
    def fixedDomain( ) :
        """
        This method does nothing.

        :returns:           (None, None).
        """

        return( None, None )

class LegendreSeries( Series1d ) :
    r"""
    This class represent a polynomial with Legendre functions as the base functions. Ergo, the function :math:`f(\mu)` as:

    ..math::

        f(\mu) = \sum_{l=0}^{\rm lMax} ( l + 0.5 ) \, C_l \, P_l(\mu)

    where :math:`C_l` is the Legendre coefficient for order :math:`l` and is defined as:

    ..math::
        C_l = \int_{-1}^1 d\mu \, P_l(\mu) \, f(\mu)

    Here lMax is the highest Legendre coefficients in the instance, and :math:`P_l(\mu)` is the Legendre polynomial of order *l*.
    """

    moniker = 'Legendre'
    dimension = 1

    def __init__(self, coefficients, domainMin = -1, domainMax = 1, lowerIndex = 0, axes = None, index = None, 
            valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):
        """
        See method :py:func:`Series1d.__init__` for definition of arguments.
        """

        if( lowerIndex != 0 ) : raise ValueError( 'lowerIndex = %s must be 0' )
        if( ( domainMin != -1 ) or ( domainMax != 1 ) ) : raise ValueError( 'domain must be [-1, 1], not [%s, %s]' % ( domainMin, domainMax ) )
        Series1d.__init__( self, coefficients, -1, 1, axes = axes, index = index, valueType = valueType, outerDomainValue = outerDomainValue, 
                label = label)

    def __iter__( self ) :
        """
        This method interates over the coefficients.

        :returns:       A python float.
        """

        length = len( self )
        for i1 in range( length ) :
            yield self[i1]

    def evaluate( self, mu ) :
        r"""
        This method returns the value of *self* at :math:`\mu`.

        :returns:       A python float.
        """

        P = 0.
        for l, c_l in enumerate( self.coefficients ) : P += ( l + 0.5 ) * c_l * Legendre( l, mu, checkXRange = False ) 
        return( P )

    def evaluateBasisFunction(self, mu, l):
        r"""
        This method returns :math:`( l + 0.5 ) \, P_l(\mu)` evaluated at *mu*.

        :param mu:      The :math:`\mu` value.
        :param l:       The order of the function to evaluate.

        :returns:       A python float.
        """

        return ( l + 0.5 ) * Legendre( l, mu, checkXRange = False )

    def integrate(self, domainMin=None, domainMax=None):
        """
        This method returns the integral of *self* from *domainMin* to *domainMax*.

        :param domainMin:   The lower limit of the domain. If None, set to -1.
        :param domainMax:   The upper limit of the domain. If None, set to 1.

        :returns:           A python float.
        """

        if domainMin is None:
            domainMin = -1.0
        if domainMax is None:
            domainMax =  1.0

        sign = 1
        if domainMin > domainMax:
            domainMin, domainMax, sign = domainMax, domainMin, -1
            
        if domainMin < -1.0:
            raise ValueError('domainMin must be greater than or equal to -1: it is %e' % domainMin )
        if domainMax >  1.0:
            raise ValueError('domainMax must be greater than or equal to 1: it is %e' % domainMax )
        if domainMin == domainMax:
            return 0.0
        if domainMin == -1.0 and domainMax == 1.0:
            return self.coefficients[0]

        maxOrder = len(self.coefficients)
        if maxOrder == 0:
            return 0.0

        integral = 0.5 * (domainMax - domainMin) * self.coefficients[0]

        if maxOrder > 1:
            integral += 0.75 * (domainMin + domainMax) * (domainMax - domainMin) * self.coefficients[1]

            P_l_m1_1 = domainMin
            P_l_m1_2 = domainMax
            P_l_1 = 0.5 * (3.0 * domainMin * domainMin - 1.0)
            P_l_2 = 0.5 * (3.0 * domainMax * domainMax - 1.0)
            for order in range(2, maxOrder):
                P_l_p1_1 = Legendre(order + 1, domainMin)
                P_l_p1_2 = Legendre(order + 1, domainMax)
                integral += 0.5 * (P_l_p1_2 - P_l_p1_1 + P_l_m1_1 - P_l_m1_2) * self.coefficients[order]
                P_l_m1_1 = P_l_1
                P_l_m1_2 = P_l_2
                P_l_1 = P_l_p1_1
                P_l_2 = P_l_p1_2
        
        return sign * integral

    def isIsotropic( self ) :
        """
        This method returns True if self is isotropic and False otherwise.

        :returns:       A python boolean.
        """

        for coefficient in self[1:] :
            if( coefficient != 0. ) : return( False )
        return( True )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        This method returns an instance of class *self* that is normalized.
        The *dimension* argument is ignored as it is only here to be compatable with other normalize methods.

        :param insitu:      If True, *self* is normalizes and returned; otherwise, a copy of *self* is normalized and returned.
        :param dimension:   This argument is ignored.

        :returns:       
        """

        norm = self.coefficients[0]
        coefficients = [ coefficient / norm for coefficient in self ]
        if( insitu ) :
            self.setData( coefficients )
            return( self )
        return( self.returnAsClass( self, coefficients ) )

    def toLinearXYsClass( self ) :
        """
        This method returns the class that is used to represent *self* as a 1d pointwise function, which is always :py:class:`XYs1dModule.XYs1d`.

        returns:        The class :py:class:`XYs1dModule.XYs1d`.
        """

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        r"""
        This method returns a representation of *self* as an :py:class:`XYs1dModule.XYs1d` instance (i.e., a pointwise :math:`P(\mu)`). 
        The accuracy of the reconstruction is determined by kwargs['accuracy'].

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : 1e-3, 'biSectionMax' : 16 } )
        accuracy = arguments['accuracy']
        biSectionMax = arguments['biSectionMax']

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1

        try :
            L = Legendre_C.Series( self.coefficients )
            P = L.toPointwiseLinear( accuracy, biSectionMax = biSectionMax, checkForRoots = True )
        except :
            P, n = [], 400
            for i in range( n ) :
                mu = -1. + ( 2. * i ) / n
                P.append( [ mu, self.evaluate( mu ) ] )
            P.append( [ 1., self.evaluate( 1. ) ] )
        axes = axesModule.Axes(2)
        unit = self.getAxisUnitSafely( 0 )
        axes[0] = axesModule.Axis( 'P(mu)', 0, unit )
        axes[1] = axesModule.Axis( 'mu', 1, '' )
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = axes, outerDomainValue = self.outerDomainValue )
        return( P.thin( accuracy = accuracy ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method returns the results of calling :py:func:`toPointwiseLinear`.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        return( self.toPointwiseLinear( **kwargs ) )

    @staticmethod
    def fixedDomain( ) :
        """
        This method returns the domain of *self* which is always [-1, 1].A

        :returns:       The tuple ( -1, 1 ).
        """

        return( -1, 1 )

    @staticmethod
    def fromXYs1d( xys1d, maxOrder ) :
        r"""
        This method returns an instance of :py:class:`LegendreSeries` that represents *xys1d* as:

        ..math::
            C_l = \int_{-1}^1 d\mu \, P_l(\mu) \, f(\mu)

        where *xys1d* is the poinwise representation of :math:`f(\mu)`.

        :param xys1d:       An instance of :py:class:`XYs1dModule.XYs1d`.
        :param maxOrder:    The maximum Legendre order of the returned instance.

        :returns:           An instance of :py:class:`LegendreSeries`
        """

        if( xys1d.domainMin != -1 ) : raise ValueError( "Domain min must be -1. It is %e" % xys1d.domainMin )
        if( xys1d.domainMax !=  1 ) : raise ValueError( "Domain max must be 1. It is %e" % xys1d.domainMax )

        axes = xys1d.axes
        return( LegendreSeries( Legendre_C.from_pointwiseXY_C( xys1d.nf_pointwiseXY, maxOrder ), axes = axes ) )

class Polynomial1d( Series1d ) :
    r"""
    This class represent a polynomial with as a simple power series in the domain variable. This is, for the
    function :math:`f(x)` the polynomial :math:`\sum_{i=0}^{\rm iMax} C_i \, x^i`.

        f(\mu) = \sum_{l=0}^{\rm lMax} ( l + 0.5 ) \, C_l \, P_l(\mu).
    """

    moniker = 'polynomial1d'
    dimension = 1

    def __init__(self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None, index = None, 
            valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):
        """
        See method :py:func:`Series1d.__init__` for definition of arguments.
        """

        Series1d.__init__( self, coefficients, domainMin, domainMax, lowerIndex = lowerIndex, axes = axes, index = index, 
                valueType = valueType, outerDomainValue = outerDomainValue, label = label)

    def __findExtrema( self ):
        """
        This method returns the list of domain values that are the roots of the polynomial.

        :returns:   A python list of numbers that may be complex.
        """

        from numpy.polynomial.polynomial import Polynomial
        p1 = Polynomial(self.coefficients, domain=(self.domainMin, self.domainMax),
            window=(self.domainMin, self.domainMax))
        pprimeroots = p1.deriv().roots()
        return [self.domainMin] + list(pprimeroots) + [self.domainMax]

    def evaluate( self, x ) :
        """
        This method returns the value of *self* evaluated at *x*.

        :param x:       The domain value to evaluate *self* at.

        :returns:       A python float.
        """

        P = 0.
        for c_l in reversed( self.coefficients ) : P = c_l  + x * P
        return( P )

    def evaluateBasisFunction(self, x, i):
        """
        This method returns that basis function evaluated at *x* which is :math:`x^i`.

        :param x:       The domain value.
        :param i:       The power the domain value.

        :returns:       A python float.
        """

        return pow(x,i)

    @property
    def rangeMin( self ):
        """
        This method returns the minimum range value for *self*.

        :returns:       A python float.
        """

        extrema = self.__findExtrema()
        return min([self.evaluate(x) for x in extrema])

    @property
    def rangeMax( self ):
        """
        This method returns the maximum range value for *self*.

        :returns:       A python float.
        """

        extrema = self.__findExtrema()
        return max([self.evaluate(x) for x in extrema])

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        This method sets the domain minimum and maximum per the arguments.

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

    def toLinearXYsClass( self ) :
        """
        This method returns the class that is used to represent *self* as a 1d pointwise function, which is always :py:class:`XYs1dModule.XYs1d`.

        returns:        The class :py:class:`XYs1dModule.XYs1d`.
        """

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        r"""
        This method returns a representation of *self* as an :py:class:`XYs1dModule.XYs1d` instance (i.e., a pointwise :math:`y(x)`).
        The accuracy of the reconstruction is determined by kwargs['accuracy'].  Currently, accuracy is not fully implemented.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : 1e-3, 'biSectionMax' : 16 } )
        accuracy = arguments['accuracy']
        biSectionMax = arguments['biSectionMax']

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1

        P, n = [], 1000
        for i in range( n + 1 ) :
            x = ( ( n - i ) * self.domainMin + self.domainMax * i ) / n
            P.append( [ x, self.evaluate( x ) ] )
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = self.axes, outerDomainValue = self.outerDomainValue )
        return( P.thin( accuracy = accuracy ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method returns the results of calling the method :py:func:`toPointwiseLinear`.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        return( self.toPointwiseLinear( **kwargs ) )

class LinearSpline1d( Series1d ) :
    r"""
    This class is a simple linear spline.  It basically wraps the XYs1d class.
    Linear interpolation uses the linear spline or the "hat" basis.  The first basis function looks like this::

        |\
        | \
        |  \
        ----

    The second basis function looks like this::

          /\
         /  \
        /    \
        ------

    and so on.  Together :math:`y(x) = \sum_i B_i(x) y_i` such that :math`y_i = y(x_i)`.
    """

    moniker = 'linearSpline1d'
    dimension = 1

    def __init__(self, xdata, ydata, axes, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue= None, label = None):
        """
        See method :py:func:`Series1d.__init__` for definition of arguments.
        """

        if( len( xdata ) != len( ydata ) ) : raise ValueError( "Number of x and y values not equal." )

        Series1d.__init__(self, ydata, xdata[0], xdata[-1], lowerIndex = 0, axes = axes, index = index, valueType = valueType, 
                outerDomainValue = outerDomainValue, label = label)

        self.axes = axes
        self.basis = XYs1dModule.XYs1d(axes=self.axes, data=list(zip(xdata, ydata)), interpolation=enumsModule.Interpolation.linlin)

    def evaluateBasisFunction( self, x, i ) :
        """
        This method returns the value of the basis function of order *i* evaluated at *x*.

        :param x:       The domain value.
        :param i:       The order of the function to evaluate.

        :returns:       A python float.
        """

        self.basis[i] = (self.basis[i][0], 1.0)
        result=self.basis.evaluate(x)
        self.basis[i] = (self.basis[i][0], 0.0)
        return result

    def toLinearXYsClass( self ) :
        """
        This method returns the class that is used to represent *self* as a 1d pointwise function, which is always :py:class:`XYs1dModule.XYs1d`.

        returns:        The class :py:class:`XYs1dModule.XYs1d`.
        """

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        """
        This method returns a representation of *self* as an :py:class:`XYs1dModule.XYs1d` instance (i.e., a pointwise :math:`y(x)`).
        The accuracy of the reconstruction is determined by kwargs['accuracy'].  
        Basically it just copy the basis function widget and put the y values back.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        result = self.basis.copy()
        for index, coefficient in enumerate( self.coefficients ) :
            result[index] = ( result[index][0], coefficient )

        return result

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method returns the restuls of Calling :py:class:`toPointwiseLinear` method.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        return( self.toPointwiseLinear( **kwargs ) )
