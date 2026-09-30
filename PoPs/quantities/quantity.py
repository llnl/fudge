# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Contains the abstract base class for storing a physical quantity (e.g., mass, spin or halflife).
Also defines specific classes for various types of quantity, i.e. integer, float, fraction or string.
"""

import sys
import abc
import fractions

from pqu import PQU as PQUModule

from .. import misc as miscModule
from .. import suite as suiteModule

from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule

class Quantity( miscModule.ClassWithLabelKey, abc.ABC ) :
    """
    This class is used to represent a physical quantity (e.g., mass, spin or halflife).
    A quantity has required members label, value and unit and optional members documentation and uncertainty.
    """

    moniker = 'quantity'
    __valueType = None
    baseUnit = None

    def __init__( self, label, value, unit, documentation = '' ) :
        """
        :param label: label for this quantity. Must be unique within the containing suite
        :param value: quantity value
        :param unit: unit (string)
        :param documentation: documentation specific to this quantity
        """

        miscModule.ClassWithLabelKey.__init__( self, label )

        self.value = value

        self.unit = unit

        if( not( isinstance( documentation, str ) ) ) : raise TypeError( 'documentation must be a string' )
        self.__documentation = documentation

        self.uncertainty = None

    def __str__( self ) :

        return( "%s %s" % ( self.value, self.unit ) )

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, value ) :

        if( not( isinstance( value, self.valueType ) ) ) : raise TypeError( 'Invalid value type must be a "%s"' % self.valueType )
        self.__value = value

    @property
    def valueType( self ) :

        return( self.__valueType )

    @property
    def unit( self ) :

        return( self.__unit )

    @unit.setter
    def unit( self, unit ) :

        if( isinstance( unit, str ) ) : unit = stringToPhysicalUnit( unit )
        if( not( isinstance( unit, PQUModule.PhysicalUnit ) ) ) : raise TypeError( 'unit must be a PQU.PhysicalUnit' )
        if( self.baseUnit is not None ) :
            if( not( self.baseUnit.isCompatible( unit ) ) ) : raise ValueError( 'unit "%s" not compatible with baseUnit "%s"' % ( unit, self.baseUnit ) )
        self.__unit = unit

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( _uncertainty is not None ) :
            if( not( isinstance( _uncertainty, uncertaintyModule.Uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )

        self.__uncertainty = _uncertainty
        if( self.__uncertainty is not None ) : self.__uncertainty.setAncestor( self )

    @property
    def documentation( self ) :

        return( self.__documentation )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        if( abs( factor - 1 ) > ( 4 * sys.float_info.epsilon ) ) :
            NotImplementedError( 'Conversion of units for quantity of "%s" not implemented: from unit "" to unit "%s"' %
                    ( self.moniker, self.unit, unit ) )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _quantity = self.__class__( self.label, self.value, self.unit, self.documentation )
        if( self.__uncertainty is not None ) : _quantity.uncertainty = self.__uncertainty.copy( )
        return( _quantity )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        if( not( self.unit.isDimensionless( ) ) ) : attributes += ' unit="%s"' % self.unit

        ending = '>'
        if( ( self.documentation == '' ) and ( self.uncertainty is None ) ) : ending = '/>'
        XMLStringList = [ '%s<%s label="%s" value="%s"%s%s' % 
                ( indent, self.moniker, self.label, self.valueToString( ), attributes, ending ) ]

        if( ending == '>' ) :
            if( self.documentation != '' ) :
                XMLStringList.append( '%s<documentation>' % indent2 )
                XMLStringList.append( '%s%s</documentation>' % ( indent3, self.documentation ) )
            if( self.uncertainty is not None ) : XMLStringList += self.uncertainty.toXML_strList( indent = indent2, **kwargs )
            XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def valueToString( self, precision = 12 ) :

        return( "%s" % self.value )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        documentation = ''
        for child in element :
            if( child.tag == 'uncertainty' ) :
                self.uncertainty = uncertaintyModule.Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)
#            elif( child.tag 'documentation' ) :
#                if( child.tag == 'documentation' ) : documentation = child.findtext( )
            else :
                raise ValueError( 'child element with tag "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( '%s[@label="%s"]' % ( element.tag, element.get( 'label' ) ) )

        attributes = ( 'label', 'value', 'unit' )
        for attributeName in element.keys() :
            if( attributeName not in attributes ) : raise ValueError( 'attribute = "%s" not allowed' % attributeName )

        value = cls.toValueType( element.get('value') )
        unit = stringToPhysicalUnit( element.get( 'unit', '' ) )

        self = cls( element.get('label'), value, unit )
        xPath.pop()

        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class String( Quantity ) :
    """
    This is an abstract base class for string quantities.
    """

    moniker = 'string'
    __valueType = str

    @property
    def valueType( self ) :

        return( self.__valueType )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class Number( Quantity, abc.ABC ) :
    """
    This is an abstract base class for numberic quantities. This class adds the pqu and float methods.
    """

    def pqu( self, unit = None ) :
        """
        Returns a PQU instance of self's value in units of unit. If unit is None, self's unit is used.
        """

        uncertainty = None
        if self.uncertainty is not None:
            uncertainty = self.uncertainty.form.value
            if uncertainty.relation == uncertaintyModule.Quantity.absolute:
                uncertainty = uncertainty.value
            elif uncertainty.relation == uncertaintyModule.Quantity.relative:
                uncertainty = uncertainty.value * self.value
            elif uncertainty.relation == uncertaintyModule.Quantity.percent:
                uncertainty = uncertainty.value * self.value / 100.0
            else:
                raise ValueError('Unsupported uncertainty relation "%s".' % uncertainty.relation)

        pqu = PQUModule.PQU( self.value, self.unit, uncertainty = uncertainty )
        if( unit is not None ) : pqu.convertToUnit( unit )
        return( pqu )

    def float( self, unit ) :
        """
        Returns a float instance of self's value in units of unit.
        """

        if( not( isinstance( unit, ( str, PQUModule.PhysicalUnit ) ) ) ) : raise TypeError( 'unit argument must be a str or a PQU.PhysicalUnit.' )
        return( float( self.pqu( unit ) ) )

class Integer( Number ) :
    """
    This class is used to represent a (physical) quantity whose value must be an integer.
    """

    moniker = 'integer'
    __valueType = int

    @property
    def valueType( self ) :

        return( self.__valueType )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class Double( Number ) :
    """
    This class is used to represent a (physical) quantity whose value must be a float.
    """

    moniker = 'double'
    __valueType = float

    def __init__( self, label, value, unit, documentation = '' ) :

        if( isinstance( value, ( int, fractions.Fraction ) ) ) : value = self.valueType( value )
        Number.__init__( self, label, value, unit, documentation = documentation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        unit = stringToPhysicalUnit( unit )
        self.value = self.value * factor
        self.unit = unit

    def valueToString( self, precision = 12 ) :

        return( PQUModule.floatToShortestString( self.value, min( max( 0, precision ), 17 ), keepPeriod = True ) )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class Fraction( Number ) :
    """
    This class is used to represent a (physical) quantity whose value must be a rational number
    (e.g., a fraction like 1/2, 3, 5/6).
    """

    moniker = 'fraction'
    __valueType = fractions.Fraction

    def __init__( self, label, value, unit, documentation = '' ) :

        if( isinstance( value, ( int, str ) ) ) : value = self.valueType( value )
        Number.__init__( self, label, value, unit, documentation = documentation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class Suite( suiteModule.Suite, abc.ABC ) :

    _allowedClasses = []

    def __init__( self ) :

        suiteModule.Suite.__init__( self, self._allowedClasses )

class NumberSuite( Suite, abc.ABC ) :
    """
    This is an abstract base class for a number suite. This class adds the pqu and float methods.
    """

    def pqu( self, unit = None ) :
        """
        Returns a PQU instance of self's recommended value in units of unit. If unit is None, self's unit is used.
        :param unit: desired unit (string)
        """

        return( self[0].pqu( unit = unit ) )

    def float( self, unit ) :
        """
        Returns a float instance of self's recommended value in units of unit.
        :param unit: desired unit (string)
        """

        return( self[0].float( unit ) )

def stringToPhysicalUnit( string ) :

    return( PQUModule._findUnit( string ) )
