# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# Notes.
# Need sphinx 1.1 or higher to use ":private-members:" parameter in PQU.rst. This is why :py:func:`PQU._getOtherAsPQU` do not display properly.

"""
---------------------------
Introduction
---------------------------

This module contains classes and functions for representing and manipulating a value with units and uncertainty, 
herein called a "Physical Quantity with Uncertainty" or :py:class:`PQU`. For example, if the distance between start and 
finish lines is measured to be '100.1 m' with an uncertainty of '0.4 m', it can be inputted to :py:class:`PQU` as the string
'100.1(4) m' or '100.1 +/- 0.4 m' (both are allowed input forms as well as several others - see :py:class:`PQU`). 

From this module, most users should only need the :py:class:`PQU` class which stores a :py:class:`PQU` object and
supports common math operations (e.g., addition, subtraction, multiplication) including the operation on 
the units and uncertainty. 
For example, if a person races between the start and finish lines in a time of '16.3 +/- 0.1 s', the :py:class:`PQU` class can 
be used to determine the person's speed as:

>>> from pqu import PQU
>>> distance = PQU.PQU( '100.1(4) m' )
>>> time = PQU.PQU( '16.3 +/- 0.1 s' )
>>> speed = distance / time
>>> print(speed)
6.14 +/- 0.04 m/s
>>> print(speed.copyToUnit( 'mi/h' ))
13.7 +/- 0.1 mi/h

As this example shows, in addition to calculating the speed and unit, the :py:class:`PQU` class also propagated the significant digits and the uncertainty. In this example
the uncertainty was propagated using Goodman's expression for uncorrelated values (see `Uncertainty propagation`_). The
propagation of significant digits explains why the printed speed has only 3 digits. The following Python code illustrates
significant digits for the speed calculation:

>>> print(100.1 / 16.3)
6.14110429448
>>> print(float( speed ))
6.14110429448
>>> print(speed.info( ))
value = 6.14110429447852724e+00, significantDigits = 3, order = 0, isPercent = False, unit = "m/s"
uncertainty = value = 4.49629906071862956e-02, significantDigits = 1, order = -2, isPercent = False

As can be seen, the internal speed value is as expected. However, since the :py:class:`PQU` class (with the help
of the :py:class:`PQU_float` class) tracks a value's significant digits, when a :py:class:`PQU` instance is printed (via the __str__
method), only at most 'significantDigits' are displayed. The allowed range for 'significantDigits' is
1 to sys.float_info.dig + 1 inclusive. For addition and subtraction, the member 'order' is also required and
is why the following output is the same for both print statements:

>>> print(distance)
100.1(4) m
>>> print(distance + '100.1(4) mum')        # note 'mum' is micrometer
100.1 +/- 0.4 m

How is 'significantDigits' determined? That depends on how :py:class:`PQU` is called. If a string without
uncertainty is entered as the only argument then 'significantDigits' is the number of digits in the string
(ignoring leading '0'). Some examples are:

>>> print(PQU.PQU( '12.345' ).info( ))
value = 1.23450000000000006e+01, significantDigits = 5, order = 1, isPercent = False, unit = ""
uncertainty = 
>>> print(PQU.PQU( '12.34500' ).info( ))
value = 1.23450000000000006e+01, significantDigits = 7, order = 1, isPercent = False, unit = ""
uncertainty = 
>>> print(PQU.PQU( '12.34500e-12' ).info( ))
value = 1.23450000000000004e-11, significantDigits = 7, order = -11, isPercent = False, unit = ""
uncertainty = 
>>> print(PQU.PQU( '0012.34500e-12' ).info( ))
value = 1.23450000000000004e-11, significantDigits = 7, order = -11, isPercent = False, unit = ""
uncertainty = 
>>> print(PQU.PQU( '00.0012' ).info( ))
value = 1.19999999999999989e-03, significantDigits = 2, order = -3, isPercent = False, unit = ""
uncertainty = 

If the string has an uncertainty, then it is also used in calculating 'significantDigits'.
Some examples are (note - these are the same as the last examples, with uncertainties added):

>>> print(PQU.PQU( '12.345 +/- 1e-8' ).info( ))
value = 1.23450000000000006e+01, significantDigits = 10, order = 1, isPercent = False, unit = ""
uncertainty = value = 1.00000000000000002e-08, significantDigits = 1, order = -8, isPercent = False
>>> print(PQU.PQU( '12.345 +/- 1e-8' ))
12.34500000 +/- 1.e-8
>>> print(PQU.PQU( '12.34500 +/- 0.12' ).info( ))
value = 1.23450000000000006e+01, significantDigits = 4, order = 1, isPercent = False, unit = ""
uncertainty = value = 1.19999999999999996e-01, significantDigits = 2, order = -1, isPercent = False
>>> print(PQU.PQU( '12.34500 +/- 0.12' ))
12.35 +/- 0.12
>>> print(PQU.PQU( '12.34500e-12(32)' ).info( ))
value = 1.23450000000000004e-11, significantDigits = 7, order = -11, isPercent = False, unit = ""
uncertainty = value = 3.20000000000000023e-16, significantDigits = 2, order = -16, isPercent = False
>>> print(PQU.PQU( '12.34500e-12(32)' ))
1.234500e-11(32)
>>> print(PQU.PQU( '0012.34500e-12 +/- 32e-15' ).info( ))
value = 1.23450000000000004e-11, significantDigits = 5, order = -11, isPercent = False, unit = ""
uncertainty = value = 3.20000000000000025e-14, significantDigits = 2, order = -14, isPercent = False
>>> print(PQU.PQU( '0012.34500e-12 +/- 32e-15' ))
1.2345e-11 +/- 3.2e-14
>>> print(PQU.PQU( '00.0012 +/- 0.000002' ).info( ))
value = 1.19999999999999989e-03, significantDigits = 4, order = -3, isPercent = False, unit = ""
uncertainty = value = 1.99999999999999991e-06, significantDigits = 1, order = -6, isPercent = False
>>> print(PQU.PQU( '00.0012 +/- 0.000002' ))
1.200e-3 +/- 2.e-6

The :py:class:`PQU` constructor (i.e., its __init__ method) allows various input options for creating an instance (see :py:class:`PQU`).

A :py:class:`PQU` object has three main members which are:

    * *value* which is stored as a :py:class:`PQU_float` instance,
    * *unit* which is stored as a :py:class:`PQU_uncertainty` instance,
    * *uncertainty* which is stored as a :py:class:`PQU_uncertainty` instance.
    
------------------------------------------------------------------------------------------
Units and conversions
------------------------------------------------------------------------------------------

This module uses the SI units along with units for angle and solid angle as its base units. The base units are:

    +-----------+-------+---------------------------+
    | Unit      | symbol| Measure                   |
    +===========+=======+===========================+
    | meter     | 'm'   | length                    |
    +-----------+-------+---------------------------+
    | kilogram  | 'kg'  | mass                      |
    +-----------+-------+---------------------------+
    | second    | 's'   | time                      |
    +-----------+-------+---------------------------+
    | Ampere    | 'A'   | electrical current        |
    +-----------+-------+---------------------------+
    | Kelvin    | 'K'   | thermodynamic temperature |
    +-----------+-------+---------------------------+
    | mole      | 'mol' | amount of substance       |
    +-----------+-------+---------------------------+
    | candela   | 'cd'  | luminous intensity        |
    +-----------+-------+---------------------------+
    | radian    | 'rad' | angle                     | 
    +-----------+-------+---------------------------+
    | steradian | 'sr'  | solid angle               |
    +-----------+-------+---------------------------+

All other units are stored with these units as a base and a conversion factor. 
For example, the unit for foot ('ft') is stored as a meter with the conversion factor of 0.3048.
Any :py:class:`PQU` instance can be represented as a combination of each of these units with an associated power for each unit.
As example, a force has the base units 'kg' and 'm' to power 1 as well as 's' to power -2 (represented as the 
string 'kg * m / s**2').  Two units are said to be compatible if they have the same power for each
base unit. Hence, 'ft' is compatible with 'm' but not 'kg' or 's'. Furthermore, the electron-Volt ('eV') is compatible
with Joule ('J') and Watt-second ('W * s') but not Watt ('W').

The :py:class:`PQU` package has many defined units with prefixes (see `Defined prefixes and units`_). In general, the
PQU methods that operate on :py:class:`PQU` objects do not attempt to simplify the units, even if the result is dimensionless.
Note, to get rid of dimensionless units use the method :py:func:`PQU.simplify`
For example, consider the division of '3.2 m' by 11.2 km:

>>> from pqu import PQU
>>> x = PQU.PQU( '3.2 m' )
>>> y = PQU.PQU( '11.2 km' )
>>> slope = x / y
>>> print(slope)
0.29 m/km
>>> dl = slope.copyToUnit( "" )
>>> print(dl)
2.9e-4
>>> slope.simplify()
PQU( "2.9e-4" )

Here the method :py:func:`PQU.copyToUnit` is used to convert the units into a dimensionless form. Here is another
example showing the use of the :py:func:`PQU.copyToUnit` method:

>>> mass = PQU.PQU( "4.321 g" )
>>> speed = PQU.PQU( "1.234 inch / mus"  )  
>>> energy = mass * speed**2
>>> print(energy)
6.580 inch**2*g/mus**2
>>> energy_joules = energy.copyToUnit( "J" )
>>> print(energy_joules)
4.245e6 J

The method :py:func:`PQU.inUnitsOf` is useful for returning a physical quantity in units of descending compatible units. 
For example, :py:func:`PQU.inUnitsOf` will convert '3123452.12 s' into days, hours, minutes and seconds, or just hours and seconds as:

>>> t = PQU.PQU( '3123452.12 s' )
>>> t.inUnitsOf( 'd', 'h', 'min', 's' )
(PQU( "36.0000000 d" ), PQU( "3.000000 h" ), PQU( "37.0000 min" ), PQU( "32.12 s" ))
>>> t.inUnitsOf( 'h',  's' )
(PQU( "867.000000 h" ), PQU( "2252.12 s" ))

Also see the methods :py:func:`PQU.convertToUnit` and :py:func:`PQU.getValueAs`.

------------------------------------------------------------------------------------------
Changing non-hardwired constants
------------------------------------------------------------------------------------------

PQU has a set of defined constants that are hardwired and a set that are not hardwired. The non-hardwired
constants reside in the module pqu_constants.py and are import by :py:class:`PQU` when it is first loaded. It is 
possible to change a non-hardwired constant by loading the pqu_constants.py module before :py:class:`PQU` is imported, redefining
the constant and then importing :py:class:`PQU`. For example, as of this writing the elementary charge is defined as '1.60217653e-19 * C':

>>> from pqu import PQU
>>> eV = PQU.PQU( 1, 'eV' )
>>> print(eV)
1. eV
>>> print(eV.copyToUnit( 'J' ))
1.60217653e-19 J

The following lines change the elementary charge to '1.6 * C':

>>> from pqu import pqu_constants
>>> pqu_constants.e = '1.6 * C'
>>> from pqu import PQU                 # This import does not work with doctest as PQU was previously imported.
>>> my_eV = PQU.PQU( 1, 'eV' )
>>> print(my_eV)
1. eV
>>> print(my_eV.copyToUnit( 'J' ))
1.6 J

The python 'reload' function can also be used as:

>>> from pqu import PQU
>>> eV = PQU.PQU( 1, 'eV' )
>>> print(eV)
1. eV
>>> print(eV.copyToUnit( 'J' ))
1.60217653e-19 J
>>> from pqu import pqu_constants
>>> pqu_constants.e = '1.6 * C'
>>> reload( PQU )
<module 'pqu.PQU' from 'pqu/PQU.pyc'>
>>> my_eV = PQU.PQU( 1, 'eV' )
>>> print(my_eV)
1. eV
>>> print(my_eV.copyToUnit( 'J' ))
1.6 J
>>> pqu_constants.e = '1.60217653e-19 * C'      # Needs to be the same as in the module pqu_constants.py for stuff below to pass doctest.
>>> reload( PQU )
<module 'pqu.PQU' from 'pqu/PQU.pyc'>

------------------------------------------------------------------------------------------
Uncertainty propagation
------------------------------------------------------------------------------------------

This section describes the propagation of uncertainties for the binary operators '+', '-', '*' and '/' as well as the power operator
(i.e., '**'). Firstly, a comment on the limits of uncertainty propagation. Uncertainties were added mainly to support the inputting
and outputting (i.e., converting to and from a string) of a physical quantity with an associated uncertainty. As can be seen from the documentation below,
the supported propagations of uncertainties are limited to simple propagation rules.

In the following discussion, Q1, Q2 and Q3 will be :py:class:`PQU` instances with values
V1, V2 and V3 respectively and uncertainties dQ1, dQ2 and dQ3 respectively. In addition, n1 will be a number (i.e., 1.23).

For the binary operators, all numbers are converted to a :py:class:`PQU` instance with uncertainty of 0. For example, the following
are all equivalent:

>>> from pqu import PQU
>>> Q1 = PQU.PQU( '12.3 +/- .2 m' )
>>> print(Q1 * 0.6)
7.4 +/- 0.1 m
>>> print(Q1 * PQU.PQU( '0.6' ))
7.4 +/- 0.1 m
>>> print(PQU.PQU( '0.6' ) * Q1)
7.4 +/- 0.1 m

All binary operators assume that the two operands are uncorrelated except when they are the same instance. When
the two operands are the same instance, 100% correlation is assumed. For example,

>>> print(Q1 - Q1)                      # Operands are the same instance
0. m
>>> print(Q1 - PQU.PQU( '12.3 +/- .2 m' ))  # Operands are not the same instance
0.0 +/- 0.3 m
>>> print(Q1 + Q1)                      # Operands are the same instance
24.6 +/- 0.4 m
>>> print(Q1 + PQU.PQU( '12.3 +/- .2 m' ))  # Operands are not the same instance
24.6 +/- 0.3 m
>>> print(Q1 * Q1)                      # Operands are the same instance
151. +/- 5. m**2
>>> print(Q1 * PQU.PQU( '12.3 +/- .2 m' ))  # Operands are not the same instance
151. +/- 3. m**2
>>> print(Q1 / Q1)                      # Operands are the same instance
1.
>>> print(Q1 / PQU.PQU( '12.3 +/- .2 m' ))  # Operands are not the same instance
1.00 +/- 0.02

There is one exception (depending on your divide-by-zero belief) to this rule of same instance and that happens when the 
value of the operand is 0. Divide-by-zero always executes a raise of 'ZeroDivisionError', even for 'a / a'.

>>> Q2 = PQU.PQU( "0.00 +/- 0.02" )
>>> print(Q2 / "2 +/- 0.5")
0.00 +/- 0.01
>>> print(Q2 / Q2)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "PQU.py", line 1189, in __div__
    extraFactor = uncertaintySelf * uncertaintyOther
ZeroDivisionError: float division by zero

For the rest of this discussion, the two binary operands will be assumed to not be the same instance.
For addition and subtraction with :py:class:`PQU` instances Q1 and Q2 as operands the uncertainty is propagated as dQ3 = sqrt( dQ1**2 + dQ2**2 ).
For multiplication with :py:class:`PQU` instances Q1 and Q2 as operands the uncertainty is propagated using Goodman's expression with
uncorrelated values (i.e., as dQ3 = sqrt( ( V2 * dQ1 )**2 + ( V1 * dQ2) **2 + ( dQ1 * dQ2 )**2 )). For division, the expression 'Q1 / Q2'
is converted to 'Q1 * Q3' with 'Q3 = PQU( 1 / V2, unit = 1 / U2, uncertainty = dQ2 / ( V2 * V2 ) )' where U2 is the unit of Q2.

For the power method, simple propagation is perform which is only valid if the uncertainty is small relative to the value. 
This is, for 'Q2 = pow( Q1, n1 )' the uncertainty for Q2 is 'dQ2 = n1 * V1**(n1-1) * dQ1'.
Note, the 'sqrt' method is equivalent to the 'pow' method with power of 0.5.

For the trigonometry methods, no propagation of uncertainty is performed.

------------------------------------------------------------------------------------------
PQU methods by functionality
------------------------------------------------------------------------------------------

This section classifies some of the methods in the :py:class:`PQU` class by their function.

For :py:class:`PQU` methods that require a second operand (e.g., '+', '/'), the second operand will be passed to the
staticmethod :py:func:`PQU._getOtherAsPQU` to determine if it is a suitable object. A suitable object includes any string that
is a valid :py:class:`PQU`. For example: :py:func:`PQU.toString`

>>> print(PQU.PQU( '12.345' ) + "3.23")
15.58
>>> print(PQU.PQU( '12.345 mm' ) + "3.23m")
3242. mm

**isCompatible**

    Returns **True** if *self*'s unit is compatible with the units of the argument and **False** otherwise.

**isPercent, isDimensionless, isLength, isMass, isTime, isCurrent, isTemperature, isMoles, isLuminousIntensity,
isAngle, isSolidAngle, isEnergy, isSpin, isCharge**

    These methods return **True** if *self* is compatible with the unit implied by the method name.

**__float__**

    This method is called when the a :py:class:`PQU` is the argument to the 'float' function (e.g., float( pqu )).
    This method returns a python float for the value of the :py:class:`PQU` instance. As example:

>>> pqu = PQU.PQU( '13.7 +/- 0.1 mi/h' )
>>> value = float( pqu )
>>> print(type( value ), value)
<type 'float'> 13.7

The :py:class:`PQU` class has methods for the following arithmetic operations. Unless stated, a new :py:class:`PQU` instance is returned.

**__abs__, __neg__**

    These are the standard unitary operators that only operate on the :py:class:`PQU`'s value.

**__add__, __radd__, __iadd__, __sub__, __rsub__, __isub__**

    These are the standard binary operators for addition and subtraction. For these operators, the other
    operand is passed to :py:func:`PQU._getOtherAsPQU`. The other operand must have units compatible with *self*.
    As expected, the __iadd__ and __isub__ methods modify *self* and return it.

**__mul__, __rmul__, __imul__, __div__, __rdiv__, __idiv__**

    These are the standard binary operators for multiplication and division. For these operators, the other
    operand is passed to :py:func:`PQU._getOtherAsPQU`. Any unit for other is allowed.
    As expected, the __imul__ and __idiv__ methods modify *self* and
    return it. All these methods call the __mul__ method; except, when a division happens, and '*self* is other'
    and the denominator is not 0.
    That is, the expression 'a / a' is handled as a special case when float( a ) != 0. In this case, a 
    dimensionless :py:class:`PQU` is returned with value 1 and no uncertainty.

    If *self* is not other, the uncertainty is propagated using Goodman's expression for uncorrelated values.
    Otherwise, 100% correlated uncertainty is assumed.

**__pow__**

    This method raises - not in the python sense but in the mathematical sense (i.e., x**y) -
    *self* to a power. Standard uncertainty propagation is performed which is only valid when
    the uncertainty is small compared to the value.  A quantity can be raised to a non-integer power 
    only if the result can be represented by integer powers of the base units.

**__eq__, __ne__, __lt__, __le__, __gt__, __ge__, compare, equivalent**

    These methods compare *self* to another object. For these operators, the other
    object is passed to :py:func:`PQU._getOtherAsPQU` and must have the same unit
    as *self*. The first 6 methods all call the :py:func:`compare` method
    with epsilonFactor = 5. All methods except 'equivalent' compare 
    *self*'s and other's values and units only. 
    The 'equivalent' method also considers *self*'s and other's uncertainty.

**convertToUnit, copyToUnit, inUnitsOf, getUnitAsString, getValueAs, getUncertaintyValueAs**

    These methods change a :py:class:`PQU`'s unit or return a new instance.

**truncate**

    This method changes a :py:class:`PQU`'s value and uncertainty to agree with its printed value as defined
    by its significant digits.

**changeUncertaintyStyle, changeUncertaintyPercent**

    These methods allow for the changing of a :py:class:`PQU` instance uncertainty.

------------------------------------------------------------------------------------------
Supporting classes
------------------------------------------------------------------------------------------
In addition to the :py:class:`PQU` class, the following classes are a part the the PQU module:
:py:class:`Parsers`, :py:class:`PQU_float`, :py:class:`PQU_uncertainty`, :py:class:`PhysicalUnit` and :py:class:`NumberDict`.
These classes are only intended as helper classes for the :py:class:`PQU` class; hence, they have limited functionality and
are only intended for internal use.

------------------------------------------------------------------------------------------
History
------------------------------------------------------------------------------------------

This module is based on the Scientific/Physics/PhysicalQuantities.py module (Last revision: 2007-5-25)
written by Konrad Hinsen <hinsen@cnrs-orleans.fr> with contributions from Greg Ward and Berthold Hoellmann.

Most of the additions to Hinsen's version were done by Nidhi R. Patel and Bret R. Beck <beck6@llnl.gov>
with feedback from Caleb M. Mattoon, Neil Summers and David Brown.

The values of physical constants are taken from the 2010 recommended values from
CODATA. Other conversion factors (e.g. for British units) come from various 
sources. The correctness of all entries in the unit table cannot be guaranteed,
so use this at your own risk.

-----------------------------------
Issues or features to be aware of
-----------------------------------

This section describes several oddities about :py:class:`PQU`.

    - Units 'oz', 'lb' and 'ton' (i.e., ounce, pound and ton respectively) are units of mass and not force.
      As example, the unit 'psi' (i.e., pounds per square inch) is equivalent to about '32.174 ft / s**2 * lb / inch**2'.

    - In the following, the first line produces a significantDigits that is small but correct per :py:class:`PQU` 
      rules as the number of significant digits is determined by the '1'.

>>> from pqu import PQU
>>> b1 = PQU.PQU( '1 Bohr' )
>>> b2 = b1.inBaseUnits( )
>>> print(b2)
5.e-11 m
>>> print(b2.value.info( ))
value = 5.29177208114537818e-11, significantDigits = 1, order = -11, isPercent = False
>>> 
>>> b3 = PQU.PQU( 1, 'Bohr' )        # this yields significantDigits = 16.
>>> b4 = b3.inBaseUnits( )
>>> print(b4)
5.291772081145378e-11 m
>>> print(b4.value.info( ))
value = 5.29177208114537818e-11, significantDigits = 16, order = -11, isPercent = False

    - The temperature units are Kelvin (i.e., 'K'), Celsius ('degC'), Rankine ('degR') and Fahrenheit ('degF').
      Two of these units (i.e., 'degC' and 'degF') cannot be used with any other unit including itself. That is,
      PQU( value, 'degC' ) and PQU( value, 'degF' ) are the only allowed forms when 'degF' or 'degC' are present.
      The two absolute temperature units (i.e., 'K' and 'degR') have no such restriction. As an example, heat 
      capacity can be expressed in units of 'W/(m * K)' or 'BTU/(hr * ft * degR)' but cannot be expressed in units
      of 'W/(m * degC)' or 'BTU/(hr * ft * degF)'.
      Any temperature unit can be converted to another temperature unit. For example,

>>> from pqu import PQU
>>> t1 = PQU.PQU( '40 degC' )
>>> print(t1.convertToUnit( 'degF' ))
104. degF
>>> print(t1.getValueAs( 'K' ))
313.15
>>> print(t1.convertToUnit( 'degR' ))
564. degR

    - :py:class:`PQU`, actually :py:class:`PhysicalUnit`, allows the unit to have values. For example 'm/s/25**3' is currently a valid
      unit and stored that way.  This should probably not be allowed and will probably be 
      deprecated (i.e., do not rely on it).
      The following illustrates the issue:

>>> from pqu import PQU
>>> c1 = PQU.PQU( '2.23 c' )
>>> print(c1)
2.23 c
>>> print(c1.getValueAs( 'm/s' ))
668537181.34
>>> print(c1.convertToUnit( 'km / s' ))
6.69e5 km/s
>>> print(c1.getValue( ))
668537.18134
>>> 
>>> c2 = PQU.PQU( '2.23 c' ) / 3
>>> print(c2)
0.743 c
>>> print(c2.getValueAs( 'm/s' ))
222845727.113
>>> print(c2.convertToUnit( 'km / s' ))
2.23e5 km/s
>>> print(c2.getValue( ))
222845.727113
>>> 
>>> c3 = PQU.PQU( '2.23 c / 3' )
>>> print(c3)
2.23 c/3
>>> print(c3.getValueAs( 'm/s' ))
222845727.113
>>> print(c3.convertToUnit( 'km / s' ))
2.23e5 km/s
>>> print(c3.getValue( ))
222845.727113
>>> 
>>> print(c1 == c2, c1 == c3, c2 == c3)
False False True
>>> 
>>> print(PQU.PQU( '2.23 2. * c / 3**2 / pi' ))
2.23 c*2.0/9/3.14159265359

----------------------------
Future plans
----------------------------

    - Use fractions.Fraction in units.
    - Added uncertainty to data in pqu_constants.py and supporting logic.
"""

import sys, math, re, copy
import fractions
from functools import reduce
from .NumberDict import NumberDict

MAJORVERSION = 1
MINORVERSION = 1
PATCHVERSION = 0
__version__ = '%s.%s.%s' % ( MAJORVERSION, MINORVERSION, PATCHVERSION )

maxSignificantDigits = sys.float_info.dig


def compare( value1, value2, epsilonFactor = 0 ) :
    """
    This function compares two floats (or objects that can be converted to floats) in a fuzzy way 
    where the fuzz factor is *epsilonFactor* * sys.float_info.epsilon. This function returns

         0          if the floats are comparable as given by *epsilonFactor* (see below),
         otherwise, it returns 1 (-1) if *value1* is greater (less) than *value2*.

    Two floats are comparable if the magnitude of the 'relative difference' between them is less than or equal to 
    *epsilonFactor* * sys.float_info.epsilon. The relative difference is defined as

         ( value1 - value2 ) / max( abs( value1 ), abs( value2 ) )

    Hence, two floats are comparable if

        ( value1 - value2 ) <= epsilonFactor * sys.float_info.epsilon * max( abs( value1 ), abs( value2 ) )

    For the default *epsilonFactor* = 0, a 0 is return (i.e., the two floats are 'equal') only if they have the same value

    :param value1: value to compare to value2
    :type value1: any object which float() accepts

    :param value2: value to compare to value1
    :type value2: any object which float() accepts

    :param epsilonFactor: The factor to scale sys.float_info.epsilon to determine the fuzz factor.
    :type epsilonFactor: any object which float() accepts

    :returns: 1 if value1 is deemed greater than value2, 0 if equal and -1 otherwise
    :rtype: `int`
    """

    valueOfSelf, valueOfOther = float( value1 ), float( value2 )
    delta = valueOfSelf - valueOfOther
    if( delta == 0 ) : return( 0 )
    _max = max( abs( valueOfSelf ), abs( valueOfOther ) )
    if( abs( delta ) <= ( float( epsilonFactor ) * sys.float_info.epsilon * _max ) ) : return( 0 )
    if( delta < 0. ) : return( -1 )
    return( 1 )

def valueOrPQ( value, unitFrom = None, unitTo = None, asPQU = False, checkOrder = True ) :
    """
    This function is designed as a convenience function for working with :py:class:`PQU`
    and float instances. That is, instead of checking the type of value, let valueOrPQ handle some of the
    details. The first argument can be either a :py:class:`PQU` or something
    that is a valid argument for the python 'float' function (e.g., 1.23, "1.23"). The returned instance
    is a :py:class:`PQU` if asPQU is **True** and a float otherwise. The unitFrom and unitTo
    arguments can be either None, a :py:class:`PhysicalUnit` or a string that represents a :py:class:`PhysicalUnit` (e.g., "m", "",
    "MeV * b"). If unitTo is not None, the returned instance will be in units of unitTo. If unitFrom is not None
    and the first argument (i.e., value) is a :py:class:`PQU` instance, it must have the
    same units.

    Note, this function was originally written when :py:class:`PQU` did not support dimensionless
    units. As it now supports dimensionless units, this function may not be of much value but should be kept anyway.

    Options for input are:
        +-----+--------+------+-----------------+---------------------------------------+
        |     |        |      |  from _getUnit  |           returned value asPQU        |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |value|unitFrom|unitTo|unitFrom|unitTo  |**False**  |**True**                   |
        +=====+========+======+========+========+===========+===========================+
        |float|None    |None  |None    |None    |float      |float                      |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |float|None    |string|None    |PU      |float in   |PQ in unitTo               |
        |     |        |or PU |        |        |unitTo     |                           |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |float|string  |None  |PU      |None    |float in   |PQ in unitFrom             |
        |     |or PU   |      |        |        |unitFrom   |                           |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |float|string  |string|PU      |PU      |float in   |PQ in unitTo               |
        |     |or PU   |or PU |        |        |unitTo     |                           |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |PQU  |None    |None  |PU      |None    |float      |PQ                         |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |PQU  |string  |None  |PU      |None    |float in   |PQU (need to check that PQ |
        |     |or PU   |      |        |        |unitFrom   |and unitFrom are the same) |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |PQU  |string  |string|PU      |PU      |float in   |PQU in unitTo (ditto)      | 
        |     |or PU   |or PU |        |        |unitTo     |                           |
        +-----+--------+------+--------+--------+-----------+---------------------------+
        |PQU  |None    |string|PU      |PU      |float in   |PQU in unitTo              |
        |     |        |or PU |        |        |unitTo     |                           |
        +-----+--------+------+--------+--------+-----------+---------------------------+

    :param value: The number object convert if requested
    :type value: :py:class:`PQU` instance or any object which float() accepts

    :param unitFrom: If not None, the unit for value
    :type unitFrom: None or instance from which units can be determined

    :param unitTo: If not None, the unit the returned value is in
    :type unitTo: None or instance from which units can be determined

    :param asPQU: If **True** returned instance is a :py:class:`PQU`; otherwise, it is a float
    :type asPQU: `bool`

    :returns: value - in units of unitTo if unitTo is not None
    :rtype: `float` or :py:class:`PQU` instance
    """

    def _getUnit( unit, default = None ) :
        """
        Returns an instance of :py:class:`PhysicalUnit` based first on the argument 'unit' and then 'default'. 'unit' 
        can only be a :py:class:`PhysicalUnit`, string (e.g., '', 'MeV', 'b * eV') or None. Default must be either 
        a :py:class:`PQU instance` or None.
        """

        if( isinstance( default, PQU ) ) :
            default = default.unit
        else :
            default = None

        if( unit is None ) :
            if( default is None ) : return( None )
            return _findUnit( default )

        return( _findUnit( unit ) )

    unitFrom = _getUnit( unitFrom, value )
    unitTo = _getUnit( unitTo )

    if( unitFrom is None ) :        # Handles first 2 cases in table above.
        if( asPQU ) : return( PQU( value, unitTo, checkOrder = checkOrder ) )
        return( float( value ) )

    if( isinstance( value, PQU ) ) :
        if( unitFrom != value.unit ) : raise Exception( "Unit '%s' not compatible with unit '%s'" % ( unitFrom, value.unit ) )
        value = copy.deepcopy(value)
    else :
        value = PQU( value, unitFrom, checkOrder = checkOrder )  # At this point unitFrom is a PhysicalUnit instance.
    if( unitTo is not None ) : value.convertToUnit( unitTo )
    if( asPQU ) : return( value )
    return( float( value ) )

def floatToShortestString( value, significantDigits = 15, trimZeros = True, keepPeriod = False, favorEFormBy = 0, includeSign = False ) :
    """
    This function returns the shortest string representation of the float 'value' that is accurate to 
    'significantDigits' digits when converted back to a float.

    The float is converted to both the E-form (i.e., '%e') and F-form (i.e., '%f') with significantDigits.
    Then, after 'trimZeros' and 'keepPeriod' options are implemented, if the length of the E-form minus 
    favorEFormBy is less than or equal to the length of the F-form, the E-form is returned. Otherwise 
    the F-form is returned.

    For example, for significantDigits = 12 the returned string for the float 1.234 will be "1.234", and
    not "1.234000000000" or "1.23400000000e+00" while the string for the float 1.234e-9 will be "1.234-9", and
    not "0.00000000123400000000" or "1.23400000000e-09".

    :param value:               The float to convert to a string.
    :type value: float

    :param significantDigits:   The number of significant digits desired. Restricted to the range 0 to 24 inclusive.
    :type significantDigits: integer

    :param trimZeros:           If **True**, unneeded zeros to the right of '.' are removed.
    :type trimZeros: bool

    :param keepPeriod:          If **False**, '.' is removed if there is no digit to its right.
    :type keepPeriod: bool

    :param favorEFormBy:        The integer subtracted from the length of the E-form before the form
                                with the shortest representation is determined.
    :type favorEFormBy: integer

    :param includeSign:         If **True**, the returned string will always start with a sign character 
                                (i.e., '+' or '-'). Otherwise, only negative values will have a sign.
    :type includeSign: bool

    :rtype: `str`
    """

    if type(value)==str: return value

    def isfinite( value ) :
        """Returns **True** if value is neither infinite nor not-a-number (NaN) and **False** otherwise."""

        if( math.isinf( value ) or math.isnan( value ) ) : return( False )
        return( True )

    signNeeded = ''
    if( includeSign ) : signNeeded = '+'

    if( not( isfinite( value ) ) ) : return( ( "%" + signNeeded + "e" ) % value )

    significantDigitsMinus1 = min( 24, max( 0, significantDigits - 1 ) )  # 24 is arbitrary.

    EForm = ( '%%%s.%de' % ( signNeeded, significantDigitsMinus1 ) ) % value

    mantissa, exponent = EForm.split( 'e' )
    if( significantDigitsMinus1 == 0 ) :
        if( mantissa[-1] != '.' ) : mantissa += '.'     # This is required, else the trimZeros will remove all characters for value = 0.
    if( trimZeros ) : mantissa = mantissa.rstrip( '0' )
    if( not( keepPeriod ) and ( mantissa[-1] == '.' ) ) : mantissa = mantissa[:-1]

    exponent = int( exponent )
    if( exponent == 0 ) : return( mantissa )
    EForm = mantissa + 'e' + "%d" % exponent            # Reconstruct "%e" string.

    digitsRightOfPeriod = significantDigitsMinus1 - exponent  # For "%f" string, determine number of digit to the right of the '.'.
    if( ( digitsRightOfPeriod > 25 ) or ( exponent > 50 ) ) : return( EForm )
    digitsRightOfPeriod = max( digitsRightOfPeriod, 0 )

    FForm = ( '%%%s.%df' % ( signNeeded, digitsRightOfPeriod ) ) % value
    if( 'e' in FForm ) : return( EForm )                # For very large numbers, the '%.?f' format returns 'e' format (e.g., .1234e300 --> "1e+299").

    if( '.' in FForm ) :
        if( trimZeros ) : FForm = FForm.rstrip( '0' )
        if( not( keepPeriod ) ) :
            if( FForm[-1] == '.' ) : FForm = FForm[:-1]
    else :
        if( keepPeriod ) : FForm += '.'
    if( ( len( FForm ) + favorEFormBy ) < len( EForm ) ) : return( FForm )

    return( EForm )

toShortestString = floatToShortestString    # For legacy use. Deprecated: last date 1-April-2014 (no joke).

class PQU_float :
    """
    This class is used by the classes :py:class:`PQU` and :py:class:`PQU_uncertainty` to store more information about a float than
    is provided by the python float class as well as methods to operator on this information. The main members are:

    - value              --- The python value of the float.
    - significantDigits  --- The number of significant digits the value possess as defined by the creator of *self*.
    - order              --- An integer that represents the order of the most significant digit of *self*. 
                             Calculated internally equivalent to int( log10( value ) ) and stored for convenience.
    - _isPercent         --- **True** if *self* represents a percent and **False** otherwise. Note, the value is always stored as
                             the non-percent value (e.g., 1% and 12% are stored as 0.01 and 0.12 respectively).
    """

    def __init__( self, value, significantDigits, isPercent = False, checkOrder = True ) :
        """
        Returns a new :py:class:`PQU_float` instance with value, significantDigits and isPercent. Order is calculated from value.
        If argument checkOrder is **True**, staticmethod *self*.fixFloatsOrder is called and its returned value and order
        are used.
        """

        value = float( value )
        if( isPercent ) : value /= 100.
        self.order = self.calculateOrder( value )
        if( checkOrder ) : value, self.order = self.fixFloatsOrder( value )
        self.value = value
        self._isPercent = bool( isPercent )
        self.setSignificantDigits( significantDigits )

    def __repr__( self ) :

        return( '%s( "%s" )' % ( self.__class__.__name__, self ) )

    def __str__( self ) :

        return( self.toString( ) )

    def __hash__( self ) :

        return( hash( self.value ) + hash( self._isPercent ) + hash( self.significantDigits ) )

    def __eq__( self, other ) :

        return( self.value == float( other ) )

    def __ne__( self, other ) :

        return( self.value != float( other ) )

    def __le__( self, other ) :

        return( self.value <= float( other ) )

    def __lt__( self, other ) :

        return( self.value < float( other ) )

    def __ge__( self, other ) :

        return( self.value >= float( other ) )

    def __gt__( self, other ) :

        return( self.value > float( other ) )

    def __abs__( self ) :

        return( PQU_float( abs( self.value ), self.significantDigits, self._isPercent, checkOrder = False ) )

    def __neg__( self ) :

        return( PQU_float( -self.value, self.significantDigits, self._isPercent, checkOrder = False ) )

    def __add__( self, other ) :

        if( isinstance( other, PQU ) ) : return( other.__radd__( self ) )
        if( self._isPercent ) : raise Exception( 'Addition of per cent is not allowed (defined)' )
        significantOrder = self.order - self.significantDigits
        if( isinstance( other, PQU_float ) ) :
            if( other._isPercent ) : raise Exception( 'Addition of per cent is not allowed (defined)' )
            value = other.value
            significantOrder = max( significantOrder, other.order - other.significantDigits )
        else :
            value = float( other )
        pqu_f = PQU_float( value + self.value, significantOrder, False )
        significantOrder = pqu_f.order - significantOrder
        pqu_f.setSignificantDigits( significantOrder )
        return( pqu_f )

    __radd__ = __add__

    def __iadd__( self, other ) :

        self._setFrom( self + other )
        return( self )

    def __sub__( self, other ) :

        if( isinstance( other, PQU ) ) : return( other.__rsub__( self ) )
        if( isinstance( other, PQU_float ) ) :
            other_ = -other
        else :
            other_ = -float( other )
        return( self.__add__( other_ ) )

    def __rsub__( self, other ) :

        return( -self + other )

    def __isub__( self, other ) :

        self._setFrom( self - other )
        return( self )

    def __mul__( self, other ) :

        if( isinstance( other, PQU ) ) : return( other.__rmul__( self ) )
        significantDigits = self.significantDigits
        if( isinstance( other, PQU_float ) ) :
            value, significantDigits = other.value, min( self.significantDigits, other.significantDigits )
        else :
            value = float( other )
        return( PQU_float( value * self.value, significantDigits, False ) )

    __rmul__ = __mul__

    def __imul__( self, other ) :

        self._setFrom( self * other )
        return( self )

    def __truediv__( self, other ) :

        if( isinstance( other, PQU ) ) : return( other.__rtruediv__( self ) )
        significantDigits = self.significantDigits
        if( isinstance( other, PQU_float ) ) :
            value, significantDigits = other.value, min( self.significantDigits, other.significantDigits )
        else :
            value = float( other )
        return( PQU_float( self.value / value, significantDigits, False ) )

    def __rtruediv__( self, other ) :

        value = float( other )
        return( PQU_float( value / self.value, self.significantDigits ) )

    def __itruediv__( self, other ) :

        self._setFrom( self / other )
        return( self )

    __div__ = __truediv__   # for Python 2.x
    __rdiv__ = __rtruediv__
    __idiv__ = __itruediv__

    def __pow__(self, other):

        power = float(other)
        return PQU_float(self.value**power, self.significantDigits)

    def __rpow__(self, other):

        base = float(other)
        return PQU_float(base**self.value, self.significantDigits)

    def __float__( self ) :

        return( self.value )

    def __bool__( self ) :

        return( self.value != 0. )

    __nonzero__ = __bool__      # for python 2.x

    def _setFrom( self, other ) :
        """Sets *self*'s members from other. Other must be a :py:class:`PQU_float` instance."""

        if( not( isinstance( other, PQU_float ) ) ) : raise TypeError( "other type %s not supported" % type( other ) )
        self.value = other.value
        self._isPercent = other._isPercent
        self.order = other.order
        self.significantDigits = other.significantDigits

    def getValue( self ) :      # Deprecated - replaced with __float__
        "Returns *self*'s value. Equivalent to float( *self* )."

        return( self.value )

    def getOrder( self ) :
        """Returns *self*'s order.

        :rtype: `int`
        """

        return( self.order )

    def getSignificantDigits( self ) :
        """Returns *self*'s significant digits.

        :rtype: `int`
        """

        return( self.significantDigits )

    def info( self, significantDigits = 17 ) :
        """
        Returns a detailed string of *self*. Mainly for debugging.

        :rtype: `str`
        """

        value = "%%.%de" % significantDigits % self.value
        return( "value = %s, significantDigits = %d, order = %d, isPercent = %s" % ( value, self.significantDigits, self.order, self._isPercent ) )

    def isPercent( self ) :
        """Returns **True** if *self* is a percent and **False** otherwise.

        :rtype: `bool`
        """

        return( self._isPercent )

    def setPercentFlag( self, isPercent, scaleValue = False ) :
        """
        If argument isPercent is **True** and *self* is not a percent, *self* is converted to a percent.
        If argument isPercent is **False** and *self* is a percent, *self* is converted to non-percent.
        Otherwise, a raise is executed.

        If scaleValue is **False**, the value is not changed; otherwise, value is scale by 0.01 (100)
        if isPercent is **True** (**False**).

        :param isPercent: see documentation.
        :type isPercent: `bool`

        :param scaleValue: see documentation.
        :type scaleValue: `bool`
        """

        if( bool( isPercent ) == self._isPercent ) : raise Exception( "Percent flag is already %s" % self._isPercent )
        if( isPercent ) :
            self._isPercent, scale, orderDelta = True, 0.01, -2
        else :
            self._isPercent, scale, orderDelta = False, 100., 2
        if( scaleValue ) :
            self.value *= scale
            self.order += orderDelta

    def setSignificantDigits( self, significantDigits ) :
        """
        Sets *self*'s significant digits to significantDigits.

        :param significantDigits: see documentation.
        :type significantDigits: `int`
        """

        self.significantDigits = max( 1, min( maxSignificantDigits + 1, int( significantDigits ) ) )

    def toString( self, trimZeros = True, keepPeriod = True, includePercent = True, favorEFormBy = 0, significantDigits = None ) :
        """
        Returns a string representation of *self* (i.e., '1.234e6' or '12%') using function :py:func:`floatToShortestString`.

        :param bool trimZeros: See function :py:func:`floatToShortestString`
        :param bool keepPeriod: See function :py:func:`floatToShortestString`
        :param int favorEFormBy: See function :py:func:`floatToShortestString`
        :param significantDigits: (*int* or *None*) Used to override self's default significantDigits
        :rtype: `str`
        """

        value = self.value
        if( self._isPercent ) : value *= 100.
        if( significantDigits is None ) : significantDigits = self.significantDigits
        str1 = floatToShortestString( value, significantDigits = significantDigits, trimZeros = trimZeros, 
            keepPeriod = keepPeriod, favorEFormBy = favorEFormBy )
        if( includePercent and self._isPercent ) : str1 += '%'
        return( str1 )

    def truncate( self, significantDigits = None ) :
        """
        Adjusts *self*'s value so that only significantDigits are non-zero. For example, if value is 3.453983207834109
        and significantDigits = 4 then the value is set to 3.454. If significantDigits is None, *self*'s significantDigits
        is used.

        :param significantDigits: 
        :type significantDigits: None or `int`
        """

        if( significantDigits is None ) : significantDigits = self.significantDigits
        significantDigits = min( maxSignificantDigits + 1, max( 1, significantDigits ) )
        self.value = float( "%%.%de" % ( significantDigits - 1 ) % self.value )
        self.setSignificantDigits( significantDigits )

    @staticmethod
    def calculateOrder( value ) :
        """
        Returns the order of value. Value must be convertible to a float. A float's order is the exponent in the
        representation 'm * 10**o' where 'm' is the mantissa and 'o' is the exponent. For example, 1.234e43 has order 43.

        :param value: any object convertible to a float.
        :rtype: integer
        """

        value = abs( float( value ) )
        if( value == 0. ) : return( 0 )
        log10 = math.log10(value)
        if log10.is_integer(): return int(log10)
        s1 = "%.17e" % value
        return( int( s1.split( 'e' )[1] ) )

    @staticmethod
    def fixFloatsOrder( value ) :
        """
        Many numbers are not represented in IEEE floating point by there true value. For example, '1e-12' is 
        stored as the float 9.9999999999999998e-13. The order of the true value is -12 while the order of the floating 
        point representation is -13.  This function is implemented to help with this issue.  For example, for 
        :py:class:`PQU` the string "2.300000000003 (1)" was originally being reprinted as "2.300000000003 (10)" because of this issue.
        To check for this, the value is multiplied by a small (of order sys.float_info.epsilon) factor to see if order 
        changes.  If it does, value is set to the smallest multiplier greater than 1 for which order is changed. 
        The fixed value and its order are returned.

        :param value: any object convertible to a float.
        :rtype: tuple of float and integer
        """

        value = float( value )
        order1 = PQU_float.calculateOrder( value )
        order2 = PQU_float.calculateOrder( value * ( 1 + 4 * sys.float_info.epsilon ) )
        if( order1 != order2 ) :
            for i1 in range( 5 ) :
                value_ = value * ( 1 + i1 * sys.float_info.epsilon )
                order2 = PQU_float.calculateOrder( value_ )
                if( order1 != order2 ) : break
            value = value_
            order1 = order2
        return( value, order1 )

    @staticmethod
    def surmiseSignificantDigits( value ) :
        """
        This factory function returns a :py:class:`PQU_float` instance given a float *value* and sets its significantDigits to
        something "reasonable". The significantDigits are determined by using PQU.floatToShortestString with trimZeros = True.

        :param value:   The float that is converted to a :py:class:`PQU_float` instance.
        """

        value_ = PQU( floatToShortestString( float( value ), trimZeros = True ) )
        return( PQU_float( value, value_.getSignificantDigits( ) ) )


class PQU_uncertainty :
    """
    This class is used by :py:class:`PQU` to store the style of uncertainty and information about its value.
    The members are:

    * style --- One of the three following uncertainty styles supported by :py:class:`PQU`:
                - pqu_uncertaintyStyleNone          - No uncertainty was given for the :py:class:`PQU` instance.
                - pqu_uncertaintyStylePlusMinus     - The uncertainty was given as '+/- value' (e.g., '123 +/- 32').
                - pqu_uncertaintyStyleParenthesis   - The uncertainty was given as '(value)' (e.g., '123(32)'.

    * value --- The value stored as a :py:class:`PQU_float`.
    """

    pqu_uncertaintyStyleNone = None
    pqu_uncertaintyStylePlusMinus = 'plusMinus'
    pqu_uncertaintyStyleParenthesis = 'parenthesis'

    def __init__( self, uncertaintyStyle, value = None, significantDigits = 2, isPercent = False, checkOrder = True ) :
        """
        Constructor method for the :py:class:`PQU_uncertainty` class. The arguments value, significantDigits, isPrecent
        and checkOrder are passed onto :py:class:`PQU_float`. significantDigits will be forced to be 1 or 2.
        """

        if( isPercent and ( uncertaintyStyle != PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) ) :
            raise Exception( 'percent only supported for style "%s"' % PQU_uncertainty.pqu_uncertaintyStylePlusMinus )
        if( value == 0. ) : uncertaintyStyle = PQU_uncertainty.pqu_uncertaintyStyleNone
        significantDigits = max( 1, min( 2, int( significantDigits ) ) )

        self.style = uncertaintyStyle
        if( uncertaintyStyle == PQU_uncertainty.pqu_uncertaintyStyleNone ) :
            self.value = PQU_float( 0., 2, checkOrder = False )
        elif( uncertaintyStyle == PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) :
            self.value = PQU_float( abs( value ), significantDigits, isPercent, checkOrder = checkOrder )
        elif( uncertaintyStyle == PQU_uncertainty.pqu_uncertaintyStyleParenthesis ) :
            self.value = PQU_float( abs( value ), significantDigits, checkOrder = checkOrder )
        else :
            raise TypeError( 'Invalid uncertainty style "%s"' % uncertaintyStyle )

    def __repr__( self ) :

        return( '%s( "%s" )' % ( self.__class__.__name__, self ) )

    def __str__( self ) :

        return( self.toString( ) )

    def __hash__( self ) :

        return( hash( self.style ) + hash( self.value ) )

    def __mul__( self, other ) :
        "For internal use only. Only works if other is a float."

        if( not( isinstance( other, float ) ) ) : TypeError( 'Other type = "%s" not supported' % type( other ) )
        return( PQU_uncertainty( self.style, self.value * other, significantDigits = self.value.significantDigits, isPercent = self.isPercent( ) ) )

    __rmul__ = __mul__

    def __float__( self ) :

        return( float( self.value ) )

    def _changeUncertaintyStyle( self, style ) :
        """For internal use."""

        if( not( self._okayToChangeUncertaintyStyleTo( style ) ) ) : 
            raise Exception( 'Cannot change uncertainty style from "%s" to "%s"' % ( self.style, style ) )
        self.style = style

    def _okayToChangeUncertaintyStyleTo( self, style, checkPercent = True ) :
        """For internal use."""

        if( self.style == style ) : return( True )
        if( self.style == PQU_uncertainty.pqu_uncertaintyStyleNone ) : return( False )
        if( style == PQU_uncertainty.pqu_uncertaintyStyleNone ) : return( False )
        if( ( style != PQU_uncertainty.pqu_uncertaintyStyleParenthesis ) and 
            ( style != PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) ) : raise TypeError( 'Invalid style "%s"' % style )
        if( checkPercent and self.value.isPercent( ) ) : raise Exception( "Cannot change pqu_uncertaintyStylePlusMinus that is a percent" )
        return( True )

    def getStyle( self ) :
        """
        Returns *self*'s style.

        :rtype: One of the three uncertainty styles
        """

        return( self.style )

    def getProperValue( self, v ) :
        """
        Returns *self*'s value. If value is a percent, returned result is multiplied by argument v first, where v
        is expected to be the number associated the the percent value.

        :param v:       Only used if *self*.value is a percent and only the float(v) is used in that case.
        :rtype: `float`
        """

        value = float( self.value )
        if( self.value.isPercent( ) ) : value *= float( v )
        return( value )

    def info( self, significantDigits = 17 ) :
        """Returns a detailed string of *self*. Mainly for debugging.

        :rtype: `str`
        """

        if( self.style == PQU_uncertainty.pqu_uncertaintyStyleNone ) : return( '' )
        return( self.value.info( significantDigits = significantDigits ) )

    def isUncertainty( self ) :
        """Returns **False** is style is pqu_uncertaintyStyleNone and **True** otherwise.

        :rtype: `bool`
        """

        return( self.style != PQU_uncertainty.pqu_uncertaintyStyleNone )

    def isPercent( self ) :
        """Returns **True** if *self* is a percent and **False** otherwise.

        :rtype: `bool`
        """

        if( self.style == PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) : return( self.value.isPercent( ) )
        return( False )

    def toString( self, prefix = ' ', significantDigits = None ) :
        """
        Returns a string representation of *self* (i.e., '+/- 1.2%').

        :param str prefix: A prefix for the '+/-' style
        :param significantDigits: Passed onto method :py:func:`PQU_float.toString`
        :rtype: `str`
        """

        if( self.style == PQU_uncertainty.pqu_uncertaintyStyleNone ) : return( '' )
        if( self.style == PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) :
            favorEFormBy = 0
            if( self.value.order - self.value.getSignificantDigits( ) == 0 ) : favorEFormBy = 1
            return( '%s+/- %s' % ( prefix, self.value.toString( trimZeros = False, favorEFormBy = favorEFormBy, significantDigits = significantDigits ) ) )
        return( '(%d)' % int( pow( 10., self.value.getSignificantDigits( ) - self.value.getOrder( ) - 1 ) * float( self.value ) + 1e-12 ) )

    def truncate( self ) :
        """
        Adjusts *self*'s value so that only the significant digits are non-zero. This is most useful when a :py:class:`PQU` instance is
        calculated from other :py:class:`PQU` instance which can create an uncertainty with many non-zero digits. As example, adding two
        :py:class:`PQU`s with uncertainties 1.3 and 3.2 each with two significant digits, produces the uncertainty sqrt( 1.3**2 + 3.2**2 ) 
        = 3.453983207834109. This uncertainty also only has two significant digits. After calling truncate the uncertainty
        become 3.5 which is the value that the method toString would return.
        """

        self.value.truncate( )

class PQU :
    """
    This class supports a float value with an optional uncertainty and unit. Many basic math operations are supported 
    (e.g., +, -, *, /).  A :py:class:`PQU` is anything that the staticmethod :py:func:`Parsers.parsePQUString` can parse.

    In this section the following notations are used:

        +---------------------------+-----------------------------------------------------------------------------------+
        | Notation                  | Denote                                                                            |
        +===========================+===================================================================================+
        | pqu                       | The argument is a pqu instance.                                                   |
        +---------------------------+-----------------------------------------------------------------------------------+
        | [],                       | The argument is optional.                                                         |
        +---------------------------+-----------------------------------------------------------------------------------+
        | int                       | A python int (e.g., 23).                                                          |
        +---------------------------+-----------------------------------------------------------------------------------+
        | float                     | A python float (e.g., 1.23e-12).                                                  |
        +---------------------------+-----------------------------------------------------------------------------------+
        | number                    | Any of the following, int, float, :py:class:`PQU_float`, a string representing a  |
        |                           | float (i.e.,                                                                      |
        |                           | the python float function must be able to convert the string - e.g., '12.3e3' but |
        |                           | not '12.3 eV') or any object having a __float__ method (hence, a numpy float      |
        |                           | will work).                                                                       |
        +---------------------------+-----------------------------------------------------------------------------------+
        | uncertainty = number_unc  | The uncertainty argument must be a valid number or a :py:class:`PQU_uncertainty`  | 
        |                           | object. If number, uncertainty is of style '+/-'.                                 |
        +---------------------------+-----------------------------------------------------------------------------------+
        | unit = string_pu          | The unit argument must be a valid unit string (e.g., "MeV", "MeV/m" but not       |
        |                           | "1 MeV") or a :py:class:`PhysicalUnit` object.                                    |
        +---------------------------+-----------------------------------------------------------------------------------+
        | string                    | String can be anything that the method :py:func:`Parsers.parsePQUString` is       |
        |                           | capable of parsing (e.g., '1.23', '1.23 m', '1.23%', '1.23(12) m',                |
        |                           | '1.23 +/- 0.12m'). If it contains a unit (uncertainty) then the unit              |
        |                           | (uncertainty) argument must be None.                                              |
        +---------------------------+-----------------------------------------------------------------------------------+
                
        Calling options are:

            - PQU( pqu )                                            # The new PQU gets all data from pqu. Unit and uncertainty arguments must be None.
            - PQU( number, [ unit = string_pu ], [ uncertainty = number_unc ] )
            - PQU( string )
            - PQU( string, unit = string_pu )                       # string must not contain a unit.
            - PQU( string, uncertainty = uncertainty )              # string must not contain an uncertainty.
            - PQU( string, unit = string_pu, uncertainty = number ) # string must not contain a unit or an uncertainty.

        For example, the following are allowed and are the same:
            - PQU( "1.23", "m", "12%" )
            - PQU( "1.23 m", uncertainty = "12%" )
            - PQU( "1.23 +/- 12%", "m" )
            - PQU( "1.23 +/- 12%", "m", None )    # Same as above as None is the default.
            - PQU( "1.23+/-12% m" )

        While the following are not allowed:
            - PQU( "1.23+/-12% m", "m" )                    # Unit already given in value.
            - PQU( "1.23+/-12% m", uncertainty = "12%" )    # Uncertainty already given in value.
            - PQU( "1.23+/-12% m", "m", "12%" )             # Unit and uncertainty already given in value.
    """

    def __init__( self, value, unit = None, uncertainty = None, checkOrder = True ) :

        if( isinstance( value, PQU ) ) :
            if( unit is not None ) : raise TypeError( 'When value is a PQU instance, unit must be None not type "%s"' % type( unit ) )
            if( uncertainty is not None ) : raise TypeError( 'When value is a PQU instance, uncertainty must be None not type "%s"' % type( uncertainty ) )
            self.value       = copy.deepcopy( value.value )
            self.unit        = copy.deepcopy( value.unit )
            self.uncertainty = copy.deepcopy( value.uncertainty )
            return

        if( isinstance( uncertainty, PQU_uncertainty ) ) : uncertainty = copy.deepcopy( uncertainty )
        if( isinstance( unit, PhysicalUnit ) ) : unit = copy.deepcopy( unit )

        if( isinstance( value, str  ) ) :
            value, unit_, uncertainty_ = Parsers.parsePQUString( value )
            if( uncertainty_.style != PQU_uncertainty.pqu_uncertaintyStyleNone ) :
                if( uncertainty is not None ) :
                    raise Exception( 'uncertainty argument = "%s" must be None when value = "%s" is string with uncertainty.' % ( uncertainty, value ) )
                uncertainty = uncertainty_
            if( not( unit_.isDimensionless( ) ) ) :
                if( unit is not None ) :
                    raise Exception( 'unit argument = "%s" must be None when value = "%s" is string with unit.' % ( unit, value ) )
            if( len( unit_.symbols ) > 0 ) : unit = unit_
        else :
            if( isinstance( value, PQU_float ) ) :
                value = copy.deepcopy( value )
            else :
                try :
                    value = PQU_float( value, maxSignificantDigits + 1, checkOrder = checkOrder )
                except :
                    raise TypeError( 'Cannot convert value = "%s" to a float.' % ( value, ) )

        if( unit is None ) : unit = ''
        unit = _findUnit( unit )           # unit can be a PhysicalUnit instance or a string convertible to PhysicalUnit object.

        if( uncertainty is None ) :
            uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleNone )
        else :
            if( not( isinstance( uncertainty, PQU_uncertainty ) ) ) :
                significantDigits, isPercent = 2, False
                uncertainty_ = uncertainty
                if( isinstance( uncertainty, str ) ) :
                    uncertainty_, significantDigitsForZero, str2 = Parsers.parseFloat( uncertainty, uncertainty )
                    if( str2[:1] == '%' ) : isPercent, str2 = True, str2[1:]
                    if( len( str2.strip( ) ) > 0 ) : raise Exception( 'Extra characters for uncertainty "%s"' % uncertainty )
                if( isinstance( uncertainty, PQU_float ) ) :
                    significantDigits, isPercent = uncertainty_.getSignificantDigits( ), uncertainty_.isPercent( )
                try :
                    uncertainty_ = float( uncertainty_ )
                except :
                    raise TypeError( 'Invalid uncertainty = "%s"' % uncertainty )
                if( uncertainty_ == 0. ) :
                    uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleNone )
                else :
                    uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStylePlusMinus, uncertainty_, 
                        significantDigits = significantDigits, isPercent = isPercent )
            if( uncertainty.getStyle( ) == PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) :
                uncertaintyLeastOrder = uncertainty.value.getOrder( ) - uncertainty.value.getSignificantDigits( )
                if( uncertainty.value.isPercent( ) ) :
                    uncertainty_ = uncertainty.value * float( value )
                    uncertaintyLeastOrder = uncertainty_.getOrder( ) - uncertainty_.getSignificantDigits( ) # As uncertainty_ is still %.
                significantDigits = value.getOrder( ) - uncertaintyLeastOrder
                value.setSignificantDigits( significantDigits )
            elif( uncertainty.getStyle( ) == PQU_uncertainty.pqu_uncertaintyStyleParenthesis ) :
                parenthesisPower = PQU_float.fixFloatsOrder( pow( 10., value.order - value.significantDigits + 1 ) )[0]
                uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleParenthesis, 
                    uncertainty.value * parenthesisPower, significantDigits = uncertainty.value.getSignificantDigits( ) )

        if( value.isPercent( ) ) :
            if( uncertainty.isPercent( ) ) : raise Exception( 'Percent for both value and uncertainty is not allowed' )
            if( not( unit.isDimensionless( ) ) ) : raise Exception( 'Percent value and unit is not allowed' )
        self.value = value                  # Private instance guaranteed by logic above.
        self.unit = unit                    # Private instance guaranteed by logic above.
        self.uncertainty = uncertainty      # Private instance guaranteed by logic above.

    def __str__( self ) :

        return( self.toString( ) )

    def __repr__( self ) :

        return( '%s( "%s" )' % ( self.__class__.__name__, self ) )

    def __hash__( self ) :

        return( hash( self.value ) + hash( self.uncertainty ) + hash( self.unit ) )

    def __abs__( self ) :

        return( self.__class__( abs( self.value ), self.unit, self.uncertainty ) )

    def __neg__( self ) :

        return( self.__class__( -self.value, self.unit, self.uncertainty ) )

    def __add__( self, other ) :

        if( self is other ) : return( 2 * self )
        other = self._getOtherAsPQU( other )
        factor, offset = other.unit.conversionTupleTo( self )
        value = self.value + ( other.value + offset ) * factor
        uncertaintySelf, uncertaintyOther = self.getUncertaintyValueAs( ), other.getUncertaintyValueAs( ) * factor
        uncertainty = math.sqrt( uncertaintySelf * uncertaintySelf + uncertaintyOther * uncertaintyOther )
                # Guess at style. If it should be pqu_uncertaintyStyleNone then uncertainty will be 0, which PQU_uncertainty will correct.
        uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStylePlusMinus, uncertainty ) 
        leastSignificantOrder = value.order - value.significantDigits
        uncertainty.value.setSignificantDigits( uncertainty.value.order - leastSignificantOrder )
        return( PQU( value, self.unit, uncertainty ) )

    __radd__ = __add__

    def __iadd__( self, other ) :

        self._setFrom( self + other )
        return( self )

    def __sub__( self, other ) :

        if( self is other ) : return( PQU( 0., self.unit, 0. ) )
        if( isinstance( other, str ) ) : other = PQU( other )
        return( self + -other )

    def __rsub__( self, other ) :

        if( self is other ) : return( PQU( 0., self.unit, 0. ) )
        return( other + -self )

    def __isub__( self, other ) :

        self._setFrom( self - other )
        return( self )

    def __mul__( self, other ) :

        other = self._getOtherAsPQU( other )

        selfValue, otherValue = float( self ), float( other )

        uncertaintySelf, uncertaintyOther = self.getUncertaintyValueAs( ), other.getUncertaintyValueAs( )
        if self is other:       # Correlated uncertainties.
            uncertainty = 2 * uncertaintySelf
        else:                   # Assume uncorrelated uncertainties.
            uncertaintySelf2, uncertaintyOther2 = uncertaintySelf * otherValue, uncertaintyOther * selfValue
            uncertainty = math.sqrt(uncertaintySelf2 * uncertaintySelf2 + uncertaintyOther2 * uncertaintyOther2)

        significantDigits = min(self.uncertainty.value.getSignificantDigits(), other.uncertainty.value.getSignificantDigits())
        uncertainty = PQU_uncertainty(PQU_uncertainty.pqu_uncertaintyStylePlusMinus, uncertainty, significantDigits)

        return PQU(self.value * other.value, unit=self.unit * other.unit, uncertainty=uncertainty)

    __rmul__ = __mul__

    def __imul__( self, other ) :

        self._setFrom( self * other )
        return( self )

    def __truediv__( self, other ) :

        if( ( self is other ) and ( float( self ) != 0 ) ) : return( PQU( 1. ) )
        other = self._getOtherAsPQU( other )
        inv_other = 1. / float( other )
        return( self * PQU( inv_other, 1 / other.unit, inv_other * inv_other * other.getUncertaintyValueAs( ) ) )

    def __rtruediv__( self, other ) :

        other = self._getOtherAsPQU( other )
        return( other / self )

    def __itruediv__( self, other ) :

        self._setFrom( self / other )
        return( self )

    __div__ = __truediv__   # for Python 2.x
    __rdiv__ = __rtruediv__
    __idiv__ = __itruediv__

    def __pow__(self, other):
        """Does not include uncertainty of other in calculation. Other must be an object convertible to a float. If it
        is a :py:class:`PQU` instance, it must be dimensionless."""

        other = self._getOtherAsPQU(other)
        if not other.isDimensionless(): raise TypeError('Power must be dimensionless. It has dimension "%s".' % other.unit)
        power = float(other)

        valueToPower = pow(self.value, power)
        value, uncertainty = float(self), None

        if value != 0:  # This is not the correct answer when value is small compared to uncertainty. Needs work?
            uncertainty = ( power * float( valueToPower ) / value ) * self.uncertainty
            if( uncertainty.getStyle( ) == PQU_uncertainty.pqu_uncertaintyStyleParenthesis ) : 
                uncertainty.style = PQU_uncertainty.pqu_uncertaintyStylePlusMinus
        valueToPower = PQU_float( valueToPower, self.value.getSignificantDigits( ) )
        return( PQU( valueToPower, unit = pow( self.unit, power ), uncertainty = uncertainty ) )

    def __rpow__(self, other):
        """Does not include uncertainty of other in calculation. Other must be an object convertible to a float. If it
        is a :py:class:`PQU` instance, it must be dimensionless."""

        if not self.isDimensionless( ): raise TypeError('Power must be dimensionless. It has dimension "%s".' % self.unit)
        other = self._getOtherAsPQU( other )
        base = float(other)

        newValue = pow(base, self.value)
        value, uncertainty = float(self), None

        if value != 0:  # This is not the correct answer when value is small compared to uncertainty. Needs work?
            uncertainty = float(newValue) * math.log(base) * self.uncertainty
            if uncertainty.getStyle() == PQU_uncertainty.pqu_uncertaintyStyleParenthesis:
                uncertainty.style = PQU_uncertainty.pqu_uncertaintyStylePlusMinus
        newValue = PQU_float(newValue, self.value.getSignificantDigits( ))
        return PQU(newValue, unit = pow(other.unit, self.value), uncertainty = uncertainty)

    def __float__( self ) :

        return( float( self.value ) )

    def __bool__( self ) :

        return( bool( self.value ) )

    __nonzero__ = __bool__      # for python 2.x

    def __eq__( self, other ) :

        return( self.compare( other, 5 ) == 0 )

    def __ne__( self, other ) :

        return( self.compare( other, 5 ) != 0 )

    def __lt__( self, other ) :

        return( self.compare( other, 5 ) < 0 )

    def __le__( self, other ) :

        return( self.compare( other, 5 ) <= 0 )

    def __gt__( self, other ) :

        return( self.compare( other, 5 ) > 0 )

    def __ge__( self, other ) :

        return( self.compare( other, 5 ) >= 0 )

    def _setFrom( self, other ) :
        """
        Sets *self*'s members from other. Other must be an instance of :py:class:`PQU`.

        :param other: The :py:class:`PQU` instance whose members are copied to *self*
        :type other: a :py:class:`PQU` instance
        """

        if( not( isinstance( other, PQU ) ) ) : raise TypeError( "other type %s not supported" % type( other ) )
        self.value = copy.deepcopy(other.value)
        self.unit = copy.deepcopy(other.unit)
        self.uncertainty = copy.deepcopy(other.uncertainty)

    def changeUncertaintyStyle( self, style ) :
        """
        Changes *self*.uncertainty's style. If style is the same as *self*.uncertainty's style nothing is done. If style
        or *self*.uncertainty's style is PQU_uncertainty.pqu_uncertaintyStyleNone as raise is executed.

        :param style: One of the three PQU_uncertainty styles
        """

        self.uncertainty._okayToChangeUncertaintyStyleTo( style, checkPercent = False )
        self.uncertainty._changeUncertaintyStyle( style )

    def changeUncertaintyPercent( self, toPercent ) :
        """
        If *self*'s style is not pqu_uncertaintyStylePlusMinus a raise is executed. Otherwise, if toPercent
        is **True** (**False**) *self*'s uncertainty will be converted to a percent (non-percent) if not already a percent
        (non-percent).

        :param toPercent: Specifies whether *self*'s uncertainty should be percent or not
        :type toPercent: `bool`
        """

        style = self.uncertainty.getStyle( )
        if( style != PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) :
            raise TypeError( 'Can only change pqu_uncertaintyStylePlusMinus style to percent and not %s' % style )
        toPercent = bool( toPercent )
        if( toPercent != self.uncertainty.isPercent( ) ) :
            if( self.isPercent( ) ) : raise Exception( 'Uncertainty cannot be percent when value is a percent' )
            if( toPercent ) :
                value = 100 * float( self.uncertainty ) / float( self )
                self.uncertainty = PQU_uncertainty( style, value, 
                    self.uncertainty.value.getSignificantDigits( ), isPercent = True )
            else :
                self.uncertainty = PQU_uncertainty( style, self.getUncertaintyValueAs( ), self.uncertainty.value.getSignificantDigits( ) )

    def compare( self, other, epsilonFactor = 0 ) :
        """
        Compares *self*'s value and unit to other's value and unit (other's value is first converted to unit of *self*) 
        using the :py:func:`compare` function.

        Also see method :py:func:`PQU.equivalent`.

        :param other: The :py:class:`PQU` instance to compare to *self*
        :type other: A :py:class:`PQU` equivalent instance

        :param epsilonFactor: See the :py:func:`compare` function

        :returns: See the :py:func:`compare` function
        """

        other = self._getOtherAsPQU( other )
        valueOfOther = other.getValueAs( self.unit )
        return( compare( self, valueOfOther, epsilonFactor = epsilonFactor ) )

    def convertToUnit( self, unit ) :
        """
        Changes *self*'s unit to unit and adjusts the value such that the new value/unit is equivalent to
        the original value/unit. The new unit must be compatible with the original unit of *self*.

        :param str unit: a unit equivalent
        :raises TypeError: if the unit string is not a known unit or is a unit incompatible with *self*'s unit
        """

        if( isinstance( unit, PQU ) ) : unit = unit.unit
        factor, offset = self.unit.conversionTupleTo( unit )
        unit = _findUnit( unit )

        self.value = ( self.value + offset ) * factor
        if( not( self.uncertainty.isPercent( ) ) ) : self.uncertainty = self.uncertainty * factor
        self.unit = unit
        return( self )

    def copyToUnit( self, unit ) :
        """
        Like convertToUnit except a copy is made and returned. The copy unit is changed to unit. Self is left unchanged.

        :param str unit: a unit equivalent
        :rtype: :py:class:`PQU`
        """

        return( copy.deepcopy(self).convertToUnit( unit ) )

    def equivalent( self, other, factor = 1. ) :
        """
        This method compares *self* to other and returns **True** if they are equivalent and **False** otherwise. The
        two are considered equivalent if their bands overlap. The band for *self* (or other) is all
        values within 'factor times uncertainty' of its value. That is, the band for *self* is the ranges from
        'value - factor * uncertainty' to 'value + factor * uncertainty'.

        Note, this is different than the :py:func:`compare` method which compares the values using sys.float_info.epsilon
        to calculate a 'band' and not their uncertainties.

        :param other: The object to compare to *self*
        :type other: a :py:class:`PQU` equivalent
        :param factor: the 1/2 width in standard deviations of the band
        :type factor: `float`
        :returns: 1 if *self* is deemed greater than other, 0 if equal and -1 otherwise
        """
# BRB. Should this be named equivalent or consistent.

        factor = abs( factor )
        other = self._getOtherAsPQU( other )
        other2 = other.copyToUnit( self.unit )
        selfValue, selfUncFactor = float( self ), factor * self.getUncertaintyValueAs( )
        otherValue, otherUncFactor = float( other2 ), factor * other2.getUncertaintyValueAs( )

        if( ( selfValue + selfUncFactor ) < ( otherValue - otherUncFactor ) ) : return( False )
        if( ( selfValue - selfUncFactor ) > ( otherValue + otherUncFactor ) ) : return( False )
        return( True )

    def inUnitsOf( self, *units ) :
        """
        Express the quantity in different units. If one unit is specified, a new
        :py:class:`PQU` object is returned that expresses the quantity
        in that unit. If several units are specified, the return value is a tuple of
        PhysicalObject instances with one element per unit such that the sum of all
        quantities in the tuple equals the the original quantity and all the values
        except for the last one are integers. This is used to convert to irregular unit
        systems like hour/minute/second.

        :param units: one or several units
        :type units: `str` or sequence of `str`
        :returns: one or more physical quantities
        :rtype: `PQU` or `tuple` of `PQU`
        :raises TypeError: if any of the specified units are not compatible with the original unit
        """

        units = list( map( _findUnit, units ) )

        if( len( units ) == 1 ) :
            elf = self.__class__( self.value, self.unit, self.uncertainty )
            elf.convertToUnit( units[0] )
            return( elf )
        else :
            units.sort( reverse = True )
            result = []
            value = self.value
            baseUnit = units[-1]
            base = value * self.unit.conversionFactorTo( baseUnit )
            for subunit in units[:-1] :
                value = base * baseUnit.conversionFactorTo( subunit )
                rounded = _round( value )
                result.append( self.__class__( rounded, subunit ) )
                base = base - rounded * subunit.conversionFactorTo( baseUnit )
            result.append( self.__class__( base, units[-1] ) )

            return( tuple( result ) )

    def getSignificantDigits( self ) :
        """
        Returns the number of significant digits for self's value.

        :rtype: `int`
        """

        return( self.value.getSignificantDigits( ) )

    def getUnitAsString( self ) :
        """
        Returns a string representation of self's unit.

        :rtype: `str`
        """

        return( self.unit.symbol( ) )

    getUnitSymbol = getUnitAsString         # To be deprecated.

    def getValueAs( self, unit, asPQU = False ) :
        """
        Returns *self*'s value in units of unit. Unit must be compatible with *self*'s unit.

        :param unit: The unit to return *self*'s value in
        :rtype: `float`
        """

        unit = _findUnit( unit )
        factor, offset = self.unit.conversionTupleTo( unit )
        value = ( float( self.value ) + offset ) * factor
        if( asPQU ) : return( PQU( value, unit ) )
        return( value )

    def getUncertaintyValueAs( self, unit = None, asPQU = False ) :
        """
        Returns a python float representing *self*'s uncertainty value. If the uncertainty is a percent, the returned
        float is *self*'s value times *self*'s uncertainty value.

        :rtype: `float`
        """

        value = PQU( self.uncertainty.getProperValue( self.value ), self.unit )
        if( unit is None ) : unit = self.unit
        return( value.getValueAs( unit, asPQU = asPQU ) )

    def getValue( self ) :      # To be deprecated - replaced by __float__
        """Returns a python float representing *self*'s value. This is equivalent to 'float( self )'."""

        return( float( self ) )

    def inBaseUnits( self ) :
        """
        :returns: A copy of *self* converted to base units, i.e. SI units in most cases
        :rtype: :py:class:`PQU`
        """

        value = ( self.value + self.unit.offset ) * self.unit.factor
        uncertainty = self.uncertainty * self.unit.factor

        num = ''
        denom = ''
        for i1 in range( len( _base_symbols ) ) :
            unit = _base_symbols[i1]
            power = self.unit.powers[i1]

            if( power < 0 ) :
                denom = denom + '/' + unit
                if( power < -1 ) : denom = denom + '**' + str( -power )
            elif( power > 0 ) :
                num = num + '*' + unit
                if( power > 1 ) : num = num + '**' + str( power )
        if( len( num ) == 0 ) :
            num = '1'
            if( len( denom ) == 0 ) : num = ''
        else:
            num = num[1:]

        return( self.__class__( value, num + denom, uncertainty ) )

    def info( self, significantDigits = 17 ) :
        """Returns a detailed string of *self*. Mainly for debugging.

        :rtype: `str`
        """

        return( '%s, unit = "%s"\nuncertainty = %s' % ( self.value.info( significantDigits = significantDigits ), \
            self.unit, self.uncertainty.info( significantDigits = significantDigits ) ) )

    def isCompatible( self, unit ) :
        """
        Returns **True** if *self*'s unit is compatible with unit.

        :param unit: a unit
        :type unit: `str`
        :rtype: `bool`
        """

        unit = _findUnit( unit )
        return( self.unit.isCompatible( unit ) )

    def isPercent( self ) :
        """Returns **True** if value is a percent and **False** otherwise.

        :rtype: `bool`
        """

        return( self.value.isPercent( ) )

    def isDimensionless( self ) :
        """Returns **True** if *self* is dimensionless and **False** otherwise.

        :rtype: `bool`
        """

        return( self.unit.isDimensionless( ) )

    def isLength( self ) :
        """Returns **True** if *self* has unit of length and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'm' ) )

    def isMass( self ) :
        """Returns **True** if *self* has unit of mass and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'g' ) )

    def isTime( self ) :
        """Returns **True** if *self* has unit of time and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 's' ) )

    def isCurrent( self ) :
        """Returns **True** if *self* has unit of current and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'A' ) )

    def isTemperature( self ) :
        """Returns **True** if *self* has unit of temperature and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'K' ) )

    def isMoles( self ) :
        """Returns **True** if *self* has unit of mol and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'mol' ) )

    def isLuminousIntensity( self ) :
        """Returns **True** if *self* has unit of luminous intensity and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'cd' ) )

    def isAngle( self ) :
        """Returns **True** if *self* has unit of angle and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'rad' ) )

    def isSolidAngle( self ) :
        """Returns **True** if *self* has unit of solid angle and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'sr' ) )

    def isEnergy( self ) :
        """Returns **True** if *self* has unit of energy and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'J' ) )

    def isSpin( self ) :
        """Returns **True** if *self* has unit of spin (i.e., same as 'hbar') and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'hbar' ) )

    def isCharge( self ) :
        """
        Returns **True** if *self* has unit of charge and **False** otherwise.

        :rtype: `bool`
        """

        return( self.isPhysicalUnitOf( self, 'e' ) )

    def setPercentFlag( self, isPercent, scaleValue = False ) :
        """
        Changes the percent flag only when *self* is dimensionless and isPercent != *self*.isPercent( ).
        Otherwise, a raise is executed.  For scaleValue argument see :py:func:`PQU_float.setPercentFlag`.

        :param isPercent: See :py:func:`PQU_float.setPercentFlag`
        :param scaleValue: See :py:func:`PQU_float.setPercentFlag`
        """

        if( not( self.isDimensionless( ) ) ) : raise Exception( 'Can only set/unset percent on a dimensionless PQU: unit = "%s"' % self.unit )
        self.value.setPercentFlag( isPercent, scaleValue )
        if( scaleValue and not( self.uncertainty.isPercent( ) ) ) :
            self.uncertainty.value.setPercentFlag( True, isPercent )
            self.uncertainty.value.setPercentFlag( False, not( isPercent ) )

    def simplify( self ) :
        """
        This method returns a :py:class:`PQU` instance that is the same as *self* but with its units simplified. This is,
        units that have the same base units are reduce to one of those units. For example, if *self* is equivalent to
        '1.2 km/m**3' then the returned instance will be equivalent to '1.2e3 / m**2' or '1.2e9 / km**2'
        (the actual units returned are not guaranteed). As another example, if *self* is equivalent to
        '1.2 ft**3 / lb / s / ( m**2 / kg )' then one returned outcome is '0.2457793723470261 ft/s'.
        That is, the 'ft' and 'm' units have the same base units and are all converted to 'ft' (in
        this case) and the 'lb' and 'kg' are all converted to 'lb' (or 'kg').
        """

        return( self * self.unit.simplifyUnitsFactor( ) )

    def toString( self, significantDigits = None, keepPeriod = True ) :
        """
        Returns a string representation of *self* (e.g., '35 +/- 1.2% J').

        :param significantDigits: Passed onto the methods :py:func:`PQU_float.toString` and :py:func:`PQU_uncertainty.toString`
        :param :param : Passed onto the method :py:func:`PQU_float.toString`.

        :rtype: `str`
        """

        trimZeros = self.value.getSignificantDigits( ) > maxSignificantDigits

        includePercent, str1, str2 = True, '', ''
        uncertainty = self.uncertainty
        if( self.isPercent( ) ) :
            if(   self.uncertainty.style == PQU_uncertainty.pqu_uncertaintyStylePlusMinus ) :
                str1, str2, includePercent = '(', ')%', False
                uncertainty = uncertainty * 100
            elif( self.uncertainty.style == PQU_uncertainty.pqu_uncertaintyStyleParenthesis ) :
                str2, includePercent = '%', False

        str1 += self.value.toString( trimZeros = trimZeros, includePercent = includePercent, \
                significantDigits = significantDigits, keepPeriod = keepPeriod ) \
                + uncertainty.toString( significantDigits = significantDigits )
        str3 = str( self.unit )
        if( str3 != '' ) : str1 += ' '
        return( str1 + str2 + str3 )

    def truncate( self, value = False, uncertainty = True ) :
        """
        If value is **True**, calls *self*.value's truncate method. If uncertainty is **True**, 
        calls *self*.uncertainty's truncate method.  See :py:func:`PQU_float.truncate`
        and :py:func:`PQU_uncertainty.truncate`

        :param value: If **True** *self*'s value is truncated
        :type value: `bool`
        :param uncertainty: If **True** *self*'s uncertainty is truncated
        :type uncertainty: `bool`
        """

        if( value ) : self.value.truncate( )
        if( uncertainty ) : self.uncertainty.truncate( )

    def sqrt( self ) :
        """Returns a :py:class:`PQU` of the square root of *self*."""

        return( pow( self, 0.5 ) )

    def sin( self ) :
        """Returns a :py:class:`PQU` of the sin of *self*. Self must have unit of angle."""

        if( self.unit.isAngle( ) ) : return( math.sin( self.value * self.unit.conversionFactorTo( _unit_table['rad'] ) ) )
        raise TypeError( 'Argument of sin must be an angle.' )

    def cos( self ) :
        """Returns a :py:class:`PQU` of the cos of *self*. Self must have unit of angle."""

        if( self.unit.isAngle( ) ) : return( math.cos( self.value * self.unit.conversionFactorTo( _unit_table['rad'] ) ) )
        raise TypeError( 'Argument of cos must be an angle.' )

    def tan( self ) :
        """Returns a :py:class:`PQU` of the tan of *self*. Self must have unit of angle."""

        if( self.unit.isAngle( ) ) : return( math.tan( self.value * self.unit.conversionFactorTo( _unit_table['rad'] ) ) )
        raise TypeError( 'Argument of tan must be an angle.' )

    @staticmethod
    def _getOtherAsPQU( other ) :
        """
        Returns a :py:class:`PQU` representation of other. If other is a :py:class:`PQU`
        instance, it is returned (i.e., no copy is made). Otherwise, the PQU.__init__ 
        method is called with only other as an argument and its return value returned.

        :param other: a :py:class:`PQU` equivalent to convert, if needed, to a :py:class:`PQU` instance
        :returns: a :py:class:`PQU` instance equivalent to other
        """

        if( isinstance( other, PQU ) ) : return( other )
        return( PQU( other ) )

    @staticmethod
    def isPhysicalUnitOf( elf, unit ) :
        """
        Returns **True** if elf has the same unit as the argument unit and **False** otherwise. Elf and unit can be an instance of string,
        :py:class:`PhysicalUnit` or :py:class:`PQU`.

        :param elf: The instance to compare to unit
        :type elf: a unit equivalent
        :param unit: The desired unit
        :type unit: a unit equivalent
        :rtype: `bool`
        """

        elfUnit = _getPhysicalUnit( elf )
        unitUnit = _getPhysicalUnit( unit )
        return( elfUnit.isCompatible( unitUnit ) )

PhysicalQuantityWithUncertainty = PQU           # This is deprecated.

class Parsers :
    """
    This is use to parse the allowed **PQU** strings (e.g., '1 m', '1.32(3) kg).
    See the table in the :py:func:`parsePQUString.parsePQUString` doc string for allowed **PQU** strings.
    This class is mainly for internal use.

    .. rubric:: Regular Expression Matching
    """

    _spaces = r'[ \t]*'
    _anything_RE = '(.*)'
    _sign = '[+-]?'                                     # Optional sign regular expression.
    _unsignedInteger = r'\d+'
    _signedInteger = _sign + _unsignedInteger
    _mantissa = r'(%s)((\d+)\.?(\d*)|\.(\d+))' % _sign   # Regular expression matching a decimal or an integer string.
    _exponent = r'([eE](%s\d+))?' % _sign                # Optional exponent regular expression.

    _floatingPoint_RE = '(%s%s)' % ( _mantissa, _exponent )   # Floating point regular expression.
    _floatingPoint_PO = re.compile( '^' + _spaces + _floatingPoint_RE + _spaces + '$' )

    _floatingPoint_andAnythingElse_RE = '^' + _floatingPoint_RE + _anything_RE
    _floatingPoint_andAnythingElse_PO = re.compile( _floatingPoint_andAnythingElse_RE )

    _uncertaintyPlusMinus_andAnythingElse_RE = '^' + _spaces + r'\+/-' + _spaces + _floatingPoint_RE + _anything_RE
    _uncertaintyPlusMinus_andAnythingElse_PO = re.compile( _uncertaintyPlusMinus_andAnythingElse_RE )

    _uncertaintyParenthesis_andAnythingElse_RE = r'^[(](\d\d?)[)]' + _anything_RE
    _uncertaintyParenthesis_andAnythingElse_PO = re.compile( _uncertaintyParenthesis_andAnythingElse_RE )

    _unitCharacters_RE = '([ a-zA-Z*/()]*)'           # This does not match a valid unit, it only match the part of a
                                                    # string that contains characters that can be present in a unit string.

    @staticmethod
    def parseFloat( str1, fullStr ) :
        """
        This method parses the beginning of str1 for a valid float. An arbitrary number of white spaces can exists before the 
        float characters. This method returns a tuple of length 3. The first item of the tuple is the float 
        value returned as a :py:class:`PQU_float` instance with significantDigits determined from the string.  The second 
        item is significantDigitsForZero which, if the string represents a 0 value, gives the number of significant 
        digits for the zero; otherwise None is returned. As example, the string '000.000' has 4 significant digits 
        (i.e., it is considered equivalent to the representation '0.000').  The last object is the string of all 
        characters after the last one used for converting the float.  For example, with str1 = '   12.345e3ABCD XYS' 
        the returned tuple is:

        ( PQU_float( 1.2345e4, 5 ), None, 'ABCD XYS' ).

        :param str1: String to parse.
        :param fullStr: String which str1 is a part of.
        :rtype: ( PQU_float, integer | None, string )
        """

        """
        The following table list the components of the groups returned by the regular expression match.

        Group |
        index | contents
        ------+-------------------------------------------------------------------------------------
          0   | Float (e.g., for '12.34e-12 +/- 6.7e-13 m / s**2' this would be '12.34e-12')
          1   | Float's sign
          2   | Mantissa part of float
          3   | Float's digits to the left of the period (e.g., for '12.34e-12' this would be '12')
          4   | Float's digits to the right of the period if there are digits to left of period (e.g., for '12.34e-12' this would be '34')
          5   | Float's digits to the right of the period if there are no digits to left of period (e.g., for '.56e-12' this would be '56')
          6   | Float's exponent (e.g., for '12.34e-12' this would be 'e-12')
          7   | Float's exponent value (e.g., for '12.34e-12' this would be '-12')
          8   | Everything after the float
        """

        if( not( isinstance( str1, str ) ) ) : raise TypeError( 'Argument must be a string: type = %s' % type( str1 ) )
        match = Parsers._floatingPoint_andAnythingElse_PO.match( str1 )
        if( match is None ) : raise TypeError( 'String does not start with a float: "%s" of "%s"' % ( str1, fullStr ) )
        groups = match.groups( )
        value, order, significantDigits, significantDigitsForZero = Parsers._parseFloatGroups( groups )
        return( PQU_float( value, significantDigits, False, checkOrder = False ), significantDigitsForZero, groups[8] )

    @staticmethod
    def parsePlusMinusUncertainty( str1, fullStr ) :
        """
        This method parses the beginning of str1 for the sub-string '+/-' followed by a float string. An arbitrary number of 
        white spaces can exists before the '+/-' sub-string and between it and the float string. This method
        returns a tuple of length 2. The first item of the tuple is float value returned as a :py:class:`PQU_uncertainty`
        of style pqu_uncertaintyStylePlusMinus. Characters after the last one used for converting the float are 
        returned as the second item. For example, with str1 = ' +/-  12. ABCD XYS' the returned tuple is:

        ( PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStylePlusMinus, 12. ), ' ABCD XYS' )

        :param str1: String to parse.
        :param fullStr: String which str1 is a part of.
        :rtype: ( PQU_uncertainty, string )
        """

        if( not( isinstance( str1, str ) ) ) : raise TypeError( 'Argument must be a string: type = %s' % type( str1 ) )
        match = Parsers._uncertaintyPlusMinus_andAnythingElse_PO.match( str1 )
        if( match is None ) : raise TypeError( 'String does not match "+/-" with a float: "%s" of "%s"' % ( str1, fullStr ) )
        groups = match.groups( )
        value, order, significantDigits, significantDigitsForZero = Parsers._parseFloatGroups( groups )
        uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStylePlusMinus, value, significantDigits, checkOrder = False )
        return( uncertainty, groups[8] )

    @staticmethod
    def parseParenthesisUncertainty( str1, fullStr ) :
        """
        This method parses the beginning of str1 for the sub-string '(#)' where '#' is a 1 or 2 digit
        number.  This method returns a tuple of length 2. The first item of the tuple is the number returned as 
        a :py:class:`PQU_uncertainty` of style pqu_uncertaintyStyleParenthesis. Characters after the ')' are returned 
        as the second item. For example, with str1 = '(34) ABCD XYS' the returned tuple is:

        ( PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleParenthesis, 34, significantDigits = 2 ), ' ABCD XYS' ).

        :param str1: String to parse.
        :param fullStr: String which str1 is a part of.
        :rtype: ( PQU_uncertainty, string )
        """

        if( not( isinstance( str1, str ) ) ) : raise TypeError( 'Argument must be a string: type = %s' % type( str1 ) )
        match = Parsers._uncertaintyParenthesis_andAnythingElse_PO.match( str1 )
        if( match is None ) : raise TypeError( 'String does not match "(#)" or "(##)" where "#" is a digit: "%s" or "%s"' % 
            ( str1, fullStr ) )
        groups = match.groups( )
        uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleParenthesis, float( groups[0] ), len( groups[0] ), checkOrder = False )
        return( uncertainty, groups[1] )

    @staticmethod
    def _parseFloatGroups( groups ) :
        """
        For internal use. This method analyzes the groups to determine the properties of a float and returns its
        value, order, significantDigits and significantDigitsForZero.
        """

        significantDigitsForZero = None
        if( groups[0] == '' ) : return( None, None, None, None )
        if( float( groups[2] ) == 0. ) :    # Mantissa is zero.
            significantDigits, index = 1, 4
            if( groups[3] is None ) : index = 5
            significantDigitsForZero = len( groups[index] )     # Needed for zero with () style uncertainty (e.g., "0.0000(4)").
        else :
            mantissa = groups[5]
            if( mantissa is None ) :       # Mantissa has digits to the left and maybe right of the period (e.g., "123.45").
                mantissa = groups[3]
                if( groups[4] is not None ) : mantissa += groups[4]
            significantDigits = len( mantissa.lstrip( '0' ) )

        value, order = PQU_float.fixFloatsOrder( float( groups[0] ) )
        return( value, order, significantDigits, significantDigitsForZero )

    @staticmethod
    def parsePQUString( str1 ) :
        """
        Parses the string str1 and returns the tuple ( value, unit, uncertainty ) where value is an instance of 
        :py:class:`PQU_float`, unit is an instance of :py:class:`PhysicalUnit` and uncertainty is an instance 
        of :py:class:`PQU_uncertainty`.

        The string can be one of the following :py:class:`PQU` forms. Here, F represents a valid float string, () represents
        the string '(#)' where '#' is a 1 or 2 digit number, +/- represents the string '+/- F', % is itself
        and u is a unit.

        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  # | Form     | Example        | Value   | Value - uncertainty | Value + uncertainty | Unit |
        +====+==========+================+=========+=====================+=====================+======+
        |  0 | F        | '234'          | 234.    | N/A                 | N/A                 | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  1 | F%       | '234%'         | 2.34    | N/A                 | N/A                 | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  2 | F u      | '234 m'        | 234.    | N/A                 | N/A                 | 'm'  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  3 | F()      | '234(5)'       | 234.    | 229.                | 239.                | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  4 | F()%     | '234(5)%'      | 2.34    | 2.29                | 2.39                | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  5 | F() u    | '234(5) m'     | 234.    | 229.                | 239.                | 'm'  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  6 | F +/-    | '234 +/- 5'    | 234.    | 229.                | 239.                | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  7 | F +/-%   | '234 +/- 5%'   | 234.    | 222.                | 246.                | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  8 | F +/- u  | '234 +/- 5 m'  | 234.    | 229.                | 239.                | 'm'  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        |  9 | F +/-% u | '234 +/- 5% m' | 234.    | 222.                | 246.                | 'm'  |
        +----+----------+----------------+---------+---------------------+---------------------+------+
        | 10 | (F +/-)% | '(234 +/- 5)%' | 2.34    | 2.29                | 2.39                | N/A  |
        +----+----------+----------------+---------+---------------------+---------------------+------+

        Note 1) 5% of 234 is 11.7 rounded to 12.

        :param str1: String to parse.
        :rtype: ( :py:class:`PQU_float`, :py:class:`PhysicalUnit`, :py:class:`PQU_uncertainty` )
        """

        str2, form10, isPercent = str1.strip( ), False, False
        if( str2[:1] == "(" ) : form10, str2 = True, str2[1:]
        value, significantDigitsForZero, str2 = Parsers.parseFloat( str2, str1 )

        if( form10 ) :
            uncertainty, str2 = Parsers.parsePlusMinusUncertainty( str2, str1 )
            uncertainty.value *= 0.01
            str2 = str2.strip( )
            if( str2 != ')%' ) : raise TypeError( 'input not of the form "(F +/- F)%' )
            isPercent, str2 = True, ''
        else :
            uncertainty = PQU_uncertainty( PQU_uncertainty.pqu_uncertaintyStyleNone )
            if( str2[:1] == '%' ) :          # The '%' must immediately follow the float string.
                isPercent, str2 = True, str2[1:]
            else :
                str2 = str2.strip( )
                if( str2[:1] == '(' ) :      # Must be '()' style uncertainty.
                    uncertainty, str2 = Parsers.parseParenthesisUncertainty( str2, str1 )
                    if( str2[:1] == '%' ) : isPercent, str2 = True, str2[1:]
                    if( significantDigitsForZero is not None ) : value.setSignificantDigits( significantDigitsForZero + 1 )
                elif( str2[:1] == '+' ) :    # Must be '+/-' style uncertainty.
                    uncertainty, str2 = Parsers.parsePlusMinusUncertainty( str2, str1 )
                    if( str2[:1] == '%' ) :          # The '%' must immediately follow the float string.
                        str2 = str2[1:]
                        uncertainty.value.setPercentFlag( True, scaleValue = True )

        str2 = str2.strip( )
        if( isPercent ) :                   # str2 must now be empty as units are not allowed with "%".
            if( len( str2 ) > 0 ) : raise TypeError( 'value with percent cannot have units: "%s"' % str1 )
            value.setPercentFlag( True, scaleValue = True )
        try :
            unit = _findUnit( str2 )
        except :
            raise TypeError( 'Invalid unit "%s", in PQU string "%s"' % ( str2, str1 ) )
        return( value, unit, uncertainty )

class PhysicalUnit :
    """
    .. rubric:: PHYSICAL UNIT

    A physical unit is defined by a symbol (possibly composite), a scaling factor, and
    the powers of each of the SI base units that enter into it. Units can be 
    multiplied, divided, and raised to integer powers. This class is mainly for internal use.

    :param symbols: a dictionary mapping each symbol component to its associated integer
                power (e.g., ``{'m': 1, 's': -1}``) for `m/s`). As a shorthand, a string may be 
                passed which is assigned an implicit power 1.
    :type symbols: `dict` or `str`

    :param factor: a scaling factor
    :type factor: `float`

    :param offset: an additive offset to the base unit (used only for temperatures)
    :type offset: `float`

    :param powers: the integer powers for each of the nine base units
    :type powers: `list` of `int`
    """
    
    def __init__( self, symbols, factor, powers, offset = 0 ) :
        '''
        :param symbols: a dictionary mapping each symbol component to its associated integer
                    power (e.g., ``{'m': 1, 's': -1}``) for `m/s`). As a shorthand, a string may be 
                    passed which is assigned an implicit power 1.
        :type symbols: `dict` or `str`

        :param factor: a scaling factor
        :type factor: `float`

        :param offset: an additive offset to the base unit (used only for temperatures)
        :type offset: `float`

        :param powers: the integer powers for each of the nine base units
        :type powers: `list` of `int`
        '''

        if( symbols is not None ) :
            if( '1' in symbols ) :
                rmOne = False
                for key in symbols:
                    if( symbols[key] > 0 ) : rmOne = True
                if( rmOne ) : del symbols['1']

        if( isinstance( symbols, str ) ) :
            self.symbols = NumberDict( )
            self.symbols[symbols] = 1
        else :
            if( symbols is None ) : symbols = {}
            self.symbols = symbols
        if( '' in self.symbols ) : del self.symbols['']

        self.factor = factor
        self.offset = offset
        self.powers = powers

    def __repr__( self ) :

        return( self.__class__.__name__ + "('" + self.symbol( ) + "')" )

    def __hash__( self ) :

        return( hash( self.factor ) + hash( self.offset ) + hash( tuple( self.powers ) ) )

    def __lt__( self, other ) :

        if( self.powers != other.powers ) : raise TypeError( 'Incompatible units' )
        return( self.factor.__lt__( other.factor ) )

    def __eq__( self, other ) :

        if( self.powers != other.powers ) : raise TypeError( 'Incompatible units' )
        return( self.factor.__eq__( other.factor ) )

    def __ne__( self, other ) :

        if( self.powers != other.powers ) : raise TypeError( 'Incompatible units' )
        return( self.factor.__ne__( other.factor ) )

    def __mul__(self, other):

        if self.offset != 0 or (isinstance(other, PhysicalUnit) and other.offset != 0):
            raise TypeError('Cannot multiply units with non-zero offset.')

        if isinstance(other, PhysicalUnit):
            pqu = PhysicalUnit(self.symbols + other.symbols, self.factor * other.factor, list(map(lambda a, b: a + b, self.powers, other.powers)))
        else:
            pqu = PhysicalUnit(self.symbols + {str(other): 1}, self.factor * other, self.powers)
        pqu.simplifyFractions()

        return pqu

    __rmul__ = __mul__

    def __truediv__(self, other):

        if self.offset != 0 or (isinstance(other, PhysicalUnit) and other.offset != 0):
            raise TypeError('Cannot divide units with non-zero offset.')

        if isinstance(other, PhysicalUnit):
            pqu = PhysicalUnit(self.symbols - other.symbols, self.factor / other.factor, list(map(lambda a, b: a - b, self.powers, other.powers)))
        else:
            pqu = PhysicalUnit(self.symbols + {str(other): -1}, self.factor / other, self.powers)
        pqu.simplifyFractions()

        return pqu

    def __rtruediv__(self, other):

        if self.offset != 0 or (isinstance(other, PhysicalUnit) and other.offset != 0):
            raise TypeError('Cannot divide units with non-zero offset.')

        if isinstance(other, PhysicalUnit):
            pqu = PhysicalUnit(other.symbols - self.symbols, other.factor / self.factor, list(map(lambda a, b: a - b, other.powers, self.powers)))
        else:
            pqu = PhysicalUnit(NumberDict({str(other): 1}) - self.symbols, other / self.factor, [-x for x in self.powers])
        pqu.simplifyFractions()

        return pqu

    __div__ = __truediv__   # for python 2.x
    __rdiv__ = __rtruediv__

    def __pow__(self, other):

        if self.offset != 0:
            raise TypeError('Cannot exponentiate units with non-zero offset.')

        if self.isDimensionless():
            return self

        power = fractions.Fraction(other).limit_denominator(10000)
        if power.denominator == 1:
            power = power.numerator
        floatOther = float(other)
        if floatOther != 0.0:
            if abs(float(power) / floatOther - 1) > 1e-14:
                raise ValueError('Could not convert "%s" to a fraction (%s)' % (other, power))
        powers = list(map(lambda x, p = power: x * p, self.powers))
        for index, p1 in enumerate(powers):
            if isinstance(p1, fractions.Fraction):
                if p1.denominator == 1:
                    powers[index] = p1.numerator

        pqu = PhysicalUnit(power * self.symbols, pow(self.factor, power), powers)
        pqu.simplifyFractions()

        return pqu

    def conversionFactorTo( self, other, ignoreTemperatureOffsets = False ) :
        """
        Returns a float factor that can be used to scale a value in unit of *self* to unit of *other*.

        :param other: another unit
        :type other: :py:class:`PhysicalUnit`
        :returns: the conversion factor from this unit to another unit
        :rtype: `float`
        :raises TypeError: if the units are not compatible
        """

        other = _getPhysicalUnit( other )
        if( self.powers != other.powers ) : raise TypeError( 'Incompatible units: cannot convert "%s" to "%s"' % ( str( self ), str( other ) ) )

        if( ( self.offset != other.offset ) and ( self.factor != other.factor ) ) :
            if ignoreTemperatureOffsets or self.powers != _unit_table['K'].powers:
                raise TypeError( 'Unit conversion (%s to %s) cannot be expressed as a simple multiplicative factor' % ( self.symbol( ), other.symbol( ) ) )

        if math.log10( self.factor ).is_integer() and math.log10( other.factor ).is_integer():
            # special treatment to reduce numeric error in order-of-magnitude conversions
            return 10**( math.log10( self.factor ) - math.log10( other.factor ) )

        return( self.factor / other.factor )

    def conversionTupleTo( self, other ) :
        """
        Returns the tuple (factor, offset) that can be used to convert a value in unit of *self* to unit of *other*.
        Offset is 0.0 except with convert some temperature units.

        :param other: another unit
        :type other: :py:class:`PhysicalUnit`
        :returns: the conversion factor and offset from this unit to another unit
        :rtype: (`float`, `float`)
        :raises TypeError: if the units are not compatible
        """

        other_ = _getPhysicalUnit( other )
        if( self.powers != other_.powers ) : raise TypeError( 'Unit "%s" not convertible with "%s"' % ( self, other ) )

        # Let (s1,d1) be the conversion tuple from 'self' to base units (i.e. (x+d1)*s1 converts a value x from 'self' to base units, 
        # and (x/s1)-d1 converts x from base to 'self' units) and (s2,d2) be the conversion tuple from 'other' to base units then 
        # we want to compute the conversion tuple (S,D) from 'self' to 'other' such that (x+D)*S converts x from 'self' units to 
        # 'other' units the formula to convert x from 'self' to 'other' units via the base units is (by definition of the conversion tuples):
        #    ( ((x+d1)*s1) / s2 ) - d2
        #  = ( (x+d1) * s1/s2) - d2
        #  = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
        #  = ( (x+d1) - (d1*s2/s1) ) * s1/s2
        #  = (x + d1 - d2*s2/s1) * s1/s2
        # thus, D = d1 - d2*s2/s1 and S = s1/s2

        factor = self.conversionFactorTo( other_ )
        offset = self.offset - ( other_.offset * other_.factor / self.factor )
        return( factor, offset )

    def isCompatible( self, other ) :
        """
        Returns **True** if *self* is compatible with other and **False** otherwise.

        :param other: another unit
        :type other: :py:class:`PhysicalUnit`
        :returns: **True** if the units are compatible, i.e. if the powers of the base units are the same
        :rtype: `bool`
        """
        return( self.powers == other.powers )

    def isDimensionless( self ) :
        """Returns **True** if *self* is dimensionless and **False** otherwise.

        :rtype: `bool`
        """

        return( not reduce( lambda a, b: a or b, self.powers ) )

    def isAngle( self ) :
        """Returns **True** if *self* is an angle and **False** otherwise.

        :rtype: `bool`
        """

        return self.powers[7] == 1 and reduce( lambda a, b: a + b, self.powers ) == 1

    def setSymbol( self, symbol ) :
        """
        Sets *self*'s symbols to symbol with power 1.

        :param symbol: 
        :type symbol: `str`
        """

        self.symbols = NumberDict( )
        self.symbols[symbol] = 1

    def simplifyFractions(self):
        '''
        This method converts any part of self's symbols and powers that are a :py:class:`fractions.Fraction` instance with 
        denominator of 1 into an integer of the numerator.
        '''

        for index, power in enumerate(self.powers):
            if isinstance(power, fractions.Fraction):
                if power.denominator == 1:
                    self.powers[index] = power.numerator
        
        symbols = self.symbols
        toDelete = []
        for symbol in symbols:
            value = self.symbols[symbol]
            if isinstance(value, fractions.Fraction):
                if value.denominator == 1:
                    self.symbols[symbol] = value.numerator
                    if value.numerator == 0:
                        toDelete.append(symbol)
        for symbol in toDelete:
            del self.symbols[symbol]

    def simplifyUnitsFactor( self ) :
        """
        This method returns a :py:class:`PQU` instance that can be used to simplify *self*'s units.
        That is, all units with the same base units are reduced to a single unit. For example,
        if *self* has units of 'km/m' then the returned :py:class:`PQU` is equivalent to '1e3 m/km'.
        """

        ignores = []
        result = PQU( 1 )
        keys = list( self.symbols.keys( ) )
        for i1, key in enumerate( keys ) :
            if( i1 in ignores ) : continue
            ignores.append( i1 )
            u1 = PQU( 1, key )
            for i2 in range( i1 + 1, len( keys ) ) :
                if( i2 in ignores ) : continue
                if( u1.isCompatible( keys[i2] ) ) :
                    ignores.append( i2 )
                    power2 = self.symbols[keys[i2]]
                    factor = PQU( 1, keys[i2] ).getValueAs( key )**power2
                    result *= factor * PQU( 1, key )**power2 / PQU( 1, keys[i2] )**power2
        return( result )

    def symbol( self ) :
        """
        Returns a string representation of *self*.

        :returns: a string representation of *self*.
        :rtype: `str`
        """

        num = ''
        denom = ''

        for unit in self.symbols :
            power = self.symbols[unit]

            leftBracket, rightBracket = '', ''
            if isinstance(power, fractions.Fraction):
                leftBracket, rightBracket = '(', ')'

            if( power < 0 ) :
                denom = denom + '/' + unit
                if power != -1:
                    denom = denom + '**' + leftBracket + str(-power) + rightBracket
            elif( power > 0 ) :
                num = num + '*' + unit
                if power != 1:
                    num = num + '**' + leftBracket + str(power) + rightBracket

        if( len( num ) == 0 ) :
            if( len( denom ) > 0 ) : num = '1'
        else:
            num = num[1:]

        return num + denom

    def __str__( self ) :

        return( self.symbol( ) )
#
# Helper functions
#
def _findUnit( unit ) :

    if( isinstance( unit, str ) ) :
        symbol = unit.strip( )

        if( symbol == '' or symbol == 'as') : # exception to handle atto-second: conflict with python keyword "as"
            unit = _unit_table[symbol]
        else:
            unit = eval( symbol, _unit_table )

        for cruft in ['__builtins__', '__args__'] :
            try:
                del _unit_table[cruft]
            except:
                pass

    if( not( isinstance( unit, PhysicalUnit ) ) ) : raise TypeError( '%s is not a unit' % type( unit ) )
    return( unit )

def _round( x ) :

    if( x > 0. ) :
        x_ = math.floor( x )
    else:
        x_ = math.ceil( x )
    if( isinstance( x, PQU_float ) ) : x_ = PQU_float( x_, x.getSignificantDigits( ) )
    return( x_ )

def _getPhysicalUnit( elf ) :

    if( isinstance( elf, PQU ) ) : return( elf.unit )
    try :
        return( _findUnit( elf ) )
    except :
        return( PQU( elf ).unit )

# SI unit definitions

_base_symbols = [ 'm', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr' ]

_base_units = [ ( 'm',   PhysicalUnit( 'm',   1.,    [ 1, 0, 0, 0, 0, 0, 0, 0, 0 ] ) ),
                ( 'g',   PhysicalUnit( 'g',   0.001, [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ] ) ),
                ( 's',   PhysicalUnit( 's',   1.,    [ 0, 0, 1, 0, 0, 0, 0, 0, 0 ] ) ),
                ( 'A',   PhysicalUnit( 'A',   1.,    [ 0, 0, 0, 1, 0, 0, 0, 0, 0 ] ) ),
                ( 'K',   PhysicalUnit( 'K',   1.,    [ 0, 0, 0, 0, 1, 0, 0, 0, 0 ] ) ),
                ( 'mol', PhysicalUnit( 'mol', 1.,    [ 0, 0, 0, 0, 0, 1, 0, 0, 0 ] ) ),
                ( 'cd',  PhysicalUnit( 'cd',  1.,    [ 0, 0, 0, 0, 0, 0, 1, 0, 0 ] ) ),
                ( 'rad', PhysicalUnit( 'rad', 1.,    [ 0, 0, 0, 0, 0, 0, 0, 1, 0 ] ) ),
                ( 'sr',  PhysicalUnit( 'sr',  1.,    [ 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ) ),
              ]

_prefixes = [( 'Y',  1.e24, 'yotta' ),
             ( 'Z',  1.e21, 'zetta' ),
             ( 'E',  1.e18, 'exa' ),
             ( 'P',  1.e15, 'peta' ),
             ( 'T',  1.e12, 'tera' ),
             ( 'G',  1.e9, 'giga' ),
             ( 'M',  1.e6, 'mega' ),
             ( 'k',  1.e3, 'kilo' ),
             ( 'h',  1.e2, 'hecto' ),
             ( 'da', 1.e1, 'deka' ),
             ( 'd',  1.e-1, 'deci' ),
             ( 'c',  1.e-2, 'centi' ),
             ( 'm',  1.e-3, 'milli' ),
             ( 'mu', 1.e-6, 'micro' ),
             ( 'n',  1.e-9, 'nano' ),
             ( 'p',  1.e-12, 'pico' ),
             ( 'f',  1.e-15, 'femto' ),
             ( 'a',  1.e-18, 'atto' ), # We support atto-second, but no operations within a single PhysicalUnit instance.
             ( 'z',  1.e-21, 'zepto' ),
             ( 'y',  1.e-24, 'yocto' ),
             ]

_help = []
_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]

def _addUnit( symbol, unit, name = '' ) :

    if( symbol in _unit_table ) : raise KeyError( 'Unit ' + symbol + ' already defined' )
    if( name ) : _help.append( ( name, unit, symbol ) )

    if( isinstance( unit, str ) ) :
        if unit == 'as': # exception to handle atto-second: conflict with python keyword "as"
            unit = _unit_table[unit]
        else:
            unit = eval( unit, _unit_table )

        for cruft in ['__builtins__', '__args__'] :
            try :
                del _unit_table[cruft]
            except :
                pass

    unit.setSymbol( symbol )
    _unit_table[symbol] = unit

def _addPrefixedUnit( unit ) :
    _prefixed_symbols = []

    for prefix in _prefixes:
        symbol = prefix[0] + unit
        _addUnit( symbol, prefix[1] * _unit_table[unit] )
        _prefixed_symbols.append( symbol )

    for entry in _help:
        if isinstance( entry, tuple ) :
            if len( entry ) == 3:
                name0, unit0, symbol0 = entry
                i = _help.index( entry )
                if unit == symbol0:
                    _help.__setitem__( i, ( name0, unit0, symbol0, _prefixed_symbols ) )

# SI derived units; these automatically get prefixes
_help.append( '----------------------------------------------------\nDefined prefixes and units\n----------------------------------------------------' )
_help.append( 'This section outlines the prefixes and unit defined by the PQU module as well as the values used for the physical constant.' )
_help.append( '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nPrefixes\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' )
_help.append( '.. highlights::\n\n::' )
_help.append( '    ----  ------ ------\n' + '    name  factor symbol\n' + '    ----  ------ ------' )
_help_prefix = []

for _prefix in _prefixes :
    _help.append( ( '%-5s %-6.0e %-s' % ( _prefix[2], _prefix[1], _prefix[0] ), '', '' ) )

_help.append( '\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nDerived SI units\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' )
_help.append( 'These automatically get prefixes.' )
_help.append( '.. highlights::\n\n::' )
_help.append( '    %-26s %-32s %-s\n    %-26s %-32s %-s\n    %-26s %-32s %-s' % ( '--------------------------', '--------------------------------', '------', 'name', 'value', 'symbol', '--------------------------', '--------------------------------', '------' ) )

_unit_table['kg'] = PhysicalUnit( 'kg',   1., [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ] )

_addUnit( 'Hz', '1/s', 'Hertz' )
_addUnit( 'N', 'm*kg/s**2', 'Newton' )
_addUnit( 'Pa', 'N/m**2', 'Pascal' )
_addUnit( 'J', 'N*m', 'Joule' )
_addUnit( 'W', 'J/s', 'Watt' )
_addUnit( 'C', 's*A', 'Coulomb' )
_addUnit( 'V', 'W/A', 'Volt' )
_addUnit( 'F', 'C/V', 'Farad' )
_addUnit( 'ohm', 'V/A', 'Ohm' )
_addUnit( 'S', 'A/V', 'Siemens' )
_addUnit( 'Wb', 'V*s', 'Weber' )
_addUnit( 'T', 'Wb/m**2', 'Tesla' )
_addUnit( 'H', 'Wb/A', 'Henry' )
_addUnit( 'lm', 'cd*sr', 'Lumen' )
_addUnit( 'lx', 'lm/m**2', 'Lux' )
_addUnit( 'Bq', '1/s', 'Becquerel' )
_addUnit( 'Gy', 'J/kg', 'Gray' )
_addUnit( 'Sv', 'J/kg', 'Sievert' )

del _unit_table['kg']

for unit in list( _unit_table.keys( ) ) : _addPrefixedUnit( unit )

# NOTE everything below does not have prefixes added to units

# Dimensionless quantity
_unit_table[''] = PhysicalUnit( '', 1., [ 0, 0, 0, 0, 0, 0, 0, 0, 0 ] )

# Fundamental constants
_help.append( '\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nFundamental constants\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' )
_help.append( '.. highlights::\n\n::' )
_help.append( '    Fundamental constants ............................................ ' )

from . import pqu_constants

_unit_table['pi'] = math.pi
_addUnit( 'c',       pqu_constants.c,               'speed of light' )
_addUnit( 'mu0',     pqu_constants.mu0,             'permeability of vacuum' )
_addUnit( 'eps0',    '1 / mu0 / c**2',             'permittivity of vacuum' )
_addUnit( 'Grav',    pqu_constants.Grav,            'gravitational constant' )
_addUnit( 'hplanck', pqu_constants.hplanck,         'Planck constant' )
_addUnit( 'hbar',    'hplanck / ( 2 * pi )',        'Planck constant / 2pi' )
_addUnit( 'e',       pqu_constants.e,               'elementary charge' )
_addUnit( 'me',      pqu_constants.me,              'electron mass' )
_addUnit( 'mp',      pqu_constants.mp,              'proton mass' )
_addUnit( 'Nav',     pqu_constants.Nav,             'Avogadro number' )
_addUnit( 'k',       pqu_constants.k,           'Boltzmann constant' )

# Time units
_help.append( '\n\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nMore units\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' )
_help.append( '.. highlights::\n\n::' )
_help.append( '    Time units ....................................................... ' )

_addUnit( 'min', '60.*s', 'minute' )
_addUnit( 'h', '60.*min', 'hour' )
_addUnit( 'd', '24.*h', 'day' )
_addUnit( 'wk', '7.*d', 'week' )
_addUnit( 'yr', '365.25*d', 'year' )

# Length units
_help.append( '    Length units ..................................................... ' )

_addUnit( 'inch', '2.54 * cm', 'inch' )
_addUnit( 'ft', '12. * inch', 'foot' )
_addUnit( 'yd', '3. * ft', 'yard' )
_addUnit( 'mi', '5280. * ft', '(British) mile' )
_addUnit( 'nmi', '1852. * m', 'Nautical mile' )
_addUnit( 'Ang', '1.e-10 * m', 'Angstrom' )
_addUnit( 'lyr', 'c * yr', 'light year' )
_addUnit( 'au', '149597870700. * m', 'Astronomical unit' )
_addUnit( 'Bohr', '4. * pi * eps0 * hbar**2/me/e**2', 'Bohr radius' )

# Area units
_help.append( '    Area units ....................................................... ' )

_addUnit( 'ha', '10000*m**2', 'hectare' )
_addUnit( 'acres', 'mi**2/640', 'acre' )
_addUnit( 'b', '1.e-28*m**2', 'barn' )

_addPrefixedUnit( 'b' )

# Volume units
_help.append( '    Volume units ..................................................... ' )

_addUnit( 'l', 'dm**3', 'liter' )
_addUnit( 'dl', '0.1 * l', 'deci liter' )
_addUnit( 'cl', '0.01 * l', 'centi liter' )
_addUnit( 'ml', '0.001 * l', 'milli liter' )
_addUnit( 'tsp', pqu_constants.tsp, 'teaspoon' )
_addUnit( 'tbsp', '3. * tsp', 'tablespoon' )
_addUnit( 'floz', '2. * tbsp', 'fluid ounce' )
_addUnit( 'cup', '8. * floz', 'cup' )
_addUnit( 'pt', '16. * floz', 'pint' )
_addUnit( 'qt', '2. * pt', 'quart' )
_addUnit( 'galUS', '4. * qt', 'US gallon' )
_addUnit( 'galUK', pqu_constants.galUK, 'British gallon' )

# Mass units
_help.append( '    Mass units ....................................................... ' )

_addUnit( 'amu', pqu_constants.amu, 'atomic mass units' )
_addUnit( 'oz', pqu_constants.oz, 'ounce' )
_addUnit( 'lb', '16. * oz', 'pound' )
_addUnit( 'ton', '2000. * lb', 'ton' )

# Force units
_help.append( '    Force units ...................................................... ' )

_addUnit( 'dyn', '1.e-5 * N', 'dyne (cgs unit)' )

# Energy units
_help.append( '    Energy units ..................................................... ' )

_addUnit( 'erg', '1.e-7*J', 'erg (cgs unit)' )
_addUnit( 'eV', 'e*V', 'electron volt' )
_addUnit( 'Hartree', 'me*e**4/16/pi**2/eps0**2/hbar**2', 'Wavenumbers/inverse cm' )
_addUnit( 'Ken', 'k*K', 'Kelvin as energy unit' )
_addUnit( 'cal', pqu_constants.cal, 'thermochemical calorie' )
_addUnit( 'kcal', '1000.*cal', 'thermochemical kilocalorie' )
_addUnit( 'cali', pqu_constants.cali, 'international calorie' )
_addUnit( 'kcali', '1000.*cali', 'international kilocalorie' )
_addUnit( 'Btu', pqu_constants.Btu, 'British thermal unit' )

_addPrefixedUnit( 'eV' )

# Power units
_help.append( '    Power units ...................................................... ' )

_addUnit( 'hp', pqu_constants.hp, 'horsepower' )

# Pressure units
_help.append( '    Pressure units ................................................... ' )

_addUnit( 'bar', '1.e5*Pa', 'bar (cgs unit)' )
_addUnit( 'atm', pqu_constants.atm, 'standard atmosphere' )
_addUnit( 'torr', 'atm/760', 'torr = mm of mercury' )
_addUnit( 'psi', pqu_constants.psi, 'pounds per square inch' )

# Angle units
_help.append( '    Angle units ...................................................... ' )

_addUnit( 'deg', 'pi*rad/180.', 'degrees' )

_help.append( '    Temperature units ................................................ ' )
# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
Kelvin = _findUnit( 'K' )
_addUnit( 'degR', '(5./9.)*K', 'degrees Rankine' )
_addUnit( 'degC', PhysicalUnit( None, 1.0, Kelvin.powers, 273.15 ), 'degrees Celsius' )
_addUnit( 'degF', PhysicalUnit( None, 5./9., Kelvin.powers, 459.67 ), 'degrees Fahrenheit' )

_help.append( '\n' )
_help.append( '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nPrefixed units\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' )
_help.append( '.. highlights::\n\n::' )

del Kelvin

def description( ) :
    """
    Return a string describing all available units.

    :rtype: `str`
    """

    s = [] # collector for description text

    for entry in _help :
        if isinstance( entry, str ) :
            if entry != None:
                s.append( '\n'  + entry + '\n' )
            else:
                pass
        elif isinstance( entry, tuple ) :
            if len( entry ) == 4:
                name, unit, symbol, prefixes = entry
            elif len( entry ) == 3:
                name, unit, symbol = entry

            s.append( '    %-26s %-32s %-s\n' % ( name, unit, symbol ) )
        else:
            # impossible
            raise TypeError( 'wrong construction of _help list' )

    entry = '    ------------- ------ --------------\n    %-13s %-6s prefixed units\n    ------------- ------ --------------' % ( 'name', 'symbol' )
    s.append( '\n'  + entry + '\n' )

    for entry in _help:
        if isinstance( entry, str ) :
            pass
        elif isinstance( entry, tuple ) :
            if len( entry ) == 4:
                name, unit, symbol, prefixes = entry
                s.append( '    %-13s %-6s %-s\n    %-13s %-6s %-s\n\n' % ( name, symbol, ', '.join( prefixes[0:10] ), ' ', ' ', ', '.join( prefixes[10:20] ) ) )
    s.append( '\n' )

    return ''.join( s )

# add the description of the units to the module's doc string:
__doc__ += '\n' + description( )

def printPredefinedUnits( ) :

    units = []
    for symbol, unit in _unit_table.items( ) :
        if( symbol in [ 'pi' ] ) : continue
        s = "    %-10s  %20.12e  %20.12e" % ( symbol, unit.factor, unit.offset )
        for power in unit.powers : s += "  %3d" % power
        units.append( s )
    units.sort( )
    for unit in units : print(unit)

def convertUnits( unit, unitMap ) :
    '''
    This function converts *unit* to a new unit based on the key/value pairs in the dictionary *unitMap*,
    and it determines the factor that converts *unit* to the new unit.
    That is, it looks up each term in *unit* to see if it is a key in *unitMap*. If it is,
    it replaces that term with the value of the key and scales a factor appropriately.
    For example, if *unit* = 'b * eV' and unitMap = {'eV': 'MeV'} then the new unit returned
    will have 'eV' replaced with 'MeV' (i.e., it will be 'b * MeV' and the factor returned will be
    1e-6. Also if *unit* = 'b * eV' and unitMap = {'eV': 'MeV', 'b': 'mb', 'm': 'km'} then the new unit returned
    it will be 'mb * MeV' and the factor returned will be 1e-3. In this latter example, 'm' is not in the 
    original unit so is not used.

    :param unit:        A valid unit as a string.
    :param unitMap:     A dictionary containing old units as keys and their new unit as the key's value.

    :returns:           (string, factor)
    '''

    PU = _findUnit( unit )
    unitString = ''
    operator = ''
    for _unit in PU.symbols :
        try :
            __unit = unitMap[_unit]
        except :
            __unit = _unit
        unitString += "%s%s**(%s)" % ( operator, __unit, PU.symbols[_unit] )
        operator = ' * '
    newUnit = _findUnit( unitString )
    factor = PU.conversionFactorTo( newUnit )
    return( newUnit.symbol( ), factor )
