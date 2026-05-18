# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains useful fudge math routines that do not fit into any other module.
"""

import math
import numpy

from pqu import PQU as PQUModule
from xData import enums as xDataEnumsModule
from fudge.core.utilities import brb


def thickenXYList(list, tester, biSectionMax=6, interpolation=xDataEnumsModule.Interpolation.linlin):
    """
    This functions takes a list of (x,y) points and a function, tester.evaluateAtX, and bi-sectionally adds points to 
    obtain linear-linear tolerance of the returned list and tester.evaluateAtX to tester.relativeTolerance. At most 
    biSectionMax bi-sections are performed between each consecutive pair of inputted points. It is assumed that the 
    initial list of points and the function tester.evaluateAtX agree to tolerance tester.relativeTolerance. The instance
    tester must contain the members relativeTolerance and absoluteTolerance as well as the method evaluateAtX. The
    method evaluateAtX takes an x-value and returns its y-value.
    """

    def thickenXYList2(interpolation, xl, yl, xu, yu, newList, tester, level):

        if level == biSectionMax: return
        level += 1
        if interpolation == xDataEnumsModule.Interpolation.linlin or interpolation == xDataEnumsModule.Interpolation.loglin:
            xMid = 0.5 * (xl + xu)
        else:
            xMid = math.sqrt(xl * xu)

        if interpolation == xDataEnumsModule.Interpolation.linlin or interpolation == xDataEnumsModule.Interpolation.linlog:
            yMid = 0.5 * (yl + yu)
        else:
            yMid = math.sqrt(yl * yu)

        y = tester.evaluateAtX(xMid)

        dy = abs(y - yMid)
        if dy > abs(y * tester.relativeTolerance) and dy > tester.absoluteTolerance:
            newList.append([xMid, y])
            thickenXYList2(interpolation, xl, yl, xMid, y, newList, tester, level)
            thickenXYList2(interpolation, xMid, y, xu, yu, newList, tester, level)

    if len(list) < 2: raise Exception("len( list ) = %d < 2" % len(list))
    newList = []
    for i1, xy in enumerate(list):
        x2, y2 = xy
        if i1 > 0: thickenXYList2(interpolation, x1, y1, x2, y2, newList, tester, 0)
        newList.append([x2, y2])
        x1, y1 = x2, y2
    newList.sort()
    return newList


def checkForNaN(v, str, printErrors=True, indentation="", messages=None):
    if v != v:
        s = "%s: value is nan" % str
        if messages is not None: messages.append(s)
        if printErrors: print(indentation + s)


def isNumber(n):
    return type(n) in (int, float, numpy.float64)


def checkNumber(v, str="", printErrors=True, indentation="", messages=None, maxAbsFloatValue=None):
    checkForNaN(v, str, printErrors=printErrors, indentation=indentation, messages=messages)
    if maxAbsFloatValue is not None:
        if abs(v) > maxAbsFloatValue:
            s = "%s: abs of number %s exceeds maxAbsFloatValue = %s" % (str, v, maxAbsFloatValue)
            if messages is not None: messages.append(s)
            if printErrors: print(indentation + s)


def getValue(n):
    if isNumber(n): return n
    if isinstance(n, PQUModule.PQU): return n.getValue()
    raise Exception('Invalid number object = %s' % brb.getType(n))


def RoundToSigFigs(x, sigfigs):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.

    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.
    """
    # The following constant was computed in maxima 5.35.1 using 64 bigfloat digits of precision
    __logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1

    if not (type(sigfigs) is int or numpy.issubdtype(sigfigs, numpy.integer)):
        raise TypeError("RoundToSigFigs: sigfigs must be an integer.")

    if not numpy.all(numpy.isreal(x)):
        raise TypeError("RoundToSigFigs: all x must be real.")

    if sigfigs <= 0:
        raise ValueError("RoundtoSigFigs: sigfigs must be positive.")

    mantissas, binaryExponents = numpy.frexp(x)

    decimalExponents = __logBase10of2 * binaryExponents
    intParts = numpy.floor(decimalExponents)

    mantissas *= 10.0 ** (decimalExponents - intParts)

    return numpy.around(mantissas, decimals=sigfigs - 1) * 10.0 ** intParts
