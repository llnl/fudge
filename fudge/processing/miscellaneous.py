# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

""" 
This module contains functions used internally by FUDGE to multi-group 1d functions.
        
This module contains the following functions: 
        
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Class                                 | Description                                                                       |
    +=======================================+===================================================================================+
    | _toLinear                             | This function returns an XYs1d linear version of its argument.                    |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | _groupFunctionsAndFluxInit            | Takes raw transfer matrices and wraps them in a                                   |
    |                                       | :py:class:`multiGroupModule.Form` instance.                                       |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | _mutualifyGrouping3Data               | Reads in a file and returns a :py:class:`Groups` instance.                        |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | groupOneFunctionAndFlux               | Reads in a file and returns a :py:class:`Groups` instance.                        |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | groupTwoFunctionsAndFlux              | Reads in a file and returns a :py:class:`Groups` instance.                        |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | groupFunctionCrossSectionAndFlux      | Reads in a file and returns a :py:class:`Groups` instance.                        |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""

accuracy = 1e-6

from xData import enums as xDataEnumsModule
from xData import link as linkModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import constant as constantModule
from xData import regions as regionsModule

from . import group as groupModule

def _toLinear(func):
    """
    This function returns an XYs1d linear version of :math:`func`.
    This function if for internal use.

    :param func:    Function to convert to a lin-lin XYs1d instance.

    :returns:       An XYs1d instace.
    """

    if isinstance(func, XYs1dModule.XYs1d):
        if func.interpolation == xDataEnumsModule.Interpolation.linlin:
            return func
        func = func.copy( )
        return( func.toPointwise_withLinearXYs( accuracy = accuracy, upperEps = 1e-8 ) )
    elif isinstance(func, (regionsModule.Regions1d, constantModule.Constant1d)):
        return( func.toPointwise_withLinearXYs( accuracy = accuracy, upperEps = 1e-8 ) )

    print(type(func))
    raise Exception( 'FIX ME' )

def _groupFunctionsAndFluxInit( style, tempInfo, f1 ) :
    """
    This function returns the flux and multi-group boundaries from *style* and *tempInfo*, and domainSlices the
    flux to the domain of :math:`f1` if it is not None.  This function if for internal use.

    :param style:           This is the multi-group style for the multi-group data.
    :param tempInfo:        This is a dictionary with needed data.
    :param f1:              A 1-d function whose domain is used to set the domain limit of the returned flux.

    :returns:               The multi-group boundaries and flux.
    """

    reactionSuite = tempInfo['reactionSuite']
    groupBoundaries = style.transportables[reactionSuite.projectile].group.boundaries
    flux = style.flux.getFluxAtLegendreOrder( 0 )
    if( f1 is not None ) :
        domainMin, domainMax = f1.domainMin, f1.domainMax
        flux = flux.domainSlice( domainMin, domainMax )
    return( groupBoundaries, flux )

def _mutualifyGrouping3Data( f1, f2, f3, printMutualDomainWarning = False ) :
    r"""
    This function returns versions of :math:`f1`, :math:`f2` and :math:`f3` whose domains are all mutual.
    This function if for internal use.

    :param f1:                          One of the functions that is mutualified.
    :param f2:                          One of the functions that is mutualified.
    :param f3:                          One of the functions that is mutualified.
    :param printMutualDomainWarning:    If True, a warning is printed if the domain of :math:`f1`, :math:`f2` and :math:`f3` are not mutual.

    :returns:                           Versions of :math:`f1`, :math:`f2` and :math:`f3` whose domains are all mutual.
    """

    domainMin1, domainMax1 = f1.domainMin, f1.domainMax
    domainMin2, domainMax2 = f2.domainMin, f2.domainMax
    if( len( f3 ) == 0 ) :
        domainMin3, domainMax3 = domainMin2, domainMax2
    else :
        domainMin3, domainMax3 = f3.domainMin, f3.domainMax
    if( ( domainMin1 != domainMin2 ) or ( domainMin1 != domainMin3 ) or ( domainMax1 != domainMax2 ) or ( domainMax1 != domainMax3 ) ) :
        domainMin, domainMax = max( domainMin1, domainMin2, domainMin3 ), min( domainMax1, domainMax2, domainMax3 )
        if( printMutualDomainWarning ) :
            print("WARNING: making domains mutual for grouping,", domainMin1, domainMin2, domainMin3, domainMax1, domainMax2, domainMax3)
        f1 = f1.domainSlice( domainMin = domainMin, domainMax = domainMax )
        f2 = f2.domainSlice( domainMin = domainMin, domainMax = domainMax )
        if( len( f3 ) > 0 ) : f3 = f3.domainSlice( domainMin = domainMin, domainMax = domainMax )
    return( f1, f2, f3 )

def groupOneFunctionAndFlux( style, tempInfo, f1, styleFilter = None ) :
    r"""
    This function mulit-groups :math:`f1`. This function if for internal use.

    :param style:           This is the multi-group style for the multi-group data.
    :param tempInfo:        This is a dictionary with needed data.
    :param f1:              The function that to multi-grouped.
    :param styleFilter:     See method findFormMatchingDerivedStyle of the class Style.

    :returns:               A list like object of the multi-group values.
    """

    if( isinstance( f1, linkModule.Link ) ) : f1 = style.findFormMatchingDerivedStyle( f1.link, styleFilter = styleFilter )
    f1 = _toLinear( f1 )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux, norm = tempInfo['groupedFlux'] ) )

def groupTwoFunctionsAndFlux( style, tempInfo, f1, f2, norm = None, printMutualDomainWarning = False ) :
    r"""
    This function mulit-groups the product of two functions, :math:`f1 \times f2`. Typically, :math:`f2` is a cross section.
    This function if for internal use.

    :param style:                       This is the multi-group style for the multi-group data.
    :param tempInfo:                    This is a dictionary with needed data.
    :param f1:                          One of the two functions that represent the product this is multi-grouped.
    :param f2:                          One of the two functions that represent the product this is multi-grouped.
    :param norm:                        A normalization to divided the multi-groups by.
    :param printMutualDomainWarning:    If True, a warning is printed if the domain of :math:`f1` and :math:`f2` are not mutual.

    :returns:                           A list like object of the multi-group values.
    """

    f1 = _toLinear( f1 )
    f2 = _toLinear( f2 )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )
    f1, f2, flux = _mutualifyGrouping3Data( f1, f2, flux, printMutualDomainWarning = printMutualDomainWarning )
    return( f1.groupThreeFunctions( groupBoundaries, flux, f2, norm = norm ) )

def groupThreeFunctionsAndFlux( style, tempInfo, f1, f2, f3, norm = None, printMutualDomainWarning = False ) :
    r"""
    This function mulit-groups the product of three functions, :math:`f1 \times f2 \times f3`.
    This function if for internal use.

    :param style:                       This is the multi-group style for the multi-group data.
    :param tempInfo:                    This is a dictionary with needed data.
    :param f1:                          One of the three functions that represent the product that is multi-grouped.
    :param f2:                          One of the three functions that represent the product that is multi-grouped.
    :param f3:                          One of the three functions that represent the product that is multi-grouped.
    :param norm:                        A normalization to divided the multi-groups by.
    :param printMutualDomainWarning:    If True, a warning is printed if the domain of :math:`f1` and :math:`f2` are not mutual.

    :returns:                           A list like object of the multi-group values.
    """

    f1 = _toLinear( f1 )
    f2 = _toLinear( f2 )
    f3 = _toLinear( f3 )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )

    minEnergy = max(f1.domainMin, f2.domainMin, f3.domainMin, flux.domainMin, groupBoundaries.domainMin)
    maxEnergy = min(f1.domainMax, f2.domainMax, f3.domainMax, flux.domainMax, groupBoundaries.domainMax)
    f1 = f1.domainSlice(domainMin=minEnergy, domainMax=maxEnergy)
    f2 = f2.domainSlice(domainMin=minEnergy, domainMax=maxEnergy)
    f3 = f3.domainSlice(domainMin=minEnergy, domainMax=maxEnergy)
    flux = flux.domainSlice(domainMin=minEnergy, domainMax=maxEnergy)

    return( f1.groupFourFunctions( groupBoundaries, flux, f2, f3, norm = norm ) )

def groupFunctionCrossSectionAndFlux( cls, style, tempInfo, f1, printMutualDomainWarning = False ) :
    """
    This function multi-groups :math:`f1` with the product of the cross section in *tempInfo*.
    This function if for internal use.

    :param cls:                         The Gridded1d class to return that has the multi-group data.
    :param style:                       This is the multi-group style for the multi-group data.
    :param tempInfo:                    This is a dictionary with needed data.
    :param f1:                          One of the two functions that represent the product that is multi-grouped.
    :param printMutualDomainWarning:    If True, a warning is printed if the domain of :math:`f1` and :math:`f2` are not mutual.

    :returns:                           An instance of *cls*.
    """

    crossSection = tempInfo['crossSection']
    norm = tempInfo['groupedFlux']
    grouped = groupTwoFunctionsAndFlux( style, tempInfo, f1, crossSection, norm = norm, printMutualDomainWarning = printMutualDomainWarning )

    unit = crossSection.axes[0].unit
    if( f1.axes[0].unit != '' ) : unit = '%s * %s' % ( f1.axes[0].unit, unit )
    axis = axesModule.Axis( label = '%s * %s' % ( f1.axes[0].label, crossSection.axes[0].label ), unit = unit, index = 0 )
    axes = f1.axes.copy( )
    axes[0] = axis
    return( groupModule.toMultiGroup1d( cls, style, tempInfo, axes, grouped, zeroPerTNSL = tempInfo['zeroPerTNSL'] ) )
