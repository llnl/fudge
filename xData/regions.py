# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the all the classes for handling regions functions. A regions container is designed to handle
a function that has 1 or more of the following::

    -) an outer dimension that is multi-value at one of more domain values (e.g., :math:`f(x)` != :math:`f(x)` for one or more values of x).
    -) the interpolation of the outer dimension of the function is different in different regions of the domain.
    -) different functions are used in different regions of the domain (e.g., Legendre functions may be used for the outer domain
        < 5 and XYs1d functions otherwise).

This module contains the following classes:
        
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | Regions           | This is the base class for the other regions classes.                                                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Regions1d         | This is a class for a 1 dimensional regions function.                                                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | RegionsMultiD     | This is the base class for the 2 and dimensional regions classes.                                     |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Regions2d         | This is a class for a 2 dimensional regions function.                                                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Regions3d         | This is a class for a 3 dimensional regions function.                                                 |
    +-------------------+-------------------------------------------------------------------------------------------------------+
"""

import abc
import math

from pqu import PQU as PQUModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from . import enums as enumsModule
from . import XYs1d as XYs1dModule
from . import series1d as series1dModule
from . import base as baseModule

domainEpsilon = 1e-15

class Regions( baseModule.XDataFunctional ) :
    """
    This is the absract base class for the other regions classes.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
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
    | __regions         | This is a python list of the regions that have been added.    |
    +-------------------+---------------------------------------------------------------+
    """

    ancestryMembers = baseModule.XDataFunctional.ancestryMembers                # See note in multiD_XYs.py about setting ancestryMembers.

    def __init__( self, axes = None, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None ) :
        """
        :param axes:                This is the axes member.
        :param index:               This is the index member.
        :param valueType:           This describes the type of data (i.e., float, int) returned by the function.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is
                                    embedded in a high dimensional function.
        :param label:               This is the label member.
        """

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)
        self.__regions = []

    def __len__( self ) :
        """
        This method returns the number of regions in self.

        :returns:       A python int.
        """

        return( len( self.__regions ) )

    @staticmethod
    @abc.abstractmethod
    def allowedSubElements( ):
        """
        This method is abstract and must be over loaded by each sub-class. This method returns a list of all functions of the 
        same dimension as the sub-class that can be added as a region of the sub-class.
        """

        pass

    def __getitem__( self, index ) :
        """
        This method returns the (i-1)^th region of self.

        :param index:   The index of the region to return.

        :returns:       A function that has the same dimension as *self*.
        """

        return( self.__regions[index] )

    def __setitem__( self, index, region ) :
        """
        This method set the region at index *index* to *region*. The *region* must be an instance of 
        base.XDataFunctional with the same dimension as self.
        The region must fix between any abutting regions. That is, if *index* is greater than 0, the domainMin of 
        *region* must be equal to the domainMax of the prior region. And, if *index* is less than the number of regions
        in *self* minus 1, then domainMax of *region* must be equal to the domainMin of the next region.

        :param index:       The index where *region* is inserted.
        :param region:      This must be an instance of a class in :py:func:`allowedSubElements`.
        """

# BRB need to check axes.
        if( not( isinstance( region, self.allowedSubElements( ) ) ) ) : raise TypeError( 'Invalid class for insertion: %s' % region.__class__ )

        n1 = len( self )
        if( index < 0 ) : index += n1
        if( not( 0 <= index <= n1 ) ) : raise IndexError( 'Index = %s not in range 0 <= index <= %d' % ( index, n1 ) )
        if( index > 0 ) :
            if( not( math.isclose( self.__regions[index-1].domainMax, region.domainMin ) ) ) :
                raise ValueError( "Prior region's domainMax %s != new region's domainMin = %s" % ( self.__regions[index-1].domainMax, region.domainMin ) )
        if( ( n1 > 0 ) and ( index < ( n1 - 1 ) ) ) :
            if( not math.isclose( self.__regions[index+1].domainMin, region.domainMax ) ) :
                raise ValueError(  "Next region's domainMin %s != new region's domainMax = %s" % ( self.__regions[index-1].domainMin, region.domainMax ) )

        if( index == n1 ) :
            self.__regions.append( region )             # Append to the end.
        else :
            self.__regions[index] = region              # Replaces the current contents of index with region.

        region.setAncestor( self )
        region.index = index

    def __add__(self, other):
        """
        """

        self2, other2 = self.copyToCommonRegions(other)
        for i1, region in enumerate(self2):
            self2[i1] = region + other2[i1]

        return self2

    __radd__ = __add__

    def __sub__( self, other ) :

        self2, other2 = self.copyToCommonRegions( other )
        for i1, region in enumerate( self2 ) : self2[i1] = region - other2[i1]
        return( self2 )

    @property
    def domainMin( self ) :
        """
        This method returns the minimum domain value for *self*.
        """

        return( self.__regions[0].domainMin )

    @property
    def domainMax( self ) :
        """
        This method returns the maximum domain value for *self*.
        """

        return( self.__regions[-1].domainMax )

    @property
    def domainGrid( self ) :
        """
        This method returns all domain values for *self* as a python list.

        :returns:           A python list.
        """

        grid = set()
        for region in self: grid.update( region.domainGrid )
        return sorted( grid )

    @property
    def domainUnit( self ) :
        """
        This method returns the domain unit for *self*.
        """

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :

        return( min( [ region.rangeMin for region in self ] ) )

    @property
    def rangeMax( self ) :

        return( max( [ region.rangeMax for region in self ] ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    @property
    def regions( self ) :
        """Returns self's __region."""

        return( self.__regions )

    @property
    def functionNdsName( self ) :
        """Returns the node name for the child "function#ds"."""

        return( "function%dds" % self.dimension )

    def append(self, curve):
        """
        Appends *curve* to the end of *self*. If *curve* is a **Regions** instance, then **append** is call for 
        each region in curve.
        """

        if isinstance(curve, Regions):
            for region in curve:
                self.append(region.copy())
        else:
            self[len(self)] = curve
            curve.keyName = 'index'

    def prepend(self, curve):
        """
        Prepends *curve* to the beginning of *self*. If *curve* is a **Regions** instance, then **append** is call 
        for each region in curve.
        """

        if isinstance(curve, Regions):
            for index in range(len(curve), 0, -1):
                self.prepend(curve[index].copy())
        else:
            if not isinstance(curve, self.allowedSubElements()):
                raise TypeError('Invalid class for insertion: %s' % region.__class__)

            if len(self) > 0:
                if not math.isclose(curve.domainMax, self[0].domainMin):
                    raise ValueError('''Prepending region's domainMax %s != first region's domainMin = %s''' % (region.domainMax, self[0].domainMin))

            self.__regions.insert(0, curve)
            curve.setAncestor(self)
            curve.index = 0

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        factors = self.axes.convertUnits( unitMap )
        for region in self : region.convertUnits( unitMap )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :
        # FIXME some of this should probably move to 'returnAsClass' method

        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        newRegions = self.__class__( axes = axes, index = self.index, valueType = self.valueType, outerDomainValue = self.outerDomainValue, label = self.label )
        for child in self : newRegions.append( child.copy( ) )
        return( newRegions )

    def hasData(self):
        """
        Returns True if one of self's regions hasData returns True and False otherwise.
        """

        for region in self:
            if region.hasData():
                return True

        return False

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) :
        """
        Splits the region containing domainValue into two regions.
        """

        for i1, region in enumerate( self ) :
            domainMin, domainMax = region.domainMin, region.domainMax
            if( domainMin < domainValue < domainMax ) :
                r1, r2 = region.splitInTwo( domainValue, epsilon = domainEpsilon )
                self.__regions[i1] = r2
                r1.setAncestor( self )
                self.__regions.insert( i1, r1 )
                r2.setAncestor( self )
                return

    def domainUnitConversionFactor( self, unitTo ) :
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A float.
        """

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def domainSlice(self, domainMin=None, domainMax=None, fill=1, dullEps=0.0, **kwargs):
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.
        Result may be a regions container, or it may be XYs, multiD_XYs, etc.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self,
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        domainMin = max( domainMin, self.domainMin )
        if( domainMax is None ) : domainMax = self.domainMax
        domainMax = min( domainMax, self.domainMax )

        for ridx1, region in enumerate( self ) :
            if( region.domainMax  > domainMin ) : break
        for ridx2, region in enumerate( self ) :
            if( region.domainMax >= domainMax ) : break

        if ridx1 == ridx2:  # only one region left after slicing, return as XYs1d or multiD_XYs
            return self[ridx1].domainSlice( domainMin=domainMin, domainMax=domainMax, fill=fill, dullEps=dullEps )
        else:
            newRegions = self.__class__( axes = self.axes, index = self.index, valueType = self.valueType, outerDomainValue = self.outerDomainValue, label = self.label )
            newRegions.append( self[ridx1].domainSlice( domainMin = domainMin, domainMax = self[ridx1].domainMax, fill = fill, dullEps = dullEps ) )
            for idx in range(ridx1+1,ridx2):
                newRegions.append(self[idx])
            newRegions.append( self[ridx2].domainSlice(
                domainMin=self[ridx2].domainMin, domainMax=domainMax, fill=fill, dullEps=dullEps ) )

            return newRegions

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

        if fixToDomain == enumsModule.FixDomain.lower:
            regions1 = self.domainSlice(domainMin = domainMin, fill = True)
        elif fixToDomain == enumsModule.FixDomain.upper:
            regions1 = self.domainSlice(domainMax = domainMax, fill = True)
        else:
            regions1 = self.domainSlice(domainMin = domainMin, domainMax = domainMax, fill = True)

        self.__regions = []
        if isinstance(regions1, Regions):
           for index, region in enumerate(regions1): self[index] = region
        else:
            if len(regions1) == 0: return 1             # Caller beware.
            self[0] = regions1

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def evaluate( self, domainValue ) :
        """
        Evaluate function at requested domain point. If at discontinuity, return upper region's value.

        :param domainValue:
        :return interpolated point at domainValue:
        """

        for region in self.__regions :
            if( domainValue < region.domainMax ) : return( region.evaluate( domainValue ) )

        return( self.__regions[-1].evaluate( domainValue ) )  # Domain value is above the last region. Let the last region determine what to do.

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :
        """
        Finds all instances of class *cls* in self's children, grand-children, etc.
        """

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for region in self :
            if( isinstance( region, cls ) ) : foundInstances.append( region )
            foundInstances += region.findInstancesOfClassInChildren( cls, level = level )

        return( foundInstances )

    def integrate( self, **limits ):
        """
        Integrate a piecewise function. Supports limits for each axis.
        Example::

            regions.integrate( energy_in = ('1e-5 eV', '10 eV'), energy_out = ('1 keV', '10 keV') )

        :param limits: dictionary containing limits for each independent axis (keyed by axis label or index).
                       If an independent axis is missing from the dictionary, integrate over the entire domain of that axis.

        :return: float
        """

        from . import multiD_XYs as multiD_XYsModule

        integral = 0
        for region in self:
            if( isinstance( region, ( XYs1dModule.XYs1d, series1dModule.Series1d ) ) ) :
                Min, Max = None, None
                if self.axes[-1].label in limits:
                    Min, Max = limits.pop( self.axes[-1].label )
                elif self.axes[-1].index in limits:
                    Min, Max = limits.pop( self.axes[-1].index )

                integral += region.integrate( domainMin = Min, domainMax = Max )
            elif isinstance( region, multiD_XYs1dModule.XYsnd ):
                integral += region.integrate( **limits )
            else:
                raise TypeError( "Unsupported class for integration: %s" % type( region ) )
        return( integral )

    def toString( self ) :

        s1 = ''
        for i1, region in enumerate( self ) : s1 += '# region %s\n' % i1 + region.toString( )
        return( s1 )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        indent3 = indent2 + kwargs.get('incrementalIndent', '  ')
        if formatVersion == GNDS_formatVersionModule.version_1_10: indent3 = indent2

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if self.isPrimaryXData():
            if self.axes is not None: XML_strList += self.axes.toXML_strList(indent2, **kwargs)

        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList.append('%s<%s>' % ( indent2, self.functionNdsName ))
        for region in self.__regions: XML_strList += region.toXML_strList(indent3, **kwargs)
        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList[-1] += "</%s>" % self.functionNdsName

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

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        regions1 = cls(**attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(regions1, node, xPath, linkData, **kwargs)

        if len(extraNodes) > 0:                                                     # The next few line support GNDS 1.10 and 2.0.
            for index, extraNode in enumerate(extraNodes):
                if extraNode.tag == regions1.functionNdsName: break
            if extraNode.tag == regions1.functionNdsName: extraNodes = extraNodes.pop(index)

        if 'axes' not in kwargs: kwargs['axes'] = regions1.axes
        allowedSubElements = cls.allowedSubElements()
        for child in extraNodes:
            childClass = None
            for allowedChildClass in allowedSubElements:
                if allowedChildClass.moniker == child.tag:
                    childClass = allowedChildClass
                    break
            if childClass is None: raise TypeError('Invalid child "%s" for node "%s"' % (child.tag, cls.moniker))
            xdata = childClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            regions1.append(xdata)

        xPath.pop()

        return regions1

class Regions1d( Regions ) :
    """
    This class supports storing a 1d function as a list of abutting 1d functions.
    """

    moniker = 'regions1d'
    dimension = 1

    def __neg__(self):
        """
        This method returns self with all regions negated.
        """

        regions1d = self.__class__(axes=self.axes)
        for region in self:
            regions1d.append(-region)

        return regions1d

    def __mul__( self, other ) :
        """
        This method multiplies *self* by *other* and returns the product. Here, *other* can be an int for a float.

        :param other:       Instance to multiply *self* by.

        :return:            A :py:class:`Regions1d` instance.
        """

        _self, _other = self.copyToCommonRegions( other )
        _regions1d = self.__class__(axes=self.axes)
        for index, region1 in enumerate( _self ) :
            region2 = _other[index]
            region = region1 * region2
            _regions1d.append( region )

        return( _regions1d )

    __rmul__ = __mul__

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance.

        :param asLinLin:    If **True**, the data have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps:    Used to dull the lower point between two regions.
        :param upperEps:    Used to dull the upper point between two regions.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xys1d = self.toPointwiseLinear(accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        return xys1d

    def copyToCommonRegions(self, other, epsilon = domainEpsilon):
        """
        This method returns two :py:class:`Regions1d` instances that yield the same function as *self* and *other*,
        but with *self* and *other* broken up so that they have the same number of regions and the regions align.
        Currenlty, *other* can only be a :py:class:`Regions1d` or :py:class:`XYs1dModule.XYs1d` instance.

        :param other:       A :py:class:`Regions1d` or :py:class:`XYs1dModule.XYs1d` instance.
        :param epsilon:     Two boundaries are considered to be the same if they are relativley within this value.

        :returns:           Two :py:class:`Regions1d` instances.
        """

        self2 = self.copy()
        if isinstance(other, (float, int)):
            other2 = self.toLinearXYsClass()(data=[[self.domainMin, other], [self.domainMax, other]])
        elif self.dimension != other.dimension:
            raise ValueError('self.dimension = %s not equal to other.dimension = %s' % (self.dimension, other.dimension))
        else:
            other2 = other.copy()

        if len(self2) == 0:
            if len(other2) == 0:
               return self2, other2

            data = [[other2.domainMin, 0.0], [other2.domainMax, 0.0]]
            xys1d = self2.toLinearXYsClass()(data=data, axes=self2.axes)
            if not isinstance(self2, Regions1d):
                self2 = Regions1d(axes=self2.axes)
            self2.append(xys1d)
        elif len(other2) == 0:
            data = [[self2.domainMin, 0.0], [self2.domainMax, 0.0]]
            xys1d = other2.toLinearXYsClass()(data=data, axes=other2.axes)
            if not isinstance(other2, Regions1d):
                other2 = Regions1d(axes=other2.axes)
            other2.append(xys1d)

        if isinstance(other2, XYs1dModule.XYs1d):
            temp = self.__class__()
            temp.append(other2)
            other2 = temp 
        elif not isinstance(other2, Regions):
            raise NotImplementedError('object of instance "%s" not implemented' % other2.__class__)

        if   self2.domainMin < other2.domainMin:
            region1 = other2[0].copy()
            region1.setData([[self2.domainMin, 0], [other2.domainMin, 0]])
            other2.prepend(region1)
        elif self2.domainMin > other2.domainMin:
            region1 = self2[0].copy()
            region1.setData([[other2.domainMin, 0], [self2.domainMin, 0]])
            self2.prepend(region1)

        if   self2.domainMax > other2.domainMax:
            region1 = other2[-1].copy()
            region1.setData([[other2.domainMax, 0], [self2.domainMax, 0]])
            other2.append(region1)
        elif self2.domainMax < other2.domainMax:
            region1 = self2[-1].copy()
            region1.setData([[self2.domainMax, 0], [other2.domainMax, 0]])
            self2.append(region1)

        boundaries = set()
        for region in self2[1:]:
            boundaries.add(region.domainMin)
        for region in other2[1:]:
            boundaries.add(region.domainMin)
        boundaries = sorted(boundaries)

        count = 0
        priorBoundary = None
        boundariesToMove = []
        for boundary in boundaries:
            if priorBoundary is not None:
                if (boundary - priorBoundary) < epsilon * max(abs(boundary), abs(priorBoundary)):
                    boundariesToMove.append([boundary, priorBoundary])
                    boundary = priorBoundary
                    count += 1
                    if count > 1:
                        raise ValueError('more than one boundary within epsilon = %s of %s' % (epsilon, priorBoundary))
                else:
                    count = 0
            priorBoundary = boundary

        for boundary, priorBoundary in boundariesToMove:
            boundaries.remove(boundary)

        def processBoundaries(regions_, boundaries, boundariesToMove):
            """
            This method adds regions to *regions_* so that its has a region boundary for each value in boundaries.

            :param regions_:            This is the regions to work on.
            :param boundaries:          This is the list of boundaries that the final regions shall have.
            :param boundariesToMove:    This is the list of (new, old) boundaries in *regions_* where old is tweaked to become new.
            """

            for region in regions_:
                domainMin, domainMax = region.domainMin, region.domainMax
                for boundary, priorBoundary in boundariesToMove:
                    if   domainMin == boundary:
                        region.tweakDomain(domainMin = priorBoundary, epsilon = epsilon)
                    elif domainMax == boundary:
                        region.tweakDomain(domainMax = priorBoundary, epsilon = epsilon)

            for boundary in boundaries:
                for region in regions_:
                    if region.domainMin < boundary < region.domainMax:
                        regions_.splitInTwo(boundary, epsilon = epsilon)
                        break

        if len(boundaries) > 0:
            processBoundaries(self2, boundaries, boundariesToMove)
            processBoundaries(other2, boundaries, boundariesToMove)

        return self2, other2

    def copyDataToXsAndYs(self):
        """
        This method copies data in *self* to a list of xs and ys. Only works if each region instance has a "copyDataToXsAndYs" method.

        :returns:   Two python lists, one with the x-values and one with the y-values.
        """

        xs = []
        ys = []
        for region in self:
            subXs, subYs = region.copyDataToXsAndYs()
            xs += subXs
            ys += subYs
        return xs, ys

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        This method normalizes the data in *self*, ergo, scaled so that its integral is 1.0.
        If *insitu* is True, *self* is normalized and returned, otherwise a copy of *self* is normalizes and the copy is returned.
        The dimension argument is currently ignored, but kept to be compatible with calling from XYsnd.normalize.

        :param insitu:          If True, *self* is normalized and returned.
        :param dimension:       The argument is ignored.

        :returns:               A *self* like instance.
        """

        factor = 1.0 / self.integrate()
        if( insitu ) :
            copy = self
        else :
            copy = self.copy()

        for region in copy : region *= factor
        return( copy )

    def plot(self, **kwargs):
        """
        Calls multiPlot with only *self* as a curve.
        """

        self.multiPlot([self], **kwargs)

    def toPointwiseLinear( self, **kwargs ) :
        """
        This method returns the results of calling :py:func:`toPointwise_withLinearXYs`.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        return( self.toPointwise_withLinearXYs( **kwargs ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method converts the regions of *self* into a single :py:class:`XYs1dModule.XYs1d` instance that has 'lin-lin' interpolation. At the
        boundary between two abutting regions, the x-values are the same, which is not allowed for an :py:class:`XYs1dModule.XYs1d` instance
        so some points may be added or removed as specified by data in the argument *kwargs* as described in the following:

            - accuracy:    This controls how many points are added when switching interpolation.
            - lowerEps:    Along with 'upperEps', used to smear the x-values at a boundary as follows. Let :math:`(x, y_l)`
                           and :math:`(x, y_u)` be the abutting points for two abutting regions. If :math:`y_l = y_u` then the 
                           point :math:`(x, y_u)` is removed.  Otherwise, if( lowerEps > 0 ) the point :math:`(x, y_l)` is 
                           moved to :math:`x' = x * ( 1 - lowerEps )` (or :math:`x' = x * ( 1 + lowerEps )` if :math:`x < 0`) 
                           and the :math:`y` value is interpolated at :math:`x'`. If :math:`x'` is less than the x-value of the point 
                           below :math:`(x', y_l)` and ``removeOverAdjustedPoints`` is True then the point :math:`(x, y_l)` is removed; 
                           otherwise, a raise is executed. Similarly for upperEps and the point :math:`(x, y_u)`.
            - upperEps:    See description of lowerEps above.

        :param kwargs:      A dictionary of data needed by *self*.

        :returns:           An instance of :py:class:`XYs1dModule.XYs1d`.
        """

        def getAdjustedX( x, eps ) :
            """
            This method returns an x-value that is shifted from *x* by *eps*. If *eps* is negative, the returned value will be
            left than *x*; otherwise, it will be greater than.  Thie method is for internal use.

            :param x:       The initial x-value.
            :param eps:     The small amount to move the new x-value from *x*.

            :returns:       A python float.
            """

            if( x == 0. ) :
                x_ = eps
            elif( x < 0 ) :
                x_ = x * ( 1. - eps )
            else :
                x_ = x * ( 1. + eps )
            return( x_ )

        pointwiseClass = self.toLinearXYsClass( )
        if( len( self.regions ) == 0 ) : return( pointwiseClass( data = [], axes = self.axes ) )

        arguments = self.getArguments( kwargs, { 'accuracy' : XYs1dModule.defaultAccuracy, 'lowerEps' : 0, 'upperEps' : 0, 
                'removeOverAdjustedPoints' : False, 'axes' : None } )
        accuracy = arguments['accuracy']
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']
        removeOverAdjustedPoints = arguments['removeOverAdjustedPoints']
        axes = arguments['axes']

        if( lowerEps < 0. ) : raise ValueError( 'lowerEps = %s must >= 0.' % lowerEps )
        if( upperEps < 0. ) : raise ValueError( 'upperEps = %s must >= 0.' % upperEps )
        if( ( lowerEps == 0. ) and ( upperEps == 0. ) ) : raise ValueError( 'lowerEps and upperEps cannot both be 0.' )

        xys = []
        for iRegion, region in enumerate( self.regions ) :
            _region = region.changeInterpolation(enumsModule.Interpolation.linlin, accuracy, lowerEps = 2 * lowerEps, upperEps = 2 * upperEps )
            _region = _region.copyDataToXYs( )
            if len(_region) < 2: continue
            if iRegion > 0 and len(xys) > 0:
                x12, y12 = xys[-1]
                x21, y21 = _region[0]
                if( y12 == y21 ) :              # Remove first point of region as it is the same as the last point.
                    del _region[0]
                else :
                    if( lowerEps != 0. ) :
                        x11, y11 = xys[-2]
                        x = getAdjustedX( x12, -lowerEps )
                        if( x <= x11 ) :
                            if( removeOverAdjustedPoints ) :
                                del xys[-1]
                            else :
                                raise ValueError( 'Adjustment at %s makes new x = %s >= prior x = %s; eps = %s' % ( x12, x, x11, lowerEps ) )
                        else :
                            xys[-1] = [x, XYs1dModule.pointwiseXY_C.interpolatePoint(str(enumsModule.Interpolation.linlin), x, x11, y11, x12, y12)]
                    if( upperEps != 0. ) :
                        x22, y22 = _region[1]
                        x = getAdjustedX( x21, upperEps )
                        if( x >= x22 ) :
                            if( removeOverAdjustedPoints ) :
                                del _region[0]
                            else :
                                raise ValueError( 'Adjustment at %s makes new x = %s >= next x = %s; eps = %s' % ( x21, x, x22, upperEps ) )
                        else :
                            _region[0] = [x, XYs1dModule.pointwiseXY_C.interpolatePoint(str(enumsModule.Interpolation.linlin), x, x21, y21, x22, y22)]
            xys += _region
        pointwise = pointwiseClass( data = xys, axes = self.axes, outerDomainValue = self.outerDomainValue ) # FIXME - need more work to insure all parameters are set properly.
        return( pointwise )

    def toLinearXYsClass( self ) :
        """
        This method always returns the class :py:class:`XYs1dModule.XYs1d`.

        :returns:       The class :py:class:`XYs1dModule.XYs1d`.
        """

        return( XYs1dModule.XYs1d )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns a list of all 1d functions class that can be added as a region to an instance of :py:class:`Regions1d`.

        :returns:       A list of 1d functions classes
        """

        return( ( XYs1dModule.XYs1d, series1dModule.Series1d ) )

    @staticmethod
    def multiPlot(curve1ds, **kwargs):
        """
        Plots a list of 1d curves on the same plot. Uses each curve's 'plotLabel' as the legend key. Each curve must have a copyDataToXsAndYs method.
        """

        XYs1dModule.XYs1d.multiPlot(curve1ds, **kwargs)

class RegionsMultiD( Regions ) :
    """
    This is the absract base class for the other regions classes of dimension 2 and 3..
    """

    def getBoundingSubFunctions( self, domainValue ) :
        """
        This method returns the two sub-functions in the region containing the domain point *domainValue*
        whose outerDomainValues bound *value*. Additonal information is also returned.
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

        for region in self.regions :
            if( domainValue < region.domainMax ) : return( region.getBoundingSubFunctions( domainValue ) )

        return( self.regions[-1].getBoundingSubFunctions( domainValue ) )  # Domain value is above the last region. Let the last region determine what to do.

class Regions2d( RegionsMultiD ) :
    """
    This class supports storing a 2d function as a list of abutting 2d functions.
    """

    moniker = 'regions2d'
    dimension = 2

    def toLinearXYsClass( self ) :
        """
        This method always returns the class :py:class:`Regions`.

        :returns:       The class :py:class:`Regions`.
        """

        return( Regions )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns a list of all 2d functions class that can be added as a region to an instance of :py:class:`Regions2d`.

        :returns:       A list of 2d functions classes
        """

        from . import multiD_XYs as multiD_XYsModule

        return( ( multiD_XYsModule.XYs2d, ) )

class Regions3d( RegionsMultiD ) :
    """
    This class supports storing a 3d function as a list of abutting 3d functions.
    """

    moniker = 'regions3d'
    dimension = 3

    def toLinearXYsClass( self ) :
        """
        This method always returns the class :py:class:`Regions`.

        :returns:       The class :py:class:`Regions`.
        """

        return( Regions )

    @staticmethod
    def allowedSubElements( ) :
        """
        This method returns a list of all 3d functions class that can be added as a region to an instance of :py:class:`Regions3d`.

        :returns:       A list of 3d functions classes
        """

        from . import multiD_XYs as multiD_XYsModule

        return( ( multiD_XYsModule.XYs3d, ) )
