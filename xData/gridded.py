# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains gridded functions. Gridded functions are functions that store the dependent data on 
a regular grid for each indenpendent axis. Each independent axis can have a grid indenpendent of the other axis grids.
These are GNDS gridded1d, gridded2d and gridded3d functions.

This module contains the following classes:

    +---------------------------+-------------------------------------------------------------------------------+
    | Class                     | Description                                                                   |
    +===========================+===============================================================================+
    | Gridded                   | This class is the base class for the gridded classes.                         |
    +---------------------------+-------------------------------------------------------------------------------+
    | Gridded1d                 | This class represents a 1-d gridded function.                                 |
    +---------------------------+-------------------------------------------------------------------------------+
    | Gridded2d                 | This class represents a 2-d gridded function.                                 |
    +---------------------------+-------------------------------------------------------------------------------+
    | Gridded3d                 | This class represents a 3-d gridded function.                                 |
    +---------------------------+-------------------------------------------------------------------------------+
"""

from . import enums as enumsModule
from . import base as baseModule
from . import axes as axesModule
from . import xDataArray as arrayModule
from . import XYs1d as XYs1dModule

class Gridded( baseModule.XDataFunctional ) :
    """
    This class is the base class for the gridded classes.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | axes              | This is the axes member.                                      |
    +-------------------+---------------------------------------------------------------+
    | array             | This member stores the dependent data.                        |
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

    ancestryMembers = baseModule.XDataFunctional.ancestryMembers + ( 'array', )

    def __init__( self, axes, array, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None ) :
        """
        :param axes:                This is the axes member.
        :param array:               This stores the dependent data.
        :param index:               This is the index member.
        :param valueType:           This describes the type of data (i.e., float, int) returned by the function.
        :param outerDomainValue:    This is the domain value for the next higher dimension for a function that is
                                    embedded in a high dimensional function.
        :param label:               This is the label member.
        """

        for axis in axes[:self.dimension+1]:
            if( axis.index == 0 ) :
                if( not( isinstance( axis, axesModule.Axis ) ) ) : raise Exception( 'dependent axis must not have a grid' )
            else :
                if( not( isinstance( axis, axesModule.Grid ) ) ) : raise Exception( 'independent axis must have a grid' )
        if( not( isinstance( array, arrayModule.ArrayBase ) ) ) : raise TypeError( 'array not an array instance' )

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

        if( self.dimension != len( array ) ) : ValueError( 'self.dimension = %d != len( array ) = %d' % 
                ( self.dimension, len( array ) ) )

        self.array = array
        self.array.setAncestor( self )

    def __len__( self ) :
        """
        This method returns the number of data in *self*'s array.

        :returns:       A python int.
        """

        return( len( self.array ) )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if( self.axes is None ) : return
        factors = self.axes.convertUnits( unitMap )
        self.array.offsetScaleValues( 0.0, factors[0] )
        self.fixValuePerUnitChange( factors )

    def copy( self ):
        """
        This method returns a copy of self.

        :returns:           A new instance of *self*.
        """

        return self.__class__( self.axes.copy( ), self.array.copy( ), index = self.index, valueType = self.valueType,
            outerDomainValue = self.outerDomainValue, label = self.label )

    __copy__ = copy

    def getGrid( self, index ) :
        """
        This method returns the values for the grid for the *index* axis.

        :param index:       The index of the axis whose grid is requested.

        :returns            A python list of numbers.
        """

        if( 0 < index <= len( self ) ) : return( self.axes[index].values )
        raise IndexError( 'index = %d out of range [1,%s]' % (index, len( self ) ) )

    def getGridUnit( self, index ) :
        """
        This method returns the unit for the *index* axis.

        :param index:       The index of the axis whose unit is requested.

        :returns            A python str.
        """

        return( self.axes[index].unit )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if self.axes is not None: XML_strList += self.axes.toXML_strList(indent2, **kwargs)
        XML_strList += self.array.toXML_strList(indent2, **kwargs)
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

        :return:            An instance of *cls* representing *node*.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        axes = kwargs.get('axes', None)
        if node.find('axes') is not None: axes = axesModule.Axes.parseNodeUsingClass(node.find('axes'), xPath, linkData, **kwargs)

        array = arrayModule.ArrayBase.parseNodeUsingClass(node[1], xPath, linkData, **kwargs)

        instance = cls(axes, array, **attributes)

        xPath.pop()

        return instance

class Gridded1d( Gridded ) :
    """
    This class stores a 1-d gridded function.
    """

    moniker = 'gridded1d'
    dimension = 1

    def asXYs1d(self, asLinLin, accuracy, lowerEps, upperEps, biSectionMax=16):
        """
        This method returns a representation of the data in *self* as an :py:class:`XYs1dModule.XYs1d` instance.

        :param asLinLin:    If **True**, the returned data will have lin-lin interpolation.
        :param accuracy:    Used to determine the accuracy if converting data to lin-lin interpolated data.
        :param lowerEps:    Used to dull the lower point for "flat" interpolation.
        :param upperEps:    Used to dull the upper point for "flat" interpolation.

        :returns:           A :py:class:`XYs1dModule.XYs1d` instance.
        """

        xs, ys = self.copyDataToXsAndYs()
        cls = self.toLinearXYsClass()
        xys1d = cls(data=[xs, ys], dataForm='xsandys', axes=self.axes, interpolation=enumsModule.Interpolation.flat)
        if asLinLin:
            xys1d = xys1d.changeInterpolation(enumsModule.Interpolation.linlin, accuracy=accuracy, lowerEps=lowerEps, upperEps=upperEps)

        return xys1d

    def constructVector(self):
        """
        Generate an xData.vector.Vector object from an xData.gridded.gridded1d instance.

        :returns: instance of xData.vector.Vector.
        """
        return self.array.constructVector()

    def copyDataToXsAndYs(self):
        """
        """

        xs = self.axes[1].values.values
        ys = list(self.constructVector())
        ys.append(ys[-1])

        return xs, ys

class Gridded2d( Gridded ) :
    """
    This class stores a 2-d gridded function.
    """

    moniker = 'gridded2d'
    dimension = 2

class Gridded3d( Gridded ) :
    """
    This class stores a 3-d gridded function.
    """

    moniker = 'gridded3d'
    dimension = 3

    def constructMatrix(self, thirdIndex):
        """
        Generate an xData.matrix.Matrix object from an xData.gridded.gridded3d instance.

        :param thirdIndex: Index along the third dimension from which to extract 2d-array.
        :returns: instance of xData.matrix.Matrix.
        """
        return self.array.constructMatrix(thirdIndex)
