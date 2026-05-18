# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Some evaluations use short-range self-scaling variance components  (ENDF LB=8 or 9)
to represent part or all of the uncertainty.
When processing these sections to generate a multi-group covariance library,
the size of the variance depends either directly or inversely on the size of the processed group.
These sections only produce a diagonal, no cross-correlations.
"""

from . import enums as covarianceEnumsModule
from . import base

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

from xData import xDataArray as arrayModule
from xData import gridded as griddedModule
from xData import link as linkModule

from pqu import PQU


class DependenceOnProcessedGroupWidth(enumsModule.Enum):
    """Defines enums of processed width dependencies."""

    inverse = enumsModule.auto()
    direct = enumsModule.auto()


class ShortRangeSelfScalingVariance(ancestryModule.AncestryIO, base.Covariance):
    moniker = 'shortRangeSelfScalingVariance'

    def __init__(self, label, type=covarianceEnumsModule.Type.absolute,
                 dependenceOnProcessedGroupWidth=DependenceOnProcessedGroupWidth.inverse, matrix=None):

        ancestryModule.AncestryIO.__init__(self)
        self.__label = label
        self.dependenceOnProcessedGroupWidth = dependenceOnProcessedGroupWidth
        self.type = type
        self.matrix = matrix

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, value):
        if value is not None:
            if not isinstance(value, str): raise TypeError('label must be a string')
        self.__label = value

    def copy(self):
        srssv = ShortRangeSelfScalingVariance(self.label, self.type, self.dependenceOnProcessedGroupWidth, self.matrix)
        return srssv

    @property
    def dependenceOnProcessedGroupWidth(self):
        return self.__dependenceOnProcessedGroupWidth

    @dependenceOnProcessedGroupWidth.setter
    def dependenceOnProcessedGroupWidth(self, value):

        self.__dependenceOnProcessedGroupWidth = DependenceOnProcessedGroupWidth.checkEnumOrString(value)

    @property
    def domainUnit(self):
        return self.matrix.axes[-1].unit

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, type):

        self.__type = covarianceEnumsModule.Type.checkEnumOrString(type)

    @property
    def matrix(self):
        return self.__matrix

    @matrix.setter
    def matrix(self, matrix):
        if not isinstance(matrix, griddedModule.Gridded2d):
            raise TypeError("Expected Gridded2d instance, got '%s'" % type(matrix))
        self.__matrix = matrix
        self.__matrix.setAncestor(self)

    @property
    def gridded2d(self):
        return self.matrix  # convenience method

    def isSymmetric(self):
        return True

    def check(self, info):
        """
        Check that matrix is diagonal and all elements >= 0
        """
        from fudge import warning
        warnings = []

        if not isinstance(self.matrix.array, arrayModule.Diagonal):
            warnings.append(warning.InvalidShortRangeVarianceData(type(self.matrix.array)))
        eigenvals = self.matrix.array.constructArray().diagonal()
        if any(eigenvals < 0):
            warnings.append(warning.NegativeEigenvalues(negativeCount=(eigenvals < 0).sum(),
                                                        worstCase=eigenvals.min(), obj=self))
        return warnings

    def convertUnits(self, unitMap):

        self.matrix.convertUnits(unitMap)

    def rowBounds(self, unit=None):
        """Get the bounds of the row.  If unit is specified, return the bounds in that unit."""
        factor = 1
        if unit:
            factor = PQU.PQU(1, self.matrix.axes[-1].unit).getValueAs(unit)
        return (self.matrix.axes[-1].values[0] * factor, self.matrix.axes[-1].values[-1] * factor)

    def columnBounds(self, unit=None):
        """Get the bounds of the column.  If unit is specified, return the bounds in that unit."""
        if isinstance(self.matrix.axes[-2].values, linkModule.Link): return self.rowBounds(unit)
        factor = 1
        if unit:
            factor = PQU.PQU(1, self.matrix.axes[-1].unit).getValueAs(unit)
        return (self.matrix.axes[-1].values[0] * factor, self.matrix.axes[-1].values[-1] * factor)

    def toCovarianceMatrix(self, domain=None):
        raise NotImplementedError("Not yet implemented for ShortRangeSelfScalingVariance")

    def getUncertaintyVector(self, theData=None, relative=True):
        raise NotImplementedError("Not yet implemented for ShortRangeSelfScalingVariance")

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attrs = ''
        if self.label is not None: attrs = ' label="%s"' % self.label
        xmlString = ['%s<%s%s type="%s" dependenceOnProcessedGroupWidth="%s">' %
                     (indent, self.moniker, attrs, self.type, self.dependenceOnProcessedGroupWidth)]
        xmlString += self.matrix.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        label = element.get('label')
        if label is not None:
            xPath.append('%s[@label="%s"]' % (element.tag, label))
        else:
            xPath.append(element.tag)
        matrix_ = griddedModule.Gridded2d.parseNodeUsingClass(element.find(griddedModule.Gridded2d.moniker), xPath,
                                                              linkData, **kwargs)
        srssv = cls(label=label, type=element.get('type'),
                    dependenceOnProcessedGroupWidth=element.get('dependenceOnProcessedGroupWidth'), matrix=matrix_)
        xPath.pop()
        return srssv
