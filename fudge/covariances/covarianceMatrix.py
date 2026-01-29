# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Base classes for covariances: matrix, axes."""

import copy
import numpy

from . import enums as covarianceEnumsModule
from . import base

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import gridded as griddedModule
from xData import xDataArray as arrayModule
from xData import XYs1d as XYs1dModule
from xData import link as linkModule
from xData import values as valuesModule

from pqu import PQU

lowerEps = 1e-8
upperEps = 1e-8


class CovarianceMatrix(ancestryModule.AncestryIO, base.Covariance):
    """
    Simplest form of covariance. covarianceMatrix contains a label, a covariance 'type' which might be 'absolute',
    'relative' or 'correlation', and an optional 'productFrame' (for covariances on outgoing distributions).
    Matrix data is stored in an :py:class:`xData.gridded.Gridded2d` class. May be diagonal, symmetric, sparse, etc
    """

    moniker = 'covarianceMatrix'
    keyName = 'label'

    def __init__(self, label, type=covarianceEnumsModule.Type.absolute, matrix=None, productFrame=None):
        """

        :param label: refers either to a style (i.e., 'eval'), or to the index inside a 'mixed' covariance
        :param type: 'relative', 'absolute' or 'correlation'
        :param productFrame: for outgoing distributions, indicates whether the covariance applies to lab or COM
        :type productFrame: str
        :param matrix: the covariance matrix itself
        :type matrix: :py:class:`xData.gridded.Gridded2d`
        :return:
        :rtype: CovarianceMatrix
        """

        ancestryModule.AncestryIO.__init__(self)
        self.__label = label
        self.__type = covarianceEnumsModule.Type.checkEnumOrString(type)
        self.productFrame = productFrame
        self.matrix = matrix  #: a :py:class:`xData.gridded.Gridded2d` instance containing the matrix
        self.matrix.setAncestor(self)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, value):

        if (value is not None):
            if (not (isinstance(value, str))): raise TypeError('label must be a string')
        self.__label = value

    @property
    def domainUnit(self):
        return self.matrix.axes[-1].unit

    @property
    def productFrame(self):
        return self.__productFrame

    @productFrame.setter
    def productFrame(self, value):

        if value is None:
            self.__productFrame = None
        else:
            self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(value)

    @property
    def type(self):
        return self.__type

    @staticmethod
    def constructFromData(energyUnit, rowGrid, columnGrid, matrix,
                          label="", type=covarianceEnumsModule.Type.relative, matrixUnit=""):
        """Construct a CovarianceMatrix from raw data.

        :param energyUnit: unit for the energy grids
        :param rowGrid: grid for the rows (as a list of values)
        :param columnGrid: grid for the columns (as a list of values). Use None if same as rowGrid
        :param matrix: 2d list or numpy array containing the matrix values
        :param label: label (e.g. 'eval')
        :param type: type of the covariance matrix (relative or absolute)
        :param matrixUnit: unit for covariance matrix elements
        :return: CovarianceMatrix instance
        """
        assert len(rowGrid) == len(matrix) + 1, "rowGrid length must be matrix size + 1"
        if columnGrid is not None:
            assert len(columnGrid) == matrix.shape[1] + 1, "columnGrid length must be matrix size + 1"

        axes = axesModule.Axes(3, labelsUnits={
            0: ('matrix_elements', matrixUnit),
            1: ('column_energy_bounds', energyUnit),
            2: ('row_energy_bounds', energyUnit)})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                style=xDataEnumsModule.GridStyle.boundaries, values=valuesModule.Values(rowGrid))
        if columnGrid is not None:
            axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                    style=xDataEnumsModule.GridStyle.boundaries,
                                    values=valuesModule.Values(columnGrid))
        else:
            axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                    style=xDataEnumsModule.GridStyle.boundaries,
                                    values=linkModule.Link(link=axes[2].values, relative=True))

        matrix = numpy.array(matrix)
        if columnGrid is not None:
            array = arrayModule.Full(shape=matrix.shape,
                                    data=matrix.flatten().tolist(),
                                    symmetry=arrayModule.Symmetry.none)
        else:
            # matrix is symmetric. Only store lower triangle
            triData = matrix[numpy.tril_indices(len(matrix))].tolist()
            array = arrayModule.Full(shape=matrix.shape,
                                    data=triData,
                                    symmetry=arrayModule.Symmetry.lower)

        gridded2d = griddedModule.Gridded2d(axes=axes, array=array)

        return CovarianceMatrix(label=label, type=type, matrix=gridded2d)

    def getValue(self, x, y):
        ix = self.matrix.axes[2].getIndexOfValue(x)
        if isinstance(self.matrix.axes[1].values, linkModule.Link):
            iy = self.matrix.axes[2].getIndexOfValue(y)
        else:
            iy = self.matrix.axes[1].getIndexOfValue(y)
        return self.matrix.array.constructArray()[ix, iy]

    def rowBounds( self, unit=None ):
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

    def isSymmetric(self):
        """
        Simple test to determine if an underlying matrix is symmetric
        :return: boolean
        """

        if (self.matrix.array.compression == arrayModule.Diagonal.moniker or
                self.matrix.array.symmetry in (arrayModule.Symmetry.lower, arrayModule.Symmetry.upper)):
            return True
        # could still be symmetric even if it doesn't use compression
        arr = self.matrix.array.constructArray()
        return numpy.all(arr == arr.T)

    def convertAxesToUnits(self, units):
        """
        Converts all the axes' units.
        The parameter ``units`` should be a list of units with the same length as self.axes
        """
        if not type(units) in [list, tuple]:
            raise TypeError("units argument must be a list or tuple of strings")
        if len(units) != len(self.matrix.axes):
            raise ValueError("requested units list has a different length than the number of axes")
        for i, a in enumerate(self.matrix.axes):
            if isinstance(a, axesModule.Grid):
                a.convertToUnit(units[i])
            elif isinstance(a, axesModule.Axis):
                a.unit = units[i]
            else:
                raise TypeError()

    def convertUnits(self, unitMap):

        self.matrix.convertUnits(unitMap)

    def domainSlice(self, rowDomainBounds, columnDomainBounds=None, label=""):
        """
        Return copy with a revised domain for rows and columns.
        If requested domain extends beyond limits of self, add rows/columns of zeros to the new matrix.

        :param rowDomainBounds: tuple (domainMin, domainMax)
        :param columnDomainBounds: tuple (domainMin, domainMax). If not supplied, columnDomainBounds = rowDomainBounds
        :param label: string to label the result
        :return: CovarianceMatrix
        """
        import bisect

        if columnDomainBounds is None:
            columnDomainBounds = rowDomainBounds

        if rowDomainBounds == self.rowBounds() and columnDomainBounds == self.columnBounds():
            return self.copy()

        def getNewBinBoundaries(bins, rowDomainBounds):
            """Get revised bin boundaries between domain min/max. Also determine whether to pad matrix with zeros"""
            matrixIndices = []
            indices = []
            for index, boundary in enumerate(rowDomainBounds):
                if boundary in bins:
                    indices.append(bins.index(boundary))
                    matrixIndices.append(bins.index(boundary))
                elif boundary > bins[-1]:
                    matrixIndices.append(len(bins) - 1)
                    indices.append(bisect.bisect(bins, boundary))
                else:
                    if index == 0:
                        matrixIndices.append(max(bisect.bisect(bins, boundary) - 1, 0))
                    else:
                        matrixIndices.append(bisect.bisect(bins, boundary))
                    indices.append(bisect.bisect(bins, boundary))

            new_bins = bins[slice(*indices)]
            slice1 = slice(*matrixIndices)
            if not new_bins:
                new_bins = [rowDomainBounds[0], rowDomainBounds[1]]
            if new_bins[0] != rowDomainBounds[0]:
                new_bins.insert(0, rowDomainBounds[0])
            if new_bins[-1] != rowDomainBounds[1]:
                new_bins.append(rowDomainBounds[1])

            indices2 = [0, matrixIndices[1] - matrixIndices[0]]
            if rowDomainBounds[0] < bins[0]:  # pad start of matrix with extra zeros
                indices2[0] += 1
                indices2[1] += 1
            slice2 = slice(*indices2)

            return slice1, slice2, new_bins

        newAxes = self.matrix.axes.copy()

        row_ebounds = list(self.matrix.axes[-1].values)
        s1, s1a, new_row_ebounds = getNewBinBoundaries(row_ebounds, rowDomainBounds)
        newAxes[2] = axesModule.Grid(newAxes[2].label, newAxes[2].index, newAxes[2].unit,
                                     style=xDataEnumsModule.GridStyle.boundaries,
                                     values=valuesModule.Values(new_row_ebounds))

        if (isinstance(self.matrix.axes[-2].values, linkModule.Link) and
                self.matrix.axes[-2].values.link is self.matrix.axes[-1].values):
            new_col_ebounds = new_row_ebounds[:]
            s2 = s1
            s2a = s1a
            newAxes[1] = axesModule.Grid(newAxes[1].label, newAxes[1].index, newAxes[1].unit,
                                         style=newAxes[1].style,
                                         values=linkModule.Link(link=newAxes[2].values, relative=True))
        else:
            col_ebounds = list(self.matrix.axes[-2].values)
            s2, s2a, new_col_ebounds = getNewBinBoundaries(col_ebounds, columnDomainBounds)
            newAxes[1] = axesModule.Grid(newAxes[1].label, newAxes[1].index, newAxes[1].unit,
                                         style=xDataEnumsModule.GridStyle.boundaries,
                                         values=valuesModule.Values(new_col_ebounds))

        rawMatrix = numpy.zeros(shape=(len(new_row_ebounds) - 1, len(new_col_ebounds) - 1))
        rawMatrix[s1a, s2a] = self.matrix.array.constructArray()[s1, s2]

        # Pack result inside a covarianceMatrix/gridded2d
        gridded = griddedModule.Gridded2d(axes=newAxes,
                                          array=arrayModule.Full(rawMatrix.shape, data=rawMatrix.flatten()))

        return CovarianceMatrix(label, type=self.type, matrix=gridded, productFrame=self.productFrame)

    def toCovarianceMatrix(self):
        """
        Return copy of self (for consistency with toCovarianceMatrix behavior in summed.py / mixed.py)
        :return:
        """
        cm = self.copy()
        cm.setAncestor(self.ancestor)
        return cm

    def getCorrelationMatrix(self):
        """
        Returns the correlation matrix generated from self's covariance matrix.  This is
        essentially a copy of self, but renormalized by the uncertainty:

            correlation[i,j] = covariance[i,j]/sqrt(covariance[i,i])/sqrt(covariance[j,j])

        We reuse the covariance matrix class so that we can do plotting, etc.  If you have
        a correlation matrix, you can safely recover it provided you have the uncertainty
        vector.

        Currently only works for a square covariance matrix and not a off-diagonal part of
        another covariance.
        """
        # Check if is full, square covariance matrix
        if not self.isSymmetric():
            raise TypeError("Can only extract correlation matrices from symmetric covariance matrices")

        # Rescale result matrix
        theCorrelationMatrix = self.matrix.array.constructArray()
        theUncertainty = copy.copy(numpy.diag(theCorrelationMatrix))
        theUncertainty[theUncertainty < 0.0] = 0.0
        theUncertainty = numpy.sqrt(theUncertainty)
        with numpy.errstate(divide='ignore', invalid='ignore'):
            theCorrelationMatrix /= theUncertainty
            theCorrelationMatrix /= theUncertainty[:,numpy.newaxis]

        # Return the result
        tridata = theCorrelationMatrix[numpy.tri(theCorrelationMatrix.shape[0], dtype=bool)].tolist()
        correlation = griddedModule.Gridded2d(
            axes=self.matrix.axes.copy(),  # FIXME: unresolvedLinks still unresolved!
            array=arrayModule.Full(shape=theCorrelationMatrix.shape, data=tridata,
                                   symmetry=arrayModule.Symmetry.lower))
        correlation.axes[0].unit = ''
        return correlation

    def toAbsolute(self, rowData=None, colData=None):
        """
        Rescales self (if it is a relative covariance) using XYs1d rowData and colData
        to convert self into an absolute covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite

        .. note:    If the column axis is a link, only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already
                    absolute because this will throw an error.

        :returns: a copy of self, but rescaled and with the type set to absolute
        """
        if self.type == covarianceEnumsModule.Type.absolute:
            return copy.copy(self)

        if hasattr(self, '_absolute'):
            return self._absolute

        # Make sure we have usable row data to rescale with
        if rowData is None:
            rowData = self.findAttributeInAncestry('rowData').link.toPointwise_withLinearXYs(lowerEps=lowerEps,
                                                                                             upperEps=upperEps)
        if not isinstance(rowData, XYs1dModule.XYs1d):
            raise TypeError('rowData must be of type XYs1d, found %s' % type(rowData))
        gRowData = rowData.group(self.matrix.axes[2].values, norm='dx')

        # Only generate the column rescaling if we need to
        if not (self.isSymmetric() and isinstance(self.matrix.axes[1].values, linkModule.Link)):
            if colData is None:
                colData = self.findAttributeInAncestry('columnData').link.toPointwise_withLinearXYs(lowerEps=lowerEps,
                                                                                                    upperEps=upperEps)
            if not isinstance(colData, XYs1dModule.XYs1d):
                raise TypeError('colData must be of type XYs1d, found %s' % type(colData))

            if isinstance(self.matrix.axes[1].values, linkModule.Link):
                columnGroups = self.matrix.axes[1].values.link.values
            else:
                columnGroups = self.matrix.axes[1].values
            gColData = colData.group(columnGroups, norm='dx')
        else:
            colData = rowData
            gColData = gRowData

        new_data = self.matrix.array.constructArray() * numpy.outer(gRowData, gColData)
        if self.matrix.array.symmetry == arrayModule.Symmetry.lower:
            new_data = new_data[numpy.tril_indices(len(new_data))]
        elif self.matrix.array.symmetry == arrayModule.Symmetry.upper:
            new_data = new_data[numpy.triu_indices(len(new_data))]
        new_array = arrayModule.Full(shape=self.matrix.array.shape,
                                     data=new_data.flatten(),
                                     symmetry=self.matrix.array.symmetry,
                                     storageOrder=self.matrix.array.storageOrder,
                                     label=self.matrix.array.label
                                     )

        newGridded2d = griddedModule.Gridded2d(axes=self.matrix.axes.copy(),
                                               array=new_array,
                                               label=self.matrix.label)

        # Set the final units
        if rowData.axes[0].unit == colData.axes[0].unit:
            if rowData.axes[0].unit != '':
                newGridded2d.axes[0].unit = rowData.axes[0].unit + '**2'
        else:
            newGridded2d.axes[0].unit = rowData.axes[0].unit + '*' + colData.axes[0].unit

        result = CovarianceMatrix(label="toAbsolute", type=covarianceEnumsModule.Type.absolute,
                                  matrix=newGridded2d, productFrame=self.productFrame)

        result.setAncestor(self.ancestor)
        result._relative = self

        return result

    def toRelative(self, rowData=None, colData=None):
        """
        Rescales self (if it is a absolute covariance) using XYs1d rowData and colData
        to convert self into a relative covariance matrix.

        :param rowData: an XYs1d instance containing data to rescale covariance in the "row direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite
        :param colData: an XYs1d instance containing data to rescale covariance in the "col direction"
                            if it isn't given, we'll compute it from the corresponding data in the reactionSuite

        .. note::   If the column axis is a link, only rowData is needed.
                    If neither rowData nor colData are specified, you'd better hope that the covariance is already
                    relative because this will throw an error.

        :returns: a copy of self, but rescaled and with the type set to relative
        """
        if self.type == covarianceEnumsModule.Type.relative:
            return copy.copy(self)

        if hasattr(self, '_relative'):
            return self._relative

        # Make sure we have usable row data to rescale with
        if rowData is None:
            rowData = self.findAttributeInAncestry('rowData').link.toPointwise_withLinearXYs(lowerEps=lowerEps,
                                                                                             upperEps=upperEps)
        if not isinstance(rowData, XYs1dModule.XYs1d):
            raise TypeError('rowData must be of type XYs1d, found %s' % type(rowData))
        gRowData = rowData.group(self.matrix.axes[2].values, norm='dx')

        # Only generate the column rescaling if we need to
        if not (self.isSymmetric() and isinstance(self.matrix.axes[1].values, linkModule.Link)):
            if colData is None:
                colData = self.findAttributeInAncestry('columnData').link.toPointwise_withLinearXYs(lowerEps=lowerEps,
                                                                                                    upperEps=upperEps)
            if not isinstance(colData, XYs1dModule.XYs1d):
                raise TypeError('colData must be of type XYs1d, found %s' % type(colData))

            if isinstance(self.matrix.axes[1].values, linkModule.Link):
                columnGroups = self.matrix.axes[1].values.link.values
            else:
                columnGroups = self.matrix.axes[1].values
            gColData = colData.group(columnGroups, norm='dx')
        else:
            gColData = gRowData

        denom = numpy.outer(gRowData, gColData)
        denom[denom == 0] = lowerEps  # remove zeros
        new_data = self.matrix.array.constructArray() / denom
        if self.matrix.array.symmetry == arrayModule.Symmetry.lower:
            new_data = new_data[numpy.tril_indices(len(new_data))]
        elif self.matrix.array.symmetry == arrayModule.Symmetry.upper:
            new_data = new_data[numpy.triu_indices(len(new_data))]
        new_array = arrayModule.Full(shape=self.matrix.array.shape,
                                     data=new_data.flatten(),
                                     symmetry=self.matrix.array.symmetry,
                                     storageOrder=self.matrix.array.storageOrder,
                                     label=self.matrix.array.label
                                     )

        newGridded2d = griddedModule.Gridded2d(axes=self.matrix.axes.copy(),
                                               array=new_array,
                                               label=self.matrix.label)
        # Set the final units
        newGridded2d.axes[0].unit = ''

        result = CovarianceMatrix(label="toRelative", type=covarianceEnumsModule.Type.relative,
                                  matrix=newGridded2d, productFrame=self.productFrame)

        result.setAncestor(self.ancestor)
        result._absolute = self

        return result

    def copy(self):

        return CovarianceMatrix(self.label, self.type, self.matrix.copy(), self.productFrame)

    def check(self, info):
        """Check if relative uncertainty is in the bounds passed into the checker,
        and checks for negative eigenvalues or other issues."""

        from fudge import warning
        warnings = []

        rowGrid = list(self.matrix.axes[2].values)
        try:
            columnData = self.findAttributeInAncestry('columnData')
            columnGrid = list(self.matrix.axes[1].values)
        except:
            columnData = None
            columnGrid = rowGrid

        if rowGrid != sorted(set(rowGrid)):
            warnings.append(warning.NonMonotonicGrid("rows", obj=self))
        if columnData and columnGrid != sorted(set(columnGrid)):
            warnings.append(warning.NonMonotonicGrid("columns", obj=self))

        A = self.matrix.array.constructArray()
        if A.shape != (len(rowGrid)-1, len(columnGrid)-1):
            warnings.append(warning.MatrixDimensionMismatch(
                A.shape, (len(rowGrid), len(columnGrid)), obj=self))

        if self.findAttributeInAncestry('crossTerm') or not info['checkUncLimits']:
            # remaining checks don't apply to cross-terms yet
            # FIXME also test FULL matrix including cross-terms
            return warnings

        # run variance and eigenvalue checks on original matrix, relative or absolute:
        variance = A.diagonal()
        if any(variance < 0):
            for badIndex in numpy.argwhere(variance < 0).flatten():
                badValue = variance[badIndex]
                warnings.append(warning.NegativeVariance(badIndex, badValue, obj=self))

            # skip additional tests if negatives variances encountered
            return warnings

        vals, vecs = numpy.linalg.eigh(A)
        if min(vals) < info['negativeEigenTolerance']:
            warnings.append(warning.NegativeEigenvalues(len(vals[vals < 0]), min(vals), self))

        # Check that the condition number of the matrix is reasonable
        minpos, maxpos = min(vals[vals >= 0]), max(vals[vals >= 0])
        if A.size != 1 and minpos / maxpos < info['eigenvalueRatioTolerance']:
            warnings.append(warning.BadEigenvalueRatio(minpos / maxpos, self))

        if self.type is not covarianceEnumsModule.Type.relative:
            # FIXME next line skips testing PFNS covariances
            rowData = self.findAttributeInAncestry('rowData')
            if not hasattr(rowData.link, 'toPointwise_withLinearXYs'):
                return warnings

        relativeMatrix = self.toRelative()

        Arr = relativeMatrix.matrix.array.constructArray()
        relativeUncertainty = numpy.sqrt(Arr.diagonal())

        minUnc = info['minRelUnc']
        maxUnc = info['maxRelUnc']
        for idx, uncert in enumerate(relativeUncertainty):
            if minUnc <= uncert <= maxUnc:
                pass  # uncertainty where it should be
            elif uncert < minUnc:
                warnings.append(warning.VarianceTooSmall(
                    idx, relativeUncertainty=uncert, absoluteUncertainty=None, obj=self))
            else:
                warnings.append(warning.VarianceTooLarge(
                    idx, relativeUncertainty=uncert, absoluteUncertainty=None, obj=self))

        return warnings

    def fix(self, **kw):
        """Fix uncertainty using the bounds passed into the fixer.
        Requires specification of the data ("theData") if the covariance is not relative.
        I was not creative when I coded this, so it will fail when theData.getValue( x )
        doesn't exist or is a function of more than one value. """

        warnings = []
        if self.isSymmetric() and kw['fixUncLimits']:
            A = numpy.array(self.matrix.data)
            relative = self.type == covarianceEnumsModule.Type.relative
            for idx in range(A.shape[0]):
                eMin = self.axes[2].grid[idx]
                eMax = self.axes[2].grid[idx + 1]
                if relative:
                    uncMin = kw['minRelUnc']
                    uncMax = kw['maxRelUnc']
                else:
                    uncMin = kw['theData'].getValue(0.5 * (eMax + eMin))
                    uncMax = kw['theData'].getValue(0.5 * (eMax + eMin))
                uncMin2 = uncMin * uncMin
                uncMax2 = uncMax * uncMax
                # eThresh = threshold.getValueAs( component.axes[0].units )
                if uncMin2 <= A[idx, idx] <= uncMax2:
                    pass  # unc is where it should be
                elif uncMin2 >= A[idx, idx]:
                    if idx + 1 < A.shape[0] and A[idx + 1, idx + 1] >= uncMin2:
                        A[idx, idx] = A[idx + 1, idx + 1]
                    else:
                        A[idx, idx] = uncMin2
                # else:                                            # above threshold and uncertainty out of bounds
                    # continue #skip this fix for now
                    # if uncMin2 > A[idx,idx]: # uncertainty too small
                        # print('    WARNING: bin', idx, 'uncertainty is too small!!!', '(', uncMin2, '>', A[idx,idx], ')')
                        # if A[idx, idx] == 0.0:
                            # for jdx in range( idx+1, A.shape[0] ):
                                # A[idx,jdx] = 0.0
                                # A[jdx,idx] = 0.0
                            # A[idx,idx] == uncMin2
                        # else:
                            # for jdx in range( idx+1, A.shape[0] ):
                                # A[idx,jdx] *= uncMin / math.sqrt( A[idx,idx] )
                                # A[jdx,idx] = A[idx,jdx]
                            # A[idx,idx] == uncMin2
                    # else:                # uncertainty too big
                        # print('    WARNING: bin', idx, 'uncertainty is too big!!!', '(', A[idx,idx], '>', uncMax2, ')')
                        # for jdx in range( idx+1, A.shape[0] ):
                            # A[idx,jdx] *= uncMax / math.sqrt( A[idx,idx] )
                            # A[jdx,idx] = A[idx,jdx]
                        # A[idx,idx] = uncMax2
            self.data = A.tolist()
        return warnings + self.matrix.fix(**kw)

    def group(self, groupBoundaries=(None, None), groupUnit=(None, None)):
        r"""
        Group the matrix in self

        :param groupBoundaries: a 2 element list containing the group boundaries for the rows
                                and columns (respectively) of the covariance to be regrouped
                                rows go in the first element, columns in the second
        :param groupUnit: a 2 element list containing the units in which group boundaries are
                          specified for the rows and columns (respectively) of the covariance
                          to be regrouped

        :returns: the regrouped matrix (an xData.array.full as the array in a gridded2d.matrix)

        .. note::  We still need to do flux weighting


        .. rubric:: Regrouping Theory

        Given a function :math:`f(E)`, we write the grouped data using fudge's ``flat`` interpolation
        scheme.  We note that we could write this scheme as an expansion over basis functions:

        .. math::
            f(E) = \sum_{i=0}^{N+1} w_i(E) * f_i

        where the weight functions :math:`w_i(E)` are

        .. math::
            w_i(E) = 1  \;\\text{for}\; E_i <= E <= E_{i+1}; \;\; 0 \;\\textrm{otherwise}

        These weights are an orthogonal (but not orthonormal) basis, with

        .. math::
            (E_{i+1}-E_i) \delta_{ij} = \int dE w_i(E) * w_j(E)

        So, to transform from basis :math:`w_i(E)` to :math:`v_i(E)` (which has group boundaries
        :math:`[ E'_0, ... ]`), do:

        .. math::
            f'_j = \sum_i m_{ji} f_i

        where :math:`f'` is the regrouped function coefficients and :math:`m_{ji}` is the matrix

        .. math::
            m_{ij} = (E'_{i+1}-E'_i)^{-1} \int dE v_i(E) w_j(E)


        .. rubric:: Applying regrouping theory to covariance matrices

        When we are given a covariance matrix :math:`c_{ij}` in ENDF, it is meant to be interpreted
        as a grouped covariance in both the direction of the matrix rows and the matrix
        columns.  Therefore, we must regroup in both the row direction and the column
        direction.  The ENDF format gives both the group boundaries for the rows and columns.
        In other words, ENDF gives us the following rule for evaluating the continuous row-
        column covariance:

        .. math::
            c( E1, E2 ) = \sum_{ij} w_i(E1) w_j(E2) c_{ij}

        Computing :math:`m_{ij}` as before,

        .. math::
            cc_{ij} = \sum_{i',j'} m_{ii'} c_{i'j'} m_{j'j}

        It is straightforward to generalize to the case where the row and column bases are
        different.

        In the routine below, we abuse :py:class:`xData.XYs1d` to specify the functions
        :math:`w_i(E)` and use the :py:func:`XYs1d.groupOneFunction()` method to perform the integrals to get
        the regrouping matrix.  We do this separately for the rows and the columns.
        The matrix multiplication that converts a covariance from one pair of bases (group
        structures) to another is accomplished using numpy.


        .. rubric:: An explanation of fudge's 'flat' interpolation

        Suppose we have a function :math:`f(E)` specified using fudge's `'flat'` interpolation.
        Then we have :math:`N` entries :math:`[f_0, f_1, ..., f_{N-1}]` and a set of group
        boundaries :math:`[E_0, E_1, ..., E_N]` and the following rule for interpolation:

            * Below :math:`E_0`, :math:`f(E)` evaluates to :math:`0.0`
            * From :math:`E_0 \\rightarrow E_1`, :math:`f(E)` evaluates to :math:`f_0`
            * From :math:`E_1 \\rightarrow E_2`, :math:`f(E)` evaluates to :math:`f_1`
            * ...
            * From :math:`E_{i} \\rightarrow E_{i+1}`, :math:`f(E)` evaluates to :math:`f_i`
            * ...
            * From :math:`E_{N-1} \\rightarrow E_N`, :math:`f(E)` evaluates to :math:`f_{N-1}`
            * Above :math:`E_N`, :math:`f(E)` evaluates to :math:`0.0`
        """
        # determine where to get the settings for the potentially mirrored second axis
        if isinstance(self.matrix.axes[1].values, linkModule.Link):
            axis1index = 2
        else:
            axis1index = 1

        # setup the old axes in a form we can (ab)use in the XYs1d class
        axes2_ = axesModule.Axes(2, labelsUnits={
            1: (self.matrix.axes[2].label, self.matrix.axes[2].unit),
            0: ('dummy', '')})
        axes1_ = axesModule.Axes(2, labelsUnits={
            1: (self.matrix.axes[axis1index].label, self.matrix.axes[axis1index].unit),
            0: ('dummy', '')})

        # define basis functions for the rows and columns
        basis2 = XYs1dModule.XYs1d(axes=axes2_,
                                   data=[(x, 0.0) for x in self.matrix.axes[2].values],
                                   interpolation=xDataEnumsModule.Interpolation.flat)
        basis1 = XYs1dModule.XYs1d(axes=axes1_,
                                   data=[(x, 0.0) for x in self.matrix.axes[axis1index].values],
                                   interpolation=xDataEnumsModule.Interpolation.flat)
        basis2 = basis2.convertAxisToUnit(1, groupUnit[0])
        basis1 = basis1.convertAxisToUnit(1, groupUnit[1])

        # build the regrouping matrices for the two bases
        w0 = []
        for idx in range(self.matrix.array.shape[0]):
            basis2[idx] = (basis2[idx][0], 1.0)
            w0.append(basis2.groupOneFunction(groupBoundaries[0], norm='dx'))
            basis2[idx] = (basis2[idx][0], 0.0)
        w0 = numpy.array(w0)
        w1 = []
        for j in range(self.matrix.array.shape[1]):
            basis1[j] = (basis1[j][0], 1.0)
            w1.append(basis1.groupOneFunction(groupBoundaries[1], norm='dx'))
            basis1[j] = (basis1[j][0], 0.0)
        w1 = numpy.array(w1)

        # set up the regrouped covariance matrix
        grouped = self.copy()
        grouped.matrix.axes[2].values.values = groupBoundaries[0]
        if groupBoundaries[1] != groupBoundaries[0]:
            grouped.matrix.axes[1].values.values = groupBoundaries[1]
        grouped.matrix.axes[2].unit = groupUnit[0]
        grouped.matrix.axes[1].unit = groupUnit[1]
        odata = self.matrix.array.constructArray()
        gdata = numpy.dot(w0.T, numpy.dot(odata, w1))
        if self.isSymmetric() and groupBoundaries[0] == groupBoundaries[1]:
            trigdata = gdata[numpy.tril_indices(gdata.shape[0])]
            grouped.matrix.array = arrayModule.Full(shape=gdata.shape, data=trigdata,
                                                    symmetry=arrayModule.Symmetry.lower)
        else:
            grouped.matrix.array = arrayModule.Full(shape=gdata.shape, data=gdata.flatten().tolist())
        return grouped

    def removeExtraZeros(self, verbose=False):
        """
        Remove all extra zeros from the underlying matrix
        :return:
        """
        theMatrix = self.matrix.array.constructArray()
        rowStart, colStart = 0, 0
        rowEnd, colEnd = theMatrix.shape

        if verbose: print('before', theMatrix)

        # Figure out which end rows/columns are full of zeros.  We can remove those
        while numpy.all(theMatrix[rowStart, :] == 0): rowStart += 1
        while numpy.all(theMatrix[:, colStart] == 0): colStart += 1
        while numpy.all(theMatrix[rowEnd - 1, :] == 0): rowEnd -= 1
        while numpy.all(theMatrix[:, colEnd - 1] == 0): colEnd -= 1

        theMatrix = theMatrix[rowStart:rowEnd, colStart:colEnd]

        if verbose: print('after', theMatrix)

        if verbose: print('before', self.matrix.axes[-1].toXML(), self.matrix.axes[-2].toXML())

        if isinstance(self.matrix.axes[-2].values, linkModule.Link):
            assert (rowStart, rowEnd) == (colStart, colEnd)
            self.matrix.axes[-1].values.values = self.matrix.axes[-1].values[rowStart:rowEnd + 1]
        else:
            self.matrix.axes[-1].values.values = self.matrix.axes[-1].values[rowStart:rowEnd + 1]
            self.matrix.axes[-2].values.values = self.matrix.axes[-2].values[colStart:colEnd + 1]

        if verbose: print('after', self.matrix.axes[-1].toXML(), self.matrix.axes[-2].toXML())

        self.matrix.array = arrayModule.Full(shape=theMatrix.shape,
                                             data=theMatrix[numpy.tri(theMatrix.shape[0]) == 1.0].tolist(),
                                             symmetry=arrayModule.Symmetry.lower)

    def getUncertaintyVector(self, theData=None, relative=True):
        """ 
        Get an XYs1d object containing uncertainty for this matrix.
        Convert relative/absolute if requested (if so, must also pass central values as theData)

        Examples:

            - if the covariance matrix is relative and we want relative uncertainty vector, just do:
            
                > matrix.getUncertaintyVector()
                
            - if we want the absolute matrix instead:
            
                > matrix.getUncertaintyVector( theData=<XYs1d instance>, relative=False )
                
        """
        if not self.isSymmetric():
            raise ValueError("getUncertaintyVector only applies to symmetric matrices!")
        energies = list(self.matrix.axes[-1].values)
        diag = numpy.diagonal(self.matrix.array.constructArray()).copy()
        diag[diag < 0.0] = 0.0
        diag = list(numpy.sqrt(diag))
        diag.append(diag[-1])  # point corresponding to final energy bin
        yunit = self.matrix.axes[0].unit
        if yunit != '':  # get square root of the unit
            yunit = PQU.PQU(1, yunit).sqrt().getUnitSymbol()
        axes_ = axesModule.Axes(2, labelsUnits={1: ('energy_in', self.matrix.axes[2].unit), 0: ('uncertainty', yunit)})
        uncert = XYs1dModule.XYs1d(list(zip(energies, copy.deepcopy(diag))), axes=axes_,
                                   interpolation=xDataEnumsModule.Interpolation.flat)
        uncert = uncert.changeInterpolation(xDataEnumsModule.Interpolation.linlin,
                                            accuracy=1e-3, lowerEps=1e-8, upperEps=1e-8)

        # do we need to convert absolute->relative or vice versa?
        if (relative and self.type == covarianceEnumsModule.Type.absolute) or (
                not relative and self.type == covarianceEnumsModule.Type.relative):
            if theData is None:
                theData = self.findAttributeInAncestry('rowData').link.toPointwise_withLinearXYs(
                    lowerEps=1e-8, upperEps=1e-8)
            try:
                if theData.domainUnit != uncert.domainUnit:
                    theData = theData.copy()
                    theData.convertUnits({theData.domainUnit: uncert.domainUnit})
                theData = theData.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
                uncert, theData = uncert.mutualify(1e-8, 1e-8, False, theData, 1e-8, 1e-8, False)
                if relative:  # convert absolute to relative
                    uncert = uncert / theData
                else:  # relative to absolute
                    uncert = uncert * theData
            except Exception as err:
                print(len(uncert), uncert.copyDataToXYs()[0], uncert.copyDataToXYs()[-1])
                print(len(theData), theData.copyDataToXYs()[0], theData.copyDataToXYs()[-1])
                raise err
        return uncert

    def plot(self, title=None, scalelabel=None, xlim=None, ylim=None, xlog=False, ylog=False):
        """

        :param title:
        :param scalelabel:
        :param xlim:
        :param ylim:
        :param xlog:
        :param ylog:
        :return:
        """
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import QuadMesh

        # determine where to get the settings for the potentially mirrored second axis
        if isinstance(self.matrix.axes[1].values, linkModule.Link):
            axis1index = 2
        else:
            axis1index = 1
        x = self.matrix.axes[2].values
        y = self.matrix.axes[axis1index].values

        X, Y = np.meshgrid(x, y)
        XY = np.hstack((X.ravel()[:, np.newaxis], Y.ravel()[:, np.newaxis]))
        Z = (self.matrix.array.constructArray()).ravel()

        ax = plt.subplot(1, 1, 1)
        if title is None: title = str(self.toXLink())
        plt.suptitle(title)

        qc = QuadMesh(
            meshWidth=len(x) - 1,
            meshHeight=len(y) - 1,
            coordinates=XY,
            #            showedges=True,
            antialiased=True,
            shading='flat',
            transOffset=ax.transData)

        qc.set_array(Z)
        ax.add_collection(qc, autolim=True)

        if xlim is None:
            ax.set_xlim(x[0], x[-1])
        else:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim is None:
            ax.set_ylim(y[0], y[-1])
        else:
            ax.set_ylim(ylim[0], ylim[1])
        if xlog:
            ax.set_xscale('log')
        if ylog:
            ax.set_yscale('log')

        xlabel = self.matrix.axes[2].label + ' (' + self.matrix.axes[2].unit + ')'
        ylabel = self.matrix.axes[axis1index].label + ' (' + self.matrix.axes[axis1index].unit + ')'

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar = plt.colorbar(qc)
        if scalelabel is not None:
            cbar.set_label(scalelabel)
        else:
            cbar.set_label(str(self.type) + ' covariance (' + str(self.matrix.axes[0].unit) + ')')
        plt.show()

    def toXML_strList(self, indent='', **kwargs):
        """

        :param indent:
        :param kwargs:
        :return:
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s' % (indent, self.moniker)]
        if self.label is not None: xmlString[0] += ' label="%s"' % self.label
        xmlString[0] += ' type="%s"' % self.type
        if self.productFrame is not None:
            xmlString[0] += ' productFrame="%s"' % self.productFrame
        xmlString[0] += '>'
        xmlString += self.matrix.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Translate <covarianceMatrix> element from xml into python class."""

        xPath.append(element.tag)

        matrix_ = griddedModule.Gridded2d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        CM = cls(label=element.get('label'), type=element.get('type'), matrix=matrix_,
                 productFrame=element.get('productFrame'))

        xPath.pop()

        return CM
