# /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


import argparse
import numpy
import copy
import math


# ------------ safer matrix and inverse construction tools, "by hand" -----------
# ------------------------- these are the core routines -------------------------
def eigenvectors_from_orthoganal_matrix(O):
    v = []
    for i in range(O.shape[0]):
        v.append(O.T[i])
    return v


def matrix_from_eigendecomposition(e, v, ndim, doInverse=False, onlyLargeEVs=True, onlyPositiveEV=True,
                                   smallEVAbsTol=1e-10, smallEVRelTol=1e-6):
    B = numpy.zeros((ndim, ndim))
    evTol = max(smallEVAbsTol, smallEVRelTol * abs(e.max()))
    for i in range(ndim):
        if abs(e[i]) < evTol and onlyLargeEVs:
            continue
        if e[i] < 0.0 and onlyPositiveEV:
            continue
        if doInverse:
            x = 1.0 / e[i]
        else:
            x = e[i]
        B = B + numpy.outer(x * v[i], v[i].T)  # += operator causes casting exception for obscure numpy reasons
    return B


def pruned_matrix_inverse(A, onlyLargeEVs=True, onlyPositiveEV=True, smallEVAbsTol=1e-10, smallEVRelTol=1e-6):
    """Build inverse of :math:`A` \"by hand\", the safe way.   :math:`A` must admit an eigenvalue decomposition."""
    e, O = numpy.linalg.eig(A)
    v = eigenvectors_from_orthoganal_matrix(O)
    return matrix_from_eigendecomposition(1.0 / e, v, A.shape[0], doInverse=True, onlyLargeEVs=onlyLargeEVs,
                                          onlyPositiveEV=onlyPositiveEV, smallEVAbsTol=smallEVAbsTol,
                                          smallEVRelTol=smallEVRelTol)


def pruned_matrix(A, onlyLargeEVs=True, onlyPositiveEV=True, smallEVAbsTol=1e-10, smallEVRelTol=1e-6):
    """Rebuild :math:`A` \"by hand\", the safe way.   :math:`A` must admit an eigenvalue decomposition."""
    e, O = numpy.linalg.eig(A)
    v = eigenvectors_from_orthoganal_matrix(O)
    return matrix_from_eigendecomposition(e, v, A.shape[0], onlyLargeEVs=onlyLargeEVs, onlyPositiveEV=onlyPositiveEV,
                                          smallEVAbsTol=smallEVAbsTol, smallEVRelTol=smallEVRelTol)


def scale_off_diagonals(A, onlyScaleThese=None, scaleFactor=0.999999):
    """Sam's trick for getting UNCOR to cooperate: shrink off diagonal elements by some (small) factor"""
    B = copy.copy(A)
    ndim = B.shape[0]
    for i in range(ndim):
        for j in range(i):
            if i == j:
                continue
            if onlyScaleThese is None or (i, j) in onlyScaleThese or (j, i) in onlyScaleThese:
                B[i, j] = scaleFactor * B[i, j]
                B[j, i] = B[i, j]
    return B


def reduce_off_diagonals(corr_mat, thresh):
    """
    Sometimes easiest way to eliminate negative eigenvalues is to
    reduce off-diagonal portion of the correlation matrix
    """
    corr_mat_new = corr_mat[:, :]
    tmp = corr_mat[off_diagonals(corr_mat)]
    tmp[tmp > thresh] = thresh
    corr_mat_new[off_diagonals(corr_mat_new)] = tmp
    return corr_mat_new


def off_diagonals(matrix):
    """
    Return indices for all off-diagonal elements.

        >>> numpy.mat[ numpy.offdiag( mat ) ] *= -1
    """
    ilen, jlen = matrix.shape
    idx = [i for i in range(ilen) for j in range(jlen) if i != j]
    jdx = [j for i in range(ilen) for j in range(jlen) if i != j]
    return idx, jdx


def dot_product(*args):
    """
    Matrix multiplication (dot product) of all arguments:

        >>> numpy.dotproduct(V.T, A, V)   # returns V.T * A * V
    """
    from functools import reduce
    return reduce(numpy.dot, args)


# ------------------------- rebinning for ENDF matrices -------------------------
def rebin_matrix(arr, N=2):
    """Rebin array by factor of N (not necessarily multiple of 2)"""
    d, r = divmod(len(arr), N)
    if r:
        d += 1
    return [sum(arr[N * i:N * i + N]) / float(N) for i in range(d)]


def hist_interp_of_matrix(supergrid, mat, erows, ecols=None):
    """
    Put 'histogram' interpolation of matrix onto a new `supergrid`:

    :param supergrid: energy list to use for new matrix
    :param mat: original matrix to be interpolated
    :param erows: energy bins for rows/x-axis of original matrix
    :param ecols: energy bins for columns/y-axis, same as ex by default

    all points in `erows` and `ecols` must also be in supergrid,
    or ValueError is raised
    """
    if not ecols:
        ecols = erows[:]

    nvals = len(supergrid) - 1
    ret_mat = numpy.zeros((nvals, nvals))
    xidx = [supergrid.index(en) for en in erows]
    yidx = [supergrid.index(en) for en in ecols]

    for i in range(1, len(xidx)):
        for j in range(1, len(yidx)):
            ret_mat[xidx[i - 1]:xidx[i], yidx[j - 1]:yidx[j]] = mat[i - 1, j - 1]
    return ret_mat


# ------------------------- symmetric matrix tools ------------------------------
def switchSymmetry(mlist, upperToLower=True):
    """
    A symmetric 2-d NxN array can be stored in memory as a list of (N*(N+1)/2) numbers. The order depends
    on whether the upper-diagonal or lower-diagonal portion of the array is being stored.  This method switches
    between the two representations.
    Note that switching from upper to lower is the same as switching the memory order from 'column-major' to 'row-major'

    :param mlist: list, tuple or 1d array representing input matrix
    :param upperToLower: boolean, True = convert upper- to lower-symmetric, False = convert lower- to upper-symmetric
    :return: list with output matrix
    """
    shape = int(math.sqrt(2 * len(mlist)))
    arrays = [[] for i in range(shape)]
    matiter = iter(mlist)
    for idx in range(shape):
        if upperToLower:
            lbound, ubound = idx, shape
        else:
            lbound, ubound = 0, idx + 1
        for jdx in range(lbound, ubound):
            arrays[jdx].append(next(matiter))

    return [val for sublist in arrays for val in sublist]


# ------------------------- covariance matrix checks -------------------------
def check_real_and_finite(A):
    """Checks that all elements in a matrix are read and finite"""
    return numpy.all(numpy.isreal(A)) and numpy.all(numpy.isfinite(A))


def check_positive_semidefinite(A, warnAll=False, verbose=False):
    """Checks that all elements in a matrix are :math:`<= 0`"""
    success = True
    e, O = numpy.linalg.eig(A)
    for i, x in enumerate(e):
        if x < 0.0:
            if warnAll:
                if verbose:
                    print("FAIL  {i}th eigenvalue < 0.0: {x}".format(i=str(i), x=str(x)))
                success = False
            elif not numpy.allclose(x, 0.0):
                if verbose:
                    print("FAIL  {i}th eigenvalue < 0.0: {x}".format(i=str(i), x=str(x)))
                success = False
    return success


def check_symmetric(A, warnAll=False, verbose=False):
    """Checks whether a matrix is symmetric, i.e. :math:`A_{ij} = A_{ji}`"""
    success = True
    ndim = A.shape[0]
    for i in range(ndim):
        for j in range(i):
            if A[i, j] != A[j, i]:
                if warnAll:
                    if verbose:
                        print("FAIL  A[{i},{j}] != A[{j},{i}]: {val1} vs. {val2}".format(i=str(i), j=str(j),
                                                                                         val1=str(A[i, j]),
                                                                                         val2=str(A[j, i])))
                    success = False
                elif not numpy.allclose(A[i, j], A[j, i]) and not numpy.allclose(A[i, j], 0.0):
                    if verbose:
                        print("FAIL  A[{i},{j}] != A[{j},{i}]: {val1} vs. {val2}".format(i=str(i), j=str(j),
                                                                                         val1=str(A[i, j]),
                                                                                         val2=str(A[j, i])))
                    success = False
    return success


def check_covariance_element_bounds(A, warnAll=False, verbose=False):
    success = True
    ndim = A.shape[0]
    for i in range(ndim):
        for j in range(i):
            if A[j, j] >= 0.0 and A[i, i] >= 0.0 and abs(A[i, j]) > math.sqrt(A[i, i] * A[j, j]):
                if warnAll:
                    if verbose:
                        print("FAIL  abs( A[{i},{j}] ) > sqrt( A[{i},{i}] * A[{j},{j}] ): {val1} vs. {val2}".format(
                            i=str(i), j=str(j), val1=str(abs(A[i, j])), val2=str(math.sqrt(A[i, i] * A[j, j]))))
                    success = False
                elif not numpy.allclose(abs(A[i, j]), math.sqrt(A[i, i] * A[j, j])):
                    if verbose:
                        print("FAIL  abs( A[{i},{j}] ) > sqrt( A[{i},{i}] * A[{j},{j}] ): {val1} vs. {val2}".format(
                            i=str(i), j=str(j), val1=str(abs(A[i, j])), val2=str(math.sqrt(A[i, i] * A[j, j]))))
                    success = False
            if A[j, j] >= 0.0 and A[i, i] >= 0.0 and abs(A[i, j]) > 0.5 * (A[i, i] + A[j, j]):
                if warnAll:
                    if verbose:
                        print("FAIL  abs( A[{i},{j}] ) > 0.5 * ( A[{i},{i}] + A[{j},{j}] ): {val1} vs. {val2}".format(
                            i=str(i), j=str(j), val1=str(abs(A[i, j])), val2=str(0.5 * (A[i, i] + A[j, j]))))
                    success = False
                elif not numpy.allclose(abs(A[i, j]), math.sqrt(A[i, i] * A[j, j])):
                    if verbose:
                        print("FAIL  abs( A[{i},{j}] ) > 0.5 * ( A[{i},{i}] + A[{j},{j}] ): {val1} vs. {val2}".format(
                            i=str(i), j=str(j), val1=str(abs(A[i, j])), val2=str(0.5 * (A[i, i] + A[j, j]))))
                    success = False
            if A[j, j] >= 0.0 and A[i, i] >= 0.0 and not numpy.allclose(abs(A[i, j]), math.sqrt(A[i, i] * A[j, j])):
                if warnAll and verbose:
                    print("WARNING  abs(A[{i},{j}]) approx. sqrt( A[{i},{i}] * A[{j},{j}] ): {val1} vs. {val2}".format(
                        i=str(i), j=str(j), val1=str(abs(A[i, j])), val2=str(math.sqrt(A[i, i] * A[j, j]))))
    return success


# ------------------------- covariance matrix utilities -------------------------
def covariance_to_correlation(matrix):
    """Convert a covariance matrix to a correlation matrix"""
    diag = numpy.sqrt(matrix.diagonal())
    corr = matrix / diag / diag[:, numpy.newaxis]
    # now fix diagonal + remove any NaN (from div/0):
    corr[[range(len(corr)), range(len(corr))]] = 1.0  # must be exactly 1
    corr[numpy.isnan(corr)] = 0
    return corr


def correlation_to_covariance(correlationMatrix, uncertainty):
    """
    Convert a correlation matrix to a covariance matrix
    """
    covMtx = numpy.zeros_like(correlationMatrix)
    for _i, u in enumerate(uncertainty[0]):
        for _j, v in enumerate(uncertainty[0]):
            covMtx[_i, _j] = u * v * correlationMatrix[_i, _j]
    return covMtx


def covariance_to_relative(matrix, rsd):
    """Convert an absolute covariance matrix to a relative covariance matrix"""
    rsd = numpy.array(rsd)
    return matrix * rsd * rsd[:, numpy.newaxis]


def relative_to_covariance(relativeMatrix, variance):
    """
    Convert a relative covariance matrix to an absolute covariance matrix

    .. warning::
        not implemented
    """
    raise NotImplementedError()


def relative_to_correlation(relativeMatrix, data):
    """
    Convert a relative covariance matrix to a correlation matrix

    .. warning::
        not implemented
    """
    raise NotImplementedError()


def correlation_to_relative(correlationMatrix, data):
    """
    Convert a correlation matrix to a relative covariance matrix

    .. warning::
        not implemented
    """
    raise NotImplementedError()


def affine_transform_covariance(covarianceMatrix, transformMatrix):
    r"""
    Perform an affine transformation to a covariance matrix, namely
    :math:`new = A^T \cdot old \cdot A`, where :math:`A` is a orthogonal matrix.
    """
    return dot_product(transformMatrix.T, covarianceMatrix, transformMatrix)


def extract_uncertainty(covarianceMatrix):
    return numpy.sqrt(extract_variance(covarianceMatrix))


def extract_variance(covarianceMatrix):
    return numpy.diag(covarianceMatrix)


# ------------------------- matrix diff -------------------------
def diff_matrices(matrixOne, matrixTwo, printDiagnostics=True, quiet=True):
    if matrixOne.shape != matrixTwo.shape:
        raise ValueError('matrices have different shapes: ' + str(matrixOne.shape) + ' vs. ' + str(matrixTwo.shape))
    varianceOne = numpy.diag(matrixOne)
    varianceTwo = numpy.diag(matrixTwo)
    minVariance = min(varianceOne[varianceOne > 0.0].min(), varianceTwo[varianceOne > 0.0].min())

    # do the diffs
    _diff = matrixOne - matrixTwo
    _reldiff = numpy.nan_to_num(_diff / numpy.matrix(matrixOne))
    perdiff = 100.0 * _reldiff

    # stats on diff
    maxPerVal = perdiff.max()
    minPerVal = perdiff.min()
    maxAbsVal = _diff.max()
    minAbsVal = _diff.min()
    absMaxAbsVal = max(map(abs, (maxAbsVal, minAbsVal)))
    absMaxPerVal = max(map(abs, (maxPerVal, minPerVal)))
    flawCount = 0
    if printDiagnostics:
        print('Matrix diff diagnostics:')
        print('    * min non-zero variance in raw covariance matrices = ' + str(minVariance))
        if absMaxPerVal > 1.0:
            flawCount += 1
            print('    * max relative difference = ' + str(absMaxPerVal) + '% > 1.0%')
        else:
            print('    * max relative difference = ' + str(absMaxPerVal) + '%')
        if absMaxAbsVal > 0.01 * minVariance:
            flawCount += 1
            print('    * max absolute difference = ' + str(absMaxAbsVal) + ' > 1.0% * ' + str(minVariance))
        else:
            print('    * max absolute difference = ' + str(absMaxAbsVal))
        if flawCount == 2:
            print('    **** check this one ****')
            try:
                if not quiet:
                    subprocess.check_call(['say', 'holy', 'guacamole'])
            except subprocess.CalledProcessError:
                pass

    # return the numerical diffs
    return _diff, perdiff


# ------------------------- print/plot matrix -------------------------
def print_matrix(M, pretty=True, elementSize=8):
    """Simple matrix printer, makes little attempt to be pretty, but does print huge matrices"""
    for _i in range(M.shape[0]):
        for _j in range(M.shape[1]):
            if pretty:
                print('{:.2g}'.format(M[_i, _j]).ljust(elementSize), end=' ')
            else:
                print(M[_i, _j], '  ', end=' ')
        print()


def print_vector(v):
    for _i in range(len(v)):
        print(v[_i], '  ', )
    print()


def print_eigenvalues(_A):
    e, _O = numpy.linalg.eig(_A)
    myEVs = [x.real for x in e.tolist()]
    myEVs.sort()
    print('\t'.join([str(x) for x in myEVs]))


def plot_matrix(m, title="a matrix", scaling=None, scalingFloor=0.0):
    """
    :param m: a numpy.mat instance
    :param title: a string to use as the plot title
    :param scaling: either None, 'log', or 'asinh'
                    this scales the value of the matrix plotted in the following ways ::

                        * None : no scaling
                        * 'log' : each element is plotted as ln(x) -- good for covariances
                          which must always be positive semidefinite.  If scalingFloor > 0.0,
                          than we do ln(max(x, scalingFloor)).
                        * 'asinh' : each element is scaled as asinh(x).  This exaggerates
                          scale for values of abs(x)<1.0.
    :param scalingFloor: the minimum value of each element that gets plotted.
    """
    import matplotlib.pyplot as plt
    # Make plot with vertical (default) colorbar
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mm = None
    if scaling is None:
        mm = m
    elif scaling == 'asinh':
        mm = numpy.asinh(m)
    elif scaling == 'log':
        mm = numpy.copy(m)
        if scalingFloor > 0.0:
            mm[mm < scalingFloor] = scalingFloor
        mm = numpy.log(mm)
    maxAbsVal = mm.max()
    minAbsVal = mm.min()
    cax = ax.imshow(mm, interpolation='nearest')
    ax.set_title(title)
    ticks = None
    if scaling is None:
        ticks = [minAbsVal, 0, maxAbsVal / 2.0, maxAbsVal]
    elif scaling == 'asinh':
        ticks = [minAbsVal + _i * (maxAbsVal - minAbsVal) / 9 for _i in range(9)] + [maxAbsVal]
    elif scaling == 'log':
        ticks = [minAbsVal + _i * (maxAbsVal - minAbsVal) / 9 for _i in range(9)] + [maxAbsVal]
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar = fig.colorbar(cax, ticks=ticks)
    cbar.ax.set_yticklabels([str(t) for t in ticks])  # vertically oriented colorbar
    plt.show()


def plot_bad_eigenspaces(_A):
    e, _O = numpy.linalg.eig(_A)
    v = eigenvectors_from_orthoganal_matrix(_O)
    ndim = _A.shape[0]
    B = numpy.zeros((ndim, ndim))
    for _i in range(ndim):
        if e[_i] >= 0.0:
            continue
        B = B + numpy.outer(v[_i], v[_i].T)  # += operator causes casting exception for obscure numpy reasons
    print("Num of elements in original:", pow(_A.shape[0], 2))
    print("Num of elements in bad-space:", len(B[B > 0.0]))
    plot_matrix(B)


# ------------------------- matrix composition through stacking -------------------------
def stackVertical(l):
    """
    :param l: a list of numpy matrices, each with same `shape[1]` (i.e. same number of columns).
              Elements equal to None are ignored

    :returns: a matrix packed as follows ::

            [  l[0]  ]
            [ ------ ]
            [  l[1]  ]
            [ ------ ]
            [    :   ]
            [ ------ ]
            [ l[n-1] ]

        where n is the number of non-None elements in l
    """
    return numpy.vstack(list(filter(lambda x: x is not None, l)))


def stackHorizontal(l):
    """
    :param l: a list of numpy matrices, each with same shape[0] (i.e. same number of rows).
              Elements equal to None are ignored

    A note about numpy.array shapes and indexing ::

        a.shape = [ nRows, nCols ]

    So, we index through ``a[ iRow, iCol ]``

    :returns: a matrix packed as follows ::

            [  l[0]  | l[1]  | ... | l[n-1] ]

        where n is the number of non-None elements in l
    """
    return numpy.hstack(list(filter(lambda x: x is not None, l)))


def stackDiagonal(l):
    """
    :param l: a list of numpy matrices, elements equal to None are ignored

    A note about numpy.array shapes and indexing ::

        a.shape = [ nRows, nCols ]

    So, we index through ``a[ iRow, iCol ]``

    :returns: a matrix packed as follows ::

            [  l[0]  |  0.0  |  0.0 | ... |   0.0  ]
            [ ------ | ----- | ---- | ... | ------ ]
            [   0.0  | l[1]  |  0.0 | ... |   0.0  ]
            [ ------ | ----- | ---- | ... | ------ ]
            [   0.0  |  0.0  | l[2] | ... |   0.0  ]
            [ ------ | ----- | ---- | ... | ------ ]
            [    :   |   :   |  :   |     |    :   ]
            [ ------ | ----- | ---- | ... | ------ ]
            [   0.0  |  0.0  |  0.0 | ... | l[n-1] ]

        where n is the number of non-None elements in l

    """
    # Determine the size of the final results matrix
    newShape = [0, 0]
    for m in l:
        if m is None:
            continue
        newShape[0] += m.shape[0]  # Keep adding to the number of rows
        newShape[1] += m.shape[1]  # Keep adding to the number of columns
    result = numpy.zeros(newShape)

    # Add matrix to results
    iStart = 0
    jStart = 0
    for m in l:
        if m is None:
            continue
        for _i in range(m.shape[0]):  # Loop over rows
            for _j in range(m.shape[1]):  # Loop over columns
                result[iStart + _i, jStart + _j] = m[_i, _j]
        iStart += m.shape[0]  # Compute new start row for next matrix
        jStart += m.shape[1]  # Compute new start column for next matrix
    return result


# ------------------------- constrained generalized least-squares -------------------------
def cglsqrSolve(data, dataUnc=None, dataCov=None,
                kernel=None,
                prior=None, priorCov=None,
                constraintVector=None, constraintMatrix=None):
    """
    Constrainted Generalized Least-Squares Solver, based on CorAL routines.
    We are minimizing the following ``chi^2`` (assuming Gaussian statistics) ::

        chi^2 = ( data - kernel * model ) * ( dataCov )^-1 * ( data - kernel * model )^T +
                ( prior - model ) * ( priorCov )^2 * ( prior - model )^T

    subject to the constraint ::

        constraintVector = constraintMatrix * model

    Here, ``^T`` means matrix transpose and ``^-1`` means a (generalized) matrix inverse.


    :returns: a tuple: ``(model, modelCovariance)``:

        - ``model`` : a M x 1 numpy.mat containing the extracted model parameter

        - ``modelCovariance`` : a M x M numpy.mat containing the covariance on the extracted model paremeters


    **Mandatory arguments:**

        - **data** : a N x 1 ``numpy.mat`` (a vector!) containing the data to fit


    **Optional arguments:**

        - **kernel** : a N x M ``numpy.mat`` that maps the model parameters into the data space.  If this
                   is not given, it is assumed that ``kernel`` = the identity matrix and N == M.

        - either

            -- **dataUnc** : a N x 1 ``numpy.mat`` (a vector) of uncertianties on the data vector.
                         This will be converted to the ``dataCov`` if the ``dataCov`` is not specified.

            -- **dataCov** : a N x N ``numpy.mat`` containing the data's covariance
                         If this is specified, the ``dataUnc`` will be ignored.

            If neither the ``dataUnc`` or ``dataCov`` are specified, we will set the covariance to the
            identity matrix. Instead of minimizing the ``chi^2``, you will minimize ::

                | data - kernel * model |^2

        - both of

            -- **prior**    : a M x 1 ``numpy.mat`` containing the fitting model's apriori values

            -- **priorCov** : a M x M ``numpy.mat`` containing the fitting model's apriori covariance

            If either one is not included, the second term in the ``chi^2`` will be ignored.

        - both of

            -- **constraintVector** : a L x 1 ``numpy.mat`` containing values to constrain the model to match

            -- **constraintMatrix** : a L x M ``numpy.mat`` relating the model parameters to the ``constraintVector``
                                      values

            If either one is not included, the constraint equation will be ignored.


    .. rubric:: HOW IT WORKS

    Because the data and the prior are independent and not correlated, we will solve the minimization problem by
    stacking the prior part and the data part to construct a new ``chi^2`` to minimize ::

                  [ data  ]             [ dataCov |     0    ]                [ kernel ]
        newData = [ ----- ],   newCov = [ --------|--------- ],   newKernel = [ ------ ]
                  [ prior ]             [    0    | priorCov ]                [    1   ]

    With these, the ``chi^2`` may be rewritten as ::

        chi^2 = ( newData - newKernel * model )^T * newCov^-1 * ( newData - newKernel * model )

    The solution to the minimization problem is the well known normal equation ::

        modelCov = ( newKernel * newCov^-1 * newKernel^T )^-1

        model = modelCov * newKernel^T * newCov^-1 * newData

    Implementing this case where no uncertainty or covariance is given on the data vector is straightforward.
    In this case, we are solving the simple least-squares problem ::

        minimize:   | data - kernel * model |^2 + ( model - prior )*( priorCov )^-1*( model - prior )^T

    In other words, we use a ``dataCov = 1``, the identity matrix.

    We comment that the equality constraints can be added in much this same way by adding the constraint
    vector as a fake data set, with infinitesimal uncertainties ::

        chi^2 += lambda * ( constraintVector - constraintMatrix * model )^2

    with lambda :math:`\\rightarrow\\infty`.  This approach works, but can lead to numerical instabilities by
    adding big number to the original covariances, then taking inverses.


    So, rather than doing this, we will do something different (using Lagrange multipliers).
    We use the "newData" etc. above and the corresponding ``chi^2`` ::

        chi^2 = ( newData - newKernel * model )^T * newCov^-1 * ( newData - newKernel * model )

    subject to the constraint ::

        constraintMatrix * model = constraintVector

    Here we are really solving the quadratic optimization problem here and the way to solve it is with
    Lagrange multipliers.  So, extend the model thusly ::

                   [   model  ]
        newModel = [ ---------]
                   [  lambda  ]

    Here lambda is the vector of Lagrange multipliers ( there are L of them ).
    Minimizing the chi^2, we find the usual normal equation ::

        ( newKernel * newCov^-1 * newKernel^T ) * model = newKernel^T * newCov^-1 * newData

    So, we'll stack the normal and constraint equations::

            [ newKernel * newCov^-1 * newKernel^T | constraintMatrix^T ]
        A = [ ------------------------------------|------------------- ]
            [          constraintMatrix           |         0          ]

        and

            [ newKernel^T * newCov^-1 * newData ]
        y = [ --------------------------------- ]
            [         constraintVector          ]

        giving,

        A * newModel = y

    We can solve this a few ways:

        1.  Performing a QR decomposition on :math:`A` gives :math:`A = QR` (using ``numpy.linalg.qr``) and ::

                R * newModel = Q^T * y

            We can then use the ``numpy.linalg.tensorsolve`` for the ``newModel`` ::

                newModel = numpy.linalg.tensorsolve( R, Q^T * y )

        2.  Just doing a Moore-Penrose inversion (uses SVD) ::

                newModel = A^-1 * y

    But now the meaning of A is clear: it is an "extended" covariance where the first
    M x M block is the covariance of the model and the rest are the covariance of the
    Lagrange multipliers (which are disposable).


    .. rubric:: TESTING

    Note on construction of dataCov with uncertainty:
        >>> unc = numpy.array( [ 1.0, 1.0, 2.0 ] )
        >>> diagCov = numpy.identity(unc.shape[0]) * unc**2
        >>> diagCov
        array([ [ 1.,  0.,  0.],
                [ 0.,  1.,  0.],
                [ 0.,  0.,  4.]])

    Also, without specifying an uncertainty or covariance, one finds:
        >>> dataCov = numpy.identity( data.shape[1] )
        >>> dataCov
        array([ [ 1.,  0.,  0.],
                [ 0.,  1.,  0.],
                [ 0.,  0.,  1.]])

    So, both alternate methods of defining the data covariance function correctly.


    **Define the problem**

    Define a test problem:
        >>> answer = numpy.array([[ 1.34883721,-0.69767442, 0.34883721, 0.1, 42.0]])
        >>> kernel = numpy.array( [ [ 1.0, 2.0, 3.0, 0.0, 0.0 ], [ 2.0, 3.0, 4.0, 0.0, 0.0 ],
        >>>                       [ 4.5, 5.4, 2.0, 0.0, 0.0 ] ] )  # Note: lower two subspaces map to 0
        >>> numpy.dot( kernel, answer.T )
            array([[1.],
                   [2.],
                   [3.]])

    Test data:
        >>> data = numpy.array([ [ 1.1, 1.89, 3.05 ] ] )
        >>> dataCov = numpy.array([ [ 1.0, 0.1, 0.1 ], [ 0.1, 1.0, 0.1 ] , [ 0.1, 0.1, 1.0 ] ] )

    Constrain the last two elements of the model to add to 42.1:
        >>> constraintVector = numpy.array( [[ 42.1 ]] )
        >>> constraintMatrix = numpy.array( [[ 0.0, 0.0, 0.0, 1.0, 1.0 ]] )

    A apriori guess to the result (good to < 5% in all 4 dimensions):
        >>> prior = numpy.array([[ 1.3, -0.7, 0.3, 0.11, 42.2 ]])
        >>> priorCov = numpy.array([
            [  0.07      ,   0.        ,   0.        ,   0.        ,   0.        ],
            [  0.        ,   0.5       ,   0.        ,   0.        ,   0.        ],
            [  0.        ,   0.        ,   0.4       ,   0.        ,   0.        ],
            [  0.        ,   0.        ,   0.        ,   0.1       ,   0.        ],
            [  0.        ,   0.        ,   0.        ,   0.        ,  20.        ]])

    This might be a little too good...


    **A data only solution**

    Using data only, we get:
        >>> modelCov = numpy.linalg.pinv( dot_product( kernel.T, numpy.linalg.pinv( dataCov ), kernel ) )
        >>> modelCov
        array([[ 19.2436993 , -17.66414278,   4.23904813,   0.        ,   0.        ],
               [-17.66414278,  16.28177393,  -3.95484045,   0.        ,   0.        ],
               [  4.23904813,  -3.95484045,   1.03439697,   0.        ,   0.        ],
               [  0.        ,   0.        ,   0.        ,   0.        ,   0.        ],
               [  0.        ,   0.        ,   0.        ,   0.        ,   0.        ]])
        >>> model = dot_product( modelCov, kernel.T, numpy.linalg.pinv( dataCov ), data.T )
        >>> model
        array([[ 0.66232558],
               [-0.05465116],
               [ 0.18232558],
               [ 0.        ],
               [ 0.        ]])

    Note: the last to elements in the model are zero because there was no way to control them given the kernel in play.
    Without a handle on the last two elements, the fitting cannot do better than this.


    **A data+prior solution**

    With data+prior, we first repack things:
        >>> newData = stackHorizontal( [ data, prior ] )
        >>> newData
        array([[  1.1 ,   1.89,   3.05,   1.3 ,  -0.7 ,   0.3 ,   0.11,  42.2 ]])
        >>> newKernel = stackVertical( [ kernel, numpy.identity( prior.shape[1] ) ] )
        >>> newKernel
        array([[ 1. ,  2. ,  3. ,  0. ,  0. ],
               [ 2. ,  3. ,  4. ,  0. ,  0. ],
               [ 4.5,  5.4,  2. ,  0. ,  0. ],
               [ 1. ,  0. ,  0. ,  0. ,  0. ],
               [ 0. ,  1. ,  0. ,  0. ,  0. ],
               [ 0. ,  0. ,  1. ,  0. ,  0. ],
               [ 0. ,  0. ,  0. ,  1. ,  0. ],
               [ 0. ,  0. ,  0. ,  0. ,  1. ]])

    Now, we can rework the algorithm above:
        >>> newCov = stackVertical([stackHorizontal([dataCov, numpy.zeros((dataCov.shape[0], priorCov.shape[1]))]),
                                    stackHorizontal([numpy.zeros((priorCov.shape[0], dataCov.shape[1])), priorCov])])
        >>> newCov
        array([[  1.  ,   0.1 ,   0.1 ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ],
               [  0.1 ,   1.  ,   0.1 ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ],
               [  0.1 ,   0.1 ,   1.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ],
               [  0.  ,   0.  ,   0.  ,   0.07,   0.  ,   0.  ,   0.  ,   0.  ],
               [  0.  ,   0.  ,   0.  ,   0.  ,   0.5 ,   0.  ,   0.  ,   0.  ],
               [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.4 ,   0.  ,   0.  ],
               [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.1 ,   0.  ],
               [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ,  20.  ]])
        # use Moore-Penrose generalized inverse (ie SVD inversion)
        >>> modelCov = numpy.linalg.pinv( dot_product( newKernel.T, numpy.linalg.pinv( newCov ), newKernel ) )
        >>> modelCov
        array([[  6.30909610e-02,  -5.01795451e-02,   5.99338297e-03,  0.00000000e+00,   0.00000000e+00],
               [ -5.01795451e-02,   9.30182090e-02,  -5.02878134e-02,  0.00000000e+00,   0.00000000e+00],
               [  5.99338297e-03,  -5.02878134e-02,   7.63220083e-02,  0.00000000e+00,   0.00000000e+00],
               [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  1.00000000e-01,   0.00000000e+00],
               [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  0.00000000e+00,   2.00000000e+01]])
        >>> model = dot_product( modelCov, newKernel.T, numpy.linalg.pinv( newCov ), newData.T )
        >>> model
        array([[  1.30359097],
               [ -0.64662084],
               [  0.32428234],
               [  0.11      ],
               [ 42.2       ]])

    A good prior gives a good result:
        >>> residual = data.T - dot_product( kernel, model )
        >>> residual
        array([[ 0.11680368],
               [-0.0744488 ],
               [ 0.02702848]])
        >>> error = answer.T - model
        >>> error
        array([[ 0.04524624],
               [-0.05105358],
               [ 0.02455487],
               [-0.01      ],
               [-0.2       ]])

    The final uncertainty is not bad:
        >>> uncertainty = [ math.sqrt(modelCov[i,i]) for i in range( modelCov.shape[0] ) ]
        >>> uncertainty
            [0.25117914115531315, 0.30498886706247513, 0.2762643810785747, 0.31622776601683794, 4.47213595499958]


    **A data+constraint+prior solution**

    Next, with everything (data+constraint+prior)...  First we stack the normal equation and the constraint equation:
        >>> reallyNewData = stackHorizontal([(newKernel.T * numpy.linalg.pinv(newCov) * newData.T).T, constraintVector])
        >>> reallyNewData
        matrix([[ 35.04920635,  19.82814815,  14.56111111, 1.1, 2.11,  42.1 ]])
        >>> reallyNewKernel = stackVertical( [ newKernel.T*numpy.linalg.pinv( newCov )*newKernel, constraintMatrix ] )
        >>> reallyNewKernel
        matrix([[ 37.13293651,  28.66666667,  15.97222222,   0.        ,   0.        ],
                [ 28.66666667,  38.82962963,  23.33333333,   0.        ,   0.        ],
                [ 15.97222222,  23.33333333,  27.22222222,   0.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,  10.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,   0.        ,   0.05      ],
                [  0.        ,   0.        ,   0.        ,   1.        ,   1.        ]])

    Now we just invert the reallyNewKernel
        >>> model = numpy.linalg.pinv( reallyNewKernel ) * reallyNewData.T
        matrix([[  1.30359097],
                [ -0.64662084],
                [  0.32428234],
                [  0.10999476],
                [ 41.99052891]])

    The fit is just a little better than with the prior alone -- but now the last element is nailed down:
        >>> residual = data.T - kernel*model
        matrix([[ 0.11680368],
                [-0.0744488 ],
                [ 0.02702848]])
        >>> error = answer.T - model
        >>> error
        matrix([[ 0.04524624],
                [-0.05105358],
                [ 0.02455487],
                [-0.00999476],
                [ 0.00947109]])
        >>> constraintMatrix*model
        matrix([[ 42.10052368]])

    We compute the covariance as the inverse of the reallyNewKernel:
        >>> modelCov = numpy.linalg.pinv( reallyNewKernel )
        >>> modelCov
        matrix([[  6.30909610e-02,  -5.01795451e-02,   5.99338297e-03,  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                [ -5.01795451e-02,   9.30182090e-02,  -5.02878134e-02,  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                [  5.99338297e-03,  -5.02878134e-02,   7.63220083e-02,  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,  9.99975063e-02,  -4.98740680e-04,   2.49370340e-05],
                [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00, -9.97481360e-02,   5.03728087e-02,   9.97481360e-01]])

    The new uncertainties are much reduced:
        >>> uncertainty = [ math.sqrt(modelCov[i,i]) for i in range( modelCov.shape[0] ) ]
        >>> uncertainty
        [0.25117914115531315, 0.30498886706247513, 0.2762643810785747, 0.316223823100982, 0.22443887510442173]


    **A data+constraint solution**

    Let's see how well we do using just the constraints and the data fitting.  As before, we stack the
    normal and constraint equation, but this time use the old kernel and old data (that is, without the prior):
        >>> reallyNewData = stackHorizontal([(kernel.T * numpy.linalg.pinv( dataCov ) * data.T ).T, constraintVector])
        >>> reallyNewData
        matrix([[ 16.47777778,  21.22814815,  13.81111111,   0.        , 0.        ,  42.1       ]])
        >>> reallyNewKernel = stackVertical( [ kernel.T*numpy.linalg.pinv( dataCov )*kernel, constraintMatrix ] )
        >>> reallyNewKernel
        matrix([[ 22.84722222,  28.66666667,  15.97222222,   0.        ,   0.        ],
                [ 28.66666667,  36.82962963,  23.33333333,   0.        ,   0.        ],
                [ 15.97222222,  23.33333333,  24.72222222,   0.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,   0.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,   0.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,   1.        ,   1.        ]])

    Now we just invert the reallyNewKernel
        >>> model = numpy.linalg.pinv( reallyNewKernel ) * reallyNewData.T
        >>> model
        matrix([[  0.66232558],
                [ -0.05465116],
                [  0.18232558],
                [ 21.05      ],
                [ 21.05      ]])

    Ack!  it made sure the constraint is obeyed, by divying up the 42.1 among the two uncontrolled model parameters!
        >>> residual = data.T - kernel*model
        >>> residual
        matrix([[ -1.07913678e-13],
                [ -2.13828955e-13],
                [ -4.26325641e-13]])
        >>> error = answer.T-model
        >>> error
        matrix([[  0.68651163],
                [ -0.64302326],
                [  0.16651163],
                [-20.95      ],
                [ 20.95      ]])

    What is telling is that, to accomodate the constraint, all the uncertainty had to be shifted to the first
    three components of the model:
        >>> modelCov = numpy.linalg.pinv( reallyNewKernel ) # use Moore-Penrose generalized inverse (ie SVD inversion)
        >>> modelCov
        matrix([[ 19.2436993 , -17.66414278,   4.23904813,   0.        ,0.        ,   0.        ],
                [-17.66414278,  16.28177393,  -3.95484045,   0.        ,0.        ,   0.        ],
                [  4.23904813,  -3.95484045,   1.03439697,   0.        ,0.        ,   0.        ],
                [  0.        ,   0.        ,   0.        ,   0.        ,0.        ,   0.5       ],
                [  0.        ,   0.        ,   0.        ,   0.        ,0.        ,   0.5       ]])
        >>> uncertainty = [ math.sqrt(modelCov[i,i]) for i in range( modelCov.shape[0] ) ]
        >>> uncertainty
        [4.386764103176406, 4.035068020722195, 1.0170530818673396, 0.0, 0.0]

    The moral is that we'd better know what spaces are constrainted by data and which ones are not!!!
    """
    # Check types of all arguments
    for x in [data, dataUnc, dataCov, kernel, prior, priorCov, constraintVector, constraintMatrix]:
        if x is not None and not isinstance(x, numpy.ndarray):
            raise TypeError("all arguments must be None or a numpy.mat, got " + str(type(x)))

    # Make sure the vector arguments are actually vectors
    if 1 not in data.shape:
        raise TypeError("data argument has wrong shape " + str(data.shape) + ', it should be a vector')
    if constraintVector is not None and 1 not in constraintVector.shape:
        raise TypeError("constraintVector argument has wrong shape " + str(constraintVector.shape) +
                        ', it should be a vector')
    if prior is not None and 1 not in prior.shape:
        raise TypeError("prior argument has wrong shape " + str(prior.shape) + ', it should be a vector')

    # Make sure all vectors are column vectors
    if data.shape[0] != 1:
        print("WARNING: data vector not a column vector, transposing it")
        data = data.T
    if constraintVector is not None and constraintVector.shape[0] != 1:
        print("WARNING: constraintVector vector not a column vector, transposing it")
        constraintVector = constraintVector.T
    if prior is not None and prior.shape[0] != 1:
        print("WARNING: prior vector not a column vector, transposing it")
        prior = prior.T

    # Construct a data covariance if we don't have one or if we have only uncertainties
    if dataCov is None and dataUnc is not None:
        dataCov = numpy.diag([dataUnc[0, i] * dataUnc[0, i] for i in range(dataUnc.shape[1])])
    if dataCov is None:
        dataCov = numpy.identity(data.shape[1])

    # Stack things to accomodate any possible prior for the model
    if prior is not None and priorCov is not None:
        newData = stackHorizontal([data, prior])
        newKernel = stackVertical([kernel, numpy.identity(prior.shape[1])])
        newCov = stackVertical([stackHorizontal([dataCov, numpy.zeros((dataCov.shape[0], priorCov.shape[1]))]),
                                stackHorizontal([numpy.zeros((priorCov.shape[0], dataCov.shape[1])),
                                                 priorCov])])  # incorrect -- 0 should be zeros()
    else:
        newData = data
        newKernel = kernel
        newCov = dataCov

    # No constraints, so perform the inversion using standard approach
    if constraintVector is None or constraintMatrix is None:
        modelCov = numpy.linalg.pinv(dot_product(newKernel.T, numpy.linalg.pinv(newCov),
                                                 newKernel))  # use Moore-Penrose generalized inverse (ie SVD inversion)
        model = dot_product(modelCov, newKernel.T, numpy.linalg.pinv(newCov), newData.T)

    # OK, there are constraints so perform the inversion using the Lagrange multiplier approach
    else:
        nConstraints = max(constraintVector.shape)
        invNewCov = numpy.linalg.pinv(newCov)
        reallyNewData = stackHorizontal([dot_product(newKernel.T, invNewCov, newData.T).T, constraintVector])
        reallyNewKernel = stackHorizontal(
            [stackVertical([dot_product(newKernel.T, invNewCov, newKernel), constraintMatrix]),
             stackVertical([constraintMatrix.T, numpy.zeros((nConstraints, nConstraints))])])
        modelCov = numpy.linalg.pinv(reallyNewKernel)  # use Moore-Penrose generalized inverse (ie SVD inversion)
        model = dot_product(modelCov, reallyNewData.T)
        # Note: extra dimensions are the fits of the Lagrange multipliers.  We don't need or want them.
        for i in range(nConstraints):
            model = numpy.delete(model, -1, 0)
            modelCov = numpy.delete(modelCov, -1, 0)
            modelCov = numpy.delete(modelCov, -1, 1)

    # How good was our inversion?
    fs = fit_statistics(data, dataCov, kernel, model)

    return model.T, modelCov, fs['residual'], fs['chi2']


def fit_statistics(data, dataCov, kernel, model):
    modelData = (dot_product(kernel, model)).T
    if data.shape != modelData.shape:
        raise ValueError('dimension mis-match:', str(data.shape), '!=', str(modelData.shape))
    residual = data - modelData
    chi2 = dot_product(residual, numpy.linalg.pinv(dataCov), residual.T)
    ndf = max(data.shape) - max(model.shape)  # an estimate, crude?
    return {'chi2': chi2, 'residual': residual, 'chi2/ndf': chi2 / ndf, 'ndf': ndf}


# ------------------------- create test matrix -------------------------
def get_test_matrix(endfFile=None, MT=None, MF=None):
    if endfFile is None:
        endfFile = 'n-099_Es_254m1.endf'
        ENDF_MFMT = "33,52"
        dimensions = "36,36"
        ndim = 36
        the_data = map(float, '''
         0.000000e+00
         0.000000e+00  0.000000e+00
         0.000000e+00  0.000000e+00  0.000000e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.690124e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.701406e+00  1.716205e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.727325e+00  1.740627e+00  1.768993e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.766112e+00  1.779705e+00  1.806902e+00  1.849737e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.805003e+00  1.818895e+00  1.846707e+00  1.888721e+00  1.932443e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.633234e+00  2.656310e+00  2.702535e+00  2.769406e+00  2.839445e+00  4.537071e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.633419e+00  2.656453e+00  2.702643e+00  2.772146e+00  2.842112e+00  4.532265e+00  4.536993e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.621032e+00  2.643952e+00  2.689979e+00  2.759142e+00  2.828787e+00  4.510953e+00  4.511160e+00  4.494515e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.294322e+00  2.314437e+00  2.354731e+00  2.415354e+00  2.476366e+00  3.932745e+00  3.936493e+00  3.918021e+00  3.422690e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.925894e+00  1.942689e+00  1.976455e+00  2.025372e+00  2.076478e+00  3.270196e+00  3.270376e+00  3.255084e+00  2.845935e+00  2.371981e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.560602e+00  1.574192e+00  1.601433e+00  1.641092e+00  1.680763e+00  2.588509e+00  2.591054e+00  2.581430e+00  2.259491e+00  1.886678e+00  1.513082e+00
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.112542e+00  1.121110e+00  1.139352e+00  1.166373e+00  1.193373e+00  1.747292e+00  1.750786e+00  1.744347e+00  1.533591e+00  1.289816e+00  1.048526e+00  7.535799e-01
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  9.208774e-01  9.279808e-01  9.422195e-01  9.642405e-01  9.854623e-01  1.388816e+00  1.392019e+00  1.387101e+00  1.223956e+00  1.035605e+00  8.521885e-01  6.280513e-01  5.345716e-01
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  6.868322e-01  6.913886e-01  7.019228e-01  7.175892e-01  7.333606e-01  9.928217e-01  9.953127e-01  9.930198e-01  8.790898e-01  7.481786e-01  6.233050e-01  4.717560e-01  4.089519e-01  3.193141e-01
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  5.577461e-01  5.614184e-01  5.699828e-01  5.822554e-01  5.944450e-01  7.851565e-01  7.871801e-01  7.855044e-01  6.967945e-01  5.954048e-01  4.998745e-01  3.846700e-01  3.370624e-01  2.655876e-01  2.228462e-01
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  3.322102e-01  3.343361e-01  3.386220e-01  3.454492e-01  3.517908e-01  4.290979e-01  4.304639e-01  4.297462e-01  3.834943e-01  3.312716e-01  2.847126e-01  2.295958e-01  2.071083e-01  1.677455e-01  1.427551e-01  9.555420e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  3.604038e-02  3.589246e-02  3.488834e-02  3.411414e-02  3.282049e-02 -3.279601e-02 -3.284671e-02 -3.271555e-02 -2.635367e-02 -1.497233e-02  1.232291e-03  2.323535e-02  3.246637e-02  3.503985e-02  3.383028e-02  3.062795e-02  2.692990e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -7.267533e-02 -7.373635e-02 -7.662602e-02 -8.038130e-02 -8.453273e-02 -1.979285e-01 -1.987893e-01 -1.986476e-01 -1.743977e-01 -1.409929e-01 -1.027300e-01 -5.400185e-02 -3.396930e-02 -1.587535e-02 -8.317105e-03  5.515628e-03  2.637916e-02  3.688719e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -8.439620e-02 -8.558441e-02 -8.881654e-02 -9.278483e-02 -9.716333e-02 -2.088966e-01 -2.102492e-01 -2.101251e-01 -1.860955e-01 -1.519483e-01 -1.129048e-01 -6.362840e-02 -4.343154e-02 -2.386577e-02 -1.535218e-02  1.041874e-03  2.649154e-02  3.981114e-02  4.414266e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -4.099937e-02 -4.161169e-02 -4.416221e-02 -4.686749e-02 -5.018817e-02 -1.297361e-01 -1.313929e-01 -1.316491e-01 -1.180110e-01 -9.560794e-02 -6.907117e-02 -3.455066e-02 -2.097080e-02 -7.999232e-03 -3.016686e-03  7.647215e-03  2.611127e-02  3.750344e-02  4.216583e-02  4.208217e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -8.692396e-03 -9.279179e-03 -1.108010e-02 -1.298799e-02 -1.528045e-02 -7.124834e-02 -7.270001e-02 -7.325542e-02 -6.765080e-02 -5.422922e-02 -3.658099e-02 -1.347648e-02 -4.733320e-03  3.182242e-03  5.588162e-03  1.224217e-02  2.591315e-02  3.602061e-02  4.106297e-02  4.208267e-02  4.331828e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.277385e-01  1.284242e-01  1.292478e-01  1.310071e-01  1.326888e-01  1.657041e-01  1.644256e-01  1.627654e-01  1.391348e-01  1.179997e-01  1.002280e-01  7.962980e-02  6.993322e-02  5.707190e-02  4.833213e-02  3.556557e-02  2.356314e-02  2.440651e-02  2.868411e-02  3.400294e-02  3.815346e-02  4.627409e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.988736e-01  2.001158e-01  2.025712e-01  2.065709e-01  2.102171e-01  2.895255e-01  2.889460e-01  2.871065e-01  2.485628e-01  2.095750e-01  1.732604e-01  1.299079e-01  1.105119e-01  8.630489e-02  7.158975e-02  4.787028e-02  2.091391e-02  1.543294e-02  1.850572e-02  2.567564e-02  3.119235e-02  4.617069e-02  5.066840e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.005723e-01  2.021068e-01  2.048516e-01  2.088673e-01  2.128834e-01  2.960696e-01  2.956199e-01  2.941737e-01  2.555868e-01  2.155465e-01  1.779475e-01  1.326661e-01  1.127036e-01  8.765646e-02  7.258235e-02  4.778240e-02  1.847589e-02  1.132459e-02  1.356651e-02  2.052769e-02  2.579265e-02  4.139233e-02  4.675601e-02  4.376726e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.682710e-01  1.693323e-01  1.714251e-01  1.747393e-01  1.780599e-01  2.444032e-01  2.436937e-01  2.425362e-01  2.109347e-01  1.783747e-01  1.478118e-01  1.113854e-01  9.527406e-02  7.458976e-02  6.200201e-02  4.128731e-02  1.688060e-02  1.109806e-02  1.307890e-02  1.888816e-02  2.330460e-02  3.623141e-02  4.049497e-02  3.776391e-02  3.280182e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.379456e-01  1.387843e-01  1.404562e-01  1.430018e-01  1.454850e-01  1.932246e-01  1.928924e-01  1.916331e-01  1.665959e-01  1.412460e-01  1.182012e-01  9.084645e-02  7.859059e-02  6.231665e-02  5.215703e-02  3.560402e-02  1.647353e-02  1.247273e-02  1.448556e-02  1.930557e-02  2.300591e-02  3.324171e-02  3.619012e-02  3.352310e-02  2.917429e-02  2.625984e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.229217e-01  1.236599e-01  1.251312e-01  1.273135e-01  1.292684e-01  1.674026e-01  1.670501e-01  1.659330e-01  1.439693e-01  1.224219e-01  1.033336e-01  8.060127e-02  7.037075e-02  5.639522e-02  4.739974e-02  3.290179e-02  1.652968e-02  1.351233e-02  1.555447e-02  1.992412e-02  2.327636e-02  3.214322e-02  3.442959e-02  3.170202e-02  2.767207e-02  2.503058e-02  2.405622e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.215488e-01  1.222745e-01  1.235435e-01  1.256653e-01  1.277727e-01  1.661270e-01  1.655886e-01  1.645113e-01  1.428780e-01  1.214350e-01  1.023434e-01  7.961931e-02  6.937311e-02  5.541419e-02  4.651601e-02  3.221869e-02  1.600271e-02  1.301406e-02  1.505617e-02  1.936928e-02  2.262628e-02  3.155131e-02  3.380264e-02  3.117874e-02  2.728977e-02  2.468726e-02  2.365762e-02  2.345737e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.209428e-01  1.216868e-01  1.229914e-01  1.251549e-01  1.273144e-01  1.678512e-01  1.675028e-01  1.664139e-01  1.447007e-01  1.228889e-01  1.030352e-01  7.963224e-02  6.891763e-02  5.464836e-02  4.576022e-02  3.136029e-02  1.486404e-02  1.154498e-02  1.337069e-02  1.760547e-02  2.089117e-02  2.985604e-02  3.244487e-02  3.002768e-02  2.626797e-02  2.368499e-02  2.268672e-02  2.244466e-02  2.164516e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.191637e-01  1.199028e-01  1.213634e-01  1.234824e-01  1.255984e-01  1.662982e-01  1.660038e-01  1.649346e-01  1.434091e-01  1.217682e-01  1.021530e-01  7.867281e-02  6.808521e-02  5.390608e-02  4.509330e-02  3.077834e-02  1.430499e-02  1.084618e-02  1.256069e-02  1.669642e-02  1.990212e-02  2.880841e-02  3.149815e-02  2.921774e-02  2.552252e-02  2.299790e-02  2.201084e-02  2.181507e-02  2.099733e-02  2.046324e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.176283e-01  1.183546e-01  1.197936e-01  1.219420e-01  1.240511e-01  1.648989e-01  1.646051e-01  1.635512e-01  1.422097e-01  1.207601e-01  1.012042e-01  7.779795e-02  6.726494e-02  5.316105e-02  4.444636e-02  3.021617e-02  1.374686e-02  1.015227e-02  1.174391e-02  1.577280e-02  1.890665e-02  2.777053e-02  3.053308e-02  2.837224e-02  2.476958e-02  2.231187e-02  2.132974e-02  2.114385e-02  2.037711e-02  1.986223e-02  1.933972e-02
         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.147770e-01  1.154951e-01  1.169225e-01  1.190371e-01  1.211243e-01  1.625308e-01  1.622074e-01  1.611758e-01  1.404440e-01  1.193101e-01  9.967764e-02  7.631129e-02  6.581783e-02  5.182580e-02  4.324281e-02  2.919151e-02  1.269126e-02  8.828116e-03  1.020002e-02  1.406360e-02  1.708787e-02  2.582078e-02  2.876437e-02  2.685574e-02  2.345379e-02  2.106032e-02  2.010601e-02  1.994271e-02  1.929665e-02  1.880384e-02  1.831317e-02  1.743624e-02
        '''.strip().split())
        import copy
        aTestMatrix = []  # copy.deepcopy( ndim*[ ndim*[ 0.0 ] ] )
        idata = 0
        for i in range(ndim):
            aTestMatrix.append([])
            for j in range(i + 1):
                aTestMatrix[i].append(the_data[idata])
                idata += 1
                if i != j:
                    aTestMatrix[j].append(aTestMatrix[i][j])
        return numpy.asmatrix(aTestMatrix)


def get_covariances_from_endf(endfFile, MT, MF=33):
    from brownies.legacy.converting import endfFileToGNDS

    rce = endfFileToGNDS.endfFileToGNDS(endfFile, toStdOut=False)
    xFileOne, cFileOne = rce['reactionSuite'], rce['covarianceSuite']
    MFMTListOne = [section.rowData.attributes['ENDF_MFMT'] for section in cFileOne]
    if not str(MF) + "," + str(MT) in MFMTListOne:
        raise ValueError("Requested MF,MT (" + str(MF) + "," + str(MT) +
                         ") not in first file, pick from " + str(MFMTListOne))
    try:
        for section in cFileOne:
            if str(MF) + ',' + str(MT) == section.rowData.attributes['ENDF_MFMT'] and \
                    section.nativeData in ['covarianceMatrix', 'mixed']:
                sectionOne = section
                break
        if sectionOne is None:
            raise KeyError('File ' + endfFile + ' missing plain old covariance matrix for MF,MT=' + str(MF) + ',' +
                           str(MT))
    except KeyError as err:
        print(err.message)
    # Sometimes the matrices are plain matrices, other times they are broken in components.
    # We'll make all of them a list of components
    componentListOne = []
    if sectionOne.nativeData == 'covarianceMatrix':
        componentListOne.append(sectionOne.forms['covarianceMatrix'])
    else:
        for component in sectionOne.forms['mixed'].components:
            componentListOne.append(component)
    return [numpy.array(x.matrix.data) for x in componentListOne]


if __name__ == "__main__":
    # ------------------------- command line parser -------------------------
    parser = argparse.ArgumentParser(description='Covariance test widget')
    parser.add_argument('--mt', dest='mt', type=int, default=53, help='MT of the covariance to check')
    parser.add_argument('--endf', dest='endf', type=str, default=None,
                        help='The endf file whose covariance you want to use')
    parser.add_argument('--mf', dest='mf', type=int, default=33,
                        help='MF of orginal data to use [default is 33, cross section covariance]')
    parser.add_argument('--tests', dest='doTests', default=False, action='store_true', help='Run unit tests')
    parser.add_argument('--plot', dest='plot', default=False, action='store_true',
                        help='Make a plot of the covariance used')
    parser.add_argument('--plotDiff', dest='plotDiff', default=False, action='store_true',
                        help='Make a plot of the differences between old & new covariances if applying fix')
    parser.add_argument('--plotBadSpace', dest='plotBadSpace', default=False, action='store_true',
                        help='Make a plot of bad eigenspaces')
    parser.add_argument('--print', dest='print', default=False, action='store_true',
                        help='Print out the entire covariance')
    parser.add_argument('--printEVs', dest='printEVs', default=False, action='store_true',
                        help='Print out the eigen values of the covariance')
    parser.add_argument('--printType', dest='printType', default=False, action='store_true',
                        help='Print type information for the covariance')
    parser.add_argument('--niter', dest='niter', default=0, type=int, help='Num iterations')
    parser.add_argument('--scaleOffDiagonals', dest='scaleOffDiagonals', default=False, action='store_true',
                        help='Scale off diagonal elements by some factor set by --scaleFactor switch')
    parser.add_argument('--scaleMatrix', dest='scaleMatrix', default=False, action='store_true',
                        help='Scale all elements by some factor set by --scaleFactor switch')
    parser.add_argument('--scaleFactor', dest='scaleFactor', default=0.999999, type=float,
                        help='Scale factor to use when scaling off diagonal elements (Default: 0.999999)')
    parser.add_argument('--pruneEVs', dest='pruneEVs', default=False, action='store_true',
                        help='Prune small/negative eigenvalues/eigenspaces')
    args = parser.parse_args()

    numpy.set_printoptions(linewidth=100)

    # Run the units tests
    if args.doTests:
        import sys
        import subprocess

        python = sys.executable
        try:
            subprocess.check_call([python, 'test_mtx.py'])
        except subprocess.CalledProcessError:
            pass

    # The test matrices
    if args.endf is not None:
        listOfA = get_covariances_from_endf(args.endf, args.mt, args.mf)
    else:
        listOfA = [get_test_matrix()]

    for icomponent, A in enumerate(listOfA):
        print()

        print("Component " + str(icomponent) + ":")

        # Check the matrix
        success = \
            check_symmetric(A) and \
            check_covariance_element_bounds(A) and \
            check_positive_semidefinite(A)

        # Any and all matrix displaying
        if args.printType:
            print(type(A))
        if args.print:
            print(A)
        if args.plot:
            plot_matrix(A, title='MT=' + str(args.mt) + ', Covariance Matrix')
        if args.printEVs:
            print_eigenvalues(A)
            plot_bad_eigenspaces(A)
        Aold = copy.copy(A)

        for i in range(args.niter):
            if success:
                print("After " + str(i) + " iterations, component " + str(icomponent) + " is OK")
                break

            # Iterating fixes
            if args.scaleOffDiagonals:
                Anew = scale_off_diagonals(Aold, scaleFactor=args.scaleFactor)
            if args.scaleMatrix:
                Anew = args.scaleFactor * Aold
            if args.pruneEVs:
                Anew = pruned_matrix(Aold)

            # Recheck everything
            success = \
                check_symmetric(Anew) and \
                check_covariance_element_bounds(Anew) and \
                check_positive_semidefinite(Anew)

            # Any and all matrix displaying (again)
            if args.print:
                print(Anew)
            if args.printType:
                print(type(Anew))
            if args.plot:
                plot_matrix(Anew, title='MT=' + str(args.mt) + ', Covariance Matrix')
            if args.printEVs:
                print_eigenvalues(Anew)
            if args.plotDiff:
                diff, reldiff = diff_matrices(A, Anew)
                plot_matrix(diff, title='Absolute difference')
                plot_matrix(reldiff, title='Relative % difference')
            Aold = Anew

    # In numpy, have around( obj[, decimals[, out] ), _round() too
    # numpy.testing.assert_approx_equal(actual, desired, significant=7, err_msg='', verbose=True)
    # numpy.nan_to_num( diff/numpy.matrix( matrixOne ) )
