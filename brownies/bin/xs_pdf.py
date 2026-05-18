#! /usr/bin/env python

import argparse, math, json, collections
import numpy as np
import scipy.special as sp

from brownies.BNL.plot_evaluation import getEXFORSets, generatePlot, getXdXYdYDataSets
from fudge.processing.resonances.reconstructResonances import URRPDFTable
from fudge.core.utilities.brb import banner
from fudge.reactionData.crossSection import XYs1d, upperEps

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule
from xData import axes as axesModule
from pqu import PQU

from brownies.legacy.endl import misc as miscENDLModule

# If numba available, use it to get some speedups
if False:
    try:
        from numba import jit, float64, int32
    except ImportError:
        def jit(func, signature=None, nopython=False, nogil=False, cache=False, forceobj=False, locals={}):
            return func
else:
    def jit(func, signature=None, nopython=False, nogil=False, cache=False, forceobj=False, locals={}):
        return func

SQRT2 = math.sqrt(2.0)
defaultEmin = 1e-11  # in MeV
defaultEmax = 20.  # in MeV
defaultXSmin = None  # if None, take from data
defaultXSmax = None  # if None, take from data
defaultNXS = 20
defaultNE = 100
exforRxnMap = {1: 'N,TOT', 2: 'N,EL', 102: 'N,G', 18: 'N,F'}
gndsRxnMap = {1: 'total', 2: 'elastic', 102: 'capture', 18: 'fission'}


# @jit
def equal_lethargy_bins(numBins, domainMin=defaultEmin, domainMax=defaultEmax):
    seq = np.logspace(start=math.log10(domainMin), stop=math.log10(domainMax), num=numBins + 1)
    return seq


# @jit
def equal_bins(numBins, domainMin=defaultEmin, domainMax=defaultEmax):
    seq = np.linspace(start=domainMin, stop=domainMax, num=numBins + 1)
    return seq


# @jit
def normal_dist_pdf(x, x0, dx):
    dx2 = dx * dx
    xx = x - x0
    xx2 = xx * xx
    return math.exp(-xx2 / 2. / dx2) / math.sqrt(2. * math.pi * dx2)


# @jit
def normal_dist_cdf(x, x0, dx):
    return 0.5 * (1. + sp.erf((x - x0) / SQRT2 / dx))


# @jit
def getPDFFromData(MT, _x4Sets, _eBins, _xsBins, skipEXFOREntries=[]):
    # The results histogram
    _binnedPDF = URRPDFTable(eBins=_eBins, xsBins=_xsBins)
    _binnedPDF[gndsRxnMap[MT]] = np.zeros((_binnedPDF.NE, _binnedPDF.NXS))

    # Create the binned PDF
    skipIt = False
    for x4Set in _x4Sets:
        for ent in skipEXFOREntries:
            if ent in x4Set.legend:
                skipIt = True
                break
        if skipIt: continue
        nPts = len(x4Set.x)
        for i in range(nPts):
            # Skip points that are off the edge of the grid
            if x4Set.x[i] >= _binnedPDF.Emax: continue
            if x4Set.x[i] <= _binnedPDF.Emin: continue
            if x4Set.y[i] >= _binnedPDF.XSmax: continue
            if x4Set.y[i] <= _binnedPDF.XSmin: continue
            if x4Set.yerr[i] == 0.0: continue
            iE = np.digitize(x4Set.x[i], _binnedPDF.eBins) - 1
            # This is stupid simple.  We made bins in the energy direction (eBins) and the
            # cross section direction (xsBins).  We then histogram the data in this 2D grid.
            # To turn it into the PDF P(xs|E), we divide all xs bins with the same energy bin
            # by the total number of points in that energy bin.
            if args.ignoreUncertainty:
                iXS = np.digitize(x4Set.y[i], _binnedPDF.xsBins) - 1
                binWeight = 1.  # / ptsInEBin[iE]
                _binnedPDF[gndsRxnMap[MT]][iE, iXS] += binWeight
            # This is more complex as we are including the experimental uncertainty in the
            # PDF generation.  As before we histogram the PDF P(xs|E), but because of the
            # experimental uncertainty, we don't know what cross section bin to add a point.
            # So, based on the experimental uncertainty we compute the probability that
            # the measured value is in each cross section bin and use that as the bin weight.
            else:
                for iXS in range(_binnedPDF.NXS):
                    binWeight = normal_dist_cdf(_binnedPDF.xsBins[iXS + 1], x4Set.y[i], x4Set.yerr[i]) - \
                                normal_dist_cdf(_binnedPDF.xsBins[iXS], x4Set.y[i], x4Set.yerr[i])
                    _binnedPDF[gndsRxnMap[MT]][iE, iXS] += binWeight
    _binnedPDF.normalize()
    return _binnedPDF


# @jit #FIXME: can't be jitted as is
def getPDFFromEvaluation(MT, _endfXS, _eBins, _xsBins, thicken=0):
    if thicken > 1: _endfXS = _endfXS.thicken(thicken)

    # The results histogram
    _binnedPDF = URRPDFTable(eBins=_eBins, xsBins=_xsBins)
    _binnedPDF[gndsRxnMap[MT]] = np.zeros((_binnedPDF.NE, _binnedPDF.NXS))
    # if args.verbose: print('binnedPDF before loading:', _binnedPDF[gndsRxnMap[MT]])

    # Build the indicator functions.  These functions work this way:
    # at a given energy, if the cross section value is within cross section
    # band iXS (that is, between _xsBins[iXS] and _xsBins[iXS+1]), than the
    # indicator function is 1.0 at that energy, otherwise it is zero.
    _junkXS = _endfXS.domainSlice(domainMin=_binnedPDF.eBins[0] * 1e6, domainMax=_binnedPDF.eBins[-1] * 1e6)
    indicatorFuncs = [_junkXS.copy() for iXS in range(_binnedPDF.NXS)]
    for iep, ep in enumerate(_junkXS):
        if ep[1] >= _binnedPDF.XSmax or ep[1] <= _binnedPDF.XSmin:
            for jXS, eta in enumerate(indicatorFuncs):
                eta[iep] = [ep[0], 0.0]
        iXS = np.digitize(ep[1], _binnedPDF.xsBins) - 1
        for jXS, eta in enumerate(indicatorFuncs):
            if jXS == iXS:
                eta[iep] = [ep[0], 1.0]
            else:
                eta[iep] = [ep[0], 0.0]

    # Group the indicator functions and read off the probabilities
    # If you group an indicator function with dx weighting, then you have
    # the probability that for a random energy in group iE, if you read the
    # cross section off, the grouped indicator function is the probability
    # that the cross section is in band iXS.
    for iXS, eta in enumerate(indicatorFuncs):
        geta = eta.group([E * 1e6 for E in _eBins], norm='dx')
        for iE, p in enumerate(geta):
            _binnedPDF[gndsRxnMap[MT]][iE, iXS] += p

    _binnedPDF.normalize()
    return _binnedPDF


def getPDFFromPURROutput(purrFile, energyBin, verbose=False):
    if energyBin is None: 
        energyBin = 0
    lines = open(purrFile).readlines()
    line = ''
    nTemps = 0
    nSigma0 = 0
    nLadders = 0
    NXS = 0
    NE = 0
    eBins = []
    xsBins = []
    probs = []
    xs = {}

    # advance to start of PURR section
    while lines and not line.startswith(" purr...probabalistic unresolved calculation"):
        line = lines.pop(0)
        if verbose: print(line)

    # do the real reading
    while lines and not line.startswith(" ******************"):

        # advance to next energy
        line = lines.pop(0)
        if verbose: print(line)
        if line.startswith(' temperatures'):
            while line.startswith('    ') or line.startswith(' temperatures'):
                nTemps += 1
                line = lines.pop(0)
            print("nTemps =", nTemps)
        if line.startswith(' sigma zero values'):
            while line.startswith('    ') or line.startswith(' sigma zero values'):
                nSigma0 += 1
                line = lines.pop(0)
            print("nSigma0 =", nSigma0)
        if line.startswith(' number of probability bins'):
            NXS = int(line.split()[-1])
            print("NXS =", NXS)
        if line.startswith(' number of resonance ladders'):
            nLadders = int(line.split()[-1])
            print("nLadders =", nLadders)
        if line.startswith(' no. of energy points (0=all)'):
            NE = int(line.split()[-1])
            print("NE =", NE)

        # now process the energy
        if line.startswith('e='):
            eBins.append(float(line.split()[1]) / 1e6)

            # we can only do one energy currently since PURR/NJOY uses a different XS grid for each incident energy
            if NE == energyBin:

                # skip ahead to the PURR table
                while lines and not line.startswith(" probability table"): 
                    line = lines.pop(0)

                # at start of table
                while lines and not line == '\n':
                    # get xs bins
                    if line.startswith(' tmin'):
                        xsBins.append(float(line.split()[-1]))
                        line = lines.pop(0)
                    elif line.startswith(' tmax'):
                        while line.startswith('    ') or line.startswith(' tmax'):
                            if line.startswith(' tmax'):
                                xsBins += [float(x) for x in line.split()[2:]]
                            else:
                                xsBins += [float(x) for x in line.split()]
                            line = lines.pop(0)
                        NXS = len(xsBins) - 1
                        print("NXS =", NXS)
                        print('xsBins:', xsBins)

                    # read in the total prob table
                    elif line.startswith(' prob'):
                        while line.startswith('    ') or line.startswith(' prob'):
                            if line.startswith(' prob'):
                                probs += [float(x) for x in line.split()[2:]]
                            else:
                                probs += [float(x) for x in line.split()]
                            line = lines.pop(0)

                    # read in the ave xs in the bin
                    elif line.startswith(' tot') or line.startswith(' els') or line.startswith(
                            ' fis') or line.startswith(' cap'):
                        rxn = line.split()[0]
                        xs[rxn] = []
                        while line.startswith('    ') or line.startswith(' ' + rxn):
                            if line.startswith(' ' + rxn):
                                xs[rxn] += [float(x) for x in line.split()[2:]]
                            else:
                                xs[rxn] += [float(x) for x in line.split()]
                            line = lines.pop(0)

                    # nothing... keep going
                    else:
                        line = lines.pop(0)

            NE += 1
    if verbose: 
        print('probs:', probs)
    if verbose: 
        print(list(xs.keys()))
    print("final NE =", NE)
    print('eBins:', eBins)
    urrPDF = URRPDFTable(eBins=np.array(eBins[energyBin:energyBin + 2]), xsBins=np.array(xsBins))
    urrPDF.NSamples = nLadders
    rxnMap = {'cap': 'capture', 'tot': 'total', 'els': 'elastic', 'fis': 'fission'}
    for rxn in xs:
        urrPDF[rxnMap[rxn]] = np.zeros(shape=(1, NXS))
        for iXS in range(NXS):
            urrPDF[rxnMap[rxn]][0, iXS] = probs[iXS] * xs[rxn][iXS] / xs['tot'][iXS]
    print(urrPDF)
    return urrPDF


def avexs_from_xspdf(key, xspdf, DEBUG=False, verbose=True):
    """

    :param key:
    :param xspdf:
    :param DEBUG:
    :param verbose:
    :return: a nested list composed of: [[E0, dE0, XS0, dXS0, skew0, kurtosis0], ...]
             E and dE are in MeV, XS and dXS are in b, skew and kurtosis are dimensionless
    """
    if not isinstance(xspdf, URRPDFTable): raise TypeError("xspdf should be a URRPDFTable, got %s" % str(type(xspdf)))
    ____EdEXSdXS = __avexs_from_xspdf(xspdf[key], xspdf.NE, xspdf.eBinCenters, xspdf.eBinWidths,
                                      xspdf.NXS, xspdf.xsBinCenters, xspdf.xsBinWidths)
    if verbose: print(xspdf.NE, xspdf.NXS)
    if DEBUG:
        for iE in range(xspdf.NE):
            if verbose: print('   ', iE)
            junk = xspdf.pdf_at_energy(key, xspdf.eBinCenters[iE])
            print("test mean:", ____EdEXSdXS[iE, 2], junk.mean)
            print("test stddev:", ____EdEXSdXS[iE, 3], junk.stddev)
            print("test skew:", ____EdEXSdXS[iE, 4], junk.skew)
            print("test kurtosis:", ____EdEXSdXS[iE, 5], junk.kurtosis)
    return ____EdEXSdXS


# @jit(nopython=True)
def __avexs_from_xspdf(xsPdf, NE, eBinCenters, eBinWidths, NXS, xsBinCenters, xsBinWidths):
    __EdEXSdXS = np.zeros(shape=(NE, 14))

    __EdEXSdXS[:, 0] = eBinCenters[:]
    __EdEXSdXS[:, 1] = eBinWidths[:] / 2.0
    for iE in np.arange(NE):

        # First pass gets the average
        for iXS in np.arange(NXS):
            __EdEXSdXS[iE, 2] += xsPdf[iE, iXS] * xsBinCenters[iXS] * xsBinWidths[iXS]

        # Second pass gets the standard deviation
        for iXS in np.arange(NXS):
            __EdEXSdXS[iE, 3] += xsPdf[iE, iXS] * xsBinWidths[iXS] * \
                                 (pow(xsBinCenters[iXS], 2.0) + pow(xsBinWidths[iXS], 3.0) / 12.)
        __EdEXSdXS[iE, 3] = math.sqrt(max(0.0, __EdEXSdXS[iE, 3] - pow(__EdEXSdXS[iE, 2], 2.0)))

        # Third pass gets the skew & kurtosis
        if __EdEXSdXS[iE, 3] > 0.0:
            for iXS in np.arange(NXS):
                __EdEXSdXS[iE, 4] += xsPdf[iE, iXS] * pow((xsBinCenters[iXS] - __EdEXSdXS[iE, 2]) / __EdEXSdXS[iE, 3],
                                                          3.0) * xsBinWidths[iXS]
                __EdEXSdXS[iE, 5] += xsPdf[iE, iXS] * pow((xsBinCenters[iXS] - __EdEXSdXS[iE, 2]) / __EdEXSdXS[iE, 3],
                                                          4.0) * xsBinWidths[iXS]

        if __EdEXSdXS[iE, 2] == 0.0: 
            continue

        # Add on log-normal columns
        mu = math.log(__EdEXSdXS[iE, 2] / math.sqrt(1.0 + pow(__EdEXSdXS[iE, 3] / __EdEXSdXS[iE, 2], 2.0)))
        sig = math.sqrt(math.log(1.0 + pow(__EdEXSdXS[iE, 3] / __EdEXSdXS[iE, 2], 2.0)))
        __EdEXSdXS[iE, 6] = mu
        __EdEXSdXS[iE, 7] = sig
        __EdEXSdXS[iE, 8] = (math.exp(sig * sig) + 2.0) * math.sqrt(math.exp(sig * sig) - 1.0)
        __EdEXSdXS[iE, 9] = math.exp(4. * sig * sig) + 2.0 * math.exp(3. * sig * sig) + 3. * math.exp(
            2. * sig * sig) - 3.0

        # Add on Poisson columns
        __EdEXSdXS[iE, 10] = __EdEXSdXS[iE, 2]  # lambda=ave
        __EdEXSdXS[iE, 11] = __EdEXSdXS[iE, 2]  # variance=lambda
        __EdEXSdXS[iE, 12] = 1.0 / math.sqrt(__EdEXSdXS[iE, 2])  # skewness
        __EdEXSdXS[iE, 13] = 1.0 / __EdEXSdXS[iE, 2] + 3.0  # kurtosis
    return __EdEXSdXS


def function_to_XYs(func, fpars,
                    Egrid=equal_bins(100),
                    domainUnit='eV', domainName='energy_in', rangeUnit='b', rangeName='crossSection',
                    accuracy=upperEps):
    """
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    """
    return XYs1dModule.XYs1d.createFromFunction(
        XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex: (rangeName, rangeUnit),
            XYs1dModule.xAxisIndex: (domainName, domainUnit)}),
        Xs=Egrid,
        func=func,
        parameters=fpars,
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def grouped_values_to_XYs(groupBdries, valueList,
                          domainUnit='eV', domainName='energy_in', rangeUnit='b', rangeName='crossSection',
                          accuracy=upperEps, domainMin=1e-5, domainMax=20.0):
    if len(groupBdries) != len(valueList) + 1:
        raise ValueError("Group boundries and value lists have incompatable lengths: len(bdries)=%i, len(vals)=%i" %
                         (len(groupBdries), len(valueList)))
    return XYs1dModule.XYs1d(
        data=[groupBdries, [valueList[0]] + valueList],
        dataForm="xsandys",
        interpolation=xDataEnumsModule.Interpolation.flat,
        axes=XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex: (rangeName, rangeUnit),
            XYs1dModule.xAxisIndex: (domainName, domainUnit)}))

def fit_XY_points(_xdata, _ydata, _yerror, form='const', convertToXYs=True):
    # Using LSQR fitting to data
    import scipy.optimize as opt
    allowedForms = ['quad', 'linear', 'const']
    if form not in allowedForms: raise ValueError('form must be one of allowed forms: ' + str(allowedForms))
    # Define the fitting function
    if form == 'quad':
        p = [0.0, 0.0, 3.0]

        def f(x, *params):
            xx = x - args.Emin
            return xx * xx * params[0] + x * params[1] + params[2]

    elif form == 'linear':
        p = [0.0, 3.0]

        def f(x, *params):
            xx = x - args.Emin
            return xx * params[0] + params[1]

    else:  # const
        p = [3.0]

        def f(x, *params):
            return params[0]

    # ... and fit!
    popt, pcov = opt.curve_fit(f, _xdata, _ydata, p0=p, sigma=_yerror)
    pars = popt.tolist()
    pars.reverse()
    dpars = [math.sqrt(pcov[i, i]) for i in range(len(pars))]
    dpars.reverse()
    if convertToXYs:
        return function_to_XYs(lambda x, y: f(x, *(popt.tolist())), [])  # awkward conversion if I ever saw one...
    return popt, pcov


def fit_series(xdata, ydata, yerr, seriesInstance, priorModel=None, priorCov=None):
    import fudge.core.math.linearAlgebra, numpy.matlib
    nData = len(ydata)
    nCoeff = len(seriesInstance.coefficients)
    if nData > 1000:
        raise ValueError("Can't handle nData > 1000")
    kernel = numpy.matlib.zeros([nData, nCoeff])
    for ic in range(nCoeff):
        for ix, x in enumerate(xdata):
            kernel[ix, ic] = seriesInstance.evaluateBasisFunction(x, ic)
    try:
        fit, cov, resid, chi2 = fudge.core.math.linearAlgebra.cglsqrSolve(
            data=numpy.asmatrix(ydata),
            dataUnc=numpy.asmatrix(yerr),
            kernel=kernel,
            prior=priorModel,
            priorCov=priorCov)
    except numpy.linalg.linalg.LinAlgError as theError:
        print('kernel', kernel)
        print('data', numpy.asmatrix(ydata))
        print('unc', numpy.asmatrix(yerr))
        raise theError
    results = {'fit': fit, 'cov': cov, 'residual': resid, 'chi2': chi2}
    return results


def get_chunk_list(N, maxChunkSize=1000):
    if N < maxChunkSize: return [(0, N)]
    nChunks = N // maxChunkSize
    leftOvers = N % maxChunkSize
    chunks = []
    for i in range(nChunks):
        chunkStart = i * maxChunkSize
        chunkStop = (i + 1) * maxChunkSize
        chunks.append((chunkStart, chunkStop))
    if leftOvers > 0: chunks.append((chunks[-1][-1], chunks[-1][-1] + leftOvers))
    return chunks


def process_args():
    # Set up command line parser
    parser = argparse.ArgumentParser(description='Get the cross section PDF given experimental data or evaluation.')

    # What to plot
    parser.add_argument('MT', type=int, help='The MT of the reaction.')

    # Decide the source of the pdf to work on
    parser.add_argument('--source', default='data', choices=['data', 'eval', 'file'],
                        help="Set where PDF comes from: data, evaluation or pre-computed file")

    # File names for ENDF data or raw URR PDF
    parser.add_argument('--endf', type=str, default=None, help="The ENDF file to use (if given)")
    parser.add_argument('--pdf', type=str, default=None, help="The file with the URR pdf (if given)")
    parser.add_argument('--pdfFormat', type=str, default='default', choices=['PURR', 'default'],
                        help="The format of file with the URR pdf (if given)")

    # Decide the nuclide to process if an ENDF file isn't given
    parser.add_argument('--sym', type=str, default=None, help='The symbol of the nucleus in question')
    parser.add_argument('--A', type=str, default=None, help='The A of the nucleus in question')

    # Cross section grid
    parser.add_argument('--XSmin', type=float, default=defaultXSmin, help="XSmin in b (default is to take from data)")
    parser.add_argument('--XSmax', type=float, default=defaultXSmax, help="XSmax in b (default is to take from data)")
    parser.add_argument('--NXS', type=int, default=defaultNXS,
                        help="Number of cross section bins (default=%i)" % defaultNXS)
    parser.add_argument('--useXSLethargy', default=False, action='store_true',
                        help='Make the cross section bins equal lethargy bins (default: False)')

    # Energy grid
    parser.add_argument('--Emin', type=float, default=defaultEmin,
                        help="Emin in MeV (default=%s MeV)" % str(defaultEmin))
    parser.add_argument('--Emax', type=float, default=defaultEmax,
                        help="Emax in MeV (default=%s MeV)" % str(defaultEmax))
    parser.add_argument('--NE', type=int, default=defaultNE,
                        help="Number of equal bins on the domain (Emin, Emax) (default=%i)" % defaultNE)
    parser.add_argument('--useELethargy', default=False, action='store_true',
                        help='Make the energy bins equal lethargy bins (default: False)')

    # Treatment of EXFOR data
    parser.add_argument('--ignoreData', default=False, action='store_true',
                        help="Ignore the experimental data unless it's what's getting averaged")
    parser.add_argument('--ignoreUncertainty', default=False, action='store_true',
                        help="Ignore the uncertainty on the experimental data")
    parser.add_argument('--printSetStats', default=False, action='store_true', help="Print EXFOR data set statistics")
    parser.add_argument('--saveSetStats', type=str, default=None, help="Store set information to this JSON file")
    parser.add_argument('--readSetModifications', type=str, default=None,
                        help='JSON file with potential file modifications')
    parser.add_argument('--editSets', default=False, action="store_true",
                        help="Edit the sets so that they are between Emin & Emax, etc")

    # Treatment of ENDF data
    parser.add_argument('--skipBadData', default=False, action='store_true',
                        help="Skip bad data in the ENDF file if needed")
    parser.add_argument('--continuumSpectraFix', default=False, action='store_true',
                        help="Apply the continuum spectrum fix to ENDF data if needed")
    parser.add_argument('--thicken', default=0, type=int,
                        help="Thicken (reconstructed) cross sections by this number of points between existing points "
                             "(default = 0, i.e. don't thicken")

    # Treatment of PDFs read from a file
    parser.add_argument('--noRenormalize', default=False, action="store_true",
                        help="Do not renormalize the PDFs read in from a file")
    parser.add_argument('--check', default=False, action="store_true",
                        help="Print out the PDf normalization inegral as a cross check and do other checks")

    # What to compute
    parser.add_argument('--getAveXS', default=False, action="store_true", help="Compute the average cross section")
    parser.add_argument('--aveXSMode', default="pdf", choices=['quad', 'linear', 'const', 'pdf', 'endf', 'mughabghab'],
                        help="Mode to compute the average")
    parser.add_argument('--getXSPDF', default=False, action="store_true",
                        help="Compute the cross section probability distribution")
    parser.add_argument('--getXSCorr', default=False, action="store_true",
                        help="Compute the cross section energy-energy correlation")

    # What is the form of the output?
    parser.add_argument('--aveOutFile', default=None, type=str,
                        help='The name of the output file for the average cross section')
    parser.add_argument('--pdfOutFile', default=None, type=str,
                        help='The name of the output file for the cross section PDF')
    parser.add_argument('--printXSPDFSlice', default=None, type=int,
                        help='Print a slice of the cross section PDF derived from the data in the bin specified')
    parser.add_argument('--plotXSCorr', default=False, action="store_true",
                        help="Plot the cross section energy-energy correlation")
    parser.add_argument('--printXSCorr', default=False, action="store_true",
                        help="Print the cross section energy-energy correlation")
    parser.add_argument('--plotXSPDF', default=False, action='store_true',
                        help='Plot cross section PDF derived from the data')
    parser.add_argument('--plotXSPDFSlice', default=None, type=int,
                        help='Plot a slice of the cross section PDF derived from the data in the bin specified')
    parser.add_argument('--plotAveXS', default=False, action='store_true',
                        help='Plot experimental data and the average cross section')

    # Other flags
    parser.add_argument('-v', dest='verbose', default=False, action='store_true',
                        help='Enable verbose output')

    return parser.parse_args()


if __name__ == "__main__":
    # Process command line options
    args = process_args()

    # Figure out sym & A
    if args.sym is None or args.A is None:
        if args.endf is None:
            print("\nERROR: No nucleus specified.  Either enter a symbol and A or enter an ENDF evaluation\n")
            exit()

    # Read the ENDF file
    endfEval = None
    endfXS = None
    aveEndfXSOverDomain = None
    if args.endf is not None:
        print(banner('Reading ENDF file ' + args.endf))
        from brownies.BNL.plot_evaluation.plotio import readEvaluation

        endfEval = \
            readEvaluation(args.endf, skipBadData=args.skipBadData, continuumSpectraFix=args.continuumSpectraFix)[0]
        targ = str(endfEval.target)
        sym, A, m = miscENDLModule.elementAFromName(targ)
        args.sym = sym
        args.A = A
        if m is not None:
            args.A += str(m)
        endfXS = endfEval.getReaction(args.MT) \
            .crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8) \
            .domainSlice(domainMin=args.Emin * 1e6, domainMax=args.Emax * 1e6)
        aveEndfXSOverDomain = endfXS.integrate() / float(PQU.PQU(endfXS.domainMax - endfXS.domainMin, 'eV'))
        print()

    # Get the EXFOR data
    x4Sets = []
    if not args.ignoreData or args.source == 'data':
        x4Sets = getEXFORSets(args.sym, args.A, reaction=exforRxnMap[args.MT], quantity="SIG",
                              nox4evals=True, nox4legend=False, forceLegend=True,
                              plotSyle={}, verbose=False)

        # Figure out if there are changes to any of the sets (particularly to exclude them)
        localMods = {}
        if args.readSetModifications is not None:
            if args.verbose:
                print("# Reading local modifications from %s" % args.readSetModifications)
            localMods = json.loads(open(args.readSetModifications, mode='r').read())

        # Edit sets
        if args.editSets:
            trimmedSets = []
            for x4Set in x4Sets:
                # Get rid of excluded sets
                if x4Set.legend in localMods and localMods[x4Set.legend]['exclude']:
                    if args.verbose:
                        print("# Excluding %s" % x4Set.legend)
                    continue
                if x4Set.legend not in localMods:
                    print("# WARNING: Set with legend %s not in localMods file" % x4Set.legend)
                    continue
                # Get rid sets with no points inside the range of interest
                if min(x4Set.x) > args.Emax or max(x4Set.x) < args.Emin:
                    continue
                # Rescale points as needed
                if localMods[x4Set.legend]["xs norm. factor"] != 1.0:
                    x4Set.y = [y * localMods[x4Set.legend]["xs norm. factor"] for y in x4Set.y]
                    x4Set.yerr = [yerr * localMods[x4Set.legend]["xs norm. factor"] for yerr in x4Set.yerr]
                # Scale and shift energies as needed
                if localMods[x4Set.legend]["E scale"] != 1.0 or localMods[x4Set.legend]["E shift (MeV)"] != 0.0:
                    x4Set.x = [x * localMods[x4Set.legend]["E scale"] + localMods[x4Set.legend]["E shift (MeV)"] for x
                               in x4Set.x]
                    x4Set.xerr = [xerr * localMods[x4Set.legend]["E scale"] for xerr in x4Set.xerr]
                # Adjust legend if rescaling happened
                if localMods[x4Set.legend]["xs norm. factor"] != 1.0:
                    x4Set.legend += ' *' + str(localMods[x4Set.legend]["xs norm. factor"])
                trimmedSets.append(x4Set)
            x4Sets = trimmedSets

        nSets = len(x4Sets)
        # Information about the sets used
        if args.printSetStats or (args.saveSetStats is not None):
            setInfo = collections.OrderedDict()
            for x4Set in x4Sets:
                # Computed suggested xs normalization based on ENDF file
                xsNormFactor = 1.0
                exclude = False
                if True and aveEndfXSOverDomain is not None:
                    xx, yy, yyerr = [], [], []
                    for ipt in range(len(x4Set.x)):
                        if x4Set.x[ipt] > args.Emax or x4Set.x[ipt] < args.Emin:
                            continue
                        if x4Set.yerr[ipt] <= 0.0:
                            continue
                        xx.append(x4Set.x[ipt])
                        yy.append(x4Set.y[ipt])
                        yyerr.append(x4Set.yerr[ipt])
                    if len(xx) > 1:
                        a = fit_XY_points(xx, yy, yyerr, form='const', convertToXYs=False)
                        xsNormFactor = (aveEndfXSOverDomain / a[0][0]).getValue()
                        if abs(math.log(xsNormFactor)) >= 1.:
                            print("# Set %s requires an unreasonable renormalization (%s) the "
                                  "'exclude' flag is being set to True" % (x4Set.legend, str(xsNormFactor)))
                            exclude = True
                    else:
                        print("# Set %s has too few usuable points (%i) in the region of "
                              "interest so the 'exclude' flag is being set to True" % (x4Set.legend, len(xx)))
                        exclude = True
                setInfo[x4Set.legend] = {
                    'legend': x4Set.legend,
                    'num. pts.': len(x4Set.x),
                    'xs mean (b)': np.mean(x4Set.y),
                    'xs stdev (b)': np.std(x4Set.y),
                    'xs min (b)': min(x4Set.y),
                    'xs max (b)': max(x4Set.y),
                    'E min (MeV)': min(x4Set.x),
                    'E max (MeV)': max(x4Set.x),
                    'xs norm. factor': xsNormFactor,
                    'E shift (MeV)': 0.0,
                    'E scale': 1.0,
                    'exclude': exclude
                }
            # Update set info with modifications
            if args.readSetModifications is not None:
                for k in setInfo:
                    if k in localMods:
                        for kk in ['xs norm. factor', 'E shift (MeV)', 'E scale', 'exclude']:
                            setInfo[k][kk] = localMods[k][kk]

            # Print set info
            if args.printSetStats:
                print(banner("Set information"))
                print(json.dumps(setInfo, sort_keys=True, indent=4, separators=(',', ': ')))

            # Save set info in file
            if args.saveSetStats is not None:
                open(args.saveSetStats, mode='w').write(json.dumps(setInfo, sort_keys=True, indent=4))

        # Get the precomputed PDF
    urrPDF = URRPDFTable()
    if args.pdf is not None:
        print(banner("Reading PDF from %s" % args.pdf))
        if args.pdfFormat == 'default':
            urrPDF.load(args.pdf)
            urrPDF.eBins /= 1e6
        elif args.pdfFormat == 'PURR':
            sliceBin = None
            if args.plotPDFSlice is not None:
                sliceBin = args.plotPDFSlice
                args.plotPDFSlice = 0  # a kludge because I can't read the full PURR table all at once, so have one bin
            elif args.printPDFSlice is not None:
                sliceBin = args.printPDFSlice
                args.printPDFSlice = 0  # a kludge because I can't read the full PURR table all at once, so have one bin
            urrPDF = getPDFFromPURROutput(args.pdf, sliceBin)
        else:
            raise ValueError("Don't know PURR format '%s'" % args.pdfFormat)
        if args.noRenormalize:
            pass
        else:
            urrPDF.normalize()

    # Set up the grid
    if args.pdf is not None:
        eBins = urrPDF.eBins
        xsBins = urrPDF.xsBins
        args.XSmin = xsBins[0]
        args.XSmax = xsBins[-1]
        args.Emin = eBins[0]
        args.Emax = eBins[-1]
        args.NE = urrPDF.NE
        args.NXS = urrPDF.NXS
    else:
        # Set min/max values of the grids
        if args.XSmin is None:
            minList = [min(x4Set.y) for x4Set in x4Sets]
            if endfXS is not None:
                minList.append(endfXS.rangeMin)
            args.XSmin = max(min(minList), 0.0)
        if args.XSmax is None:
            maxList = [max(x4Set.y) for x4Set in x4Sets]
            if endfXS is not None:
                maxList.append(endfXS.rangeMax)
            args.XSmax = max(maxList)
        if args.useELethargy:
            eBins = equal_lethargy_bins(args.NE, domainMin=args.Emin, domainMax=args.Emax)
        else:
            eBins = equal_bins(args.NE, domainMin=args.Emin, domainMax=args.Emax)
        if args.useXSLethargy:
            xsBins = equal_lethargy_bins(args.NXS, domainMin=max(args.XSmin, 1e-15), domainMax=args.XSmax)
        else:
            xsBins = equal_bins(args.NXS, domainMin=args.XSmin, domainMax=args.XSmax)
    if args.verbose:
        print('energy bins:', eBins)
        print('cross section bins:', xsBins)

    # Compute the PDF
    if args.getXSPDF or (args.getAveXS and args.aveXSMode == 'pdf'):
        print(banner("Getting PDF from " + args.source))
        if args.source == 'data' and x4Sets != []:
            binnedPDF = getPDFFromData(args.MT, x4Sets, eBins, xsBins)
        elif args.source == 'eval' and endfXS is not None:
            binnedPDF = getPDFFromEvaluation(args.MT, endfXS, eBins, xsBins, args.thicken)
        elif args.source == 'file' and urrPDF is not None:
            binnedPDF = urrPDF
        else:
            exit()

    # Check the pdf normalization
    if args.check:
        print(banner("Checking PDF norm from " + args.source))
        for iE in range(binnedPDF.NE):
            norm = 0.0
            for iXS in range(binnedPDF.NXS):
                norm += binnedPDF[gndsRxnMap[args.MT]][iE, iXS] * binnedPDF.xsBinWidths[iXS]
            print(iE, binnedPDF.eBinCenters[iE], binnedPDF.eBinWidths[iE], norm)

    # Compute the average cross section
    # aveXSMode: default="pdf", choices=['quad', 'linear', 'const', 'endf']
    if args.getAveXS or args.getXSCorr:
        if args.aveXSMode == 'pdf':
            # Using the PDF
            EdEXSdXS = avexs_from_xspdf(gndsRxnMap[args.MT], binnedPDF, DEBUG=args.check)
            # aveXS = grouped_values_to_XYs(EdEXSdXS, popt) # FIXME
        elif args.aveXSMode == 'endf':
            aveXS = endfXS
        elif args.aveXSMode == "mughabghab":
            # Use Said Mughaghab's systematics
            raise NotImplementedError()
            # Said gives <\sigma_tot> in S.F. Mughabghab, Nucl. Data Sheets 118, 287-291 (2014), eqs. (1), (2)
            # Said gives <\sigma_g> in S.F. Muchabghab, Atlas of Neutron Resonances, 5th Ed., Elsevier (2006), eq.(2.19)
        else:
            # Build the list of data to fit
            xdata = []
            ydata = []
            yerror = []
            for x4Set in x4Sets:
                for ipt in range(len(x4Set.x)):
                    if x4Set.x[ipt] > args.Emax or x4Set.x[ipt] < args.Emin:
                        continue
                    if x4Set.yerr[ipt] <= 0.0:
                        continue
                    xdata.append(x4Set.x[ipt])
                    ydata.append(x4Set.y[ipt])
                    yerror.append(x4Set.yerr[ipt])
            # Now fit!
            aveXS = fit_XY_points(xdata, ydata, yerror, form=args.aveXSMode, convertToXYs=True)
            if False:
                aveXS.plot()

    # Compute the correlation function
    if args.getXSCorr:
        expR = collections.OrderedDict()
        expRTemp = collections.OrderedDict()
        maxPairs = 1000000000
        for x4Set in x4Sets:
            iPair = 0
            stopSet = False
            try:
                binSizeRatio = (x4Set.x[1] - x4Set.x[0]) / (args.Emax - args.Emin)
            except:
                continue
            expR[x4Set.legend] = []
            expRTemp[x4Set.legend] = collections.OrderedDict()

            # Build the integrand by looping over pairs
            for i in range(len(x4Set.x)):
                if stopSet:
                    break
                if x4Set.x[i] > args.Emax or x4Set.x[i] < args.Emin:
                    continue
                if x4Set.yerr[i] <= 0.0:
                    continue
                for j in range(i + 1, len(x4Set.x)):

                    # Basic point pair quality control
                    if x4Set.x[j] > args.Emax or x4Set.x[j] < args.Emin:
                        continue
                    if x4Set.yerr[j] <= 0.0:
                        continue
                    if x4Set.x[i] == x4Set.x[j]:
                        continue

                    # Compute the point pair quantities
                    yiAve = aveXS.evaluate(x4Set.x[i])
                    yjAve = aveXS.evaluate(x4Set.x[j])
                    EDiff = abs(x4Set.x[i] - x4Set.x[j])  # really epsilon
                    dEDiff = x4Set.xerr[i] + x4Set.xerr[j]
                    R = (x4Set.y[i] - yiAve) * (x4Set.y[j] - yjAve) / (yiAve * yjAve)
                    relerri = x4Set.yerr[i] / x4Set.y[i]
                    relerrj = x4Set.yerr[j] / x4Set.y[j]
                    dR = abs(R) * math.sqrt(relerri * relerri + relerrj * relerrj)
                    RR = R * R

                    # Log the results from this pair of points
                    epsilonKey = int(round(EDiff, 6) * 1e6)
                    if epsilonKey not in expRTemp[x4Set.legend]:
                        expRTemp[x4Set.legend][epsilonKey] = []
                    expRTemp[x4Set.legend][epsilonKey].append([EDiff, R, dEDiff, dR, RR])
                    iPair += 1
                    if iPair > maxPairs:
                        stopSet = True

            # Make a sorted list of epsilons so we can prepare the results dict
            epsilonList = sorted(expRTemp[x4Set.legend].keys())

            # Integrate to get the correlator
            minNpts = 10
            for epsilonKey in epsilonList:
                if False:
                    Npts = 1. / binSizeRatio  # Simpson's rule suggests this
                else:
                    Npts = len(expRTemp[x4Set.legend][epsilonKey])  # Monte-Carlo approach suggests this
                if Npts < minNpts:
                    continue
                epsilon = np.mean([x[0] for x in expRTemp[x4Set.legend][epsilonKey]])
                R = np.sum([x[1] for x in expRTemp[x4Set.legend][epsilonKey]]) / Npts
                dEpsilon = np.mean([x[2] for x in expRTemp[x4Set.legend][epsilonKey]])
                RR = np.sum([x[4] for x in expRTemp[x4Set.legend][epsilonKey]]) / Npts
                dR = math.sqrt(max((RR - R * R), 0.0) / Npts)
                expR[x4Set.legend].append([epsilon, R, dEpsilon, dR])

            # Print the correlator from this set
            saveXSCorr = True
            if len(expR[x4Set.legend]) > 2:
                print('    Writing R(eps) for ', x4Set.legend)
                if True:  # Latex summary table of data sets used
                    year = x4Set.legend.split(')')[0].split('(')[1].strip()
                    author = x4Set.legend.split(')')[1].split('(')[0].strip()
                    exfor = x4Set.legend.split('(')[-1].split(')')[0].strip()
                    xsScale = x4Set.legend.split('(')[-1].split(')')[-1].strip()
                    Npts = len(x4Set.x)
                    NptsInRange = len([x for x in x4Set.x if args.Emax >= x >= args.Emin])
                    Npairs = len(expR[x4Set.legend])
                    print(' & '.join([year, author, exfor, xsScale, str(Npts), str(NptsInRange), str(Npairs)]) + '\\\\')
                if saveXSCorr:
                    setOut = open(x4Set.legend.split('(')[-1].split(')')[0] + '_Reps.dat', mode='w')
                    setOut.write('#   ' + x4Set.legend + '\n')
                    for x in expR[x4Set.legend]:
                        setOut.write('  '.join([str(xx) for xx in x]) + '\n')
            else:
                print('    Did not write R(eps) for ', x4Set.legend,
                      'it has too few points (' + str(len(expR[x4Set.legend])) + ')')
                pass

        # Spline fit to correlation
        if True:
            from xData import series1d as series1dModule

            axes = axesModule.Axes(2)
            axes[0] = axesModule.Axis('R(epsilon)', 0, '')
            axes[1] = axesModule.Axis('epsilon', 1, 'MeV')
            nModel = 100  # args.NE
            epsilonMin = 0.0
            epsilonMax = 0.8 * (args.Emax - args.Emin)
            bins = equal_bins(nModel, domainMin=epsilonMin, domainMax=epsilonMax)
            Rfunc = series1dModule.LinearSpline1d(
                xdata=bins,
                ydata=[0.0 for i in range(nModel + 1)],
                axes=axes)
            priorModel = None
            priorCov = None

            # Chop expR into 500 point or less chunks
            maxChunkSize = 500
            expRChunk = collections.OrderedDict()
            for k in expR:
                #                if "Weston" in k: continue
                if len(expR[k]) <= 2:
                    continue  # not worth using in fit
                elif len(expR[k]) > maxChunkSize:
                    for ichunk, chunk in enumerate(get_chunk_list(len(expR[k]), maxChunkSize=maxChunkSize)):
                        kk = k + ', chunk #' + str(ichunk)
                        expRChunk[kk] = expR[k][chunk[0]:chunk[1]]
                else:
                    expRChunk[k] = expR[k]

            # Bayesian update
            cumChi2 = 0.
            cumNPts = 0
            for k in expRChunk:
                print('    Bayesian updating with set %s' % k)
                xdataR = []
                ydataR = []
                yerrR = []
                for x in expRChunk[k]:
                    if x[0] < epsilonMin or x[0] > epsilonMax:
                        continue  # off the grid, we can't use the point
                    xdataR.append(x[0])
                    ydataR.append(x[1])
                    yerrR.append(x[3])
                if len(xdataR) == 0:
                    continue
                results = fit_series(xdataR, ydataR, yerrR, seriesInstance=Rfunc, priorModel=priorModel,
                                     priorCov=priorCov)
                priorModel = results['fit']
                priorCov = results['cov']
                cumChi2 += results['chi2']
                cumNPts += len(ydataR)
                print('       ', 'Chi2=', results['chi2'][0, 0], 'Chi2/ndf=', cumChi2[0, 0] / (cumNPts - nModel))
            Rfunc.setData(priorModel.tolist()[0])

            # Plot the cross section correlation
            if args.plotXSCorr:
                Rfunc.toPointwise_withLinearXYs().plot()

            # Print the correlator
            for i in range(len(priorModel.tolist()[0])):
                print(i, bins[i],
                      priorModel[0, i],
                      math.sqrt(priorCov[i, i]))

    # Output the PDF
    if args.getXSPDF:

        # Save to file
        if args.pdfOutFile is not None and args.source != 'file':
            binnedPDF.save(args.pdfOutFile)

        # Plot the PDF of the cross section
        if args.plotXSPDF:
            from fudge.vis.matplotlib.plot_matrix import plot_matrix
            import matplotlib.pyplot as plt

            title = '$^{' + args.A + '}$' + args.sym + '(' + exforRxnMap[args.MT].lower() + ') cross section PDF'
            if args.ignoreUncertainty:
                title += ', ignoring experimental uncertainty'
            plot_matrix(binnedPDF[gndsRxnMap[args.MT]].T, energyBoundariesX=eBins, energyBoundariesY=xsBins,
                        title=title,
                        xyTitle=('energy (MeV)', 'cross section (b)'),
                        switchY=False)
            plt.show()

        # Plot a slice of the PDF of the cross section
        if args.plotXSPDFSlice is not None:
            E = binnedPDF.eBinCenters[args.plotXSPDFSlice]
            binnedPDF.pdf_at_energy(gndsRxnMap[args.MT], E).plot()

        # Plot a slice of the PDF of the cross section
        if args.printXSPDFSlice is not None:
            E = binnedPDF.eBinCenters[args.printXSPDFSlice]
            print('\n'.join(binnedPDF.pdf_at_energy(gndsRxnMap[args.MT], E).toXML_strList()))

    # Output the average
    if args.getAveXS:

        # Print to stdout or a file
        def formatField(x, outputFieldSize=20):
            return str(x).ljust(outputFieldSize)


        def formatLine(l):
            return '  '.join([formatField(x) for x in l])


        colHeaders = ["# iE", "E (MeV)", "dE (MeV)", "XS (b)", "dXS (b)", 'XS skew', 'XS kurtosis', 'LN mu', 'LN sigma',
                      'LN skew', 'LN kurtosis', 'Poi lambda', 'Poi ave', 'Poi skew', 'Poi kurtosis']
        if args.aveOutFile is not None:
            out = open(args.aveOutFile, mode='w')
            out.write(formatLine(colHeaders) + '\n')
            for iE in np.arange(args.NE):
                out.write(formatLine([iE] + list(EdEXSdXS[iE])) + '\n')
            out.close()
        else:
            if args.ignoreUncertainty:
                print(banner("average cross section, ignoring experimental uncertainty"))
            else:
                print(banner("average cross section"))
            print(formatLine(colHeaders))
            for iE in range(args.NE):
                print(formatLine([iE] + list(EdEXSdXS[iE])))

        # Plot the average cross section
        # It uses plot_evaluation coding, so it is kinda kludgy.
        if args.plotAveXS:
            title = '$^{' + args.A + '}$' + args.sym + '(' + exforRxnMap[args.MT].lower() + ')'
            if args.ignoreUncertainty:
                title += ', ignoring experimental uncertainty'
            theAve = getXdXYdYDataSets({'aveXS': EdEXSdXS}, {})
            theAve[0].legend = 'average cross section'
            theAve[0].symbol = 'o'
            theAve[0].symbolSyize = 10
            theAve[0].dataType = 'scatter'
            theAve[0].errorbarLineWidth = 2
            theAve[0].errorbarColor = 'k'
            theAve[0].lineWidth = 0
            theAve[0].xUnit = 'MeV'
            endfSet = []
            if endfXS is not None:
                from fudge.vis.matplotlib import DataSet2d

                endfSet.append(
                    DataSet2d(
                        endfXS,
                        legend="ENDF",
                        lineWidth=4,
                        lineStyle="-",
                        color='blue'))
            whatToPlot = endfSet + theAve
            if not args.ignoreData:
                whatToPlot += x4Sets
            generatePlot(
                'crossSection',
                whatToPlot,
                plotStyle={
                    'plotStyles': {'crossSection': {
                        'xAxis': {'min': args.Emin * 1e3, 'max': args.Emax * 1e3, 'unit': 'keV'},
                        'yAxis': {'min': args.XSmin, 'max': args.XSmax, 'unit': 'b'}}}},
                suggestTitle=title,
                figsize=(15, 6),
                suggestYLog=args.MT in (1, 2, 18, 102),
                suggestXLog=args.MT in (1, 2, 18, 102))
