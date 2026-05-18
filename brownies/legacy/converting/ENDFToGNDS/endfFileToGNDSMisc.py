# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Helper functions to streamline reading ENDF-6 formatted data
"""

import sys

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule
from xData import date as dateModule
from xData.Documentation import dates as datesModule
from xData.Documentation import author as authorModule

from fudge import reactionSuite as reactionSuiteModule

FUDGE_EPS = 1e-8


def parseENDFByMT_MF(fileName, stripMATMFMTCount=True, logFile=sys.stderr):
    """This function reads an endf file and creates a python dictionary with the keys being MT
    numbers and the values being the data for the MT key. For each MT, the data are stored in a python
    dictionary with the keys being MF numbers and the values being the data for the MF key."""

    MTDatas = {}
    with open(fileName) as fin:
        header = next(fin)[:-1]           # Occasionally, the first line does not have MT or MF values.
        secondLine = None
        for lineNum, line in enumerate(fin):
            line = line.rstrip()
            if secondLine is None:
                secondLine = line
            if line == '':
                break
            try:
                MF = int(line[70:72])
                MT = int(line[72:75])
            except Exception as ex:
                logFile.write("Cannot extract MF/MT from line #%s\n" % lineNum)
                raise ex
            if MT not in MTDatas:
                MTDatas[MT] = {}
            MTData = MTDatas[MT]
            if MF not in MTData:
                MTData[MF] = []
            if stripMATMFMTCount:
                MTData[MF].append(line[:66].rstrip())
            else:
                MTData[MF].append(line)

    try:
        MAT = int(secondLine[66:70])
    except ValueError:
        MAT = int(MTDatas[451][1][6][31:35])

    if 0 in MTDatas:
        del MTDatas[0]

    return header, MAT, MTDatas


def ENDFInterpolationToGNDS1d(interpolation):

    if interpolation == 1:
        return xDataEnumsModule.Interpolation.flat
    if interpolation == 2:
        return xDataEnumsModule.Interpolation.linlin
    if interpolation == 3:
        return xDataEnumsModule.Interpolation.linlog
    if interpolation == 4:
        return xDataEnumsModule.Interpolation.loglin
    if interpolation == 5:
        return xDataEnumsModule.Interpolation.loglog
    if interpolation == 6:
        return xDataEnumsModule.Interpolation.chargedParticle
    raise Exception('Unsupported 2d interpolation = %d' % interpolation)


def ENDFInterpolationToGNDS2plusd(interpolation):

    if    0 < interpolation < 10:
        interpolationQualifier = xDataEnumsModule.InterpolationQualifier.none
        normalInterpolation = interpolation
    elif 10 < interpolation < 20:
        interpolationQualifier = xDataEnumsModule.InterpolationQualifier.correspondingPoints
        normalInterpolation = interpolation - 10
    elif 20 < interpolation < 30:
        interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase
        normalInterpolation = interpolation - 20
    else:
        raise Exception('Unsupported 2+d interpolation = %d' % interpolation)
    return interpolationQualifier, ENDFInterpolationToGNDS1d(normalInterpolation)


def funkyFloatStringToFloat(index, s, logFile=sys.stderr):

    index6 = index % 6
    value = s[11 * index6: 11 * (index6 + 1)]
    try:
        f = float(value)
    except ValueError:
        value = value.lower().replace('d', '')
        value = value.replace(' ', '')
        i1 = value.find('.')
        if i1 != -1:
            j = value[i1:].find('+')
            if j == -1: j = value[i1:].find('-')
            if j == -1: raise Exception("Float value at index = %d (%d) is not funky enough <%s>" % (index, index6, s))
            value = value[:i1+j] + 'e' + value[i1+j:]
        try:
            f = float(value)
        except ValueError:
            if value == len(value) * ' ': return 0
            if logFile is not None:
                logFile.write('Could not convert value "%s" at index %s (%s) of string "%s"\n' %
                              (value, index6, index, s))
            raise
    return f


def sixFunkyFloatStringsToFloats(s, logFile=sys.stderr):

    return [funkyFloatStringToFloat(i1, s, logFile=logFile) for i1 in range(6)]


def sixFunkyFloatStringsToIntsAndFloats(s, intIndices=(), logFile=sys.stderr):

    values = sixFunkyFloatStringsToFloats(s, logFile=logFile)
    for i1 in intIndices: values[i1] = int(values[i1])
    return values


def nFunkyFloatStringsToFloats(n1, startLine, lines, dimension=1, logFile=sys.stderr):

    floats = []
    nd = dimension * n1
    if len(lines[startLine:]) < nd // 6:
        logFile.write("WARNING: trying to read %d values, but only %d remain in the section!\n" %
                      (nd, len(lines[startLine:]) * 6))
    while True:
        line = lines[startLine]
        n6 = min(nd, 6)
        for i1 in range(n6):
            floats.append(funkyFloatStringToFloat(i1, line, logFile=logFile))
        nd -= n6
        if nd == 0: break
        startLine += 1
    if dimension > 1:
        floats1D, floats = floats, []
        for i1 in range(0, dimension * n1, dimension):
            floats.append(floats1D[i1:i1+dimension])
    return floats


def readEndfINTG(string, ndigit):
    """Special case in ENDF: line of all integers, used to represent sparse correlation matrix."""

    ix = int(string[:5])
    iy = int(string[5:10])

    try:
        i1, i2 = {2: [11, 65], 3: [11, 63], 4: [11, 66], 5: [11, 65], 6: [10, 66]}[ndigit]
    except KeyError:
        raise ValueError("Encountered illegal ndigit (%d) for INTG format" % ndigit)
    sub = string[i1:i2]
    n1 = ndigit + 1
    vals = [sub[n1*i1:n1*(i1+1)] for i1 in range(len(sub) // n1)]

    def parseInt(intString):

        return 0 if not intString.strip() else int(intString)

    arr = [parseInt(a) for a in vals]
    return ix, iy, arr


def readDiscreteAndLegendre(numDiscrete, numContinua, startLine, lines, dimension, logFile=sys.stderr):

    floats = []
    nd = dimension * numDiscrete + dimension * numContinua   # total amount of data
    while True:
        line = lines[startLine]
        n6 = min(nd, 6)
        for i1 in range(n6):
            floats.append(funkyFloatStringToFloat(i1, line, logFile=logFile))
        nd -= n6
        if nd == 0: break
        startLine += 1

    discretes = []              # Data for discrete gammas is [ gammaEnergy, multiplicity ]
    duplicateCounter = 0
    EgPrior = None
    for i1 in range(0, dimension * numDiscrete, dimension):
        Eg = floats[i1]
        if i1 > 0:
            if Eg != EgPrior: duplicateCounter = 0
        duplicateCounter += 1
        discretes.append([Eg, duplicateCounter, floats[i1+1:i1+dimension]])
        EgPrior = Eg
        
    continua = []
    offset = dimension * numDiscrete
    for i1 in range(offset, offset + dimension * numContinua, dimension):
        continua.append(floats[i1:i1+dimension])

    return discretes, continua


def nStringsToInts(n_, startLine, lines, dimension=1):

    ints = []
    nd = dimension * n_
    while True:
        line = lines[startLine]
        n6 = min(nd, 6)
        for i1 in range(n6): ints.append(int(line[11*i1: 11*(i1+1)]))
        nd -= n6
        if nd == 0: break
        startLine += 1
    if dimension > 1:
        ints1D, ints = ints, []
        for i1 in range(0, dimension * n_, dimension): ints.append(ints1D[i1:i1+dimension])
    return ints


def getENDFDate(date):

    try:
        y, m, d = date[4:6], date[6:8], date[8:]
        if int(y) < 50:
            y = '20' + y
        else:
            y = '19' + y
        date = '%s-%s-%s' % (y, m, d)
    except Exception:
        sMonth = date[5:8]
        try:
            month = {"JAN": '01', "FEB": '02', "MAR": '03', "APR": '04', "MAY": '05', "JUN": '06',
                     "JUL": '07', "AUG": '08', "SEP": '09', "OCT": '10', "NOV": '11', "DEC": '12'
                     }[sMonth.upper()]
        except KeyError:
            raise Exception('Bad evaluation date: could not convert month = "%s" to integer '
                            'for "%s". See option --ignoreBadDate' % (sMonth, date))
        try:
            year = 1900 + int(date[8:10])
            if year < 1950: year += 100
        except ValueError:
            raise Exception('Bad evaluation date: could not convert year = "%s" to integer '
                            'for "%s". See option --ignoreBadDate' % (date[8:10], date))
        day = '01'
        date = '%s-%s-%s' % (year, month, day)
    return date


def getTAB1(startLine, dataLines, logFile=sys.stderr):

    lineNumber = startLine
    C1, C2, L1, L2, NR, NP = sixFunkyFloatStringsToFloats(dataLines[startLine], logFile=logFile)
    NR = int(NR)                  # number of interpolation flags
    NP = int(NP)                  # number of pairs of data
    lineNumber += 1
    interpolationInfo = nStringsToInts(NR, lineNumber, dataLines, dimension=2)
    if interpolationInfo[0][0] == NP:
        interpolationInfo = interpolationInfo[:1]  # There is really only one region, so fix it.

    lineNumber += (NR + 2) // 3
    data = nFunkyFloatStringsToFloats(NP, lineNumber, dataLines, dimension=2)

    lineNumber += (NP + 2) // 3
    return lineNumber, {'C1': C1, 'C2': C2, 'L1': L1, 'L2': L2, 'NR': NR, 'NP': NP,
                        'interpolationInfo': interpolationInfo, 'data': data}


def getList(startLine, dataLines, logFile=sys.stderr):

    lineNumber = startLine
    C1, C2, L1, L2, NPL, N2 = sixFunkyFloatStringsToFloats(dataLines[startLine], logFile=logFile)
    NPL = int(NPL)                # number of items in list
    N2 = int(N2)                  #

    lineNumber += 1
    data = nFunkyFloatStringsToFloats(NPL, lineNumber, dataLines, dimension=1)

    lineNumber += (NPL + 5) // 6
    return lineNumber, {'C1': C1, 'C2': C2, 'L1': L1, 'L2': L2, 'NPL': NPL, 'N2': N2, 'data': data}


def getTAB2Header(startLine, dataLines, logFile=sys.stderr):

    lineNumber = startLine
    C1, C2, L1, L2, NR, NZ = sixFunkyFloatStringsToFloats(dataLines[startLine], logFile=logFile)
    L1, L2 = int(L1), int(L2)
    NR = int(NR)                  # number of interpolation flags
    NZ = int(NZ)                  # number of TAB1's

    lineNumber += 1
    interpolationInfo = nStringsToInts(NR, lineNumber, dataLines, dimension=2)

    lineNumber += (NR + 2) // 3
    return lineNumber, {'C1': C1, 'C2': C2, 'L1': L1, 'L2': L2, 'NR': NR, 'NZ': NZ,
                        'interpolationInfo': interpolationInfo}


def getTAB2_TAB1s(startLine, dataLines, logFile=sys.stderr, axes=None):
    """
    This function currently does not support multiple interpolation regions or discontinuous functions vs. C2.
    If this is changed, some calling functions will need to be changed.
    """

    lineNumber, TAB2 = getTAB2Header(startLine, dataLines, logFile=logFile)
    TAB2s = []
    iStart, C2Prior = 0, -1e99          # C2 may be mu that starts at -1.
    piecewise = None
    for i1, (NZSub, interpolation) in enumerate(TAB2['interpolationInfo']):
        TAB1s = []
        if piecewise is not None: TAB1s.append(piecewise)
        for i2 in range(iStart, NZSub):
            lineNumber, TAB1, piecewise = getTAB1Regions(lineNumber, dataLines, logFile=logFile, dimension=2, axes=axes)
            if C2Prior == TAB1['C2']: raise NotImplementedError('hell - see document')
            TAB1s.append(piecewise)
        iStart = NZSub
        TAB2s.append([interpolation, TAB1s])
    TAB2['TAB2s'] = TAB2s
    return lineNumber, TAB2


def getTAB2_Lists(startLine, dataLines, logFile=sys.stderr):

    lineNumber, TAB2 = getTAB2Header(startLine, dataLines, logFile=logFile)
    Lists = []
    for i1 in range(TAB2['NZ']):
        lineNumber, List = getList(lineNumber, dataLines, logFile=logFile)
        Lists.append(List)
    TAB2['Lists'] = Lists
    return lineNumber, TAB2


def getTAB1Regions(startLine, dataLines, allowInterpolation6=False, logFile=sys.stderr,
                   dimension=1, axes=None, cls=XYs1dModule.XYs1d):

    endLine, TAB1 = getTAB1(startLine, dataLines, logFile=logFile)
    n1, i1, mode, data = TAB1['NR'], 0, 0, TAB1['data']

    for regionsEnd, interpolation in TAB1['interpolationInfo']:        # Search for 'flat' region next to linear region.
        if interpolation not in [1, 2]: break
    if interpolation in [1, 2]:                                # Only consider if all regions are flat and/or linear.
        regionsStart, interpolationRegions2 = 0, []            # We are now going to combine as many adjacent flat and linear regions as possible.
        for index, interpolationInfo in enumerate(TAB1['interpolationInfo']):
            regionsEnd, interpolation = interpolationInfo               # This loop converts all flat regions to linear if possible
            interpolationRegions2.append([regionsEnd, interpolation])   # (i.e. x_i < x_{i+1} and all y values are the same).
            if interpolation == 1:
                convertToLinear = True
                x1, y1 = data[regionsStart]
                for x2, y2 in data[regionsStart+1: regionsEnd]:
                    if x1 == x2 or y1 != y2:                # x_i < x_{i+1} and all y values must be the same.
                        convertToLinear = False
                        break
                    x1 = x2
                if convertToLinear: interpolationRegions2[-1][1] = 2
            regionsStart = regionsEnd
        x1, i1 = interpolationRegions2[0]
        interpolationRegions = []
        for x2, i2 in interpolationRegions2[1:]:                       # Now combine adjacent linear regions if possible
            if i1 == 2 and i2 == 2:
                if data[x1-1][0] == data[x1][0]: interpolationRegions.append([x1, i1])
            else:
                interpolationRegions.append([x1, i1])
            x1, i1 = x2, i2
        interpolationRegions.append([x1, i1])
    else:
        interpolationRegions = TAB1['interpolationInfo']

    if len(interpolationRegions) != len(TAB1['interpolationInfo']):     # Print a warning message if number of regions were reduced.
# BRB, check this with real data.
        if not (len(TAB1['interpolationInfo']) == 1 and
                TAB1['interpolationInfo'][0][1] == 1 and
                interpolationRegions[0][1] == 2):
            logFile.write(' reduced regions from %s to %s' %
                          (len(TAB1['interpolationInfo']), len(interpolationRegions)))

    if data[-2][0] == data[-1][0] and data[-1][1] == 0:     # discontinuity at final outgoing energy used to drop spectrum to 0
        if interpolation == 1:
            data.pop(-2)                                    # flat interpolation: omit 2nd-to-last point
        else:
            data[-1][0] *= 1.000000001                      # lin-lin interpolation: blur edge of discontinuity

    n1, regions = 0, []
    for interpolationInfo in interpolationRegions:
        n2, interpolation = interpolationInfo
        if n2 == n1:
            logFile.write("encountered interpolation region with no points")
            continue
        if interpolation == 6:
            threshold = 0.
            if not allowInterpolation6:
                raise Exception("charged-particle interpolation not allowed here")
        np1 = n1
        if n1 > 0:
            if data[n1-1][0] != data[n1][0]: np1 -= 1

        regionData2, regionData, x1, y1 = data[np1:n2], [], None, None
        for x2, y2 in regionData2:                                     # Loop to remove duplicate points.
            if x1 is not None:
                if (x1 == x2) and (y1 == y2): continue
            regionData.append([x2, y2])
            x1, y1 = x2, y2

        while len(regionData) > 0:
            i3, n3 = 1, len(regionData)
            x1, y1 = regionData[0]
            x2, y2 = regionData[1]
            while (i3 < n3) and (x1 != x2):
                x1, y1 = x2, y2
                i3 += 1
                if i3 < n3: x2, y2 = regionData[i3]
            value = None
            if dimension > 1: value = TAB1['C2']
            regions.append(cls(data=regionData[:i3], outerDomainValue=value,
                               interpolation=ENDFInterpolationToGNDS1d(interpolation),
                               axes=axes))
            regionData = regionData[i3:]
            for i3 in range(1, len(regionData)):
                if regionData[i3][0] != x1: break
            i3 -= 1
            regionData = regionData[i3:]
        n1 = n2
    return endLine, TAB1, regions


niceSortOfMTs = reactionSuiteModule.niceSortOfMTs


def getMFDataInMFList(MFs, MFData):

    Data = {}
    for MF in MFData:
        if MF in MFs: Data[MF] = MFData[MF]
    return Data


def toEnergyFunctionalData(info, dataLine, MF5Data, LF, moniker, unit, xLabel='energy_in',
                           xUnit='eV'):

    import fudge.productData.distributions.energy as energyModule

    axes = axesModule.Axes(2, labelsUnits={0: (moniker, unit), 1: (xLabel, xUnit)})
    dataLine, TAB1, data = getTAB1Regions(dataLine, MF5Data, logFile=info.logs, axes=axes)
    EFclass = None
    for tmp in (energyModule.A, energyModule.B, energyModule.Theta, energyModule.G, energyModule.T_M):
        if tmp.moniker == moniker:
            EFclass = tmp

    if len(data) == 1:
        oneD = data[0]
    else:
        oneD = regionsModule.Regions1d(axes=data[0].axes)
        for datum in data: oneD.append(datum)

    return dataLine, TAB1, EFclass(oneD)


def completeDocumentation(info, documentation):

    documentation.body.body = 'See the endfCompatible section.'
    # Save the evaluation date
    try:
        date = datesModule.Date(dateModule.Date.parse(info.Date), datesModule.DateType.created)
    except Exception as e:
        if not info.ignoreBadDate:
            info.doRaise.append(str(e))
        import datetime
        dateString = datetime.datetime.today().strftime("%Y-%m-%d")
        date = datesModule.Date(dateModule.Date.parse(dateString), datesModule.DateType.created)
    documentation.dates.add(date)
    documentation.title.body = ('ENDF-6 file from %s version %s translated to GNDS by FUDGE.' %
                                (info.library, info.libraryVersion))
    if '&' in info.author:
        info.author = info.author.replace('&', 'and')  # make XML-friendly.  Could also split into two authors
    documentation.authors.add(authorModule.Author(info.author, '', ''))
