import fudge

from xData import XYs1d as XYs1dModule
from xData import axes as axesModule
from xData import regions as regionsModule
from fudge.vis.matplotlib import plot2d, DataSet2d, DataSet3d
import fudge.sums
from brownies.legacy.endl import misc as miscENDLModule
from . import plotstyles

observableNameMap = {
    3: "Cross Section ( Sigma(E) )",
    23: "Cross Section ( Sigma(E) )",
    4: "Angular Distribution ( dSigma(E)/d(cos Theta) )",
    5: "Energy Distribution ( dSigma(E)/d(E') )",
    6: "Energy-Angle Distribution ( dSigma(E)/d(E')d(cos Theta) )"}

ENDFStandardsReactions = (
    ('H1', 2), ('C0', 2), ('AU197', 102), ('U235', 18), ('U238', 18), ('U238', 102), ('PU239', 18))


def isStandards(projectile, target, mt):
    if projectile != 'n':
        return False
    if mt not in [2, 102, 18]:
        return False
    if target not in ['H1', 'C0', 'AU197', 'U235', 'U238', 'PU239']:
        return False
    for x in ENDFStandardsReactions:
        if (target, mt) == x:
            return True
    return False


def hasMT(reactionSuite, mt):
    return mt in getEvaluationMTs(reactionSuite)


def getPlottableReactions(reactionSuite, observable='crossSection'):
    """
    Retrieve references to all reactions with the observable specified.
    For now, the observable should be an attribute of the reaction.
    :param reactionSuite:
    :param observable:
    :return:
    """
    result = []
    for reactionList in reactionSuite.reactions, reactionSuite.sums.crossSectionSums, reactionSuite.sums.multiplicitySums, \
                        reactionSuite.fissionComponents, reactionSuite.productions:
        for r in reactionList:
            if hasattr(r, observable):
                result.append(r)
    return result


def getEvaluationMTs(reactionSuite, mtFilter=None):
    reactionList = getPlottableReactions(reactionSuite)
    if mtFilter is None:
        return [int(r.ENDF_MT) for r in reactionList]
    else:
        return [int(r.ENDF_MT) for r in reactionList if int(r.ENDF_MT) in mtFilter]


def getReactions(reactionSuite, MT):
    return [r for r in getPlottableReactions(reactionSuite) if r.ENDF_MT == MT]


# ---------------------------------------------------
# ENDF-specific plotting helpers
# ---------------------------------------------------

# Simplified MT - EXFOR reaction mapping
def getEXFORRxn(MT, projectile='N'):
    try:
        from . import exfor2endf
        query_results = exfor2endf.get_exfor_reactions(projectile, mt=MT)
    except ImportError:
        # oops, no sqlalchemy, default to old heuristic
        query_results = None

    if query_results:
        # Got it from x4toc4 database, take first one (we're ignoring the observable)
        return query_results[0].split(')').lstrip('(').split(',')[-1]

    else:
        # Use older heuristic
        if MT == 2:
            return 'EL'
        if MT in range(50, 92) or MT == 4:
            if projectile == 'N':
                return 'INEL'
            else:
                return 'N'
        if MT in range(600, 650) or MT == 103:
            if projectile == 'P':
                return 'INEL'
            else:
                return 'P'
        if MT in range(650, 700) or MT == 104:
            if projectile == 'D':
                return 'INEL'
            else:
                return 'D'
        if MT in range(700, 750) or MT == 105:
            if projectile == 'T':
                return 'INEL'
            else:
                return 'T'
        if MT in range(750, 800) or MT == 106:
            if projectile == 'HE3':
                return 'INEL'
            else:
                return 'HE3'
        if MT in range(800, 850) or MT == 107:
            if projectile == 'A':
                return 'INEL'
            else:
                return 'A'
        if MT in range(875, 891):
            return '2N'
        if MT in [19, 20, 21, 38]:
            return 'F'

        # Or, do everything for this projectile
        if projectile == 'N':
            return {
                1: 'TOT', 2: 'EL', 3: 'NON', 4: 'N', 5: 'X', 11: '2N+D', 16: '2N', 17: '3N', 18: 'F',
                22: 'N+A', 23: 'N+3A', 24: '2N+A', 25: '3N+A', 27: 'ABS', 28: 'N+P', 29: 'N+2A',
                30: '2N+2A', 32: 'N+D', 33: 'N+T', 34: 'N+HE3', 35: 'N+D+2A', 36: 'N+T+2A', 37: '4N',
                41: '2N+P', 42: '3N+P', 44: 'N+2P', 45: 'N+P+A', 102: 'G', 103: 'P', 104: 'D', 105: 'T',
                106: 'HE3', 107: 'A', 108: '2A', 109: '3A', 111: '2P', 112: 'P+A', 113: 'T+2A',
                114: 'D+2A', 115: 'P+D', 116: 'P+T', 117: 'D+A'}.get(MT, "MT=" + str(MT))
        elif projectile in ['P', 'D', 'T', 'A', 'HE3']:
            raise NotImplementedError('projectile %s' % projectile)
        elif projectile in ['G']:
            raise NotImplementedError('projectile %s' % projectile)
        elif projectile in ['E']:
            raise NotImplementedError('projectile %s' % projectile)


MTShellAssignments = {534: "K",
                      535: "L1", 536: "L2", 537: "L3",
                      538: "M1", 539: "M2", 540: "M3", 541: "M4", 542: "M5",
                      543: "N1", 544: "N2", 545: "N3", 546: "N4", 547: "N5", 548: "N6", 549: "N7",
                      550: "O1", 551: "O2", 552: "O3", 553: "O4", 554: "O5", 555: "O6", 556: "O7", 557: "O8", 558: "O9",
                      559: "P1", 560: "P2", 561: "P3", 562: "P4", 563: "P5", 564: "P6", 565: "P7", 566: "P8", 567: "P9",
                      568: "P10", 569: "P11",
                      570: "Q1", 571: "Q2", 572: "Q3"}


def getSuggestTitle(target, projectile, reaction, mt):
    if reaction is None:
        # Generic production data
        return target.capitalize() + '(' + projectile + ',X)'

    if reaction.lower() == 'inel' or mt >= 600:
        # Inelastic scattering
        if mt in range(50, 91):
            reactionString = "n[" + str(mt % 50) + "]"
        elif mt == 91:
            reactionString = "n[c]"
        elif mt in range(600, 649):
            reactionString = "p[" + str((mt - 550) % 50) + "]"
        elif mt == 649:
            reactionString = "p[c]"
        elif mt in range(650, 691):
            reactionString = "d[" + str((mt - 600) % 50) + "]"
        elif mt == 691:
            reactionString = "d[c]"
        elif mt in range(700, 749):
            reactionString = "t[" + str((mt - 650) % 50) + "]"
        elif mt == 749:
            reactionString = "t[c]"
        elif mt in range(750, 791):
            reactionString = "3He[" + str((mt - 700) % 50) + "]"
        elif mt == 791:
            reactionString = "3He[c]"
        elif mt in range(800, 849):
            reactionString = "a[" + str((mt - 750) % 50) + "]"
        elif mt == 849:
            reactionString = "a[c]"
        elif mt in range(850, 891):
            reactionString = "2n[" + str((mt - 800) % 50) + "]"
        elif mt == 891:
            reactionString = "2n[c]"
        else:
            reactionString = 'inel'

    elif mt in range(500, 573):
        # Atomic reaction data
        if mt == 500:
            reactionString = "Total charged particle stopping power"
        elif mt == 501:
            reactionString = "Total photon interaction"
        elif mt == 502:
            reactionString = "Photon coherent scattering"
        elif mt == 504:
            reactionString = "Photon incoherent scattering"
        elif mt in [505, 506]:
            reactionString = "Photon scattering factor"
        elif mt == 515:
            reactionString = "Pair production, electron field"
        elif mt == 516:
            reactionString = "Total pair production"
        elif mt == 517:
            reactionString = "Pair production, nuclear field"
        elif mt == 522:
            reactionString = "Photo-electric absorption"
        elif mt == 523:
            reactionString = "Photo-excitation"
        elif mt == 526:
            reactionString = "Electro-atomic scattering"
        elif mt == 527:
            reactionString = "Bremstrahlung"
        elif mt == 528:
            reactionString = "Electro-atomic excitation"
        elif mt == 533:
            reactionString = "Atomic relaxation"
        else:
            reactionString = MTShellAssignments.get(mt, 'MT=%i' % mt) + ' shell ionization'
        return projectile + '+' + target.capitalize() + ", " + reactionString + " (MT=%i)" % mt
    else:
        # The reaction string better be correct
        reactionString = reaction.lower()

    return target.capitalize() + '(' + projectile + ',' + reactionString + ')'


def generatePlot(observable, dataSets, xyData=None, xydyData=None, xdxydyData=None, plotStyle=None, suggestTitle='',
                 suggestXLog=None, suggestYLog=None, suggestFrame=None, figsize=(20, 10),
                 useBokeh=False, outFile=None):
    # Set up the XY data
    xyDataSets = []
    if xyData is not None:
        xyDataSets = getXYDataSets(xyData, plotStyle)

    # Set up the XYdY data
    xydyDataSets = []
    if xydyData is not None:
        xydyDataSets = getXYdYDataSets(xydyData, plotStyle)

    # Set up the XdXYdY data
    xdxydyDataSets = []
    if xdxydyData is not None:
        xdxydyDataSets = getXdXYdYDataSets(xdxydyData, plotStyle)

    # Set up this plot styles
    thisPlotStyle = plotstyles.getThisPlotStyle(plotStyle, observable)

    # Save the original setting of the logx & logy switches
    originalLogX = thisPlotStyle["xAxis"]["log"]
    originalLogY = thisPlotStyle["yAxis"]["log"]

    # Whether to do lin-lin or log-log based on suggestLog flag
    if thisPlotStyle["xAxis"]["log"] is None:
        if suggestXLog is not None:
            thisPlotStyle["xAxis"]["log"] = suggestXLog
    if thisPlotStyle["yAxis"]["log"] is None:
        if suggestYLog is not None:
            thisPlotStyle["yAxis"]["log"] = suggestYLog

    # Check the reference frame
    if thisPlotStyle['referenceFrame'] is not None:
        if suggestFrame != thisPlotStyle['referenceFrame']:
            raise ValueError("Suggested frame from data of %s does not match required plot style %s" % (
                suggestFrame, thisPlotStyle['referenceFrame']))

    # load the axis style information into the plotting widget class instances
    xAxisSettings = plot2d.AxisSettings(axisMin=thisPlotStyle["xAxis"]["min"],
                                        axisMax=thisPlotStyle["xAxis"]["max"],
                                        label=thisPlotStyle["xAxis"]["label"],
                                        isLog=thisPlotStyle["xAxis"]["log"],
                                        unit=thisPlotStyle["xAxis"]["unit"])
    yAxisSettings = plot2d.AxisSettings(axisMin=thisPlotStyle["yAxis"]["min"],
                                        axisMax=thisPlotStyle["yAxis"]["max"],
                                        label=thisPlotStyle["yAxis"]["label"],
                                        isLog=thisPlotStyle["yAxis"]["log"],
                                        unit=thisPlotStyle["yAxis"]["unit"])
    # Set the plot title, if needed
    if thisPlotStyle['title'] is None:
        thisPlotStyle['title'] = suggestTitle

    # Finally make the plot
    plot2d.makePlot2d(
        dataSets + xyDataSets + xydyDataSets + xdxydyDataSets,
        xAxisSettings=xAxisSettings,
        yAxisSettings=yAxisSettings,
        title=thisPlotStyle['title'],
        legendOn=True,
        outFile=outFile,
        legendXY=(thisPlotStyle['legendX'], thisPlotStyle['legendY']),
        figsize=figsize,
        useBokeh=useBokeh)

    # Put logx & logy back the way they were
    thisPlotStyle["xAxis"]["log"] = originalLogX
    thisPlotStyle["yAxis"]["log"] = originalLogY


# ---------------------------------------------------
# stuff for interacting with various data sources
# ---------------------------------------------------

def getEXFORSets(sym, A, metastable, reaction=None, quantity="SIG", nox4evals=True, nox4legend=False,
                 forceLegend=False, plotSyle=None, verbose=True):
    exforData = []
    try:
        from x4i import exfor_manager, exfor_entry
    except ImportError:
        print('WARNING: x4i not successfully imported (check your PYTHONPATH?), so EXFOR data not plotted')
        return exforData

    db = exfor_manager.X4DBManagerPlainFS()
    i = 0
    M = ""
    if metastable:
        metastableIndex = int(metastable.replace('m',''))
        if metastableIndex == 1:
            M = "-M"  # reformat for exfor retrieval
        else:
            print("***WARNING: EXFOR data retrieval not supported for metastable %s***" % metastable)
            return exforData

    subents = db.retrieve(target=sym + '-' + str(A) + M, reaction=reaction, quantity=quantity)
    if verbose:
        print(fudge.core.utilities.brb.banner(
            "Preparing EXFOR data for " + sym + '-' + str(A) + '(' + reaction.upper() + ')' + ', ' + quantity))
    if verbose:
        print('Retrieving entries:' + str(list(subents.keys())))
    suppressEXFORLegend = (nox4legend or len(list(subents.keys())) > 20) and not forceLegend
    if verbose:
        print('Search parameters:', sym + '-' + str(A), reaction, quantity)
    if verbose:
        print('Found following (sub)entries:')
    for e in subents:
        if nox4evals and e.startswith('V'):
            continue
        if verbose:
            print('    Entry:', e)
        try:
            # New version of x4i returns an X4Entry directly
            if isinstance(subents[e], exfor_entry.X4Entry):
                ds = subents[e].getSimplifiedDataSets(makeAllColumns=True)
            # Old versions of x4i return a string that we have to convert to an X4Entry
            else:
                ds = exfor_entry.X4Entry(subents[e]).getSimplifiedDataSets(makeAllColumns=True)
        except KeyError:
            continue
        for d in ds:
            # TO DO: skip if data ratio data or SPA or MXW ave'd
            legend = ds[d].legend()
            if verbose:
                print('       ', d, legend)
            if 'Energy' in ds[d].data and 'Data' in ds[d].data:
                dat = ds[d].data[['Energy', 'Data']].pint.dequantify().values.tolist()
                try:
                    unc = ds[d].data[['d(Energy)', 'd(Data)']].pint.dequantify().values.tolist()
                except KeyError as err:  # There was a problem getting uncertainties, so skip it
                    unc = None
                xUnit = str(ds[d].data['Energy'].pint.units)
                yUnit = str(ds[d].data['Data'].pint.units)
                if yUnit == 'b/sr':
                    yUnit='b'  # newish strange choice in EXFOR management
                if e.startswith('V') or not suppressEXFORLegend:
                    theLegend = legend + ' (' + str(d[0]) + '.' + str(d[1][-3:]) + ')'
                else:
                    theLegend = '_noLegend_'
                exforData.append(
                    DataSet2d(data=dat, uncertainty=unc, 
                              xUnit=xUnit, 
                              yUnit=yUnit,
                              legend=theLegend, lineStyle=' ', symbol=plotstyles.getPlotSymbol(i),
                              color=plotstyles.getPlotColor(theLegend, False)))
                i += 1
    if verbose:
        print()
    return exforData


def getC4DataSets(c4File, mt, mf, target, projectile, plotSyle=None):
    try:
        from brownies.BNL import c4
    except ImportError:
        print('WARNING: c4 not successfully imported (check your PYTHONPATH?), so X4TOC4 data not plotted')
        return []
    print(fudge.core.utilities.brb.winged_banner('Retrieving C4 data from: ' + str(c4File)))
    flist = [x for x in c4.readC4File(open(c4File).readlines(), asPointList=True) if x.MT == mt and x.MF == mf
             and x.target == miscENDLModule.getZAFromName(target)
             and x.projectile == miscENDLModule.getZAFromName(projectile)]
    i = -1
    dat = []
    unc = []
    c4Data = []
    lastSet = (None, None)
    thisSet = (None, None)
    for p in flist:
        thisSet = (p.reference, p.exforEntry + str(p.exforSubEntry).zfill(3))
        if lastSet != thisSet:
            if lastSet != (None, None):
                print('    Set:', lastSet)
                theLegend = lastSet[0] + ' (' + lastSet[1] + ')'
                if mf == 3:
                    c4Data.append(
                        DataSet2d(data=dat, uncertainty=unc, xUnit='eV', yUnit='b', legend=theLegend, lineStyle=' ',
                                  symbol=plotstyles.getPlotSymbol(i), color=plotstyles.getPlotColor(theLegend, False)))
                elif mf == 4:
                    c4Data.append(
                        DataSet3d(data=dat, uncertainty=unc, xUnit='eV', yUnit='', zUnit='b', legend=theLegend,
                                  lineStyle=' ', symbol=plotstyles.getPlotSymbol(i),
                                  color=plotstyles.getPlotColor(theLegend, False)))
                dat = []
                unc = []
            i += 1
            lastSet = thisSet
        dat.append([p.energy, p.data])
        if p.dEnergy is not None:
            unc.append([p.dEnergy])
        else:
            unc.append([0.0])
        if p.dData is not None:
            unc[-1].append(p.dData)
        else:
            unc[-1].append(0.0)
        if mf == 4:
            dat[-1].append(p.cosMuOrLegendreOrder)
            if p.dCosMuOrLegendreOrder is not None:
                unc[-1].append(p.dCosMuOrLegendreOrder)
            else:
                unc[-1].append(0.0)
    lastSet = thisSet
    if lastSet == (None, None):
        return []
    print('    Set:', lastSet)
    theLegend = lastSet[0] + ' (' + lastSet[1] + ')'
    if mf == 3:
        c4Data.append(DataSet2d(data=dat, uncertainty=unc, xUnit='eV', yUnit='b', legend=theLegend, lineStyle=' ',
                                symbol=plotstyles.getPlotSymbol(i), color=plotstyles.getPlotColor(theLegend, False)))
    elif mf == 4:
        c4Data.append(
            DataSet3d(data=dat, uncertainty=unc, xUnit='eV', yUnit='', zUnit='b', legend=theLegend, lineStyle=' ',
                      symbol=plotstyles.getPlotSymbol(i), color=plotstyles.getPlotColor(theLegend, False)))
    print()
    return c4Data


def getXYDataSets(xyData, plotStyle):
    xyDataSets = []
    for k in xyData:
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'xyCurve', k)
        xyDataSets.append(DataSet2d(data=xyData[k],
                                    legend=thisSetStyle['legend'],
                                    lineWidth=thisSetStyle['lineWidth'],
                                    lineStyle=thisSetStyle['lineStyle'],
                                    color=thisSetStyle['lineColor'],
                                    symbol=thisSetStyle['symbol'],
                                    dataType=thisSetStyle['dataType'],
                                    xUnit=thisSetStyle['xUnit'],
                                    yUnit=thisSetStyle['yUnit']))
    return xyDataSets


def getXYdYDataSets(xydyData, plotStyle):
    # Set up theXYdY data
    xydyDataSets = []
    for k in xydyData:
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'xydyCurve', k)
        d = [[x[0], x[1]] for x in xydyData[k]]
        u = [[0.0, x[2]] for x in xydyData[k]]
        xydyDataSets.append(DataSet2d(data=d,
                                      uncertainty=u,
                                      legend=thisSetStyle['legend'],
                                      lineWidth=thisSetStyle['lineWidth'],
                                      lineStyle=thisSetStyle['lineStyle'],
                                      color=thisSetStyle['lineColor'],
                                      symbol=thisSetStyle['symbol'],
                                      dataType=thisSetStyle['dataType'],
                                      xUnit=thisSetStyle['xUnit'],
                                      yUnit=thisSetStyle['yUnit'],
                                      errorbarColor=thisSetStyle['errorColor']))
    return xydyDataSets


def getXdXYdYDataSets(xdxydyData, plotStyle):
    # Set up the XdXYdY data
    xdxydyDataSets = []
    for k in xdxydyData:
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'xdxydyCurve', k)
        d = [[x[0], x[2]] for x in xdxydyData[k]]
        u = [[x[1], x[3]] for x in xdxydyData[k]]
        xdxydyDataSets.append(DataSet2d(data=d,
                                        uncertainty=u,
                                        legend=thisSetStyle['legend'],
                                        lineWidth=thisSetStyle['lineWidth'],
                                        lineStyle=thisSetStyle['lineStyle'],
                                        color=thisSetStyle['lineColor'],
                                        symbol=thisSetStyle['symbol'],
                                        dataType=thisSetStyle['dataType'],
                                        xUnit=thisSetStyle['xUnit'],
                                        yUnit=thisSetStyle['yUnit'],
                                        errorbarColor=thisSetStyle['errorColor']))
    return xdxydyDataSets


# ---------------------------------------------------
# cross section plots
# ---------------------------------------------------
def makeCrossSectionPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None, c4File=None, mt=0,
                         projectile='n', target='1H', nounc=False, nox4=False, showparts=False, nox4evals=True,
                         nox4legend=False, evaluationStyle='recon', logX=None, logY=None,
                         outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):

    # Deal with default values
    if gndsMap is None:
        gndsMap = {}
    if xyData is None:
        xyData = []
    if xydyData is None:
        xydyData = []
    if plotStyle is None:
        plotStyle = {}
    isoFiles = []
    mixtureSum = None
    eMin = None
    eMax = None
    reactionName = ''

    # Declare the main reaction type
    try:
        reactionName = getEXFORRxn(mt)
    except KeyError:
        if mt in range(850, 871, 1):
            print("Got lumped covariance for MT = " + str(mt))
        else:
            raise KeyError("Unknown MT: " + str(mt))

    # A place to store the sum of isotopes if we are dealing with isotopic mixtures
    if 'mixture' in gndsMap:
        import pqu.PQU as PQU
        eMin = PQU.PQU(str(gndsMap['mixture']['eMin']))
        eMax = PQU.PQU(str(gndsMap['mixture']['eMax']))
        isoFiles = [gndsMap['mixture']['isotopes'][iso]['pathToFile'] for iso in gndsMap['mixture']['isotopes']]

    # Get the ENDF data for plotting
    endfData = []
    rawEndfData = []
    rawEndfUnc = []
    for endf in gndsMap:
        if endf == 'mixture':
            continue
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'evaluation', endf)

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        weight = 1.0
        if 'mixture' in gndsMap:
            for iso in gndsMap['mixture']['isotopes']:
                if gndsMap['mixture']['isotopes'][iso]['pathToFile'] == endf:
                    weight = gndsMap['mixture']['isotopes'][iso]['atomicFraction']
                    thisSetStyle['legend'] = gndsMap['mixture']['isotopes'][iso].get('legend', thisSetStyle['legend'])
                    thisSetStyle['lineWidth'] = gndsMap['mixture']['isotopes'][iso].get('lineWidth',
                                                                                        thisSetStyle['lineWidth'])
                    thisSetStyle['lineStyle'] = gndsMap['mixture']['isotopes'][iso].get('lineStyle',
                                                                                        thisSetStyle['lineStyle'])
                    thisSetStyle['lineColor'] = gndsMap['mixture']['isotopes'][iso].get('lineColor',
                                                                                        thisSetStyle['lineColor'])
                    break

        try:
            try:
                thisReaction = reactionSuite.getReaction(mt)
            except KeyError:
                thisReaction = None

            if thisReaction is None:
                print("Reaction with MT=%i not present in evaluation %s" % (mt, endf))
            else:
                try:
                    # Extract the cross section data from the endf file
                    csData = thisReaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8) * weight
                    if 'mixture' in gndsMap and endf in isoFiles:
                        if mixtureSum is None:
                            mixtureSum = csData.domainSlice(eMin.getValueAs(csData.domainUnit),
                                                            eMax.getValueAs(csData.domainUnit))
                        else:
                            mixtureSum += csData.domainSlice(eMin.getValueAs(csData.domainUnit),
                                                             eMax.getValueAs(csData.domainUnit))
                    rawEndfData.append(csData)

                    # Extract the uncertainty on the cross section
                    if nounc:
                        rawEndfUnc.append(None)
                    else:
                        covariance = thisReaction.crossSection.getMatchingCovariance(covarianceSuite=covarianceSuite)
                        if covariance is None:
                            print("MT = %i (%s) covariance not present in the file" % (mt, reactionName))
                            rawEndfUnc.append(None)
                        else:
                            try:
                                rawEndfUnc.append(covariance.getUncertaintyVector(relative=False))
                            except ValueError as err:
                                print(
                                    "WARNING: When adding uncertainty to plot from %s for MT %i failed "
                                    "with error ValueError(\"%s\")" % (endf, mt, str(err)))
                                rawEndfUnc.append(None)
                            except IndexError as err:
                                print(
                                    "WARNING: When adding uncertainty to plot from %s for MT %i failed "
                                    "with error IndexError(\"%s\")" % (endf, mt, str(err)))
                                rawEndfUnc.append(None)
                            except TypeError as err:
                                print(
                                    "WARNING: When adding uncertainty to plot from %s for MT %i failed "
                                    "with error TypeError(\"%s\")" % (endf, mt, str(err)))
                                rawEndfUnc.append(None)

                    # Add it to the plot
                    endfData.append(
                        DataSet2d(
                            rawEndfData[-1],
                            uncertainty=rawEndfUnc[-1],
                            legend=thisSetStyle['legend'],
                            lineWidth=thisSetStyle['lineWidth'],
                            lineStyle=thisSetStyle['lineStyle'],
                            color=thisSetStyle['lineColor'],
                            errorbarColor=thisSetStyle['errorColor']))

                except TypeError as err:
                    print("WARNING: When adding data to plot from %s for MT %i failed with error TypeError(\"%s\")" %
                          (endf, mt, str(err)))

            # and plot all the rest, they don't get covariances...
            if showparts:
                if isinstance(thisReaction, fudge.sums.CrossSectionSum):
                    for summand in thisReaction.summands:
                        # FIXME: I shouldn't have to do this, but apparently I need to
                        summand.setAncestor(thisReaction.summands)
                        summand.updateXPath()
                        # Now add it to the plot
                        otherCS = summand.follow(reactionSuite)
                        otherReaction = otherCS.ancestor
                        csData = otherCS.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
                        if csData is not None:
                            rawEndfData.append(csData)
                        else:
                            continue
                        if len(gndsMap) == 1:
                            endfData.append(DataSet2d(rawEndfData[-1], legend='MT=%i' % otherReaction.ENDF_MT))
                        else:
                            endfData.append(
                                DataSet2d(rawEndfData[-1], legend=endf + ' (MT=%i)' % otherReaction.ENDF_MT))
            print()

        except TypeError as err:
            print("WARNING: Adding data from %s for MT %i failed with error \"%s\"" % (endf, mt, str(err)))

    # Add in the sum of any isotopic mixture
    if 'mixture' in gndsMap:
        endfData.append(
            DataSet2d(
                mixtureSum,
                legend=gndsMap['mixture']['legend'],
                lineWidth=gndsMap['mixture']['lineWidth'],
                lineStyle=gndsMap['mixture']['lineStyle'],
                color=gndsMap['mixture']['lineColor']))

    # Get the C4 data for plotting
    c4Data = []
    if c4File is not None:
        c4Data = getC4DataSets(c4File, mt, 3, target, projectile, plotSyle=plotStyle)

    # Get the EXFOR data for plotting
    exforData = []
    if (not nox4) and (reactionName is not None):
        sym, A, m = miscENDLModule.elementAFromName(target)
        exforData = getEXFORSets(sym, A, m, reaction=projectile + ',' + reactionName, quantity="SIG", nox4evals=nox4evals,
                                 nox4legend=nox4legend, plotSyle=plotStyle)

    # Actually make the plot
    if endfData + exforData + c4Data != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        if logX is not None:
            suggestXLog = logX
        else:
            suggestXLog = mt in [1, 2, 18, 102] + list(range(501, 574))
        if logY is not None:
            suggestYLog = logY
        else:
            suggestYLog = mt in [1, 2, 18, 102] + list(range(501, 574))

        generatePlot(observable='crossSection',
                     dataSets=endfData + exforData + c4Data,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, reactionName, mt),
                     suggestXLog=suggestXLog,
                     suggestYLog=suggestYLog,
                     figsize=figsize,
                     useBokeh=useBokeh,
                     outFile=outFile)


def makeCrossSectionUncertaintyPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                                    mt=0, projectile='n', target='1H',
                                    evaluationStyle='eval',
                                    outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("Must refactor")

    # Assemble the cross section data into something matplotlib can plot
    if not cov:
        raise TypeError("File has no covariance")
    ratioDomain = (max(rawEndfUnc[-1].domainMin, rawEndfData[-1].domainMin),
                   min(rawEndfUnc[-1].domainMax, rawEndfData[-1].domainMax))
    ratio = rawEndfUnc[-1].xSlice(*ratioDomain) / rawEndfData[-1].xSlice(*ratioDomain)
    endfData.append(
        DataSet2d(ratio,
                  uncertainty=None,
                  legend=getStyle(plotStyle, 'evaluations', endf, 'legend',
                                  default=endf + ': (' + projectile + ',' + reaction.lower() + ')'),
                  lineWidth=getStyle(plotStyle, 'evaluations', endf, 'lineWidth'),
                  lineStyle=getStyle(plotStyle, 'evaluations', endf, 'lineStyle'),
                  color=getStyle(plotStyle, 'evaluations', endf, 'lineColor'),
                  errorbarColor=getStyle(plotStyle, 'evaluations', endf, 'errorColor')))


def makeCrossSectionIntegralsPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                                  mt=0, projectile='n', target='1H',
                                  evaluationStyle='eval',
                                  outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError()


# ---------------------------------------------------
# special atomic data plots
# ---------------------------------------------------
def makeEnergyTransferPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                           mt=0, projectile='e-', target='1H', product='e-',
                           evaluationStyle='eval',
                           outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("Ask Bret how to get at the energy transfer for LAW=8 data in MF=26, MT=528 or 527")

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture':
            continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        # collect complete list of matching reactions
        reactionList = []
        for MT in mtList:
            try:
                reactionList += endf_utils.getReactions(reactionSuite, mt)
            except:
                pass
        print("Exit Channels:", ', '.join(map(str, reactionList)))

        # Extract the cross section data from the endf file
        dist = reactionList[0].outputChannel.getProductWithName(
            product).distributions.energyTransfer.getNativeData().toPointwise_withLinearXYs(accuracy=1e-3,
                                                                                            upperEps=1e-8)
        thisSetStyle = plot_utils.getThisSetStyle(plotStyle, 'evaluation', endf)
        endfData.append(
            DataSet2d(
                dist,
                legend=thisSetStyle['legend'],
                lineWidth=thisSetStyle['lineWidth'],
                lineStyle=thisSetStyle['lineStyle'],
                color=thisSetStyle['lineColor'],
                errorbarColor=thisSetStyle['errorColor']))

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='energyTransfer',
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=projectile + '+' + target + '->' + str(reactionList[0]),
                     suggestXLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     suggestYLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     figsize=figsize,
                     useBokeh=useBokeh)


def makeFormFactorPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                       mt=0, projectile='gamma', target='1H', product='gamma',
                       evaluationStyle='eval',
                       outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("FIXME: I need updating for new FUDGE")

    # collect complete list of matching reactions
    if mt not in [502, 504]:
        mtList = [502, 504]
    else:
        mtList = [mt]
    mtPrefix = {502: "coherent", 504: "incoherent"}
    mtLineStyle = {502: "-", 504: ":"}

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture': continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        # Extract the cross section data from the endf file
        for MT in mtList:
            dist = getReactions(reactionSuite, MT)[0].outputChannel.getProductWithName(
                product).distributions.getNativeData()
            if hasattr(dist, 'formFactor'): dist = dist.formFactor
            dist = dist.toPointwise_withLinearXYs(accuracy=1e-3, upperEps=1e-8)
            thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'evaluation', endf)
            endfData.append(
                DataSet2d(
                    dist,
                    legend=thisSetStyle['legend'] + " (" + mtPrefix[MT] + ")",
                    lineWidth=thisSetStyle['lineWidth'],
                    lineStyle=mtLineStyle[MT],
                    # thisSetStyle['lineStyle'],
                    color=thisSetStyle['lineColor'],
                    errorbarColor=thisSetStyle['errorColor']))

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='formFactor',
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=projectile + '+' + target.capitalize() +
                                  ', Photon (in)coherent scattering form factors (MT=502,504)',
                     suggestXLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     suggestYLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     figsize=figsize,
                     useBokeh=useBokeh)


def makeAnomolousScatteringFactorPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                                      mt=502, projectile='gamma', target='1H', product='gamma',
                                      evaluationStyle='eval',
                                      outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("FIXME: I need updating for new FUDGE")

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture': continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        # Extract the cross section data from the endf file
        dist = getReactions(reactionSuite, mt)[0].outputChannel.getProductWithName(
            product).distributions.getNativeData()
        for x in ['anomalousScatteringFactor_imaginaryPart', 'anomalousScatteringFactor_realPart']:
            if not hasattr(dist, x): continue
            thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'evaluation', endf)
            endfData.append(
                DataSet2d(
                    getattr(dist, x),
                    legend=thisSetStyle['legend'] + ', ' + x.split('_')[-1],
                    lineWidth=thisSetStyle['lineWidth'],
                    lineStyle={'imaginaryPart': ':', 'realPart': '-'}[x.split('_')[-1]],
                    color=thisSetStyle['lineColor'],
                    errorbarColor=thisSetStyle['errorColor']))

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='anomolousScatteringFactor',
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, reaction, mt) + ', anomolous scattering factor',
                     suggestXLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     suggestYLog=len(set(mtList).intersection(list(range(501, 574)))) != 0,
                     figsize=figsize,
                     useBokeh=useBokeh)


# ---------------------------------------------------
# multiplicity plots
# ---------------------------------------------------
def makeMultiplicityPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                         mt=0, projectile='n', target='H1', product='n',
                         c4File=None,
                         showparts=True,
                         nox4evals=True,
                         nox4legend=False,
                         evaluationStyle='eval',
                         logX=False,
                         logY=False,
                         figsize=(20, 10),
                         outFile=None, plotStyle=None, useBokeh=False):
    import collections
    from fudge.outputChannelData.fissionFragmentData import delayedNeutron

    # Defaults
    if gndsMap is None:
        gndsMap = {}
    if xyData is None:
        xyData = []
    if xydyData is None:
        xydyData = []
    if xdxydyData is None:
        xdxydyData = []
    if plotStyle is None:
        plotStyle = {}
    showDelayed = False
    showPrompt = True
    nox4 = False
    reactionName = ''

    # Declare the main reaction type
    try:
        reactionName = getEXFORRxn(mt)
    except KeyError:
        if mt in range(850, 871, 1):
            print("Got lumped covariance for MT = " + str(mt))
        else:
            raise KeyError("Unknown MT: " + str(mt))

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture':
            continue
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'evaluation', endf)
        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        # Extract the data from the endf file & give the data rational names (esp. for nubar)
        thisReaction = getReactions(reactionSuite, mt)[0]
        thisChannel = thisReaction.outputChannel
        products = collections.OrderedDict()
        delayedTimeGroupIndex = 0
        for p in thisChannel.getProductsWithName(product):
            if thisReaction.isFission():
                if isinstance(p, delayedNeutron.Product):
                    rate = p.findAttributeInAncestry("rate")[evaluationStyle].float('1/s')
                    label = 'delayed[timegroup #%i (%s)]' % (delayedTimeGroupIndex, rate)
                    delayedTimeGroupIndex += 1
                else:
                    label = 'prompt'
            else:
                label = p.label
            products[label] = p.multiplicity.evaluated.toPointwise_withLinearXYs()

        # Compute total multiplicity
        totalMult = sum([products[p] for p in products])

        # The ENDF data to plot
        totalLegend = endf  # thisSetStyle['legend']
        if thisReaction.isFission():
            totalLegend += ' (total nubar)'
        endfData.append(
            DataSet2d(
                totalMult,
                uncertainty=None,
                legend=totalLegend,
                lineWidth=thisSetStyle['lineWidth'],
                lineStyle=thisSetStyle['lineStyle'],
                color=thisSetStyle['lineColor'],
                errorbarColor=thisSetStyle['errorColor']))

        # Partial multiplicities (if requested)
        if showparts:
            for p in products:
                thisLegend = endf  # thisSetStyle['legend']
                if thisChannel.isFission():
                    if p == 'prompt' and not showPrompt:
                        continue
                    if p != 'prompt' and not showDelayed:
                        continue
                thisLegend += ' (%s)' % p
                endfData.append(
                    DataSet2d(
                        products[p],
                        uncertainty=None,
                        legend=thisLegend,
                        lineWidth=thisSetStyle['lineWidth'] / 2,
                        lineStyle=thisSetStyle['lineStyle'],
                        color=thisSetStyle['lineColor'],
                        errorbarColor=thisSetStyle['errorColor']))

    # Get the C4 data for plotting
    c4Data = []
    if c4File is not None:
        c4Data = getC4DataSets(c4File, mt, 3, target, projectile, plotSyle=plotStyle)

    # Get the EXFOR data for plotting
    exforData = []
    if (not nox4) and (reactionName is not None):
        sym, A, m = miscENDLModule.elementAFromName(target)
        exforData = getEXFORSets(sym, A, m, reaction=projectile + ',' + reactionName, quantity="NU",
                                 nox4evals=nox4evals, nox4legend=nox4legend, plotSyle=plotStyle)

    # Actually make the plot
    if endfData + exforData + c4Data != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='nubar',
                     dataSets=endfData + exforData + c4Data,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=("multiplicity of %s from " % product) + getSuggestTitle(target, projectile,
                                                                                           reactionName, mt),
                     suggestXLog=logX,
                     suggestYLog=logY,
                     figsize=figsize,
                     useBokeh=useBokeh,
                     outFile=outFile)


# ---------------------------------------------------
# energy/momentum deposition plots
# ---------------------------------------------------
def makeEnergyDepositPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                          mt=0, projectile='n', target='1H', product='n',
                          evaluationStyle='eval', observable=None,
                          outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("FIXME: I need updating for new FUDGE")

    from fudge.productData import averageProductEnergy as averageProductEnergyModule

    endfData = []

    # Get the ENDF data for plotting
    for endf in gndsMap:
        if endf == 'mixture': continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        # Set up MT list.
        # This is the list of all MT that we will search for data to sum into the plot.
        # Usually, it is one element long, but there are some special cases:
        #     summed channels that were not specified in the ENDF file and lumped channels....
        if not hasMT(reactionSuite, mt): continue

        # We may need resonance reconstruction
        if mt in [1, 2, 18, 102] and args.doResonanceReconstruction:
            reactionSuite.reconstructResonances()

        # We need to compute energy depositions
        reactionSuite.calculateDepositionData({'verbosity': 0})

        # Get the right reaction
        reaction = getReactions(reactionSuite, mt)[0]
        print("Exit Channel:", str(reaction))

        # Compute the energy deposits and available energy
        Q = reaction.outputChannel.Q.toPointwise_withLinearXYs(accuracy=1e-3, upperEps=1e-8)
        energyDep = [
            [prod.getLabel(), prod.data[averageProductEnergyModule.Component.genre].getNativeData()]
            for prod in reaction.outputChannel.particles if averageProductEnergyModule.Component.genre in prod.data]
        if energyDep:
            totalEDep = energyDep[0][1].copy()
            for idx in range(1, len(energyDep)):
                if totalEDep.getDomain() != energyDep[idx][1].getDomain():
                    totalEDep, energyDep[idx][1] = totalEDep.mutualify(1e-8, 0, 0, energyDep[idx][1], 1e-8, 0, 0)
                totalEDep += energyDep[idx][1]
        else:
            totalEDep = []
        expectedTotalEnergy = XYs1dModule.XYs1d(axes_=Q.axes, data=[[x[0], x[1] + x[0]] for x in Q.copyDataToXYs()],
                                        accuracy=1e-8)

        # Add sets to plot list
        thisSetStyle = plotstyles.getThisSetStyle(plotStyle, 'evaluation', endf, verbose=False)
        if observable == "energyDeposit":

            for i, p in enumerate(energyDep):
                endfData.append(DataSet2d(
                    p[1],
                    legend='Energy for product %s (%s)' % (str(p[0]), thisSetStyle['legend']),
                    lineWidth=max(thisSetStyle['lineWidth'] - 2, 2),
                    lineStyle=plotstyles.getLineStyle(i + 1),
                    color=thisSetStyle['lineColor']))
            endfData.append(DataSet2d(
                totalEDep,
                legend='Actual total energy (%s)' % thisSetStyle['legend'],
                lineWidth=thisSetStyle['lineWidth'] + 2,
                lineStyle=thisSetStyle['lineStyle'],
                color=thisSetStyle['lineColor']))
            endfData.append(DataSet2d(
                expectedTotalEnergy,
                legend='Expected total energy (%s)' % thisSetStyle['legend'],
                lineWidth=thisSetStyle['lineWidth'] + 5,
                lineStyle=thisSetStyle['lineStyle'],
                color='#cccccc'))

        else:
            raise NotImplementedError(observable)

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable=observable,
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, str(reaction), mt),
                     suggestXLog=len(
                         set(mtList).intersection([1, 2, 4, 18, 102, 103, 105, 106, 107] + range(501, 574))) != 0,
                     suggestYLog=len(
                         set(mtList).intersection([1, 2, 4, 18, 102, 103, 105, 106, 107] + range(501, 574))) != 0,
                     suggestFrame='lab',
                     figsize=figsize,
                     useBokeh=useBokeh)


def makeMomentumDepositPlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                            mt=0, projectile='n', target='1H', product='n',
                            evaluationStyle='eval', observable=None,
                            outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError("FIXME: I need updating for new FUDGE")

    from fudge.productData import averageProductMomentum as averageProductMomentumModule
    import math

    projectileMass = None
    endfData = []

    # Get the ENDF data for plotting
    for endf in gndsMap:
        if endf == 'mixture':
            continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]

        if projectileMass is None:
            projectileMass = gndsMap[endf][0].projectile.getMass('eV/c/c')

        # Set up MT list.
        # This is the list of all MT that we will search for data to sum into the plot.
        # Usually, it is one element long, but there are some special cases:
        #     summed channels that were not specified in the ENDF file and lumped channels....
        if not endf_utils.hasMT(reactionSuite, mt):
            continue

        # We may need resonance reconstruction
        if mt in [1, 2, 18, 102] and args.doResonanceReconstruction:
            reactionSuite.reconstructResonances()

        # We need to compute energy depositions
        reactionSuite.calculateDepositionData({'verbosity': 0})

        # Get the right reaction
        reaction = endf_utils.getReactions(reactionSuite, mt)[0]
        print("Exit Channel:", str(reaction))

        # Compute the energy deposits and available energy
        #        Q = reaction.outputChannel.Q.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
        momDep = [
            [prod.getLabel(), prod.data[averageProductMomentumModule.Component.genre].getNativeData()]
            for prod in reaction.outputChannel.particles if averageProductMomentumModule.Component.genre in prod.data]
        if momDep:
            totalMomDep = momDep[0][1].copy()
            for idx in range(1, len(momDep)):
                if totalMomDep.getDomain() != momDep[idx][1].getDomain():
                    totalMomDep, momDep[idx][1] = totalMomDep.mutualify(1e-8, 0, 0, momDep[idx][1], 1e-8, 0, 0)
                totalMomDep += momDep[idx][1]
        else:
            totalMomDep = []
        expectedTotalMom = XYs1dModule.XYs1d(axes_=totalMomDep.axes,
                                     data=[[x[0], math.sqrt(2.0 * x[0] * projectileMass)] for x in
                                           totalMomDep.copyDataToXYs()], accuracy=1e-8)

        # Add sets to plot list
        thisSetStyle = plot_utils.getThisSetStyle(plotStyle, 'evaluation', endf)
        if observable == "momentumDeposit":

            for i, p in enumerate(momDep):
                endfData.append(DataSet2d(
                    p[1],
                    legend='Forward momentum of product %s (%s)' % (str(p[0]), thisSetStyle['legend']),
                    lineWidth=max(thisSetStyle['lineWidth'] - 2, 2),
                    lineStyle=plot_utils.getLineStyle(i + 1),
                    color=thisSetStyle['lineColor']))
            endfData.append(DataSet2d(
                totalMomDep,
                legend='Actual total forward momentum (%s)' % thisSetStyle['legend'],
                lineWidth=thisSetStyle['lineWidth'] + 2,
                lineStyle=thisSetStyle['lineStyle'],
                color=thisSetStyle['lineColor']))
            endfData.append(DataSet2d(
                expectedTotalMom,
                legend='Expected total forward momentum (%s)' % thisSetStyle['legend'],
                lineWidth=thisSetStyle['lineWidth'] + 5,
                lineStyle=thisSetStyle['lineStyle'],
                color='#cccccc'))

        else:
            raise NotImplementedError(observable)

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable=observable,
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, str(reaction), mt),
                     suggestXLog=len(
                         set(mtList).intersection([1, 2, 4, 18, 102, 103, 105, 106, 107] + list(range(501, 574)))) != 0,
                     suggestYLog=len(
                         set(mtList).intersection([1, 2, 4, 18, 102, 103, 105, 106, 107] + list(range(501, 574)))) != 0,
                     suggestFrame='lab',
                     figsize=figsize,
                     useBokeh=useBokeh)


def makeFissionEnergyReleasePlot(gndsMap=None, xyData=None, xydyData=None, xdxydyData=None,
                                 mt=0, projectile='n', target='1H', product='n',
                                 evaluationStyle='eval',
                                 outFile=None, plotStyle=None, figsize=(20, 10), useBokeh=False):
    raise NotImplementedError()


# ---------------------------------------------------
# angular distribution plots
# ---------------------------------------------------
def makeAngDistMubarPlot(gndsMap={}, xyData={}, xydyData={}, xdxydyData={},
                         mt=0, projectile='n', target='1H', product='n',
                         evaluationStyle='eval', referenceFrame=None,
                         outFile=None, plotStyle={}, figsize=(20, 10), useBokeh=False):
    # Declare the main reaction type
    try:
        reactionName = getEXFORRxn(mt)
    except KeyError:
        if mt in range(850, 871, 1):
            print("Skipping MT = " + str(mt) + ", it is lumped covariance")
            return
        else:
            raise KeyError("Unknown MT: " + str(mt))

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture': continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]
        try:
            thisReaction = reactionSuite.getReaction(mt)
        except KeyError:
            thisReaction = None  # no data to plot 'cuz MT not in the reactionSuite
        if thisReaction is not None:
            print("Retrieving data for" + str(thisReaction))
            angDistEvaluated = thisReaction.outputChannel.getProductWithName(product).distribution.evaluated
            angDistPointwise = angDistEvaluated.toPointwise_withLinearXYs(upperEps=1e-8)
            suggestFrame = angDistEvaluated.productFrame
            csData = thisReaction.crossSection.toPointwise_withLinearXYs()
            legend = product + ' given as ' + str(angDistEvaluated) + " in " + str(
                angDistEvaluated.productFrame) + ' frame (' + endf + ')'
            endfData.append(
                DataSet2d(
                    [[E, angDistPointwise.averageMu(E)] for E in angDistPointwise.getEnergyArray()],
                    xUnit=angDistPointwise.axes[0].unit,
                    yUnit='',
                    lineWidth=3,
                    legend=legend))
            print()

    # figure out reference frame
    if referenceFrame is None:
        if suggestFrame is None:
            suggestFrame = 'lab'
    else:
        suggestFrame = referenceFrame

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='mubar',
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, reactionName, mt),
                     suggestXLog=(mt in [1, 2, 18, 102] + list(range(501, 574))),
                     suggestYLog=(mt in [1, 2, 18, 102] + list(range(501, 574))),
                     suggestFrame=suggestFrame,
                     figsize=figsize,
                     useBokeh=useBokeh)


def makeAngDistLegendreMomentPlot(gndsMap={}, xyData={}, xydyData={}, xdxydyData={},
                                  mt=0, L=1, projectile='n', target='1H', product='n',
                                  evaluationStyle='eval', referenceFrame='centerOfMass',
                                  outFile=None, plotStyle={}, figsize=(20, 10), useBokeh=False):
    def get_LegendreMoments(theAngDist):

        if isinstance(theAngDist, regionsModule.Regions2d):
            myRegions = theAngDist.regions
        else:
            myRegions = [theAngDist]
        print('        Num. regions:', len(myRegions))
        Lmax = len(myRegions[-1].functionals[-1])
        theLegMoments = {L: [] for L in range(Lmax + 1)}

        # Reorder the distributions so we have a table of energies for each L
        for ireg, reg in enumerate(myRegions):
            eList = reg.getEnergyArray()
            print('            Regions #%i with %i points, ' % (ireg, len(reg)), reg.domainMin, '-', reg.domainMax,
                  reg.domainUnit)
            for iE, funcs in enumerate(reg):  # loop over functional containers in this region

                # If in Legendre moments already, just add them to the plot
                if isinstance(funcs, fudge.productData.distributions.angular.Legendre):
                    legs = funcs
                    for L, coeff in enumerate(legs):
                        theLegMoments[L].append((eList[iE], coeff))

                # If is pointwise, we'll fit each energy's XYs1d with Legendre moments before plotting
                else:
                    import numpy.polynomial.legendre as legendreModule
                    xys = (funcs - 0.5).copyDataToXsAndYs()  # subtract L=0 term, that better give 1/2
                    LMax = min(max(1, len(xys[1]) - 3), 50)
                    minResidual = 100.e12
                    for l in range(LMax):
                        fit_result = legendreModule.Legendre.fit(xys[0], xys[1], l, full=True)
                        legs = fit_result[0].coef
                        residuals = fit_result[1][0]
                        if minResidual < residuals:
                            break
                        minResidual = min(minResidual, residuals)
                    legs[0] = 0.5  # set the L=0 term since it had better be 1/2 (before fixing normalization)
                    for L, coeff in enumerate(legs):
                        theLegMoments[L].append((eList[iE], coeff * 2. / (
                                2. * L + 1.)))  # numpy & ENDF normalization a little different, so fix it here

        # Put the points in order
        for L in theLegMoments: theLegMoments[L].sort(key=lambda x: x[0])

        # Remove duplicates
        for L in theLegMoments:
            newList = []
            for pair in theLegMoments[L]:
                if len(newList) == 0 or pair[0] > newList[-1][0]:
                    newList.append(pair)
            theLegMoments[L] = newList

        # Construct the axes for all the XY's we'll need below
        coeffAxes = axesModule.Axes(2)
        coeffAxes[0] = theAngDist.axes[0]
        coeffAxes[1] = theAngDist.axes[-1]

        # Convert them to XYs
        return {L: XYs1dModule.XYs1d(data=theLegMoments[L], axes=coeffAxes) for L in range(Lmax + 1)}

    # Declare the main reaction type
    try:
        reactionName = getEXFORRxn(mt)
    except KeyError:
        if mt in range(850, 871, 1):
            print("Skipping MT = " + str(mt) + ", it is lumped covariance")
            return
        else:
            raise KeyError("Unknown MT: " + str(mt))

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndsMap:
        if endf == 'mixture':
            continue

        print(fudge.core.utilities.brb.winged_banner("Preparing data for " + endf))
        reactionSuite, covarianceSuite = gndsMap[endf]
        try:
            thisReaction = reactionSuite.getReaction(mt)
        except KeyError:
            thisReaction = None  # no data to plot 'cuz MT not in the reactionSuite
        if thisReaction is not None:
            print("Retrieving data for", thisReaction)
            angDistData = thisReaction.outputChannel.getProductWithName(product).distribution.evaluated.subforms[0]
            legend = 'L=%i Legendre moment for %s given in %s frame (%s)' % (
                L, product, angDistData.ancestor.productFrame, endf)
            print('   ', legend)
            LMomentData = get_LegendreMoments(angDistData)[L]
            endfData.append(
                DataSet2d(
                    LMomentData,
                    xUnit=LMomentData.axes[1].unit,
                    yUnit=LMomentData.axes[0].unit,
                    lineWidth=3,
                    legend=legend))
        print()

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
        generatePlot(observable='LegendreMoment',
                     dataSets=endfData,
                     xyData=xyData,
                     xydyData=xydyData,
                     xdxydyData=xdxydyData,
                     plotStyle=plotStyle,
                     suggestTitle=getSuggestTitle(target, projectile, reactionName, mt),
                     suggestXLog=(mt in [1, 2, 18, 102] or (501 <= mt <= 574)),
                     suggestYLog=(mt in [1, 2, 18, 102] or (501 <= mt <= 574)),
                     suggestFrame=referenceFrame,
                     figsize=figsize,
                     useBokeh=useBokeh,
                     outFile=outFile)
