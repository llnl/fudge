# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.core.utilities.brb import uniquify
from pqu import PQU

INCLASSAXISSETTING = False

# --------------------------------------------------------
#
# global data
#
# --------------------------------------------------------
"""
A few words about MatPlotLib symbols ("markers" in MatPlotLib-speak)...

.    point
,    pixel
o    circle
v    triangle_down 
^    triangle_up
<    triangle_left 
>    triangle_right 
1    tri_down -- DAB 11/17/2014: doesn't show up 
2    tri_up -- DAB 11/17/2014: doesn't show up 
3    tri_left -- DAB 11/17/2014: doesn't show up 
4    tri_right -- DAB 11/17/2014: doesn't show up 
8    octagon
s    square
p    pentagon
*    star
h    hexagon1
H    hexagon2
+    plus -- DAB 11/17/2014: barely visable
x    x -- DAB 11/17/2014: barely visable
D    diamond
d    thin_diamond
|    vline -- DAB 11/17/2014: doesn't show up 
_    hline -- DAB 11/17/2014: doesn't show up 
"""

defaultSymbols = ['o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']
nSymbols = len(defaultSymbols)

defaultLineStyles = ['-', '--', '-.', ':']
nLineStyles = len(defaultLineStyles)

"""
A few words about MatPlotLib colors....

Commands which take color arguments can use several formats to specify the colors. For the basic built-in colors, you 
can use a single letter

b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
w: white
Gray shades can be given as a string encoding a float in the 0-1 range, e.g.:

color = '0.75'
For a greater range of colors, you have two options. You can specify the color using an html hex string, as in:

color = '#eeefff'
or you can pass an R , G , B tuple, where each of R , G , B are in the range [0,1].

Finally, legal html names for colors, like 'red', 'burlywood' and 'chartreuse' are supported."""

defaultColors = ['k', 'r', 'g', 'b', 'y', 'brown', 'gray', 'violet', 'cyan',
                 'magenta', 'orange', 'indigo', 'maroon', 'turquoise',
                 'darkgreen', '0.20', '0.40', '0.60', '0.80']

nColors = len(defaultColors)


def convertUnit(unitFrom, unitTo, xs, xErrs, legend):
    if unitTo is None:
        return None
    if unitFrom is None:
        return None

    unitFrom = unitFrom.replace('EV', 'eV')
    if 'PER-CENT' in unitFrom:
        unitFrom = 'PERCENT'
    elif 'no-dim' in unitFrom:
        unitFrom = ''
    elif 'DEG-K' in unitFrom:
        unitFrom = 'K'
    elif 'KeV' in unitFrom:
        unitFrom = 'keV'
    unit = unitFrom
    if (unitFrom is not None) and (unitTo is not None) and (unitFrom != unitTo):
        try:
            conversionFactor = PQU.valueOrPQ(1.0, unitFrom=unitFrom, unitTo=unitTo, asPQU=False)
        except TypeError as err:
            if legend is None:
                legend = "name suppressed"
            raise TypeError("In plot2d.convertUnit, " + err.args[0] + ', error from PQU: ' + unitFrom +
                            ' -> ' + unitTo + ' for dataset ' + str(legend))
        except NameError as err:
            if legend is None:
                legend = "name suppressed"
            raise NameError("In plot2d.convertUnit, " +
                            err.args[0] + ', error from PQU: ' + unitFrom + ' -> ' + unitTo + ' for dataset ' +
                            str(legend))
        if xErrs is not None:
            for i1, x1 in enumerate(xErrs):
                xErrs[i1] = conversionFactor * x1
        for i1, x1 in enumerate(xs):
            xs[i1] = conversionFactor * x1
        unit = unitTo
    return unit


# --------------------------------------------------------
#
# Classes for new, matplotlib based plotting
#
# --------------------------------------------------------
class TextSettings:
    """
    Base class for text objects in plots
    
    Member data (set these with __init__ function):
    
    * **text** the text itself
    * **font** the font (e.g. Helvetica) (optional)
    * **fontSize** font size in points (optional, defaults to 14)
    * **fontWeight** font weight (e.g. bold) (optional)
    
    """

    def __init__(self, text=None, font=None, fontSize=12, fontWeight='bold'):
        self.text = text
        self.font = font
        self.fontSize = fontSize
        self.fontWeight = fontWeight


class AxisSettings(TextSettings):
    """
    Member data (set these with __init__ function) ::
    
        - `axisMin` :  min value of the axis 
            This is optional and defaults to `None` causing autoscaling.  
            If given as float, we assume in units of `unit`.  
            If given as a PQU instance, we'll 
            attempt to convert it to `unit`.
        - `axisMax` : max value of the axis.  See axisMin for usage
        - `isLog` : (optional, defaults to False) False == linear, True == semilog
        - `unit` : (optional, defaults to None), must be None or something that is usable in a 
             PQU instance
        - `gridOn` : (optional)
    
    """

    def __init__(self, axisMin=None, axisMax=None, autoscale=True, isLog=False, unit=None,
                 gridOn=False, label=None, font=None, fontSize=20, fontWeight='normal'):

        TextSettings.__init__(self, text=label, font=font, fontSize=fontSize, fontWeight=fontWeight)
        if isinstance(axisMin, PQU.PQU):
            if unit is not None:
                self.axisMin = axisMin.getValueAs(unit)
            else:
                self.axisMin = axisMin.getValue()
        else:
            self.axisMin = axisMin
        if isinstance(axisMax, PQU.PQU):
            if unit is not None:
                self.axisMax = axisMax.getValueAs(unit)
            else:
                self.axisMax = axisMax.getValue()
        else:
            self.axisMax = axisMax
        self.isLog = isLog
        self.gridOn = gridOn
        self.autoscale = autoscale
        if axisMin is not None and axisMax is not None:
            self.autoscale = False
        self.unit = unit
        self.label = self.text
        if self.unit is not None:
            if self.unit == '':
                self.label += ' (dimensionless)'
            else:
                self.label += ' (' + self.unit + ')'

    def __str__(self):

        units = ''
        if self.unit is not None:
            units = " %s" % self.unit
        mode = 'linear'
        if self.isLog:
            mode = 'log'
        if self.autoscale:
            s = 'Axis "%s" is autoscaled' % self.label
        else:
            s = 'Axis "%s" spans %s to %s%s' % (self.label, self.axisMin, self.axisMax, units)
        s += ' as %s' % mode
        if self.gridOn:
            s += ' with grids'
        font = 'default'
        if self.font is None:  # Need to fix this
            pass
        s += '. Font is "%s" with size %s and weight "%s".' % (font, self.fontSize, self.fontWeight)
        return s

    def setupAxis(self, theAxis):
        kwargs = {}
        if self.gridOn:
            kwargs['linewidth'] = 1.5
        theAxis.grid(self.gridOn, **kwargs)
        if INCLASSAXISSETTING:  # this doesn't' work, should I remove it?
            if self.isLog:
                theAxis.set_scale('log', linewidth=2)
            else:
                theAxis.set_scale('linear', linewidth=2)
            if not self.autoscale:
                theAxis.set_data_interval(self.axisMin, self.axisMax)
                theAxis.set_view_interval(self.axisMin, self.axisMax)


class DataSetParameters:
    """
    Base class for all DataSet classes (DataSet2d, DataSet3d, ...)
    
    Member data (set these with __init__ function):
    
    * legend (optional, defaults to '')
    * symbol
    * symbolColor
    * symbolSize
    * lineStyle
    * lineColor
    * lineWidth

    Commands which take color arguments can use several formats to specify the colors. 
    For the basic builtin colors, you can use a single letter
    
      =========   ====================
      character   description
      =========   ====================
              b   blue
              g   green
              r   red
              c   cyan
              m   magenta
              y   yellow
              k   black
              w   white
      =========   ====================
    
    Gray shades can be given as a string encoding a float in the 0-1 range, e.g.:
    
    >>> color = '0.75'
    
    For a greater range of colors, you have two options. You can specify the color 
    using an html hex string, as in:
    
    >>> color = '#eeefff'
    
    or you can pass an R , G , B tuple, where each of R , G , B are in the range [0,1].
    Finally, legal html names for colors, like "red", "burlywood" and "chartreuse" are 
    supported.

    The following format string characters are accepted to control the line style or 
    marker:
    
        =========   ====================
        character   description
        =========   ====================
        '-'         solid line style
        '--'        dashed line style
        '-.'        dash-dot line style
        ':'         dotted line style
        '.'         point marker
        ','         pixel marker
        'o'         circle marker
        'v'         triangle_down marker
        '^'         triangle_up marker
        '<'         triangle_left marker
        '>'         triangle_right marker
        '1'         tri_down marker
        '2'         tri_up marker
        '3'         tri_left marker
        '4'         tri_right marker
        's'         square marker
        'p'         pentagon marker
        '*'         star marker
        'h'         hexagon1 marker
        'H'         hexagon2 marker
        '+'         plus marker
        'x'         x marker
        'D'         diamond marker
        'd'         thin_diamond marker
        '|'         vline marker
        '_'         hline marker
        =========   ====================
    """

    def __init__(self, legend='', symbol=None, symbolSize=5, color=None, lineStyle=None,
                 lineWidth=None, errorbarCapSize=3, errorbarLineWidth=1.5, errorbarColor=None,
                 formatString=None):
        self.legend = legend
        self.symbol = symbol
        self.color = color
        self.symbolSize = symbolSize  # not yet used
        self.lineStyle = lineStyle
        self.lineWidth = lineWidth
        self.errorbarCapSize = errorbarCapSize
        self.errorbarLineWidth = errorbarLineWidth
        self.errorbarColor = errorbarColor
        self.formatString = formatString

    def getFormatString(self):
        """Attempt to make simple format string for matplotlib use."""
        if self.formatString is not None:
            return self.formatString
        r = ''
        if self.color is not None and len(self.color) == 1:
            r += self.color
        if self.symbol is not None and len(self.symbol) == 1:
            r += self.symbol
        if self.lineStyle is not None and len(self.lineStyle) == 1:
            r += self.lineStyle
        return r

    def getFormatMap(self, doErrors=False):
        def updateMap(theMap, key, theValue, defaultValue):
            if theValue is not None:
                theMap[key] = theValue
            elif defaultValue is not None:
                theMap[key] = defaultValue

        kw = {}
        updateMap(kw, 'color', self.color, None)
        updateMap(kw, 'marker', self.symbol, None)
        if self.symbol in ['+', 'x', '1', '2', '3', '4', '|', '_']:
            kw['markeredgewidth'] = 2
        else:
            kw['markeredgewidth'] = 1
        updateMap(kw, 'markersize', self.symbolSize, 5)
        updateMap(kw, 'linestyle', self.lineStyle, '-')
        updateMap(kw, 'linewidth', self.lineWidth, 1)
        updateMap(kw, 'label', self.legend, None)
        if doErrors:
            updateMap(kw, 'elinewidth', self.errorbarLineWidth, 1.5)
            updateMap(kw, 'capsize', self.errorbarCapSize, 3)
            updateMap(kw, 'ecolor', self.errorbarColor, self.color)
        return kw

    def addToPlot(self, thePlot):
        raise NotImplementedError("Implement in derived classes")


class DataSet2d(DataSetParameters):
    """
    Member data (set these with __init__ function):
    
    * data: the data as either a list of [ x, y ] pairs or as an object that has the method copyDataToXsAndYs
      (and also toPointwise_withLinearXYs if infill is True).
    * uncertainty: the uncertainty on the data as either a list of [ x, y ] pairs or as an object that has
      the method evaluate.
    * dataType ('scatter', 'line', or None): by default lists of points are rendered as 'scatter' and the
      others as 'line'.  Override the default behavior with this keyword.
    * xUnits: If the data object has an axes member, units are taken from it. Otherwise this argument is used
      to for the x units. When the units are given as None and the data/uncertainty do not have an axes member
      units are assumed to match those of the plot axis. The default is None.
    * yUnits: Units for the y-axis. See xUnits for more information.
    
    type help(DataSetParameters) to find out about how to set the line properties, etc.  
    All of this information is passed on though the __init__ function of DataSet2d to 
    DataSetParameters.
    """

    def __init__(self, data, uncertainty=None, xUnit=None, yUnit=None, dataType=None, infill=False, **kw):

        self.xerr = None
        self.yerr = None
        self.xUnit = xUnit
        self.yUnit = yUnit

        if isinstance(data, list) or isinstance(data, tuple):
            self.x = [p[0] for p in data]
            self.y = [p[1] for p in data]
            if dataType is None:
                self.dataType = 'scatter'
            else:
                self.dataType = dataType
            if 'lineStyle' not in kw:
                kw['lineStyle'] = '-'
            if isinstance(uncertainty, list) or isinstance(uncertainty, tuple):
                try:
                    self.xerr = [p[0] for p in uncertainty]
                except IndexError:
                    pass
                try:
                    self.yerr = [p[1] for p in uncertainty]
                except IndexError:
                    pass
        else:
            if dataType is None:
                self.dataType = 'line'
            else:
                self.dataType = dataType
            if hasattr(data, 'axes'):
                self.xUnit = data.axes[1].unit
                self.yUnit = data.axes[0].unit
            if infill:
                if not hasattr(data, 'toPointwise_withLinearXYs'):
                    raise TypeError("Unsupported type for infilling: %s" % type(data))
                data = data.toPointwise_withLinearXYs(accuracy=1e-3, upperEps=1e-6)
            if hasattr(data, 'copyDataToXsAndYs'):
                self.x, self.y = data.copyDataToXsAndYs()
                if uncertainty is not None:
                    self.yerr = [uncertainty.evaluate(xx) for xx in self.x]
            else:
                raise TypeError("Unsupported type: " + str(type(data)))

        if 'legend' not in kw:
            if hasattr(data, 'legend'):
                kw_ = {}
                for key in kw:
                    kw_[key] = kw[key]
                kw_['legend'] = data.legend
                kw = kw_

        if dataType is not None:
            self.dataType = dataType
        DataSetParameters.__init__(self, **kw)

    def convertUnits(self, xUnit=None, yUnit=None):
        """Convert self's units to xUnit and yUnit"""

        self.xUnit = convertUnit(self.xUnit, xUnit, self.x, self.xerr, self.legend)
        self.yUnit = convertUnit(self.yUnit, yUnit, self.y, self.yerr, self.legend)

    def addToPlot(self, thePlot, logY=False, verbose=False):
        if self.dataType == 'scatter':
            thePlot.errorbar(self.x, self.y, yerr=self.yerr, xerr=self.xerr,
                             **self.getFormatMap(doErrors=self.xerr is not None or self.yerr is not None))
        elif self.dataType == 'line':
            nRange = range(len(self.y))
            if not logY:
                thePlot.plot(self.x, self.y, **self.getFormatMap())
                if self.yerr is not None:
                    lowerbound = [self.y[_i] - abs(self.yerr[_i]) for _i in nRange]
                    upperbound = [self.y[_i] + abs(self.yerr[_i]) for _i in nRange]
                    thePlot.fill_between(self.x, lowerbound, upperbound, facecolor=self.color, alpha=0.25)
            else:
                theYs = []
                if self.yerr is not None:
                    theYLowBounds = []
                    theYUpBounds = []
                for _i in range(len(self.y)):
                    if verbose and self.y[_i] <= 0.0:
                        print("WARNING: y value <= 0.0 (", self.y[_i], ") at index", _i, "on a log plot")
                    theYs.append(self.y[_i])

                    if self.yerr is not None:
                        lb = self.y[_i] - abs(self.yerr[_i])
                        if verbose and lb <= 0.0:
                            print("WARNING: y-yerr value <= 0.0 (", lb, ") at index", _i, "on a log plot")
                        theYLowBounds.append(lb)

                        ub = self.y[_i] + abs(self.yerr[_i])
                        if verbose and ub <= 0.0:
                            print("WARNING: y+yerr value <= 0.0 (", ub, ") at index", _i, "on a log plot")
                        theYUpBounds.append(ub)
                thePlot.plot(self.x, theYs, **self.getFormatMap())
                if self.yerr is not None:
                    thePlot.fill_between(self.x, theYLowBounds, theYUpBounds, facecolor=self.color, alpha=0.25)
        else:
            raise TypeError('dataType must be "scatter" or "line", found ' + str(self.dataType))

    @property
    def dimension(self):

        return 2


class DataSet3d(DataSetParameters):
    """
    Member data (set these with __init__ function):
    
    * data
    * dataType (None)
    
    type help(DataSetParameters) to find out about how to set the line properties, etc.  
    All of this information is passed on though the __init__ function of DataSet2d to 
    DataSetParameters.
    """

    def __init__(self, data, dataType=None, xUnit=None, yUnit=None, zUnit=None, **kw):

        data_ = data
        if hasattr(data, 'toPointwise_withLinearXYs'):
            data_ = data.toPointwise_withLinearXYs(accuracy=1e-5, upperEps=1e-8)
        self.x, self.y, self.z = data_.copyDataToGridWsAndXsAndYs()

        self.xUnit, self.yUnit, self.zUnit = xUnit, yUnit, zUnit
        if hasattr(data, 'axes'):
            self.xUnit = data.axes[2].unit
            self.yUnit = data.axes[1].unit
            self.zUnit = data.axes[0].unit
        self.dataType = 'line'
        if dataType is not None:
            self.dataType = dataType
        DataSetParameters.__init__(self, **kw)

    def convertUnits(self, xUnit=None, yUnit=None, zUnit=None):
        """Convert self's units to xUnit and yUnit and zUnit"""

        self.xUnit = convertUnit(self.xUnit, xUnit, self.x, None, self.legend)
        self.yUnit = convertUnit(self.yUnit, yUnit, self.y, None, self.legend)
        self.zUnit = convertUnit(self.zUnit, zUnit, self.z, None, self.legend)

    def addToPlot(self, thePlot, plotType='contour', numContours=10, logX=False, logY=False):
        """Default plotType set to 'contour' so that can use the __makePlot2d function to drive both the regular 
        2d plotting and contour plotting"""

        if self.dataType == 'line':
            if plotType == 'contour':
                import matplotlib.pyplot as plt
                plt.clabel(thePlot.contour(self.x, self.y, self.z, numContours, **self.getFormatMap()), inline=1,
                           fontsize=10)
            else:
                raise TypeError('plotType ' + plotType + ' not supported')
        else:
            raise TypeError('dataType must be "scatter" or "line", found ' + str(self.dataType))

    def getW_XYs(self):

        from xData import axes as axesModule
        from xData import XYs1d as XYs1dModule
        from xData import multiD_XYs as multiD_XYsModule

        xUnit, yUnit, zUnit = '', '', ''
        if self.xUnit is not None:
            xUnit = self.xUnit
        if self.yUnit is not None:
            yUnit = self.yUnit
        if self.zUnit is not None:
            zUnit = self.zUnit
        axes_3d = axesModule.Axes(3)
        axes_3d[0] = axesModule.Axis('z', 0, zUnit)
        axes_3d[1] = axesModule.Axis('y', 0, yUnit)
        axes_3d[2] = axesModule.Axis('x', 0, xUnit)
        w_xys = multiD_XYsModule.XYs2d(axes=axes_3d)

        axes_2d = axesModule.Axes(2)
        axes_2d[0] = axes_3d[0]
        axes_2d[1] = axes_3d[1]

        for ix, x in enumerate(self.x):
            xys = [[y, self.z[iy][ix]] for iy, y in enumerate(self.y)]
            w_xys[ix] = XYs1dModule.XYs1d(xys, axes=axes_2d, outerDomainValue=x)
        return w_xys

    @property
    def dimension(self):

        return 3


def __makePlot2d(datasets, xAxisSettings=None, yAxisSettings=None, theTitle=None,
                 legendOn=False, legendXY=(0.05, 0.95), figsize=(20, 10), outFile=None,
                 thePlot=None, useBokeh=False, **kw):
    """
    Main driver routine for all 2d plots (regular, contour, slices, ...)
    """

    # Sometimes we'll adjust the layout of the plot external to this routine so the backend may already be set
    if outFile is not None:
        defaultBackend = 'Agg'
        backendMap = {'png': 'AGG', 'ps': 'PS', 'eps': 'PS', 'pdf': 'PDF', 'svg': 'SVG'}
        import matplotlib
        if matplotlib.get_backend() != defaultBackend:
            matplotlib.use(defaultBackend)
    import matplotlib.pyplot as plt

    # Check the axis settings.  If they are set to None and we can safely initialize them using information in the
    # appropriate XYs1d data structure, then do so
    if xAxisSettings is None:
        xUnit = None
        xLabel = None
        for ds in datasets:
            if hasattr(ds, 'axes'):
                if xUnit is None:
                    xUnit = ds.axes[1].unit
                if xUnit != ds.axes[1].unit:
                    raise ValueError("Incompatible units found in x axes of datasets, %s vs. %s" %
                                     (xUnit, ds.axes[1].unit))
                if xLabel is None:
                    xLabel = ds.axes[1].label
        xAxisSettings = AxisSettings(unit=xUnit, label=xLabel)
    if yAxisSettings is None:
        yUnit = None
        yLabel = None
        for ds in datasets:
            if hasattr(ds, 'axes'):
                if yUnit is None:
                    yUnit = ds.axes[0].unit
                if yUnit != ds.axes[0].unit:
                    raise ValueError("Incompatible units found in y axes of datasets, %s vs. %s" %
                                     (yUnit, ds.axes[0].unit))
                if yLabel is None:
                    yLabel = ds.axes[0].label
        yAxisSettings = AxisSettings(unit=yUnit, label=yLabel)

    # Create a plotting instance in a figure if one is needed.  
    # If the user passed in a plot, we'll return that.
    # Otherwise, all control of the plot will be handled in here.
    if thePlot is None:
        returnPlotInstance = False
        # Set up the plot canvas
        fig = plt.figure(figsize=figsize)
        # Add the plot title
        if theTitle is not None:
            fig.suptitle(theTitle.text, fontsize=theTitle.fontSize, fontweight=theTitle.fontWeight)
        # Set up the plot region
        thePlot = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
    else:
        returnPlotInstance = True

    # Add the data to the plot
    for dataset in datasets:
        try:
            dataset.convertUnits(xUnit=xAxisSettings.unit, yUnit=yAxisSettings.unit)
            dataset.addToPlot(thePlot, logY=yAxisSettings.isLog, **kw)
        except Exception as err:
            if 'error from PQU' in err.args[0]:
                pass  # harmless, data in incompatable units so we can't plot it
            else:
                raise err
    #            print("WARNING: could not add dataset "+str(dataset.legend)+", found error:")
    #            print(err.__class__, err)
    # Put on the plot legend
    if legendOn:
        thePlot.legend(bbox_to_anchor=legendXY, loc='upper left', borderaxespad=0., markerscale=1.0,
                       ncol=max(1, len(datasets) // 20), fancybox=True,
                       framealpha=0.5)  # , scatterpoints = 1, frameon = False ) # these last 3 options don't work?

    # Set up the axes
    xAxisSettings.setupAxis(thePlot.get_xaxis())
    yAxisSettings.setupAxis(thePlot.get_yaxis())

    if not INCLASSAXISSETTING:
        # Set the scales using the SubPlot calls (direct calls to XAxis and YAxis don't work because the scale
        # information is contained in the Axes class instance that contains both the XAxis and YAxis instances.
        # The Axis.set_limit functions actually call the stuff in the Axes instance and there is no equivalent
        # set_scale functions).
        if xAxisSettings.isLog:
            thePlot.set_xscale('log')
        else:
            thePlot.set_xscale('linear')
        if yAxisSettings.isLog:
            thePlot.set_yscale('log')
        else:
            thePlot.set_yscale('linear')
        if xAxisSettings.autoscale:
            thePlot.set_autoscalex_on(True)
        else:
            thePlot.set_xlim((xAxisSettings.axisMin, xAxisSettings.axisMax))
            thePlot.set_xbound((xAxisSettings.axisMin, xAxisSettings.axisMax))
        if yAxisSettings.autoscale:
            thePlot.set_autoscaley_on(True)
        else:
            thePlot.set_ylim((yAxisSettings.axisMin, yAxisSettings.axisMax))
            thePlot.set_ybound((yAxisSettings.axisMin, yAxisSettings.axisMax))
        thePlot.set_xlabel(xAxisSettings.label, size=xAxisSettings.fontSize, weight=xAxisSettings.fontWeight,
                           family=xAxisSettings.font)
        thePlot.set_ylabel(yAxisSettings.label, size=yAxisSettings.fontSize, weight=yAxisSettings.fontWeight,
                           family=yAxisSettings.font)
        for label in thePlot.get_xticklabels():
            label.set_size(xAxisSettings.fontSize)
            label.set_weight(xAxisSettings.fontWeight)
            label.set_family(xAxisSettings.font)
        for label in thePlot.get_yticklabels():
            label.set_size(yAxisSettings.fontSize)
            label.set_weight(yAxisSettings.fontWeight)
            label.set_family(yAxisSettings.font)
        # the following worked, but when you zoom in on the plot, the axes don't go too:
        #thePlot.set_xticklabels(thePlot.get_xticks(), size=xAxisSettings.fontSize, weight=xAxisSettings.fontWeight,
        #                        family=xAxisSettings.font )
        #thePlot.set_yticklabels(thePlot.get_yticks(), size=yAxisSettings.fontSize, weight=yAxisSettings.fontWeight,
        #                        family=yAxisSettings.font )
        thePlot.minorticks_on()

    if returnPlotInstance:
        return thePlot
    else:
        # Generate the output
        if outFile is not None:
            if useBokeh:
                #
                # Experimental bokeh usage
                #
                import bokeh.mpl
                import bokeh.plotting
                bokeh.plotting.output_file(outFile.replace('.png', '.html'))
                bokeh.plotting.show(bokeh.mpl.to_bokeh())
            else:
                plt.savefig(outFile)
                plt.clf()
        else:
            plt.show()


def makePlot2d(datasets, xAxisSettings=None, yAxisSettings=None, title='', legendOn=False, outFile=None,
               legendXY=(0.05, 0.95), figsize=(20, 10), infill=False, useBokeh=False):
    """
    Plain, vanilla, 2d plots.
    
    Arguments ::
    
        * datasets (required), list of DataSet objects to plot
        * xAxisSettings (optional) 
        * yAxisSettings (optional)
        * title (optional, defaults to '' )
        * legendOn (optional, defaults to True)
        * outFile, if not None, use this to output plot instead of generating a matplotlib window
        * legendXY, upper left corner of legend (optional, defaults to ( 0.05, 0.95 ))
        * figsize (optional, defaults to ( 20, 10 ))
        * infill
    """
    # Process the function arguments
    if isinstance(title, TextSettings):
        theTitle = title
    else:
        theTitle = TextSettings(title, fontSize=24)
    if not isinstance(xAxisSettings, AxisSettings):
        xAxisSettings = AxisSettings()
    if not isinstance(yAxisSettings, AxisSettings):
        yAxisSettings = AxisSettings()
    if not isinstance(datasets, (list, tuple)):
        datasets = [datasets]
    _datasets = []
    for dataset in datasets:
        if not isinstance(dataset, DataSet2d):
            dataset = DataSet2d(dataset, infill=infill)
        _datasets.append(dataset)
    __makePlot2d(_datasets, xAxisSettings, yAxisSettings, theTitle=theTitle, legendOn=legendOn, legendXY=legendXY,
                 figsize=figsize, outFile=outFile, useBokeh=useBokeh)


def makePlot2dContour(datasets, xAxisSettings=None, yAxisSettings=None, numContours=10, title='', figsize=(20, 10),
                      legendOn=False, outFile=None, useBokeh=False):
    """
    Contour plots:
    
    * datasets (required): list of DataSet objects to plot, one works better than many for a contour plot
    * xAxisSettings (optional) 
    * yAxisSettings (optional)
    * numContours (optional, defaults to 10): number of contours in the plot assuming linear spacing between intervals
                                              from the min to max data values
    * title (optional, defaults to '' )
    * legendOn (optional, defaults to True)
    * outFile (optional)
    
    """
    # Process the function arguments
    if isinstance(title, TextSettings):
        theTitle = title
    else:
        theTitle = TextSettings(title, fontSize=24)
    if not isinstance(xAxisSettings, AxisSettings):
        xAxisSettings = AxisSettings()
    if not isinstance(yAxisSettings, AxisSettings):
        yAxisSettings = AxisSettings()
    if not isinstance(datasets, (list, tuple)):
        datasets = [datasets]
    _datasets = []
    for dataset in datasets:
        if isinstance(dataset, DataSet3d):
            _datasets.append(dataset)
        else:
            _datasets.append(DataSet3d(dataset))
    if len(_datasets) != 1:
        print('Warning:  Contour plots with more than one dataset may not be very readable')
    __makePlot2d(_datasets, xAxisSettings, yAxisSettings, theTitle=theTitle, legendOn=legendOn, outFile=outFile,
                 figsize=figsize, numContours=numContours, useBokeh=useBokeh)


def makePlot2dSlice(datasets, xyAxisSettings=None, zAxisSettings=None, sliceXCoords=None, sliceYCoords=None,
                    sliceUnits='', title='', legendOn=True, legendXY=(0.05, 0.95), figsize=(20, 10), outFile=None,
                    useBokeh=False):
    """
        Make a slice of a dataset of 3d data objects along a fixed x line or a fixed y line.  
        
        Any 2d objects get plotted as if they are xy data in the correct plane for plotting.
        
        In other words, at fixed sliceXCoord = 2, all 3d functions F( x, y ) get plotted 
        as y vs. F( 2.0, y ) and all 2d functions G( x ) are plotted as x vs. G( x ) even 
        if you meant something different.
        
        * datasets (required): list of DataSet objects to plot
        * xyAxisSettings (optional) 
        * zAxisSettings (optional)
        * sliceXCoords: a list or value
        * sliceYCoords: a list or value
        * title (optional, defaults to '' )
        * legendOn (optional, defaults to True)
        * outFile (optional)
        
    """

    # Process the function arguments
    if isinstance(title, TextSettings):
        theTitle = title
    else:
        theTitle = TextSettings(title, fontSize=24)
    if not isinstance(xyAxisSettings, AxisSettings):
        xyAxisSettings = AxisSettings()
    if not isinstance(zAxisSettings, AxisSettings):
        zAxisSettings = AxisSettings()
    if sliceXCoords is None and sliceYCoords is None:
        raise ValueError("sliceXCoords or sliceYCoords must have a list or value, otherwise where do I slice?")
    if sliceXCoords is not None and sliceYCoords is not None:
        raise ValueError("sliceXCoords or sliceYCoords both have a list or value, which one do I slice?")

    if not isinstance(datasets, (list, tuple)):
        datasets = [datasets]
    _datasets = []
    for dataset in datasets:
        dimension = dataset.dimension
        if dimension == 2:
            _datasets.append(dataset)
        elif dimension == 3:
            if isinstance(dataset, DataSet3d):
                asDataSet3d = dataset
                likeW_XYs = dataset.getW_XYs()
            else:
                asDataSet3d = DataSet3d(dataset)
                likeW_XYs = dataset
            if sliceXCoords is not None:
                for x in sliceXCoords:
                    _datasets.append(
                        DataSet2d(likeW_XYs.getValue(x), legend=asDataSet3d.legend + ' @ ' + str(x) + ' ' + sliceUnits))
            if sliceYCoords is not None:
                for y in sliceYCoords:
                    _datasets.append(DataSet2d([[x, likeW_XYs.evaluate(x).evaluate(y)] for x in asDataSet3d.x],
                                               legend=asDataSet3d.legend + ' @ ' + str(y) + ' ' + sliceUnits))
        else:
            raise TypeError('Can only put 2d or 3d objects in slice plots')
    __makePlot2d(_datasets, xyAxisSettings, zAxisSettings, theTitle=theTitle, legendOn=legendOn, outFile=outFile,
                 legendXY=legendXY, figsize=figsize, useBokeh=useBokeh)


# --------------------------------------------------------
#
#  testing
#
# --------------------------------------------------------
def plotTests(_tests=(False, False, False, False, False, False, False, False, False, False, False )):
    from brownies.legacy.endl import fudgeParameters
    fudgeParameters.VerboseMode = 0
    from brownies.legacy.endl.endlProject import endlProject
    from brownies.legacy.endl.endl3dmathClasses import endl3dmath
    from fudge import __path__
    from xData import axes as axesModule
    from xData import XYs1d as XYs1dModule
    from xData import multiD_XYs as multiD_XYsModule

    testData = """
        #  Authors:   T.S.Suzuki, Y.Nagai, T.Shima, T.Kikuchi, H.Sato, T.Kii, M.Igashira
        #  Title:     First Measurement Of A P(N,Gamma)D Reaction Cross Section Between 10 And 80 Kev
        #  Year:      1995
        #  Institute: Tokyo Inst.of Technology, Tokyo
        #  Reference: Astrophysical Journal 439, (L), 59 (1995)
        #  Subent:    22310002
        #  Reaction:  Cross section for 1H(n,gamma)2H 
        # Note: the d(Energy) errorbars are fake and are used for plot testing
        #   Energy        Data          d(Energy)     d(Data)       
        #   MeV           barns         MeV           barns         
            0.02          0.000353      0.001         2.6e-05       
            0.02          0.000329      0.001         2.6e-05       
            0.02          0.000287      0.0           2.2e-05       
            0.02          0.000304      0.0           1.8e-05       
            0.04          0.00023       0.001         1.5e-05       
            0.04          0.000198      0.001         1.2e-05      
            0.04          0.000177      0.0015        1e-05         
            0.04          0.000207      0.0           1.3e-05       
            0.064         0.000156      0.0           1.1e-05       
            0.064         0.00015       0.0           7e-06         
            0.064         0.000158      0.0           1.1e-05       
            0.064         0.00014       0.0           9e-06
    """

    e = endlProject(__path__[0] + "/legacy/endl/test/testdb")
    za = e.readZA(1001)
    za.read()
    xAxis = AxisSettings(isLog=True, label='$E_n$', axisMin=0.5e-2, axisMax=1.5e-1, gridOn=True, autoscale=False,
                         unit='MeV')
    yAxis = AxisSettings(isLog=True, label='$\\sigma(E_n)$', gridOn=True, unit='b')
    d = []
    u = []
    for line in testData.split('\n'):
        if line.strip().startswith('#'):
            continue
        sline = [float(x) for x in line.split()] #list(map(float, line.split()))
        if len(sline) != 0:
            d.append(sline[0:2])
            u.append(sline[2:4])

    xyAxes = axesModule.Axes(labelsUnits={1: (xAxis.label, xAxis.unit), 0: (yAxis.label, yAxis.unit)})

    # Simple test, here we make a plot of one set

    if _tests[0]:
        xSec = za.findData(I=0, C=46)
        makePlot2d([xSec], xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,\\gamma)$ Cross Section',
                   outFile=None)
        xys = XYs1dModule.XYs1d(xSec.data, axes=xyAxes)
        makePlot2d([xys], xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,\\gamma)$ Cross Section',
                   outFile=None)
        dataset = DataSet2d(xys, xUnit=xAxis.unit, yUnit=yAxis.unit)
        dataset.convertUnits('eV', 'mb')
        makePlot2d([dataset], xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,\\gamma)$ Cross Section',
                   outFile=None)

    # Plot all the cross section data in the fudge2 test library, but unthemed!
    if _tests[1]:
        makePlot2d(za.findDatas(I=0), outFile=None)

    # Plot all the cross section data in the fudge2 test library
    if _tests[2]:
        xSecs = za.findDatas(I=0)
        makePlot2d(xSecs, xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,*)$ Cross Sections', outFile=None)
        xySecs = [XYs1dModule.XYs1d(xSec.data, axes=xyAxes) for xSec in xSecs]
        xySecs = (xySecs[0], xySecs[1], xySecs[2])
        makePlot2d(xySecs, xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,*)$ Cross Sections', outFile=None)

        # Fancy test, here we make a plot of one dataset (the testData above)
    if _tests[3]:
        endfData = za.findData(I=0, C=46)  # .slicex(domainMin=1e-2,domainMax=1e-1)
        endfUnc = 0.1 * endfData  # 10% error bars
        makePlot2d([
            DataSet2d(data=endfData, uncertainty=endfUnc, legend='ENDF/B-VII.0', color='g', lineStyle='-'),
            DataSet2d(data=d, uncertainty=u, legend='T.S.Suzuki, et al. (1995) EXFOR entry # 22310002', color='g',
                      symbol='o'),
        ], xAxisSettings=xAxis, yAxisSettings=yAxis, title='$^1$H$(n,\\gamma)$ Cross Section', legendOn=True,
            outFile=None)

    # Contour test
    if _tests[4]:  # Contour plots are meaningful if there is only one dataset plotted, how do we enforce this?
        endfData = za.findData(yo=1, I=1, C=10)
        EAxis = AxisSettings(isLog=False, label='$E_n$ (MeV)', axisMin=1.0, axisMax=20.0, gridOn=True, autoscale=False,
                             unit='MeV')
        muAxis = AxisSettings(isLog=False, label='$\\mu$ = cos( $\\theta$ )', axisMin=-1.0, axisMax=1.0,
                              autoscale=False,gridOn=True)
        makePlot2dContour(DataSet3d(data=endfData, legend='ENDF/B-VII.0', xUnit=EAxis.unit, yUnit=muAxis.unit),
                          xAxisSettings=EAxis, yAxisSettings=muAxis, title='$^1$H$(n,el)$ Angular Distribution',
                          outFile=None)
        w_xys = multiD_XYsModule.XYs2d(axes=axesModule.Axes(3, labelsUnits={2: ('$E_n$', 'MeV')}))
        for w, xy in endfData.data:
            w_xys.append(XYs1dModule.XYs1d(xy, axes=axesModule.Axes(2), outerDomainValue=w))
        dataset = DataSet3d(data=w_xys, legend='ENDF/B-VII.0')
        dataset.convertUnits('eV', None, None)
        makePlot2dContour(dataset, xAxisSettings=EAxis, yAxisSettings=muAxis,
                          title='$^1$H$(n,el)$ Angular Distribution', outFile=None)

    # Slice tests
    if _tests[5] or _tests[6] or _tests[7] or _tests[8] or _tests[9]:
        endfDataI0 = za.findData(yo=0, I=0, C=10).slicex(1.0, 20.0)  # simplify the plot
        endfDataI1 = za.findData(yo=1, I=1, C=10)
        Es = [p[0] for p in endfDataI0.data]  # in [MeV]
        mus = []  # in [mu]
        for t in endfDataI1.data:
            mus += [p[0] for p in t[1]]
        mus = sorted(uniquify(mus))
        table = []
        for E in Es:
            muDist = []
            for mu in mus:
                muDist.append([mu, endfDataI0.getValue(E) * endfDataI1.getValue(E, mu)])
            table.append([E, muDist])
        endfDataCombined = endl3dmath(data=table)
        EAxis = AxisSettings(isLog=True, label='$E_n$ (MeV)', axisMin=1.0, axisMax=20.0, gridOn=True, autoscale=False)
        muAxis = AxisSettings(isLog=False, label='$\\mu = \\cos{( \\theta )}$', axisMin=-1.0, axisMax=1.0,
                              autoscale=False, gridOn=True)
        zAxis = AxisSettings(isLog=True, label='$d\\sigma(E)/d\\mu$ (b)', axisMin=0.0, axisMax=5.0, autoscale=False,
                             gridOn=True)

        # unthemed slice tests
        if _tests[5]:
            makePlot2dSlice(endfDataCombined, xyAxisSettings=EAxis, zAxisSettings=zAxis, sliceXCoords=None,
                            sliceYCoords=[-1.0, -0.75, -0.5, -.25, 0.0, .25, .5, .75, 1.0], title='', outFile=None)
            dataSet3d = DataSet3d(data=endfDataCombined, legend='ENDF/B-VII.0')
            makePlot2dSlice(dataSet3d.getW_XYs(), xyAxisSettings=EAxis, zAxisSettings=zAxis, sliceXCoords=None,
                            sliceYCoords=[-1.0, -0.75, -0.5, -.25, 0.0, .25, .5, .75, 1.0], title='', outFile=None)
        if _tests[6]:
            makePlot2dSlice(endfDataCombined, xyAxisSettings=muAxis, zAxisSettings=zAxis,
                            sliceXCoords=[1.0, 5.0, 10.0, 15.0, 20.0], sliceYCoords=None, title='',  outFile=None)

        # themed slice test and contour test
        if _tests[7]:
            makePlot2dContour(DataSet3d(data=endfDataCombined, legend='ENDF/B-VII.0'), xAxisSettings=EAxis,
                                        yAxisSettings=muAxis, numContours=10,
                                        title='$d\\sigma(E)/d\\mu$ for $^1$H$(n,el)$', outFile=None)
        if _tests[8]:
            makePlot2dSlice(DataSet3d(data=endfDataCombined, legend='ENDF/B-VII.0'), xyAxisSettings=EAxis,
                                      zAxisSettings=zAxis, sliceXCoords=None,
                                      sliceYCoords=[-1.0, -0.75, -0.5, -.25, 0.0, .25, .5, .75, 1.0], title='',
                                      outFile=None)

        # themed slice test, with test 2d experimental data
        if _tests[9]:
            makePlot2dSlice([
                DataSet3d(data=endfDataCombined, legend='ENDF/B-VII.0'),
                DataSet2d(data=d, uncertainty=u, legend='T.S.Suzuki, et al. (1995) EXFOR entry # 22310002', color='g',
                          symbol='o'),
                ], xyAxisSettings=muAxis, zAxisSettings=zAxis, sliceXCoords=[1.0, 5.0, 10.0, 15.0, 20.0],
                sliceYCoords=None, sliceUnits='MeV', title='Slice test', outFile=None)


# --------------------------------------------------------
#
# main!
#
# --------------------------------------------------------
if __name__ == "__main__":

    for i in range(10):
        #        if( i != 8 ) : continue
        tests = 10 * [False]
        tests[i] = True
        if i in [4, 5, 9]:
            print('Test %d needs to be fixed.' % i)
            continue  # Dave, fix me.
        plotTests(tests=tests)
