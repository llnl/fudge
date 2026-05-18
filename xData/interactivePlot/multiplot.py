# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Interactive plotting with PyQt5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main entry is the class MultiPlotWithPyQt5 which requires input arguments for the plot title, x-axis label, y-axis
label, and a dictionary which contains the plot data. The latter dictionary is expected to have the plot labels as keys
and each value should be a tuple of two lists. The list at index 0 is expected to the plot's independent axis data,
while the list at index 1 is expected to contain the dependent axis data.

Sample Use:
MultiPlotWithPyQt5(title, xLabel, yLabel, plotData)
"""

import matplotlib # Weird bug on LC with Python 3.8.2 if matplotlib is imported after PyQt5

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLineEdit, QLabel, QCheckBox, QVBoxLayout, QTabWidget, \
    QDialog, QScrollArea, QComboBox, QGridLayout, QGroupBox
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QColor, QPalette, QFont, QIntValidator, QDoubleValidator

from collections import OrderedDict

import re
import sys
import numpy
import warnings
import math

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import pyplot
from matplotlib.transforms import Bbox
matplotlib.use('Qt5Agg')

from pqu import PQU as PQUModule


# ignore matplotlib warning for loc="best" since it is only used for cases with less than 10 plots
warnings.filterwarnings("ignore", message=".*Creating legend with loc.*")

def niceInitialLimits(left, right):

    left1, right1 = min(right, left), max(right, left)
    delta = right1 - left1
    if delta == 0:
        if right == 0:
            return -1, 1
        right = 1.01 * right
        left1 = 0.99 * left1
        left1, right1 = min(right1, left1), max(right, left1)        # Needed for negative value.
        delta = right1 - left1

    numberOfDigits = 3 + math.floor(math.log10(max(abs(right1), abs(left1)) / delta))

    left2 = float(PQUModule.floatToShortestString(left1, numberOfDigits))
    if left2 > left1:
        adjust = '%.0e' % (left2 - left1)
        left2 = left1 - float('9.999e' + adjust.split('e')[-1])
        left2 = float(PQUModule.floatToShortestString(left2, numberOfDigits))

    right2 = float(PQUModule.floatToShortestString(right1, numberOfDigits))
    if right2 < right1:
        adjust = '%.0e' % (right1 - right2)
        right2 = right1 + float('5e' + adjust.split('e')[-1])
        right2 = float(PQUModule.floatToShortestString(right2, numberOfDigits))

    return left2, right2

class MultiPlotWithPyQt5:
    """
    Main class to plot cross sections with PyQt5 which in turn calls the main PyQt5 application window class

    :param plotAttributes: Dictionary with the plot title and axis labels as entries.
    :param plotData: Dictionary with plot labels as keys and each value as a tuple of x-value list and y-value list.
    :param plotType: Variable indicating the plot type (i.e. '2d' or '3d') (string)
    """
    def __init__(self, plotAttributes, plotData, plotType='2d'):

        plotWindows = MainWindow(plotAttributes, plotData, plotType)
        plotWindows.applicationObject.exec()

        # explicitly delete the reference to QMainWindow to avoid segmentation fault 11
        plotWindows.deleteLater()


class MainWindow(QMainWindow):
    """
    The main PyQt5 application window class which is derived from PyQt5.QtWidgets.QMainWindow.

    this class initializes the PyQt5.QtWidgets.QApplication which manages the GUI application's control flow and main
    settings. It initially sets the application style to "fusion", i.e. QApplication::setStyle("fusion") which is
    required to address OSX idiosyncrasies with the use of colors with PyQt5.QtWidgets.QComboBox.

    After setting the geometry it launches a main window (which embeds the matplotlib plot) and an auxiliary window
    which contains controls for plot attributes.

    :param plotAttributes: Dictionary with the plot title and axis labels as entries.
    :param plotData: Dictionary with plot labels as keys and each value as a tuple of x-value list and y-value list.
    :param plotType: Variable indicating the plot type (i.e. '2d' or '3d') (string)
    """

    def __init__(self, plotAttributes, plotData, plotType):
        self.applicationObject = QApplication([])
        self.applicationObject.setStyle("fusion")  # required to have QCombobox colors on osx

        super().__init__()

        # initialize plot window
        self.setWindowTitle(plotAttributes['title'])

        self.plotType = plotType
        magnification = 1 if plotType == '2d' else 1.5

        # window dimensions from default pyplot parameters
        dpi = 100
        widthInches, heightInches = pyplot.rcParams["figure.figsize"]
        width = int(widthInches*dpi*magnification)
        height = int(heightInches*dpi*magnification)
        self.setGeometry(0, 0, width, height)

        # generate plot
        self.plotInstance = PlotCanvas(self, widthInches, heightInches, plotAttributes, plotData, dpi, plotType)
        self.setCentralWidget(self.plotInstance)

        self.addToolBar(NavigationToolbar(self.plotInstance, self))

        self.show()

        x0, y0 = width, 0
        width = 800
        height = 300
        self.dialogWindow = DialogWindow(width, height, plotAttributes, plotObject=self.plotInstance, plotType=plotType,
                                         parent=self)

        self.dialogWindow.setGeometry(x0, y0, width, height)
        self.dialogWindow.show()

    def closeEvent(self, event):
        # ensure dialog window is closed
        if self.plotType == '2d':
            self.dialogWindow.accept()


class PlotCanvas(FigureCanvas):
    """
    The canvas the plot figure renders to.

    Class derived from matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg to embed the the plot. It loops through the
    entries in the plot data dictionary and uses the keys as plot labels. The values of the plot data dictionary are
    expected to be a tuple withthe list at index 0 represeting the independent ("x") data and the list at index 1, the
    dependent data ("y").

    For cases with less than 10 plots the legend is placed at the matplotlib "best" location and is allowed to be
    draggable. For cases with more than 10 plots the legend is fixed to the left of the plot window and is allowed to be
    scrollable.

    :param parent: The parent class which in this case is the PyQt5 application window class.
    :param width: Window width in inches, derived from pyplot.rcParams["figure.figsize"] in the parent class (float).
    :param height: Window height in inches, derived from pyplot.rcParams["figure.figsize"] in the parent class (float).
    :param plotAttributes: Dictionary with the plot title and axis labels as entries.
    :param plotData: Dictionary with plot labels as keys and each value as a tuple of x-value list and y-value list.
    :param dpi: Figure resolution (float).
    :param plotType: Variable indicating the plot type (i.e. '2d' or '3d') (string)
    """

    def __init__(self, parent, width, height, plotAttributes, plotData, dpi, plotType):

        self.plotFigure = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.plotFigure)
        self.setParent(parent)

        # not sure what this does ... keep it around anyway
        # FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        # FigureCanvas.updateGeometry(self)

        availablePlots = ['2d', '3d']
        assert plotType in availablePlots, \
            'Unrecognized plot type ... expected one of the following: %s' % ', '.join(availablePlots)

        self.plotAxis = None
        self.gridOn = None
        self.dataLimits = None
        if plotType == '2d':
            self.plot2d(plotData, plotAttributes)
            assert self.plotAxis is not None
        
            self.legendDraggableNotScrollable = len(self.plotAxis.lines) < 10
            self.percentPlotWidth = 1.0 if self.legendDraggableNotScrollable else 0.8
            self.legend = self.addLegend(self.legendDraggableNotScrollable)

        elif plotType == '3d':
            self.plot3d(plotData, plotAttributes)

        self.draw()

    def plot2d(self, plotData, plotAttributes):

        self.plotAxis = self.figure.add_subplot(111)

        floatInfo = numpy.finfo(numpy.float64())
        self.dataLimits = {'x': [floatInfo.max, floatInfo.min], 'y': [floatInfo.max, floatInfo.min]}
        for plotLegend in plotData.keys():
            if (not isinstance(plotData[plotLegend], (list, tuple))) or (len(plotData[plotLegend]) != 2):
                warnings.warn('No plot for legend "%s" due to incorrect dictionary value object type' % plotLegend)
                continue

            xValues = plotData[plotLegend][0]
            self.dataLimits['x'] = [min(self.dataLimits['x'][0], min(xValues)),
                                    max(self.dataLimits['x'][1], max(xValues))]
            if (not isinstance(xValues, (list, tuple))) or (len(xValues) == 0):
                warnings.warn('No plot for legend "%s" due to zero length xValue array' % plotLegend)
                continue

            yValues = plotData[plotLegend][1]
            self.dataLimits['y'] = [min(self.dataLimits['y'][0], min(yValues)),
                                    max(self.dataLimits['y'][1], max(yValues))]
            if len(yValues) == 0:
                warnings.warn('No plot for legend "%s" due to zero length yValue array' % plotLegend)
                continue

            self.plotAxis.plot(xValues, yValues, label=plotLegend)

        xleft = self.dataLimits['x'][0] if plotAttributes.get('xMin', '') == '' else float(plotAttributes['xMin'])
        xright = self.dataLimits['x'][1] if plotAttributes.get('xMax', '') == '' else float(plotAttributes['xMax'])
        xleft, xright = niceInitialLimits(xleft, xright)
        self.plotAxis.set_xlim(xleft, xright)
        plotAttributes['xMin'] = str(xleft)
        plotAttributes['xMax'] = str(xright)

        yleft = self.dataLimits['y'][0] if plotAttributes.get('yMin', '') == '' else float(plotAttributes['yMin'])
        yright = self.dataLimits['y'][1] if plotAttributes.get('yMax', '') == '' else float(plotAttributes['yMax'])
        yleft, yright = niceInitialLimits(yleft, yright)
        self.plotAxis.set_ylim(yleft, yright)
        plotAttributes['yMin'] = str(yleft)
        plotAttributes['yMax'] = str(yright)

        if plotAttributes.get('title', '') != '':
            self.plotAxis.set_title(plotAttributes['title'])

        self.plotAxis.set_xlabel(plotAttributes['xLabel'])
        self.plotAxis.set_ylabel(plotAttributes['yLabel'])

        logs = plotAttributes.get('logs', 0)
        xlog = logs % 2
        ylog = ( logs // 2 ) % 2
        linlog = {0: 'linear', 1: 'log'}
        self.plotAxis.set_xscale(linlog[xlog])
        self.plotAxis.set_yscale(linlog[ylog])

        self.plotAxis.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
        self.plotAxis.tick_params(which='minor', bottom=False, top=False, left=False, right=False)
        self.gridOn = True
        self.plotAxis.grid(self.gridOn, linestyle='dashed')

    def plot3d(self, plotData, plotAttributes):
        from matplotlib import cm as colorMap

        self.plotAxis = self.figure.add_subplot(111, projection='3d')

        assert len(plotData) == 1, '3d Plots limited to a single plot'
        plotLegend = [key for key in plotData.keys()][0]
        if (not isinstance(plotData[plotLegend], (list, tuple))) or (len(plotData[plotLegend]) != 3):
            warnings.warn('No plot for legend "%s" due to incorrect dictionary value object type' % plotLegend)

        xValues = plotData[plotLegend][0]
        self.dataLimits = {'x': [xValues.min(), xValues.max()]}
        if len(xValues) == 0:
            warnings.warn('No plot for legend "%s" due to zero length xValue array' % plotLegend)

        yValues = plotData[plotLegend][1]
        self.dataLimits['y'] = [yValues.min(), yValues.max()]
        if len(yValues) == 0:
            warnings.warn('No plot for legend "%s" due to zero length yValue array' % plotLegend)

        zValues = plotData[plotLegend][2]
        self.dataLimits['z'] = [zValues.min(), zValues.max()]
        if len(zValues) == 0:
            warnings.warn('No plot for legend "%s" due to zero length zValue array' % plotLegend)

        surfacePlot = self.plotAxis.plot_surface(xValues, yValues, zValues, cmap=colorMap.jet, linewidth=0,
                                                 antialiased=False)

        self.plotAxis.set_xlabel(plotAttributes['xLabel'])
        self.plotAxis.set_ylabel(plotAttributes['yLabel'])
        self.plotAxis.set_zlabel(plotAttributes['zLabel'])

        self.gridOn = True
        self.plotAxis.grid(self.gridOn, linestyle='dashed')

        self.plotAxis.view_init(25, 15)
        self.plotFigure.colorbar(surfacePlot, shrink=0.5, aspect=10)

    def addLegend(self, draggableLegend):
        """
        Method to create plot legend.

        The variable self.legendDraggableNotScrollable is used to decide whether the legend should be inside or outside
        the plot region.

        If self.legendDraggableNotScrollable == True:
          - Legend is placed inside the plot region;
          - Plot region maintains the default size;
          - Depending on the value of the draggableLegend variable, the legend location is fixed
            (draggableLegend == False) or draggable (draggableLegend == True); and
          - Legend is never scrollable

        If self.legendDraggableNotScrollable == False:
          - Plot region is shrinked and legend is placed on the outside; and
          - Legend is allowed to be scrollable

        :param draggableLegend: Variable to set legend draggable (boolean) if self.legendDraggableNotScrollable is True
        :return: matplotlib.axes.legend object
        """

        if self.legendDraggableNotScrollable:
            if draggableLegend:
                # draggable legend
                legend = self.plotAxis.legend(loc='best')
                legend.set_draggable(True)

            else:
                transformation = self.legend.axes.transAxes.inverted()
                boundingBox = self.legend.get_bbox_to_anchor().transformed(transformation)
                location = self.legend.get_window_extent().transformed(transformation)
                legend = self.plotAxis.legend(loc=(location.x0, location.y0), bbox_to_anchor=boundingBox)

        else:
            # Shrink current axis by 20%
            box = self.plotAxis.get_position()
            self.plotAxis.set_position([box.x0, box.y0, box.width * self.percentPlotWidth, box.height])
            self.percentPlotWidth = 1.0  # ensure that axis is only shrinked once

            legend = self.plotAxis.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
            self.mpl_connect("scroll_event", self.legendScroll)

        return legend

    def redrawPlot(self):
        """
        Method to redraw plot and is normally called after a plot property is changed.
        """

        draggableLegend = self.legend.get_draggable()
        self.legend = self.addLegend(draggableLegend)
        self.plotFigure.canvas.draw_idle()

    def legendScroll(self, event):
        """
        Method to enable scrolling of the plot legend.

        Useful for cases with a large number of plots (currently implemented for cases with more than 10 plots) and
        based on a response to a stackoverflow question (
        see https://stackoverflow.com/questions/55863590/adding-scroll-button-to-matlibplot-axes-legend accessed on
        03/24/2020).


        :param event: Mouse scroll event (type matplotlib.backend_bases.MouseEvent).
        """
        # pixels to scroll per mousewheel event
        d = {"down": 30, "up": -30}

        if self.legend.contains(event):
            boundingBox = self.legend.get_bbox_to_anchor()
            boundingBox = Bbox.from_bounds(boundingBox.x0, boundingBox.y0 + d[event.button], boundingBox.width,
                                           boundingBox.height)
            transformation = self.legend.axes.transAxes.inverted()

            # limit scrolling to keep legend visible
            windowExtent = self.legend.get_window_extent().transformed(transformation)
            legendScrollable = (event.button == 'up' and windowExtent.y1 >= 1) or \
                               (event.button == 'down' and windowExtent.y0 < 0)
            if legendScrollable:
                self.legend.set_bbox_to_anchor(boundingBox.transformed(transformation))
                self.plotFigure.canvas.draw_idle()


class DialogWindow(QDialog):
    """
    The dialogue window with user controls for plot attributes.

    Dialogue window derived from PyQt5.QtWidgets.QDialog. It generates separate tabs for overall plot options and
    individual plot data options.

    The overall options include controls for the font size, labels, data ranges and the ability to toggle the legend
    display.

    The individual plot data options allows for the control of the plot display; plot label; plot line color, type and
    width; as well as the plot marker type and size. This tab is allowed to be scrollable to allow for the control of a
    large number of individual plots.

    :param width: Window width (float).
    :param height: Window height (float).
    :param plotObject: Plot object to be passed on to the children classes to access the plot figure and axis handles.
    :param parent: The parent class which in this case is the MainWindow class.
    """

    def __init__(self, width, height, plotAttributes, plotObject, plotType, parent=None):
        super().__init__(parent)

        self.parentWindow = parent

        self.setMinimumSize(QSize(width, height))
        self.setWindowTitle('Fudge Interactive Plot Controls')

        tabWidget = QTabWidget()
        tabWidget.addTab(OverallPlotOptions(plotObject, plotAttributes, plotType), 'Overall Plot Options')
        if plotType == '2d':
            tabWidget.addTab(IndividualPlotOptions(plotObject), 'Individual DataSet Options')

        vboxLayout = QVBoxLayout()
        vboxLayout.addWidget(tabWidget)

        self.setLayout(vboxLayout)

    def closeEvent(self, event):
        # close main window
        self.parentWindow.applicationObject.quit()


class OverallPlotOptions(QWidget):
    """
    Class derived from PyQt5.QtWidgets.QWidget to generate the controls for the overall plot options.

    Include controls for the font size, labels, data ranges and the ability to toggle the legend display.

    :param plotObject: Plot object to be passed on to the children classes to access the plot figure and axis handles.
    """

    def __init__(self, plotObject, plotAttributes, plotType):
        self.plotInstance = plotObject
        self.axisHandle = plotObject.plotAxis
        self.figureHandle = plotObject.plotFigure

        super().__init__()

        self.userUpdatedRange = False

        # Font size
        yCoordinate = 20
        self.fontSizeLabel = self.createTextBoxLabel('     font size = ', yCoordinate)
        self.fontSize = self.createTextBox('%d' % self.axisHandle.title.get_fontsize(), yCoordinate, 100, self.updateFontSize,
                                           setValidator=QIntValidator(2, 60))

        # Title
        deltaYCoordinate = 30
        yCoordinate += deltaYCoordinate
        self.plotTitleLabel = self.createTextBoxLabel('              title = ', yCoordinate)
        self.plotTitle = self.createTextBox(self.axisHandle.get_title(), yCoordinate, 300, self.updatePlotTitle)

        # X axislabel
        yCoordinate += deltaYCoordinate
        self.xAxisLabelLabel = self.createTextBoxLabel('         x label = ', yCoordinate)
        self.xAxisLabel = self.createTextBox(self.axisHandle.get_xlabel(), yCoordinate, 300, self.updateXAxisLabel)

        # Y axislabel
        deltaYCoordinate = 25
        yCoordinate += deltaYCoordinate
        self.yAxisLabelLabel = self.createTextBoxLabel('         y label = ', yCoordinate)
        self.yAxisLabel = self.createTextBox(self.axisHandle.get_ylabel(), yCoordinate, 300, self.updateYAxisLabel)
        deltaYCoordinateZLabel = deltaYCoordinate
        yCoordinateZLabel = yCoordinate + deltaYCoordinateZLabel

        # X range
        deltaYCoordinate = 35
        maxSystemFloat = sys.float_info.max
        maxNumberDigits = 16
        yCoordinate += deltaYCoordinate if plotType == '2d' else deltaYCoordinate + deltaYCoordinateZLabel
        self.xRangeLabel = self.createTextBoxLabel('        x range = ', yCoordinate)
        self.xRangeLower = self.createTextBox(plotAttributes['xMin'], yCoordinate, 100, self.updateXRange, \
            setValidator=QDoubleValidator(-maxSystemFloat, maxSystemFloat, maxNumberDigits) )
        self.xRangeUpper = self.createTextBox(plotAttributes['xMax'], yCoordinate, 100, self.updateXRange, \
            setValidator=QDoubleValidator(-maxSystemFloat, maxSystemFloat, maxNumberDigits), moveX=205 )
        self.xLog = self.createCheckBox('xlog', 350, yCoordinate, self.axisHandle.get_xscale() == 'log', self.toggleXScale)
        self.xrange = plotObject.dataLimits['x']

        # update textbox if zoom changes
        self.axisHandle.callbacks.connect('xlim_changed', self.onXLimitChange)

        # Y range
        deltaYCoordinate = 25
        yCoordinate += deltaYCoordinate
        self.yRangeLabel = self.createTextBoxLabel('        y range = ', yCoordinate)
        self.yRangeLower = self.createTextBox(plotAttributes['yMin'], yCoordinate, 100, self.updateYRange, \
            setValidator=QDoubleValidator(-maxSystemFloat, maxSystemFloat, maxNumberDigits) )
        self.yRangeUpper = self.createTextBox(plotAttributes['yMax'], yCoordinate, 100, self.updateYRange, \
            setValidator=QDoubleValidator(-maxSystemFloat, maxSystemFloat, maxNumberDigits), moveX=205 )
        self.yLog = self.createCheckBox('ylog', 350, yCoordinate, self.axisHandle.get_yscale() == 'log',
                                        self.toggleYScale)

        self.yrange = plotObject.dataLimits['y']

        # update textbox if zoom changes
        self.axisHandle.callbacks.connect('ylim_changed', self.onYLimitChange)

        # Grid checkbox
        self.gridCheckBox = self.createCheckBox('Grid Lines', 600, 80, self.plotInstance.gridOn, self.toggleGrid)

        if plotType == '2d':
            self.legendDraggableNotScrollable = plotObject.legendDraggableNotScrollable

            # Legend checkbox
            self.legendVisibleCheckBox = \
                self.createCheckBox('Show Legend', 600, 20, self.plotInstance.legend.get_visible(), self.toggleLegend)

            # checkbox to fix legend location
            if self.legendDraggableNotScrollable:
                self.legendLocationFixedCheckBox = \
                    self.createCheckBox('Fix Legend Location', 600, 40, not self.legendDraggableNotScrollable,
                                        self.fixLegendLocation)

        else:
            # Z axislabel
            self.zAxisLabelLabel = self.createTextBoxLabel('         z label = ', yCoordinateZLabel)
            self.zAxisLabel = self.createTextBox(self.axisHandle.get_zlabel(), yCoordinateZLabel, 300,
                                                 self.updateZAxisLabel)

            # dropdown list for plot type
            # self.addDropDownList(100, 'Plot Type', Surface Plot', )



    def createTextBoxLabel(self, labelText, moveY, moveX=10):
        """
        General method to create a PyQt5 label, i.e. PyQt5.QtWidgets.QLabel

        :param labelText: Text to be placed in the label (string).
        :param moveY: Vertical label location (float).
        :param moveX: Horizontal label location (float).
        """

        labelBox = QLabel(self)
        labelBox.move(moveX, moveY)
        labelBox.setText(labelText)

        return labelBox

    def createTextBox(self, defaultText, moveY, sizeX, updateMethod, setValidator=None, moveX=100, sizeY=20):
        """
        General method to create a PyQt5 text box, i.e. PyQt5.QtWidgets.QLineEdit

        :param defaultText: Initial text to be placed in text box (string).
        :param moveY: Vertical textbox location (float).
        :param sizeX: Horizontal textbox size (float).
        :param updateMethod: Method to use after ENTER is pressed in textbox.
        :param setValidator: Input to validate input text.
        :param moveX: Horizontal textbox location (float).
        :param sizeY: Vertical textbox size (float).
        """

        textBox = QLineEdit(self)

        if setValidator is not None:
            textBox.setValidator(setValidator)

        textBox.move(moveX, moveY)
        textBox.resize(sizeX, sizeY)
        textBox.setText(defaultText)
        textBox.returnPressed.connect(updateMethod)

        return textBox

    def createCheckBox(self, label, moveX, moveY, checkRequired, updateMethod):
        """
        General method to create a PyQt5 checkbox, i.e. PyQt5.QtWidgets.QCheckBox

        :param label: Checkbox label (string).
        :param moveX: Horizontal checkbox size (float).
        :param moveY: Vertical checkbox size (float).
        :param checkRequired: Default check state (bool).
        :param updateMethod: Method to use after checkbox state change.
        """
        checkBox = QCheckBox(label, self)
        checkBox.move(moveX, moveY)

        if checkRequired:
            if not checkBox.isChecked():
                checkBox.toggle()
        else:
            if checkBox.isChecked():
                checkBox.toggle()

        checkBox.stateChanged.connect(updateMethod)

        return checkBox

    def updateFontSize(self):
        """
        Method to update font size for the plot title, x-axis label, y-axis labels, x-tick labels, and y-tick labels.
        """

        newFontSize = self.fontSize.text()
        axisList = [self.axisHandle.title]
        axisList += [getattr(self.axisHandle, axisObject).label for axisObject in dir(self.axisHandle)
                     if re.match('^[xyz]axis$', axisObject)]
        for axisObject in dir(self.axisHandle):
            if re.match('^get_.ticklabels$', axisObject):
                axisList += getattr(self.axisHandle, axisObject)()

        for axisItem in axisList:
            axisItem.set_fontsize(newFontSize)

        self.axisHandle.axes.tick_params(labelsize=newFontSize)
        self.axisHandle.xaxis.get_offset_text().set_fontsize(newFontSize)
        self.axisHandle.yaxis.get_offset_text().set_fontsize(newFontSize)
        if hasattr(self.axisHandle, 'zaxis'):
            self.axisHandle.zaxis.get_offset_text().set_fontsize(newFontSize)

        self.figureHandle.canvas.draw_idle()

    def updateXAxisLabel(self):
        """
        Method x-axis label text.
        """

        newLabel = self.xAxisLabel.text()
        self.axisHandle.set_xlabel(newLabel)
        self.axisHandle.xaxis.label.set_fontsize(self.fontSize.text())
        self.figureHandle.canvas.draw_idle()

    def updateYAxisLabel(self):
        """
        Method y-axis label text
        """

        newLabel = self.yAxisLabel.text()
        self.axisHandle.set_ylabel(newLabel)
        self.axisHandle.yaxis.label.set_fontsize(self.fontSize.text())
        self.figureHandle.canvas.draw_idle()

    def updateZAxisLabel(self):
        """
        Method y-axis label text
        """

        newLabel = self.zAxisLabel.text()
        self.axisHandle.set_zlabel(newLabel)
        self.axisHandle.zaxis.label.set_fontsize(self.fontSize.text())
        self.figureHandle.canvas.draw_idle()

    def updatePlotTitle(self):
        """
        Method to update main plot title.
        """

        newTitle = self.plotTitle.text()
        self.axisHandle.set_title(newTitle)
        self.axisHandle.title.set_fontsize(self.fontSize.text())
        self.figureHandle.canvas.draw_idle()

    def onXLimitChange(self, axes):
        """
        Method to update the x-axis limits checkbox when the actual axis limit changes via panning or zooming.

        :param axes: The matplotlib axes variable.
        """

        if self.userUpdatedRange:
            return
            self.userUpdatedRange = False

        xLimits = axes.get_xlim()
        self.xRangeLower.setText('%G' % xLimits[0])
        self.xRangeUpper.setText('%G' % xLimits[1])

    def updateXRange(self):
        """
        Method to update x-axis range.
        """

        self.userUpdatedRange = True

        xrange = [self.xRangeLower.text(), self.xRangeUpper.text()]
        xrange[0] = float(xrange[0]) if xrange[0] != '*' else self.xrange[0]
        xrange[1] = float(xrange[1]) if xrange[1] != '*' else self.xrange[1]

        self.axisHandle.set_xlim(xrange)
        self.figureHandle.canvas.draw_idle()

    def onYLimitChange(self, axes):
        """
        Method to update the y-axis limits checkbox when the actual axis limit changes via panning or zooming.

        :param axes: The matplotlib axes variable.
        """

        if self.userUpdatedRange:
            self.userUpdatedRange = False
            return

        yLimits = axes.get_ylim()
        self.yRangeLower.setText('%G' % yLimits[0])
        self.yRangeUpper.setText('%G' % yLimits[1])

    def updateYRange(self):
        """
        Method to update y-axis range.
        """

        self.userUpdatedRange = True

        yrange = [self.yRangeLower.text(), self.yRangeUpper.text()]
        yrange[0] = float(yrange[0]) if yrange[0] != '*' else self.yrange[0]
        yrange[1] = float(yrange[1]) if yrange[1] != '*' else self.yrange[1]

        self.axisHandle.set_ylim(yrange)
        self.figureHandle.canvas.draw_idle()

    def toggleXScale(self):
        """
        Method to toggle the x-axis between a linear and log scale.
        """

        if self.xLog.isChecked():
            if self.axisHandle.get_xscale() != 'log':
                self.axisHandle.set_xscale('log')
        elif self.axisHandle.get_xscale() != 'linear':
            self.axisHandle.set_xscale('linear')

        self.figureHandle.canvas.draw_idle()

    def toggleYScale(self):
        """
        Method to toggle the y-axis between a linear and log scale.
        """

        if self.yLog.isChecked():
            if self.axisHandle.get_yscale() != 'log':
                self.axisHandle.set_yscale('log')
        elif self.axisHandle.get_yscale() != 'linear':
            self.axisHandle.set_yscale('linear')

        self.figureHandle.canvas.draw_idle()

    def toggleLegend(self):
        """
        Method to toggle legend visibility.
        """

        self.plotInstance.legend.set_visible(self.legendVisibleCheckBox.isChecked())
        self.figureHandle.canvas.draw_idle()

    def redrawPlot(self):
        """
        Local method that calls the PlotCanvas.redrawPlot method.
        """

        self.plotInstance.redrawPlot()

    def fixLegendLocation(self):
        """
        Method to fix draggable legend location
        """

        checkBoxState = self.legendLocationFixedCheckBox.isChecked()
        self.plotInstance.legend.set_draggable(not checkBoxState)

    def toggleGrid(self):
        """
        Method to toggle plot axes grid.
        """

        if self.gridCheckBox.isChecked():
            self.plotInstance.gridOn = True
            self.axisHandle.grid(True)

        else:
            self.plotInstance.gridOn = False
            self.axisHandle.grid(False)

        self.figureHandle.canvas.draw_idle()


class IndividualPlotOptions(QWidget):
    """
    Class derived from PyQt5.QtWidgets.QWidget to generate the controls for the individual plot options.

    Include controls for the individual plot display; label; line color, type and width; as well as the marker type and
    size. This tab is allowed to be scrollable to allow for the control of a large number of plots.

    :param plotObject: Plot object to be passed on to the children classes to access the plot figure and axis handles.
    """

    def __init__(self, plotObject):
        self.plotInstance = plotObject
        self.axisHandle = plotObject.plotAxis
        self.figureHandle = plotObject.plotFigure

        super().__init__()

        self.selectionNotAvailable = 'N/A'

        self.colorDictionary = dict([(colorHex.lower(), colorName.split(':')[-1])
                                     for colorName, colorHex in matplotlib.colors.TABLEAU_COLORS.items()])
        self.markerStyleDictionary = \
            OrderedDict([('None', 'None'), ('.', 'Point'), ('o', 'Circle'), ('v', 'Triangle Down'),
                         ('^', 'Triangle Up'), ('<', 'Triangle Left'), ('>', 'Triangle Right'), ('s', 'Square'),
                         ('p', 'Pentagon'), ('*', 'Star'), ('+', 'Plus'), ('x', 'X'), ('D', 'Diamond')])
        self.lineStyleDictionary = OrderedDict([('-', 'Solid'), ('--', 'Dashed'), ('-.', 'Dash Dot'), (':', 'Dotted'),
                                                ('', 'None')])
        availableColors = [self.colorDictionary[hexKey] for hexKey in
                           pyplot.rcParams['axes.prop_cycle'].by_key()['color']]

        # headings
        self.xSize = 100
        self.ySize = 20
        self.xSizeLabel = 2 * self.xSize
        self.xSizeLineStyle = 1.5 * self.xSize
        self.xSizeLineWidth = 0.25 * self.xSize
        formLayout = self.headings()

        # rows of widgets for each individual plot
        self.widgetList = []
        for i, (line, legendLabel) in enumerate(zip(*self.axisHandle.get_legend_handles_labels())):
            irow = i + 1
            # row variables
            rowName = '%3.3d%s' % (i, legendLabel)

            # checkbox: toggle plot selection
            checkBox = self.addCheckBox(rowName, self.toggleSinglePlotVisibility)
            formLayout.addWidget(checkBox, irow, 0)
            self.widgetList.append(checkBox)

            # label: plot label used in legend
            widget = self.addTextBox(self.xSizeLabel, rowName, legendLabel, self.changePlotLabel)
            formLayout.addWidget(widget, irow, 1)
            self.widgetList.append(widget)

            # dropdown list for plot line color selection
            comboBox = self.addLineColorDropDown(rowName, line.get_color(), availableColors)
            formLayout.addWidget(comboBox, irow, 2)
            self.widgetList.append(comboBox)

            # dropdown list for plot line style selection
            comboBox = self.addDropDownList(self.xSizeLineStyle, rowName, line.get_linestyle(),
                                            self.lineStyleDictionary, self.changeLineStyle, irow, 'line stye')
            formLayout.addWidget(comboBox, irow, 3)
            self.widgetList.append(comboBox)

            # textbox for plot linewidth
            textBox = self.addTextBox(self.xSizeLineWidth, rowName, '%.2f' % line.get_linewidth(),
                                      self.changeLineWidth, setValidator=QDoubleValidator(0.0, 20.0, 3))
            formLayout.addWidget(textBox, irow, 4)
            self.widgetList.append(textBox)

            # dropdown list for marker style
            comboBox = self.addDropDownList(self.xSizeLineStyle, rowName, line.get_marker(),
                                            self.markerStyleDictionary, self.changeMarkerStyle, irow, 'marker style')
            formLayout.addWidget(comboBox, irow, 5)
            self.widgetList.append(comboBox)

            # textbox for plot marker size
            textBox = self.addTextBox(self.xSizeLineWidth, rowName, '%.2f' % line.get_markersize(),
                                      self.changeMarkerSize, setValidator=QDoubleValidator(0.0, 20.0, 3))
            formLayout.addWidget(textBox, irow, 6)
            self.widgetList.append(textBox)

        groupBox = QGroupBox()
        formLayout.setColumnStretch(1, 4)
        formLayout.setColumnStretch(2, 1)
        formLayout.setColumnStretch(3, 1)
        formLayout.setColumnStretch(4, 1)
        formLayout.setColumnStretch(5, 1)
        formLayout.setColumnStretch(6, 1)
        groupBox.setLayout(formLayout)
        scroll = QScrollArea()
        scroll.setWidget(groupBox)
        scroll.setWidgetResizable(True)

        layout = QVBoxLayout(self)
        layout.addWidget(scroll)

    def redrawPlot(self):
        """
        Local method that calls the PlotCanvas.redrawPlot method.
        """

        self.plotInstance.redrawPlot()

    def headings(self, skip=False):
        """
        Method to add headings above individual plot options ... currently not used since skip == True

        :param skip: Option to add heading above individual plot options (bool)
        """

        formLayout = QGridLayout()

        if skip:
            return formLayout

        formLayout.addWidget(self.addCheckBox('Label', self.toggleAllPlotVisibility), 0, 0)
        formLayout.addWidget(self.addLabel('Line Labels'), 0, 1)
        formLayout.addWidget(self.addLabel('Line Color'), 0, 2)
        formLayout.addWidget(self.addLabel('Line Style'), 0, 3)
        formLayout.addWidget(self.addLabel('Line Width'), 0, 4)
        formLayout.addWidget(self.addLabel('Marker Style'), 0, 5)
        formLayout.addWidget(self.addLabel('Marker Size'), 0, 6)

        return formLayout

    def addLabel(self, labelText):
        """
        General method to create a PyQt5 label, i.e. PyQt5.QtWidgets.QLabel

        :param labelText: Text to be placed in the label (string).
        """

        labelBox = QLabel(self)
        labelBox.setText(labelText)
        labelBox.setAlignment(Qt.AlignCenter)
        labeFont = QFont()
        labeFont.setBold(True)
        labelBox.setFont(labeFont)

        return labelBox

    def addCheckBox(self, objectName, connectMethod):
        """
        General method to create a PyQt5 checkbox, i.e. PyQt5.QtWidgets.QCheckBox

        :param objectName: Checkbox object name (string).
        :param connectMethod: Method to use after checkbox state change.
        """

        # noinspection PyArgumentList
        checkBox = QCheckBox('', self, objectName=objectName)
        if not checkBox.isChecked():
            checkBox.toggle()

        checkBox.stateChanged.connect(connectMethod)

        return checkBox

    def addTextBox(self, xSize, objectName, defaultText, connectMethod=None, setValidator=None):
        """
        General method to create a PyQt5 textbox, i.e. PyQt5.QtWidgets.QLineEdit

        :param xSize: Horizontal textbox size (float).
        :param objectName: Textbox object name (string).
        :param defaultText: Initial text to be placed in textbox (string).
        :param connectMethod: Method to use after ENTER is pressed in textbox.
        :param setValidator: Input to validate input text.
        """

        # noinspection PyArgumentList
        textBox = QLineEdit(self, objectName=objectName)
        textBox.resize(int(xSize), int(self.ySize))
        textBox.setText(defaultText)

        if setValidator is not None:
            textBox.setValidator(setValidator)

        if connectMethod is not None:
            textBox.returnPressed.connect(connectMethod)

        return textBox

    def addDropDownList(self, xSize, objectName, defaultValue, itemDictionary, connectMethod, plotIndex, plotAttribute):
        """
        General method to create a PyQt5 drop-down list, i.e. PyQt5.QtWidgets.QComboBox

        :param xSize: Horizontal drop-down list size (float).
        :param objectName: Drop-down list object name (string).
        :param defaultValue: Default selected value for drop-down list (string).
        :param itemDictionary: Ordered Dictionary with possible drop-down values.
        :param connectMethod: Method to use after selected drop-down value changes.
        :param plotIndex: Index of plot associated with the drop-down list.
        :param plotAttribute: Plot attribute associated with the drop-down list.
        """

        # noinspection PyArgumentList
        comboBox = QComboBox(self, objectName=objectName)
        for dropDownItem in itemDictionary.values():
            comboBox.addItem(dropDownItem)

        indexComboBox = comboBox.findText(itemDictionary[defaultValue], Qt.MatchFixedString)
        if indexComboBox >= 0:
            comboBox.setCurrentIndex(indexComboBox)

        else:
            comboBox.addItem(self.selectionNotAvailable)
            comboBox.setCurrentIndex(-1)
            warnings.warn('Default %s for plot number %d may be used but is intentionally not made available in the '
                          'drop-down list.' % (plotAttribute, plotIndex))

        comboBox.resize(int(xSize), int(self.ySize))
        # noinspection PyUnresolvedReferences
        comboBox.currentIndexChanged.connect(connectMethod)

        return comboBox

    def addLineColorDropDown(self, objectName, newColor, availableColors):
        """
        Method to add drop-down list for plot line colors.

        In addition to the color names, the corresponding backgrounds are also changed to the corresponding colors.

        :param objectName: Drop-down list object name (string).
        :param newColor: New selected color for drop-down list (string).
        :param availableColors: List of available colors.
        """

        # noinspection PyArgumentList
        comboBox = QComboBox(self, objectName=objectName)
        for colorIndex in range(len(availableColors)):
            comboBox.addItem(availableColors[colorIndex])

            # selection background color: requires self.applicationObject.setStyle("fusion")
            comboBox.model().item(colorIndex).setBackground(QColor(availableColors[colorIndex]))

        indexComboBox = comboBox.findText(self.colorDictionary[newColor], Qt.MatchFixedString)
        if indexComboBox >= 0:
            comboBox.setCurrentIndex(indexComboBox)

            colors = comboBox.model().item(indexComboBox).background().color()
            palette = comboBox.palette()
            palette.setColor(QPalette.Button, colors)
            comboBox.setPalette(palette)

        comboBox.resize(self.xSize, self.ySize)
        # noinspection PyUnresolvedReferences
        comboBox.currentIndexChanged.connect(self.changeLineColor)

        return comboBox

    def toggleAllPlotVisibility(self):
        """
        Method to toggle the visibility of all plots (either all visible or all invisible)

        When all plots are made invisible the warning 'No handles with labels found to put in legend.'
        is written to the terminal.
        """

        checkBox = self.sender()
        checkBoxName = checkBox.objectName()
        checkBoxState = Qt.Checked if checkBox.isChecked() else Qt.Unchecked
        for widget in self.widgetList:
            if isinstance(widget, QCheckBox) and widget.objectName() != checkBoxName:
                widget.setCheckState(checkBoxState)

    def toggleSinglePlotVisibility(self):
        """
        Method to toggle the visibility of individual plots.
        """

        checkBox = self.sender()
        checkBoxName = checkBox.objectName()
        plotIndex = int(checkBoxName[:3])
        if checkBox.isChecked():
            self.axisHandle.lines[plotIndex].set_visible(True)
            self.axisHandle.lines[plotIndex].set_label(checkBoxName[3:])
        else:
            self.axisHandle.lines[plotIndex].set_visible(False)
            self.axisHandle.lines[plotIndex].set_label('_nolegend_')

        self.redrawPlot()

    def changePlotLabel(self):
        """
        Method to update the label of an individual plot after triggering of a UI event.

        This method is typically called after ENTER is pressed in the corresponding textbox.
        """

        textBox = self.sender()
        textBoxContent = textBox.text()
        textBoxName = textBox.objectName()
        newTextBoxName = textBoxName[:3] + textBoxContent
        for widget in self.widgetList:
            if widget.objectName() == textBoxName:
                widget.setObjectName(newTextBoxName)

        plotIndex = int(textBoxName[:3])
        if self.axisHandle.lines[plotIndex].get_label() != '_nolegend_':
            self.axisHandle.lines[plotIndex].set_label(textBoxContent)
            self.redrawPlot()

    def changeLineColor(self):
        """
        Method to update the line color of an individual plot after triggering of a UI event.

        This method is called after a new color is selected on the corresponding drop-down list.
        """

        comboBox = self.sender()
        objectName = comboBox.objectName()
        plotIndex = int(objectName[:3])
        colorName = comboBox.currentText()
        colorHex = [x[0] for x in self.colorDictionary.items() if x[1] == colorName][0]
        self.axisHandle.lines[plotIndex].set_color(colorHex)
        self.redrawPlot()

        # selection background color: requires self.applicationObject.setStyle("fusion")
        selectedIndex = comboBox.currentIndex()
        colors = comboBox.model().item(selectedIndex).background().color()
        palette = comboBox.palette()
        palette.setColor(QPalette.Button, colors)
        comboBox.setPalette(palette)

    def changeLineStyle(self):
        """
        Method to update the line style of an individual plot after triggering of a UI event.

        This method is called after a new line style is selected on the corresponding drop-down list.
        """

        comboBox = self.sender()
        objectName = comboBox.objectName()
        plotIndex = int(objectName[:3])
        lineStyle = [styleEntry[0] for styleEntry in self.lineStyleDictionary.items()
                     if styleEntry[1] == comboBox.currentText()][0]
        self.axisHandle.lines[plotIndex].set_linestyle(lineStyle)
        self.redrawPlot()

    def changeLineWidth(self):
        """
        Method to update the line width of an individual plot after triggering of a UI event.

        This method is typically called after ENTER is pressed in the corresponding textbox.
        """

        textBox = self.sender()
        objectName = textBox.objectName()
        plotIndex = int(objectName[:3])
        self.axisHandle.lines[plotIndex].set_linewidth(textBox.text())
        self.redrawPlot()

    def changeMarkerStyle(self):
        """
        Method to update the marker style of an individual plot after triggering of a UI event.

        This method is called after a new marker style is selected on the corresponding drop-down list.
        """

        comboBox = self.sender()
        objectName = comboBox.objectName()
        plotIndex = int(objectName[:3])
        markerStyle = [styleEntry[0] for styleEntry in self.markerStyleDictionary.items()
                       if styleEntry[1] == comboBox.currentText()][0]
        self.axisHandle.lines[plotIndex].set_marker(markerStyle)
        self.redrawPlot()

    def changeMarkerSize(self):
        """
        Method to update the marker size of an individual plot after triggering of a UI event.

        This method is typically called after ENTER is pressed in the corresponding textbox.
        """

        textBox = self.sender()
        objectName = textBox.objectName()
        plotIndex = int(objectName[:3])
        self.axisHandle.lines[plotIndex].set_markersize(textBox.text())
        self.redrawPlot()
