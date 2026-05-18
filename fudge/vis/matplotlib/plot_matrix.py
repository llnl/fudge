# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy
from matplotlib import pyplot
from matplotlib import colors
from matplotlib import ticker

def plot_matrix( matrix, energyBoundariesX=None, energyBoundariesY=None, title="Matrix", xyTitle='Energy (MeV)',
        xylog=False, zlog=False, zRange=(), switchY=False, colorMap='jet', subplot=None ):
    """Plot a matrix. Arguments are:
    :matrix (required): 2-dimensional numpy array (or equivalent)
    :energyBoundariesX: list of group boundary edges. If omitted, use matrix indices for axis
    :energyBoundariesY: "". If omitted, use same boundaries as x-axis
    :title: specify plot title
    :xyTitle: specify x (and optionally y) titles. May be string or tuple of two strings
    :xylog: if True, energyBoundariesX (and Y if different) must also be specified
    :zlog: if True, use log scale for z-axis colorbar. See note below for symmetric log scale.
    :zRange: for plotting a sub-range of z
    :switchY: if True, plot matrix[0,0] at upper left rather than lower left
    :colorMap: select color range for z-axis. Options include 'jet', 'jet_r', 'correlation', etc.
    :subplot: optional matplotlib.axes.AxesSubplot instance (for drawing matrix inside larger figure)

    This prepares a 2-d matrix plot. Use pyplot.show() to view the result.

    plot_matrix supports symmetric log-scale for the z-axis, using these arguments:
        zlog=True, zRange=(-z,z), colorMap=plot_matrix.generateColorMap()"""

    matrix = numpy.array( matrix )
    if matrix.ndim != 2:
        raise Exception("plot_matrix can only handle 2-d arrays!")

    if not subplot:
        fig = pyplot.figure()
        subplot = fig.add_subplot(111)

    if type(xyTitle) is str: xtitle = ytitle = xyTitle
    elif len(xyTitle)==2: xtitle, ytitle = xyTitle
    else:
        raise Exception("Can't understand xyTitle argument: %s" % xyTitle)
    pyplot.xlabel(xtitle)
    pyplot.ylabel(ytitle)

    if energyBoundariesY is None: energyBoundariesY = energyBoundariesX

    if energyBoundariesX is not None:
        lowX, highX = energyBoundariesX[0], energyBoundariesX[-1]
        lowY, highY = energyBoundariesY[0], energyBoundariesY[-1]
    else:
        lowX = lowY = 0
        highX, highY = matrix.shape[::-1]

    pyplot.xlim( lowX, highX )
    if switchY: # switch y-axis: put low energy at top
        pyplot.ylim( highY, lowY )
    else:
        pyplot.ylim( lowY, highY )

    if zRange: vmin, vmax = zRange
    else:
        vmin, vmax = matrix.min(), matrix.max()
        if zlog:
            try: vmin = min( matrix[matrix>0] )
            except ValueError:
                print("WARNING: no matrix elements > 0. Switching to linear scale for z-axis")
                zlog=False

    colorBarOpts = {}
    if zlog:    # z-axis log scale
        if vmin < 0:
            zopts = {'norm': colors.SymLogNorm( linthresh=vmax*1e-9, vmin=vmin, vmax=vmax )}

            ticks = []  # FIXME future versions of matplotlib may do this automatically for SymLogNorm
            zmax = 10**numpy.floor( numpy.log10( vmax ) )
            for i in range(4): ticks.append( -zmax * 10**(-3*i) )
            ticks.append(0)
            for i in range(3,-1,-1): ticks.append( zmax * 10**(-3*i) )
            colorBarOpts['ticks'] = ticks
            colorBarOpts['format'] = ticker.LogFormatterMathtext()
        else:
            zopts = {'norm': colors.LogNorm( vmin=vmin, vmax=vmax ) }
    else: zopts = {'vmin':vmin, 'vmax':vmax }

    # color map:
    zopts['cmap'] = colorMap
    if colorMap=='correlation': zopts['cmap'] = generateColorMap()

    if xylog:
        if energyBoundariesX is None: raise Exception("xylog option requires energy group boundary info")
        X,Y = numpy.meshgrid(energyBoundariesX,energyBoundariesY)
        pyplot.pcolormesh( X,Y, matrix, **zopts )
        pyplot.setp( subplot, xscale='log', yscale='log' )
    elif energyBoundariesX is not None:
        pyplot.pcolormesh( energyBoundariesX, energyBoundariesY, matrix, **zopts )
    else:
        pyplot.imshow( matrix, interpolation='none', **zopts )

    pyplot.title(title)
    pyplot.colorbar( **colorBarOpts )
    pyplot.matplotlib.cm.jet.set_under('w')

def generateColorMap():
    # based on jet_r, but the center is white (for 0 correlation)
    import matplotlib
    cdict = {'blue': [(0.0, 0, 0),
            (0.34999999999999998, 0, 0),
            (0.4995, 0.4828709677, 1),
            (0.5005, 1, 0.4828709677),
            (0.65999999999999992, 1, 1),
            (0.89000000000000001, 1, 1),
            (1.0, 0.5, 0.5)],
        'green': [(0.0, 0, 0),
            (0.089999999999999969, 0, 0),
            (0.35999999999999999, 1, 1),
            (0.4995, 1, 1),
            (0.5005, 1, 1),
            (0.625, 1, 1),
            (0.875, 0, 0),
            (1.0, 0, 0)],
        'red': [(0.0, 0.5, 0.5),
            (0.10999999999999999, 1, 1),
            (0.33999999999999997, 1, 1),
            (0.4995, 0.516129, 1),
            (0.5005, 1, 0.516129),
            (0.65000000000000002, 0, 0),
            (1.0, 0, 0)]} 
    return matplotlib.colors.LinearSegmentedColormap('correlation',cdict,1024)

