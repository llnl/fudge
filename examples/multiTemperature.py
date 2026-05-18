# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

""" Basic plotting example: reconstructs resonances (if necessary), and plots cross sections for requested MTs.
 The original file can be either in GNDS or ENDF format. This example also requires matplotlib

 Sample use:
>python plotCrossSection.py <filename> 102  293 600 1200 --temperatureUnit K   # plot capture at three temperatures """

import os
import sys
import argparse

from fudge import GNDS_file
from fudge import styles as stylesModule
from fudge import physicalQuantity as physicalQuantityModule

from fudge.vis.matplotlib import plot2d

from pqu import PQU

exampleDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, exampleDir)

parser = argparse.ArgumentParser(description="""
This example plots a cross section heated to one or more temperatures. The 'filename' argument can be either
a GNDS or ENDF-formatted file """)
parser.add_argument('filename', type=str, help='ENDF or GNDS file with cross section data')
parser.add_argument('mt', type=int, help='MT number to plot')
parser.add_argument('temps', type=float, nargs="+", help="Temperature(s) for heating.")
parser.add_argument("--temperatureUnit", default="K", help="Unit for temperatures (default='K')")
parser.add_argument("--energyUnit", help="Incident energy unit, e.g. 'eV' or 'MeV'")
parser.add_argument("--title", help="Plot title")
args = parser.parse_args()

filename = args.filename
if open(filename).readline().startswith("<?xml"):
    RS = GNDS_file.read(filename)
else:
    from brownies.legacy.converting import endfFileToGNDS
    rce = endfFileToGNDS.endfFileToGNDS(filename)
    RS, CS = rce['reactionSuite'], rce['covarianceSuite']

xUnit = RS.domainUnit
if args.energyUnit:
    RS.convertUnits({RS.domainUnit: args.energyUnit})
    xUnit = args.energyUnit

reaction = RS.getReaction(args.mt)

data = {}
for temp in args.temps:

    temp = PQU.PQU( temp, args.temperatureUnit )
    print("Heating to %s" % temp)
    heatedStyle = stylesModule.Heated(
            'heated', derivedFrom=RS.styles.getEvaluatedStyle().label,
            temperature=physicalQuantityModule.Temperature(temp.value, temp.unit))
    data[temp] = reaction.crossSection.heat( heatedStyle, EMin=PQU.PQU(1e-5, 'eV'))


# plotting:
lows, highs = zip(*[d.domain() for d in data.values()])
low, high = min(lows), max(highs)
title = "%s MT=%s" % (os.path.split(filename)[-1], args.mt)
if args.title is not None:
    title = args.title

datasets = []
for key in sorted( data.keys() ):   # sort by temperature
    datasets.append( plot2d.DataSet2d( data[key], legend=key ) )

plot2d.makePlot2d(
        datasets,
        xAxisSettings=plot2d.AxisSettings(label=f"$E_n$ ({xUnit})", isLog=True, axisMin=low, axisMax=high),
        yAxisSettings=plot2d.AxisSettings(label="Cross Section (barn)", isLog=True, autoscale=True),
        title=title, legendOn=True)
