#! /usr/bin/env python3
  
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib
import shutil

from xData import XYs1d as XYs1dModule

summaryDocString__xData = '''Reads 1d data from each file listed into an XYs1d instance and plots all using XYs1d.multiPlot.'''

description = '''
Reads 1d data from each file listed into an XYs1d instance and plots all on a single plot using XYs1d.multiPlot.
The data in each file must be in columns and currently only columns 1 and 2 of each file are plotted.
The '#' character is treated as a comment character and, for each line, it and all remaining character 
on the line are ignored. Also, blank lines, after the comment is removed, are ignored.
'''

parser = argparse.ArgumentParser(description=description)

parser.add_argument('files', type=pathlib.Path, nargs='*',                  help='The list of files whose data are plotted.')
parser.add_argument('--xMin', type=float, default=None,                     help='The x-axis minimum value.')
parser.add_argument('--xMax', type=float, default=None,                     help='The x-axis maximum value.')
parser.add_argument('--yMin', type=float, default=None,                     help='The y-axis minimum value.')
parser.add_argument('--yMax', type=float, default=None,                     help='The y-axis maximum value.')
parser.add_argument('--logs', type=int, default=0,                          help='Set the x (y)-axis to log is value is odd (value/2 is oadd).')
parser.add_argument('--saveToFile', type=pathlib.Path,                      help='Set the x (y)-axis to log is value is odd (value/2 is oadd).')

def readCurve(file):

    with open(file) as fIn:
        lines = fIn.readlines()

    data = []
    for line in lines:
        line = line.split('#')[0]
        values = list(map(float, line.split()))
        if len(values) > 0:
            data.append([values[0], values[1]])

    curve = XYs1dModule.XYs1d(data=data)
    curve.plotLabel = file

    return curve

if __name__ == '__main__':
    args = parser.parse_args()

    curves = []
    for file in args.files:
        try:
            curve = readCurve(file)
            if len(curve) == 0:
                continue
        except:
            print('Reading %s failed.' % file)
            continue
        curves.append(curve)

    if len(curves) > 0:
        XYs1dModule.XYs1d.multiPlot(curves, domainMin=args.xMin, domainMax=args.xMax, 
                rangeMin=args.yMin, rangeMax=args.yMax, logs=args.logs)
