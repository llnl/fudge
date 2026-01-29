# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Compare two Merced output files, keep track of maximum relative / absolute diffs
"""

import sys
import numpy

f1,f2 = sys.argv[1:]

print('  ',f1,'output differs from baseline: ', end='')

def parse( utfil ):
    xsc, tm1, tm2 = [], [], []

    with open(utfil) as fin:
        for line in fin:
            if line.startswith('outputLegendreOrder'): break

            if line.startswith('Cross section'):
                npoints = int(line.split()[-1])
                next(fin)  # interpolation line

                for idx in range(npoints):
                    xsc.append( list(map(float, next(fin).split())) )

        for i in range(4): next(fin)
        for line in fin:
            if line.startswith('EinBin'): continue
            if line.startswith('Integrals'): break
            tm1.append( list(map(float, line.split())) )

        for line in fin:
            if line.startswith('EinBin'): continue
            tm2.append( list(map(float, line.split())) )

    return numpy.array(xsc), numpy.array(tm1), numpy.array(tm2)

xs1, tm1_1, tm1_2 = parse(f1)
xs2, tm2_1, tm2_2 = parse(f2)

if xs1.shape != xs2.shape:
    print("cross section size mismatch!")
    sys.exit()
if tm1_1.shape != tm2_1.shape or tm1_2.shape != tm2_2.shape:
    print("transfer matrix size mismatch!")
    sys.exit()

rtol = 1e-5
atol = 1e-8 # use numpy defaults for now

xscs = numpy.allclose(xs1,xs2, rtol, atol)
tm1s = numpy.allclose(tm1_1, tm2_1, rtol, atol)
tm2s = numpy.allclose(tm1_2, tm2_2, rtol, atol)

def summarize_diffs( arr1, arr2 ):
    absdiff = numpy.abs(arr1-arr2)
    mask = absdiff > (atol + rtol * numpy.abs(arr2))
    vals = absdiff[mask]
    # FIXME not complete yet

status = ''
if not xscs: status += ' xsc diffs'
if not tm1s: status += ' tm1 diffs'
if not tm2s: status += ' tm2 diffs'
if not status: status = ' all within tolerance'

print(status)
