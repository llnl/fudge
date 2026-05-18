#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# This script takes a gnds/XML file, reads it in and rewrites it to a file as an ENDF6 file
# (in the currently directory) with the extension '.endf' added.

import os
import argparse
import pathlib

from fudge import styles as stylesModule, GNDS_file as GNDS_fileModule

# add 'toENDF6' methods to GNDS classes:
import brownies.legacy.toENDF6.toENDF6

description = """Translate a GNDS file to ENDF-6.
Sample use: python gnds2endf.py n-001_H_001.xml n-001_H_001.endf
If n-001_H_001.xml links to a covarianceSuite (in the 'externalFiles' section), covariances will automatically be loaded.
The output file name is optional, defaults to the input file with '.endf' appended."""

__doc__ = description

parser = argparse.ArgumentParser(description)
parser.add_argument('-l', '--lineNumbers', action='store_true',             help='Add line numbers')
stylesGroup = parser.add_mutually_exclusive_group()
stylesGroup.add_argument('--style', default=None,                           help='Label of style to translate. If not supplied, choose the latest evaluated style')
stylesGroup.add_argument('--writeReconstructed', action='store_true',       help='If present, the resonance reconstructed cross section data are written, if present, instead of the background cross section data.')
parser.add_argument('-v', '--verbose', action='count', default=0,           help='Determines verbosity level, the more the merrier.')
parser.add_argument('--skipCovariances', action='store_true',               help='If present, any covariance files in are not written.')
parser.add_argument('--NLIB', type=int, default=-1,                         help='Specifies the NLIB value to set in the output ENDF-6 file.')
parser.add_argument('--useRedsFloatFormat', action='store_true',            help='''If present, add float are written using Red Cullen's floating point format.''')
parser.add_argument('gnds',                                                 help='GNDS file to translate.')
parser.add_argument('output', nargs='?', default=None,                      help='The translated ENDF file.')

if __name__ == '__main__':
    args = parser.parse_args( )

    gndsCov = None
    gnds = GNDS_fileModule.read(args.gnds)

    if hasattr(gnds, 'domainUnit') and gnds.domainUnit != 'eV':
        gnds.convertUnits({gnds.domainUnit: 'eV'})

    if not args.skipCovariances:
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()
            if len(covariances) == 1:
                gndsCov = covariances[0]
            elif len(covariances) > 1:
                raise NotImplementedError("Converting multiple covarianceSuites back to ENDF-6")

    if gndsCov is not None and gndsCov.domainUnit != 'eV':
        gndsCov.convertUnits({gndsCov.domainUnit: 'eV'})

    styleLabel = args.style
    if styleLabel is None:
        styleLabel = gnds.styles[0].label
        if hasattr(gnds.styles, 'preProcessingChains'):
            preProcessingChains = gnds.styles.preProcessingChains( ends = True )
            if( len( preProcessingChains ) != 1 ) : raise Exception( 'Currently, only one chain is supported.' )
            for style in reversed(preProcessingChains[0]) :
                if( isinstance( style, stylesModule.Evaluated ) ) : styleLabel = style.label
                if args.writeReconstructed:
                    if( isinstance( style, stylesModule.CrossSectionReconstructed ) ) : styleLabel = style.label
            if len(gnds.styles.getStylesOfClass(stylesModule.Evaluated)) > 1:
                print("INFO: translating style '%s' to ENDF-6" % styleLabel)

    output = args.output
    if output is None:
        output = pathlib.Path(args.gnds).with_suffix('.endf')
    with open(output, 'w') as fout:
        fout.write(gnds.toENDF6(styleLabel, flags={'verbosity': args.verbose * 10}, covarianceSuite=gndsCov, lineNumbers=args.lineNumbers, NLIB=args.NLIB,
                useRedsFloatFormat=args.useRedsFloatFormat))
