#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from pqu import PQU as PQUModule

from fudge.processing import flux as fluxModule

energyUnitDefault = "MeV"

from argparse import ArgumentParser

summaryDocString__FUDGE = '''Adds a flux definition (label and f(T,E,mu) data) to a fluxes file (e.g., fluxes.xml).'''

description = """Read a list of energies and fluxes from a file or command-line and adds it to a flux file. For each energy, the fluxes must be stored as Legendre coefficients.

This script does not currently support creating a temperature-dependent flux."""

parser = ArgumentParser( description = description )
parser.add_argument( "label",                                           help = """The label to assign to the flux.""" )
parser.add_argument( "flux",                                            help = """A file containing a list of energies and fluxes for each energy. For each energy, the flux is represented as Legendre coefficients. Each energy and associated Legendre coefficient must be on one line.""" )
parser.add_argument( "output",                                          help = """Name of the outputted flux file. If the "input" argument is not present, then the output file, if present, is used for the input file. Otherwise, no input file is used.""" )
parser.add_argument( "input", nargs = "?", default = None,              help = """If present, the file to read existing flux data from. The new fluxes are appended to the fluxes in this file and all are written to the output file.""" )
parser.add_argument( "--override", action = "store_true",               help = """If label exists in flux file and option present, replace flux; otherwise, execute a raise.""" )
parser.add_argument( "-u", "--unit", default = energyUnitDefault,       help = """The unit for the energies. Default is %s.""" % energyUnitDefault )
parser.add_argument( "--values", action = "store_true",                 help = """If option is present, the "flux" argument should be a quoted string like "0.0 1.0; 200.0 1.0" rather than a file name.""" )

args = parser.parse_args( )

unit = PQUModule.PQU( 1, args.unit )
if not unit.isEnergy( ): raise TypeError( "Unit must be an energy unit." )

XYs2d = fluxModule.XYs2d( outerDomainValue = 0.0 )
if args.values:
    for xys1d in args.flux.split( ';' ) :
        energy_coefficients = list( map( float, xys1d.split( ) ) )
        XYs2d.append( fluxModule.LegendreSeries( energy_coefficients[1:], outerDomainValue = energy_coefficients[0] ) )
else :
    fIn = open( args.flux )
    lines = fIn.readlines( )
    fIn.close( )

    for line in lines :
        data = list( map( float, line.split( ) ) )
        energy = data.pop( 0 )
        XYs2d.append( fluxModule.LegendreSeries( data, outerDomainValue = energy ) )

assert len(XYs2d) > 1, "Flux definition must contain at least two incident energies"

flux = fluxModule.XYs3d( label = args.label, axes = fluxModule.axes( args.unit ) )
flux.append( XYs2d )

input = args.input
if input is None:
    if os.path.exists( args.output ): input = args.output

if input is not None:
    fluxes = fluxModule.Fluxes.readXML_file(input)
else :
    fluxes = fluxModule.Fluxes( )

if flux.label in fluxes:
    if not args.override: raise ValueError( """Label "%s" already in flux file.""" % flux.label )
    fluxes.replace( flux )
else :
    fluxes.add( flux )

fluxes.saveToFile( args.output )
