#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import traceback
import argparse

binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from fudge import externalFile as externalFileModule
from LUPY import checksums as checksumsModule

from brownies.legacy.converting import endfFileToGNDS

summaryDocString__Brownies = '''Converts an ENDF6 file into GNDS.'''

def process_args( ) :

    parser = argparse.ArgumentParser( description = "Translates an ENDF file to the GNDS format." )
    parser.add_argument( "inFile",                                                              help = "input file" )
    parser.add_argument( "outputFile", nargs = "?", default = None,                             help = "output file" )
    parser.add_argument( "-v", "--verbose", action = "count", default = 0,                      help = "enable verbose output. The more the gabbier." )
    parser.add_argument( "--skipBadData", action = "store_true", default = False,               help = "skip bad data, rather than throw an exception, when reading an ENDF file" )
    parser.add_argument( "--skipCovariances", action = "store_true", default = False,           help = "skip the covariance, if present" )
    parser.add_argument( "--verboseWarnings", action = "store_true", default = False,           help = "print verbose warnings" )
    parser.add_argument( "--printBadNK14", action = "store_true", default = False,              help = "print bad NK's if found" )
    parser.add_argument( "--continuumSpectraFix", action = "store_true", default = False,       help = "fix continuous spectra on read, if foobar" )
    parser.add_argument( "--ignoreBadDate", action = "store_true", default = False,             help = "ignore malformed ENDF dates" )
    parser.add_argument( "--acceptBadMF10FissionZAP", action = "store_true", default = False,   help = "allow MF=10 MT=18 IZAP=0" )
    parser.add_argument( "--traceback", action = "store_true", default = False,                 help = "print traceback on exception" )
    parser.add_argument( "--formatVersion", default = GNDS_formatVersionModule.default, choices = GNDS_formatVersionModule.allowed,
                                                                                                help = "Specifies the format for the outputted GNDS file. " )
    parser.add_argument('--energyUnit', type=str, default='eV',                                 help='Convert all energies in the gnds file to this unit.')
    parser.add_argument("--JENDL_stylePrimarygammas", action="store_true",                      help="If prsent, treat MF=6 primary gamma energies as binding energy, otherwise as gamma energy.")
    parser.add_argument('--IDs', choices = ('familiar', 'nucleus', 'nuclide'), default='nuclide',
            help='Choose between light charged particle naming conversion: "nuclide" (i.e., H1, H2, H3, He3 and He4), "nucleus" (i.e., h1, h2, h3, he3 and he4), or "familiar" (i.e., p, d, t, h and a).')
    return parser.parse_args()

args = process_args()

specialNuclearParticleID = specialNuclearParticleIDPoPsModule.Mode(args.IDs)

inFile = args.inFile
outFile = args.outputFile
if outFile is None:
    outFile = os.path.basename(inFile) + '.gnds.xml'
    outCovFile = os.path.basename(inFile) + '.gnds-covar.xml'
elif '.xml' in outFile:
    outCovFile = outFile.replace('.xml', '-covar.xml')
else:
    outCovFile = outFile + "-covar.xml"

try:
    kwargs = {}
    if args.JENDL_stylePrimarygammas:
        kwargs['JENDL_stylePrimarygammas'] = True
    results = endfFileToGNDS.endfFileToGNDS(inFile, toStdOut=args.verbose, verbose=args.verbose, skipBadData=args.skipBadData, 
            doCovariances=not args.skipCovariances, 
            verboseWarnings=args.verboseWarnings, printBadNK14=args.printBadNK14, continuumSpectraFix=args.continuumSpectraFix,
            ignoreBadDate=args.ignoreBadDate, acceptBadMF10FissionZAP=args.acceptBadMF10FissionZAP, formatVersion=args.formatVersion,
            specialNuclearParticleID=specialNuclearParticleID, **kwargs)

    reactionSuite = results.get('reactionSuite', None)
    covariance = results.get('covarianceSuite', None)
    fissionFragmentData = results.get('fissionFragmentData', None)
    pops = results.get('PoPs', None)
    errors = results['errors']
except Exception as err:
    sys.stderr.write('WARNING: read ENDF error: %s\n' % err)
    if args.traceback:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)
    exit(1)

if covariance:
    if args.energyUnit != covariance.domainUnit:
        covariance.convertUnits({covariance.domainUnit: args.energyUnit})

try:
    for topLevel in (fissionFragmentData, reactionSuite, pops):
        if topLevel is not None:
            if hasattr(topLevel, 'convertUnits') and args.energyUnit != 'eV':
                topLevel.convertUnits({'eV': args.energyUnit})
            if hasattr(topLevel, 'saveAllToFile'):
                topLevel.saveAllToFile(outFile, formatVersion=args.formatVersion)
            else:
                topLevel.saveToFile(outFile, formatVersion=args.formatVersion)
            break
except Exception as err:
    sys.stderr.write('WARNING: ENDF write error: %s\n' % err)
    if args.traceback:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback)

exit(len(errors))
