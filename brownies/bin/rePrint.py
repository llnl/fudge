#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

usage = """
>rePrint.py ENDF-formatted-file.endf [MT]

rePrint translates an ENDF file to the new GNDS format. If MT is given, only translate the specified MT.
"""

import sys, os, shutil

binDir = os.path.dirname( os.path.abspath( __file__ ) )

from LUPY import subprocessing
from LUPY import checksums as checksumsModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import externalFile as externalFileModule

import brownies.legacy.toENDF6.toENDF6     # this import adds 'toENDF6' methods to many GNDS classes
import brownies.legacy.toENDF6.endfFormats as endfFormatsModule

from brownies.legacy.converting import endfFileToGNDS
from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

def process_args( ) :       # see https://docs.python.org/2/howto/argparse.html

    from argparse import ArgumentParser

    parser = ArgumentParser( description=usage )
    parser.add_argument( 'inputFile', type=str,                                                         help='ENDF-6 file to translate to/from GNDS')
    parser.add_argument( 'MT', type=int, nargs="?", default=None )
    parser.add_argument( '-v', '--verbose', default = 0, action = "count",                              help = "enable verbose output")
    parser.add_argument( '--verboseWarnings', default = False, action = "store_true",                   help = "increase detail for warnings" )
    parser.add_argument( '-x', '--xmlOnly', default=False, action="store_true",                         help = "only translate ENDF to GNDS/Python and write xml files. Do not re-write endf file" )
    parser.add_argument( '-t', '--translateOnly', default=False, action="store_true",                   help = "only translate ENDF to GNDS/Python, do not write any files" )
    parser.add_argument( "-c", '--checkCovars', default = False, action = "store_true",                 help = "enable covariance checking" )
    parser.add_argument( "-s", "--skipBadData", default = False, action = "store_true",                 help = "Recover from format errors if possible" )
    parser.add_argument( "-o", "--outline", default = False, action = "store_true",                     help = "Set outline option to True for GNDS/XML output" )
    parser.add_argument( "--skipCovariances", default = False, action = "store_true",                   help = "Do not translate covariance data" )
    parser.add_argument( "--skipReconstruction", default = False, action = "store_true",                help = "Do not reconstruct resonances" )

    parser.add_argument( "--printBadNK14", default = False, action = "store_true",                      help = "If true print warning for NK mismatch between MF 12 and 14" )
    parser.add_argument( "--continuumSpectraFix", default = False, action = "store_true",               help = "Skip unnormalizeable continuum gamma distributions" )
    parser.add_argument( "--ignoreBadDate", default = False, action = "store_true",                     help = "If true ignore malformed date in MF=1 MT=")
    parser.add_argument( "--acceptBadMF10FissionZAP", default = False, action = "store_true",           help = "allow MF=10 MT=18 IZAP=0" )
    parser.add_argument( "--ignoreMF32Errors", default = False, action = "store_true",                  help = "Skip MF32 sections with internal consistency errors and continue processing" )
    parser.add_argument( "--output", default = 'test', dest='outputFile',                               help = "Prefix for resulting endf6 output file names" )
    parser.add_argument( "--pythonInterpreter", default = sys.executable, dest='pythonInterpreter',     help = "Python interpreter to use " )
    parser.add_argument( "--formatVersion", default = GNDS_formatVersionModule.default, choices = GNDS_formatVersionModule.allowed,
                                                                                                        help = "Specifies the format for the outputted GNDS file. " )
    parser.add_argument('--IDs', choices = ('familiar', 'nucleus', 'nuclide'), default='nuclide',       help='Choose between light charged particle naming conversion: "nuclide" (i.e., H1, H2, H3, He3 and He4), "nucleus" (i.e., h1, h2, h3, he3 and he4), or "familiar" (i.e., p, d, t, h and a).')
    parser.add_argument('--printMassHistory', action='store_true',                                      help='If present, information on the loading of masses from the ENDF-6 file is printed.')

    return parser.parse_args()

args = process_args()

specialNuclearParticleID = specialNuclearParticleIDPoPsModule.Mode(args.IDs)

for fname in ('.endf6.xml', '.endf6-covar.xml', '.endf6.orig.noLineNumbers', '.endf6.noLineNumbers'):
    if os.path.exists( args.outputFile + fname ):
        os.remove( args.outputFile + fname )

if( args.MT is None ) :
    shutil.copy2( args.inputFile, args.outputFile + '.endf6.orig2' )
else :
    header, MAT, MTDatas = endfFileToGNDSMisc.parseENDFByMT_MF(args.inputFile, stripMATMFMTCount = False)
    f = open( args.outputFile + '.endf6.orig2', 'w' )
    f.write( header + '\n' )
    MTs = [ 451, args.MT ]
    if( args.MT == 18 ) :
        MTsp = [ 452, 455, 456, 458 ]
        for MTp in MTsp :
            if( MTp in MTDatas ) : MTs.append( MTp )
    MFLines = {}
    for MTp in MTs :
        if( MTp not in MTDatas ) : continue
        MFs = sorted( MTDatas[MTp] )
        for MF in MFs : 
            if( MF != 1 ) : continue
            if( MF not in MFLines ) : MFLines[MF] = []
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormatsModule.endfSENDLine( MAT, MF ) ]
    MTs.sort( )
    for MTp in MTs :
        if( MTp not in MTDatas ) : continue
        MFs = sorted( MTDatas[MTp] )
        for MF in MFs : 
            if( MF == 1 ) : continue
            if( MF not in MFLines ) : MFLines[MF] = []
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormatsModule.endfSENDLine( MAT, MF ) ]
    MFs = sorted( MFLines.keys( ) )
    for MF in MFs :
        f.write( '\n'.join( MFLines[MF] ) + '\n' )
        f.write( endfFormatsModule.endfFENDLine( MAT ) + '\n' )
    f.write( endfFormatsModule.endfMENDLine( ) + '\n' )
    f.write( '%s' % endfFormatsModule.endfTENDLine( ) + '\n' )
    f.close( )

flags = { 'verbosity' : 0 }
subprocessing.executeCommand([args.pythonInterpreter, os.path.join(binDir, 'reForm.py'), args.outputFile + '.endf6.orig2', args.outputFile + '.endf6.orig'])
os.remove( args.outputFile + '.endf6.orig2' )
style = 'eval'
rce = endfFileToGNDS.endfFileToGNDS(args.inputFile, singleMTOnly=args.MT, toStdOut=args.verbose,
                                    skipBadData=args.skipBadData, doCovariances=not args.skipCovariances,
                                    verboseWarnings=args.verboseWarnings, verbose=args.verbose, formatVersion=args.formatVersion,
                                    printBadNK14=args.printBadNK14, continuumSpectraFix=args.continuumSpectraFix,
                                    ignoreBadDate=args.ignoreBadDate, acceptBadMF10FissionZAP=args.acceptBadMF10FissionZAP,
                                    ignoreMF32Errors=args.ignoreMF32Errors,
                                    reconstructResonances=not args.skipReconstruction, specialNuclearParticleID=specialNuclearParticleID,
                                    printMassHistory=args.printMassHistory)

errs = rce['errors']
if( 'reactionSuite' in rce ) :
    reactions, covariances = rce['reactionSuite'], rce['covarianceSuite']
    if( args.checkCovars ) :
        sys.stderr.write( ''.join(
            [ "%s: %s" % ( os.path.split( args.inputFile )[-1], warning ) for warning in covariances.check( ) ]
        ) )
    if( args.translateOnly ) : sys.exit( len( errs ) )

    RSfile = args.outputFile + ".endf6.xml"
    CSfile = args.outputFile + ".endf6-covar.xml"
    if( covariances is not None ) :
        covariances.externalFiles.add(externalFileModule.ExternalFile("reactions", RSfile))
        covariances.saveToFile( CSfile, outline = args.outline, formatVersion = args.formatVersion )

        sha1sum = checksumsModule.Sha1sum.from_file(CSfile)
        reactions.externalFiles.add(externalFileModule.ExternalFile("covariances", path=CSfile, checksum=sha1sum))

    reactions.saveToFile( RSfile, outline = args.outline, formatVersion = args.formatVersion )
    if( args.xmlOnly ) : sys.exit( len( errs ) )

    if( args.verbose ) : flags['verbosity'] = 31
    with open( args.outputFile + '.endf6', 'w' ) as fout:
        fout.write( reactions.toENDF6( style, flags, covarianceSuite = covariances ) )

elif( 'fissionFragmentData' in rce ) :
    fissionFragmentData = rce['fissionFragmentData']
    fissionFragmentData.saveToFile( args.outputFile + '.endf6.xml', formatVersion=args.formatVersion)
    sys.exit( 0 )

elif( 'PoPs' in rce ) :
    pops = rce['PoPs']
    pops.saveToFile( args.outputFile + '.endf6.xml', outline = args.outline, formatVersion = args.formatVersion )

    if( args.verbose ) : flags['verbosity'] = 31
    with open( args.outputFile + '.endf6', 'w' ) as fout:
        fout.write( pops.toENDF6( style, flags ) )

else :
    raise Exception( 'Unsupported return from endfFileToGNDS.endfFileToGNDS' )

subprocessing.executeCommand([args.pythonInterpreter, os.path.join(binDir, 'noLineNumbers.py'), args.outputFile + '.endf6', args.outputFile + '.endf6.noLineNumbers'])
subprocessing.executeCommand([args.pythonInterpreter, os.path.join(binDir, 'noLineNumbers.py'), args.outputFile + '.endf6.orig', args.outputFile + '.endf6.orig.noLineNumbers'])
if( True ) : subprocessing.executeCommand( [ 'rm', '-f', args.outputFile + '.endf6', args.outputFile + '.endf6.orig' ] )

sys.exit( len( errs ) )
