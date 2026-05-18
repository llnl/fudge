#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

parser = argparse.ArgumentParser( '' )
parser.add_argument( '-v', '--verbose', action = 'store_true', 
    help = 'print each file searched and, for each breakup reaction found, print MT and LR' )
parser.add_argument( 'paths', nargs='+', help = 'path for each ENDF to search for breakup reactions' )

args = parser.parse_args( )

filesWithBreakup = 0
for file in args.paths :
    header, MAT, MTDatas = endfFileToGNDSMisc.parseENDFByMT_MF(file)
    if( args.verbose ) : print(file)
    doPrintHeader, counter = True, 0
    for MT in sorted( MTDatas.keys( ) ) :
        MTData = MTDatas[MT]
        if( 3 in MTData ) :
            MF3 = MTData[3]
            dummy, dummy, dummy, LR, dummy, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MF3[1])
            if( LR != 0 ) :
                counter += 1
                if( args.verbose ) : 
                    if( doPrintHeader ) : print( '     MT  LR' )
                    doPrintHeader = False
                    print( '    %3d  %d' % ( MT, LR )  )
                elif( counter == 1 ) :
                    print( file )
    if( counter > 0 ) : filesWithBreakup += 1
if( not( args.verbose ) ) : print( '# %d files checked, %d with breakup reactions' % ( len( args.paths ), filesWithBreakup ) )
