# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities as utilitiesModule

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/UnitTesting/thinning'

os.system( 'cd %s; ./thinDomain -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.14e' % v1, '%.14e' % v2
    if( sv1 != sv2 ) :
        raise Exception( '%s: values %.17g and %.17g diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def checkThin( ls ) :

    ls = utilitiesModule.skipBlankLines( ls )
    if( ls == [] ) : return( ls )

    ls, label = utilitiesModule.getStringValue( 'label', ls )
    ls, epsilon = utilitiesModule.getDoubleValue( 'epsilon', ls )
    ls, original = utilitiesModule.getXYData( ls )
    ls, answer = utilitiesModule.getXYData( ls )
    ls, dummy = utilitiesModule.getXYData( ls )
    thinned = original.thinDomain( epsilon )
    if( len( answer ) != len( thinned ) ) : raise Exception( '%s: len( answer ) = %d != len( thinned ) = %d for label = "%s"' % \
        ( __file__, len( answer ), len( thinned ), label ) )
    for i, xy in enumerate( answer ) :
        xc, yc = xy
        xp, yp = thinned[i]
        compareValues( label, i, xc, xp )
        compareValues( label, i, yc, yp )
    return( ls )

while( 1 ) :
    if( ls == [] ) : break
    ls = checkThin( ls )
