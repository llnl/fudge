# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/UnitTesting/thinning'

os.system( 'cd %s; ./thin -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls, hasLabel ) :

    i = 0
    for l in ls :
        if( l.strip( ) != '' ) : break
        i = i + 1
    ls = ls[i:]
    if( len( ls ) == 0 ) : return( None, None, None, None )
    label, accuracy = None, None
    if( hasLabel ) :
        label, ls = ls[0].strip( ), ls[1:]
        accuracy, ls = ls[0].strip( ), ls[1:]
        if( accuracy[:13] != '# accuracy = ' ) : raise Exception( '%s: line does not contain accuracy info: "%s"' % ( __file__, accuracy[:-1] ) )
        accuracy = float( accuracy.split( '=' )[1] )
    length, ls = ls[0], ls[1:]
    if( '# length = ' != length[:11] ) : raise Exception( '%s: line does not contain length info: "%s"' % ( __file__, length.strip( ) ) )
    length = int( length.split( '=' )[1] )
    data = [ list( map( float, ls[i].split( )[:2] ) ) for i in range( length ) ]
    return( ls[length:], label, accuracy, pointwiseXY_C.pointwiseXY_C( data, initialSize = 100, overflowSize = 10 ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def thin( label, accuracy, original, data ) :

    thin = original.thin( accuracy )
    if( len( data ) != len( thin ) ) : raise Exception( '%s: len( data ) = %d != len( thin ) = %d for label = "%s"' % \
        ( __file__, len( data ), len( thin ), label ) )
    for i, xy in enumerate( data ) :
        xc, yc = xy
        xp, yp = thin[i]
        compareValues( label, i, xc, xp )
        compareValues( label, i, yc, yp )

while( 1 ) :
    ls, label, accuracy, original = getData( ls, True )
    if( ls is None ) : break
    ls, label, dummy, data = getData( ls, False )
    if( ls is None ) : raise Exception( '%s: missing thinned data for label = "%s"' % ( __file__, label ) )
    thin( label, accuracy, original, data )
