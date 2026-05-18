# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import math
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../../Test/UnitTesting/toFunction'

os.system( 'cd %s; ./createFromFunction -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls ) :

    ls = utilities.skipBlankLines( ls )
    if( len( ls ) == 0 ) : return( None, None, None, None, None, None )
    label, accuracy = None, None
    label, ls = ls[0].strip( ), ls[1:]
    if( label != '# Xs' ) : raise Exception( '%s: missing line "# Xs": found "%s"' % ( __file__, label ) )
    ls, Xs = utilities.getXData( ls, '=' )
    ls, accuracy = utilities.getDoubleValue( 'accuracy', ls )
    ls, biSectionMax = utilities.getDoubleValue( 'biSectionMax', ls )
    ls, data = utilities.getXYData( ls )
    return( ls, label, accuracy, biSectionMax, Xs, pointwiseXY_C.pointwiseXY_C( data, accuracy = accuracy ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def f_x_Sin_xx( x, args ) :

       return( x * math.sin( x * x ) )

def createFromFunction( label, accuracy, biSectionMax, Xs, original ) :

    data = pointwiseXY_C.createFromFunction( Xs, f_x_Sin_xx, None, accuracy, biSectionMax, checkForRoots = 1 )
    if( len( data ) != len( original ) ) : raise Exception( '%s: len( data ) = %d != len( original ) = %d for label = "%s"' % \
        ( __file__, len( data ), len( original ), label ) )

while( 1 ) :
    ls, label, accuracy, biSectionMax, Xs, original = getData( ls )
    if( ls is None ) : break
    createFromFunction( label, accuracy, biSectionMax, Xs, original )
