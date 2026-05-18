/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtX = "%18.11e\n", *fmtXY = "%18.11e %18.11e\n";

static void compareResults( statusMessageReporting *smr, ptwXPoints *groupedData, ptwXPoints *groupedData4, char const *label );
static ptwXPoints *getGroupBoundaries( statusMessageReporting *smr );
static ptwXYPoints *getFluxData( statusMessageReporting *smr );
static ptwXYPoints *getCrossSectionData( statusMessageReporting *smr );
static ptwXYPoints *getMultiplicityData( statusMessageReporting *smr );
static ptwXYPoints *getConstantOneData( statusMessageReporting *smr );
static void writeXYDataOnVerbosity( ptwXYPoints *data, const char * const fileName );
static void writeXDataOnVerbosity( statusMessageReporting *smr, ptwXPoints *data, const char * const fileName );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    ptwXYPoints *flux, *crossSection, *multiplicity, *constantOne;
    ptwXPoints *groupBoundaries, *groupedData3, *groupedData4, *fluxGrouped_None;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    groupBoundaries = getGroupBoundaries( &smr );
    writeXDataOnVerbosity( &smr, groupBoundaries, "groupBoundaries.dat" );

    flux = getFluxData( &smr );
    writeXYDataOnVerbosity( flux, "flux.dat" );

    crossSection = getCrossSectionData( &smr );
    writeXYDataOnVerbosity( crossSection, "crossSection.dat" );

    multiplicity = getMultiplicityData( &smr );
    writeXYDataOnVerbosity( multiplicity, "multiplicity.dat" );

    constantOne = getConstantOneData( &smr );
    writeXYDataOnVerbosity( constantOne, "constantOne.dat" );

    if( ( fluxGrouped_None = ptwXY_groupOneFunction( &smr, flux, groupBoundaries, ptwXY_group_normType_none, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

/*
*   ptwXY_groupFourFunctions testing.
*/
    if( ( groupedData3 = ptwXY_groupThreeFunctions( &smr, crossSection, multiplicity, flux,             groupBoundaries, ptwXY_group_normType_norm,
            fluxGrouped_None ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    writeXDataOnVerbosity( &smr, groupedData3, "multiplicity_norm.dat" );

    if( ( groupedData4 = ptwXY_groupFourFunctions( &smr, constantOne, crossSection, multiplicity, flux,        groupBoundaries, 
            ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    compareResults( &smr, groupedData3, groupedData4, "norm-1" );

    if( ( groupedData4 = ptwXY_groupFourFunctions( &smr, crossSection, constantOne, multiplicity, flux,        groupBoundaries, 
            ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    compareResults( &smr, groupedData3, groupedData4, "norm-2" );

    if( ( groupedData4 = ptwXY_groupFourFunctions( &smr, crossSection, multiplicity, constantOne, flux,        groupBoundaries, 
            ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    compareResults( &smr, groupedData3, groupedData4, "norm-3" );

    if( ( groupedData4 = ptwXY_groupFourFunctions( &smr, crossSection, multiplicity, flux,        constantOne, groupBoundaries, 
            ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    compareResults( &smr, groupedData3, groupedData4, "norm-4" );

    ptwX_free( groupBoundaries );
    ptwXY_free( flux );
    ptwX_free( fluxGrouped_None );
    ptwXY_free( crossSection );
    ptwXY_free( multiplicity );

    exit( errCount );
}

/*
************************************************************
*/
static void compareResults( statusMessageReporting *smr, ptwXPoints *groupedData3, ptwXPoints *groupedData4, char const *label ) {

    int i1, n = (int) ptwX_length( smr, groupedData3 );
    double x1, x2, s, d, r;

    for( i1 = 0; i1 < n; ++i1 ) {
        x1 = groupedData3->points[i1];
        x2 = groupedData4->points[i1];
        s = fabs( x1 ) + fabs( x2 );
        d = x1 - x2;
        r = d;
        if( s != 0 ) r /= s;
        if( fabs( r ) > 1e-12 )
            printf( "ERROR %s: compare for '%s' at group %4d, %e %e %e %e %e\n", __FILE__, label, i1, x1, x2, s, d, r );
    }

    ptwX_free( groupedData4 );
}

/*
************************************************************
*/
static ptwXPoints *getGroupBoundaries( statusMessageReporting *smr ) {

    double data[] = {
        1.30680e-09, 2.09080e-08, 1.30680e-07, 3.34530e-07, 1.17610e-06, 2.09080e-06, 5.65780e-06, 1.30680e-05, 2.07460e-05, 5.12300e-05,
        1.02450e-04, 2.09080e-04, 3.81050e-04, 5.65780e-04, 7.15580e-04, 1.05850e-03, 1.30680e-03, 1.88170e-03, 2.94020e-03, 3.34530e-03,
        4.23390e-03, 5.76280e-03, 7.52700e-03, 1.02450e-02, 1.51060e-02, 2.09080e-02, 2.64620e-02, 3.26690e-02, 3.95300e-02, 7.00200e-02,
        9.89090e-02, 1.30680e-01, 1.81950e-01, 2.07460e-01, 2.41700e-01, 2.70970e-01, 2.94020e-01, 3.34530e-01, 3.77650e-01, 5.12300e-01,
        6.32470e-01, 7.52700e-01, 8.83370e-01, 1.02450e+00, 1.17610e+00, 1.33810e+00, 1.51060e+00, 1.69360e+00, 2.09080e+00, 2.30510e+00,
        2.52990e+00, 2.74110e+00, 3.01080e+00, 3.53350e+00, 4.06880e+00, 4.39600e+00, 4.70440e+00, 4.99080e+00, 5.35250e+00, 5.65780e+00,
        6.04250e+00, 6.36660e+00, 6.73670e+00, 7.15580e+00, 7.54790e+00, 7.90960e+00, 8.32150e+00, 8.78670e+00, 9.17670e+00, 9.66480e+00,
        1.01200e+01, 1.05850e+01, 1.10120e+01, 1.15470e+01, 1.19930e+01, 1.24990e+01, 1.30680e+01, 1.35420e+01, 1.38630e+01, 1.41340e+01,
        1.44070e+01, 1.46830e+01, 1.51860e+01, 1.57540e+01, 1.63340e+01, 1.69230e+01, 1.81340e+01, 2.00000e+01 };
    int nData = sizeof( data ) / sizeof( double );
    ptwXPoints *groupBoundaries;

    if( ( groupBoundaries = ptwX_create( smr, nData, nData, data ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( groupBoundaries );
}

/*
************************************************************
*/
static ptwXYPoints *getFluxData( statusMessageReporting *smr ) {

    double data[] = {  
        0.0000000e+00, 5.4598150e+01, 5.9634465e+00, 5.6305711e+00, 8.2124320e+00, 2.3903851e+00, 9.6171596e+00, 1.3997871e+00, 1.0649959e+01, 9.4447384e-01,
        1.1453293e+01, 6.9547558e-01, 1.2120111e+01, 5.3946075e-01, 1.2689944e+01, 4.3419372e-01, 1.3187622e+01, 3.5920720e-01, 1.3620180e+01, 3.0463531e-01,
        1.4009880e+01, 2.6260684e-01, 1.4364428e+01, 2.2942870e-01, 1.4689668e+01, 2.0269278e-01, 1.4989942e+01, 1.8078369e-01, 1.5268929e+01, 1.6255565e-01,
        1.5529441e+01, 1.4719790e-01, 1.5773780e+01, 1.3411483e-01, 1.5998524e+01, 1.2311016e-01, 1.6211129e+01, 1.1353226e-01, 1.6412837e+01, 1.0513506e-01,
        1.6604714e+01, 9.7724282e-02, 1.6787655e+01, 9.1145600e-02, 1.6962472e+01, 8.5273285e-02, 1.7129855e+01, 8.0005574e-02, 1.7290411e+01, 7.5258733e-02,
        1.7444603e+01, 7.0965397e-02, 1.7592983e+01, 6.7065295e-02, 1.7735973e+01, 6.3509804e-02, 1.7873951e+01, 6.0257736e-02, 1.8007253e+01, 5.7274143e-02,
        1.8136190e+01, 5.4528889e-02, 1.8261037e+01, 5.1996154e-02, 1.8382047e+01, 4.9653589e-02, 1.8607903e+01, 4.5560014e-02, 1.8821502e+01, 4.1999560e-02,
        1.9024105e+01, 3.8879888e-02, 1.9216791e+01, 3.6128174e-02, 1.9400467e+01, 3.3686629e-02, 1.9575956e+01, 3.1508211e-02, 1.9743954e+01, 2.9554878e-02,
        1.9905077e+01, 2.7795351e-02, 2.0059790e+01, 2.6204483e-02, 2.0208653e+01, 2.4759787e-02, 2.0352091e+01, 2.3443133e-02, 2.0490488e+01, 2.2239172e-02,
        2.0624179e+01, 2.1134888e-02, 2.0753480e+01, 2.0119061e-02, 2.0878669e+01, 1.9182082e-02, 2.1000000e+01, 1.8315639e-02 };
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *flux;

    if( ( flux = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( flux );
}

/*
************************************************************
*/
static ptwXYPoints *getCrossSectionData( statusMessageReporting *smr ) {

#include "Data/crossSection.dat.h"     /* This lines include data. */
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *crossSection;

    if( ( crossSection = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( crossSection);
}

/*
************************************************************
*/
static ptwXYPoints *getMultiplicityData( statusMessageReporting *smr ) {

    double data[] = {  
            1.00000000e-11,  3.95456000e+00, 2.00000000e-02,  3.95456000e+00, 3.00000000e-02,  3.95403000e+00, 4.00000000e-02,  3.95466000e+00,
            5.00000000e-02,  3.95591000e+00, 6.00000000e-02,  3.95750000e+00, 8.00000000e-02,  3.96126000e+00, 9.00000000e-02,  3.96508000e+00,
            1.00000000e-01,  3.96818000e+00, 1.20000000e-01,  3.97448000e+00, 1.50000000e-01,  3.98320000e+00, 2.00000000e-01,  3.99771000e+00,
            3.00000000e-01,  4.02781000e+00, 4.00000000e-01,  4.05829000e+00, 5.00000000e-01,  4.08896000e+00, 6.00000000e-01,  4.11958000e+00,
            7.00000000e-01,  4.15106000e+00, 8.00000000e-01,  4.17956000e+00, 9.00000000e-01,  4.20731000e+00, 1.00000000e+00,  4.23581000e+00,
            1.10000000e+00,  4.26337000e+00, 1.20000000e+00,  4.28905000e+00, 1.30000000e+00,  4.31651000e+00, 1.40000000e+00,  4.34344000e+00,
            1.50000000e+00,  4.36066000e+00, 1.60000000e+00,  4.38702000e+00, 1.70000000e+00,  4.41008000e+00, 1.80000000e+00,  4.43265000e+00,
            1.90000000e+00,  4.45567000e+00, 2.00000000e+00,  4.47668000e+00, 2.50000000e+00,  4.58012000e+00, 3.00000000e+00,  4.67677000e+00,
            3.50000000e+00,  4.76295000e+00, 4.00000000e+00,  4.84127000e+00, 5.00000000e+00,  4.97391000e+00, 6.00000000e+00,  5.03102000e+00,
            7.00000000e+00,  4.96053000e+00, 8.00000000e+00,  4.75777000e+00, 9.00000000e+00,  4.50265000e+00, 1.00000000e+01,  4.32857000e+00,
            1.10000000e+01,  4.34556000e+00, 1.20000000e+01,  4.45396000e+00, 1.30000000e+01,  4.59275000e+00, 1.40000000e+01,  4.71241000e+00,
            1.45000000e+01,  4.76865000e+00, 1.50000000e+01,  4.82063000e+00, 1.60000000e+01,  4.89104000e+00, 1.70000000e+01,  4.92432000e+00,
            1.80000000e+01,  4.90964000e+00, 1.90000000e+01,  4.87693000e+00, 2.00000000e+01,  4.86958000e+00 };
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *multiplicity;

    if( ( multiplicity = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( multiplicity );
}

/*
************************************************************
*/
static ptwXYPoints *getConstantOneData( statusMessageReporting *smr ) {

    double data[] = { 1.00000000e-11,  1.0, 2.00000000e+01, 1.0 };
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *constantOne;

    if( ( constantOne = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( constantOne );
}

/*
************************************************************
*/
static void writeXYDataOnVerbosity( ptwXYPoints *data, const char * const fileName ) {

    FILE *f;

    if( !verbose ) return;
    if( ( f = fopen( fileName, "w" ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: could not open file %s\n", __FILE__, fileName );
    fprintf( f, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, f, fmtXY );
    fclose( f );
}

/*
************************************************************
*/
static void writeXDataOnVerbosity( statusMessageReporting *smr, ptwXPoints *data, const char * const fileName ) {

    int64_t i;
    FILE *f;

    if( !verbose ) return;
    if( ( f = fopen( fileName, "w" ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: could not open file %s\n", __FILE__, fileName );
    fprintf( f, "# length = %d\n", (int) ptwX_length( smr, data ) );
    for( i = 0; i < ptwX_length( smr, data ); i++ ) fprintf( f, fmtX, ptwX_getPointAtIndex_Unsafely( data, i ) );
    fclose( f );
}
