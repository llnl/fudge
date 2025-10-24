/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*  
    todo:

        -) Should check for nan or infty?
*/

#include <Python.h>
#include <pyport.h>
#include "structmember.h"
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#define POINTWISEXY_C_MODULE 1
#include <ptwXY.h>
#include <pointwiseXY_C.h>

#if PY_MAJOR_VERSION < 3
    #define PYUNICODE_CHECK PyString_Check
    #define PYUNICODE_ASSTRING PyString_AsString
    #define PYUNICODE_FROMSTRING PyString_FromString

    #define STATICFORWARD staticforward

#else
    #define PYUNICODE_CHECK PyUnicode_Check
    #define PYUNICODE_ASSTRING PyUnicode_AsUTF8
    #define PYUNICODE_FROMSTRING PyUnicode_FromString

    #define STATICFORWARD static
    #define Py_TPFLAGS_CHECKTYPES 0
#endif

#define pyObject_Unsupported 0
#define pyObject_Number 1
#define pyObject_List 2
#define pyObject_Tuple 3
#define pyObject_ptwXY 4

enum e_interpolationType { e_interpolationTypeInvalid = -1, e_interpolationTypeLinear, e_interpolationTypeLog, e_interpolationTypeFlat, e_interpolationTypeOther };

typedef nfu_status (*ptwXY_ptwXY_d_func)( statusMessageReporting *smr, ptwXYPoints *, double );
typedef ptwXYPoints *(*ptwXY_ptwXY_ptwXY_func)( statusMessageReporting *smr, ptwXYPoints *, ptwXYPoints * );

static double _defaultAccuracy = 1e-3;

STATICFORWARD PyTypeObject pointwiseXY_CPyType;

typedef struct pointwiseXY_CPy_s {
    PyObject_HEAD
    statusMessageReporting smr;
    ptwXYPoints *ptwXY;
    int infill;
    int safeDivide;
} pointwiseXY_CPy;

#define is_pointwiseXY_CPyObject( v ) ( Py_TYPE( v ) == &pointwiseXY_CPyType )

typedef struct GnG_parameters_s {
    PyObject *func;
    PyObject *argList;
} GnG_parameters;

static char pointwiseXY_C__doc__[] = 
    "The pointwiseXY_C class stores and manipulates a list of XY points (i.e., [x, y] pairs).\n" \
    "Methods to add, substract, multiply and divide a pointwiseXY_C object with a scaler\n" \
    "(i.e., a number) or other pointwiseXY_C object are provided.\n" \
    "\n" \
    "Constructor:\n" \
    "pointwiseXY_C( data = [], dataForm = 'xys', initialSize = 100, overflowSize = 10, accuracy = defaultAccuracy( ), biSectionMax = 3." \
    ", interpolation = 'lin-lin', infill = True, safeDivide = False, userFlag = 0 )\n\n" \
    "Constructor arguments are ([o] implies optional argument):\n" \
    "   data            [o] the [x_i, y_i] pairs given as described by the dataForm argument,\n" \
    "   dataForm        [o] can be one of three strings (case is ignored) that describes the form of data:\n" \
    "                          'XYs'      data are a sequence of [x_i, y_i] pairs with x_i ascending (i.e., x_i < x_{i+1})\n" \
    "                                          (e.g., data = [ [ 1, 1 ], [ 2.3, 2 ], [ 3.4, 6 ], [ 5.1, 4.3 ] ])\n" \
    "                          'XsAndYs'  data are given as [xs, ys] where xs the list of x-values and ys the matching list of y-values,\n" \
    "                                          (e.g., data = [ [ 1, 2.3, 3.4, 5.1 ], [ 1, 2, 6, 4.3 ] ]),\n" \
    "                          'list'     data are given as [x_0, y_0, x_1, y_1, ..., x_n, y_n],\n" \
    "                                          (e.g., data = [ 1, 1, 2.3, 2, 3.4, 6, 5.1, 4.3 ]),\n" \
    "   initialSize     [o] the initial size of the primary data cache (default = 100),\n" \
    "   overflowSize    [o] the initial size of the secondary (overflow) data cache (default = 10),\n" \
    "   accuracy        [o] the accuracy of the data for the given interpolation (default is value returned by function defaultAccuracy)\n" \
    "   biSectionMax    [o] at times (e.g., multiplication), points may need to be added to maintain the give accuracy. In this case, a region\n" \
    "                       is continuously divided into two until the accuracy is met or biSectionMax divisions have occurred (default = 3),\n" \
    "   interpolation   [o] can be one of the following strings:\n" \
    "                          'flat'              for the domain [x_i,x_{i+1}), the y-value is y_i\n" \
    "                          'lin-lin'\n" \
    "                          'lin-log'\n" \
    "                          'log-lin'\n" \
    "                          'log-log'\n" \
    "                          'other'\n" \
    "   infill          [o] if True, multiplication will continuously divide a region until accuracy or biSectionMax is met (default = True),\n" \
    "   safeDivide      [o] if True, safe division is used (default = False)\n" \
    "   userFlag        [o] an integer (of type C int) used to store a user defined flag. (see getUserFlag and setUserFlag).\n";

static pointwiseXY_CPy *pointwiseXY_CNewInitialize( int infill, int safeDivide );
static int pointwiseXY_C__init__( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static void pointwiseXY_C_dealloc( PyObject *self );
static PyObject *pointwiseXY_C__repr__( pointwiseXY_CPy *self );
static Py_ssize_t pointwiseXY_C__len__( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C__getitem__( pointwiseXY_CPy *self, Py_ssize_t index );
static PyObject *pointwiseXY_C__getslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_ );
static PyObject *pointwiseXY_C_getslice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static int pointwiseXY_C__setitem__( pointwiseXY_CPy *self, Py_ssize_t index, PyObject *value );
static int pointwiseXY_C__setslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *value );
static PyObject *pointwiseXY_C__add__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__iadd__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__sub__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__isub__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__mul__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__imul__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__div__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C__idiv__( PyObject *self, PyObject *other );
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( statusMessageReporting *smr,
        ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( statusMessageReporting *smr,
        ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
static PyObject *pointwiseXY_C__mod__( PyObject *self, PyObject *other );
static PyObject *pointwiseXY_C_add_sub_mul_div_number( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func );
static PyObject *pointwiseXY_C_add_sub_mul_div_number_insitu( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func );
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( pointwiseXY_CPy *self, pointwiseXY_CPy *other, 
        ptwXY_ptwXY_ptwXY_func func );
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( pointwiseXY_CPy *self, pointwiseXY_CPy *other, 
        ptwXY_ptwXY_ptwXY_func func );
static PyObject *pointwiseXY_C__pow__( PyObject *self, PyObject *other, PyObject *dummy );
static PyObject *pointwiseXY_C__neg__( PyObject *self );
static PyObject *pointwiseXY_C__abs__( PyObject *self );
static PyObject *pointwiseXY_C_pop( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_allocatedSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_applyFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static nfu_status pointwiseXY_C_applyFunction2( statusMessageReporting *smr, ptwXYPoint *ptwXY, void *argList );
static PyObject *pointwiseXY_C_changeInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_changeInterpolationIfNeeded( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_changeInterpolation2( pointwiseXY_CPy *self, ptwXY_interpolation interpolation, double accuracy, double lowerEps, double upperEps );
static PyObject *pointwiseXY_C_clip( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_coalescePoints( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_convolute( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_copy( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_cloneToInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_copyDataToXYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_copyDataToXsAndYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_dullEdges( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_exp( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_getAccuracy( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getBiSectionMax( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getInfill( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getInterpolation( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getNFStatus( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_lowerIndexBoundingX( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_isInterpolationLinear( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_isInterpolationOther( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getSecondaryCacheSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_getSafeDivide( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_evaluate( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_getUserFlag( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_integrate( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_integrateWithFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static nfu_status pointwiseXY_C_integrateWithFunction_callback( statusMessageReporting *smr, double x, double *y, void *argList );
static PyObject *pointwiseXY_C_integrateWithWeight_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_integrateWithWeight_sqrt_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_inverse( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_normalize( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupOneFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupTwoFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupThreeFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_groupFunctionsCommon( pointwiseXY_CPy *f1, PyObject *f2, PyObject *f3, PyObject *f4, 
        PyObject *groupBoundariesPy, PyObject *normPy );
static PyObject *pointwiseXY_C_areDomainsMutual( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_mutualify( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_overflowAllocatedSize( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_overflowLength( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_plot( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_reallocatePoints( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_reallocateOverflowPoints( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setData( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setData2( pointwiseXY_CPy *self, ptwXYPoints *ptwXY, PyObject *PyXYList );
static int pointwiseXY_C_setDataFromPtwXY( ptwXYPoints *ptwXY, pointwiseXY_CPy *otherPY );
static PyObject *pointwiseXY_C_setAccuracy( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setBiSectionMax( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setSecondaryCacheSize( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setDataFromList( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setDataFromList2( ptwXYPoints *ptwXY, PyObject *list );
static PyObject *pointwiseXY_C_setDataFromXsAndYs( pointwiseXY_CPy *self, PyObject *args );
static int pointwiseXY_C_setDataFromXsAndYs2( ptwXYPoints *ptwXY, PyObject *PyXs, PyObject *PyYs );
static PyObject *pointwiseXY_C_setInfill( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setSafeDivide( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setValue( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_setUserFlag( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_showInteralStructure( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_thicken( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_thin( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_thinDomain( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_trim( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_toString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_toString2( pointwiseXY_CPy *self, int pairsPerLine, char *format, char *pairSeparator );
static char *pointwiseXY_C_toString_isFormatForDouble( char *format, int returnFormat );
static PyObject *pointwiseXY_C_union(  pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_scaleOffsetXAndY(  pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_mergeClosePoints(  pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_domainGrid( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_domain( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_domainMin( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_domainMax( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_domainSlice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_range( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_rangeMin( pointwiseXY_CPy *self );
static PyObject *pointwiseXY_C_rangeMax( pointwiseXY_CPy *self );
static int pointwiseXY_C_domainMinMax( pointwiseXY_CPy *self, double *domainMin, double *domainMax );
static int pointwiseXY_C_rangeMinMax( pointwiseXY_CPy *self, double *rangeMin, double *rangeMax );
static PyObject *pointwiseXY_C_ysMappedToXs( pointwiseXY_CPy *self, PyObject *args );

static PyObject *pointwiseXY_C_emptyXYs( void );
static PyObject *pointwiseXY_C_createFromFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static nfu_status pointwiseXY_C_createFromFunction2( statusMessageReporting *smr, double x, double *y, void *argList );
static PyObject *pointwiseXY_C_createFromString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );

static PyObject *pointwiseXY_C_ptwXY_interpolatePoint( PyObject *self, PyObject *args );

static PyObject *pointwiseXY_C_gaussian( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords );
static PyObject *pointwiseXY_C_basicGaussian( pointwiseXY_CPy *self, PyObject *args );
static PyObject *pointwiseXY_C_unitbaseInterpolate( pointwiseXY_CPy *self, PyObject *args );

static PyObject *floatToShortestString_C( PyObject *self, PyObject *args, PyObject *keywords );

static void pointwiseXY_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 );
static int pointwiseXY_C_PythonXYPairToCPair( PyObject *XYPairPy, double *x, double *y, int64_t index );
static int64_t pointwiseXY_C_PythonXYListToCList( PyObject *PyXYList, double **xys );
static int64_t pointwiseXY_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds, int ascending );
static ptwXPoints *pointwiseXY_C_PyFloatList_to_ptwXPoints( PyObject *PyFloatList );
static PyObject *pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwXPoints *ptwX );
static int pointwiseXY_C_PyNumberToFloat( PyObject *n, double *d );
static int pyObject_NumberOrPtwXY( PyObject *other );
#if 0
static int pyObject_getNumberOrPtwXYAsPtwXY( statusMessageReporting *smr, pointwiseXY_CPy *self, PyObject *object, ptwXYPoints **ptwXY, char *variableName );
#endif
static int pointwiseXY_C_Get_pointwiseXY_CAsSelf( PyObject *self, PyObject *other, PyObject **self2, PyObject **other2,
    pointwiseXY_CPy **self3, pointwiseXY_CPy **other3 );
static int pointwiseXY_C_addedItemToPythonList( PyObject *list, PyObject *item );
static int isOkayAndHasData( pointwiseXY_CPy *self );
static int pointwiseXY_C_checkInterpolationString( char const *interpolationStr, ptwXY_interpolation *interpolation, int allowOther );
static PyObject *pointwiseXY_C_GetNone( void );
static void pointwiseXY_C_SetPyErrorExceptionFromSMR( PyObject *type, statusMessageReporting *smr );
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( const char *s, ... );
static int pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... );
static int pointwiseXY_C_checkStatus( pointwiseXY_CPy *self );
static int pointwiseXY_C_checkStatus2( pointwiseXY_CPy *self, char const *name );

/*
******************** pointwiseXY_CNewInitialize ************************
*/
static pointwiseXY_CPy *pointwiseXY_CNewInitialize( int infill, int safeDivide ) {

    pointwiseXY_CPy *self = (pointwiseXY_CPy *) PyObject_New( pointwiseXY_CPy, &pointwiseXY_CPyType );

    if( self != NULL ) {
        smr_initialize( &(self->smr), smr_status_Ok );
        self->ptwXY = NULL;
        self->infill = infill;
        self->safeDivide = safeDivide;
    }
    return( self );
}
/*
************************ __init__ ************************************
*/
static int pointwiseXY_C__init__( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int infill = 1, safeDivide = 0, userFlag = 0, status;
    int initialSize = 100, overflowSize = 10;
    double accuracy = _defaultAccuracy, biSectionMax = 3.;
    static char *kwlist[] = { "data", "dataForm", "initialSize", "overflowSize", "accuracy", "biSectionMax", "interpolation", "infill", 
        "safeDivide", "userFlag", NULL };
    ptwXYPoints *ptwXY = NULL;
    PyObject *dataPy = NULL, *dataFormPy = NULL, *xsPy = NULL, *ysPy = NULL, *theEnd = NULL, *iterator;
    char *interpolationStr = NULL, dataFormXYs[] = "xys", dataFormXsAndYs[] = "xsandys", dataFormList[] = "list", dataFormToLower[12], *c;
    char const *dataForm;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;

    self->ptwXY = NULL;
    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|OOiiddsiiiOO", kwlist, &dataPy, &dataFormPy, &initialSize, &overflowSize, &accuracy, 
        &biSectionMax, &interpolationStr, &infill, &safeDivide, &userFlag ) ) return( -1 );
    self->infill = infill;
    self->safeDivide = safeDivide;

    if( dataFormPy == NULL ) {
        dataForm = dataFormXYs; }
    else {
        if( !PYUNICODE_CHECK( dataFormPy ) ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "dataForm must be a string" ) );
        dataForm = PYUNICODE_ASSTRING( dataFormPy );
        if( strlen( dataForm ) > strlen( dataFormXsAndYs ) )
            return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid dataForm = '%s'", dataForm ) );
        strcpy( dataFormToLower, dataForm );
        for( c = dataFormToLower; *c != 0; c++ ) *c = (unsigned char ) tolower( *c );
        if( strcmp( dataFormToLower, dataFormXYs ) == 0 ) {
            dataForm = dataFormXYs; }
        else if( strcmp( dataFormToLower, dataFormXsAndYs ) == 0 ) {
            dataForm = dataFormXsAndYs ; }
        else if( strcmp( dataFormToLower, dataFormList ) == 0 ) {
            dataForm = dataFormList; }
        else {
            return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid dataForm = '%s'", dataForm ) );
        }
    }

    if( pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 1 ) != 0 ) return( -1 );

    if( ( ptwXY = ptwXY_new( NULL, interpolation, interpolationStr, biSectionMax, accuracy, initialSize, overflowSize, userFlag ) ) == NULL ) {
        PyErr_NoMemory( );
        goto err;
    }

    if( dataPy != NULL ) {
        if( dataForm == dataFormXYs ) {
            if( ( status = PyObject_IsInstance( dataPy, (PyObject* ) &pointwiseXY_CPyType ) ) == 1 ) {
                if( pointwiseXY_C_setDataFromPtwXY( ptwXY, (pointwiseXY_CPy *) dataPy ) != 0 ) goto err; }
            else {
                if( status < 0 ) goto err;
                if( pointwiseXY_C_setData2( self, ptwXY, dataPy ) != 0 ) goto err;
            } }
        else if( dataForm == dataFormXsAndYs ) {
            if( ( iterator = PyObject_GetIter( dataPy ) ) == NULL ) {
                pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs', data must be a list of length 2" );
                goto err;
            }
            xsPy = PyIter_Next( iterator );
            if( xsPy != NULL ) ysPy = PyIter_Next( iterator );
            if( ysPy != NULL ) theEnd = PyIter_Next( iterator );
            if( theEnd != NULL ) {
                pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs', data must be a list length 2, it is longer" );
                goto err; }
            else {
                if( ysPy == NULL ) {
                    pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "for dataForm = 'XsAndYs', data must be a list length 2, not 1" );
                    goto err; }
                else {
                    if( pointwiseXY_C_setDataFromXsAndYs2( ptwXY, xsPy, ysPy ) != 0 ) goto err;
                }
            }
            Py_DECREF( iterator );
            if( xsPy != NULL ) { Py_DECREF( xsPy ); }
            if( ysPy != NULL ) { Py_DECREF( ysPy ); }
            if( theEnd != NULL ) { Py_DECREF( theEnd ); } }
        else {
            if( pointwiseXY_C_setDataFromList2( ptwXY, dataPy ) != 0 ) goto err;
        }
    }
    self->ptwXY = ptwXY;
    return( 0 );

err:
    if( ptwXY != NULL ) free( ptwXY );
    return( -1 );
}
/*
************************************************************
*/
static void pointwiseXY_C_dealloc( PyObject *self ) {

    pointwiseXY_CPy *ptwXY = (pointwiseXY_CPy *) self;

    smr_release( &(ptwXY->smr) );

    if( ptwXY->ptwXY != NULL ) {
        ptwXY_free( ptwXY->ptwXY );
    }

    Py_TYPE( self )->tp_free( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__repr__( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    return( pointwiseXY_C_toString2( self, 1, " %16.8e %16.8e", "" ) );
}
/*
************************************************************
*/
static Py_ssize_t pointwiseXY_C__len__( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( -1 );
    return( (Py_ssize_t) ptwXY_length( NULL, self->ptwXY ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__getitem__( pointwiseXY_CPy *self, Py_ssize_t _index ) {

    int64_t index = (int64_t) _index;
    ptwXYPoint *point;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( ( point = ptwXY_getPointAtIndex( smr, self->ptwXY, index ) ) == NULL ) {
        if( smr_isOk( smr ) ) {
            PyErr_SetString( PyExc_IndexError, "index out of range" ); }
        else {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        }
        return( NULL );
    }
    return( Py_BuildValue( "[dd]", point->x, point->y ) );
}
/*
************************************************************
*/
static void pointwiseXY_C_getSliceIndices( int64_t length, int64_t *index1, int64_t *index2 ) {

    if( length == 0 ) {
        *index1 = 0;
        *index2 = 0; }
    else {
        if( *index1 < 0 ) {
            *index1 += length;
            if( *index1 < 0 ) *index1 = 0;
        }
        if( *index1 > length ) *index1 = length;
        if( *index2 < 0 ) {
            *index2 = *index2 + length;
            if( *index2 < 0 ) *index2 = 0;
        }
        if( *index2 > length ) *index2 = length;
        if( *index1 > *index2 ) *index2 = *index1;
    }
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__getslice__( pointwiseXY_CPy *self, Py_ssize_t _index1, Py_ssize_t _index2 ) {

    int64_t index1 = (int64_t) _index1, index2 = (int64_t) _index2, overflowSize = self->ptwXY->overflowAllocatedSize;
    pointwiseXY_CPy *newPy;
    ptwXYPoints *n;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    pointwiseXY_C_getSliceIndices( self->ptwXY->length, &index1, &index2 );
    if( ( newPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) return( NULL );
    if( overflowSize > 10 ) overflowSize = 10;
    if( ( n = ptwXY_slice( smr, self->ptwXY, index1, index2, overflowSize ) ) == NULL ) {
        Py_DECREF( newPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_IndexError, smr );
        return( NULL );
    }
    newPy->ptwXY = n;
    return( (PyObject *) newPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getslice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    long lower = 0, upper = (long) self->ptwXY->length;
    static char *kwlist[] = { "lower", "upper", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ll", kwlist, &lower, &upper ) ) return( NULL );
    return( pointwiseXY_C__getslice__( self, lower, upper ) );
}
/*
************************************************************
*/
static int pointwiseXY_C__setitem__( pointwiseXY_CPy *self, Py_ssize_t _index, PyObject *value ) {
/*
*   This routine can be called either when "self[_index] = number" or "del self[_index]".
*/
    int status = 0;
    int64_t index = (int64_t) _index;
    ptwXYPoint *point;
    pointwiseXY_CPy *other = (pointwiseXY_CPy *) value;
    double x, y;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( -1 );

    if( ( index < 0 ) || ( index >= self->ptwXY->length ) ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( -1 );
    }
    if( value == NULL ) {
        if( ptwXY_deletePoints( smr, self->ptwXY, index, index + 1 ) != nfu_Okay ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            status = -1;
        } }
    else {
        switch( pyObject_NumberOrPtwXY( value ) ) {
        case pyObject_Unsupported :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid object for set method" );
            break;
        case pyObject_Number :
            point = ptwXY_getPointAtIndex( NULL, self->ptwXY, index );
            pointwiseXY_C_PyNumberToFloat( value, &y );
            point->y = y;
            break;
        case pyObject_ptwXY :
            if( other->ptwXY->length != 1 ) {
                status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "invalid length = %ld for set object; must be 1", other->ptwXY->length ); }
            else {
                point = ptwXY_getPointAtIndex( NULL, other->ptwXY, 0 );
                if( ptwXY_setXYPairAtIndex( smr, self->ptwXY, index, point->x, point->y ) != nfu_Okay ) {
                    pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                    status = -1;
                }
            }
            break;
        default :       /* Some kind of sequence object. */
            status = pointwiseXY_C_PythonXYPairToCPair( value, &x, &y, 0 );
            if( status == 0 ) {
                 if( ptwXY_setXYPairAtIndex( smr, self->ptwXY, index, x, y ) != nfu_Okay ) {
                    pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                    status = -1; } }
            else {
                status = -1;
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C__setslice__( pointwiseXY_CPy *self, Py_ssize_t index1_, Py_ssize_t index2_, PyObject *value ) {
/*
*   This routine can be called either when "self[index1_:index2_] = number" or "del self[index1_:index2_]".
*/
    int status = 0;
    int64_t i, index1 = (int64_t) index1_, index2 = (int64_t) index2_, length;
    double x, domainMin, domainMax, *xys = NULL;
    ptwXYPoints *n = NULL;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( -1 );

    pointwiseXY_C_getSliceIndices( self->ptwXY->length, &index1, &index2 );

    if( value == NULL ) {
        if( ptwXY_deletePoints( smr, self->ptwXY, index1, index2 ) != nfu_Okay ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            status = -1;
        } }
    else {
        if( ( length = pointwiseXY_C_PythonXYListToCList( value, &xys ) ) < (int64_t) 0 ) return( -1 );
        if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) goto Err;
        if( ( self->ptwXY->length > 0 ) && ( length > 0 ) ) {
            x = domainMin = *xys;
            domainMax = xys[2 * ( length - 1 )];
            for( i = 1; i < length; i++ ) { 
                if( xys[2 * i] <= x ) {
                    smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, "x[%d] = %.17e >= x[%d] = %.17e", 
                            (int) i - 1, x, (int) i, xys[2 * i] );
                    goto Err;
                }
                x = xys[2 * i];
            }
            if( ptwXY_deletePoints( smr, n, index1, index2 ) != nfu_Okay ) goto Err;
            if( index1 > 0 ) {                                  /* Logic here requires that ptwXY_deletePoints coalesced points. */
                i = index1 - 1;
                if( index1 >= n->length ) i = n->length - 1;
                if( domainMin <= n->points[i].x ) {
                    smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, "self of x[%d] = %.17e >= slice's x[0] = %.17e", 
                            (int) i, n->points[i].x, domainMin );
                    goto Err;
                }
            }
            if( index1 < n->length ) {
                if( domainMax >= n->points[index1].x ) {
                    smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XNotAscending, "slice's x[-1] = %.17e >= self x[%d] = %.17e", 
                            domainMax, (int) index1, n->points[index1].x );
                    goto Err;
                }
            }
        }
        for( i = 0; i < length; i++ ) {
            if( ptwXY_setValueAtX( smr, n, xys[2 * i], xys[2 * i + 1] ) != nfu_Okay ) goto Err;
        }
        ptwXY_free( self->ptwXY );
        self->ptwXY = n;
    }
    return( status );

Err:
    if( n != NULL ) ptwXY_free( n );
    if( xys != NULL ) free( xys );
    self->ptwXY->status = nfu_badSelf;
    pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
    return( -1 );
}
/*
************************************************************
*/
static PySequenceMethods pointwiseXY_CPy_sequence = {
    (lenfunc) pointwiseXY_C__len__,                     /* sq_length */
    0,                                                  /* sq_concat */
    0,                                                  /* sq_repeat */
    (ssizeargfunc) pointwiseXY_C__getitem__,            /* sq_item */
#if PY_MAJOR_VERSION < 3
    (ssizessizeargfunc) pointwiseXY_C__getslice__,      /* sq_slice */
    (ssizeobjargproc) pointwiseXY_C__setitem__,         /* sq_ass_item */
    (ssizessizeobjargproc) pointwiseXY_C__setslice__,   /* sq_ass_slice */
#else
    0,                                                  /* sq_slice */
    (ssizeobjargproc) pointwiseXY_C__setitem__,         /* sq_ass_item */
    0,                                                  /* sq_ass_slice */
#endif
    0,                                                  /* sq_contains */
    0,                                                  /* sq_inplace_concat */
    0                                                   /* sq_inplace_repeat */
};
#if PY_MAJOR_VERSION < 3
    #define pointwiseXY_CPy_mappingPointer 0
#else
/*
************************************************************
*/
static PyObject *pointwiseXY_C_subscript( pointwiseXY_CPy *self, PyObject *item ) {

    if( PyIndex_Check( item ) ) {
        Py_ssize_t index;

        index = PyNumber_AsSsize_t( item, PyExc_IndexError );
        if( ( index == -1 ) && PyErr_Occurred( ) ) return( NULL );
        if( index < 0 ) index += pointwiseXY_C__len__( self );
        return( pointwiseXY_C__getitem__( self, index ) ); }
    else if ( PySlice_Check( item ) ) {
        Py_ssize_t start, stop, step, slicelength;

        if( PySlice_Unpack( item, &start, &stop, &step ) < 0 ) return( NULL );
        slicelength = PySlice_AdjustIndices( pointwiseXY_C__len__( self ), &start, &stop, step );
        if( slicelength <= 0 ) stop = start;
        if( step == 1 ) return( pointwiseXY_C__getslice__( self, start, stop ) );

        PyErr_Format( PyExc_ValueError, "step can only be 1" );
        return( NULL );
    }
    PyErr_Format( PyExc_TypeError, "indices must be integers or slices, not %.200s", item->ob_type->tp_name );
    return( NULL );
}
/*
************************************************************
*/
static int pointwiseXY_ass_subscript( pointwiseXY_CPy *self, PyObject *item, PyObject *value ) {

    if( PyIndex_Check( item ) ) {
        Py_ssize_t index;

        index = PyNumber_AsSsize_t( item, PyExc_IndexError );
        if( ( index == -1 ) && PyErr_Occurred( ) ) return( -1 );
        if( index < 0 ) index += pointwiseXY_C__len__( self );
        return( pointwiseXY_C__setitem__( self, index, value ) ); }
    else if ( PySlice_Check( item ) ) {
        Py_ssize_t start, stop, step, slicelength;

        if( PySlice_Unpack( item, &start, &stop, &step ) < 0 ) return( -1 );
        slicelength = PySlice_AdjustIndices( pointwiseXY_C__len__( self ), &start, &stop, step );
        if( slicelength <= 0 ) stop = start;
        if( step == 1 ) return( pointwiseXY_C__setslice__( self, start, stop, value ) );

        PyErr_Format( PyExc_ValueError, "step can only be 1" );
        return( -1 );
    }
    PyErr_Format( PyExc_TypeError, "indices must be integers or slices, not %.200s", item->ob_type->tp_name );
    return( -1 );
}
/*
************************************************************
*/
static PyMappingMethods pointwiseXY_CPy_mapping = {
    (lenfunc) pointwiseXY_C__len__,
    (binaryfunc) pointwiseXY_C_subscript,
    (objobjargproc) pointwiseXY_ass_subscript
};

#define pointwiseXY_CPy_mappingPointer (&pointwiseXY_CPy_mapping)
#endif
/*
************************************************************
*/
static PyObject *pointwiseXY_C__add__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s2, *o2;
    pointwiseXY_CPy *s3, *o3;

    if( pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s2, &o2, &s3, &o3 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( o2 ) ) {
    case pyObject_Number :
        n = pointwiseXY_C_add_sub_mul_div_number( s3, o2, ptwXY_add_double );
        break;
    case pyObject_ptwXY :
        n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, ptwXY_add_ptwXY );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__iadd__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    if( pointwiseXY_C_checkStatus( s1 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_add_double );
        break;
    case pyObject_ptwXY :
        if( pointwiseXY_C_checkStatus2( o1, "other" ) != 0 ) return( NULL );
        self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_add_ptwXY );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__sub__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s2, *o2;
    pointwiseXY_CPy *s3, *o3;

    if( pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s2, &o2, &s3, &o3 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( o2 ) ) {
    case pyObject_Number :
        if( s2 == self ) {
            n = pointwiseXY_C_add_sub_mul_div_number( s3, o2, ptwXY_sub_doubleFrom ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_number( s3, o2, ptwXY_sub_fromDouble );
        }
        break;
    case pyObject_ptwXY :
        n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, ptwXY_sub_ptwXY );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__isub__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    if( pointwiseXY_C_checkStatus( s1 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_sub_doubleFrom );
        break;
    case pyObject_ptwXY :
        if( pointwiseXY_C_checkStatus2( o1, "other" ) != 0 ) return( NULL );
        self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_sub_ptwXY );
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__mul__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s2, *o2;
    pointwiseXY_CPy *s3, *o3;

    if( pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s2, &o2, &s3, &o3 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( o2 ) ) {
    case pyObject_Number :
        n = pointwiseXY_C_add_sub_mul_div_number( s3, o2, ptwXY_mul_double );
        break;
    case pyObject_ptwXY :
        if( s3->infill || o3->infill ) {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, ptwXY_mul2_ptwXY ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, ptwXY_mul_ptwXY );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__imul__( PyObject *self, PyObject *other ) {

    pointwiseXY_CPy *s1 = (pointwiseXY_CPy *) self, *o1 = (pointwiseXY_CPy *) other;

    if( pointwiseXY_C_checkStatus( s1 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s1, other, ptwXY_mul_double );
        break;
    case pyObject_ptwXY :
        if( pointwiseXY_C_checkStatus2( o1, "other" ) != 0 ) return( NULL );
        if( s1->infill || o1->infill ) {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_mul2_ptwXY ); }
        else {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s1, o1, ptwXY_mul_ptwXY );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == (PyObject *) s1 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__div__( PyObject *self, PyObject *other ) {

    PyObject *n = Py_NotImplemented, *s2, *o2;
    pointwiseXY_CPy *s3, *o3;
    ptwXYPoints *numerator, *ptwXYresult;
    statusMessageReporting *smr;

    if( pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s2, &o2, &s3, &o3 ) != 0 ) return( NULL );

    smr = &(s3->smr);
    switch( pyObject_NumberOrPtwXY( o2 ) ) {
    case pyObject_Number :
        if( s2 == self ) {
            n = pointwiseXY_C_add_sub_mul_div_number( s3, other, ptwXY_div_doubleFrom ); }
        else {
            double number = 0, domainMin, domainMax;

            switch( pointwiseXY_C_domainMinMax( s3, &domainMin, &domainMax ) ) {
            case 1 :
                return( pointwiseXY_C_emptyXYs( ) );
            case -1 :
                return( NULL );
            }

            pointwiseXY_C_PyNumberToFloat( o2, &number );
            if( ( numerator = ptwXY_valueTo_ptwXY( smr, domainMin, domainMax, number ) ) == NULL ) {
                pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                return( NULL );
            }
            if( s3->safeDivide ) {
                ptwXYresult = pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( smr, numerator, s3->ptwXY ); }
            else {
                ptwXYresult = pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( smr, numerator, s3->ptwXY );
            }
            ptwXY_free( numerator );
            if( ptwXYresult == NULL ) {
                pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                n = NULL; }
            else {
                if( ( n = (PyObject *) pointwiseXY_CNewInitialize( s3->infill, s3->safeDivide ) ) == NULL ) {
                    ptwXY_free( ptwXYresult );
                    return( NULL );
                }
                ((pointwiseXY_CPy *) n)->ptwXY = ptwXYresult;
            }
        }
        break;
    case pyObject_ptwXY :
        if( s3->safeDivide || o3->safeDivide ) {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide ); }
        else {
            n = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( s3, o3, pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__idiv__( PyObject *self, PyObject *other ) {

    PyObject *s2, *o2;
    pointwiseXY_CPy *s3, *o3;

    if( pointwiseXY_C_Get_pointwiseXY_CAsSelf( self, other, &s2, &o2, &s3, &o3 ) != 0 ) return( NULL );

    switch( pyObject_NumberOrPtwXY( other ) ) {
    case pyObject_Number :
        self = pointwiseXY_C_add_sub_mul_div_number_insitu( s3, other, ptwXY_div_doubleFrom );
        break;
    case pyObject_ptwXY :
        if( s3->safeDivide || o3->safeDivide ) {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s3, o3, pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide ); }
        else {
            self = pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( s3, o3, pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide );
        }
        break;
    default :
        Py_INCREF( Py_NotImplemented );
        self = Py_NotImplemented;
    }

    if( self == s2 )  Py_INCREF( self );
    return( self );
}
/*
************************************************************
*/
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withoutSafeDivide( statusMessageReporting *smr, 
        ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    return( ptwXY_div_ptwXY( smr, ptwXY1, ptwXY2, 0 ) );
}
/*
************************************************************
*/
static ptwXYPoints *pointwiseXY_C_div_pointwiseXY_C_Py_withSafeDivide( statusMessageReporting *smr,
        ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    return( ptwXY_div_ptwXY( smr, ptwXY1, ptwXY2, 1 ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__mod__( PyObject *self, PyObject *other ) {

    double m;
    PyObject *n, *o; 
    pointwiseXY_CPy *self2 = (pointwiseXY_CPy *) self;
    statusMessageReporting *smr = &(self2->smr);

    if( pointwiseXY_C_checkStatus( self2 ) != 0 ) return( NULL );

    if( ( o = PyNumber_Float( other ) ) == NULL ) return( NULL );
    m = PyFloat_AsDouble( o );
    if( ( n = pointwiseXY_C_copy( self2 ) ) == NULL ) return( NULL );
    if( ptwXY_mod( smr, ((pointwiseXY_CPy *) n)->ptwXY, m, 1 ) ) {
        Py_DECREF( n );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_number( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func ) {

    ptwXYPoints *n = NULL;
    pointwiseXY_CPy *nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide );
    double number = 0.;
    statusMessageReporting *smr = &(self->smr);

    if( nPy != NULL ) {
        pointwiseXY_C_PyNumberToFloat( other, &number );       /* Assume calling routine has checked that other is a number. */
        n = ptwXY_clone( smr, self->ptwXY );
        if( n == NULL ) {
            Py_DECREF( nPy );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL ); }
        else {
            if( func( smr, n, number ) == nfu_Okay ) {
                nPy->ptwXY = n; }
            else {
                pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                n = ptwXY_free( n );
            }
        }
    }
    if( ( n == NULL ) && ( nPy != NULL ) ) {
        Py_DECREF( nPy );
        nPy = NULL;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_number_insitu( pointwiseXY_CPy *self, PyObject *other, ptwXY_ptwXY_d_func func ) {

    ptwXYPoints *n1 = NULL;
    double number = 0.;
    statusMessageReporting *smr = &(self->smr);

    pointwiseXY_C_PyNumberToFloat( other, &number );       /* Assume calling routine has checked that other is a number. */
    if( ( n1 = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) goto Err;
    if( func( smr, n1, number ) != nfu_Okay ) goto Err;
    ptwXY_free( self->ptwXY );
    self->ptwXY = n1;
    return( (PyObject *) self );

Err:
    if( n1 != NULL ) ptwXY_free( n1 );
    pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
    self->ptwXY->status = nfu_badSelf;
    return( NULL );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C( pointwiseXY_CPy *self, pointwiseXY_CPy *other, 
        ptwXY_ptwXY_ptwXY_func func ) {

    ptwXYPoints *n = NULL;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide );
    if( nPy != NULL ) {
        if( ( n = func( smr, self->ptwXY, other->ptwXY ) )  == NULL ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr ); }
        else {
            nPy->ptwXY = n;
        }
    }
    if( ( n == NULL ) && ( nPy != NULL ) ) {
        Py_DECREF( nPy );
        nPy = NULL;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_add_sub_mul_div_pointwiseXY_C_insitu( pointwiseXY_CPy *self, pointwiseXY_CPy *other, 
        ptwXY_ptwXY_ptwXY_func func ) {

    ptwXYPoints *n1 = NULL;
    statusMessageReporting *smr = &(self->smr);

    if( ( n1 = func( smr, self->ptwXY, other->ptwXY ) )  == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        self->ptwXY->status = nfu_badSelf;
        self = NULL; }
    else {
        ptwXY_free( self->ptwXY );
        self->ptwXY = n1;
    }
    return( (PyObject *) self );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__pow__( PyObject *self, PyObject *other, PyObject *dummy ) {

    double p;
    PyObject *o, *n;
    pointwiseXY_CPy *self2 = (pointwiseXY_CPy *) self;
    statusMessageReporting *smr = &(self2->smr);

    if( pointwiseXY_C_checkStatus( self2 ) != 0 ) return( NULL );

    if( dummy != Py_None) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "pow() 3rd argument not allowed unless all arguments are integers" ) );
    if( ( o = PyNumber_Float( other ) ) == NULL ) return( NULL );
    p = PyFloat_AsDouble( o );
    if( ( n = pointwiseXY_C_copy( self2 ) ) == NULL ) return( NULL );
    if( ptwXY_pow( smr, ((pointwiseXY_CPy *) n)->ptwXY, p ) ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        Py_DECREF( n );
        n = NULL;
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__neg__( PyObject *self ) {

    PyObject *n;
    pointwiseXY_CPy *self2 = (pointwiseXY_CPy *) self;
    statusMessageReporting *smr = &(self2->smr);

    if( pointwiseXY_C_checkStatus( self2 ) != 0 ) return( NULL );

    if( ( n = pointwiseXY_C_copy( self2 ) ) == NULL ) return( NULL );
    if( ptwXY_neg( smr, ((pointwiseXY_CPy *) n)->ptwXY ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        Py_DECREF( n );
        n = NULL;
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C__abs__( PyObject *self ) {

    PyObject *n;
    pointwiseXY_CPy *self2 = (pointwiseXY_CPy *) self;
    statusMessageReporting *smr = &(self2->smr);

    if( pointwiseXY_C_checkStatus( self2 ) != 0 ) return( NULL );

    if( ( n = pointwiseXY_C_copy( self2 ) ) == NULL ) return( NULL );
    if( ptwXY_abs( smr, ((pointwiseXY_CPy *) n)->ptwXY ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        Py_DECREF( n );
        n = NULL;
    }
    return( n );
}
/*
************************************************************
*/
static PyNumberMethods pointwiseXY_CPy_number = {
#if PY_MAJOR_VERSION < 3
    (binaryfunc) pointwiseXY_C__add__,      /* nb_add */
    (binaryfunc) pointwiseXY_C__sub__,      /* nb_subtract */
    (binaryfunc) pointwiseXY_C__mul__,      /* nb_multiply */
    (binaryfunc) pointwiseXY_C__div__,      /* nb_divide */
    (binaryfunc) pointwiseXY_C__mod__,      /* nb_remainder */
    0,                                      /* binaryfunc nb_divmod */
    (ternaryfunc) pointwiseXY_C__pow__,     /* nb_power */
    (unaryfunc) pointwiseXY_C__neg__,       /* nb_negative */
    0,                                      /* unaryfunc nb_positive */
    (unaryfunc) pointwiseXY_C__abs__,       /* nb_absolute */
    0,                                      /* inquiry nb_nonzero */
    0,                                      /* unaryfunc nb_invert */
    0,                                      /* binaryfunc nb_lshift */
    0,                                      /* binaryfunc nb_rshift */
    0,                                      /* binaryfunc nb_and */
    0,                                      /* binaryfunc nb_xor */
    0,                                      /* binaryfunc nb_or */
    0,                                      /* nb_coerce */
    0,                                      /* nb_int */
    0,                                      /* nb_long */
    0,                                      /* nb_float */
    0,                                      /* unaryfunc nb_oct */
    0,                                      /* unaryfunc nb_hex */
    (binaryfunc) pointwiseXY_C__iadd__,     /* nb_inplace_add */
    (binaryfunc) pointwiseXY_C__isub__,     /* nb_inplace_subtract */
    (binaryfunc) pointwiseXY_C__imul__,     /* nb_inplace_multiply */
    (binaryfunc) pointwiseXY_C__idiv__,     /* nb_inplace_divide */
    (binaryfunc) 0,                         /* nb_inplace_remainder */
    (ternaryfunc) 0,                        /* nb_inplace_power */

#else
    (binaryfunc) pointwiseXY_C__add__,      /* nb_add */
    (binaryfunc) pointwiseXY_C__sub__,      /* nb_subtract */
    (binaryfunc) pointwiseXY_C__mul__,      /* nb_multiply */
    (binaryfunc) pointwiseXY_C__mod__,      /* nb_remainder */
    0,                                      /* nb_divmod */
    (ternaryfunc) pointwiseXY_C__pow__,     /* nb_power */
    (unaryfunc) pointwiseXY_C__neg__,       /* nb_negative */
    0,                                      /* nb_positive */
    (unaryfunc) pointwiseXY_C__abs__,       /* nb_absolute */
    0,                                      /* nb_bool */
    0,                                      /* nb_invert */
    0,                                      /* nb_lshift */
    0,                                      /* nb_rshift */
    0,                                      /* nb_and */
    0,                                      /* nb_xor */
    0,                                      /* nb_or */
    0,                                      /* nb_int */
    0,                                      /* nb_reserved */
    0,                                      /* nb_float */

    (binaryfunc) pointwiseXY_C__iadd__,     /* nb_inplace_add */
    (binaryfunc) pointwiseXY_C__isub__,     /* nb_inplace_subtract */
    (binaryfunc) pointwiseXY_C__imul__,     /* nb_inplace_multiply */
     0,                                     /* nb_inplace_remainder */
     0,                                     /* nb_inplace_power */
     0,                                     /* nb_inplace_lshift */
     0,                                     /* nb_inplace_rshift */
     0,                                     /* nb_inplace_and */
     0,                                     /* nb_inplace_xor */
     0,                                     /* nb_inplace_or */

     0,                                     /* nb_floor_divide */
     (binaryfunc) pointwiseXY_C__div__,     /* nb_true_divide */
     0,                                     /* nb_inplace_floor_divide */
     (binaryfunc) pointwiseXY_C__idiv__,    /* nb_inplace_true_divide */
#endif
};

/*
************************************************************
*/
static PyObject *pointwiseXY_C_pop( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int index = (int) self->ptwXY->length - 1;
    ptwXYPoint point;
    static char *kwlist[] = { "index", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &index ) ) return( NULL );

    if( index < 0 ) index += self->ptwXY->length;
    if( ( index < 0 ) || ( index >= self->ptwXY->length ) ) {
        PyErr_SetString( PyExc_IndexError, "index out of range" );
        return( NULL );
    }

    point = *ptwXY_getPointAtIndex_Unsafely( self->ptwXY, index );
    if( ptwXY_deletePoints( smr, self->ptwXY, index, index + 1 ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( Py_BuildValue( "[dd]", point.x, point.y ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_allocatedSize( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "l", self->ptwXY->allocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_applyFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    nfu_status status_nf;
    PyObject *f_args[2];
    PyObject *f, *parameters;
    int biSectionMax = -1, checkForRoots = 0;
    double accuracy = -1.;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "f", "parameters", "accuracy", "biSectionMax", "checkForRoots", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OO|dii", kwlist, &f, &parameters, &accuracy, &biSectionMax, &checkForRoots ) ) return( NULL );
    f_args[0] = f;
    f_args[1] = parameters;

    if( biSectionMax < 0 ) biSectionMax = ptwXY_getBiSectionMax( self->ptwXY );
    if( accuracy < 0 ) accuracy = ptwXY_getAccuracy( self->ptwXY );

    if( !PyCallable_Check( f ) ) {
        PyErr_SetString( PyExc_TypeError, "First argument must be callable" );
        return( NULL );
    }

    if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    ptwXY_setBiSectionMax( n, biSectionMax );
    ptwXY_setAccuracy( n, accuracy );

    Py_INCREF( f );
    status_nf = ptwXY_applyFunction( smr, n, pointwiseXY_C_applyFunction2, (void *) f_args, checkForRoots );
    Py_DECREF( f );

    if( status_nf != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        ptwXY_free( n );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static nfu_status pointwiseXY_C_applyFunction2( statusMessageReporting *smr, ptwXYPoint *ptwXY, void *argList ) {

    PyObject **f_args = (PyObject **) argList, *result;
    nfu_status status_nf = nfu_Okay;

    result = PyObject_CallFunction( (PyObject *) f_args[0], "(d,O)", ptwXY->x, f_args[1] );
    if( result == NULL ) return( nfu_badInput );
    if( pointwiseXY_C_PyNumberToFloat( result, &(ptwXY->y) ) != 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badInput, "could not convert returned value to float" );
        status_nf = nfu_badInput;
    }
    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    char *interpolationStr = NULL;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;
    double accuracy = self->ptwXY->accuracy, lowerEps = 0., upperEps = 0.;
    static char *kwlist[] = { "interpolation", "accuracy", "lowerEps", "upperEps", NULL };

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|sddd", kwlist, &interpolationStr, &accuracy, &lowerEps, &upperEps ) ) return( NULL );

    if( pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 0 ) != 0 ) return( NULL );
    return( pointwiseXY_C_changeInterpolation2( self, interpolation, accuracy, lowerEps, upperEps ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolationIfNeeded( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int status, length = 0, i1;
    double accuracy = self->ptwXY->accuracy, lowerEps = 0., upperEps = 0.;
    char const *interpolationStr;
    static char *kwlist[] = { "interpolation", "accuracy", "lowerEps", "upperEps", NULL };
    PyObject *allowedInterpolations, *iterator, *interpolationItem = NULL;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin, firstInterpolation = ptwXY_interpolationLinLin;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|ddd", kwlist, &allowedInterpolations, &accuracy, &lowerEps, &upperEps ) ) return( NULL );

    if( ( iterator = PyObject_GetIter( allowedInterpolations ) ) == NULL ) return( NULL );
    for( interpolationItem = PyIter_Next( iterator ), i1 = 0; interpolationItem != NULL; 
            interpolationItem = PyIter_Next( iterator ), ++i1 ) {
        if( ( interpolationStr = PYUNICODE_ASSTRING( interpolationItem ) ) == NULL ) {
            Py_DECREF( interpolationItem );
            Py_DECREF( iterator );
            smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                    "Interpolation item as index %d is not a string", i1 );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
        status = pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 0 );
        Py_DECREF( interpolationItem );
        if( status < 0 ) {
            Py_DECREF( iterator );
            return( NULL );
        }
        if( interpolation == self->ptwXY->interpolation ) {
            Py_DECREF( iterator );
            Py_INCREF( self );
            return( (PyObject *) self );
        }
        if( length == 0 ) firstInterpolation = interpolation;
        ++length;
    }
    Py_DECREF( iterator );
    if( length == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Allowed interpolation sequence is empty" ) );
    return( pointwiseXY_C_changeInterpolation2( self, firstInterpolation, accuracy, lowerEps, upperEps ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_changeInterpolation2( pointwiseXY_CPy *self, ptwXY_interpolation interpolation, double accuracy, double lowerEps, double upperEps ) {

    ptwXYPoints *n1;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( self->ptwXY->interpolation == ptwXY_interpolationFlat ) {
        if( interpolation != ptwXY_interpolationLinLin )
            return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: can only change 'flat' to 'lin-lin'" ) );
        if( ( n1 = ptwXY_flatInterpolationToLinear( smr, self->ptwXY, lowerEps, upperEps ) ) == NULL ) {
            if( ( lowerEps == 0 ) && ( upperEps == 0 ) )
                return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from changeInterpolation: both lowerEps and upperEps are 0." ) );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        } }
    else {
        if( ( n1 = ptwXY_toOtherInterpolation( smr, self->ptwXY, interpolation, accuracy ) ) == NULL ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
    }

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n1 );
        return( NULL );
    }
    nPy->ptwXY = n1;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_clip( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    double rangeMin, rangeMax;
    static char *kwlist[] = { "rangeMin", "rangeMax", NULL };
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( pointwiseXY_C_rangeMinMax( self, &rangeMin, &rangeMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &rangeMin, &rangeMax ) ) return( NULL );

    if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( self->ptwXY->length != 0 ) {
        if( ptwXY_clip( smr, n, rangeMin, rangeMax ) != nfu_Okay ) {
            ptwXY_free( n );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
    }

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_coalescePoints( pointwiseXY_CPy *self ) {

    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ptwXY_coalescePoints( smr, self->ptwXY, 0, NULL, 0 ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_convolute( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int mode = 0;
    ptwXYPoints *n;
    PyObject *otherPy;
    pointwiseXY_CPy *nPy, *other;
    static char *kwlist[] = { "other", "mode", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|i", kwlist, &otherPy, &mode ) ) return( NULL );

    other = (pointwiseXY_CPy *) otherPy;
    if( pointwiseXY_C_checkStatus2( other, "other" ) != 0 ) return( NULL );

    if( ( n = ptwXY_convolution( smr, self->ptwXY, ((pointwiseXY_CPy *) otherPy)->ptwXY, mode ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copy( pointwiseXY_CPy *self ) {

    ptwXYPoints *ptwXY;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ( ptwXY = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( ptwXY );
        return( NULL );
    }
    nPy->ptwXY = ptwXY;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_cloneToInterpolation( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    char *interpolationStr = NULL;
    ptwXY_interpolation interpolation;
    static char *kwlist[] = { "interpolation", NULL };
    statusMessageReporting *smr = &(self->smr);
    ptwXYPoints *cloned;
    pointwiseXY_CPy *clonedPy;

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|s", kwlist, &interpolationStr ) ) return( NULL );
    if( pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 1 ) != 0 ) return( NULL );

    if( ( cloned = ptwXY_cloneToInterpolation( smr, self->ptwXY, interpolation ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    if( ( clonedPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( cloned );
        return( NULL );
    }
    clonedPy->ptwXY = cloned;
    return( (PyObject *) clonedPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copyDataToXYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int64_t i1;
    PyObject *l = PyList_New( 0 ), *ni;
    double xScale = 1.0, yScale = 1.0, *xs, *ys;
    static char *kwlist[] = { "xScale", "yScale", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xScale, &yScale ) ) return( NULL );

    if( l == NULL ) return( NULL );
    if( ptwXY_coalescePoints( smr, self->ptwXY, 0, NULL, 0 ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ptwXY_valuesToC_XsAndYs( smr, self->ptwXY, &xs, &ys ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    for( i1 = 0; i1 < self->ptwXY->length; i1++ ) {
        ni = Py_BuildValue( "[dd]", xScale * xs[i1], yScale * ys[i1] );
        if( pointwiseXY_C_addedItemToPythonList( l, ni ) != 0 ) {
            free( xs );
            free( ys );
            return( NULL );
        }
    }
    free( xs );
    free( ys );
    return( l );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_copyDataToXsAndYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int64_t i;
    PyObject *lx, *ly, *item;
    double xScale = 1.0, yScale = 1.0;
    static char *kwlist[] = { "xScale", "yScale", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &xScale, &yScale ) ) return( NULL );

    if( ptwXY_coalescePoints( smr, self->ptwXY, 0, NULL, 0 ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    if( ( lx = PyList_New( 0 ) ) == NULL ) return( NULL );
    if( ( ly = PyList_New( 0 ) ) == NULL ) {
        Py_DECREF( lx );
        return( NULL );
    }

    for( i = 0; i < self->ptwXY->length; i++ ) {
        item = Py_BuildValue( "d", xScale * self->ptwXY->points[i].x );
        if( pointwiseXY_C_addedItemToPythonList( lx, item ) != 0 ) {
            Py_DECREF( ly );
            return( NULL );
        }
        item = Py_BuildValue( "d", yScale * self->ptwXY->points[i].y );
        if( pointwiseXY_C_addedItemToPythonList( ly, item ) != 0 ) {
            Py_DECREF( lx );
            return( NULL );
        }
    }

    item = Py_BuildValue( "OO", lx, ly );
    Py_DECREF( lx );
    Py_DECREF( ly );

    return( item );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_dullEdges( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "lowerEps", "upperEps", "positiveXOnly", NULL };
    int positiveXOnly = 0;
    double lowerEps = 0., upperEps = 0.;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ddi", kwlist, &lowerEps, &upperEps, &positiveXOnly ) ) return( NULL );

    if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    if( ptwXY_dullEdges( smr, n, lowerEps, upperEps, positiveXOnly ) != nfu_Okay ) {
        ptwXY_free( n );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_exp( pointwiseXY_CPy *self, PyObject *args ) {

    double a;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &a ) ) return( NULL );

    if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    if( ptwXY_exp( smr, n, a ) != nfu_Okay ) {
        ptwXY_free( n );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getAccuracy( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "d", ptwXY_getAccuracy( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getBiSectionMax( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "d", (double) ptwXY_getBiSectionMax( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getInfill( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "i", self->infill ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getInterpolation( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "s", ptwXY_getInterpolationString( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setInterpolation( pointwiseXY_CPy *self, PyObject *args ) {

    char *interpolationString;

    if( !PyArg_ParseTuple( args, "s", &interpolationString ) ) return( NULL );

    if( ptwXY_setInterpolationString( self->ptwXY, interpolationString ) == nfu_Okay ) return( pointwiseXY_C_GetNone( ) );

    return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Other interpolation not support for setInterpolation: interpolation = '%s'", interpolationString ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getNFStatus( pointwiseXY_CPy *self ) {

    return( Py_BuildValue( "i", ptwXY_getStatus( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_lowerIndexBoundingX( pointwiseXY_CPy *self, PyObject *args ) {

    int64_t index;
    double x;
    statusMessageReporting *smr = &(self->smr);

    if( !PyArg_ParseTuple( args, "d", &x ) ) return( NULL );

    if( ptwXY_getLowerIndexBoundingX( smr, self->ptwXY, x, &index ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    return( Py_BuildValue( "i", (int) index ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_isInterpolationLinear( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "i", self->ptwXY->interpolation == ptwXY_interpolationLinLin ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_isInterpolationOther( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "i", self->ptwXY->interpolation == ptwXY_interpolationOther ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getSecondaryCacheSize( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "l", (long int) self->ptwXY->overflowAllocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getSafeDivide( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "i", self->safeDivide ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_evaluate( pointwiseXY_CPy *self, PyObject *args ) {

    double x, y;
    nfu_status status_nf;
    PyObject *value = NULL;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &x ) ) return( NULL );

    switch( status_nf = ptwXY_getValueAtX( smr, self->ptwXY, x, &y ) ) {
    case nfu_Okay :
        value = Py_BuildValue( "d", y );
        break;
    case nfu_XOutsideDomain :
        value = pointwiseXY_C_GetNone( );
        break;
    case nfu_invalidInterpolation :
        value = pointwiseXY_C_SetPyErrorExceptionReturnNull( "unsupported interpolation = %d", self->ptwXY->interpolation );
        break;
    default :
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
    }

    return( value );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_getUserFlag( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "i", ptwXY_getUserFlag( self->ptwXY ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrate( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "domainMin", "domainMax", NULL };
    double domainMin, domainMax, integral;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &domainMin, &domainMax ) ) return( NULL );

    if( ptwXY_integrate( smr, self->ptwXY, domainMin, domainMax, &integral ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", integral ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrateWithFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int degree = 4, recursionLimit = 10;
    double tolerance, domainMin, domainMax, integral;
    PyObject *func_Py, *parameters_Py = Py_None;
    static char *kwlist[] = { "f", "tolerance", "parameters", "domainMin", "domainMax", "degree", "recursionLimit", NULL };
    GnG_parameters argList;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "Od|Oddii", kwlist, &func_Py, &tolerance, &parameters_Py, 
        &domainMin, &domainMax, &degree, &recursionLimit ) ) return( NULL );

    argList.func = func_Py;
    argList.argList = parameters_Py;

    if( ptwXY_integrateWithFunction( smr, self->ptwXY, pointwiseXY_C_integrateWithFunction_callback, &argList,
            domainMin, domainMax, degree, recursionLimit, tolerance, &integral ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", integral ) );
}
/*
************************************************************
*/
static nfu_status pointwiseXY_C_integrateWithFunction_callback( statusMessageReporting *smr, double x, double *y, void *argList ) {
/*
BRB FIXME , check if smr is needed.
*/

    nfu_status status_nf = nfu_Okay;
    GnG_parameters *parameters = (GnG_parameters *) argList;
    PyObject *result;

    result = PyObject_CallFunction( parameters->func, "(d,O)", x, parameters->argList );
    if( result == NULL ) return( nfu_badInput );

    if( PyFloat_Check( result ) ) {
        *y = PyFloat_AsDouble( result ); }
    else {
        status_nf = nfu_badInput;
    }

    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrateWithWeight_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "domainMin", "domainMax", NULL };
    double domainMin, domainMax, integral;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &domainMin, &domainMax ) ) return( NULL );

    if( ptwXY_integrateWithWeight_x( smr, self->ptwXY, domainMin, domainMax, &integral ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", integral ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_integrateWithWeight_sqrt_x( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    static char *kwlist[] = { "domainMin", "domainMax", NULL };
    double domainMin, domainMax, integral;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( self->ptwXY->length == 0 ) return( Py_BuildValue( "d", 0. ) );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|dd", kwlist, &domainMin, &domainMax ) ) return( NULL );

    if( ptwXY_integrateWithWeight_sqrt_x( smr, self->ptwXY, domainMin, domainMax, &integral ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", integral ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_inverse( pointwiseXY_CPy *self ) {

    ptwXYPoints *ptwXY;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ( ptwXY = ptwXY_inverse( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( ptwXY );
        return( NULL );
    }
    nPy->ptwXY = ptwXY;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_normalize( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int insitu = 0;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy = NULL;
    static char *kwlist[] = { "insitu", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &insitu ) ) return( NULL );

    if( insitu == 0 ) {
        if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
        if( ptwXY_normalize( smr, n ) != nfu_Okay ) {
            ptwXY_free( n );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
        if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
            ptwXY_free( n );
            return( NULL );
        }
        nPy->ptwXY = n; }
    else {
        if( ptwXY_normalize( smr, self->ptwXY ) != nfu_Okay ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        }
        nPy = self;
        Py_INCREF( self );      /* FIXME: Why is this needed when insitu if True? Or is it? */
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_runningIntegral( pointwiseXY_CPy *self ) {

    ptwXPoints *ptwX;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ( ptwX = ptwXY_runningIntegral( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwX ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupOneFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL;
    static char *kwlist[] = { "groupBoundaries", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|O", kwlist, &groupBoundariesPy, &normPy ) ) return( NULL );
    return( pointwiseXY_C_groupFunctionsCommon( self, NULL, NULL, NULL, groupBoundariesPy, normPy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupTwoFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL, *f2;
    static char *kwlist[] = { "groupBoundaries", "f2", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OO|O", kwlist, &groupBoundariesPy, &f2, &normPy ) ) return( NULL );
    return( pointwiseXY_C_groupFunctionsCommon( self, f2, NULL, NULL, groupBoundariesPy, normPy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupThreeFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL, *f2, *f3;
    static char *kwlist[] = { "groupBoundaries", "f2", "f3", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OOO|O", kwlist, &groupBoundariesPy, &f2, &f3, &normPy ) ) return( NULL );
    return( pointwiseXY_C_groupFunctionsCommon( self, f2, f3, NULL, groupBoundariesPy, normPy ) );
}

/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupFourFunctions( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *groupBoundariesPy, *normPy = NULL, *f2, *f3, *f4;
    static char *kwlist[] = { "groupBoundaries", "f2", "f3", "f4", "norm", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OOOO|O", kwlist, &groupBoundariesPy, &f2, &f3, &f4, &normPy ) ) return( NULL );
    return( pointwiseXY_C_groupFunctionsCommon( self, f2, f3, f4, groupBoundariesPy, normPy ) );
}

/*
************************************************************
*/
static PyObject *pointwiseXY_C_groupFunctionsCommon( pointwiseXY_CPy *f1, PyObject *of2, PyObject *of3, PyObject *of4,
                PyObject *groupBoundariesPy, PyObject *normPy ) {

    ptwXPoints *ptwXGBs, *ptwX_norm = NULL, *groups = NULL;
    pointwiseXY_CPy *f2 = (pointwiseXY_CPy *) of2, *f3 = (pointwiseXY_CPy *) of3, *f4 = (pointwiseXY_CPy *) of4;
    PyObject *newPy = NULL;
    ptwXY_group_normType norm = ptwXY_group_normType_none;
    char const *normChars;
    int status;
    statusMessageReporting *smr = &(f1->smr);

    if( pointwiseXY_C_checkStatus( f1 ) != 0 ) return( NULL );
    if( f2 != NULL ) {
        if( ( status = PyObject_IsInstance( of2, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
        if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "f2 must be a pointwiseXY_C instance" ) );
        if( pointwiseXY_C_checkStatus2( f2, "f2" ) != 0 ) return( NULL );
        if( f3 != NULL ) {
            if( ( status = PyObject_IsInstance( of3, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
            if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "f3 must be a pointwiseXY_C instance" ) );
            if( pointwiseXY_C_checkStatus2( f3, "f3" ) != 0 ) return( NULL );
            if( f4 != NULL ) {
                if( ( status = PyObject_IsInstance( of4, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
                if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "f4 must be a pointwiseXY_C instance" ) );
                if( pointwiseXY_C_checkStatus2( f4, "f4" ) != 0 ) return( NULL );
            }
        }
    }
/*
BRB FIXME, need to check status of normPy
*/
    if( ( normPy != NULL ) && ( normPy != Py_None ) ) {
        if( PYUNICODE_CHECK( normPy ) ) {
            if( ( normChars = PYUNICODE_ASSTRING( normPy ) ) == NULL ) {
                return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid string norm" ) ); }
            else {
                if( strcmp( "dx", normChars ) != 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid string norm identifier, must be 'dx'" ) );
                norm = ptwXY_group_normType_dx;
            } }
        else {
            if( ( ptwX_norm = pointwiseXY_C_PyFloatList_to_ptwXPoints( normPy ) ) == NULL ) return( NULL );
            norm = ptwXY_group_normType_norm;
        }
    }

    if( ( ptwXGBs = pointwiseXY_C_PyFloatList_to_ptwXPoints( groupBoundariesPy ) ) == NULL ) {
        if( ptwX_norm != NULL ) ptwX_free( ptwX_norm );
        return( NULL );
    }

    if( f2 == NULL ) {
        groups = ptwXY_groupOneFunction( smr, f1->ptwXY, ptwXGBs, norm, ptwX_norm ); }
    else if( f3 == NULL ) {
        groups = ptwXY_groupTwoFunctions( smr, f1->ptwXY, ((pointwiseXY_CPy *) f2)->ptwXY, ptwXGBs, norm, ptwX_norm ); }
    else if( f4 == NULL ) {
        groups = ptwXY_groupThreeFunctions( smr, f1->ptwXY, ((pointwiseXY_CPy *) f2)->ptwXY, ((pointwiseXY_CPy *) f3)->ptwXY, ptwXGBs, norm, ptwX_norm ); }
    else {
        groups = ptwXY_groupFourFunctions( smr, f1->ptwXY, f2->ptwXY, f3->ptwXY, f4->ptwXY, ptwXGBs, norm, ptwX_norm );
    }

    ptwX_free( ptwXGBs );
    if( ptwX_norm != NULL ) ptwX_free( ptwX_norm );

    if( groups == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr ); }
    else {
        newPy = pointwiseXY_C_ptwXPoints_to_PyFloatList( groups );
        ptwX_free( groups );
    }

    return( newPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_areDomainsMutual( pointwiseXY_CPy *self, PyObject *args ) {

    nfu_status status_nf;
    PyObject *otherPy, *rPy = NULL;
    pointwiseXY_CPy *other;
    int status;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "O", &otherPy ) ) return( NULL );
    other = (pointwiseXY_CPy *) otherPy;

    if( ( status = PyObject_IsInstance( otherPy, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "other must be a pointwiseXY_C instance" ) );
    if( pointwiseXY_C_checkStatus2( other, "other" ) != 0 ) return( NULL );

    switch( status_nf = ptwXY_areDomainsMutual( NULL, self->ptwXY, other->ptwXY ) ) {
    case nfu_Okay :
    case nfu_empty :
        if( ( rPy = Py_BuildValue( "i", 1 ) ) == NULL ) return( NULL );
        break;
    case nfu_domainsNotMutual :
    case nfu_tooFewPoints :
        if( ( rPy = Py_BuildValue( "i", 0 ) ) == NULL ) return( NULL );
        break;
    default :
        rPy = pointwiseXY_C_SetPyErrorExceptionReturnNull( "Error from ptwXY_areDomainsMutual: %s", nfu_statusMessage( status_nf ) );
    }

    return( rPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_mutualify( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int positiveXOnly1, positiveXOnly2, status;
    double lowerEps1, upperEps1, lowerEps2, upperEps2;
    static char *kwlist[] = { "lowerEps1", "upperEps1", "positiveXOnly1", "other", "lowerEps2", "upperEps2", "positiveXOnly2", NULL };
    PyObject *other, *lPy;
    ptwXYPoints *n1 = NULL, *n2 = NULL;
    pointwiseXY_CPy *n1Py = NULL, *n2Py = NULL, *other2;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "ddiOddi", kwlist, &lowerEps1, &upperEps1, &positiveXOnly1, &other, &lowerEps2, 
        &upperEps2, &positiveXOnly2 ) ) return( NULL );

    if( ( status = PyObject_IsInstance( other, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "other must be a pointwiseXY_C instance" ) );
    other2 = (pointwiseXY_CPy *) other;
    if( pointwiseXY_C_checkStatus2( other2, "other" ) != 0 ) return( NULL );

    if( ( n1 = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) goto ErrSMR;
    if( ( n2 = ptwXY_clone( smr, other2->ptwXY ) ) == NULL ) goto ErrSMR;

    if( ptwXY_mutualifyDomains( smr, n1, lowerEps1, upperEps1, positiveXOnly1, n2, lowerEps2, upperEps2, positiveXOnly2 ) != nfu_Okay )
        goto ErrSMR;

    if( ( n1Py = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) goto Err;
    if( ( n2Py = pointwiseXY_CNewInitialize( other2->infill, other2->safeDivide ) ) == NULL ) goto Err;
    if( ( lPy = Py_BuildValue( "OO", n1Py, n2Py ) ) == NULL ) goto Err;
    n1Py->ptwXY = n1;
    n2Py->ptwXY = n2;
    Py_DECREF( n1Py );          /* This is needed as the tuple in Py_BuildValue is adding a reference counts, making the count for the object 2 and not 1. */
    Py_DECREF( n2Py );

    return( lPy );

ErrSMR:
    pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
Err:
    if( n1 != NULL ) ptwXY_free( n1 );
    if( n2 != NULL ) ptwXY_free( n2 );
    if( n1Py != NULL ) { Py_DECREF( n1Py ); }
    if( n2Py != NULL ) { Py_DECREF( n2Py ); }
    return( NULL );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_overflowAllocatedSize( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "l", self->ptwXY->overflowAllocatedSize ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_overflowLength( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    return( Py_BuildValue( "l", self->ptwXY->overflowLength ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_plot( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *moduleName, *GnuplotModule, *Gnuplot = NULL, *Data = NULL, *g = NULL, *data = NULL, *status = NULL;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ( moduleName = PYUNICODE_FROMSTRING( "Gnuplot" ) ) == NULL ) return( NULL );
    GnuplotModule = PyImport_Import( moduleName );
    Py_DECREF( moduleName );
    if( GnuplotModule == NULL ) return( NULL );

    if( ( Gnuplot = PyObject_GetAttrString( GnuplotModule, "Gnuplot" ) ) == NULL ) goto err;
    if( ( g = PyObject_CallFunction( Gnuplot, "()" ) ) == NULL ) goto err;
    if( PyObject_CallFunction( g, "(s)", "set style data linespoints" ) == NULL ) goto err;

    if( ( Data = PyObject_GetAttrString( GnuplotModule, "Data" ) ) == NULL ) goto err;

    if( ptwXY_length( NULL, self->ptwXY ) > 0 ) {
        if( ( data = PyObject_CallFunction( Data, "(O)", self ) ) == NULL ) goto err; }
    else {
        if( ( data = PyObject_CallFunction( Data, "([[d,d]])", -.9999, -.9999 ) ) == NULL ) goto err;
        PyObject_CallFunction( g, "(s)", "set xrange [-1:1]" );
        PyObject_CallFunction( g, "(s)", "set yrange [-1:1]" );
        PyObject_CallFunction( g, "(s)", "set label 'Instance has no data' center at 0,0" );
    }
    if( PyObject_CallMethod( g, "plot", "(O)", data ) == NULL ) goto err;

    status = g;
    goto theEnd;

err:
    status = NULL;

theEnd:
    Py_DECREF( GnuplotModule );
    if( Data != NULL ) { Py_DECREF( Data ); }
    if( data != NULL ) { Py_DECREF( data ); }
    if( Gnuplot != NULL ) { Py_DECREF( Gnuplot ); }
    if( ( status == NULL ) && ( g != NULL ) ) { Py_DECREF( g ); }
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_reallocatePoints( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int size;
    int forceSmallerResize = 1;
    static char *kwlist[] = { "size", "forceSmaller", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "i|i", kwlist, &size, &forceSmallerResize ) ) return( NULL );

    if( ptwXY_reallocatePoints( smr, self->ptwXY, size, forceSmallerResize ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_reallocateOverflowPoints( pointwiseXY_CPy *self, PyObject *args ) {

    int size;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "i", &size ) ) return( NULL );

    if( ptwXY_reallocateOverflowPoints( smr, self->ptwXY, size ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setData( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXYList;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "O", &PyXYList ) ) return( NULL );
    if( pointwiseXY_C_setData2( self, self->ptwXY, PyXYList ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setData2( pointwiseXY_CPy *self, ptwXYPoints *ptwXY, PyObject *PyXYList ) {

    int status = 0;
    int64_t length;
    double *xys;
    statusMessageReporting *smr = &(self->smr);

    length = pointwiseXY_C_PythonXYListToCList( PyXYList, &xys );
    if( length == -1 ) return( -1 );
    if( ptwXY_setXYData( smr, ptwXY, length, xys ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( -1 );
    }
    free( xys );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromPtwXY( ptwXYPoints *ptwXY, pointwiseXY_CPy *otherPy ) {

    statusMessageReporting *smr = &(otherPy->smr);

    if( ptwXY_copyPointsOnly( smr, ptwXY, otherPy->ptwXY ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setAccuracy( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );

    return( Py_BuildValue( "d", ptwXY_setAccuracy( self->ptwXY, accuracy ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setBiSectionMax( pointwiseXY_CPy *self, PyObject *args ) {

    double biSectionMax;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &biSectionMax ) ) return( NULL );

    return( Py_BuildValue( "d", ptwXY_setBiSectionMax( self->ptwXY, biSectionMax ) ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setSecondaryCacheSize( pointwiseXY_CPy *self, PyObject *args ) {

    long int length;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "l", &length ) ) return( NULL );

    if( ptwXY_reallocateOverflowPoints( smr, self->ptwXY, length ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setDataFromList( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXYs;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "O", &PyXYs) ) return( NULL );
    if( pointwiseXY_C_setDataFromList2( self->ptwXY, PyXYs ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromList2( ptwXYPoints *ptwXY, PyObject *list ) {

    int status = 0;
    int64_t i, j, length;
    double *xys = NULL, x, xPrior = 0.;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    ptwXY->length = 0;
    if( ( length = pointwiseXY_C_pythonDoubleListToCList( list, &xys, 0 ) ) == -1 ) return( -1 );
    if( 2 * ( length / 2 ) != length ) {
        status = -1;
        free( xys );
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "length = %d is not even", length ); }
    else {
        if( ptwXY_reallocatePoints( &smr, ptwXY, length / 2, 1 ) == nfu_Okay ) {
            for( i = 0, j = 0; i < length; i += 2, j++ ) {
                x = xys[i];
                if( i > 0 ) {
                    if( x <= xPrior ) {
                        pointwiseXY_C_SetPyErrorExceptionReturnNull( "data not in ascending order at index %d and %d", (int) i - 1, (int) i );
                        status = -1;
                        break;
                    }
                }
                ptwXY->points[j].x = x;
                ptwXY->points[j].y = xys[i+1];
                xPrior = x;
            }
            ptwXY->length = length / 2; }
        else {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
            status = -1;
        }
        free( xys );
    }
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setDataFromXsAndYs( pointwiseXY_CPy *self, PyObject *args ) {

    PyObject *status = NULL;
    PyObject *PyXs, *PyYs;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "OO", &PyXs, &PyYs) ) return( NULL );
    if( pointwiseXY_C_setDataFromXsAndYs2( self->ptwXY, PyXs, PyYs ) == 0 ) status = pointwiseXY_C_GetNone( );
    return( status );
}
/*
************************************************************
*/
static int pointwiseXY_C_setDataFromXsAndYs2( ptwXYPoints *ptwXY, PyObject *PyXs, PyObject *PyYs ) {

    int status = 0;
    int64_t i, xLength, yLength;
    double *xs = NULL, *ys = NULL;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    ptwXY->length = 0;
    if( ( xLength = pointwiseXY_C_pythonDoubleListToCList( PyXs, &xs, 1 ) ) == -1 ) return( -1 );
    if( ( yLength = pointwiseXY_C_pythonDoubleListToCList( PyYs, &ys, 0 ) ) == -1 ) {
        free( xs );
        return( -1 );
    }
    if( xLength != yLength ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "x and y lengths are not the same" );
        status = -1;
    }
    if( status == 0 ) {
        if( ptwXY_reallocatePoints( &smr, ptwXY, xLength, 1 ) == nfu_Okay ) {
            for( i = 0; i < xLength; i++ ) {
                ptwXY->points[i].x = xs[i];
                ptwXY->points[i].y = ys[i];
            }
            ptwXY->length = xLength; }
        else {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
            status = -1;
        }
    }
    if( xs != NULL ) free( xs );
    if( ys != NULL ) free( ys );
    return( status );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setInfill( pointwiseXY_CPy *self, PyObject *args ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "i", &(self->infill) ) ) return( NULL );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setSafeDivide( pointwiseXY_CPy *self, PyObject *args ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "i", &(self->safeDivide) ) ) return( NULL );

    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setValue( pointwiseXY_CPy *self, PyObject *args ) {

    double x, y;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "dd", &x, &y ) ) return( NULL );

    if( ptwXY_setValueAtX( smr, self->ptwXY, x, y ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_setUserFlag( pointwiseXY_CPy *self, PyObject *args ) {

    int flag;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "i", &flag ) ) return( NULL );

    ptwXY_setUserFlag( self->ptwXY, flag );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_showInteralStructure( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int printPointersAsNull = 0;
    static char *kwlist[] = { "printPointersAsNull", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|i", kwlist, &printPointersAsNull ) ) return( NULL );

    ptwXY_showInteralStructure( self->ptwXY, stdout, printPointersAsNull );
    return( pointwiseXY_C_GetNone( ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_thicken( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int sectionSubdivideMax = 1;
    double dDomainMax = 0, fDomainMax = 1.;
    static char *kwlist[] = { "sectionSubdivideMax", "dDomainMax", "fDomainMax", NULL };
    PyObject *n;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|idd", kwlist, &sectionSubdivideMax, &dDomainMax, &fDomainMax ) ) return( NULL );
    if( ( n = pointwiseXY_C_copy( self ) ) == NULL ) return( NULL );
    if( ptwXY_thicken( smr, ((pointwiseXY_CPy *) n)->ptwXY, sectionSubdivideMax, dDomainMax, fDomainMax ) ) {
        Py_DECREF( n );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( n );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_thin( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );

    if( ( n = ptwXY_thin( smr, self->ptwXY, accuracy ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_thinDomain( pointwiseXY_CPy *self, PyObject *args ) {

    double epsilon;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &epsilon ) ) return( NULL );

    if( ( n = ptwXY_thinDomain( smr, self->ptwXY, epsilon ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_toString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int pairsPerLine = 1;
    char *format = " %16.8e %16.8e", *pairSeparator = "";
    static char *kwlist[] = { "pairsPerLine", "format", "pairSeparator", NULL };

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|iss", kwlist, &pairsPerLine, &format, &pairSeparator ) ) return( NULL );
    return( pointwiseXY_C_toString2( self, pairsPerLine, format, pairSeparator ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_toString2( pointwiseXY_CPy *self, int pairsPerLine, char *format, char *pairSeparator ) {

    int lineFeed, extraSpace = 4, dummyLen;     /* '\n' plus three more for reserve. */
    ptwXYPoints *ptwXY = self->ptwXY;
    int64_t index, strLength, pairSeparatorLength = strlen( pairSeparator );
    char *s, *e, *p, dummy[1024];
    ptwXYPoint *point;
    PyObject *str;

    s = pointwiseXY_C_toString_isFormatForDouble( format, 0 );
    if( ( s = pointwiseXY_C_toString_isFormatForDouble( s, 0 ) ) != NULL ) {
        s = pointwiseXY_C_toString_isFormatForDouble( s, 1 );
/*
        snprintf( dummy, 0, "'%s' %% ( -3.1238972718972345e-21, -6.1432987987397648974e+33 )", format );
        if( PyRun_SimpleString( dummy ) < 0 ) s = NULL;
*/
    }
    if( s == NULL ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "invalid format = '%s' for converting two doubles to a string", format ) );

    dummyLen = snprintf( dummy, 0, format, -3.1238972718972345e-21, -6.1432987987397648974e-33 );
    if( dummyLen < 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "Invalid format string = '%s'", format ) );
    strLength = ptwXY->length * ( dummyLen + pairSeparatorLength + extraSpace );

    if( ( s = (char *) malloc( (size_t) strLength * sizeof( char ) ) ) == NULL ) {
        PyErr_NoMemory( );
        return( NULL );
    }
    s[0] = 0;
    for( index = 0, e = s, lineFeed = pairsPerLine; index < ptwXY->length; index++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( ptwXY, index );
        sprintf( e, format, point->x, point->y );
        while( *e ) e++;
        if( index < ptwXY->length - 1 ) {
            for( p = pairSeparator; *p; p++, e++ ) *e = *p;
            *e = 0;
        }
        lineFeed--;
        if( lineFeed <= 0 ) {
            *e = '\n';
            e++;
            *e = 0;
            lineFeed = pairsPerLine;
        }
    }
    if( lineFeed != pairsPerLine ) {
        *e = '\n';
        e++;
        *e = 0;
    }
    str = PYUNICODE_FROMSTRING( s );
    free( s );
    return( str );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_trim( pointwiseXY_CPy *self ) {

    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( ( n = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ptwXY_trim( smr, n ) != nfu_Okay ) {
        ptwXY_free( n );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static char *pointwiseXY_C_toString_isFormatForDouble( char *format, int returnFormat ) {

    size_t i;
    char *s;

    if( format == NULL ) return( format );
    for( s = format; *s; s++ ) {
        if( *s == '%' ) {
            if( s[1] != '%' ) break;
            s++;
        }
    }
    if( returnFormat ) return( ( *s == 0 ) ? format : NULL );
    if( *s == 0 ) return( NULL );
    s++;
    i = strspn( s, "+-0123456789." );
    s = &(s[i]);
    if( ( *s != 'e' ) && ( *s != 'E' ) && ( *s != 'f' ) && ( *s != 'F' ) && ( *s != 'g' ) && ( *s != 'G' ) ) return( NULL );
    s++;
    return( s );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_union(  pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int fillWithSelf = 0, trim = 0, unionOptions = 0, status;
    static char *kwlist[] = { "other", "fillWithSelf", "trim", NULL };
    PyObject *otherPy;
    pointwiseXY_CPy *other, *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "O|ii", kwlist, &otherPy, &fillWithSelf, &trim ) ) return( NULL );
    other = (pointwiseXY_CPy *) otherPy;

    if( ( status = PyObject_IsInstance( otherPy, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "other instance is not an pointwiseXY_C instance" ) );
    if( pointwiseXY_C_checkStatus2( other, "other" ) != 0 ) return( NULL );

    if( fillWithSelf ) unionOptions |= ptwXY_union_fill;
    if( trim ) unionOptions |= ptwXY_union_trim;

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill || other->infill, self->safeDivide || other->safeDivide ) ) == NULL ) return( NULL );
    if( ( nPy->ptwXY = ptwXY_union( smr, self->ptwXY, other->ptwXY, unionOptions ) ) == NULL ) {
        Py_DECREF( (PyObject *) nPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_scaleOffsetXAndY( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int insitu = 0;
    double xScale = 1., xOffset = 0., yScale = 1., yOffset = 0.;
    pointwiseXY_CPy *nPy = NULL;
    static char *kwlist[] = { "xScale", "xOffset", "yScale", "yOffset", "insitu", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ddddi", kwlist, &xScale, &xOffset, &yScale, &yOffset, &insitu ) ) return( NULL );

    if( insitu == 0 ) {
        if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) return( NULL );

        if( ( nPy->ptwXY = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
            Py_DECREF( (PyObject *) nPy );
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( NULL );
        } }
    else {
        nPy = self;
    }

    if( ptwXY_scaleOffsetXAndY( smr, nPy->ptwXY, xScale, xOffset, yScale, yOffset ) != nfu_Okay ) {
        Py_DECREF( (PyObject *) nPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( insitu ) Py_INCREF( (PyObject *) self );    /* FIXME: Why is this needed when insitu if True? Or is it? */
    return( (PyObject *) nPy );
}
#if 0
/*
************************************************************
*/
static PyObject *pointwiseXY_C_scaleAndOffsetDomainWithPointwiseXYs( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {
/* This does not work. */

    int offsetStatus = 0, slopeStatus = 0, skipLastPoint = 0;
    pointwiseXY_CPy *nPy = NULL;
    PyObject *offsetPy, *slopePy;
    ptwXYPoints *offset_ptwXY = NULL, *slope_ptwXY = NULL;
    static char *kwlist[] = { "offset", "slope", "skipLastPoint", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|OOi", kwlist, &offsetPy, &slopePy, &skipLastPoint ) ) return( NULL );

    if( ( offsetStatus = pyObject_getNumberOrPtwXYAsPtwXY( NULL, self, offsetPy, &offset_ptwXY, "offset" ) ) < 0 ) goto error;
    if( ( slopeStatus  = pyObject_getNumberOrPtwXYAsPtwXY( NULL, self, slopePy,  &slope_ptwXY, "slope"  ) ) < 0 ) goto error;

    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) goto error;

    if( ( nPy->ptwXY = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        Py_DECREF( (PyObject *) nPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        goto error;
    }

error:
    if( offsetStatus > 0 ) ptwXY_free( offset_ptwXY );
    if( slopeStatus > 0 ) ptwXY_free( slope_ptwXY );
    return( NULL );
}
#endif
/*
************************************************************
*/
static PyObject *pointwiseXY_C_mergeClosePoints(  pointwiseXY_CPy *self, PyObject *args ) {

    double epsilon;
    pointwiseXY_CPy *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "d", &epsilon ) ) return( NULL );
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) return( NULL );

    if( ( nPy->ptwXY = ptwXY_clone( smr, self->ptwXY ) ) == NULL ) {
        Py_DECREF( (PyObject *) nPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    if( ptwXY_mergeClosePoints( smr, nPy->ptwXY, epsilon ) != nfu_Okay ) {
        Py_DECREF( (PyObject *) nPy );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_domainGrid( pointwiseXY_CPy *self, PyObject *args ) {

    int i;
    double scale = 1.0;
    ptwXPoints *ptwX;
    PyObject *nPy;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "|d", &scale ) ) return( NULL );
    if( ( ptwX = ptwXY_getXArray( smr, self->ptwXY ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    for( i = 0; i < ptwX_length( NULL, ptwX ); i++ ) ptwX_setPointAtIndex( NULL, ptwX, i, scale * ptwX_getPointAtIndex_Unsafely( ptwX, i ) );
    nPy = pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwX );
    ptwX_free( ptwX );
    return( nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_domain( pointwiseXY_CPy *self ) {

    double domainMin, domainMax;

    if( isOkayAndHasData( self ) == -1 ) return( NULL );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    return( Py_BuildValue( "[dd]", domainMin, domainMax ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_domainMin( pointwiseXY_CPy *self ) {

    double domainMin;

    if( isOkayAndHasData( self ) == -1 ) return( NULL );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, NULL ) == -1 ) return( NULL );

    return( Py_BuildValue( "d", domainMin ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_domainMax( pointwiseXY_CPy *self ) {

    double domainMax;

    if( isOkayAndHasData( self ) == -1 ) return( NULL );
    if( pointwiseXY_C_domainMinMax( self, NULL, &domainMax ) == -1 ) return( NULL );

    return( Py_BuildValue( "d", domainMax ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_domainSlice( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int fill = 1;
    double domainMin, domainMax, dullEps = 0;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "domainMin", "domainMax", "fill", "dullEps", NULL };
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );
    if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) return( NULL );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "|ddid", kwlist, &domainMin, &domainMax, &fill, &dullEps) ) return( NULL );

    if( ( n = ptwXY_domainSlice( smr, self->ptwXY, domainMin, domainMax, 10, fill ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( self->infill, self->safeDivide ) ) == NULL ) {
        ptwXY_free( n ); }
    else {
        nPy->ptwXY = n;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_range( pointwiseXY_CPy *self ) {

    double rangeMin, rangeMax;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( pointwiseXY_C_rangeMinMax( self, &rangeMin, &rangeMax ) == -1 ) return( NULL );
    return( Py_BuildValue( "[dd]", rangeMin, rangeMax ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_rangeMin( pointwiseXY_CPy *self ) {

    double rangeMin;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( pointwiseXY_C_rangeMinMax( self, &rangeMin, NULL ) == -1 ) return( NULL );
    return( Py_BuildValue( "d", rangeMin ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_rangeMax( pointwiseXY_CPy *self ) {

    double rangeMax;

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( pointwiseXY_C_rangeMinMax( self, NULL, &rangeMax ) == -1 ) return( NULL );
    return( Py_BuildValue( "d", rangeMax ) );
}
/*
************************************************************
*/
static int pointwiseXY_C_domainMinMax( pointwiseXY_CPy *self, double *domainMin, double *domainMax ) {

    statusMessageReporting *smr = &(self->smr);

    if( isOkayAndHasData( self ) == -1 ) return( 1 );
    if( domainMin != NULL ) {
        if( ptwXY_domainMin( smr, self->ptwXY, domainMin ) != nfu_Okay ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( -1 );
        }
    }
    if( domainMax != NULL ) {
        if( ptwXY_domainMax( smr, self->ptwXY, domainMax ) != nfu_Okay ) {
            pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
            return( -1 );
        }
    }

    return( 0 );
}
/*
************************************************************
*/
static int pointwiseXY_C_rangeMinMax( pointwiseXY_CPy *self, double *rangeMin, double *rangeMax ) {

    statusMessageReporting *smr = &(self->smr);
    double _rangeMin, _rangeMax;

    if( isOkayAndHasData( self ) == -1 ) return( -1 );

    if( ptwXY_range( smr, self->ptwXY, &_rangeMin, &_rangeMax ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( -1 );
    }
    if( rangeMin != NULL ) *rangeMin = _rangeMin;
    if( rangeMax != NULL ) *rangeMax = _rangeMax;

    return( 0 );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_ysMappedToXs( pointwiseXY_CPy *self, PyObject *args ) {

    int64_t offset;
    ptwXPoints *ptwX, *Xs = NULL;
    PyObject *nPy, *Xs_Py;
    statusMessageReporting *smr = &(self->smr);

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( NULL );

    if( !PyArg_ParseTuple( args, "O", &Xs_Py ) ) return( NULL );

    if( ( Xs = pointwiseXY_C_PyFloatList_to_ptwXPoints( Xs_Py ) ) == NULL ) return( NULL );

    ptwX = ptwXY_ysMappedToXs( smr, self->ptwXY, Xs, &offset );
    ptwX_free( Xs );

    if( ptwX == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
        return( NULL );
    }

    nPy = pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwX );
    ptwX_free( ptwX );

    Xs_Py = Py_BuildValue( "lO", (long) offset, nPy );
    Py_DECREF( nPy );

    return( Xs_Py );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_defaultAccuracy( pointwiseXY_CPy *self, PyObject *args ) {

    return( Py_BuildValue( "d", _defaultAccuracy ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_emptyXYs( void ) {

    ptwXYPoints *ptwXY;
    pointwiseXY_CPy *XYsPy;

    if( ( XYsPy = pointwiseXY_CNewInitialize( 1, 1 ) ) != NULL ) {
        if( ( ptwXY = ptwXY_new( NULL, ptwXY_interpolationLinLin, NULL, 12, 1e-3, 0, 0, 0 ) ) == NULL ) {
            Py_DECREF( (PyObject *) XYsPy );
            PyErr_NoMemory( );
            return( NULL );
        }
        XYsPy->ptwXY = ptwXY;
    }
    return( (PyObject *) XYsPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_createFromFunction( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    PyObject *f_args[3], *xs_Py;
    PyObject *f, *parameters;
    int checkForRoots = 0, infill = 1, safeDivide = 1;
    double accuracy, biSectionMax;
    ptwXPoints *xs = NULL;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "xs", "f", "parameters", "accuracy", "biSectionMax", "checkForRoots", "infill", "safeDivide", NULL };
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "OOOdd|iii", kwlist, &xs_Py, &f, &parameters, &accuracy, &biSectionMax, &checkForRoots,
        &infill, &safeDivide ) ) return( NULL );
    f_args[0] = f;
    f_args[1] = parameters;
    f_args[2] = NULL;

    if( !PyCallable_Check( f ) ) {
        PyErr_SetString( PyExc_TypeError, "Second argument must be callable" );
        return( NULL );
    }

    if( ( xs = pointwiseXY_C_PyFloatList_to_ptwXPoints( xs_Py ) ) == NULL ) return( NULL );

    Py_INCREF( f );
    n = ptwXY_createFromFunction2( &smr, xs, pointwiseXY_C_createFromFunction2, (void *) f_args, accuracy, 
        checkForRoots, (int) biSectionMax );
    Py_DECREF( f );

    ptwX_free( xs );
    if( n == NULL ) {
        if( f_args[2] != 0 ) return( NULL );
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }

    if( ( nPy = pointwiseXY_CNewInitialize( infill, safeDivide ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static nfu_status pointwiseXY_C_createFromFunction2( statusMessageReporting *smr, double x, double *y, void *argList ) {

    PyObject **f_args = (PyObject **) argList, *result;
    nfu_status status_nf = nfu_Okay;

    result = PyObject_CallFunction( (PyObject *) f_args[0], "(d,O)", x, f_args[1] );
    if( result == NULL ) {
        f_args[2] = f_args[0];
        return( nfu_badInput );
    }
    if( pointwiseXY_C_PyNumberToFloat( result, y ) != 0 ) {
        pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert returned value from createFromFunction function to float" );
        status_nf = nfu_badInput;
        f_args[2] = f_args[0];
    }

    Py_DECREF( result );
    return( status_nf );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_createFromString( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    int infill = 1, safeDivide = 1;
    double accuracy, biSectionMax;
    ptwXYPoints *ptwXY = NULL;
    pointwiseXY_CPy *ptwXYPy;
    char *str, *interpolationStr = NULL, *endCharacter, sep = ' ';
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;
    static char *kwlist[] = { "str", "accuracy", "biSectionMax", "interpolation", "infill", "safeDivide", "sep", NULL };
    statusMessageReporting smr;
    PyObject *xys_Py;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "sdd|siiOOc", kwlist, 
        &str, &accuracy, &biSectionMax, &interpolationStr, &infill, &safeDivide, &sep ) ) return( NULL );
    if( pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 1 ) != 0 ) return( NULL );

    if( ( ptwXY = ptwXY_fromString( &smr, str, sep, interpolation, interpolationStr, biSectionMax, accuracy, &endCharacter, 1 ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        goto err;
    }

    if( ( ptwXYPy = pointwiseXY_CNewInitialize( infill, safeDivide ) ) == NULL ) goto err;
    ptwXYPy->ptwXY = ptwXY;

    xys_Py = Py_BuildValue( "Os", ptwXYPy, endCharacter );
    Py_DECREF( ptwXYPy );

    return( xys_Py );

err:
    if( ptwXY != NULL ) ptwXY_free( ptwXY );
    return( NULL );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_ptwXY_interpolatePoint( PyObject *self, PyObject *args ) {

    double x, y, x1, y1, x2, y2;
    ptwXY_interpolation interpolation;
    char *interpolationStr;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    PyArg_ParseTuple( args, "sddddd", &interpolationStr, &x, &x1, &y1, &x2, &y2 );
    if( pointwiseXY_C_checkInterpolationString( interpolationStr, &interpolation, 0 ) != 0 ) return( NULL );
    if( ptwXY_interpolatePoint( &smr, interpolation, x, &y, x1, y1, x2, y2 ) != nfu_Okay ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    return( Py_BuildValue( "d", y ) );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_gaussian( pointwiseXY_CPy *self, PyObject *args, PyObject *keywords ) {

    double accuracy, domainMin, domainMax, offset = 0., sigma = 1., amplitude = 1., dullEps = 0.;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    static char *kwlist[] = { "accuracy", "domainMin", "domainMax", "offset", "sigma", "amplitude", "dullEps", NULL };
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "ddd|dddd", kwlist, &accuracy, &domainMin, &domainMax, &offset, &sigma, 
        &amplitude, &dullEps) ) return( NULL );

    if( ( n = ptwXY_createGaussian( &smr, accuracy, offset, sigma, amplitude, domainMin, domainMax, dullEps ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_basicGaussian( pointwiseXY_CPy *self, PyObject *args ) {

    double accuracy;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTuple( args, "d", &accuracy ) ) return( NULL );

    if( ( n = ptwXY_createGaussianCenteredSigma1( &smr, accuracy ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_unitbaseInterpolate( pointwiseXY_CPy *self, PyObject *args ) {

    double w, lw, uw;
    ptwXYPoints *n;
    pointwiseXY_CPy *nPy, *lXY, *uXY;
    int status, scaleRange;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( !PyArg_ParseTuple( args, "ddOdOi", &w, &lw, &lXY, &uw, &uXY, &scaleRange ) ) return( NULL );

    if( ( status = PyObject_IsInstance( (PyObject* ) lXY, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "lowerXYs instance is not an pointwiseXY_C instance" ) );
    if( pointwiseXY_C_checkStatus2( lXY, "lowerXYs" ) != 0 ) return( NULL );

    if( ( status = PyObject_IsInstance( (PyObject* ) uXY, (PyObject* ) &pointwiseXY_CPyType ) ) < 0 ) return( NULL );
    if( status == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnNull( "upperXYs instance is not an pointwiseXY_C instance" ) );
    if( pointwiseXY_C_checkStatus2( uXY, "upperXYs" ) != 0 ) return( NULL );

    if( ( n = ptwXY_unitbaseInterpolate( &smr, w, lw, lXY->ptwXY, uw, uXY->ptwXY, scaleRange ) ) == NULL ) {
        pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, &smr );
        return( NULL );
    }
    if( ( nPy = pointwiseXY_CNewInitialize( 1, 1 ) ) == NULL ) {
        ptwXY_free( n );
        return( NULL );
    }
    nPy->ptwXY = n;
    return( (PyObject *) nPy );
}
/*
************************************************************
*/
static PyObject *floatToShortestString_C( PyObject *self, PyObject *args, PyObject *keywords ) {

    int significantDigits = 15, trimZeros = 1, keepPeriod = 0, favorEFormBy = 0, includeSign = 0, flags = 0;
    double value;
    char *string;
    PyObject *StringPy;
    static char *kwlist[] = { "value", "significantDigits", "trimZeros", "keepPeriod", "favorEFormBy", "includeSign", NULL };

    if( !PyArg_ParseTupleAndKeywords( args, keywords, "d|iiiii", kwlist, &value, &significantDigits, &trimZeros, &keepPeriod,
            &favorEFormBy, &includeSign ) ) return( NULL );

    if( trimZeros ) flags += nf_floatToShortestString_trimZeros;
    if( keepPeriod ) flags += nf_floatToShortestString_keepPeriod;
    if( includeSign ) flags += nf_floatToShortestString_includeSign;

    string = nf_floatToShortestString( value, significantDigits, favorEFormBy, flags );
    StringPy = Py_BuildValue( "s", string );
    free( string );
    return( StringPy );
}
/*
************************************************************
*/
static int pointwiseXY_C_PythonXYPairToCPair( PyObject *XYPairPy, double *x, double *y, int64_t index ) {
/*
*   Status is 0) if all is ok.
*             1) if XYPairPy does not support iteration
*             2) if the list is empty (i.e., no x object).
*             3) if the list contains only one entry (i.e., no y object).
*             4) if the list contains more than 2 entries.
*             5) if x entry could not be converted to a float.
*             6) if y entry could not be converted to a float.
*/
    int status = 1;
    PyObject *iterator = PyObject_GetIter( XYPairPy ), *xPy = NULL, *yPy = NULL, *theEnd = NULL;

    if( iterator != NULL ) {
        status = 2;
        xPy = PyIter_Next( iterator );
    }
    if( xPy != NULL ) {
        status = 3;
        yPy = PyIter_Next( iterator );
        if( yPy != NULL ) {
            status = 4;
            theEnd = PyIter_Next( iterator );
        }
    }
    if( ( xPy != NULL ) && ( yPy != NULL ) && ( theEnd == NULL ) ) {
        status = 0;
        if( pointwiseXY_C_PyNumberToFloat( xPy, x ) != 0 ) status = 5;
        if( ( status == 0 ) && ( pointwiseXY_C_PyNumberToFloat( yPy, y ) != 0 ) ) status = 6;
    }
    if( iterator != NULL ) { Py_DECREF( iterator ); }
    if( xPy != NULL ) { Py_DECREF( xPy ); }
    if( yPy != NULL ) { Py_DECREF( yPy ); }
    if( theEnd != NULL ) { Py_DECREF( theEnd ); }
    if( status != 0 ) {
        switch( status ) {
        case 1 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld does not support iteration", index );
            break;
        case 2 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld is empty", index );
            break;
        case 3 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld only has one entry", index );
            break;
        case 4 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "element at index = %ld has more than two entries", index );
            break;
        case 5 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "x value for element at index = %ld is not a number", index );
            break;
        case 6 :
            status = pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "y value for element at index = %ld is not a number", index );
            break;
        }
    }
    return( status );
}
/*
************************************************************
*/
static int64_t pointwiseXY_C_PythonXYListToCList( PyObject *XYListPy, double **xys ) {

    int status = 0;
    int64_t i = 0, length, n;
    double *d;
    PyObject *item, *iterator;
    pointwiseXY_CPy *ptwXY = (pointwiseXY_CPy *) XYListPy;
    statusMessageReporting *smr = &(ptwXY->smr);

    *xys = NULL;
    if( PyObject_TypeCheck( XYListPy, &pointwiseXY_CPyType ) ) {
        length = ptwXY->ptwXY->length;
        if( length != (int64_t) 0 ) {
            if( ( *xys = (double *) malloc( 2 * (size_t) length * sizeof( double ) ) ) == NULL ) {
                PyErr_NoMemory( );
                return( -1 );
            }
            if( ptwXY_copyToC_XY( smr, ptwXY->ptwXY, 0, length, length, &n, *xys ) != nfu_Okay ) {
                free( *xys );
                pointwiseXY_C_SetPyErrorExceptionFromSMR( PyExc_Exception, smr );
                return( -1 );
            }
        } }
    else {
        if( ( iterator = PyObject_GetIter( XYListPy ) ) == NULL ) return( -1 );
        length = (int64_t) PySequence_Size( XYListPy );
        if( length != (int64_t) 0 ) {
            *xys = d = (double *) malloc( 2 * (size_t) length * sizeof( double ) );
            if( d == NULL ) {
                PyErr_NoMemory( );
                Py_DECREF( iterator );
                return( -1 );
            }

            while( ( status == 0 ) && ( ( item = PyIter_Next( iterator ) ) != NULL ) ) {
                status = pointwiseXY_C_PythonXYPairToCPair( item, d, &(d[1]), i );
                Py_DECREF( item );
                i++;
                d += 2;
            }

            if( status != 0 ) {
                free( *xys );
                length = -1;
            }
        }
        Py_DECREF( iterator );
    }

    return( length );
}
/*
************************************************************
*/
static int64_t pointwiseXY_C_pythonDoubleListToCList( PyObject *PyDoubleList, double **ds, int ascending ) {

    int status = 0;
    int64_t length, i;
    double *d, dPrior = 0.;
    PyObject *item, *iterator;

    *ds = NULL;
    if( ( iterator = PyObject_GetIter( PyDoubleList ) ) == NULL ) return( -1 );
    if( ( length = (int64_t) PySequence_Size( PyDoubleList ) ) != (int64_t) 0 ) {
        if( ( *ds = (double *) malloc( (size_t) length * sizeof( double ) ) ) == NULL ) {
            PyErr_NoMemory( );
            Py_DECREF( iterator );
            return( -1 );
        }
        for( i = 0, d = *ds, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i++, d++ ) {
            if( ( status = pointwiseXY_C_PyNumberToFloat( item, d ) ) != 0 )
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i );
            if( ( status == 0 ) && ascending && ( i != 0 ) ) {
                if( *d <= dPrior ) {
                    status = -1;
                    pointwiseXY_C_SetPyErrorExceptionReturnNull( "data not in ascending order at index %d and %d", (int) i - 1, (int) i );
                }
            }
            dPrior = *d;
            Py_DECREF( item );
        }
    }
    Py_DECREF( iterator );
    if( status != 0 ) {
        length = -1;
        free( *ds );
        *ds = NULL;
    }

    return( length );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_ptwXPoints_to_PyFloatList( ptwXPoints *ptwX ) {

    int64_t i;
    PyObject *l = PyList_New( 0 ), *ni;

    if( l == NULL ) return( NULL );
    for( i = 0; i < ptwX->length; i++ ) {
        ni = Py_BuildValue( "d", ptwX->points[i] );
        if( pointwiseXY_C_addedItemToPythonList( l, ni ) != 0 ) return( NULL );
    }
    return( l );
}
/*
************************************************************
*/
static ptwXPoints *pointwiseXY_C_PyFloatList_to_ptwXPoints( PyObject *PyFloatList ) {

    int status = 0;
    int64_t i, length;
    double d;
    PyObject *item, *iterator;
    ptwXPoints *ptwX;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    if( ( iterator = PyObject_GetIter( PyFloatList ) ) == NULL ) return( NULL );
    length = (int64_t) PySequence_Size( PyFloatList );
    if( ( ptwX = ptwX_new( &smr, length ) ) == NULL ) {
        PyErr_NoMemory( ); }
    else {
        ptwX->length = length;
        for( i = 0, item = PyIter_Next( iterator ); ( status == 0 ) && ( item != NULL ); item = PyIter_Next( iterator ), i++ ) {
            if( ( status = pointwiseXY_C_PyNumberToFloat( item, &d ) ) == 0 ) {
                ptwX->points[i] = d; }
            else {
                pointwiseXY_C_SetPyErrorExceptionReturnNull( "could not convert item at index = %d to float", i );
            }
            Py_DECREF( item );
        }
    }
    Py_DECREF( iterator );
    if( status != 0 ) ptwX = ptwX_free( ptwX );

    return( ptwX );
}
/*
************************************************************
*/
static int pointwiseXY_C_PyNumberToFloat( PyObject *n, double *d ) {

    if( PyFloat_Check( n ) ) {
        *d = PyFloat_AsDouble( n ); }
#if PY_MAJOR_VERSION < 3
    else if( PyInt_Check( n ) ) {
        *d = (double) PyInt_AsLong( n ); }
#endif
    else if( PyLong_Check( n ) ) {
        *d = (double) PyLong_AsLongLong( n ); }
    else {
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
static int pyObject_NumberOrPtwXY(  PyObject *other ) {

    if( PyLong_Check( other ) ) {
        return( pyObject_Number ); }
#if PY_MAJOR_VERSION < 3
    else if( PyInt_Check( other ) ) {
        return( pyObject_Number ); }
#endif
    else if( PyFloat_Check( other ) ) {
        return( pyObject_Number ); }
    else if( PyList_Check( other ) ) {
        return( pyObject_List ); }
    else if( PyTuple_Check( other ) ) {
        return( pyObject_Tuple ); }
    else if( PyObject_TypeCheck( other, &pointwiseXY_CPyType ) ) {
        return( pyObject_ptwXY );
    }
    return( pyObject_Unsupported );
}

#if 0
/*
************************************************************
*/
static int pyObject_getNumberOrPtwXYAsPtwXY( statusMessageReporting *smr, pointwiseXY_CPy *self, PyObject *object, ptwXYPoints **ptwXY, char *variableName ) {
/*
    If return value is 
            -1 - an error occured and nothing created by this function that needs to be freed,
             0 - object is a ptwXYPoints instance and nothing created by this function that needs to be freed,
             1 - object is a number and a ptwXYPoints instance is created that the caller needs to free.
*/

    int status = -1;
    double domainMin, domainMax, xys[4];

    *ptwXY = NULL;

    switch( pyObject_NumberOrPtwXY( object ) ) {
    case pyObject_Number :
        if( pointwiseXY_C_domainMinMax( self, &domainMin, &domainMax ) == -1 ) {
            return( -1 );
        }
        xys[0] = domainMin;
        xys[2] = domainMax;

        if( pointwiseXY_C_PyNumberToFloat( object, &xys[1] ) < 0 ) return( -1 );
        xys[3] = xys[1];

        *ptwXY = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 12, 1e-3, 4, 0, 4, xys, 0 );
        status = 1;
        break;
    case pyObject_ptwXY :
        *ptwXY = ((pointwiseXY_CPy *)object)->ptwXY;
        status = 0;
        break;
    default :
        pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "%s is not a number or a pointwiseXY_C.", variableName );
    }

    return( status );
}
#endif

/*
************************************************************
*/
static int pointwiseXY_C_Get_pointwiseXY_CAsSelf( PyObject *self, PyObject *other, PyObject **self2, PyObject **other2,
    pointwiseXY_CPy **self3, pointwiseXY_CPy **other3 ) {
/*
    BRB FIXME. Should PyObject_TypeCheck test use PyObject_IsInstance.
*/
    *self2 = self;
    *other2 = other;
    if( !PyObject_TypeCheck( self, &pointwiseXY_CPyType ) ) {
        *self2 = other;
        *other2 = self;
    }

    *self3 = (pointwiseXY_CPy *) *self2;
    *other3 = (pointwiseXY_CPy *) *other2;

    if( pointwiseXY_C_checkStatus( *self3 ) != 0 ) return( -1 );
    if( PyObject_TypeCheck( *other2, &pointwiseXY_CPyType ) ) {
        if( pointwiseXY_C_checkStatus2( *other3, "other" ) != 0 ) return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
static int pointwiseXY_C_addedItemToPythonList( PyObject *list, PyObject *item ) {

    if( item == NULL ) {
        Py_DECREF( list );
        return( -1 ); }
    else {
        if( PyList_Append( list, item ) != 0 ) {
            Py_DECREF( item );          /* append created a new ref */
            Py_DECREF( list );
            return( -1 );
        }
        Py_DECREF( item );              /* append created a new ref */
    }
    return( 0 );
}
/*
************************************************************
*/
static int isOkayAndHasData( pointwiseXY_CPy *self ) {

    if( pointwiseXY_C_checkStatus( self ) != 0 ) return( -1 );
    if( self->ptwXY->length == 0 ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "pointwiseXY_C object has no data" ) );
    return( 0 );
}
/*
************************************************************
*/
static int pointwiseXY_C_checkInterpolationString( char const *interpolationStr, ptwXY_interpolation *interpolation, int allowOther ) {

    if( interpolationStr == NULL ) interpolationStr = "lin-lin";
    if( strcmp( interpolationStr, "lin-lin" ) == 0 ) {
        *interpolation = ptwXY_interpolationLinLin; }
    else if( strcmp( interpolationStr, "log-lin" ) == 0 ) {
        *interpolation = ptwXY_interpolationLogLin; }
    else if( strcmp( interpolationStr, "lin-log" ) == 0 ) {
        *interpolation =  ptwXY_interpolationLinLog; }
    else if( strcmp( interpolationStr, "log-log" ) == 0 ) {
        *interpolation = ptwXY_interpolationLogLog; }
    else if( strcmp( interpolationStr, "flat" ) == 0 ) {
        *interpolation = ptwXY_interpolationFlat; }
    else {
        *interpolation = ptwXY_interpolationOther;
        if( !allowOther ) return( pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( "interpolation = '%s' not allowed", interpolationStr ) );
    }
    return( 0 );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_GetNone( void ) {

    Py_INCREF( Py_None );
    return( Py_None );
}
/*
************************************************************
*/
static void pointwiseXY_C_SetPyErrorExceptionFromSMR( PyObject *type, statusMessageReporting *smr ) {

    if( smr_isOk( smr ) ) return;               /* This and the next line are probably the same. But will check anyway. */
    if( PyErr_Occurred() == NULL ) {            /* Do not set if exception if on is already set. */
        PyErr_SetString( type, smr_getMessage( smr_firstReport( smr ) ) );
    }
    smr_release( smr );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_SetPyErrorExceptionReturnNull( const char *s, ... ) {

    va_list args;
    char Str[1024];

    va_start( args, s );
    vsnprintf( Str, sizeof( Str ), s, args );
    Str[sizeof( Str ) - 1] = 0;
    PyErr_SetString( PyExc_Exception, Str );
    va_end( args );
    return( NULL );
}
/*
************************************************************
*/
static int pointwiseXY_C_SetPyErrorExceptionReturnMinusOne( const char *s, ... ) {

    va_list args;
    char Str[1024];

    va_start( args, s );
    vsnprintf( Str, sizeof( Str ), s, args );
    Str[sizeof( Str ) - 1] = 0;
    PyErr_SetString( PyExc_Exception, Str );
    va_end( args );
    return( -1 );
}
/*
************************************************************
*/
static int pointwiseXY_C_checkStatus( pointwiseXY_CPy *self ) {

    return( pointwiseXY_C_checkStatus2( self, "self" ) );
}
/*
************************************************************
*/
static int pointwiseXY_C_checkStatus2( pointwiseXY_CPy *self, char const *name ) {

    nfu_status status_nf = self->ptwXY->status;

    if( status_nf == nfu_Okay ) return( 0 );
    pointwiseXY_C_SetPyErrorExceptionReturnMinusOne(
            "pointwiseXY_C object '%s' had a prior error and is not usable: status = %d", name, status_nf,
            nfu_statusMessage( status_nf ) );
    return( -1 );
}
/*
************************************************************
*/
static PyMethodDef pointwiseXY_CPyMethods[] = {

    { "allocatedSize", (PyCFunction) pointwiseXY_C_allocatedSize, METH_NOARGS, "Returns the size of memory allocated in the points region." },
    { "applyFunction", (PyCFunction) pointwiseXY_C_applyFunction, METH_VARARGS | METH_KEYWORDS, \
        "Returns a new instance which is y(x) = f(y_s(x)) where f is a function of one variable and y_s(x) is self's.\n" \
        "y-value at x. The function f must take two arguments. The first is y_s(x) and the second (not including self)\n" \
        "is the second argument to applyFunction. The returned object may contain more points then self as applyFunction\n" \
        "does infilling. That is, points are added recursively between points until accuracy or biSectionMax is reached.\n" \
        "\nArguments are:\n" \
        "   f             a function that returns the new y(x) given the old y(x),\n" \
        "   parameters    any python object. Passed as the second argument to f,\n" \
        "   accuracy      [o] The accuracy for infilling. If not present, taken from self,\n" \
        "   biSectionMax  [o] The maximum number of bi-sections for infilling. If not present, taken from self.\n" \
        "   checkForRoots [o] If true, and biSectionMax > 0, an addition point is added whenever two consecutive points cross the y-axis\n" \
        "                     at (or close to) the root (i.e., the y crossing point).\n" \
        "\nTo return f(x), set self's y-values to x (i.e., y(x) = x for self)." },
    { "changeInterpolation", (PyCFunction) pointwiseXY_C_changeInterpolation, METH_VARARGS | METH_KEYWORDS,
        "Returns a new instance that is equivalent to self, but with the interpolation changed to interpolation.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   interpolation  [o] the new interpolation (default is 'lin-lin'),\n" \
        "   accuracy       [o] the accuracy of the conversion from the old to the new interpolation (default is self's accuracy),\n" \
        "   lowerEps       [o] for flat to linear, the lower eps,\n" \
        "   upperEps       [o] the flat to linear, the upper eps. lowerEps and upperEps cannot both be 0.\n" },
    { "changeInterpolationIfNeeded", (PyCFunction) pointwiseXY_C_changeInterpolationIfNeeded, METH_VARARGS | METH_KEYWORDS,
        "Returns self if it is one of the allowed interpolations. Otherwise, a new instance is returned that is equivalent to self,\n" \
        "but with the interpolation changed to the interpolation of the first allowed interpolations.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   allowedInterpolations the list of allowed interpolations,\n" \
        "   accuracy              [o] the accuracy of the conversion from the old to the new interpolation (default is self's accuracy),\n" \
        "   lowerEps              [o] for flat to linear, the lower eps,\n" \
        "   upperEps              [o] the flat to linear, the upper eps. lowerEps and upperEps cannot both be 0.\n" },
    { "clip", (PyCFunction) pointwiseXY_C_clip, METH_VARARGS | METH_KEYWORDS,
        "Returns a new instance which is the same as self, but with the y-values clipped between rangeMin and rangeMax\n" \
        "clip may add points, to insure that the return instance has the same shape as self between rangeMin and rangeMax\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   rangeMin  [o]   the lower y-value for clipping,\n" \
        "   rangeMax  [o]   the upper y-value for clipping.\n" },
    { "coalescePoints", (PyCFunction) pointwiseXY_C_coalescePoints, METH_VARARGS, "Moves all points in overflow region to points region." },
    { "convolute", (PyCFunction) pointwiseXY_C_convolute, METH_VARARGS | METH_KEYWORDS, \
        "Returns a new instance which is the convolution of self with the first argument, that must also be a pointwiseXY_C instance." },
    { "copy", (PyCFunction) pointwiseXY_C_copy, METH_NOARGS, "Returns a copy of self." },
    { "cloneToInterpolation", (PyCFunction) pointwiseXY_C_cloneToInterpolation, METH_VARARGS | METH_KEYWORDS, 
        "A clone of self is returned with its interpolation changed, but no points altered or added." \
        "\nArguments are:\n" \
        "   interpolation   the interpolation for the cloned instance." },
    { "copyDataToXYs", (PyCFunction) pointwiseXY_C_copyDataToXYs, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of self's data as [ x1, y1 ], ... [ xn, yn ] ].\n" \
        "\nArguments are:\n" \
        "   xScale [o]  a float to scale all x-values by,\n" \
        "   yScale [o]  a float to scale all y-values by." },
    { "copyDataToXsAndYs", (PyCFunction) pointwiseXY_C_copyDataToXsAndYs, METH_VARARGS | METH_KEYWORDS, 
        "Returns python list of length two. The first element is a python list of\n" \
        "self's x-values and the second element is a python list of self's y-values.\n" \
        "\nArguments are:\n" \
        "   xScale [o]  a float to scale all x-values by,\n" \
        "   yScale [o]  a float to scale all y-values by." },
    { "dullEdges", (PyCFunction) pointwiseXY_C_dullEdges, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new instance that is a copy of self, except the endpoints are guaranteed to have 0's for y-values\n" \
        "as long as lowerEps and/or upperEps are none zero (see below). The following algorithm is used.\n\n" \
        "If lowerEps is zero or the y-value at domainMin is 0, then nothing is added at the lower end. Otherwise, the following is\n"\
        "done at the lower end.\n" \
        "   If lowerEps is positive, a point 'abs( domainMin * lowerEps )' above domainMin is added with its interpolated y-value; provided,\n" \
        "       this x value is 'abs( domainMin * lowerEps )' less than the next x-value. In the prior sentence, if domainMin is 0, then\n" \
        "       replace 'abs( domainMin * lowerEps )' with 'abs( lowerEps )'. Independent of whether a point is added, domainMin's y-value is\n" \
        "       set to 0.\n" \
        "   If lowerEps is negative, the logic for adding the point above domainMin for positive lowerEps is followed. In addition, a\n" \
        "       point is added 'abs( domainMin * lowerEps )' below domainMin with a value of zero and the point at domainMin is reset by\n" \
        "       interpolating the new surrounding values. However, if positiveXOnly is True and the point below domainMin would cause\n" \
        "       a negative x-value (and domainMin is not negative) then the logic for positive lowerEps is implemented instead.\n" \
        "\n" \
        "   The logic for upperEps is similar to lowerEps except, replace domainMin with domainMax, below with above and above with below.\n" \
        "       Also positiveXOnly is ignored.\n" \
        "\nArguments are:\n" \
        "   lowerEps       [o] a point (or two if lowerEps is negative) is (are) added a distance domainMin * lowerEps from domainMin,\n" \
        "   upperEps       [o] a point (or two if upperEps is negative) is (are) added a distance domainMax * upperEps from domainMax,\n" \
        "   positiveXOnly  [o] this only applies to lowerEps and only if an added point would be negative when domainMin is non-negative." },
    { "exp", (PyCFunction) pointwiseXY_C_exp, METH_VARARGS,
        "self.exp( a )\n\n" \
        "Returns a new instance with ts y-values set to exp( a * [self's y-values] ). x-values are added to meet required accuaracy."},
    { "getAccuracy", (PyCFunction) pointwiseXY_C_getAccuracy, METH_NOARGS, "Returns self's accuracy value." },
    { "getBiSectionMax", (PyCFunction) pointwiseXY_C_getBiSectionMax, METH_NOARGS, "Returns self's biSectionMax value." },
    { "getInfill", (PyCFunction) pointwiseXY_C_getInfill, METH_NOARGS, "Returns self's infill flag." },
    { "getInterpolation", (PyCFunction) pointwiseXY_C_getInterpolation, METH_NOARGS, "Returns self's interpolation as a string." },
    { "setInterpolation", (PyCFunction) pointwiseXY_C_setInterpolation, METH_VARARGS, "Set self's interpolation. Other interpolation not supported." },
    { "getNFStatus", (PyCFunction) pointwiseXY_C_getNFStatus, METH_NOARGS, "Returns self's numerical functions status flag (an integer)." },
    { "lowerIndexBoundingX", (PyCFunction) pointwiseXY_C_lowerIndexBoundingX, METH_VARARGS,
        "Returns the lower index in self that bounds the x-value. -1 is returned if x is outside domain.\n" \
        "\nArguments are:\n" \
        "   x       the x-value whose lower bounding index is return." },
    { "isInterpolationLinear", (PyCFunction) pointwiseXY_C_isInterpolationLinear, METH_NOARGS,
        "Returns True if self's interpolation is 'lin-lin', otherwise returns False." },
    { "isInterpolationOther", (PyCFunction) pointwiseXY_C_isInterpolationOther, METH_NOARGS,
        "Returns True if self's interpolation is 'other', otherwise returns False." },
    { "getslice", (PyCFunction) pointwiseXY_C_getslice, METH_VARARGS | METH_KEYWORDS,
        "getslice( lowerIndex, upperIndex )\n\n" \
        "Returns the slice of self from lowerIndex to upperIndex with a step of 1. This is redundent with __getslice__ and will be deprecated.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   lowerIndex [o]      Value of the lower index. Default is 0.\n" \
        "   upperIndex [o]      Value of the upper index. Default is length of pointwiseXY_C.\n" \
        "\nExample: Return a pointwise representation of 'x * sin( x**2 )' in the domain [1, 10]." },
    { "getSecondaryCacheSize", (PyCFunction) pointwiseXY_C_getSecondaryCacheSize, METH_NOARGS, "Returns the size of self's secondary cache." },
    { "getSafeDivide", (PyCFunction) pointwiseXY_C_getSafeDivide, METH_NOARGS, "Returns self's safeDivide flag." },
    { "evaluate", (PyCFunction) pointwiseXY_C_evaluate, METH_VARARGS, "Gets the y value at x." },
    { "getUserFlag", (PyCFunction) pointwiseXY_C_getUserFlag, METH_NOARGS, "Gets the user flag." },
    { "groupOneFunction", (PyCFunction) pointwiseXY_C_groupOneFunction, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of self between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "groupTwoFunctions", (PyCFunction) pointwiseXY_C_groupTwoFunctions, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of the product of self and f2 between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   f2                  the second pointwiseXY_C function with the integrand being the product of self * f2\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "groupThreeFunctions", (PyCFunction) pointwiseXY_C_groupThreeFunctions, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of the product of self, f2 and f3 between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   f2                  the second pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   f3                  the third pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "groupFourFunctions", (PyCFunction) pointwiseXY_C_groupFourFunctions, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list of float values. Each value is the integral of the product of self, f2, f3 and f4 between two consecutive group boundaries.\n" \
        "\nArguments are:\n" \
        "   f2                  the second pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   f3                  the third pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   f4                  the fourth pointwiseXY_C function with the integrand being the product of self * f2 * f3\n" \
        "   groupBoundaries     the list of group boundaries,\n" \
        "   norm                each value returned can be normalized as directed by one of the following allowed values\n" \
        "       None                no normalization is applied,\n" \
        "       'dx'                each value is normalized by the width of its interval,\n" \
        "       list                a list of floats, one each for each group which the group is normalized by." },
    { "integrate", (PyCFunction) pointwiseXY_C_integrate, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self from domainMin to domainMax.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   domainMin    [o] the lower limit of the integral, default is domainMin of self,\n" \
        "   domainMax    [o] the upper limit of the integral, default is domainMax of self." },
    { "integrateWithFunction", (PyCFunction) pointwiseXY_C_integrateWithFunction, METH_VARARGS | METH_KEYWORDS, 
        "Returns a float value that is the integral of self times f(x) for x from domainMin to domainMax (i.e. integral dx self(x) * f(x))." \
        " This method uses an adaptive method with a Gauss-Legendre quadrature to calculate the integral to tolerance." \
        " The first argument, func, must return the value of f(x) at x. func is called as func( x, parameters ) where" \
        " parameters is the third argument to this method. The Gauss-Legendre quadrature is exact for any polynomial" \
        " of degree n for which n < degree where degree is one of the arguments. If recursionLimit is needed by the" \
        " adaptive method, tolerance will, most likely, not be met.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   func                a python function taking x and parameters as arguments and returning the value f(x),\n" \
        "   tolerance           the desired tolerance for the integral,\n" \
        "   parameters      [o] a python object passed passed to func [defalut=None],\n" \
        "   domainMin            [o] the lower limit of the integral [default is domainMin of self],\n" \
        "   domainMax            [o] the upper limit of the integral [default is domainMax of self],\n" \
        "   degree          [o] the polynomial degree for which the Gauss-Legendre quadrature is exact [default=4],\n" \
        "   recursionLimit  [o] the maximum recursion depth used in the adaptive quadrature [default=10]." },
    { "integrateWithWeight_x", (PyCFunction) pointwiseXY_C_integrateWithWeight_x, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self weighted by x from domainMin to domainMax (i.e. integral dx x * self(x)).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   domainMin    [o] the lower limit of the integral, default is domainMin of self,\n" \
        "   domainMax    [o] the upper limit of the integral, default is domainMax of self." },
    { "integrateWithWeight_sqrt_x", (PyCFunction) pointwiseXY_C_integrateWithWeight_sqrt_x, METH_VARARGS | METH_KEYWORDS, 
        "Returns float a value that is the integral of self weighted by sqrt( x ) from domainMin to domainMax (i.e. integral dx sqrt( x ) * self(x)).\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   domainMin    [o] the lower limit of the integral, default is domainMin of self,\n" \
        "   domainMax    [o] the upper limit of the integral, default is domainMax of self." },
    { "inverse", (PyCFunction) pointwiseXY_C_inverse, METH_NOARGS, 
        "Returns a new instance with the x- and y-values of self being the y- and x-values of the returned instance, respectively." },
    { "normalize", (PyCFunction) pointwiseXY_C_normalize, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new instance with the same x-values as self but with the y-values scaled so that the area of the curve is 1.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   insitu  [o] If True, self is normalize and returned. Otherwise, a new instances is created and normalize and self is unchanged." },
    { "pop", (PyCFunction) pointwiseXY_C_pop, METH_VARARGS | METH_KEYWORDS, \
        "Removes and returns point at index (default last).\n" },
    { "runningIntegral", (PyCFunction) pointwiseXY_C_runningIntegral, METH_NOARGS, 
        "Returns, as a python list, the integrals between successive x-values." },
    { "areDomainsMutual", (PyCFunction) pointwiseXY_C_areDomainsMutual, METH_VARARGS, 
        "This routine returns True if self and the first argument have a mutual domain, and false otherwise." },
    { "mutualify", (PyCFunction) pointwiseXY_C_mutualify, METH_VARARGS | METH_KEYWORDS, 
        "Returns a python list containing two pointwiseXY_C instances that are the mutualification of self and other.\n" \
        "Mutualification is the act of modifying, if needed, two pointwiseXY_C instances so that their domains are mutual.\n" \
        "Self and other are not altered.\n" \
        "\nArguments are:\n" \
        "   lowerEps1       the lowerEps applied to self if needed, see dullEdges for meaning,\n" \
        "   upperEps1       the upperEps applied to self if needed, see dullEdges for meaning,\n" \
        "   positiveXOnly1  the positiveXOnly applied to self if needed, see dullEdges for meaning,\n" \
        "   other           the second pointwiseXY_C instance to mutualify self with,\n" \
        "   lowerEps2       the lowerEps applied to other if needed, see dullEdges for meaning,\n" \
        "   upperEps2       the upperEps applied to other if needed, see dullEdges for meaning,\n" \
        "   positiveXOnly2  the positiveXOnly applied to other if needed, see dullEdges for meaning,\n" },
    { "overflowAllocatedSize", (PyCFunction) pointwiseXY_C_overflowAllocatedSize, METH_NOARGS, 
        "Returns the size of memory allocated for the overflow region." },
    { "overflowLength", (PyCFunction) pointwiseXY_C_overflowLength, METH_NOARGS, "Returns of number of points in the overflow region." },
    { "plot", (PyCFunction) pointwiseXY_C_plot, METH_NOARGS, 
        "Calls Gnuplot, if it exists, to plot self. This is a simple method, mainly for debugging.\n" \
        "\nArguments: None." },
    { "reallocatePoints", (PyCFunction) pointwiseXY_C_reallocatePoints, METH_VARARGS | METH_KEYWORDS, \
        "self.reallocatePoints( size, forceSmaller = True )\n\n" \
        "Adjusts the memory allocated for primary points to the maximum of size and the current length of self.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   size              the desired allocated size of self (actual size will be larger if length is greater than size),\n" \
        "   forceSmaller  [o] if False action is only taken if allocated is significantly greater than size (default is True).\n" },
    { "reallocateOverflowPoints", (PyCFunction) pointwiseXY_C_reallocateOverflowPoints, METH_VARARGS, \
        "self.reallocateOverflowPoints( size )\n\n" \
        "Adjusts the memory allocated for the overflow points to size (the first and only argument)." },
    { "setData", (PyCFunction) pointwiseXY_C_setData, METH_VARARGS, "Replaces the data in self with the frist argument. This argument must be a list. " \
        "Each item of the list must contain two floats, or objects that can be convert to float (e.g., [ [ 1, 2 ], [ 2, 4 ], [ 4, 0.5 ] ]" },
    { "setAccuracy", (PyCFunction) pointwiseXY_C_setAccuracy, METH_VARARGS, "Sets self's accuracy value and returns the actual value set." },
    { "setBiSectionMax", (PyCFunction) pointwiseXY_C_setBiSectionMax, METH_VARARGS, "Sets self's biSectionMax value and returns the actual value set." },
    { "setSecondaryCacheSize", (PyCFunction) pointwiseXY_C_setSecondaryCacheSize, METH_VARARGS, "Sets the size of self's secondary cache to size." },
    { "setDataFromList", (PyCFunction) pointwiseXY_C_setDataFromList, METH_VARARGS, 
        "Replaces the data in self with the python list given by the first argument.\n" \
        "\nArguments are:\n" \
        "   xys      the list of xy-values as a single python list (e.g., [ 1, 2, 4, 3, 5, 6]\n" },
    { "setDataFromXsAndYs", (PyCFunction) pointwiseXY_C_setDataFromXsAndYs, METH_VARARGS, 
        "Replaces the data in self with the 2 python lists given by the first 2 arguments.\n" \
        "\nArguments are:\n" \
        "   xs      the list of x-values\n" \
        "   ys      the list of y-values\n" },
    { "setInfill", (PyCFunction) pointwiseXY_C_setInfill, METH_VARARGS, "Sets self's infill flag to the first argument." },
    { "setSafeDivide", (PyCFunction) pointwiseXY_C_setSafeDivide, METH_VARARGS, "Sets self's safeDivide flag to the first argument." },
    { "setValue", (PyCFunction) pointwiseXY_C_setValue, METH_VARARGS, "Sets the y value at x." },
    { "setUserFlag", (PyCFunction) pointwiseXY_C_setUserFlag, METH_VARARGS, "Sets the users flag via the first argument." },
    { "showInteralStructure", (PyCFunction) pointwiseXY_C_showInteralStructure, METH_VARARGS | METH_KEYWORDS, 
        "For debbuging only. Dumps information on internal data." },
    { "thicken", (PyCFunction) pointwiseXY_C_thicken, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new instance with denser points than self filled in using self's interpolation. The number of points added are determined by the following arguments\n" \
            "sectionSubdivideMax    maximum number of points to insert between consecutive points (default 1),\n" \
            "dDomainMax             minimum dx step (default 0),\n" \
            "fDomainMax             minimum fractional step (default 1)." },
    { "thin", (PyCFunction) pointwiseXY_C_thin, METH_VARARGS, 
        "Returns a new instance of self thinned to accuracy (i.e., points are removed if possible).\n" \
        "\nArguments are:\n" \
        "   accuracy        The accuracy to maintain for the function when thinning.\n" },
    { "thinDomain", (PyCFunction) pointwiseXY_C_thinDomain, METH_VARARGS, 
        "Returns a new instance of self with fractional x-spacing between points of at least\n" \
        "epsilon. Points are removed if needed. Fractional x-spacing is defined as\n" \
        "       ( x_{i+1} - x_i ) / ( fabs( x_{i+1} ) - fabs( x_i ) ).\n" \
        "\nArguments are:\n" \
        "   epsilon      The fraction minimum distance between sequential x-value.\n" },
    { "toString", (PyCFunction) pointwiseXY_C_toString, METH_VARARGS | METH_KEYWORDS, "Returns a string representation of self." \
        " This method has three keyword parameters:\npairsPerLine, format and pairSeparator which are defined as,\n" \
        "    pairsPerLine    the number of pairs to put on each line\n" \
        "    format          a valid format to convert an (x,y) pair (i.e., two floats) into a string (e.g. format = ' %.3f %12.5e')\n" \
        "    pairSeparator   a string to put between every pair (e.g, to put a comma to separate pairs use pairSeparator = ',')" },
    { "trim", (PyCFunction) pointwiseXY_C_trim, METH_NOARGS, "Returns a new instance with excessive 0. y-value points at beginning and end of self removed." },
    { "union", (PyCFunction) pointwiseXY_C_union, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new pointwiseXY_C object whose x values are the union of self and other." },
    { "scaleOffsetXAndY", (PyCFunction) pointwiseXY_C_scaleOffsetXAndY, METH_VARARGS | METH_KEYWORDS, 
        "Returns a new pointwiseXY_C object whose x and y values are scaled and offset.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xScale   [o] the scale for the x-axis,\n" \
        "   xOffset  [o] the offset for the x-axis,\n" \
        "   yScale   [o] the scale for the y-axis,\n" \
        "   yOffset  [o] the offset for the y-axis,\n" \
        "   insitu   [o] If True, self is scaled and offset, and returned.\n Otherwise, a new instances is created and returned, and self is unchanged." },
    { "mergeClosePoints", (PyCFunction) pointwiseXY_C_mergeClosePoints, METH_VARARGS, 
        "Returns a new pointwiseXY_C object whose x values are thinned. The first and last point of *self* are not altered but intermediate points are average so that the fractional distance between points is greater than *epsilon*." \
        "\nArguments are: ([o] implies optional argument)\n" \
        "epsilon         the minimum fractional distance between x values in the return instance.\n" },
    { "domainGrid", (PyCFunction) pointwiseXY_C_domainGrid, METH_VARARGS, "Returns a list of x-values for self." },
    { "domain", (PyCFunction) pointwiseXY_C_domain, METH_NOARGS, "Returns the x-value of the first and last points as a tuple." },
    { "domainMin", (PyCFunction) pointwiseXY_C_domainMin, METH_NOARGS, "Returns the x-value of the first point." },
    { "domainMax", (PyCFunction) pointwiseXY_C_domainMax, METH_NOARGS, "Returns the x-value of the last point." },
    { "domainSlice", (PyCFunction) pointwiseXY_C_domainSlice, METH_VARARGS | METH_KEYWORDS, "Returns a new instance with self sliced between domainMin and domainMax.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   domainMin    [o] the lower x-value of the slice, default is domainMin of self,\n" \
        "   domainMax    [o] the upper x-value of the slice, default is domainMax of self,\n" \
        "   fill    [o] if True, points are added at domainMin and domainMax if they are not in self,\n" \
        "               else only existing points in the range [domainMin, domainMax] are included.\n" \
        "   dullEps [o] (Currently not implemented) the lower and upper points are dulled, default is 0." },
    { "range",    (PyCFunction)    pointwiseXY_C_range, METH_NOARGS, "Returns the minimum and maximum y-values in self or 0 if self is empty." },
    { "rangeMin", (PyCFunction) pointwiseXY_C_rangeMin, METH_NOARGS, "Returns the minimum y-value in self or 0 if self is empty." },
    { "rangeMax", (PyCFunction) pointwiseXY_C_rangeMax, METH_NOARGS, "Returns the maximum y-value in self or 0 if self is empty." },
    { "ysMappedToXs", (PyCFunction) pointwiseXY_C_ysMappedToXs, METH_VARARGS, "Returns the tuple (offset, y-values) where y-values is self's y-values mapped to Xs.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   Xs  [o] List of x-values to map self to." },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
static PyTypeObject pointwiseXY_CPyType = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    "pointwiseXY_C.pointwiseXY_C",          /* tp_name        */
    sizeof( pointwiseXY_CPy ),              /* tp_basicsize   */
    0,                                      /* tp_itemsize    */
    /* methods */ 
    (destructor) pointwiseXY_C_dealloc,     /* tp_dealloc     */
    0,                                      /* tp_print       */
    0,                                      /* tp_getattr     */
    0,                                      /* tp_setattr     */
    0,                                      /* tp_compare     */
    (reprfunc) pointwiseXY_C__repr__,       /* tp_repr        */
    &pointwiseXY_CPy_number,                /* tp_as_number   */
    &pointwiseXY_CPy_sequence,              /* tp_as_sequence */
    pointwiseXY_CPy_mappingPointer,         /* tp_as_mapping  */
    0,                                      /* tp_hash        */
    0,                                      /* tp_call        */
    0,                                      /* tp_str         */
    0,                                      /* tp_getattro    */
    0,                                      /* tp_setattro    */
    0,                                      /* tp_as_buffer   */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES, /* tp_flags     */
    pointwiseXY_C__doc__,                   /* tp_doc         */
    0,                                      /* tp_traverse    */
    0,                                      /* tp_clear       */
    0,                                      /* tp_richcompare */
    0,                                      /* tp_weaklistoffset */
    0,                                      /* tp_iter        */
    0,                                      /* tp_iternext    */
    pointwiseXY_CPyMethods,                 /* tp_methods     */
    0,                                      /* tp_members     */
    0,                                      /* tp_getset      */
    0,                                      /* tp_base        */
    0,                                      /* tp_dict        */
    0,                                      /* tp_descr_get   */
    0,                                      /* tp_descr_set   */
    0,                                      /* tp_dictoffset  */
    (initproc) pointwiseXY_C__init__,       /* tp_init        */
    0,                                      /* tp_alloc       */
    0                                       /* tp_new         */
};
/*
************************************************************
*/

static PyMethodDef pointwiseXY_CMiscPyMethods[] = {

    { "defaultAccuracy", (PyCFunction) pointwiseXY_C_defaultAccuracy, METH_NOARGS,
        "defaultAccuracy( )\n\nReturns the value of the default accuracy for the __init_ method." },
    { "createFromFunction", (PyCFunction) pointwiseXY_C_createFromFunction, METH_VARARGS | METH_KEYWORDS,
        "createFromFunction( xs, f, parameters, accuracy, biSectionMax, checkForRoots = False, infill = True, saveDivide = True )\n\n" \
        "Returns a pointwiseXY_C instance which represents the function f(x) where f is the second argument and must be a\n" \
        "function reference. The function f must take two arguments. The first arugment to f is x and the second is the third\n" \
        "argument to createFromFunction. The function f must return, as a float, the y-value at x. The returned pointwiseXY_C\n" \
        "instance may contain more points then xs (the first argument) as createFromFunction does infilling. That is, points\n" \
        "are added recursively between the points in xs until accuracy or biSectionMax is reached.\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   xs            An ascending list of x-values to use as an initial guess for the refinement of the returned pointwiseXY_C instance.\n" \
        "   f             A function that returns the y-value given an x-value.\n" \
        "   parameters    Any python object which may be needed by f. Passed as the second argument to f.\n" \
        "   accuracy      The accuracy for infilling.\n" \
        "   biSectionMax  The maximum number of bi-sections for infilling.\n" \
        "   checkForRoots [o] If true, and biSectionMax > 0, an additional point is added whenever two consecutive points cross the y-axis\n" \
        "                     at (or close to) the root (i.e., the y crossing point) (default is False).\n" \
        "   infill        [o] Infill value used for the returned pointwiseXY_C instance (default is True).\n" \
        "   safeDivide    [o] SafeDivide value used for the returned pointwiseXY_C instance (default is True).\n" \
        "\nExample: Return a pointwise representation of 'x * sin( x**2 )' in the domain [1, 10].\n\n" \
        "import math\n" \
        "import pointwiseXY_C\n" \
        "def f( x, args ) :\n" \
        "   return( x * math.sin( x * x ) )\n" \
        "\n" \
        "xSin_xx = pointwiseXY_C.createFromFunction( ( 1, 10 ), f, None, 1e-3, 12 )\n" \
        "\n\n# A better solution may be to input the known zeros of f and force them to be 0., as\n" \
        "\n" \
        "xs = [ 1 ] + [ math.sqrt( i * math.pi ) for i in xrange( 1, int( math.sqrt( 100 / math.pi ) ) ) ] + [ 10 ]\n" \
        "xSin_xx = pointwiseXY_C.createFromFunction( xs, f, None, 1e-3, 12, checkForRoots = True )" },
    { "createFromString", (PyCFunction) pointwiseXY_C_createFromString, METH_VARARGS | METH_KEYWORDS,
        "createFromString( str, accuracy, biSectionMax, interpolation = 'lin-lin', infill = True, safeDivide = True )\n\n" \
        "Returns a tuple of two elements. The first element is a pointwiseXY_C instance representing the float values\n" \
        "translated from 'str'. The second element is the portion of 'str' not translated\n" \
        "\nArguments are: ([o] implies optional argument)\n" \
        "   str           The containing a list of floats to be converted.\n" \
        "   accuracy      The accuracy for infilling.\n" \
        "   biSectionMax  The maximum number of bi-sections for infilling.\n" \
        "   interpolation [o] the interpolation string (default is 'lin-lin').\n" \
        "   infill        [o] Infill value used for the returned pointwiseXY_C instance (default is True).\n" \
        "   safeDivide    [o] SafeDivide value used for the returned pointwiseXY_C instance (default is True).\n" },
    { "interpolatePoint", (PyCFunction) pointwiseXY_C_ptwXY_interpolatePoint, METH_VARARGS,
        "interpolatePoint( interpolation, x, x1, y1, x2, y2 )\n\n" \
        "Returns interpolation of x for x1 <= x <= x2 given x1, y1, x2, y2 and interpolation law." \
        "\n" \
        "Arguments are:\n" \
        "   interpolation   a string representing the interpolation law (e.g., 'log-lin'; see constructor's docstring),\n" \
        "   x               x point of interpolated y-value,\n" \
        "   x1              lower x-value,\n" \
        "   y1              y(x1),\n" \
        "   x2              upper x-value,\n" \
        "   y2              y(x2)." },
    { "gaussian", (PyCFunction) pointwiseXY_C_gaussian, METH_VARARGS | METH_KEYWORDS, 
        "gaussian( accuracy, domainMin, domainMax, offset = 0., sigma = 1., amplitude = 1., dullEps = False )\n\n" \
        "Returns a new pointwiseXY_C instance constructed from the following equation\n\n" \
        "       amplitude * exp( ( ( x - offset ) / sigma )^2 / 2 )        for domainMin <= x <= domainMax\n" \
        "\n" \
        "Arguments are:  ([o] implies optional argument)\n" \
        "   accuracy        the accuracy of linear interpolation,\n" \
        "   domainMin            the lower x point generated,\n" \
        "   domainMax            the upper x point generated,\n" \
        "   offset     [o]  the x offset of the center of the Gaussian,\n" \
        "   sigma      [o]  width of the Gaussian,\n" \
        "   amplitude  [o]  the Gaussian's amplitude,\n" \
        "   dullEps    [o]  currently not used." },
    { "basicGaussian", (PyCFunction) pointwiseXY_C_basicGaussian, METH_VARARGS, 
        "basicGaussian( accuracy )\n\n" \
        "Returns a new pointwiseXY_C instance constructed from the following equation\n\n" \
        "       exp( ( x^2 / 2 )\n" \
        "\n" \
        "Arguments are:\n" \
        "   accuracy        the accuracy of linear interpolation.\n" },
    { "unitbaseInterpolate", (PyCFunction) pointwiseXY_C_unitbaseInterpolate, METH_VARARGS,
        "unitbaseInterpolate( x, lw, lXY, uw, uXY )\n\n"
        "Returns the unitbase interpolation of two XYs objects at w where the axes are labeled (w, x, y).\n" \
        "\n" \
        "Arguments are:\n" \
        "   w           the point between lw and uw to return the unitbase interpolation of lXY and uXY,\n" \
        "   lw          the w point where lXY resides,\n" \
        "   lXY         a pointwiseXY_C instance for a function y(x),\n" \
        "   uw          the w point where uXY resides,\n" \
        "   uXY         a pointwiseXY_C instance for a function y(x),\n" \
        " scaleRange    if True range values are scaled, otherwise they are unchanged.\n" },
    { "floatToShortestString", (PyCFunction) floatToShortestString_C, METH_VARARGS | METH_KEYWORDS,
        "floatToShortestString( value, significantDigits = 15, trimZeros = True, keepPeriod = False,\n" \
        "        favorEFormBy = 0, includeSign = False )\n\n" \
        "Returns the float converted to either the E-form (i.e., '%e') and F-form (i.e., '%f') with significantDigits.\n" \
        "The form with the shortest string representation of the float value is returned.\n" \
        "\n" \
        "Arguments are:\n" \
        "   value               the float value to convert to a string,\n" \
        "   significantDigits   The number of significant digits desired. Restricted to the range 0 to 24 inclusive,\n" \
        "   trimZeros           If True, unneeded zeros to the right of '.' are removed,\n" \
        "   keepPeriod          If False, '.' is removed if there is no digit to its right,\n" \
        "   favorEFormBy        The integer subtracted from the length of the E-form before the form\n" \
        "                       with the shortest representation is determined,\n" \
        "   includeSign         If True, the returned string will always start with a sign character\n" \
        "                       (i.e., '+' or '-'). Otherwise, only negative values will have a sign.\n" },
    { NULL, NULL, 0, NULL }        /* Sentinel (i.e., the end of the list) */
};
/*
************************************************************
*/
ptwXYPoints *pointwiseXY_C_factory_get_ptwXYPoints( PyObject *ptwXY_Py ) {

    int status = PyObject_IsInstance( ptwXY_Py, (PyObject* ) &pointwiseXY_CPyType );

    if( status < 0 ) return( NULL );
    if( status == 0 ) return( (ptwXYPoints *) pointwiseXY_C_SetPyErrorExceptionReturnNull( "Object is not a pointwiseXY_C instance" ) );
    if( pointwiseXY_C_checkStatus( (pointwiseXY_CPy *) ptwXY_Py ) != 0 ) return( NULL );
    return( ((pointwiseXY_CPy *) ptwXY_Py)->ptwXY );
}
/*
************************************************************
*/
static PyObject *pointwiseXY_C_factory_create( ptwXYPoints *ptwXY, int infill, int safeDivide ) {

    pointwiseXY_CPy *nPy = pointwiseXY_CNewInitialize( infill, safeDivide );

    if( nPy == NULL ) {
        ptwXY_free( ptwXY ); }
    else {
        nPy->ptwXY = ptwXY;
    }
    return( (PyObject *) nPy );
}
/*
************************************************************
*/

static char const doc[] = "A module that contains the class pointwiseXY_C.";

#if PY_MAJOR_VERSION < 3

    PyMODINIT_FUNC initpointwiseXY_C( void ) {

        PyObject *module;
        static void *Py_pointwiseXY_C_API[Py_pointwiseXY_C_API_numberOfPointers];
        PyObject *c_api_object;

        pointwiseXY_CPyType.tp_new = PyType_GenericNew;
        if( PyType_Ready( &pointwiseXY_CPyType ) < 0 ) return;

        if( ( module = Py_InitModule3( "pointwiseXY_C", pointwiseXY_CMiscPyMethods, doc ) ) == NULL ) return;

        Py_INCREF( &pointwiseXY_CPyType );
        PyModule_AddObject( module, "pointwiseXY_C", (PyObject *) &pointwiseXY_CPyType );

        Py_pointwiseXY_C_API[Py_pointwiseXY_C_get_ptwXYPoints_NUM] = (void *) pointwiseXY_C_factory_get_ptwXYPoints;
        Py_pointwiseXY_C_API[Py_pointwiseXY_C_NUM] = (void *) pointwiseXY_C_factory_create;

/* Cannot get the PyCapsule_Import to work for python 2.7.  */
        c_api_object = PyCObject_FromVoidPtr( (void *) Py_pointwiseXY_C_API, NULL );

        if( c_api_object != NULL ) PyModule_AddObject( module, "_C_API", c_api_object );
    }

#else

    static struct PyModuleDef pointwiseXY_CModule = {
        PyModuleDef_HEAD_INIT,
        "pointwiseXY_C",            /* name of module */
        doc,                        /* module documentation, may be NULL */
        -1,                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        pointwiseXY_CMiscPyMethods
    };

    PyMODINIT_FUNC PyInit_pointwiseXY_C( void ) {

        PyObject *module;
        static void *Py_pointwiseXY_C_API[Py_pointwiseXY_C_API_numberOfPointers];
        PyObject *c_api_object;

        if( ( module = PyModule_Create( &pointwiseXY_CModule ) ) == NULL ) return( NULL );

        pointwiseXY_CPyType.tp_new = PyType_GenericNew;
        if( PyType_Ready( &pointwiseXY_CPyType ) < 0 ) return( NULL );
        Py_INCREF( &pointwiseXY_CPyType );
        PyModule_AddObject( module, "pointwiseXY_C", (PyObject *) &pointwiseXY_CPyType );

        Py_pointwiseXY_C_API[Py_pointwiseXY_C_get_ptwXYPoints_NUM] = (void *) pointwiseXY_C_factory_get_ptwXYPoints;
        Py_pointwiseXY_C_API[Py_pointwiseXY_C_NUM] = (void *) pointwiseXY_C_factory_create;

        c_api_object = PyCapsule_New( (void *) Py_pointwiseXY_C_API, "numericalFunctions.pointwiseXY_C._C_API", NULL );
        if( c_api_object != NULL ) PyModule_AddObject( module, "_C_API", c_api_object );

        return( module );
    }

#endif
