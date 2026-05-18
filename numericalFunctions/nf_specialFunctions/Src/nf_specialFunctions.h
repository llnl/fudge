/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#ifndef specialFunctions_h_included
#define specialFunctions_h_included

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

double nf_polevl( double x, double coef[], int N );
double nf_p1evl( double x, double coef[], int N );
nfu_status nf_exponentialIntegral( statusMessageReporting *smr, int n, double x, double *value );
nfu_status nf_gammaFunction( statusMessageReporting *smr, double x, double *value );
nfu_status nf_logGammaFunction( statusMessageReporting *smr, double x, double *value );
nfu_status nf_incompleteGammaFunction( statusMessageReporting *smr, double a, double x, double *value );
nfu_status nf_incompleteGammaFunctionComplementary( statusMessageReporting *smr, double a, double x, double *value );

double  nf_amc_log_factorial( int );
double  nf_amc_factorial( int );
double  nf_amc_wigner_3j( int, int, int, int, int, int );
double  nf_amc_wigner_6j( int, int, int, int, int, int );
double  nf_amc_wigner_9j( int, int, int, int, int, int, int, int, int );
double  nf_amc_racah( int, int, int, int, int, int );
double  nf_amc_clebsh_gordan( int, int, int, int, int );
double  nf_amc_z_coefficient( int, int, int, int, int, int );
double  nf_amc_zbar_coefficient( int, int, int, int, int, int );
double  nf_amc_reduced_matrix_element( int, int, int, int, int, int, int );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwXY_h_included. */
