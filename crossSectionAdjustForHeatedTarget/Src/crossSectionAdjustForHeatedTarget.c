/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "crossSectionAdjustForHeatedTarget.h"

#define sqrtpi 1.7724538509055160273
static const double vCutoffRatio = 6., fInterpolationMin = 1e-6, fInterpolationMax = 0.1, xsecMin1 = 1e-20, xsecMin2 = 1e-40;

#define FreeReturn( error ) { Free_E_cs_point( &E_cs_Info ); return( error ); }

static double xMaxForExp = .98 * DBL_MAX_10_EXP * 2.302585;         /* 2.302585 = log( 10 ), the 0.98 is a fudge factor */
static double xMaxForErfc = -1;

static double crossSectionAdjustForHeatedTarget_next_point( double mass_Ratio, double T, double E, double fStep, int *index, 
    double prior_dE, E_cs_heated_point_Info *E_cs_Info );
static int crossSectionAdjustForHeatedTarget_thicken( E_cs_heated_point *E_cs_point, E_cs_heated_point_Info *E_cs_Info, int checkPoint );
static int crossSectionAdjustForHeatedTarget_add_point( E_cs_heated_point *E_cs_point, E_cs_heated_point *E_cs_point_prior,
    E_cs_heated_point_Info *E_cs_Info );
#if 0
static double crossSectionAdjustForHeatedTarget_sum_factors( double s, double E, double Ei, double csi, double SEi, double K[5], double dxnerf[5] );
#endif
static double oneOverVStepSize( double f );
static double myExp( double x );
static double myErf( double x );
static double myErfc( double x );
static double mySinh( double x );
static double crossSectionAdjustForHeatedTarget_interpolate_cs( double E, int index, double *E_cs_in );
/*
****************************************************
*/
static void Free_E_cs_point( E_cs_heated_point_Info *p ) {

    free( p->E_cs_in );
    if( p->E_cs != NULL ) free( p->E_cs );
}
/*
****************************************************
*/
int crossSectionAdjustForHeatedTarget_init( crossSectionAdjustForHeatedTarget_limit lowerlimit, crossSectionAdjustForHeatedTarget_limit upperlimit,
    crossSectionAdjustForHeatedTarget_info *info, double mass_Ratio, double T, double fInterpolation, int n_pairs, double *E_cs_in, 
    E_cs_heated_point_Info *E_cs_Info ) {
/*
*    This routine initializes E_cs_Info, info and xMaxForErfc.
*/
    int i, j, nAddedPoints = 0;
    double dE, *p;

    E_cs_Info->mass_Ratio = mass_Ratio;
    E_cs_Info->T = T;
    E_cs_Info->csMax = 0.;
    E_cs_Info->lowerlimit = lowerlimit;
    E_cs_Info->upperlimit = upperlimit;
    E_cs_Info->info = info;
    E_cs_Info->n_pairs = n_pairs;
    E_cs_Info->E_cs_in_orig = E_cs_in;
    E_cs_Info->E_cs_in = NULL;
    E_cs_Info->n_done = 0;
    E_cs_Info->n_allocated = 0;
    E_cs_Info->dEMin = 0.;
    E_cs_Info->E_cs = NULL;

    if( ( info->mode & crossSectionAdjustForHeatedTarget_mode_all ) || 
            ( info->mode & crossSectionAdjustForHeatedTarget_mode_allEDomain ) ) {
        E_cs_Info->iStart = 0;
        E_cs_Info->iEnd = n_pairs; }
    else {
        for( i = 0; i < n_pairs; i++ ) if( E_cs_in[2*i+1] != 0. ) break;
        if( i > 0 ) i--;
        E_cs_Info->iStart = 1;
        for( i = n_pairs; i > 0; i-- ) if( E_cs_in[2*i-1] != 0. ) break;
        if( i < n_pairs ) i++;
        E_cs_Info->iEnd = i;
    }
    E_cs_Info->iStart_orig = E_cs_Info->iStart;
    E_cs_Info->iEnd_orig = E_cs_Info->iEnd;

    fInterpolation *= 0.5;
    info->bytes = 0;
    info->evaluations = 0;
    info->InfoStats = 0;
    info->WarningStats = 0;
    info->ErrorStats = 0;
    info->resolutionStats = 0;
    if( fInterpolation < fInterpolationMin ) fInterpolation = fInterpolationMin;
    if( fInterpolation > fInterpolationMax ) fInterpolation = fInterpolationMax;
    (E_cs_Info->info)->fInterpolation = fInterpolation;
    E_cs_Info->fInterpolationOffset = log( fInterpolation );
    E_cs_Info->fInterpolationSlope = log( fInterpolationMax / fInterpolation ) / log( xsecMin1 / xsecMin2 );

    if( ( E_cs_Info->iEnd - E_cs_Info->iStart ) < 2 ) return( -1 );
    if( mass_Ratio  <= 0. ) return( -2 );
    if( E_cs_in[0] <= 0. ) return( -3 );
    if( T <= 0. ) return( -4 );
    for( i = 2; i < 2 * n_pairs; i += 2 ) {
        if( E_cs_in[i-2] > E_cs_in[i] ) return( -5 );
        dE = 4. * sqrt( ( E_cs_in[i-2] + E_cs_in[i] ) * T / mass_Ratio );
        if( ( E_cs_in[i] - E_cs_in[i-2] ) > 5 * dE ) nAddedPoints++;
    }
    info->bytes = ( 4 * n_pairs + 8 * nAddedPoints ) * sizeof( double );
    E_cs_Info->E_cs_in = (double *) malloc( info->bytes );
    if( E_cs_Info->E_cs_in == NULL ) return( -6 );
    p = E_cs_Info->E_cs_in;
    *(p++) = E_cs_in[2 * E_cs_Info->iStart_orig];
    *(p++) = E_cs_in[2 * E_cs_Info->iStart_orig + 1];
    for( i = 2 * E_cs_Info->iStart_orig + 2, E_cs_Info->iEnd = 1; i < 2 * E_cs_Info->iEnd_orig; i += 2, (E_cs_Info->iEnd)++ ) {
        dE = 4. * sqrt( ( E_cs_in[i-2] + E_cs_in[i] ) * T / mass_Ratio );
        if( ( E_cs_in[i] - E_cs_in[i-2] ) > 5 * dE ) {
            *(p++) = E_cs_in[i-2] + dE;
            *(p++) = crossSectionAdjustForHeatedTarget_interpolate_cs( E_cs_in[i-2] + dE, i - 2, E_cs_in );
            *(p++) = E_cs_in[i] - dE;
            *(p++) = crossSectionAdjustForHeatedTarget_interpolate_cs( E_cs_in[i] - dE, i - 2, E_cs_in );
            E_cs_Info->iEnd += 2;
        }
        *(p++) = E_cs_in[i];
        *(p++) = E_cs_in[i + 1];
    }
    E_cs_Info->iStart = 0;
    E_cs_Info->sqrtEi = &(E_cs_Info->E_cs_in[2 * E_cs_Info->iEnd] );
    E_cs_Info->SEi = &(E_cs_Info->sqrtEi[E_cs_Info->iEnd]);
    E_cs_Info->E_cs = NULL;

    E_cs_Info->sqrtEi[0] = sqrt( E_cs_Info->E_cs_in[0] );                    /* Precalculate v and cross-section slope. */
    for( i = 1, j = 2; i < E_cs_Info->iEnd; i++, j += 2 ) {
        E_cs_Info->sqrtEi[i] = sqrt( E_cs_Info->E_cs_in[j] );
        E_cs_Info->SEi[i-1] = 0.;
        if( E_cs_Info->E_cs_in[j] != E_cs_Info->E_cs_in[j - 2] ) E_cs_Info->SEi[i-1] = ( E_cs_Info->E_cs_in[j + 1] - E_cs_Info->E_cs_in[j - 1] ) / 
            ( E_cs_Info->E_cs_in[j] - E_cs_Info->E_cs_in[j - 2] );
    }
    return( 0 );
}
/*
****************************************************
*/
int crossSectionAdjustForHeatedTarget( crossSectionAdjustForHeatedTarget_limit lowerlimit, crossSectionAdjustForHeatedTarget_limit upperlimit, 
    crossSectionAdjustForHeatedTarget_info *info, double EMin, double mass_Ratio, double T, double fInterpolation, int n_pairs, 
    double *E_cs_in, double **E_cs_out ) {
/*
*   info              See "Info data" below. If info == NULL then default values are used.
*   EMin              The lowest energy to heat to, example, in ENDL 1e-10 (MeVs) is the lowest energy (must be in same units as E in E_cs_in).
*   mass_Ratio        Ratio of target mass to incident particle mass.
*   T                 Temperature of target (must be in same units as E in E_cs_in).
*   fInterpolation   Data is thicken so that Lin-lin interpolation is accurate to this faction
*                         (e.g., fInterpolation = 0.01 is 1% lin-lin interpolation accuracy).
*   n_pairs           Number of E, cross-section pairs in E_cs_in array.
*   E_cs_in           Array containing E, cross-section pair data
*
* Info data:
*   mode              Determines the spacing of the initially heated points.
*                         if bits 0 is off then a judicial choice is made based on the Temperature, otherwise all points of E_cs_in are heated.
*                         if bits 1 is off then data is thinned, else it is not.
*   verbose           flag used to determine message to print (input):
*                         0 = no printing, 1 = print errors, 2 = 1 + print warning, or 3 = 2 + print info (e.g., print all).
*   InfoStats         Currently not used, set to 0 (output).
*   WarningStats      Number of warnings detected (output).
*   ErrorStats        Currently not used, set to 0 (output).
*   resolutionStats   The number of time interpolation fails for given resolution (output).
*/
    int i, err, iPrior = -1, index, iNotHeated, n, debugFlag = 0;
    double E, EPrior, *d, fStep, prior_dE;
    E_cs_heated_point E_cs_point, E_cs_point_init = { -1, -1, -1, 0., 0. }, *E_cs_point_prior, *E_cs_point_current = NULL, *E_cs_point_next;
    E_cs_heated_point_Info E_cs_Info;
    double x, y, x1, y1, x2, y2, cs, EStart;

    if( ( err = crossSectionAdjustForHeatedTarget_init( lowerlimit, upperlimit, info, mass_Ratio, T, fInterpolation, n_pairs, E_cs_in, &E_cs_Info ) )
        < 0 ) FreeReturn( err );
    fInterpolation = (E_cs_Info.info)->fInterpolation;

    E_cs_point_prior = &E_cs_point_init;

    cs = 0.;                                 /* The following is done to set csMax. */
    for( i = E_cs_Info.iStart; i < E_cs_Info.iEnd; i++ ) if( cs < E_cs_Info.E_cs_in[2 * i + 1] ) cs = E_cs_Info.E_cs_in[2 * i + 1];
    E_cs_Info.csMax = cs;

    EPrior = -1.;
    if( lowerlimit == crossSectionAdjustForHeatedTarget_limit_threshold ) {    /* May need to add energy point below inputted EMin. */
        EStart = E_cs_Info.E_cs_in[2 * E_cs_Info.iStart];
        double lastNonZeroEnergy = EStart;
        if( EMin < EStart ) {
            double ELower = EMin;
            double EUpper = EStart;
            EPrior = EMin;
            i = 0;
            do {
                E = sqrt( ELower * EUpper );
                crossSectionAdjustForHeatedTarget_heat_at_E( E, &E_cs_Info, &E_cs_point );
                if( E_cs_point.cs < xsecMin2 * E_cs_Info.csMax ) {
                    ELower = E; }
                else {
                    if( i == 0 ) {          /* If first E is acceptable, just use EMin. */
                        E = EMin;
                        crossSectionAdjustForHeatedTarget_heat_at_E( E, &E_cs_Info, &E_cs_point );
                        break;
                    }
                    EUpper = E;
                }
                if( ( E < lastNonZeroEnergy ) && ( E_cs_point.cs != 0 ) ) {
                    lastNonZeroEnergy = E;
                }
                i += 1;
            } while( i < 40 );
            if( ( i == 40 ) && ( E_cs_point.cs == 0.0 ) ) { 
                E = lastNonZeroEnergy;
                crossSectionAdjustForHeatedTarget_heat_at_E( E, &E_cs_Info, &E_cs_point );
            }
            if( ( err = crossSectionAdjustForHeatedTarget_add_point( &E_cs_point, E_cs_point_prior, &E_cs_Info ) ) < 0 ) FreeReturn( err );
            E_cs_point_prior = &(E_cs_Info.E_cs[err]);
            EPrior = E;
        }
    }

    fStep = oneOverVStepSize( fInterpolation );
    for( i = E_cs_Info.iStart; i < E_cs_Info.iEnd; ) {                                        /* Loop to heat cross-sections at each input energy. */
        index = i;
        prior_dE = 0.;
        if( iPrior > -1 ) {
            E_cs_point_prior = &(E_cs_Info.E_cs[iPrior]);
            if( E_cs_point_prior->prior != -1 ) prior_dE = E_cs_point_prior->E - E_cs_Info.E_cs[E_cs_point_prior->prior].E;
        }
        if( info->mode & crossSectionAdjustForHeatedTarget_mode_all ) {        /* Heat all inputted points. */
            E = E_cs_Info.E_cs_in[2 * i]; }
        else {                                                            /* Pick points judicially. */
            E = crossSectionAdjustForHeatedTarget_next_point( mass_Ratio, T, EPrior, fStep, &index, prior_dE, &E_cs_Info );
        }
        crossSectionAdjustForHeatedTarget_heat_at_E( E, &E_cs_Info, &E_cs_point );
        if( ( err = crossSectionAdjustForHeatedTarget_add_point( &E_cs_point, E_cs_point_prior, &E_cs_Info ) ) < 0 ) FreeReturn( err );
        iPrior = err;
        if( ( E_cs_Info.n_done > 1 ) && ( E != EPrior ) ) {                /* Check if we need to thicken Energy points. */
            if( ( err = crossSectionAdjustForHeatedTarget_thicken( &E_cs_point, &E_cs_Info, 0 ) ) < 0 ) FreeReturn( err );
        }
        EPrior = E;
        if( info->mode & crossSectionAdjustForHeatedTarget_mode_all ) {
            i++; }
        else {
            i = index;
        }
        if( ( !(info->mode & crossSectionAdjustForHeatedTarget_mode_allEDomain) ) && ( E > 1e7 * T / mass_Ratio ) ) {
            if( ( index > E_cs_Info.iStart ) && ( index < E_cs_Info.iEnd ) ) {
                cs = E_cs_Info.SEi[index - 1] * ( E - E_cs_Info.E_cs_in[2 * index - 2] ) + E_cs_Info.E_cs_in[2 * index - 1];
                if( fabs( cs - E_cs_point.cs ) < E_cs_point.cs * fInterpolation ) break;
            }
        }
    }

    iNotHeated = i;

    if( debugFlag && ( E_cs_Info.n_done > 0 ) ) {       /* The code in this if statement is only turned on in a debugger. */
        int i = 0;
        E_cs_heated_point *E_cs;
        info->bytes = E_cs_Info.n_done * sizeof( E_cs_heated_point );
        E_cs = (E_cs_heated_point *) malloc( info->bytes );
        E_cs_point_current = &(E_cs_Info.E_cs[0]);

        do {
            E_cs[i] = *E_cs_point_current;
            E_cs[i].index = i;
            E_cs[i].prior = i - 1;
            E_cs[i].next = i + 1;
            E_cs_point_prior = E_cs_point_current;
            if( E_cs_point_prior->next != -1 ) {
                i++;
                E_cs_point_current = &(E_cs_Info.E_cs[E_cs_point_prior->next]);
            }
        } while( E_cs_point_prior->next != -1 );
        E_cs[i].next = -1;
        E_cs_Info.n_allocated = E_cs_Info.n_done;
        free( E_cs_Info.E_cs );
        E_cs_Info.E_cs = E_cs;
    }

    if( !( info->mode & crossSectionAdjustForHeatedTarget_mode_do_not_thin ) && ( E_cs_Info.n_done > 0 ) ) {    /* Now, thin the data. */

        E_cs_point_prior = &(E_cs_Info.E_cs[0]);            /* This loop removes the middle points in a region where all the cs values are the same. */
        for( i = 0; i < E_cs_Info.n_done - 1; i++ ) {
            E_cs_point_current = &(E_cs_Info.E_cs[E_cs_point_prior->next]);
            E_cs_point_next = &(E_cs_Info.E_cs[E_cs_point_current->next]);
            while( ( E_cs_point_prior->cs == E_cs_point_current->cs ) && ( E_cs_point_prior->cs == E_cs_point_next->cs ) 
                    && ( E_cs_point_next->next != -1 ) ) {
                E_cs_point_prior->next = E_cs_point_current->next;
                E_cs_point_current = E_cs_point_next;
                E_cs_point_next = &(E_cs_Info.E_cs[E_cs_point_next->next]);
                E_cs_Info.n_done--;
            }
            E_cs_point_prior = E_cs_point_current;
        }

        E_cs_point_prior = &(E_cs_Info.E_cs[0]);
        E_cs_point_prior->thinnable = 0;
        if( E_cs_Info.n_done > 1 ) {
            E_cs_point_current = &(E_cs_Info.E_cs[E_cs_point_prior->next]);
            E_cs_point_current->thinnable = 0;                                /* Second point is not thinnable. */
        }
        if( ( E_cs_Info.n_done > 2 ) ) {
            do {                                                            /* Do not thin local minimums and maximums. */
                E_cs_point_next =  &(E_cs_Info.E_cs[E_cs_point_current->next]);
                if( (E_cs_point_next->cs - E_cs_point_current->cs) * (E_cs_point_current->cs - E_cs_point_prior->cs) < 0. ) E_cs_point_current->thinnable = 0;
                E_cs_point_prior = E_cs_point_current;
                E_cs_point_current = E_cs_point_next;
            } while( E_cs_point_current->next != -1 );
            E_cs_point_prior->thinnable = 0;                                /* Second to last point is not thinnable. */
            E_cs_point_current->thinnable = 0;                                /* Last point is not thinnable. */

            E_cs_point_prior = &(E_cs_Info.E_cs[0]);
            do {                                                            /* Now thin the other points. */
                x1 = E_cs_point_prior->E;
                y1 = E_cs_point_prior->cs;
                E_cs_point_next = &(E_cs_Info.E_cs[E_cs_point_prior->next]);
                if( E_cs_point_next->next == -1 ) break;
                do {                                                        /* Increment E_cs_point_next until end of list or unthinnable point found. */
                    E_cs_point_current =  &(E_cs_Info.E_cs[E_cs_point_prior->next]);
                    E_cs_point_next =  &(E_cs_Info.E_cs[E_cs_point_next->next]);
                    x2 = E_cs_point_next->E;
                    y2 = E_cs_point_next->cs;
                    do {                                                    /* Increment E_cs_point_current until E_cs_point_next or unthinnable point found. */
                        x = E_cs_point_current->E;
                        y = E_cs_point_current->cs;
                        if( fabs( ( y2 - y1 ) * ( x - x1 ) + ( y1 - y ) * ( x2 - x1 ) ) > fInterpolation * y * ( x2 - x1 ) ) E_cs_point_current->thinnable = 0;
                        if( !E_cs_point_current->thinnable ) break;
                        E_cs_point_current =  &(E_cs_Info.E_cs[E_cs_point_current->next]);
                    } while( E_cs_point_current != E_cs_point_next );
                } while( ( E_cs_point_next->next != -1 ) && ( E_cs_point_current == E_cs_point_next ) );
                E_cs_point_prior = E_cs_point_current;
            } while( E_cs_point_prior->next != -1 );
        }
    }

    n = 1;                                                                    /* Start at 1 since last point is not checked. */
    if( E_cs_Info.n_done == 0 ) n = 0;
    for( E_cs_point_current = &(E_cs_Info.E_cs[0]); E_cs_point_current->next != -1; E_cs_point_current = &(E_cs_Info.E_cs[E_cs_point_current->next]) ) {
        if( !E_cs_point_current->thinnable ) n++;
    }
    info->bytes = 2 * ( n + ( E_cs_Info.iEnd - iNotHeated ) ) * sizeof( double );
    *E_cs_out = (double *) malloc( info->bytes );
    if( E_cs_out == NULL ) FreeReturn( -7 );
    d = *E_cs_out;
    E_cs_point_current = E_cs_Info.E_cs;
    for( i = 0; i < E_cs_Info.n_done; i++ ) {
        if( !E_cs_point_current->thinnable ) {
            *(d++) = E_cs_point_current->E;
            *(d++) = E_cs_point_current->cs;
        }
        E_cs_point_current = &(E_cs_Info.E_cs[E_cs_point_current->next]);
    }
    for( i = iNotHeated; i < E_cs_Info.iEnd; i++ ) {
        *(d++) = E_cs_Info.E_cs_in[2 * i];
        *(d++) = E_cs_Info.E_cs_in[2 * i + 1];
    }

    FreeReturn( n + ( E_cs_Info.iEnd - iNotHeated ) );
}
/*
****************************************************
*/
static double crossSectionAdjustForHeatedTarget_next_point( double mass_Ratio, double T, double E, double fStep, int *index,
    double prior_dE, E_cs_heated_point_Info *E_cs_Info ) {

    int dEMode = 0;
    double dE, EPrior = E, dEMin = sqrt( E * T / mass_Ratio ), *E_cs_in = E_cs_Info->E_cs_in;

    if( E < 0. ) {                                            /* The first point. */
        *index = E_cs_Info->iStart + 1;
        E = E_cs_in[2 * E_cs_Info->iStart]; }
    else if( E < T / mass_Ratio ) {
        E *= fStep; }
    else {
        dE = 0.2 * dEMin;
        prior_dE *= 2.0;
        if( dE < prior_dE ) {
            dEMode = 1;
            dE = prior_dE;
        }
        E += dE;
        if( *index < E_cs_Info->iEnd ) {
            if( ( *index + 1 == E_cs_Info->iEnd ) && ( E > E_cs_in[2 * *index] ) ) E = E_cs_in[2 * *index];
            if( E > E_cs_in[2 * ( E_cs_Info->iEnd - 1 )] ) E = E_cs_in[2 * ( E_cs_Info->iEnd - 1 )];
            while( ( *index < E_cs_Info->iEnd ) && ( E_cs_in[2 * *index] <= ( E + 0.01 * dE ) ) ) {
                if( ( dEMode == 1 ) && ( EPrior < E_cs_in[2 * *index] ) ) {
                    E = E_cs_in[2 * *index];
                    (*index)++;
                    break;
                }
                (*index)++;
            } }
        else {
            E = E_cs_in[2 * ( E_cs_Info->iEnd - 1 )];
        }
    }
    if( E >= E_cs_in[2 * ( E_cs_Info->iEnd - 1 )] ) {
        *index = E_cs_Info->iEnd;
        E = E_cs_in[2 * ( E_cs_Info->iEnd - 1 )];
    }
    E_cs_Info->dEMin = (E_cs_Info->info)->fInterpolation * dEMin;
    if( E_cs_Info->dEMin > (E_cs_Info->info)->fInterpolation * E ) E_cs_Info->dEMin = (E_cs_Info->info)->fInterpolation * E;
    E_cs_Info->dEMin *= 0.1;
    return( E );
}
/*
****************************************************
*/
static int crossSectionAdjustForHeatedTarget_thicken( E_cs_heated_point *E_cs_point, E_cs_heated_point_Info *E_cs_Info, int checkPoint ) {

    int err;
    E_cs_heated_point E_cs_point_mid, *E_cs_point_prior = &(E_cs_Info->E_cs[E_cs_point->prior]);
    double fInterpolation = (E_cs_Info->info)->fInterpolation, epsE = 1e-8;
    double E_mid  = 0.5 * ( E_cs_point->E  + E_cs_point_prior->E  );
    double cs_mid = 0.5 * ( E_cs_point->cs + E_cs_point_prior->cs );

    if( checkPoint && ( ( E_mid -  E_cs_point_prior->E ) < E_cs_Info->dEMin ) ) {
        if( ( (E_cs_Info->info)->verbose > 1 ) && ( (E_cs_Info->info)->resolutionStats < 256 ) ) 
            fprintf( stderr, "crossSectionAdjustForHeatedTarget: resolution accuracy could not be met at E = %.12e\n", E_mid );
        (E_cs_Info->info)->resolutionStats++;
        return( 0 );
    }
    if( checkPoint && ( E_mid -  E_cs_point_prior->E ) < epsE * E_mid ) {
        if( ( (E_cs_Info->info)->verbose > 1 ) && ( (E_cs_Info->info)->WarningStats < 256 ) )
            fprintf( stderr, "crossSectionAdjustForHeatedTarget: possible infinite loop detected at E = %.12e\n", E_mid );
        (E_cs_Info->info)->WarningStats++;
        return( 0 );
    }
    crossSectionAdjustForHeatedTarget_heat_at_E( E_mid, E_cs_Info, &E_cs_point_mid );
    if( E_cs_point_mid.cs < xsecMin1 * E_cs_Info->csMax ) {
        if( E_cs_point_mid.cs < xsecMin2 * E_cs_Info->csMax ) return( 0 );
        fInterpolation = exp( E_cs_Info->fInterpolationSlope * log( xsecMin2 * E_cs_Info->csMax / ( xsecMin1 * E_cs_point_mid.cs ) ) 
            + E_cs_Info->fInterpolationOffset );
    }
    if( fabs( cs_mid - E_cs_point_mid.cs ) <= cs_mid * fInterpolation ) return( 0 );
    if( ( err = crossSectionAdjustForHeatedTarget_add_point( &E_cs_point_mid, E_cs_point_prior, E_cs_Info ) ) < 0 ) return( err );
    if( ( err = crossSectionAdjustForHeatedTarget_thicken( &E_cs_point_mid, E_cs_Info, 1 ) ) < 0 ) return( err );
    E_cs_point->prior = E_cs_Info->E_cs[E_cs_point->index].prior;
    err = crossSectionAdjustForHeatedTarget_thicken( E_cs_point, E_cs_Info, 1 );
    return( err );
}
/*
****************************************************
*/
static int crossSectionAdjustForHeatedTarget_add_point( E_cs_heated_point *E_cs_point, E_cs_heated_point *E_cs_point_prior,
    E_cs_heated_point_Info *E_cs_Info ) {

    int i, n_want;
    E_cs_heated_point *p, E_cs_point_prior_Saved;

    if( ( E_cs_Info->n_done + 1 ) > E_cs_Info->n_allocated ) {
        E_cs_point_prior_Saved = *E_cs_point_prior;
        E_cs_point_prior = &E_cs_point_prior_Saved;
        n_want = E_cs_Info->n_allocated + 1000;
        E_cs_Info->E_cs = (E_cs_heated_point *) realloc( E_cs_Info->E_cs, n_want * sizeof( E_cs_heated_point ) );
        if( E_cs_Info->E_cs == NULL ) return( -11 );
        p = &(E_cs_Info->E_cs[E_cs_Info->n_allocated]);
        for( i = E_cs_Info->n_allocated; i < n_want; i++, p++ ) {
            p->index = i;
            p->next = -1;
            p->prior = -1;
            p->thinnable = !((E_cs_Info->info)->mode & crossSectionAdjustForHeatedTarget_mode_do_not_thin);
            p->E = 0.;
            p->cs = 0.;
        }
        E_cs_Info->n_allocated = n_want;
    }
    E_cs_point->prior = E_cs_point_prior->index;
    E_cs_point->next = E_cs_point_prior->next;
    E_cs_point->thinnable = !((E_cs_Info->info)->mode & crossSectionAdjustForHeatedTarget_mode_do_not_thin);
    E_cs_point->index = E_cs_Info->E_cs[E_cs_Info->n_done].index;
    E_cs_Info->E_cs[E_cs_Info->n_done] = *E_cs_point;
    if( E_cs_point->prior != -1 ) E_cs_Info->E_cs[E_cs_point->prior].next = E_cs_point->index;
    if( E_cs_point->next  != -1 ) E_cs_Info->E_cs[E_cs_point->next].prior = E_cs_point->index;
    E_cs_Info->n_done++;
    return( E_cs_Info->n_done - 1 );
}
/*
****************************************************
*/
void crossSectionAdjustForHeatedTarget_heat_at_E( double E, E_cs_heated_point_Info *E_cs_Info, E_cs_heated_point *E_cs_point ) {

    int i, doAgain = 1;
    double a, b, SqrtE = sqrt( E ), InvSqrtE, K[5], *E_cs_in = &(E_cs_Info->E_cs_in[2 * E_cs_Info->iStart]), Ei, csi, pm, scale;
    double *sqrtEi = &(E_cs_Info->sqrtEi[E_cs_Info->iStart]), *SEi = &(E_cs_Info->SEi[E_cs_Info->iStart]), Cutoff, cs, cs_m = 0, cs_p = 0, cs_l = 0, dcs;
    double dxnerf[5];

    (E_cs_Info->info)->evaluations++;
    InvSqrtE = 1. / SqrtE;
    K[1] = E_cs_Info->T / ( E * E_cs_Info->mass_Ratio );        /* K      */
    K[0] = sqrt( K[1] );                                        /* K^1/2  */
    K[2] = K[0] * K[1];                                         /* K^3/2  */
    K[3] = K[1] * K[1];                                         /* K^2    */
    K[4] = 1. / K[0];                                           /* K^-1/2 */
    Cutoff = ( vCutoffRatio * K[0] - 1. ) * SqrtE;

    for( i = E_cs_Info->iStart; i < E_cs_Info->iEnd - 1; i++, sqrtEi++, SEi++ ) {    /* Calculate the contribution from sigma+. */
        if( *sqrtEi > Cutoff ) break;
        Ei = *(E_cs_in++);
        csi = *(E_cs_in++);
        a = K[4] * ( InvSqrtE * sqrtEi[0] + 1. );
        b = K[4] * ( InvSqrtE * sqrtEi[1] + 1. );
        crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default, a, b, dxnerf );
        dcs = *SEi * E * K[3] * dxnerf[4] 
            - 4. * *SEi * E * K[2] * dxnerf[3]
            +      ( csi + *SEi * ( 6. * E - Ei ) ) * K[1] * dxnerf[2]
            - 2. * ( csi + *SEi * ( 2. * E - Ei ) ) * K[0] * dxnerf[1]
            +      ( csi + *SEi * (      E - Ei ) ) *        dxnerf[0];
        cs_p -= dcs;
    }

    Cutoff = ( 1. - vCutoffRatio * K[0] ) * SqrtE;
    for( i =  E_cs_Info->iStart, sqrtEi = &(E_cs_Info->sqrtEi[E_cs_Info->iStart]); i < E_cs_Info->iEnd; i++, sqrtEi++ ) if( *sqrtEi > Cutoff ) break;
    if( i > E_cs_Info->iStart ) i--;
    sqrtEi = &(E_cs_Info->sqrtEi[i]);
    SEi = &(E_cs_Info->SEi[i]);
    E_cs_in = &(E_cs_Info->E_cs_in[2 * i]);
    Cutoff = ( 1. + vCutoffRatio * K[0] ) * SqrtE;
    for( ; i < E_cs_Info->iEnd - 1; i++, sqrtEi++, SEi++ ) {
        if( *sqrtEi > Cutoff ) {
            if( doAgain == 0 ) break;           /* doAgain is needed to solve issue that can arise with E is below threshold. */
            doAgain = 0;
        }
        Ei = *(E_cs_in++);
        csi = *(E_cs_in++);
        a = K[4] * ( InvSqrtE * sqrtEi[0] - 1. );
        b = K[4] * ( InvSqrtE * sqrtEi[1] - 1. );
        pm = 1.;
        if( a < 0. ) {
            if( b < 0. ) {
                dcs = -a;
                a = -b;
                b = dcs;
                pm = -1.; }
            else {
                crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default, 0., -a, dxnerf );
                dcs = *SEi * E * K[3] * dxnerf[4]
                    - 4. * *SEi * E * K[2] * dxnerf[3]
                    +      ( csi + *SEi * ( 6. * E - Ei ) ) * K[1] * dxnerf[2]
                    - 2. * ( csi + *SEi * ( 2. * E - Ei ) ) * K[0] * dxnerf[1]
                    +      ( csi + *SEi * (      E - Ei ) ) *        dxnerf[0];
                cs_m += dcs;
                a = 0.;
            }
        }
        crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default, a, b, dxnerf );
        dcs = *SEi * E * K[3] * dxnerf[4] 
            + pm * 4. * *SEi * E * K[2] * dxnerf[3]
            +           ( csi + *SEi * ( 6. * E - Ei ) ) * K[1] * dxnerf[2]
            + pm * 2. * ( csi + *SEi * ( 2. * E - Ei ) ) * K[0] * dxnerf[1]
            +           ( csi + *SEi * (      E - Ei ) ) *        dxnerf[0];
        cs_m += dcs;
        if( fabs( cs_m * (E_cs_Info->info)->fInterpolation ) < fabs( dcs ) ) doAgain = 1;
    }
    if( *sqrtEi < Cutoff ) {
        switch( E_cs_Info->upperlimit ) {
        case crossSectionAdjustForHeatedTarget_limit_constant :
        case crossSectionAdjustForHeatedTarget_limit_one_over_v :
            Ei = *(E_cs_in++);
            csi = *(E_cs_in++);
            a = K[4] * ( InvSqrtE * sqrtEi[0] - 1. );
            b = K[4] * ( InvSqrtE * sqrtEi[0] + 1. );
            crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default, a, b, dxnerf );
            pm = K[0] * myExp( -b * b );
            switch( E_cs_Info->upperlimit ) {
            case crossSectionAdjustForHeatedTarget_limit_constant :
                a = 1.;
                b = 2.;
                scale = 1.;
                break;
            default :
                a = 0.;
                b = 1.;
                scale = sqrt( Ei / E );                 /* v_N / v_n */
            }
            dcs = scale * csi * ( a * K[1] * dxnerf[2] + b * K[0] * dxnerf[1] + dxnerf[0] + b * pm );
            cs_l += dcs;
            break;
        case crossSectionAdjustForHeatedTarget_limit_threshold :
        default :
            break;
        }
    }

    cs = cs_p + cs_m + cs_l;
    E_cs_point->E = E;
    if( cs < xsecMin2 * E_cs_Info->csMax ) cs = 0.;
    E_cs_point->cs = cs;
    if( E_cs_Info->csMax < cs ) E_cs_Info->csMax = cs;
}
#if 0
/*
****************************************************
*/
static double crossSectionAdjustForHeatedTarget_sum_factors( double pm, double E, double Ei, double csi, double SEi, double K[5], double dxnerf[5] ) {

    double dcs = 0, dcs1, dcs2;

    dcs1 = ( ( ( K[0] * dxnerf[4] + pm * 4. * dxnerf[3] ) * K[0] + 6. * dxnerf[2] ) * K[0] + pm * 4. * dxnerf[1] ) * K[0] + dxnerf[0];
    dcs2 = ( dxnerf[2] * K[0] + pm * 2. * dxnerf[1] ) * K[0] + dxnerf[0];
    dcs = SEi * ( E * dcs1 - Ei * dcs2 ) + csi * dcs2;

    dcs1 = E * K[3] * dxnerf[4] + pm * 4. * E * K[2] * dxnerf[3] + ( 6. * E - Ei ) * K[1] * dxnerf[2]
            + pm * 2. * ( 2. * E - Ei ) * K[0] * dxnerf[1] + ( E - Ei ) * dxnerf[0];
    dcs = SEi * dcs1 + csi * dcs2;
    return( dcs );
}
#endif
/*
****************************************************
*/
int crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode Mode , double a, double b, double dxnerf[5] ) {
/*
*   Returns the integrals dxnerf[i] = int_a^b x^i e^(-x^2) dx for 0 <= i <= 4.  The parameter Mode is
*   intended for debugging, and the value crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default is intended for normal use.
*/
    double x, x2, eps, eps2, eps3, eps5, eps7, exp_two_sqrtpi;
    double a2, b2, expa2, expb2;

    x = ( a + b ) / 2.;
    x2 = x * x;

    eps = ( b - a ) / 2.;    /* Actually, eps / 2. */
    eps2 = eps  * eps;

    a2 = a * a;
    b2 = b * b;
    if( ( eps < 10. ) && ( x < 25. ) ) {
        dxnerf[1] = myExp( -x2 - eps2 ) * mySinh( 2. * eps * x ) / sqrtpi; }
    else {
        dxnerf[1] = ( myExp( -a2 ) - myExp( -b2 ) ) / ( 2. * sqrtpi );
    }
    if( ( Mode == crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_taylor ) || 
        ( ( Mode == crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default ) && ( eps < 0.001 ) ) ) {
        eps3 = eps  * eps2;
        eps5 = eps3 * eps2;
        eps7 = eps5 * eps2;

        exp_two_sqrtpi = myExp( -x2 ) * ( 2. / sqrtpi );

        dxnerf[0] = exp_two_sqrtpi * ( eps + ( 2. * x2 - 1. ) / 3. * eps3 + ( 4. * ( x2 - 3 ) * x2 + 3. ) / 30. * eps5 +
            ( ( ( 8. * x2 - 60. ) * x2 + 90. ) * x2 - 15. ) / 630. * eps7 );
        dxnerf[2] = exp_two_sqrtpi * ( x2 * eps + ( ( 2. * x2 - 5. ) * x2 + 1. ) / 3. * eps3 +
            ( ( ( 4. * x2 - 28. ) * x2 + 39. ) - 6. ) / 30. * eps5 + 
            ( ( ( ( 8. * x2 - 108. ) * x2 + 390. ) * x2 - 375. ) * x2 + 45 ) / 630. * eps7 );
        dxnerf[3] = exp_two_sqrtpi * x * ( x2 * eps + ( ( 2. * x2 - 7. ) * x2 + 3. ) / 3. * eps3 +
            ( ( ( 4. * x2 - 36. ) * x2 + 75. ) - 30. ) / 30. * eps5 + 
            ( ( ( ( 8. * x2 - 132. ) * x2 + 630. ) * x2 - 945. ) * x2 + 315 ) / 630. * eps7 );
        dxnerf[4] = exp_two_sqrtpi * ( x2 * x2 * eps + x2 * ( ( 2. * x2 - 9. ) * x2 + 6. ) / 3. * eps3 +
            ( ( ( ( 4. * x2 - 44. ) * x2 + 123. ) - 84. ) * x2 + 6. ) / 30. * eps5 +
            ( ( ( ( ( 8. * x2 - 156. ) * x2 + 930. ) * x2 - 1935. ) * x2 + 1170 ) * x2 - 90. ) / 630. * eps7 ); }
    else {
        expa2 = myExp( - a2 ) / ( 2. * sqrtpi );
        expb2 = myExp( - b2 ) / ( 2. * sqrtpi );
        if( a < 2. ) {
            dxnerf[0] = 0.5 * ( myErf( b ) - myErf( a ) ); }
        else {
            dxnerf[0] = 0.5 * ( myErfc( a ) - myErfc( b ) );
        }
        dxnerf[2] = 0.5 * dxnerf[0] + a * expa2 - b * expb2;
        dxnerf[3] = ( a2 + 1. ) * expa2 - ( b2 + 1. ) * expb2;
        dxnerf[4] = 0.75 * dxnerf[0] + 0.5 * ( a * ( 2. * a2 + 3. ) * expa2 - b * ( 2. * b2 + 3. ) * expb2 );
    }
    return( 0 );
}
/*
****************************************************
*/
static double oneOverVStepSize( double f ) {

    int i;
    double xMin = 1.01, xMax = 10., x, xMid;

    for( i = 0; i < 20; i++ ) {
        xMid = sqrt( xMax * xMin );
        x = pow( xMid * sqrt( 0.5 * ( 1 + xMid ) ), 1. / 3. );
        if( ( sqrt( x ) * ( ( 1. - x ) / ( ( sqrt( xMid ) + 1. ) * sqrt( xMid ) ) + 1. ) - 1. ) < f ) {
            xMin = xMid; }
        else {
            xMax = xMid;
        }
    }
    xMid = sqrt( xMax * xMin );
    return( xMid );
}
/*
****************************************************
*/
static double myExp( double x ) {

    if( x < -xMaxForExp ) {
        x = -xMaxForExp; }
    else if( x > xMaxForExp ) {
        x = xMaxForExp;
    }
    return( exp( x ) );
}
/*
****************************************************
*/
static double myErf( double x ) {

    if( xMaxForErfc < 0 ) xMaxForErfc = sqrt( xMaxForExp );
    if( x < -xMaxForErfc ) {
        x = -xMaxForErfc; }
    else if( x > xMaxForErfc ) {
        x = xMaxForErfc;
    }
    return( erf( x ) );
}
/*
****************************************************
*/
static double myErfc( double x ) {

    if( xMaxForErfc < 0 ) xMaxForErfc = sqrt( xMaxForExp );
    if( x < -xMaxForErfc ) {
        x = -xMaxForErfc ; }
    else if( x > xMaxForErfc ) {
        x = xMaxForErfc;
    }
    return( erfc( x ) );
}
/*
****************************************************
*/
static double mySinh( double x ) {

    if( x < -xMaxForExp ) {
        x = -xMaxForExp; }
    else if( x > xMaxForExp ) {
        x = xMaxForExp;
    }
    return( sinh( x ) );
}
/*
****************************************************
*/
static double crossSectionAdjustForHeatedTarget_interpolate_cs( double E, int index, double *E_cs_in ) {

    double *p = &(E_cs_in[index]), dE = p[2] - p[0], cs;

    if( dE == 0. ) {
        cs = 0.5 * ( p[1] + p[3] ); }
    else {
        cs = ( p[1] * ( p[2] - E ) + p[3] * ( E - p[0] ) ) / dE;
    }
    return( cs );
}
