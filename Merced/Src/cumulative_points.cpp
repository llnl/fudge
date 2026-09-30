/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-28 (Fri, Jan 28, 2011) $
 * $Author: hedstrom $
 * $Id: cumulative_points.cpp 1 2011-01-28 hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// classes used to handle interpolation by cumulative points

#include <cmath>
#include <cstdlib>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "cumulative_points.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class Cum::cumulative_prob_entry *****************
// ----------- Cum::cumulative_prob_entry::get_prob --------------
// Gets the probability density at outgoing energy E_in
double Cum::cumulative_prob_entry::get_prob( double E ) const
{
  double this_prob = Prob + slope*( E - E_out );
  return this_prob;
}
// ----------- Cum::cumulative_prob_entry::get_cum_prob --------------
// Gets the cumulative probability at outgoing energy E_in
double Cum::cumulative_prob_entry::get_cum_prob( double E ) const
{
  double dE = E - E_out;
  double this_cum = cum_prob + dE*( Prob + 0.5*slope*dE );
  return this_cum;
}
// ----------- Cum::cumulative_prob_entry::get_cum_inv --------------
// Gets the energy corresponding to cumulative probability A
double Cum::cumulative_prob_entry::get_cum_inv( double A ) const
{
  double dA = A - cum_prob;
  double this_E;
  
  static double abs_tol = Global.Value( "tight_tol" );
  if( dA <= abs_tol )
  {
    this_E = E_out;
  }
  else if( slope == 0.0 )
  {
    this_E = E_out + dA/ Prob;
  }
  else
  {
    double root = Prob*Prob + 2.0*slope*dA;
    if( root >= 0.0 )
    {
      this_E = E_out + 2.0*dA/( Prob + std::sqrt( root ) );
    }
      else
    {
      this_E = E_out - Prob/slope;
    }
  }
  
  return this_E;
}
// ----------- Cum::cumulative_prob_entry::set_slope --------------
// Sets the slope for this interval
void Cum::cumulative_prob_entry::set_slope( double next_E, double next_A )
{
  double dE = next_E - E_out;
  if( dE <= 0.0 )
  {
    Msg::FatalError( "Cum::cumulative_prob_entry::set_slope",
		"improper energy interval" );
  }
  
  slope = 2.0 * ( next_A - cum_prob - dE* Prob ) / ( dE * dE );
}
// ----------- Cum::cumulative_prob_entry::CP_interpolate --------------
// Interpolates between two cumulative_prob_entrys
void Cum::cumulative_prob_entry::CP_interpolate( double alpha, double this_A,
		       const Cum::cumulative_prob_entry &prev_CP,
		       const Cum::cumulative_prob_entry &next_CP )
{
  double prev_E = prev_CP.get_cum_inv( this_A );
  double prev_prob = prev_CP.get_prob( prev_E );

  double next_E = next_CP.get_cum_inv( this_A );
  double next_prob = next_CP.get_prob( next_E );

  // interpolate
  cum_prob = this_A;
  E_out = ( 1.0 - alpha ) * prev_E + alpha * next_E;
  Prob = ( 1.0 - alpha ) * prev_prob + alpha * next_prob;

  // the slope is calculated later
  
}
// ----------- Cum::cumulative_prob_entry::print --------------
// Prints the entry
void Cum::cumulative_prob_entry::print( )
{
  std::cout << "E: " << E_out << " P: " << Prob <<
    " S: " << slope << " A: " << cum_prob << std::endl;
}

// ************* class Cum::cumulative_prob_list *****************
// ----------- Cum::cumulative_prob_list::get_cum_prob_flat --------------
// Computes the cumulative probabilities and slopes for histogram data
void Cum::cumulative_prob_list::get_cum_prob_flat( )
{
  if( Eout_interp != Terp::HISTOGRAM )
  {
    Msg::FatalError( "Cum::cumulative_prob_list::get_cum_prob_flat",
		"wrong interpolation type" );
  }

  // pointers to the data
  static bool zero_found = false;
  Cum::cumulative_prob_list::iterator this_link = begin( );
  Cum::cumulative_prob_list::iterator next_link = this_link;
  ++next_link;

  // set up the list of cumulative probabilities
  double sum = 0.0;
  this_link->cum_prob = 0.0;
  this_link->slope = 0.0;

  for( ; next_link != end( ); this_link = next_link, ++next_link )
  {
    // norm for histogram data
    double to_add = this_link->Prob * ( next_link->E_out - this_link->E_out );
    if( ( to_add <= 0.0 ) && ( !zero_found ) )
    {
      Msg::Warning( "Cum::cumulative_prob_list::get_cum_prob_flat",
	       "data have an interval with zero probability" );
      zero_found = true;
    }
    sum += to_add;
    next_link->cum_prob = sum;
    next_link->slope = 0.0;
  }
  // ensure that cum_prob <= 1
  for( this_link = begin( ); this_link != end( ); ++this_link )
  {
    if( this_link->cum_prob > 1.0 )
    {
      this_link->cum_prob = 1.0;
    }
  }
  this_link = end( );
  --this_link;
  this_link->cum_prob = 1.0;
}
// ----------- Cum::cumulative_prob_list::get_cum_prob_linlin --------------
// Computes the cumulative probabilities and slopes for lin-lin data
void Cum::cumulative_prob_list::get_cum_prob_linlin( )
{
  if( Eout_interp != Terp::LINLIN )
  {
    Msg::FatalError( "Cum::cumulative_prob_list::get_cum_prob_linlin",
		"wrong interpolation type" );
  }

  // pointers to the data
  static bool zero_found = false;
  Cum::cumulative_prob_list::iterator this_link = begin( );
  Cum::cumulative_prob_list::iterator next_link = this_link;
  ++next_link;
  double E0 = this_link->E_out;
  double P0 = this_link->Prob;

  // set up the list of cumulative probabilities
  double sum = 0.0;
  this_link->cum_prob = 0.0;

  for( ; next_link != end( ); this_link = next_link, ++next_link )
  {
    double E1 = next_link->E_out;
    double P1 = next_link->Prob;
    double dE = E1 - E0;
    // norm for linlin data
    double to_add = 0.5*dE*( P1 + P0 );
    if( ( to_add <= 0.0 ) && ( !zero_found ) )
    {
      Msg::Warning( "Cum::cumulative_prob_list::get_cum_prob_linlin",
	       "data have an interval with zero probability" );
      zero_found = true;
    }
    sum += to_add;
    next_link->cum_prob = sum;
    if( dE == 0.0 )
    {
      this_link->slope = 0.0;
    }
    else
    {
      this_link->slope = ( P1 - P0 )/dE;
    }
    E0 = E1;
    P0 = P1;
  }
  // ensure that cum_prob <= 1
  for( this_link = begin( ); this_link != end( ); ++this_link )
  {
    if( this_link->cum_prob > 1.0 )
    {
      this_link->cum_prob = 1.0;
    }
  }
  this_link = end( );
  --this_link;
  this_link->cum_prob = 1.0;
  this_link->slope = 0.0;
}
// ----------- Cum::cumulative_prob_list::set_slopes --------------
// Sets the slopes for all of the entries
void Cum::cumulative_prob_list::set_slopes( )
{
  Cum::cumulative_prob_list::iterator this_entry = begin( );
  Cum::cumulative_prob_list::iterator next_entry = this_entry;
  ++next_entry;

  for( ; next_entry != end( ); this_entry = next_entry, ++next_entry )
  {
    this_entry->set_slope( next_entry->E_out, next_entry->cum_prob );
  }

  // the last slope is zero
  this_entry->slope = 0.0;
}
