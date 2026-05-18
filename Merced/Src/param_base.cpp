/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: param_base.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle uncorrelated energy-angle distributions

#include <cmath>

#include "param_base.hpp"
#include "messaging.hpp"
#include "global_params.hpp"


// ************* class Pbase::param_base *****************
// ----------- Pbase::param_base::setup --------------
// Sets up the quadrature parameters
void Pbase::param_base::setup( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
  const Ddvec::dd_vector& weight_,
  const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups )
{
  // save the interpolation type for the model weight
  weight_interp = weight_.interp_type;
  
  // pointers to the cross section
  first_ladder_sigma = sigma_.begin( );
  sigma_end = sigma_.end( );
  if( first_ladder_sigma->y == 0.0 )
  {
    next_sigma = first_ladder_sigma;
    for( ++next_sigma; next_sigma != sigma_end;
	 first_ladder_sigma = next_sigma, ++next_sigma )
    {
      if( next_sigma->y != 0.0 )
      {
	break;
      }
    }
  }

  // pointers to the multiplicity
  this_mult = mult_.begin( );
  next_mult = this_mult;
  ++next_mult;
  mult_end = mult_.end( );

  // pointers to the model weight
  this_weight = weight_.begin( );
  next_weight = this_weight;
  ++next_weight;
  weight_end = weight_.end( );

  // pointers to the flux data
  flux_ptr = e_flux_.begin( );
  next_flux = flux_ptr;
  ++next_flux;
  flux_end = e_flux_.end( );
  order = flux_ptr->order;
  current_weight.initialize( order );

  // pointers to the incident energy boundaries
  Ein_ptr = Ein_groups.begin( );
  next_Ein = Ein_ptr;
  ++next_Ein;
  Ein_end = Ein_groups.end( );
  Ein_count = 0;

  // Set the first common incident energy
  common_E0( );
}
// ----------- Pbase::param_base::setup_bin --------------
// Sets up the quadrature parameters for integration on a bin
void Pbase::param_base::setup_bin( int Ein_bin, const Ddvec::dd_vector& sigma_,
  const Ddvec::dd_vector& mult_, const Ddvec::dd_vector& weight_,
  const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups )
{
  // the incident energy bin
  Ein_count = Ein_bin;
  Ein_ptr = Ein_groups.begin( ) + Ein_count;
  next_Ein = Ein_ptr;
  ++next_Ein;

  // save the interpolation type for the model weight
  weight_interp = weight_.interp_type;

  // pointers to the cross section
  first_ladder_sigma = sigma_.begin( );
  next_sigma = first_ladder_sigma;
  ++next_sigma;
  sigma_end = sigma_.end( );
  while( ( next_sigma != sigma_end ) && ( next_sigma->x <= *Ein_ptr ) )
  {
    first_ladder_sigma = next_sigma;
    ++next_sigma;
  }

  // pointers to the multiplicity
  this_mult = mult_.begin( );
  next_mult = this_mult;
  ++next_mult;
  mult_end = mult_.end( );
  while( ( next_mult->x <= *Ein_ptr ) && ( next_mult != mult_end ) )
  {
    this_mult = next_mult;
    ++next_mult;
  }

  // pointers to the model weight
  this_weight = weight_.begin( );
  next_weight = this_weight;
  ++next_weight;
  weight_end = weight_.end( );
  while( ( next_weight->x <= *Ein_ptr ) && ( next_weight != weight_end ) )
  {
    this_weight = next_weight;
    ++next_weight;
  }

  // pointers to the flux data
  flux_ptr = e_flux_.begin( );
  next_flux = flux_ptr;
  ++next_flux;
  flux_end = e_flux_.end( );
  while( ( next_flux->get_E_in( ) <= *Ein_ptr ) && ( next_flux != flux_end ) )
  {
    flux_ptr = next_flux;
    ++next_flux;
  }

  order = flux_ptr->order;
  current_weight.initialize( order );

  // Set the first common incident energy
  common_E0( );
}
// ----------- Pbase::param_base::common_E0 --------------
// Sets the first common incident energy
void Pbase::param_base::common_E0( )
{
  data_E_0 = ( first_ladder_sigma->x < *Ein_ptr ) ? *Ein_ptr :
    first_ladder_sigma->x;
  if( this_mult->x > data_E_0 ) data_E_0 = this_mult->x;
  if( this_weight->x > data_E_0 ) data_E_0 = this_weight->x;
  if( flux_ptr->get_E_in( ) > data_E_0 ) data_E_0 = flux_ptr->get_E_in( );
}
// ----------- Pbase::param_base::set_Ein_range --------------
// Sets the range of interpolagtion over incident energy
void Pbase::param_base::set_Ein_range( )
{
  data_E_0 = ( this_mult->x < *Ein_ptr ) ? *Ein_ptr :
    this_mult->x;
  if( this_weight->x > data_E_0 ) data_E_0 = this_weight->x;
  if( flux_ptr->get_E_in( ) > data_E_0 ) data_E_0 = flux_ptr->get_E_in( );

  // Ensure that energies are in order
  bool data_bad = update_pointers( data_E_0 );
  if( data_bad )
  {
    Msg::FatalError( "Pbase::param_base::set_Ein_range",
		     "energies inconsistent" );
  }

  data_E_1 = ( next_mult->x > *next_Ein ) ? *next_Ein :
    next_mult->x;
  if( next_weight->x < data_E_1 ) data_E_1 = next_weight->x;
  if( next_flux->get_E_in( ) < data_E_1 ) data_E_1 = next_flux->get_E_in( );
}
// ----------- Pbase::param_base::set_sigma_range --------------
// Sets the range of pointers to cross sections for this set of data
void Pbase::param_base::set_sigma_range( )
{
  static double skip_tol = Global.Value( "tight_tol" );

  Ddvec::dd_vector::const_iterator sigma_ptr = first_ladder_sigma;
  while( sigma_ptr->x < data_E_0 * ( 1.0 + skip_tol ) )
  {
    ++sigma_ptr;
    if( sigma_ptr == sigma_end )
    {
      Msg::FatalError( "Pbase::param_base::set_sigma_range",
		       "we ran out of cross section data" );
    }
  }
  if( first_ladder_sigma != sigma_ptr )
  {
    first_ladder_sigma = sigma_ptr;
    --first_ladder_sigma;
  }

  while( sigma_ptr->x < data_E_1 )
  {
    ++sigma_ptr;
    if( sigma_ptr == sigma_end )
    {
      break;
    }
  }
  last_ladder_sigma = sigma_ptr;
}
// ----------- Pbase::param_base::get_Ein_range --------------
// Gets the range of nontrivial incident energy bins
// returns true if the threshold is too high for the energy bins
bool Pbase::param_base::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
  const Ddvec::dd_vector& weight,
  const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups,
  double *first_Ein, double *last_Ein )
{
  bool done = false;
  // The cross section sometimes starts with several zeros.
  Ddvec::dd_vector::const_iterator this_sigma = sigma.begin( );
  Ddvec::dd_vector::const_iterator next_sigma = this_sigma;
  ++next_sigma;
  if( this_sigma->y == 0.0 )
  {
    for( ; next_sigma != sigma_end;
	 this_sigma = next_sigma, ++next_sigma )
    {
      if( next_sigma->y != 0.0 )
      {
	break;
      }
    }
  }
  double E_sigma = this_sigma->x;
  double E_mult = mult.begin( )->x;
  double E_weight = weight.begin( )->x;
  double E_flux = e_flux.begin( )->get_E_in( );

  double E_first = ( E_sigma < E_mult ) ? E_mult : E_sigma;
  if( E_weight > E_first ) E_first = E_weight;
  if( E_flux > E_first ) E_first = E_flux;

  double top_Ein = Ein_groups.get_top_E( );
  *last_Ein = top_Ein;
  if( E_first >= top_Ein )
  {
    *first_Ein = 0.0;
    *last_Ein = 0.0;
    return true;
  }
  *first_Ein = E_first;

  Ddvec::dd_vector::const_iterator sigma_ptr = sigma.end( );
  --sigma_ptr;
  E_sigma = sigma_ptr->x;
  Ddvec::dd_vector::const_iterator mult_ptr = mult.end( );
  --mult_ptr;
  E_mult = mult_ptr->x;
  Ddvec::dd_vector::const_iterator weight_ptr = weight.end( );
  --weight_ptr;
  E_weight = weight_ptr->x;
  Lgdata::Flux_List::const_iterator Flux_ptr = e_flux.end( );
  --Flux_ptr;
  E_flux = Flux_ptr->get_E_in( );

  double E_last = ( E_sigma > E_mult ) ? E_mult : E_sigma;
  if( E_weight < E_last ) E_last = E_weight;
  if( E_flux < E_last ) E_last = E_flux;
  *last_Ein = E_last;
  return done;
}
// ----------- Pbase::param_base::update_pointers --------------
// Sets the data pointers for a new incident energy interval.
// Returns "true" at the end of the data
bool Pbase::param_base::update_pointers( double E_in )
{
  static double etol = Global.Value( "tight_tol" );
  double E_tol = E_in * etol;
  //  double E_tol = 0.0;
    while( E_in + E_tol >= next_mult->x )
    {
      // increment multiplicity
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == mult_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_weight->x )
    {
      // increment model weight
      this_weight = next_weight;
      ++next_weight;
      if( next_weight == weight_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= *next_Ein )
    {
      // go to the next incident energy bin
      ++Ein_count;
      Ein_ptr = next_Ein;
      ++next_Ein;
      if( next_Ein == Ein_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_flux->get_E_in( ) )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == flux_end )
      {
        return true;
      }
    }

  return false;
}
// ----------- Pbase::param_base::update_bin_pointers --------------
// Sets the data pointers for a new incident energy interval in one bin.
// Returns "true" at the end of the data on this bin.
bool Pbase::param_base::update_bin_pointers( double E_in )
{
  static double etol = Global.Value( "tight_tol" );
  double E_tol = E_in * etol;
  //  double E_tol = 0.0;
    if( E_in + E_tol >= *next_Ein )
    {
      return true;
    }
    while( E_in + E_tol >= next_mult->x )
    {
      // increment multiplicity
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == mult_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_weight->x )
    {
      // increment model weight
      this_weight = next_weight;
      ++next_weight;
      if( next_weight == weight_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_flux->get_E_in( ) )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == flux_end )
      {
        return true;
      }
    }

  return false;
}
// ----------- Pbase::param_base::set_weight --------------
// Calculates the weight (cross section) * flux * multiplicity * model weight
void Pbase::param_base::set_weight( double E_in )
{
  // Ignore this bool
  bool is_OK;

  // interpolate the weight
  double vikt;
  if( weight_interp == Terp::HISTOGRAM )
  {
    vikt = this_weight->y;
  }
  else
  {
    vikt = this_weight->linlin_interp( E_in, *next_weight, &is_OK );
  }
  
  is_OK = current_weight.linlin_interp( E_in, *flux_ptr, *next_flux );
  double scale_by = this_sigma->linlin_interp( E_in, *next_sigma, &is_OK ) *
    this_mult->linlin_interp( E_in, *next_mult, &is_OK ) * vikt;
  current_weight *= scale_by;
}
// ----------- Pbase::param_base::flux_weight --------------
// Calculates the weight for gammas: flux * multiplicity * model weight
void Pbase::param_base::flux_weight( double E_in )
{
  // Ignore this bool
  bool is_OK;

  is_OK = current_weight.linlin_interp( E_in, *flux_ptr, *next_flux );
  double scale_by = this_mult->linlin_interp( E_in, *next_mult, &is_OK ) * this_weight->y;
  current_weight *= scale_by;
}
