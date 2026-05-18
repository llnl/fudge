/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-28 (Fri, Jan 28, 2011) $
 * $Author: hedstrom $
 * $Id: standard_Legendre.cpp 1 2011-01-28 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle Legendre expansions of energy probability density

#include <cmath>
#include <cstdlib>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "standard_Legendre.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class StdLg::standard_Legendre_vector *****************
// ----------- StdLg::standard_Legendre_vector::read_coef --------------
// Reads the Legendre coefficients of the energy probability density
void StdLg::standard_Legendre_vector::read_coef( Dpar::data_parser& infile, int num_Eout, int max_order )
{
  StdLg::standard_Legendre_vector::iterator new_Eout_ptr;
  // read the data
  for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
  {
    // make a new set of Legendre coefficients
    new_Eout_ptr = insert( end( ), Lgdata::Legendre_coefs( ) );
    new_Eout_ptr->initialize( max_order );
    new_Eout_ptr->set_E_out( infile.get_next_double( ) );  // energy of outgoing particle
    int file_order = infile.get_next_int( ) - 1;  // Legendre order of input data
    int save_coefs = ( file_order > max_order ) ? max_order : file_order;
    int coef_count = 0;
    for( ; coef_count <= save_coefs; ++coef_count )
    {
      ( *new_Eout_ptr )[ coef_count ] = infile.get_next_double( );
    }
    // we may need to discard high-order coefficients
    for( coef_count = save_coefs+1; coef_count <= file_order; ++coef_count )
    {
      infile.get_next_double( );
    }
  }
  // ensure proper normalization, data not truncated
  renorm( false );

  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    to_unit_base( );
  }
  else if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    form_cum_prob( );
  }
}
// ----------- StdLg::standard_Legendre_vector::append_data --------------
// Appends a copy of the data to the list
void StdLg::standard_Legendre_vector::append_data( double E_out, const Lgdata::Legendre_coefs &to_copy )
{
  // make a new set of Legendre coefficients
  StdLg::standard_Legendre_vector::iterator new_Eout_ptr;
  new_Eout_ptr = insert( end( ), Lgdata::Legendre_coefs( ) );
  new_Eout_ptr->initialize( to_copy.order );
  new_Eout_ptr->set_E_out( E_out );
  for( int L_count = 0; L_count <= to_copy.order; ++L_count )
  {
    (*new_Eout_ptr)[ L_count ] = to_copy.value( L_count );
  }
}
// ----------- StdLg::standard_Legendre_vector::copy --------------
// Copies a vector
void StdLg::standard_Legendre_vector::copy( const StdLg::standard_Legendre_vector& vector_from )
{
  set_E_in( vector_from.get_E_in( ) );
  Ein_interp = vector_from.Ein_interp;
  Eout_interp = vector_from.Eout_interp;

  for( StdLg::standard_Legendre_vector::const_iterator entry = vector_from.begin( );
       entry != vector_from.end( ); ++entry )
  {
    append_data( entry->get_E_out( ), *entry );
  }
}
// ----------- StdLg::standard_Legendre_vector::truncate_copy --------------
// Copies a vector with truncation
bool StdLg::standard_Legendre_vector::truncate_copy(
     const StdLg::standard_Legendre_vector& vector_from,
     double min_E, double max_E )
{
  if( min_E >= max_E )
  {
    return false;
  }
  
  set_E_in( vector_from.get_E_in( ) );
  Ein_interp = vector_from.Ein_interp;
  Eout_interp = vector_from.Eout_interp;

  static double abs_tol = Global.Value( "tight_tol" );
  StdLg::standard_Legendre_vector::const_iterator prev_entry = vector_from.begin( );
  StdLg::standard_Legendre_vector::const_iterator next_entry = prev_entry;
  ++next_entry;

  Lgdata::Legendre_coefs new_entry;  // for creating entries
  int prev_order = prev_entry->order;
  int next_order = next_entry->order;
  int use_order = ( prev_order > next_order ) ?
    prev_order : next_order;
  new_entry.initialize( use_order );

  // find the start
  while( next_entry->get_E_out( ) < ( 1.0 + abs_tol ) * min_E )
  {
    prev_entry = next_entry;
    ++next_entry;
  }
  if( Eout_interp == Terp::HISTOGRAM )
  {
    append_data( min_E, *prev_entry );
  }
  else  // LTerp::INLIN
  {
    if( prev_entry->get_E_out( ) < ( 1.0 - abs_tol ) * min_E )
    {
      new_entry.linlin_interp( min_E, *prev_entry, *next_entry );
      append_data( min_E, new_entry );
    }
    else
    {
      append_data( min_E, *prev_entry );
    }
  }

  // do the middle ones
  while( next_entry->get_E_out( ) < ( 1.0 - abs_tol ) * max_E )
  {
    append_data( next_entry->get_E_out( ), *next_entry );
    prev_entry = next_entry;
    ++next_entry;
  }

  // do the last one
  if( Eout_interp == Terp::HISTOGRAM )
  {
    append_data( max_E, *prev_entry );
  }
  else  // Terp::LINLIN
  {
    if( prev_entry->get_E_out( ) < ( 1.0 - abs_tol ) * max_E )
    {
      prev_order = prev_entry->order;
      next_order = next_entry->order;
      use_order = ( prev_order > next_order ) ?
          prev_order : next_order;
      new_entry.initialize( use_order );
      new_entry.linlin_interp( max_E, *prev_entry, *next_entry );
      append_data( max_E, new_entry );
    }
    else
    {
      append_data( max_E, *prev_entry );
    }
  }

  // renorm truncated data
  bool norm_OK = renorm( true );

  return norm_OK;
}
// ----------- StdLg::standard_Legendre_vector::extrapolate_copy --------------
// Copies a vector with extrapolation
void StdLg::standard_Legendre_vector::extrapolate_copy(
     const StdLg::standard_Legendre_vector& vector_from,
     double min_E, double max_E )
{
  set_E_in( vector_from.get_E_in( ) );
  Ein_interp = vector_from.Ein_interp;
  Eout_interp = vector_from.Eout_interp;

  Lgdata::Legendre_coefs null_entry;  // for creating entries, initially all zero
  null_entry.initialize( 0 );

  static double abs_tol = Global.Value( "tight_tol" );
  StdLg::standard_Legendre_vector::const_iterator this_entry = vector_from.begin( );

  double E_value = this_entry->get_E_out( );
  if( min_E < ( 1.0 - abs_tol )*E_value )
  {
    // extrapolate a zero head of the list
    append_data( min_E, null_entry );
    if( this_entry->value( 0 ) != 0.0 )
    {
      append_data( ( 1.0 - abs_tol )*E_value, null_entry );
    }
    append_data( E_value, *this_entry );
  }
  else if( min_E < E_value )
  {
    append_data( min_E, *this_entry );
  }
  else
  {
    append_data( E_value, *this_entry );
  }
  for( ++this_entry; this_entry != vector_from.end( ); ++this_entry )
  {
    append_data( this_entry->get_E_out( ), *this_entry );
  }
  StdLg::standard_Legendre_vector::iterator Lptr = end( );
  --Lptr;
  E_value = Lptr->get_E_out( );
  if( max_E > ( 1.0 + abs_tol )*E_value )
  {
    // extrapolate a zero tail
    if( Lptr->value( 0 ) != 0.0 )
    {
      append_data( ( 1.0 + abs_tol )*E_value, null_entry );
    }
    append_data( max_E, null_entry );
  }
  else if( max_E > E_value )
  {
    Lptr->set_E_out( max_E );
  }
}
// ----------- StdLg::standard_Legendre_vector::form_cum_prob --------------
// Forms the list of cumulative probabilities
void StdLg::standard_Legendre_vector::form_cum_prob( )
{
  // copy the data
  cum_prob.Eout_interp = Eout_interp;
  for( StdLg::standard_Legendre_vector::const_iterator Eout_ptr = begin( );
       Eout_ptr != end( ); ++Eout_ptr )
  {
    Cum::cumulative_prob_list::iterator cum_prob_ptr = cum_prob.insert(
      cum_prob.end( ), Cum::cumulative_prob_entry( ) );
    cum_prob_ptr->E_out = Eout_ptr->get_E_out( );
    cum_prob_ptr->Prob = Eout_ptr->value( 0 );
  }
  // now form the slopes and cumulative probabilities
  if( Eout_interp == Terp::HISTOGRAM )
  {
    cum_prob.get_cum_prob_flat( );
  }
  else // lin-lin
  {
    cum_prob.get_cum_prob_linlin( );
  }
}

// *************** class StdLg::standard_Legendre_param *************************
// ------------------ StdLg::standard_Legendre_param::initialize ----------------
// Allocates space
void StdLg::standard_Legendre_param::initialize( int Order )
{
  mid_lower_Eout.initialize( Order );
  mid_upper_Eout.initialize( Order );
  use_prev_Eout.initialize( Order );
  use_next_Eout.initialize( Order );
  Ein0_data.prev_data.initialize( Order );
  Ein0_data.next_data.initialize( Order );
  Ein1_data.prev_data.initialize( Order );
  Ein1_data.next_data.initialize( Order );
}
// ------------------ StdLg::standard_Legendre_param::reset_start ----------------
// Initializes the data pointers for one incident energy range
void StdLg::standard_Legendre_param::reset_start( )
{
  Ein0_data.set_E_in( this_Ein->get_E_in( ) );
  Ein1_data.set_E_in( next_Ein->get_E_in( ) );
  left_Ein = this_Ein->get_E_in( );
  right_Ein = next_Ein->get_E_in( );

  // initialize the pointers
  if( Ein_interp.qualifier == Terp::DIRECT )
  {
    if( Ein_interp.flag == Terp::LINLIN )
    {
      setup_Ein_linlin( );
    }
    else if( Ein_interp.flag == Terp::HISTOGRAM )
    {
      setup_Ein_flat( );
    }
  }
  else if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    setup_Ein_ubase( );
  }
  else if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    setup_Ein_cum_prob( );
  }

  if( Ein_interp.qualifier != Terp::CUMULATIVE_POINTS )
  {
    // for the range of Eout values
    double lower_Eout;
    double higher_Eout;

    lower_Eout = ( left_ptr->get_E_out( ) > right_ptr->get_E_out( ) )?
        left_ptr->get_E_out( ) : right_ptr->get_E_out( );
    higher_Eout = ( next_left_ptr->get_E_out( ) < next_right_ptr->get_E_out( ) )?
        next_left_ptr->get_E_out( ) : next_right_ptr->get_E_out( );
    if( higher_Eout <= lower_Eout )
    {
      Msg::FatalError( "StdLg::standard_Legendre_param::reset_start",
		       "Check the Eout values." );
    }

    // Interpolate to the common Eout values
    common_low_Eout( lower_Eout );
    common_high_Eout( higher_Eout );
  }

  // physical outgoing energy ranges
  if( Ein_interp.qualifier == Terp::DIRECT )
  {
    Eout_0_range.x = Ein0_data.prev_data.get_E_out( );
    Eout_0_range.y = Ein0_data.next_data.get_E_out( );
    Eout_1_range.x = Ein1_data.prev_data.get_E_out( );
    Eout_1_range.y = Ein1_data.next_data.get_E_out( );
  }
  else
  {
    Eout_0_range.x = Ein0_data.ubase_map.Eout_min;
    // Don't test for good interpolation
    Eout_0_range.y = Ein0_data.ubase_map.un_unit_base( Ein0_data.next_data.get_E_out( ) );
    Eout_1_range.x = Ein1_data.ubase_map.Eout_min;
    // Don't test for good interpolation
    Eout_1_range.y = Ein1_data.ubase_map.un_unit_base( Ein1_data.next_data.get_E_out( ) );
  }
}
// ---------------- StdLg::standard_Legendre_param::setup_Ein_ubase ------------------
// Sets up the data for unit-base interpolation in incident energy
void StdLg::standard_Legendre_param::setup_Ein_ubase( )
{
  // lower incident energy
  left_ptr = this_Ein->begin( );
  next_left_ptr = left_ptr;
  ++next_left_ptr;
  last_left_ptr = this_Ein->end( );

  // higher incident energy
  right_ptr = next_Ein->begin( );
  next_right_ptr = right_ptr;
  ++next_right_ptr;
  last_right_ptr = next_Ein->end( );

  // save the unit-base interpolation
  Ein0_data.ubase_map.copy( this_Ein->ubase_map );
  Ein1_data.ubase_map.copy( next_Ein->ubase_map );
}
// ---------------- StdLg::standard_Legendre_param::setup_Ein_cum_prob ------------------
// Sets up the data for cumulative points interpolation in incident energy
void StdLg::standard_Legendre_param::setup_Ein_cum_prob( )
{
  // lower incident energy
  left_ptr = this_Ein->begin( );
  next_left_ptr = left_ptr;
  ++next_left_ptr;
  last_left_ptr = this_Ein->end( );
  left_cum_prob = this_Ein->cum_prob.begin( );
  next_left_cum_prob = left_cum_prob;
  ++next_left_cum_prob;

  // skip zero probability intervals
  while( ( left_cum_prob->Prob == 0.0 ) &&
         ( left_cum_prob->slope == 0.0 ) )
  {
    left_cum_prob = next_left_cum_prob;
    ++next_left_cum_prob;
    left_ptr = next_left_ptr;
    ++next_left_ptr;
  }

  // higher incident energy
  right_ptr = next_Ein->begin( );
  next_right_ptr = right_ptr;
  ++next_right_ptr;
  last_right_ptr = next_Ein->end( );
  right_cum_prob = next_Ein->cum_prob.begin( );
  next_right_cum_prob = right_cum_prob;
  ++next_right_cum_prob;

  // skip zero probability intervals
  while( ( right_cum_prob->Prob == 0.0 ) &&
         ( right_cum_prob->slope == 0.0 ) )
  {
    right_cum_prob = next_right_cum_prob;
    ++next_right_cum_prob;
    right_ptr = next_right_ptr;
    ++next_right_ptr;
  }

  // for the range of cumulative probabilities A
  //  double lower_A = 0.0;
  double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
      next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;

  // set up Ein0_data and Ein1_data
  setup_low_A( );
  setup_high_A( higher_A );

  if( Ein0_data.prev_data.get_E_out( ) < Ein0_data.next_data.get_E_out( ) )
  {
    Ein0_data.to_unit_base( );
  }
  else
  {
    Ein0_data.short_to_unit_base( higher_A );
  }

  if( Ein1_data.prev_data.get_E_out( ) < Ein1_data.next_data.get_E_out( ) )
  {
    Ein1_data.to_unit_base( );
  }
  else
  {
    Ein1_data.short_to_unit_base( higher_A );
  }
}
// ---------------- StdLg::standard_Legendre_param::setup_Ein_linlin ------------------
// Sets up the data for direct linlin interpolation in incident energy
void StdLg::standard_Legendre_param::setup_Ein_linlin( )
{
  // remove previous data
  if( !low_linlin.empty( ) )
  {
    low_linlin.erase( low_linlin.begin( ), low_linlin.end( ) );
    high_linlin.erase( high_linlin.begin( ), high_linlin.end( ) );
  }

  // trancate or extrapolate data?
  static int truncate = Global.Value( "truncate_direct" );
  bool use_truncate = ( truncate > 0 );

  // get the outgoing energy range
  left_ptr = this_Ein->begin( );
  right_ptr = next_Ein->begin( );

  double Eout_min;
  double left_Eout = left_ptr->get_E_out( );
  double right_Eout = right_ptr->get_E_out( );
  if( use_truncate )
  {
    Eout_min = ( left_Eout > right_Eout ) ? left_Eout : right_Eout;
  }
  else
  {
    Eout_min = ( left_Eout < right_Eout ) ? left_Eout : right_Eout;
  }

  left_ptr = this_Ein->end( );
  --left_ptr;
  right_ptr = next_Ein->end( );
  --right_ptr;
  double Eout_max;
  left_Eout = left_ptr->get_E_out( );
  right_Eout = right_ptr->get_E_out( );
  if( use_truncate )
  {
    Eout_max = ( left_Eout < right_Eout ) ? left_Eout : right_Eout;
  }
  else
  {
    Eout_max = ( left_Eout > right_Eout ) ? left_Eout : right_Eout;
  }

  // make copies, truncated or extrapolated
  if( use_truncate )
  {
    bool low_OK = low_linlin.truncate_copy( *this_Ein, Eout_min, Eout_max );
    bool high_OK = high_linlin.truncate_copy( *next_Ein, Eout_min, Eout_max );
    if( !low_OK || !high_OK )
    {
      Msg::Warning( "StdLg::standard_Legendre_param::setup_Ein_linlin",
	       "truncation gave norm 0, using histogram" );
      // we got norm zero; use histogram in incident energy
      if( !low_linlin.empty( ) )
      {
        low_linlin.erase( low_linlin.begin( ), low_linlin.end( ) );
        high_linlin.erase( high_linlin.begin( ), high_linlin.end( ) );
      }

      low_linlin.copy( *this_Ein );
      high_linlin.copy( *this_Ein );
      high_linlin.set_E_in( next_Ein->get_E_in( ) );
    }
  }
  else
  {
    low_linlin.extrapolate_copy( *this_Ein, Eout_min, Eout_max );
    high_linlin.extrapolate_copy( *next_Ein, Eout_min, Eout_max );
  }

  // set pointers at lower incident energy
  left_ptr = low_linlin.begin( );
  next_left_ptr = left_ptr;
  ++next_left_ptr;
  last_left_ptr = low_linlin.end( );

  // higher incident energy
  right_ptr = high_linlin.begin( );
  next_right_ptr = right_ptr;
  ++next_right_ptr;
  last_right_ptr = high_linlin.end( );
}
// ---------------- StdLg::standard_Legendre_param::setup_Ein_flat ------------------
// Sets up the data for histogram interpolation in incident energy
void StdLg::standard_Legendre_param::setup_Ein_flat( )
{
  // lower incident energy
  left_ptr = this_Ein->begin( );
  next_left_ptr = left_ptr;
  ++next_left_ptr;
  last_left_ptr = this_Ein->end( );
  // we need to set up the right pointers
  right_ptr = this_Ein->begin( );
  next_right_ptr = next_left_ptr;
  last_right_ptr = this_Ein->end( );
}
// ------------------ StdLg::standard_Legendre_param::get_next_Eout ----------------
// Increments the data pointers for one incident energy range
bool StdLg::standard_Legendre_param::get_next_Eout( )
{
  // ignore intervals with probability less than skip_tol
  static double skip_tol = Global.Value( "tight_tol" );

  bool done = false;
  if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    // undo the unit-base map before copying data
    // Don't test for good interpolation
    Ein0_data.un_unit_base( );
    Ein1_data.un_unit_base( );
  }

  Ein0_data.prev_data.set_E_out( Ein0_data.next_data.get_E_out( ) );
  Ein0_data.prev_data.copy_coef( Ein0_data.next_data );
  Ein1_data.prev_data.set_E_out( Ein1_data.next_data.get_E_out( ) );
  Ein1_data.prev_data.copy_coef( Ein1_data.next_data );

  // update the pointers
  if( next_left_ptr->get_E_out( ) < Ein0_data.prev_data.get_E_out( ) * 
      ( 1.0 + skip_tol ) )
  {
    left_ptr = next_left_ptr;
    ++next_left_ptr;
    if( next_left_ptr == last_left_ptr )
    {
      return true;
    }
    if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      left_cum_prob = next_left_cum_prob;
      ++next_left_cum_prob;

      bool do_skip = false;
      // skip intervals with essentially zero probability
      while( next_left_cum_prob->cum_prob - left_cum_prob->cum_prob <= skip_tol )
      {
	do_skip = true;
        left_cum_prob = next_left_cum_prob;
        ++next_left_cum_prob;
        left_ptr = next_left_ptr;
        ++next_left_ptr;
        if( next_left_ptr == last_left_ptr )
        {
          return true;
        }
      }
      if( do_skip )
      {
        Ein0_data.prev_data.set_E_out( left_ptr->get_E_out( ) );
        Ein0_data.prev_data.copy_coef( *left_ptr );
      }
    }
  }
  if( next_right_ptr->get_E_out( ) < Ein1_data.prev_data.get_E_out( ) * 
      ( 1.0 + skip_tol ) )
  {
    right_ptr = next_right_ptr;
    ++next_right_ptr;
    if( next_right_ptr == last_right_ptr )
    {
      return true;
    }
    if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      right_cum_prob = next_right_cum_prob;
      ++next_right_cum_prob;

      bool do_skip = false;
      // skip intervals with essentially zero probability
      while( next_right_cum_prob->cum_prob - right_cum_prob->cum_prob <= skip_tol )
      {
	do_skip = true;
        right_cum_prob = next_right_cum_prob;
        ++next_right_cum_prob;
        right_ptr = next_right_ptr;
        ++next_right_ptr;
        if( next_right_ptr == last_right_ptr )
        {
          return true;
        }
      }
      if( do_skip )
      {
        Ein1_data.prev_data.set_E_out( right_ptr->get_E_out( ) );
        Ein1_data.prev_data.copy_coef( *right_ptr );
      }
    }
  }

  if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    // Interpolate to the common higher cumulative probability
    double lower_A = ( left_cum_prob->cum_prob > right_cum_prob->cum_prob )?
      left_cum_prob->cum_prob : right_cum_prob->cum_prob;
    double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
        next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;
    setup_high_A( higher_A );
    double dA = higher_A - lower_A;

    if( Ein0_data.prev_data.get_E_out( ) < Ein0_data.next_data.get_E_out( ) )
    {
      Ein0_data.to_unit_base( );
    }
    else
    {
      Ein0_data.short_to_unit_base( dA );
    }
    
    if( Ein1_data.prev_data.get_E_out( ) < Ein1_data.next_data.get_E_out( ) )
    {
      Ein1_data.to_unit_base( );
    }
    else
    {
      Ein1_data.short_to_unit_base( dA );
    }
  }
  else
  {
    // Interpolate to the common higher Eout value
    double higher_Eout;
    higher_Eout = ( next_left_ptr->get_E_out( ) < next_right_ptr->get_E_out( ) )?
        next_left_ptr->get_E_out( ) : next_right_ptr->get_E_out( );
    common_high_Eout( higher_Eout );
  }

  // Reset the physical E_out ranges
  Eout_0_range.x = Eout_0_range.y;
  Eout_1_range.x = Eout_1_range.y;
  if( ( Ein_interp.qualifier == Terp::UNITBASE ) ||
      ( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) )
  {
    // Don't test for good interpolation
    Eout_0_range.y = Ein0_data.ubase_map.un_unit_base( Ein0_data.next_data.get_E_out( ) );
    Eout_1_range.y = Ein1_data.ubase_map.un_unit_base( Ein1_data.next_data.get_E_out( ) );
    // We may have skipped an interval with zero probability
    Eout_0_range.x = Ein0_data.ubase_map.un_unit_base( Ein0_data.prev_data.get_E_out( ) );
    Eout_1_range.x = Ein1_data.ubase_map.un_unit_base( Ein1_data.prev_data.get_E_out( ) );
  }
  else if( Ein_interp.flag == Terp::LINLIN )
  {
    Eout_0_range.y = Ein0_data.next_data.get_E_out( );
    Eout_1_range.y = Ein1_data.next_data.get_E_out( );
  }
  else // Ein_interp == Terp::HISTOGRAM
  {
    Eout_0_range.y = Ein0_data.next_data.get_E_out( );
    Eout_1_range.y = Eout_0_range.y;
  }
  return done;
}
// ---------------- StdLg::standard_Legendre_param::common_low_Eout ------------------
// Interpolates (Eout, probability) data to the lower common Eout value
void StdLg::standard_Legendre_param::common_low_Eout( double lower_Eout )
{
  Ein0_data.prev_data.set_E_out( lower_Eout );
  Ein1_data.prev_data.set_E_out( lower_Eout );

  if( ( left_ptr->get_E_out( ) == lower_Eout ) || ( Ein0_data.Eout_interp == Terp::HISTOGRAM ) )
  {
    Ein0_data.prev_data.copy_coef( *left_ptr );
  }
  else
  {
    Ein0_data.prev_data.linlin_interp( lower_Eout, *left_ptr, *next_left_ptr );
  }

  if( ( right_ptr->get_E_out( ) == lower_Eout ) || ( Ein1_data.Eout_interp == Terp::HISTOGRAM ) )
  {
    Ein1_data.prev_data.copy_coef( *right_ptr );
  }
  else
  {
    Ein1_data.prev_data.linlin_interp( lower_Eout, *right_ptr, *next_right_ptr );
  }
}
// ---------------- StdLg::standard_Legendre_param::common_high_Eout ------------------
// Interpolates (Eout, probability) data to the higher common Eout value
void StdLg::standard_Legendre_param::common_high_Eout( double higher_Eout )
{
  Ein0_data.next_data.set_E_out( higher_Eout );
  Ein1_data.next_data.set_E_out( higher_Eout );

  static double abs_tol = Global.Value( "tight_tol" );

  if( next_left_ptr->get_E_out( ) < higher_Eout * ( 1 + abs_tol ) )
  {
    Ein0_data.next_data.copy_coef( *next_left_ptr );
  }
  else if( Ein0_data.Eout_interp == Terp::HISTOGRAM )
  {
    Ein0_data.next_data.copy_coef( *left_ptr );
  }
  else
  {
    Ein0_data.next_data.linlin_interp( higher_Eout, *left_ptr, *next_left_ptr );
  }

  if( next_right_ptr->get_E_out( ) < higher_Eout * ( 1 + abs_tol ) )
  {
    Ein1_data.next_data.copy_coef( *next_right_ptr );
  }
  else if( Ein1_data.Eout_interp == Terp::HISTOGRAM )
  {
    Ein1_data.next_data.copy_coef( *right_ptr );
  }
  else
  {
    Ein1_data.next_data.linlin_interp( higher_Eout, *right_ptr, *next_right_ptr );
  }
}
// ---------------- StdLg::standard_Legendre_param::setup_low_A ------------------
// Sets (Eout, probability) data to the lower zero cumulative probability
void StdLg::standard_Legendre_param::setup_low_A( )
{
  Ein0_data.prev_data.set_E_out( left_ptr->get_E_out( ) );
  Ein1_data.prev_data.set_E_out( right_ptr->get_E_out( ) );

  Ein0_data.prev_data.copy_coef( *left_ptr );
  Ein1_data.prev_data.copy_coef( *right_ptr );
}
// ---------------- StdLg::standard_Legendre_param::setup_high_A ------------------
// Interpolates (Eout, probability) data to the higher common cumulative probability
void StdLg::standard_Legendre_param::setup_high_A( double higher_A )
{
  double higher_Eout;

  if( next_left_cum_prob->cum_prob == higher_A )
  {
    Ein0_data.next_data.set_E_out( next_left_ptr->get_E_out( ) );
    Ein0_data.next_data.copy_coef( *next_left_ptr );
  }
  else
  {
    higher_Eout = left_cum_prob->get_cum_inv( higher_A );
    Ein0_data.next_data.set_E_out( higher_Eout );
    if( Ein0_data.Eout_interp == Terp::HISTOGRAM )
    {
      Ein0_data.next_data.copy_coef( *left_ptr );
    }
    else
    {
      Ein0_data.next_data.linlin_interp( higher_Eout, *left_ptr, *next_left_ptr );
    }
  }

  if( next_right_cum_prob->cum_prob == higher_A )
  {
    Ein1_data.next_data.set_E_out( next_right_ptr->get_E_out( ) );
    Ein1_data.next_data.copy_coef( *next_right_ptr );
  }
  else
  {
    higher_Eout = right_cum_prob->get_cum_inv( higher_A );
    Ein1_data.next_data.set_E_out( higher_Eout );
    if( Ein1_data.Eout_interp == Terp::HISTOGRAM )
    {
      Ein1_data.next_data.copy_coef( *right_ptr );
    }
    else
    {
      Ein1_data.next_data.linlin_interp( higher_Eout, *right_ptr, *next_right_ptr );
    }
  }
}
// ------------------ StdLg::standard_Legendre_param::unitbase_interpolate ----------------
// Does unit-base interpolation between two incident energies
bool StdLg::standard_Legendre_param::unitbase_interpolate( double Ein )
{
  bool is_OK = true;
  
  double dEin = right_Ein - left_Ein;
  if( dEin <= 0.0 )
  {
    Msg::DebugInfo( "StdLg::standard_Legendre_param::unitbase_interpolate",
		"Incident energies out of order" );
    return false;
  }
  // lin-lin interpolate the physical energy range
  double alpha = ( Ein - left_Ein )/dEin;
  Eout_range.x = ( 1.0 - alpha )*Eout_0_range.x + alpha*Eout_1_range.x;
  Eout_range.y = ( 1.0 - alpha )*Eout_0_range.y + alpha*Eout_1_range.y;
  // interpolate the unit-base map
  mid_ubase_map.interpolate( alpha, Ein0_data.ubase_map, Ein1_data.ubase_map );

  // interpolate data at the lower and upper unit-base outgoing energies
  if(  Ein_interp.flag == Terp::LINLOG )
  {
    alpha = log( Ein / left_Ein )/log( right_Ein / left_Ein );
  }
  is_OK = mid_lower_Eout.unitbase_interp( Ein, alpha, Ein0_data.prev_data,
     Ein1_data.prev_data );
  if( !is_OK ) return false;
  is_OK = mid_upper_Eout.unitbase_interp( Ein, alpha, Ein0_data.next_data,
     Ein1_data.next_data );
  if( !is_OK ) return false;
  
  double Eout_min_ubase;
  double Eout_max_ubase;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    // save the unit-base outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_ubase = mid_ubase_map.to_unit_base( Eout_min, &is_OK );
      if( !is_OK ) return false;
      use_prev_Eout.set_E_out( Eout_min_ubase );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
    }
    use_prev_Eout.copy_coef( mid_lower_Eout );
    if( use_Eout_max )
    {
      Eout_max_ubase = mid_ubase_map.to_unit_base( Eout_max, &is_OK );
      if( !is_OK ) return false;
      use_next_Eout.set_E_out( Eout_max_ubase );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
    }
    use_next_Eout.copy_coef( mid_lower_Eout );
  }
  else
  {
    // interpolate data to the unit-base outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_ubase = mid_ubase_map.to_unit_base( Eout_min, &is_OK );
      if( !is_OK ) return false;
      use_prev_Eout.linlin_interp( Eout_min_ubase, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
      use_prev_Eout.copy_coef( mid_lower_Eout );
    }
    if( use_Eout_max )
    {
      Eout_max_ubase = mid_ubase_map.to_unit_base( Eout_max, &is_OK );
      if( !is_OK ) return false;
      use_next_Eout.linlin_interp( Eout_max_ubase, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
      use_next_Eout.copy_coef( mid_upper_Eout );
    }
  }

  return is_OK;
}
// ------------------ StdLg::standard_Legendre_param::direct_linlin_interpolate -----------
// Does direct lin-lin interpolation between two incident energies
void StdLg::standard_Legendre_param::direct_linlin_interpolate( double Ein )
{
  double dEin = right_Ein - left_Ein;
  if( dEin <= 0.0 )
  {
    Msg::FatalError( "StdLg::standard_Legendre_param::linlin_interpolate",
		"Incident energies out of order" );
  }
  double alpha = ( Ein - left_Ein )/dEin;
  // interpolate data at the lower and upper unit-base outgoing energies
  mid_lower_Eout.Ein_linlin_interp( Ein, left_Ein, Ein0_data.prev_data, right_Ein,
     Ein1_data.prev_data );
  mid_upper_Eout.Ein_linlin_interp( Ein, left_Ein, Ein0_data.next_data, right_Ein,
     Ein1_data.next_data );

  double Eout_min_linlin;
  double Eout_max_linlin;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    // save the outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_linlin = Eout_min;
      use_prev_Eout.set_E_out( Eout_min_linlin );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
    }
    use_prev_Eout.copy_coef( mid_lower_Eout );
    if( use_Eout_max )
    {
      Eout_max_linlin = Eout_max;
      use_next_Eout.set_E_out( Eout_max_linlin );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
    }
    use_next_Eout.copy_coef( mid_lower_Eout );
  }
  else
  {
    // interpolate data to the outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_linlin = Eout_min;
      use_prev_Eout.linlin_interp( Eout_min_linlin, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
      use_prev_Eout.copy_coef( mid_lower_Eout );
    }
    if( use_Eout_max )
    {
      Eout_max_linlin = Eout_max;
      use_next_Eout.linlin_interp( Eout_max_linlin, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
      use_next_Eout.copy_coef( mid_upper_Eout );
    }
  }
  // interpolate the physical energy range
  Eout_range.x = ( 1.0 - alpha )*Eout_0_range.x + alpha*Eout_1_range.x;
  Eout_range.y = ( 1.0 - alpha )*Eout_0_range.y + alpha*Eout_1_range.y;
}
// ------------------ StdLg::standard_Legendre_param::flat_interpolate ----------------
// Does histogram interpolation between two incident energies
void StdLg::standard_Legendre_param::flat_interpolate( )
{
  // copy the data at the lower incident energy
  mid_lower_Eout.set_E_out( Ein0_data.prev_data.get_E_out( ) );
  mid_lower_Eout.copy_coef( Ein0_data.prev_data );
  mid_upper_Eout.set_E_out( Ein0_data.next_data.get_E_out( ) );
  mid_upper_Eout.copy_coef( Ein0_data.next_data );

  double Eout_min_flat;
  double Eout_max_flat;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    // save the outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_flat = Eout_min;
      use_prev_Eout.set_E_out( Eout_min_flat );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
    }
    use_prev_Eout.copy_coef( mid_lower_Eout );
    if( use_Eout_max )
    {
      Eout_max_flat = Eout_max;
      use_next_Eout.set_E_out( Eout_max_flat );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
    }
    use_next_Eout.copy_coef( mid_lower_Eout );
  }
  else
  {
    // interpolate data to the outgoing energies that are to be used
    if( use_Eout_min )
    {
      Eout_min_flat = Eout_min;
      use_prev_Eout.linlin_interp( Eout_min_flat, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_prev_Eout.set_E_out( mid_lower_Eout.get_E_out( ) );
      use_prev_Eout.copy_coef( mid_lower_Eout );
    }
    if( use_Eout_max )
    {
      Eout_max_flat = Eout_max;
      use_next_Eout.linlin_interp( Eout_max_flat, mid_lower_Eout, mid_upper_Eout );
    }
    else
    {
      use_next_Eout.set_E_out( mid_upper_Eout.get_E_out( ) );
      use_next_Eout.copy_coef( mid_upper_Eout );
    }
  }
  // the physical energy range
  Eout_range.x = Eout_0_range.x;
  Eout_range.y = Eout_0_range.y;
}

// ************* class StdLg::standard_Legendre *****************
// ----------- StdLg::standard_Legendre::standard_Legendre ------------------
// constructor
StdLg::standard_Legendre::standard_Legendre( )
{
  Ein_interp.flag = Terp::NOTSET;  // Terp::LINLIN or Terp::HISTOGRAM
  Eout_interp = Terp::NOTSET;   // Terp::LINLIN or Terp::HISTOGRAM
}
// ----------- StdLg::standard_Legendre::~standard_Legendre ------------------
// destructor
StdLg::standard_Legendre::~standard_Legendre( )
{
}
// ----------- StdLg::standard_Legendre::read_data ------------------
// Reads the Legendre data
void StdLg::standard_Legendre::read_data( Dpar::data_parser& infile, int num_Ein )
{
  order = Global.Value( "outputLegendreOrder" );
  StdLg::standard_Legendre::iterator new_Ein_ptr;
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new StdLg::standard_Legendre_vector
    new_Ein_ptr = insert( end( ), StdLg::standard_Legendre_vector( ) );
    new_Ein_ptr->set_E_in( infile.get_next_double( ) );  // energy of incident particle
    new_Ein_ptr->Ein_interp = Ein_interp;
    new_Ein_ptr->Eout_interp = Eout_interp;
    int num_Eout = infile.get_next_int( );  // how many outgoing energies
    new_Ein_ptr->read_coef( infile, num_Eout, order );
  }
}
// -----------  StdLg::standard_Legendre::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool StdLg::standard_Legendre::get_Ein_range( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
    const Ddvec::dd_vector& weight_,
    const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  StdLg::standard_Legendre_param initial_param;
  bool done = initial_param.get_Ein_range( sigma_, mult_, weight_, e_flux_,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  StdLg::standard_Legendre::const_iterator this_ptr = begin( );
  double E_data = this_ptr->get_E_in( );
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  this_ptr = end( );
  --this_ptr;
  E_data = this_ptr->get_E_in( );
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ----------- StdLg::standard_Legendre::get_T ------------------
// Computes the transfer matrix
void StdLg::standard_Legendre::get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple, 
  const Ddvec::dd_vector& weight, Trf::T_matrix& transfer )
{
  bool interp_OK = ( ( Ein_interp.qualifier == Terp::UNITBASE ) &&
		     ( ( Ein_interp.flag == Terp::LINLIN ) ||
		       ( Ein_interp.flag == Terp::LINLOG ) ) ) ||
    ( ( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) &&
      ( ( Ein_interp.flag == Terp::LINLIN ) ||
        ( Ein_interp.flag == Terp::HISTOGRAM ) ) ) ||
    ( ( Ein_interp.qualifier == Terp::DIRECT ) &&
      ( ( Ein_interp.flag == Terp::LINLIN ) ||
        ( Ein_interp.flag == Terp::HISTOGRAM ) ) );

  if( !interp_OK )
  {
    Msg::FatalError( "StdLg::standard_Legendre::get_T",
		     "Incident interpolation type not implemented" );
  }
  interp_OK = ( Eout_interp == Terp::LINLIN ) || ( Eout_interp == Terp::HISTOGRAM );
  if( !interp_OK )
  {
    Msg::FatalError( "StdLg::standard_Legendre::get_T",
		     "Outgoing interpolation type not implemented" );
  }
  if( Ein_interp.qualifier == Terp::DIRECT )
  {
    Msg::Warning( "standard_Legend::get_T",
      "direct interpolation may violate energy conservation" );
  }

  // Is there data inside the incident energy range?
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }
  transfer.threshold = sigma.begin( )->x;

  long int quad_count = 0;  // number of 3-d quadratures
  long int Ein_F_count= 0;  // number of calls to standard_Legendre_F::Ein_F

  // now do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    StdLg::standard_Legendre_param Ein_param;
    Ein_param.Ein_interp = Ein_interp;
    Ein_param.Eout_interp = Eout_interp;
    Ein_param.initialize( transfer.order );

    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( &Ein_param );
    // work on this bin
    for( ; ; )
    {
      // get the incident energy interval common to all data
      set_Ein_range( Ein_bin, Ein_param );
      E_data_ladder( transfer, &Ein_param );  // loop over the outgoing energies
      // go to the next interval
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  std::cout << "2d quadratures: " << quad_count << std::endl;
  std::cout << "standard_Legendre_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "average standard_Legendre_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
}
// ----------- StdLg::standard_Legendre::setup_data --------------
// Initializes the quadrature parameters
void StdLg::standard_Legendre::setup_data( StdLg::standard_Legendre_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );

  Ein_param->Ein0_data.Ein_interp = Ein_interp;
  Ein_param->Ein0_data.Eout_interp = Eout_interp;
  Ein_param->Ein1_data.Ein_interp = Ein_interp;
  Ein_param->Ein1_data.Eout_interp = Eout_interp;

  Ein_param->this_Ein = begin( );
  Ein_param->next_Ein = Ein_param->this_Ein;
  ++Ein_param->next_Ein;

  while( ( Ein_param->next_Ein->get_E_in( ) < E_first * ( 1.0 + skip_tol ) ) ||
	 ( Ein_param->next_Ein->get_E_in( ) < (*Ein_param->Ein_ptr) *
           ( 1.0 + skip_tol ) ) )
  {
    Ein_param->this_Ein = Ein_param->next_Ein;
    ++Ein_param->next_Ein;
  }

  double first_E = Ein_param->this_Ein->get_E_in( );
  if( first_E > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_E;
    bool data_bad = Ein_param->update_pointers( first_E );
    if( data_bad )
    {
      Msg::FatalError( "StdLg::standard_Legendre::setup_data",
		       "energies inconsistent" );
    }
  }
}
// -----------  StdLg::standard_Legendre::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void StdLg::standard_Legendre::set_Ein_range( int Ein_bin, StdLg::standard_Legendre_param &Ein_param )
{
  Ein_param.set_Ein_range( );
  double this_E = Ein_param.this_Ein->get_E_in( );
  if( this_E > Ein_param.data_E_0 ) Ein_param.data_E_0 = this_E;
  this_E = Ein_param.next_Ein->get_E_in( );
  if( this_E < Ein_param.data_E_1 ) Ein_param.data_E_1 = this_E;

  if( Ein_param.data_E_1 < Ein_param.data_E_0 )
  {
    Msg::FatalError( "StdLg::standard_Legendre::set_Ein_range",
		     "check the incident energies" );
  }
  Ein_param.set_sigma_range( );
}
// ----------- StdLg::standard_Legendre::E_data_ladder --------------
// Loops through the energy data
 void StdLg::standard_Legendre::E_data_ladder( Trf::T_matrix& transfer,
    StdLg::standard_Legendre_param *Ein_param )
{
  Ein_param->reset_start( );
  // loop through the energy data
  for( ; ; )
  {
    // the physical (E, E') values for the data
    Ein_param->lower_hits.E_Eout.first.x = Ein_param->this_Ein->get_E_in( );
    Ein_param->lower_hits.E_Eout.first.y = Ein_param->Eout_0_range.x;  // physical lower E'
    Ein_param->lower_hits.E_Eout.second.x = Ein_param->next_Ein->get_E_in( );
    Ein_param->lower_hits.E_Eout.second.y = Ein_param->Eout_1_range.x;  // physical lower E'
    Ein_param->upper_hits.E_Eout.first.x = Ein_param->this_Ein->get_E_in( );     // Ein
    Ein_param->upper_hits.E_Eout.first.y = Ein_param->Eout_0_range.y;  // physical upper E'
    Ein_param->upper_hits.E_Eout.second.x = Ein_param->next_Ein->get_E_in( );    // Ein
    Ein_param->upper_hits.E_Eout.second.y = Ein_param->Eout_1_range.y;  // physical upper E'

    Eout_ladder( transfer, Ein_param );
    bool done = Ein_param->get_next_Eout( );
    if( done ) break;
  }
}
// ----------- StdLg::standard_Legendre::Eout_ladder --------------
// Loops through the outgoing energy bins
void StdLg::standard_Legendre::Eout_ladder( Trf::T_matrix& transfer,
   StdLg::standard_Legendre_param *Ein_param )
{
  double dummy = 0.0;
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    // how does the lowest interpolation line meet this E-E' box?
    Ein_param->lower_hits.hit_box( dummy, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
        ( Ein_param->lower_hits.is_above( ) ) )
    {
      // go on to the next E-E' box
      continue;
    }
    // how does the next interpolation line meet this E-E' box?
    Ein_param->upper_hits.hit_box( dummy, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( ( Eout_count > 0 ) && ( Ein_param->upper_hits.is_below( ) ) )
    {
      // we are done with this pair of E values
      break;
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
  }
}
// ----------- StdLg::standard_Legendre::one_Ebox --------------
// Integrate over one E-E' box
void StdLg::standard_Legendre::one_Ebox( Trf::T_matrix& transfer, int Eout_count,
   StdLg::standard_Legendre_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the hyperbolas eta = const meet the box
  Box::energy_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Box::energy_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Box::energy_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Box::energy_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( low_hit_ptr->hit_edge == Box::ABOVE ) ||
        ( low_hit_ptr->hit_edge == Box::TOP_OUT ) )
    {
      // do nothing---we are above the E-E' box
      continue;
    }
    else if( ( low_hit_ptr->hit_edge == Box::TOP_IN ) ||
             ( low_hit_ptr->hit_edge == Box::BOTTOM_IN ) ||
             ( low_hit_ptr->hit_edge == Box::INSIDE ) )
    {
      // the lower eta = const hyperbola is inside the E-E' box
      Ein_param->use_Eout_min = false;
      // where is the upper hyperbola?
      if( ( high_hit_ptr->hit_edge == Box::ABOVE ) ||
          ( high_hit_ptr->hit_edge == Box::TOP_OUT ) )
      {
        // integrate up to the top of the E-E' bin
        Ein_param->use_Eout_max = true;
      }
      else
      {
        // integrate up to the next eta = const hyperbola
        Ein_param->use_Eout_max = false;
      }
    }
    else
    {
      // the lower eta = const hyperbola is below the E-E' box;
      // integrate from Eout_min
      Ein_param->use_Eout_min = true;
      // where is the upper eta = const hyperbola?
      if( ( high_hit_ptr->hit_edge == Box::BOTTOM_OUT ) ||
          ( high_hit_ptr->hit_edge == Box::BELOW ) )
      {
        // do nothing---we are below the E-E' box
        continue;
      }
      else if( ( high_hit_ptr->hit_edge == Box::TOP_IN ) ||
               ( high_hit_ptr->hit_edge == Box::BOTTOM_IN ) ||
               ( high_hit_ptr->hit_edge == Box::INSIDE ) )
      {
        // the upper eta = const hyperbola is inside the E-E' box
        Ein_param->use_Eout_max = false;
      }
      else
      {
        // the upper eta = const hyperbola is above the E-E' box
        Ein_param->use_Eout_max = true;
      }
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- StdLg::standard_Legendre::next_ladder --------------
bool StdLg::standard_Legendre::next_ladder( double E_in, StdLg::standard_Legendre_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  if( !done )
  {
    static double etol = Global.Value( "tight_tol" );
    double E_tol = E_in * etol;
    //    double E_tol = 0.0;
    if( E_in + E_tol >= Ein_param->next_Ein->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->next_Ein->get_E_in( ) )
      {
        // get the next energy data
        Ein_param->this_Ein = Ein_param->next_Ein;
        ++Ein_param->next_Ein;
        if( Ein_param->next_Ein == end ( ) )
        {
          return true;
        }
      }
    }
  }

  return done;
}
// ----------- StdLg::standard_Legendre::update_T --------------
// Increments the transfer matrix
void StdLg::standard_Legendre::update_T( Trf::T_matrix &transfer, int Eout_count,
   StdLg::standard_Legendre_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // differences of nearly-equal numbers can cause problems; when to skip an interval
  static double skip_tol = Global.Value( "tight_tol" );
  
  // a vector to store the integrals, one Legendre order
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );
  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  double Ein_0 = Ein_param->Ein_0;
  double Ein_1 = Ein_param->Ein_1;
  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x <  Ein_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    Ein_param->Ein_0 = ( Ein_param->this_sigma->x < Ein_0 ) ? Ein_0 :
      Ein_param->this_sigma->x;
    Ein_param->Ein_1 = ( Ein_param->next_sigma->x > Ein_1 ) ? Ein_1 :
      Ein_param->next_sigma->x;
    // differences of nearly-equal numbers can cause problems
    if( ( Ein_param->Ein_1 - Ein_param->Ein_0 <= 
	  Ein_param->Ein_1 * skip_tol ) ||
	( Ein_param->Ein0_data.next_data.get_E_out( ) - Ein_param->Ein0_data.prev_data.get_E_out( ) <= 
	  Ein_param->Ein0_data.next_data.get_E_out( ) * skip_tol ) )
    {
      continue;  // skip this interval
    }
    else
    {
      quad_F::integrate( standard_Legendre_F::Ein_F, transfer.Ein_quad_rule,
                         Ein_param->Ein_0,
			 Ein_param->Ein_1, params, tol, &value );
      // add this integral
      transfer( Ein_param->Ein_count, Eout_count ) += value;
      // increment the function counts
      Ein_param->Ein_F_count += Ein_param->func_count;
      ++Ein_param->quad_count;
    }
  }
}
// ----------- StdLg::standard_Legendre::print --------------
void StdLg::standard_Legendre::print( )
{
  for( StdLg::standard_Legendre::iterator energy_ptr = begin( );
       energy_ptr != end( ); ++energy_ptr )
  {
    energy_ptr->print( );
  }
}

// *************** standard_Legendre_F::Ein_F **********************************
// Gets the energy probability density at given incident and outgpoing energies
bool standard_Legendre_F::Ein_F( double E_in, Qparam::QuadParamBase *Ein_param,
   Coef::coef_vector *value )
{
  bool is_OK = true;
  
  // the parameters are really StdLg::standard_Legendre_param *
  StdLg::standard_Legendre_param *e_params =
    static_cast<StdLg::standard_Legendre_param *>( Ein_param );
  e_params->func_count += 1;
  double scale = 1;  // jacobian for the unit-base map
  if( ( e_params->Ein_interp.qualifier == Terp::UNITBASE ) ||
      ( e_params->Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) )
  {
    is_OK = e_params->unitbase_interpolate( E_in );
    if( !is_OK ) return false;
    scale = e_params->mid_ubase_map.Eout_max - e_params->mid_ubase_map.Eout_min;
  }
  else if( e_params->Ein_interp.flag == Terp::HISTOGRAM )
  {
    e_params->flat_interpolate( );
  }
  else  // direct lin-lin
  {
    e_params->direct_linlin_interpolate( E_in );
  }
  int use_order = e_params->use_prev_Eout.order;
  double Eout_high = e_params->use_next_Eout.get_E_out( );
  double Eout_low = e_params->use_prev_Eout.get_E_out( );
  double d_Eout = Eout_high - Eout_low;  // unit-base outgoing energy difference

  // omit really small intervals
  static double etol = Global.Value( "tight_tol" );
  if( std::abs( d_Eout ) <= etol*Eout_high )
  {
    value->set_zero( );  // return zero
    return false;
  }
  // test for bad input
  if( ( d_Eout < 0.0 ) || ( scale <= 0.0 ) )
  {
    Msg::DebugInfo( "standard_Legendre_F::Ein_F", "energies out of order" );
    if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) )
    {
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        value->weight_1[ L_count ] = 0.0;
      }
    }
    if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
    {
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        value->weight_E[ L_count ] = 0.0;
      }
    }
    return false;
  }
  // input OK
  if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) )
  {
    if( e_params->Eout_interp == Terp::HISTOGRAM )
    {
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        value->weight_1[ L_count ] = d_Eout*e_params->mid_lower_Eout.data[ L_count ];
      }
    }
    else  // Terp::LINLIN
    {
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        value->weight_1[ L_count ] = 0.5*d_Eout*
	  ( e_params->use_prev_Eout.data[ L_count ] +
	    e_params->use_next_Eout.data[ L_count ] );
      }
    }
  }
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    double av_Eout = 0.5*( Eout_high + Eout_low );  // (unit-base)
    if( ( e_params->Ein_interp.qualifier == Terp::UNITBASE ) ||
	( e_params->Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) )
    {
      av_Eout = e_params->mid_ubase_map.un_unit_base( av_Eout );  // physical
      if( !is_OK )
      {
	Msg::DebugInfo( "standard_Legendre_F::Ein_F", "bad interpolation" );
	return false;
      }
    }
    if( e_params->Eout_interp == Terp::HISTOGRAM )
    {
      double integral = d_Eout*av_Eout;
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        value->weight_E[ L_count ] = integral*e_params->mid_lower_Eout.data[ L_count ];
      }
    }
    else  // Terp::LINLIN
    {
      for( int L_count = 0; L_count <= use_order; ++L_count )
      {
        double av_Prob = 0.5*( e_params->use_prev_Eout.data[ L_count ] +
            e_params->use_next_Eout.data[ L_count ] );
        double Prob_slope = ( e_params->use_next_Eout.data[ L_count ] -
          e_params->use_prev_Eout.data[ L_count ] )/d_Eout;
        value->weight_E[ L_count ] = d_Eout*( av_Prob*av_Eout +
          scale*Prob_slope*d_Eout*d_Eout/12.0 );
      }
    }
  }
  // weight it by flux * cross section * multiplicity * model weight
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;

  return true;
}
