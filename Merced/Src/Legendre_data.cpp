/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Legendre_data.cpp 1 2006-02-02 03:06:56Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// routines for the list of (incident energy, flux)

#include <cmath>

#include "Legendre_data.hpp"
#include "adapt_quad.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"


// *************** class Lgdata::Legendre_coefs *********************
// --------------- Lgdata::Legendre_coefs::initialize ----------------
// Allocates space
void Lgdata::Legendre_coefs::initialize( int order )
{
  LgBase::Legendre_base::initialize( order );
}
// --------------- Lgdata::Legendre_coefs::zero_data ----------------
// Sets all coefficients to zero
void Lgdata::Legendre_coefs::zero_data( )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] = 0.0;
  }
}
// --------------- Lgdata::Legendre_coefs::copy_coef ----------------
// Copies the Legendre coefficients
void Lgdata::Legendre_coefs::copy_coef( const Lgdata::Legendre_coefs& to_copy )
{
  // reset the order if necessary
  if( order != to_copy.order )
  {
    if( order >= 0 )
    {
      clean_data( );
    }
    order = to_copy.order;
    data = new double[ order + 1 ];
  }
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] = to_copy.data[ L_count ];
  }
}
// -------------- Lgdata::Legendre_coefs::only_copy_coef ------------
// Copies the Legendre coefficients; do not change the order
void Lgdata::Legendre_coefs::only_copy_coef( const Lgdata::Legendre_coefs& to_copy )
{
  // make sure that space is allocated
  if( order < 0 )
  {
    Msg::FatalError( "Lgdata::Legendre_coefs::only_copy_coef",
		     "no space allocated" );
  }
  int use_order = ( order <= to_copy.order ) ? order : to_copy.order;
  for( int L_count = 0; L_count <= use_order; ++L_count )
  {
    data[ L_count ] = to_copy.data[ L_count ];
  }
}
// -------------- Lgdata::Legendre_coefs::set_max_order ------------
// Sets the order for interpolated data
void Lgdata::Legendre_coefs::set_max_order( int left_order, int right_order )
{
  int max_order = ( left_order >= right_order ) ? left_order : right_order;
  if( order != max_order )
  {
    if( order >= 0 )
    {
      clean_data( );
    }
    order = max_order;
    data = new double[ max_order + 1 ];
  }
}
// ----------- Lgdata::Legendre_coefs::basic_linlin_interp ----------
// Interpolates the flux with weight alpha
void Lgdata::Legendre_coefs::basic_linlin_interp( double alpha,
  const Lgdata::Legendre_coefs& prev_flux, const Lgdata::Legendre_coefs& next_flux )
{
  int max_order;
  int min_order;
  bool prev_longer;
  if( prev_flux.order > next_flux.order )
  {
    max_order = prev_flux.order;
    min_order = next_flux.order;
    prev_longer = true;
  }
  else
  {
    min_order = prev_flux.order;
    max_order = next_flux.order;
    prev_longer = false;
  }
  if( order < max_order )
  {
    max_order = order;
  }
  if( order < min_order )
  {
    min_order = order;
  }
  int L_count;
  for( L_count = 0; L_count <= min_order; ++L_count )
  {
    (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ] +
      alpha * next_flux.data[ L_count ];
  }
  if( prev_longer )
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ];
    }
  }
  else
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = alpha * next_flux.data[ L_count ];
    }
  }
  for( L_count = max_order+1; L_count <= order; ++L_count )
  {
    (*this)[ L_count ] = 0.0;
  }
}
// -------------- Lgdata::Legendre_coefs::linlin_interp ------------
// Interpolates the flux at energy E_in
bool Lgdata::Legendre_coefs::linlin_interp( double E_in, const Lgdata::Legendre_coefs& prev_flux,
  const Lgdata::Legendre_coefs& next_flux )
{
  if( order < 0 )
  {
    Msg::FatalError( "Lgdata::Legendre_coefs::linlin_interp",
		     "no space for data allocated" );
  }
  set_E_in( E_in );
  // linearly interpolate
  double denom = next_flux.get_E_in( ) - prev_flux.get_E_in( );
  if( denom <= 0.0 )
  {
    Msg::DebugInfo( "Lgdata::Legendre_coefs::linlin_interp",
		    "incident energies out of order" );
    return false;
  }
  double alpha = ( E_in - prev_flux.get_E_in( ) ) / denom;
  basic_linlin_interp( alpha, prev_flux, next_flux );

  return true;
}
// ----------- Lgdata::Legendre_coefs::Ein_linlin_interp ------------
// Interpolates Legendre-coefficient data at energy E_in
bool Lgdata::Legendre_coefs::Ein_linlin_interp( double E_in, double left_Ein, 
    const Lgdata::Legendre_coefs& left_flux, double right_Ein,
    const Lgdata::Legendre_coefs& right_flux )
{
  set_E_in( left_flux.get_E_in( ) );
  // linearly interpolate
  double denom = right_Ein - left_Ein;
  if( denom <= 0.0 )
  {
    Msg::DebugInfo( "Lgdata::Legendre_coefs::Ein_linlin_interp",
		    "incident energies out of order" );
    return false;
  }
  double alpha = ( E_in - left_Ein ) / denom;
  basic_linlin_interp( alpha, left_flux, right_flux );

  return true;
}
// ----------- Lgdata::Legendre_coefs::unitbase_interp ------------
// Interpolates unit-base Legendre-coefficient data
bool Lgdata::Legendre_coefs::unitbase_interp( double E_in, double alpha, 
    const Lgdata::Legendre_coefs& left_flux,
    const Lgdata::Legendre_coefs& right_flux )
{
  static double E_tol = Global.Value( "looser_tol" );
  if( ( alpha < 0.0 ) || ( alpha > 1.0 + E_tol ) )
  {
    Msg::DebugInfo( "Lgdata::Legendre_coefs::unitbase_interp",
		    "energies out of order" );
    return false;
  }
  set_E_in( left_flux.get_E_in( ) );
  basic_linlin_interp( alpha, left_flux, right_flux );

  return true;
}
// ------------- Lgdata::Legendre_coefs::linlog_interp -------------
// Interpolates the flux linearly with respect to the logarithm of the energy
bool Lgdata::Legendre_coefs::linlog_interp( double E_in, const Lgdata::Legendre_coefs& prev_flux,
  const Lgdata::Legendre_coefs& next_flux )
{
  if( order < 0 )
  {
    Msg::FatalError( "Lgdata::Legendre_coefs::linlog_interp",
		     "no space for data allocated" );
  }
  set_E_in( E_in );
  // interpolate
  double denom = log( next_flux.get_E_in( ) / prev_flux.get_E_in( ) );
  if( denom <= 0.0 )
  {
    Msg::DebugInfo( "Lgdata::Legendre_coefs::linlog_interp",
		    "incident energies out of order" );
    return false;
  }
  int max_order;
  int min_order;
  bool prev_longer;
  if( prev_flux.order > next_flux.order )
  {
    max_order = prev_flux.order;
    min_order = next_flux.order;
    prev_longer = true;
  }
  else
  {
    min_order = prev_flux.order;
    max_order = next_flux.order;
    prev_longer = false;
  }
  if( order < max_order )
  {
    max_order = order;
  }
  if( order < min_order )
  {
    min_order = order;
  }
  int L_count;
  double alpha = log( E_in / prev_flux.get_E_in( ) ) / denom;
  for( L_count = 0; L_count <= min_order; ++L_count )
  {
    (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ] +
      alpha * next_flux.data[ L_count ];
  }
  if( prev_longer )
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ];
    }
  }
  else
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = alpha * next_flux.data[ L_count ];
    }
  }
  for( L_count = max_order+1; L_count <= order; ++L_count )
  {
    (*this)[ L_count ] = 0.0;
  }

  return true;
}

// ------------------ Lgdata::Legendre_coefs::operator+= --------------
Lgdata::Legendre_coefs& Lgdata::Legendre_coefs::operator+=( const Coef::coef_vector &to_add )
{
  if( to_add.conserve != Coef::NUMBER )
  {
    Msg::FatalError( "Lgdata::Legendre_coefs::operator+=",
		     "incompatible conserve" );
  }
  for( int ell = 0; ell <= order; ++ell )
  {
    data[ ell ] += to_add.weight_1[ ell ];
  }
  return *this;
}

// ****************** class Lgdata::mu_param ***************
// ---------------- Lgdata::mu_param::value ------------------
double Lgdata::mu_param::value( double mu, bool *is_OK )
{
  double p = left_data->linlin_interp( mu, *right_data, is_OK );
  return p;
};

// ************ class Lgdata::Legendre_data_range *******************
// ----------- Lgdata::Legendre_data_range::new_Ein --------------
// Sets up a new incident energy
void Lgdata::Legendre_data_range::new_Ein( double Ein,
					   const Ddvec::unit_base_map &ubasemap,
  Terp::Interp_Type Eoutinterp )
{
  E_in = Ein;
  Eout_interp = Eoutinterp;
  ubase_map.copy( ubasemap );
}
// ----------- Lgdata::Legendre_data_range::set_data --------------
// Sets up the data for a given range of outgloing energies
void Lgdata::Legendre_data_range::set_data(  const Lgdata::Legendre_coefs &prevdata,
  const Lgdata::Legendre_coefs &nextdata, double Eout_min, double Eout_max )
{
  if( Eout_min >= Eout_max )
  {
    Msg::FatalError( "Lgdata::Legendre_data_range::set_data",
		     "improper energy range" );
  }
  static double E_tol = Global.Value( "looser_tol" );
  if( Eout_min < ( 1.0 - E_tol )*prevdata.get_E_out( ) )
  {
    Msg::FatalError( "Lgdata::Legendre_data_range::set_data",
		     "Eout_min too low" );
  }
  if( Eout_max > ( 1.0 + E_tol )*nextdata.get_E_out( ) )
  {
    Msg::FatalError( "Lgdata::Legendre_data_range::set_data",
		     "Eout_max too high" );
  }

  if( Eout_interp == Terp::HISTOGRAM )
  {
    prev_data.set_E_out( Eout_min );
    prev_data.copy_coef( prevdata );
    next_data.set_E_out( Eout_max );
  }
  else
  {
    if( Eout_min < ( 1.0 + E_tol )*prevdata.get_E_out( ) )
    {
      prev_data.set_E_out( Eout_min );
      prev_data.copy_coef( prevdata );
    }
    else
    {
      prev_data.set_max_order( prevdata.order, nextdata.order );
      prev_data.linlin_interp( Eout_min, prevdata, nextdata );
    }
    if( Eout_max > ( 1.0 - E_tol )*nextdata.get_E_out( ) )
    {
      next_data.set_E_out( Eout_max );
      next_data.copy_coef( nextdata );
    }
    else
    {
      next_data.set_max_order( prevdata.order, nextdata.order );
      next_data.linlin_interp( Eout_max, prevdata, nextdata );
    }
  }
}
// ----------- Lgdata::Legendre_data_range::ubase_interpolate --------------
// Do unit-base interpolation between incident energies
bool Lgdata::Legendre_data_range::ubase_interpolate( double E_in,
  const Lgdata::Legendre_data_range &left_data, const Lgdata::Legendre_data_range &right_data )
{
  bool is_OK = true;
  
  Eout_interp = left_data.Eout_interp;
  double left_Ein = left_data.get_E_in( );
  double right_Ein = right_data.get_E_in( );
  double denom = right_Ein - left_Ein;
  if( denom <= 0.0 )
  {
    Msg::FatalError( "Lgdata::Legendre_data_range::ubase_interpolate",
                "incident energies out of order" );
    return false;
  }
  double alpha = ( E_in - left_Ein )/denom;
  if( ( alpha < 0.0 ) || ( alpha > 1.0 ) )
  {
    Msg::DebugInfo( "Lgdata::Legendre_data_range::ubase_interpolate",
                "extrapolation" );
    return false;
  }
  prev_data.set_max_order( left_data.prev_data.order, right_data.prev_data.order );
  is_OK = prev_data.unitbase_interp( E_in, alpha, left_data.prev_data,
			     right_data.prev_data );
  if( !is_OK ) return false;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    next_data.set_E_out( left_data.next_data.get_E_out( ) );
  }
  else
  {
    next_data.set_max_order( left_data.next_data.order, right_data.next_data.order );
    is_OK = next_data.unitbase_interp( E_in, alpha, left_data.next_data,
			       right_data.next_data );
    if( !is_OK ) return false;
  }
  set_E_in( E_in );
  ubase_map.interpolate( alpha, left_data.ubase_map, right_data.ubase_map );

  return true;
}
// ----------- Lgdata::Legendre_data_range::to_unit_base --------------
//  Maps from physical variables to unit-base
bool Lgdata::Legendre_data_range::to_unit_base( )
{
  double Eout_min = prev_data.get_E_out( );
  double Eout_max = next_data.get_E_out( );
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    char Str[1024];

    snprintf( Str, sizeof( Str ), "bad energy range: (Eout_min = %.17e, Eout_max = %.17e: diff = %e)",  Eout_min, Eout_max, scale );
    Msg::DebugInfo( "Lgdata::Legendre_data_range::to_unit_base", Str );
    return false;
  }
  ubase_map.Eout_min = Eout_min;
  ubase_map.Eout_max = Eout_max;

  prev_data.set_E_out( 0.0 );
  prev_data *= scale;
  next_data.set_E_out( 1.0 );
  next_data *= scale;

  return true;
}
// ----------- Lgdata::Legendre_data_range::short_to_unit_base --------------
// Used by cumulative points interpolation for intervals of length zero
void Lgdata::Legendre_data_range::short_to_unit_base( double dA )
{
  double Eout_max = next_data.get_E_out( );
  double scale = dA;
  ubase_map.Eout_min = Eout_max;
  ubase_map.Eout_max = Eout_max;

  prev_data.set_E_out( 0.0 );
  prev_data *= scale;
  next_data.set_E_out( 1.0 );
  next_data *= scale;
}
// ----------- Lgdata::Legendre_data_range::un_unit_base --------------
// Maps from unit-base to physical variables
bool Lgdata::Legendre_data_range::un_unit_base( )
{
  double scale = ubase_map.Eout_max - ubase_map.Eout_min;
  if( scale <= 0.0 )
  {
    Msg::DebugInfo( "Lgdata::Legendre_data_range::un_unit_base",
                 "bad unit-base map" );
    return false;
  }

  bool interp_OKL;
  double phys_Eout = ubase_map.un_unit_base( prev_data.get_E_out( ), &interp_OKL );
  prev_data.set_E_out( phys_Eout );
  prev_data *= 1.0/scale;

  bool interp_OKR;
  phys_Eout = ubase_map.un_unit_base( next_data.get_E_out( ), &interp_OKR );
  next_data.set_E_out( phys_Eout );
  next_data *= 1.0/scale;

  return( interp_OKL && interp_OKR );
}
/*
// ----------- Lgdata::Legendre_data_range::Eout_interpolate --------------
// Returns the Legendre coefficients for this outgoing energy
Lgdata::Legendre_coefs Lgdata::Legendre_data_range::Eout_interpolate( double E_out )
{ 
  Lgdata::Legendre_coefs interpolated_data;
  interpolated_data.set_max_order( prev_data.order, next_data.order );
  interpolated_data.linlin_interp( E_out, prev_data, next_data );
  return interpolated_data;
}
*/

// *************** class Lgdata::Legendre_list_base *************************
// ----------- Lgdata::Legendre_list_base::to_unit_base --------------
// Maps the data to unit base
void Lgdata::Legendre_list_base::to_unit_base( )
{
  // find the energy range
  Lgdata::Legendre_list_base::iterator Eprob_ptr = begin( );
  double Eout_min = Eprob_ptr->get_E_out( );
  ubase_map.Eout_min = Eout_min;
  Eprob_ptr = end( );
  --Eprob_ptr;
  double Eout_max = Eprob_ptr->get_E_out( );
  ubase_map.Eout_max = Eout_max;
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    Msg::FatalError( "Lgdata::Legendre_list_base::to_unit_base",
                 "improper energy range" );
  }

  for( Eprob_ptr = begin( ); Eprob_ptr != end( ); ++Eprob_ptr )
  {
    Eprob_ptr->set_E_out( ( Eprob_ptr->get_E_out( ) - Eout_min ) / scale );
    *Eprob_ptr *= scale;
  }
}
// ----------- Lgdata::Legendre_list_base::get_norm --------------
// Finds the total probability
double Lgdata::Legendre_list_base::get_norm( ) const
{
  // coding based on dd_vector::get_norm
  Lgdata::Legendre_list_base::const_iterator entry = begin();
  
  double norm = 0.0;
  double last_E = entry->get_E_out( );
  double last_P = entry->value( 0 );
  double next_E;
  double next_P;
  
// loop through the vector
  for(++entry; entry != end(); ++entry)
  {
    // get the next energy and probability density
    next_E = entry->get_E_out( );
    next_P = entry->value( 0 );

    // accumulate the probability
    if( Eout_interp == Terp::LINLIN )
    {
      norm += 0.5*(next_P + last_P)*(next_E - last_E);
    }
    else // if( Eout_interp == Terp::HISTOGRAM )
    {
      norm += last_P*(next_E - last_E);
    }

    // save the values for the next interval
    last_E = next_E;
    last_P = next_P;
  }

  return norm;
}
// ----------- Lgdata::Legendre_list_base::renorm --------------
// Normalizes the total probability
bool Lgdata::Legendre_list_base::renorm( bool truncated )
{
  // coding based on dd_vector::renorm
  double norm = get_norm( );
  Lgdata::Legendre_list_base::iterator entry = begin( );

  // error check
  if( norm == 0.0 )
  {
    if( truncated )
    {
      // direct interpolation with truncation gave norm 0
      return false;
    }
    else
    {
      Msg::Warning("Lgdata::Legendre_list_base::renorm",
 	    Msg::pastenum( "zero norm in renorm routine for E_in: ", get_E_in( ) ) );
      // fudge something
      entry->data[ 0 ] = 1.0;
      norm = get_norm( );
    }
  }

  static int norm_warn = Global.Value( "norm_warn" );
  static double norm_tol = Global.Value( "norm_tol" );
  static bool norm_good = true;
  if( norm_warn && norm_good && ( std::abs( norm - 1.0 ) > norm_tol ) )
  {
    Msg::Warning("Lgdata::Legendre_list_base::renorm",
            Msg::pastenum( "bad norm in renorm routine: ", norm ) );
    norm_good = false;
  }
  // Now, renormalize
  for( ; entry != end( ); ++entry )
  {
    *entry *= 1.0/norm;
  }
  // this is not truncated data with zero norm
  return true;
}
// ----------- Lgdata::Legendre_list_base::print --------------
// For debugging
void Lgdata::Legendre_list_base::print( )
{
  std::cout << "E_in: " << E_in << std::endl;
  Lgdata::Legendre_list_base::const_iterator this_entry = begin( );
  for( ; this_entry != end( ); ++this_entry )
  {
    this_entry->print( );
  }
  std::cout << std::endl;
}

// *************** class Lgdata::Flux_List *************************
// ------------------ Lgdata::Flux_List::read_flux ----------------
// Constructs the list from the Python data
void Lgdata::Flux_List::read_flux( Dpar::data_parser &infile, int num_Ein )
{
  // Read the interpolation
  interp = interp_flag_F::read_1d_interpolation( infile );

  order = Global.Value( "outputLegendreOrder" );
  if( order < 0 )
  {
    Msg::FatalError( "Lgdata::Flux_List::convert_flux",
		     "Desired Legendre order not set" );
  }

  Lgdata::Flux_List::iterator next_ptr;  // point to the next data
  // read the flux data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // what is the next incident energy?
    double this_E = infile.get_next_double( );

    // append a link
    next_ptr = insert( end( ), Legendre_coefs( ) );
    next_ptr->initialize( order );
    next_ptr->set_E_in( this_E );

    int file_order = infile.get_next_int( ) - 1;  // Legendre order of input data
    int save_coefs = ( file_order > order ) ? order : file_order;
    int coef_count;
    // A coefficient-less flux is zero; this also makes padding well-defined.
    double next_flux = 0.0;
    for( coef_count = 0; coef_count <= save_coefs; ++coef_count )
    {
      next_flux = infile.get_next_double( );
      ( *next_ptr )[ coef_count ] = next_flux;
    }
    // we may need to fill in the higher-order Legendre coefficients
    if( file_order < order )
    {
      for( coef_count = save_coefs+1; coef_count <= order; ++coef_count )
      {
        ( *next_ptr )[ coef_count ] = next_flux;
      }
    }
    // we may need to discard high-order coefficients
    else if( file_order > order )
    {
      for( coef_count = save_coefs+1; coef_count <= file_order; ++coef_count )
      {
        next_flux = infile.get_next_double( );
      }
    }
  }
}
// ------------------ Lgdata::Flux_List::value ----------------
Lgdata::Legendre_coefs Lgdata::Flux_List::value( double E_in, Lgdata::Flux_List::const_iterator &ptr ) const
{
  Lgdata::Flux_List::const_iterator next_ptr = ptr;
  ++next_ptr;  // point to the next data
  if( next_ptr == end( ) )
  {
    Msg::FatalError( "Lgdata::Flux_List::value", Msg::pastenum( "energy ", E_in) +
		 " out of range" );
  }
  double denom = next_ptr->get_E_in( ) - ptr->get_E_in( );
  if( denom <= 0.0 )
  {
    Msg::FatalError( "Lgdata::Flux_List::value", "energies out of order" );
  }

  // linear interpolate
  double alpha = ( E_in - ptr->get_E_in( ) ) / denom;
  Lgdata::Legendre_coefs answer;
  answer.initialize( order );
  answer.set_E_in( E_in );
  for( int L_order = 0; L_order <= order; ++L_order )
  {
    answer[ L_order ] = ( 1 - alpha ) * ptr->data[ L_order ] +
      alpha * next_ptr->data[ L_order ];
  }
  return answer;
}

// *************** class Lgdata::weight_vector *************************
// ------------------ Lgdata::weight_vector::increment ----------------
// Adds the integrals over ( E_left, E_right )
void Lgdata::weight_vector::increment( Lgdata::Flux_List &e_flux, 
  Lgdata::Flux_List::const_iterator this_flux, double E_left, double E_right )
{
  Lgdata::Legendre_coefs left_value = e_flux.value( E_left, this_flux );
  Lgdata::Legendre_coefs right_value = e_flux.value( E_right, this_flux );
  // use the trapezoid rule
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] += 0.5 * ( E_right - E_left ) *
      ( left_value[ L_count ] + right_value[ L_count ] );
  }
}
// ------------------ Lgdata::weight_vector::invert ----------------
// take the reciprocals
void Lgdata::weight_vector::invert( )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    if( data[ L_count ] <= 0.0 )
    {
      Msg::FatalError( "Lgdata::weight_vector::invert",
                   "negative flux weight" );
    }
    data[ L_count ] = 1.0 / data[ L_count ];
  }
}

// **************** functions to integrate ******************
// ----------------- to_Legendre_F::mu_F ----------------------
// Function for the quadrature over mu: Legendre * probability density table
bool to_Legendre_F::mu_F( double mu_in, Qparam::QuadParamBase *void_param,
   Coef::coef_vector *value )
{
  // the parameters are really Lgdata::mu_param
  Lgdata::mu_param *params = static_cast< Lgdata::mu_param* >( void_param );
  params->func_count += 1;

  math_F::Legendre( mu_in, value );   // the Legendre polynomials
  bool is_OK = true;
  *value *= params->value( mu_in, &is_OK );

  return is_OK;
}

// ------------------ to_Legendre_F::from_table --------------
// Evaluates the Legendre moments
int to_Legendre_F::from_table( const Ddvec::dd_vector& mu_table,
			       Lgdata::Legendre_coefs *coefs )
{
  int mu_quad_count = 0;
  coefs->zero_data( );
  
  Coef::coef_vector integral( coefs->order, Coef::NUMBER );
  integral.set_zero( );
  coefs->set_E_in( mu_table.get_E_in( ) );
  Lgdata::mu_param mu_params;   // the quadrature parameters
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &mu_params );
  mu_params.left_data = mu_table.begin( );
  mu_params.right_data = mu_params.left_data;
  ++mu_params.right_data;
  for( ; mu_params.right_data != mu_table.end( );
       mu_params.left_data = mu_params.right_data, ++mu_params.right_data )
  {
    static double tol = Global.Value( "quad_tol" );
    Qmeth::Quadrature_Rule quad_rule;
    quad_rule.adaptive = false;   // We usually don't need adaptive quadrature
    integral.set_zero( );
    if( coefs->order < 2 )
    {
      quad_rule.quad_method = Qmeth::GAUSS2;
    }
    else if( coefs->order < 7 )
    {
      quad_rule.quad_method = Qmeth::GAUSS4;
    }
    else if( coefs->order < 11 )
    {
      quad_rule.quad_method = Qmeth::GAUSS6;
    }
    else if( coefs->order < 19 )
    {
      quad_rule.quad_method = Qmeth::GAUSS10;
    }
    else
    {
      // use adaptive quadrature
      quad_rule.adaptive = true;
      quad_rule.quad_method = Qmeth::GAUSS4;
    }
    quad_F::integrate( to_Legendre_F::mu_F, quad_rule, mu_params.left_data->x,
		       mu_params.right_data->x, params, tol, &integral );
    *coefs += integral;
    mu_quad_count += mu_params.func_count;
  }
  return mu_quad_count;
}
