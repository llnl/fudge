/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: reaction.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Implement the classes used for identification of the reaction data

#include <iostream>
#include <string>

#include "reaction.hpp"
#include "messaging.hpp"
#include "global_params.hpp"


//****************** class React::num_check_param *****************
// --------------- React::num_check_param::value -----------------
// Computes the product sigma*multiplicity*flux*weight
double React::num_check_param::value( double E_in )
{
  // ignore the checks on interpolation
  bool is_OK;
  
  if( E_in < model_weight.first.x )
  {
    return 0.0;
  }
  double vaerde = sigma.value( E_in, &is_OK ) *
    mult.value( E_in, &is_OK ) * e_flux.value( E_in, &is_OK ) *
    model_weight.value( E_in, &is_OK );

  return vaerde;
}

//****************** class React::reaction  *****************
// --------------- React::reaction::read_quadrature -----------------
// Interprets the quadrature rule
Qmeth::Quadrature_Method React::reaction::read_quadrature( Dpar::data_parser &input_file )
{
  std::string quadrature = input_file.get_text( );
  string_F::Tolower( quadrature );
  Qmeth::Quadrature_Method use_quad = Qmeth::GAUSS2;
  if( quadrature == "gauss1" )
  {
    use_quad = Qmeth::GAUSS1;  // midpoint rule
  }
  else if( quadrature == "gauss2" )
  {
    use_quad = Qmeth::GAUSS2;  // Gauss 2-nd order
  }
  else if( quadrature == "gauss3" )
  {
    use_quad = Qmeth::GAUSS3;  // Gauss 3-rd order
  }
  else if( quadrature == "gauss4" )
  {
    use_quad = Qmeth::GAUSS4; // Gauss 4-th order
  }
  else if( quadrature == "gauss6" )
  {
    use_quad = Qmeth::GAUSS6; // Gauss 6-th order
  }
  else if( quadrature == "square root" )
  {
    use_quad = Qmeth::WEIGHT_L1; // 1st-order Gauss with sqrt(x) singularity
  }
  else
  {
    Msg::FatalError( "React::reaction::read_quadrature",
		     "quadrature rule " +
		quadrature + " undefined." );
  }
  return use_quad;
}
// --------------- React::reaction::read_Global -----------------
// Reads Global parameters from the input file
bool React::reaction::read_Global( const std::string &dataID, Dpar::data_parser &input_file )
{
  std::string testCopy = dataID;
  string_F::Tolower( testCopy );
  Glb::ss_list::iterator SSlist_ptr;
  std::string pvalue;

  bool found_it = Global.find( testCopy, SSlist_ptr );
  if( found_it )
  {
    pvalue = input_file.get_text( );
    if( !SSlist_ptr->set_in_command )
    {
      string_F::Tolower( pvalue );
      SSlist_ptr->set_value( pvalue );
    }
  }
  else  // some special flags
  {
    if( testCopy == "inversewavelengthtoenergyfactor" )
    {
      pvalue = input_file.get_text( );
      Global.set( "x_to_energy", pvalue );
      found_it = true;
    }
    else if( testCopy == "electron mass" )
    {
      pvalue = input_file.get_text( );
      Global.set( "m_electron", pvalue );
      found_it = true;
    }
    else if( testCopy == "thompsonscattering" )
    {
      pvalue = input_file.get_text( );
      Global.set( "thompson", pvalue );
      found_it = true; 
    }
  }

  return found_it;
}
// --------------- React::reaction::process_data -----------------
// process data
void React::reaction::process_data( Dpar::data_parser &input_file,
				    std::ofstream *output_file )
{
  // the default quadrature methods
  transfer.Ein_quad_rule.quad_method = Qmeth::GAUSS2;  // for incident energy
  transfer.Eout_quad_rule.quad_method = Qmeth::GAUSS2;  // for outgoing energy
  transfer.mu_quad_rule.quad_method = Qmeth::GAUSS2;  // for outgoing cosine
  transfer.mucm2_quad_rule.quad_method = Qmeth::GAUSS2;  // for second cosine, 2-step reaction
  transfer.w_quad_rule.quad_method = Qmeth::GAUSS2;  // for w parameter, 2-step reaction

  transfer.output_file = output_file;
  std::string File_type = input_file.get_dataID( );
  // Test the file identity and the version
  if( File_type == "DONE" )
  {
    Msg::FatalError( "React::reaction::process_data", "The input file is empty" );
  }
  else if( ( File_type != "xndfgenTransferMatrix" ) &&
      ( File_type != "xendlTransferMatrix" ) )
  {
    Msg::FatalError( "React::reaction::process_data", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": File type: " + File_type + " not implemented" );
  }
  version = input_file.get_next_double(  );
  if( version != 1.0 )
  {
    Msg::FatalError( "React::reaction::process_data",
      Msg::pastenum( "Implement version ", version ) );
  }
  *output_file << "merced: version " << version << std::endl;

  // Get the type of data
  std::string processID = input_file.get_dataID( );
  string_F::Tolower( processID );
  while( processID == "comment" )
  {
    input_file.get_comment( *(transfer.output_file) );
    processID = input_file.get_dataID( );
    string_F::Tolower( processID );
  }

  std::string process_orig = input_file.get_text( );
  std::string process = process_orig;
  string_F::Tolower( process );
  if( process == "two body transfer matrix" )
  {
    two_body( input_file );
  }
  else if( process == "legendre two body transfer matrix" )
  {
    Legendre_two_body( input_file );
  }
  else if( process == "isotropic table" )
  {
    do_isotropic( input_file );
  }
  else if( process == "legendre eepp data transfer matrix" )
  {
    Legendre( input_file );
  }
  else if( process == "legendre energy-angle data" )
  {
    do_ENDFLegendre( input_file );
  }
  else if( process == "double differential emuepp data transfer matrix" )
  {
    joint_dist_table.ENDL_data = true;
    do_joint_dist( input_file );
  }
  else if( process == "endf double differential emuepp data" )
  {
    joint_dist_table.ENDL_data = false;
    do_joint_dist( input_file );
  }
  else if( process == "pointwise energy-angle data" )
  {
    do_energy_angle( input_file );
  }
  else if( process == "uncorrelated energy-angle data transfer matrix" )
  {
    do_uncorr( input_file );
  }
  else if( process == "compton scattering" )
  {
    do_Compton( input_file );
  }
  else if( process == "coherent scattering" )
  {
    do_coherent( input_file );
  }
  else if( process == "evaporation spectrum" )
  {
    do_evaporation( input_file );
  }
  else if( process == "maxwell spectrum" )
  {
    do_Maxwell( input_file );
  }
  else if( process == "watt spectrum" )
  {
    do_Watt( input_file );
  }
  else if( process == "madland-nix spectrum" )
  {
    do_MadlandNix( input_file );
  }
  else if( process == "kalbach spectrum" )
  {
    do_Kalbach( input_file );
  }
  else if( process == "phase space spectrum" )
  {
    do_phase_space( input_file );
  }
  else if( process == "general evaporation" )
  {
    do_gen_evap( input_file );
  }
  else if( process == "two step two body reaction" )
  {
    do_two_step( input_file );
  }
  else
  {
    Msg::FatalError( "React::reaction::process_data", 
		Msg::pastenum( "line ", input_file.line_count ) +
             ": process " + process_orig + " not implemented" );
  }

  if( transfer.order >= 0 )
  {
    // Print the input parameters
    Global.print();
    write_transfer( );
  }
}
// --------------- React::reaction::common_input -----------------
// Reads and processes input data common to all reactions
bool React::reaction::common_input( const std::string &dataID, Dpar::data_parser &input_file )
{
  bool found_it = false;
  // what type of data is it?
  std::string tempID = dataID;
  string_F::Tolower( tempID );
  if( tempID == "comment" )
  {
    input_file.get_comment( *(transfer.output_file) );
    found_it = true;
  }
  else if( tempID == "projectile's group boundaries" )
  {
    int num_bd = input_file.get_next_int( );
    transfer.num_Ein_bins = num_bd - 1;
    // read the bin boundaries
    transfer.in_groups.read_bd( input_file, num_bd );
    found_it = true;
  }
  else if( tempID == "product's group boundaries" )
  {
    int num_bd = input_file.get_next_int( );
    transfer.num_Eout_bins = num_bd - 1;
    // read the bin boundaries
    transfer.out_groups.read_bd( input_file, num_bd );
    //! Put all low-energy output in the bottom bin
    transfer.out_groups[ 0 ] = 0.0;
    //! Put all high-energy output in the top bin
    transfer.out_groups[ transfer.num_Eout_bins ] *= 5.0;
    found_it = true;
  }
  else if( tempID == "fluxes" )
  {
    int num_flux = input_file.get_next_int( );
    // read the approximate flux
    transfer.e_flux.read_flux( input_file, num_flux );
    if( transfer.e_flux.interp != Terp::LINLIN )
    {
      Msg::FatalError( "React::reaction::common_input",
		  "only linlear-linear flux is implemented" );
    }
    // Validate that flux spans at least the same range as projectile group boundaries
    if( transfer.in_groups.size( ) > 0 )
    {
      double flux_min = transfer.e_flux.begin( )->get_E_in( );
      double flux_max = ( --transfer.e_flux.end( ) )->get_E_in( );
      double group_min = transfer.in_groups.front( );
      double group_max = transfer.in_groups.back( );

      if( flux_min > group_min || flux_max < group_max )
      {
        Msg::FatalError( "React::reaction::common_input",
          Msg::pastenum( "Flux energy range [", flux_min ) +
          Msg::pastenum( ", ", flux_max ) +
          "] does not span the projectile group boundaries [" +
          Msg::pastenum( "", group_min ) +
          Msg::pastenum( ", ", group_max ) +
          "]. The flux must be defined over at least the same energy range as the group boundaries." );
      }
    }
    found_it = true;
  }
  else if( tempID == "cross section" )
  {
    int num_sigma = input_file.get_next_int( );
    // read the cross section
    cross_section.read_data_interp( input_file, num_sigma );
    if( cross_section.interp_type != Terp::LINLIN )
    {
      Msg::FatalError( "React::reaction::common_input", 
		  "only linlear-linear cross section is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "multiplicity" )
  {
    int num_mult = input_file.get_next_int( );
    // read the multiplicities
    multiple.read_data_interp( input_file, num_mult );
    if( multiple.interp_type != Terp::LINLIN )
    {
      Msg::FatalError( "React::reaction::common_input", 
		  "only linlear-linear multiplicity is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "weight" )
  {
    int num_weight = input_file.get_next_int( );
    // read the weights
    model_weight.read_data_interp( input_file, num_weight );
    if( ( model_weight.interp_type != Terp::HISTOGRAM ) &&
	( model_weight.interp_type != Terp::LINLIN ) )
    {
      Msg::FatalError( "React::reaction::common_input", 
		  "only histogram and linear-linear weight are implemented" );
    }
    found_it = true;
  }
  else if( tempID == "projectile's mass" )
  {
    particle_info.mProj = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "target's mass" )
  {
    particle_info.mTarg = input_file.get_next_double( );
    found_it = true;
  }
  else if( ( tempID == "product's mass" ) ||
	   ( tempID == "first product's mass" ) )
  {
    particle_info.mProd = input_file.get_next_double( );
    found_it = true;
  }
  else if( ( tempID == "residual's mass" ) ||
	   ( tempID == "first residual's mass" ) )
  {
    particle_info.mRes = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "projectile frame" )
  {
    std::string frame = input_file.get_text( );
    string_F::Tolower( frame );
    if( frame == "lab" )
    {
      frame_in = React::LAB;
    }
    else
    {
      Msg::FatalError( "React::reaction::common_input",
		       "invalid incident frame: " + frame );
    }
    found_it = true;
  }
  else if( tempID == "product frame" )
  {
    std::string frame = input_file.get_text( );
    string_F::Tolower( frame );
    if( frame == "lab" )
    {
      frame_out = React::LAB;
    }
    else if( frame == "centerofmass" )
    {
      frame_out = React::CM;
    }
    else
    {
      Msg::FatalError( "React::reaction::common_input",
		       "invalid outgoing frame: " + frame );
    }
    found_it = true;
  }
  else if( tempID == "conserve" )
  {
    std::string cons = input_file.get_text( );
    string_F::Tolower( cons );
    if( cons == "number" )
    {
      transfer.conserve = Coef::NUMBER;
    }
    else if( cons == "energy" )
    {
      transfer.conserve = Coef::ENERGY;
    }
    else if( cons == "both" )
    {
      transfer.conserve = Coef::BOTH;
    }
    else
    {
      Msg::Warning( "React::reaction::common_input",
		    "invalid conservation type: " + cons );
    }
    found_it = true;
  }
  else if( tempID == "do adaptive quadrature" )
  {
    std::string adapt = input_file.get_text( );
    string_F::Tolower( adapt );
    if( adapt == "false" )
    {
      // turn off all quadratures
      transfer.Ein_quad_rule.adaptive = false;
      transfer.Eout_quad_rule.adaptive = false;
      transfer.mu_quad_rule.adaptive = false;
      transfer.mucm2_quad_rule.adaptive = false;
      transfer.w_quad_rule.adaptive = false;
      found_it = true;
    }
  }
  else if( tempID == "quadrature method" )
  {
    // all quadratures use the same method
    transfer.Ein_quad_rule.quad_method = read_quadrature( input_file );
    transfer.Ein_quad_rule.input_set = true;
    transfer.Eout_quad_rule = transfer.Ein_quad_rule;
    transfer.mu_quad_rule = transfer.Ein_quad_rule;
    transfer.mucm2_quad_rule = transfer.Ein_quad_rule;
    transfer.w_quad_rule = transfer.Ein_quad_rule;
    found_it = true;
  }
  else if( tempID == "ein quadrature method" )
  {
    transfer.Ein_quad_rule.quad_method = read_quadrature( input_file );
    transfer.Ein_quad_rule.input_set = true;
    found_it = true;
  }
  else if( tempID == "eout quadrature method" )
  {
    transfer.Eout_quad_rule.quad_method = read_quadrature( input_file );
    transfer.Eout_quad_rule.input_set = true;
    found_it = true;
  }
  else if( tempID == "mu quadrature method" )
  {
    transfer.mu_quad_rule.quad_method = read_quadrature( input_file );
    transfer.mu_quad_rule.input_set = true;
    found_it = true;
  }
  else if( tempID == "mucm2 quadrature method" )
  {
    transfer.mucm2_quad_rule.quad_method = read_quadrature( input_file );
    transfer.mucm2_quad_rule.input_set = true;
    found_it = true;
  }
  else if( tempID == "w quadrature method" )
  {
    transfer.w_quad_rule.quad_method = read_quadrature( input_file );
    transfer.w_quad_rule.input_set = true;
    found_it = true;
  }
  else if( tempID == "interpolate eout integrals" )
  {
    Msg::Warning( "React::reaction::common_input",
	     "interpolate eout integrals option obsolete" );
    input_file.get_text( );  // read and ignore the option
    found_it = true;
  }
  else if( read_Global( tempID, input_file ) )
  {
    found_it = true;
  }
  else if( dataID == "Bad dataID" )
  {
    Msg::FatalError( "React::reaction::common_input", "Expected dataID in " +
		Msg::pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
  }

  return found_it;
}
// --------------- React::reaction::two_body -----------------
// Reads and processes input data for discrete 2-body reactions
void React::reaction::two_body( Dpar::data_parser &input_file )
{
  frame_out = React::CM;  // default
  for ( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
    dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "reaction's q value" )
    {
      angdist.Q =  input_file.get_next_double( );
    }
    else if( tempID == "angular data" )
    {
      // set a temporary threshold based on the cross section data
      // this is reset later, based on the kinetics
      angdist.threshold = cross_section.begin( )->x;
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      angdist.read_data( input_file, num_Ein );
    }
    else
    {
      Msg::FatalError( "React::reaction::two_body", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::CM )
  {
    Msg::FatalError( "React::reaction::two_body",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  angdist.particle_info.copy( particle_info );
  angdist.get_T( cross_section, model_weight, transfer );
}
// --------------- React::reaction::Legendre_two_body ----------------------
// Gets the transfer matrix from discrete 2-body Legendre polynomial angular data
void React::reaction::Legendre_two_body( Dpar::data_parser &input_file )
{
  frame_out = React::CM;  // default
  for ( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
    dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "reaction's q value" )
    {
      // not used by gammas from neutron capture
      LegendreAngle.Q =  input_file.get_next_double( );
    }
    else if( tempID == "legendre coefficients" )
    {
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      LegendreAngle.read_data( input_file, num_Ein );
    }
    else
    {
      Msg::FatalError( "React::reaction::Legendre_two_body", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::CM )
  {
    Msg::FatalError( "React::reaction::Legendre_two_body",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }

  LegendreAngle.particle_info.copy( particle_info );
  LegendreAngle.get_T( cross_section, model_weight, transfer );
}
// --------------- React::reaction::do_two_step ----------------------
// Gets the transfer matrix for a 2-step 2-body reaction
void React::reaction::do_two_step( Dpar::data_parser &input_file )
{
  frame_out = React::CM;  // default
  for ( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
    dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "first step's q value" )
    {
      twoStep.first_Q = input_file.get_next_double( );
    }
    else if( tempID == "second step's q value" )
    {
      twoStep.second_Q = input_file.get_next_double( );
    }
    else if( tempID == "second product's mass" )
    {
      twoStep.step2_particles.mProd = input_file.get_next_double( );
    }
    else if( tempID == "second residual's mass" )
    {
      twoStep.step2_particles.mRes = input_file.get_next_double( );
    }
    else if( tempID == "legendre coefficients" )
    {
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      twoStep.read_data( input_file, num_Ein );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_two_step", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::CM )
  {
    Msg::FatalError( "React::reaction::do_two_step",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  twoStep.step1_particles.copy( particle_info );
  twoStep.get_T( cross_section, model_weight, transfer );
}
// --------------- React::reaction::Legendre ----------------------
// Reads and processes Legendre expansions of double differential data
void React::reaction::Legendre( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "leeppdata" )
    {
      int num_moments = input_file.get_next_int( );
      // read the Legendre coefficients of the energy probability density
      energyMoments.read_data( input_file, num_moments );
    }
    else
    {
      Msg::FatalError( "React::reaction::Legendre", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::Legendre",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  energyMoments.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_ENDFLegendre ----------------------
// Reads and processes Legendre expansions of double differential data
void React::reaction::do_ENDFLegendre( Dpar::data_parser &input_file )
{
  Terp::two_d_interp Ein_interp;
  Terp::Interp_Type Eout_interp;
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "legendre data by incident energy" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      // read the Legendre coefficients of the energy probability density
      if( frame_out == React::LAB )
      {
        standard_legendre.Ein_interp = Ein_interp;
        standard_legendre.Eout_interp = Eout_interp;
	standard_legendre.read_data( input_file, num_Ein );
      }
      else
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_data( input_file, num_Ein );
      }
    }
    else
    {
      Msg::FatalError( "React::reaction::do_ENDFLegendre", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
                  " not implemented" );
    }
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == React::LAB )
  {
    standard_legendre.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- React::reaction::do_isotropic ----------------------
// Reads and processes tables of isotropic energy probability density
void React::reaction::do_isotropic( Dpar::data_parser &input_file )
{
  Terp::two_d_interp Ein_interp;
  Terp::Interp_Type Eout_interp;
  frame_out = React::LAB;  // default

  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "eeppdata" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      if( frame_out == React::CM )
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_isotropic( input_file, num_Ein );
      }
      else
      {
        isotrop.Ein_interp = Ein_interp;
        isotrop.Eout_interp = Eout_interp;
        isotrop.read_data( input_file, num_Ein );
      }
    }
    else
    {
      Msg::FatalError( "React::reaction::do_isotropic", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }

  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == React::LAB )
  {
    isotrop.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- React::reaction::do_uncorr ----------------------
// Reads and processes uncorrelated expansions of double differential data
void React::reaction::do_uncorr( Dpar::data_parser &input_file )
{
  Terp::two_d_interp Ein_interp;
  Terp::Interp_Type Eout_interp;
  frame_out = React::LAB;  // default
  bool is_isotropic = false;    // Is the angular data isotropic?
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "angular data" )
    {
      uncorr.mu_table = true;
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      angdist.threshold = cross_section.begin( )->x;
      angdist.read_data( input_file, num_Ein );
      is_isotropic = angdist.isotropic( );
    }
    else if( tempID == "legendre coefficients" )
    {
      uncorr.mu_table = false;
      int num_Ein = input_file.get_next_int( );
      // read the Legendre coefficients
      uncorr.read_Legendre( input_file, num_Ein );
    }
    else if( tempID == "eeppdata" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      // read the energy probability density
      if( is_isotropic && ( frame_out == React::CM ) )
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_isotropic( input_file, num_Ein );
      }
      else
      {
        uncorr.Ein_interp = Ein_interp;
        uncorr.Eout_interp = Eout_interp;
	uncorr.read_data( input_file, num_Ein, &angdist );
      }
    }
    else
    {
      Msg::FatalError( "React::reaction::do_uncorr", 
		Msg::pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }

  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == React::LAB )
  {
    uncorr.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    if( !is_isotropic )
    {
      Msg::FatalError( "React::reaction::do_uncorr",
        "center-of-mass data is implemented only for isotropic data" );
    }
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- React::reaction::do_joint_dist ----------------------
// Reads and processes double-differential tabular data
void React::reaction::do_joint_dist( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "angular data" )
    {
      int num_Ein = input_file.get_next_int( );
      joint_dist_table.angle_data.threshold = cross_section.begin( )->x;
      // read the angular probability density
      joint_dist_table.angle_data.read_data( input_file, num_Ein );
    }
    else if( tempID == "emueppdata" )
    {
      int num_Eout = input_file.get_next_int( );
      // read the double-differential data tables
      joint_dist_table.read_data( input_file, num_Eout );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_joint_dist", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_joint_dist",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  joint_dist_table.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_energy_angle -----------------
// Reads and processes double-differential tabular cm data
void React::reaction::do_energy_angle( Dpar::data_parser &input_file )
{
  frame_out = React::CM;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "eepmupdata" )
    {
      int num_Eout = input_file.get_next_int( );
      // read the double-differential data tables
      if( frame_out == React::CM )
      {
        cm_joint_table.read_data( input_file, num_Eout );
      }
      else
      {
	lab_joint_table.read_data( input_file, num_Eout );
      }
    }
    else
    {
      Msg::FatalError( "React::reaction::do_cm_joint", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }

  if( frame_out == React::CM )
  {
    cm_joint_table.particles.copy( particle_info );
    cm_joint_table.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    lab_joint_table.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- React::reaction::do_Compton ----------------------
// Reads and processes scattering factor data for Compton scattering
void React::reaction::do_Compton( Dpar::data_parser &input_file )
{
  // set defaults
  frame_out = React::LAB;

  Global.set( "scale_rows", 0 );   // the check is invalid---sigma not lin-lin

  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "scatteringfactordata" )
    {
      int num_factor = input_file.get_next_int( );
      // read the scattering factor data
      compton.file_data.read_data_interp( input_file, num_factor );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_Compton", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_Compton",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  compton.get_T( transfer, cross_section );
  write_xsec( );
}
// --------------- React::reaction::do_coherent ----------------------
// Reads and processes scattering factor data for coherent scattering
void React::reaction::do_coherent( Dpar::data_parser &input_file )
{
  // set defaults
  frame_out = React::LAB;

  Global.set( "scale_rows", 0 );   // the check is invalid---sigma not lin-lin

  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "formfactordata" )
    {
      int num_factor = input_file.get_next_int( );
      // read the form factors
      coherent_model.file_data.read_data_interp( input_file, num_factor );
    }
    else if( tempID == "realanomalousfactor" )
    {
      int num_factor = input_file.get_next_int( );
      // read the real anomalous factor
      coherent_model.realAnomalous.read_data_interp( input_file, num_factor );
    }
    else if( tempID == "imaginaryanomalousfactor" )
    {
      int num_factor = input_file.get_next_int( );
      // read the imaginary anomalous factor
      coherent_model.imaginaryAnomalous.read_data_interp( input_file, num_factor );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_coherent", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_coherent",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  coherent_model.get_T( transfer, cross_section );
  write_xsec( );
}
// --------------- React::reaction::do_evaporation ----------------------
// Reads and processes data for evaporation spectra
void React::reaction::do_evaporation( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      evaporation_model.U = input_file.get_next_double( );
    }
    else if( tempID == "theta" )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta values
      evaporation_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_evaporation", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_evaporation",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  transfer.allocate( );
  evaporation_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_Maxwell ----------------------
// Reads and processes data for Maxwell spectra
void React::reaction::do_Maxwell( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      Maxwell_model.U = input_file.get_next_double( );
    }
    else if( tempID == "theta" )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta values
      Maxwell_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_Maxwell", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_Maxwell",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Maxwell_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_Watt ----------------------
// Reads and processes data for Watt spectra
void React::reaction::do_Watt( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      Watt_model.U = input_file.get_next_double( );
    }
    else if( tempID == "b" )
    {
      int num_b = input_file.get_next_int( );
      // read the b values
      Watt_model.b_data.read_data_interp( input_file, num_b );
    }
    else if( ( tempID == "theta" ) || ( tempID == "a" ) )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta, i.e., a values
      Watt_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_Watt", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_Watt",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Watt_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_MadlandNix ----------------------
// Reads and processes data for MadlandNix spectra
void React::reaction::do_MadlandNix( Dpar::data_parser &input_file )
{
  double Efl = 0.0;  // kinetic energy of the lighter fragment
  double Efh = 0.0;  // kinetic energy of the heavier fragment
  frame_out = React::LAB;  // default
  // The Madland-Nix model does not specify a maximum outgoing energy
  MadlandNix_model.maxEout = 0.0;  // not set, default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "maxeout" )
    {
      MadlandNix_model.maxEout = input_file.get_next_double( );
    }
    else if( tempID == "efl" )
    {
      Efl = input_file.get_next_double( );
      MadlandNix_model.Efl = Efl;
    }
    else if( tempID == "efh" )
    {
      Efh = input_file.get_next_double( );
      MadlandNix_model.Efh = Efh;
    }
    else if( tempID == "u" )
    {
      MadlandNix_model.U = input_file.get_next_double( );
      Msg::Warning( "React::reaction::do_MadlandNix",
        "U is not used in the Madland-Nix model." );
    }
    else if( tempID == "tm" )
    {
      int num_TM = input_file.get_next_int( );
      // read the TM values
      MadlandNix_model.read_data_interp( input_file, num_TM );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_MadlandNix", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( ( Efl <= 0.0 ) || ( Efh <= 0.0 ) )
  {
    Msg::FatalError( "React::reaction::do_MadlandNix",
      "The kinetic energies of the fission fragments should be positive." );
  }
  if( transfer.conserve != Coef::NUMBER )
  {
    Msg::Warning( "React::reaction::do_MadlandNix", 
	     "Energy-conserving transfer matrices not implemented." );
    transfer.conserve = Coef::NUMBER;
  }

  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_MadlandNix",
		"only emission in the laboratory frame is supported" );
  }
  // set up the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  MadlandNix_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_Kalbach ----------------------
// Reads and processes data for Kalbach spectra
void React::reaction::do_Kalbach( Dpar::data_parser &input_file )
{
  frame_out = React::CM;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "projectile's za" )
    {
      Kalbach_model.kalbach_a.projectile.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "target's za" )
    {
      Kalbach_model.kalbach_a.target.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "compound's mass" )
    {
      Kalbach_model.kalbach_a.compound.mass = input_file.get_next_double( );
    }
    else if( tempID == "product's za" )
    {
      Kalbach_model.kalbach_a.eject.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "kalbach probabilities" )
    {
      int num_Ein = input_file.get_next_int( );
      Kalbach_model.read_probability( input_file, num_Ein );
    }
    else if( tempID == "kalbach r parameter" )
    {
      int num_Ein = input_file.get_next_int( );
      Kalbach_model.read_r( input_file, num_Ein );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_Kalbach", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::CM )
  {
    Msg::FatalError( "React::reaction::do_Kalbach",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  Kalbach_model.kalbach_a.copy_masses( particle_info );
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Kalbach_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_phase_space ----------------------
// Reads and processes data for phase_space spectra
void React::reaction::do_phase_space( Dpar::data_parser &input_file )
{
  transfer.mu_quad_rule.quad_method = Qmeth::GAUSS4;  // default mu quadrature

  frame_out = React::CM;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "total mass" )
    {
      phase_space_model.totalMass = input_file.get_next_double( );
    }
    else if( tempID == "number of particles" )
    {
      phase_space_model.numParticles = input_file.get_next_int( );
    }
    else if( tempID == "q value" )
    {
      phase_space_model.Q_value = input_file.get_next_double( );
    }
    else
    {
      Msg::FatalError( "React::reaction::do_phase_space", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }

  if( frame_out != React::CM )
  {
    Msg::FatalError( "React::reaction::do_phase_space",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  phase_space_model.copy_masses( particle_info );
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  phase_space_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::do_gen_evap ----------------------
// Reads and processes data for gen_evap spectra
void React::reaction::do_gen_evap( Dpar::data_parser &input_file )
{
  frame_out = React::LAB;  // default
  for( std::string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    std::string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "g" )
    {
      int num_g = input_file.get_next_int( );
      // read the g values
      gen_evap_model.read_data_interp( input_file, num_g );
    }
    else if( tempID == "theta" )
    {
      int num_theta = input_file.get_next_int( );
      // read the theta values
      gen_evap_model.theta.read_data_interp( input_file, num_theta );
    }
    else if( tempID == "u" )
    {
      double U = input_file.get_next_double( );  // not used
      if( U != 0.0 )
      {
	Msg::Info( "React::reaction::do_gen_evap",
	      "The parameter U is not used in this model" );
      }
    }
    else
    {
      Msg::FatalError( "React::reaction::do_gen_evap", 
		Msg::pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != React::LAB )
  {
    Msg::FatalError( "React::reaction::do_gen_evap",
		"only emission in the laboratory frame is supported" );
  }
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  // now get the transfer matrix
  transfer.allocate( );
  gen_evap_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- React::reaction::write_transfer -----------------
// Writes the output data
void React::reaction::write_transfer( )
{
  // do the proper scaling
  transfer.get_flux_weight( );
  transfer.use_weight( );

  if( cross_section.size( ) > 0 )
  {
    // check the row sums against the integrals
    number_check( transfer.threshold );
    transfer.check_ell0( );
  }
  else
  {
    Msg::Warning( "React::reaction::write_transfer",
		  "Compute the cross section for this reaction" );
  }
  transfer.write_transfer( );
}
// --------------- React::reaction::write_xsec -----------------
// Writes the cross section for gamma data
void React::reaction::write_xsec( )
{
  int num_xsec = cross_section.size( );
  *transfer.output_file << "Cross section: n = " << num_xsec << std::endl;
  *transfer.output_file << "Interpolation: linear-linear" << std::endl;
  transfer.output_file->setf( std::ios::scientific, std::ios::floatfield);
  static int data_precision = Global.Value( "datafield_precision" );
  static int field_width = Global.get_field_width( );
  for( Ddvec::dd_vector::const_iterator xsec = cross_section.begin( ); 
       xsec != cross_section.end( ); ++xsec )
  {
    *transfer.output_file << std::setw(field_width) <<
            std::setprecision(data_precision) <<
      xsec->x << " " << xsec->y << std::endl;
  }
}
// --------------- React::reaction::number_check -----------------
// Evaluates the integrals of multiplicity * cross section * flux
void React::reaction::number_check( double threshold )
{
  React::num_check_param num_check;

  // pointers to the cross section
  num_check.sigma.Eout_interp = cross_section.interp_type;
  Ddvec::dd_vector::const_iterator this_sigma = cross_section.begin( );
  Ddvec::dd_vector::const_iterator next_sigma = this_sigma;
  ++next_sigma;

  // pointers to the multiplicity
  if( multiple.size( ) == 0 )
  {
    multiple.make_flat( cross_section, 1.0 );
  }
  num_check.mult.Eout_interp = multiple.interp_type;
  Ddvec::dd_vector::const_iterator this_mult = multiple.begin( );
  Ddvec::dd_vector::const_iterator next_mult = this_mult;
  ++next_mult;

  // pointers to the flux data
  num_check.e_flux.Eout_interp = transfer.e_flux.interp;
  Lgdata::Flux_List::const_iterator flux_ptr = transfer.e_flux.begin( );
  Lgdata::Flux_List::const_iterator next_flux = flux_ptr;
  ++next_flux;

  // pointers to the model weight
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  num_check.model_weight.Eout_interp = model_weight.interp_type;
  Ddvec::dd_vector::const_iterator this_model_weight = model_weight.begin( );
  Ddvec::dd_vector::const_iterator next_model_weight = this_model_weight;
  ++next_model_weight;

  // count through the incident energy groups
  int Ein_count = 0;
  // which energy group
  std::vector< double >::const_iterator Ein_ptr = transfer.in_groups.begin( );
  std::vector< double >::const_iterator next_Ein = Ein_ptr;
  ++next_Ein;

  // synchronize the first energies
  double E_in0 = ( threshold > *Ein_ptr ) ? threshold :
    *Ein_ptr;
  E_in0 = ( this_sigma->x < E_in0 ) ? E_in0 :
      this_sigma->x;
  E_in0 = ( this_mult->x < E_in0 ) ? E_in0 :
      this_mult->x;
  E_in0 = ( this_model_weight->x < E_in0 ) ? E_in0 :
      this_model_weight->x;
  E_in0 = ( flux_ptr->get_E_in( ) < E_in0 ) ? E_in0 :
      flux_ptr->get_E_in( );

    // we may need to increment num_check.sigma
  for( ; next_sigma->x <= E_in0; 
    this_sigma = next_sigma, ++next_sigma )
  {
    if( next_sigma == cross_section.end( ) )
    {
      Msg::Warning( "React::reaction::number_check",
        "omit reaction: all energy groups too high" );
      transfer.zero_transfer( );
    }
  }
  num_check.sigma.set_pair( *this_sigma, *next_sigma );
  
  // we may need to increment num_check.mult
  for( ; next_mult->x <= E_in0; 
    this_mult = next_mult, ++next_mult )
  {
    if( next_mult == multiple.end( ) )
    {
      Msg::Warning( "React::reaction::number_check",
        "omit reaction: no particle multiplicity given" );
      transfer.zero_transfer( );
    }
  }
  num_check.mult.set_pair( *this_mult, *next_mult );
  
  // we may need to increment the energy bin
  for( ; *next_Ein <= E_in0;
       ++Ein_count, Ein_ptr = next_Ein, ++next_Ein )
  {
    if( next_Ein == transfer.in_groups.end( ) )
    {
      Msg::Warning( "React::reaction::number_check", 
        "omit reaction: threshold too high" );
      transfer.zero_transfer( );
    }
  }
  
  // we may need to increment the incident flux
  for( ; next_flux->get_E_in( ) <= E_in0;
       flux_ptr = next_flux, ++next_flux )
  {
    if( next_flux == transfer.e_flux.end( ) )
    {
      Msg::Warning( "React::reaction::number_check",
        "omit reaction: ran out of flux data" );
      transfer.zero_transfer( );
    }
  }
  num_check.e_flux.first.x = flux_ptr->get_E_in( );
  num_check.e_flux.first.y = flux_ptr->data[ 0 ];
  num_check.e_flux.second.x = next_flux->get_E_in( );
  num_check.e_flux.second.y = next_flux->data[ 0 ];

  // we may need to increment the model weight
  for( ; next_model_weight->x <= E_in0; 
    this_model_weight = next_model_weight, ++next_model_weight )
  {
    if( next_model_weight == model_weight.end( ) )
    {
      Msg::Warning( "React::reaction::number_check",
        "omit reaction: no model weight given" );
      transfer.zero_transfer( );
    }
  }
  num_check.model_weight.set_pair( *this_model_weight, *next_model_weight );

  static double abs_tol = Global.Value( "tight_tol" );
  // now do the integrals
  for( ; ; )
  {
    // the range of integration
    double E_in1;
    E_in1 = ( next_sigma->x > *next_Ein ) ? *next_Ein :
      next_sigma->x;
    E_in1 = ( next_mult->x > E_in1 ) ? E_in1 :
      next_mult->x;
    E_in1 = ( next_model_weight->x > E_in1 ) ? E_in1 :
      next_model_weight->x;
    E_in1 = ( next_flux->get_E_in( ) > E_in1 ) ? E_in1 :
      next_flux->get_E_in( );

    // Use Simpson's rule
    double Simpson = ( ( E_in1 - E_in0 ) / 6.0 ) *
      ( num_check.value( E_in0 ) +
	4.0 * num_check.value( 0.5 * ( E_in0 + E_in1 ) ) +
	num_check.value( E_in1 ) );
    transfer.row_checks[ Ein_count ].weight_1[ 0 ] += Simpson;

    // got to the next interval
    while( next_sigma->x < ( 1.0 + abs_tol ) * E_in1 )
    {
      // increment num_check.sigma
      this_sigma = next_sigma;
      ++next_sigma;
      if( next_sigma == cross_section.end( ) )
      {
        goto slut;
      }
      num_check.sigma.set_pair( *this_sigma, *next_sigma );
    }
    while( next_mult->x < ( 1.0 + abs_tol ) * E_in1 )
    {
      // increment num_check.mult
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == multiple.end( ) )
      {
        goto slut;
      }
      num_check.mult.set_pair( *this_mult, *next_mult );
    }
    while( next_model_weight->x < ( 1.0 + abs_tol ) * E_in1 )
    {
      // increment num_check.model_weight
      this_model_weight = next_model_weight;
      ++next_model_weight;
      if( next_model_weight == model_weight.end( ) )
      {
        goto slut;
      }
      num_check.model_weight.set_pair( *this_model_weight, *next_model_weight );
    }
    while( *next_Ein < ( 1.0 + abs_tol ) * E_in1 )
    {
      // go to the next incident energy bin
      ++Ein_count;
      Ein_ptr = next_Ein;
      ++next_Ein;
      if( next_Ein == transfer.in_groups.end( ) )
      {
        goto slut;
      }
    }
    while( next_flux->get_E_in( ) < ( 1.0 + abs_tol ) * E_in1 )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == transfer.e_flux.end( ) )
      {
        goto slut;
      }
      num_check.e_flux.first.x = flux_ptr->get_E_in( );
      num_check.e_flux.first.y = flux_ptr->data[ 0 ];
      num_check.e_flux.second.x = next_flux->get_E_in( );
      num_check.e_flux.second.y = next_flux->data[ 0 ];
    }
    E_in0 = E_in1;
  }
 slut:
  // scale the row sums by the flux weights
  transfer.scale_row_check( );
}
