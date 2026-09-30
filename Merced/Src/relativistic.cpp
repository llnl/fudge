/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: relativistic.cpp 1 2015-05-20 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

#include <fstream>  // standard file stream package
#include <iostream>
#include <cmath>

#include "relativistic.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"


// implementation of the classes used in relativistic mechanics

// ********* class Rel::Ep_vector ******************************
// ------------- Rel::Ep_vector::get_p --------------------
// gets the length of the momentum
double Rel::Ep_vector::get_p( )
{
  double laengd = Ep[ 1 ] * Ep[ 1 ] + Ep[ 2 ] * Ep[ 2 ] +
    Ep[ 3 ] * Ep[ 3 ];
  return std::sqrt( laengd );
}
// ********* class Rel::relativistic_masses **********************
// ------------- Rel::relativistic_masses::setup_masses --------------------
// Saves the rest masses
void Rel::relativistic_masses::setup_masses( Maps::particleInfo *to_save, double *file_Q )
{
  rest_masses = to_save;
  Q_value = *file_Q;

  static double etol = Global.Value( "looser_tol" );

  // check for photon as residual
  if( rest_masses->mRes == 0.0 )
  {
    double true_Q = ( rest_masses->mProj + rest_masses->mTarg ) -
      rest_masses->mProd;
    if( std::abs( true_Q - Q_value ) > etol*std::abs( true_Q ) )
    {
      Msg::Warning( "Rel::relativistic_masses::setup_masses",
		    "Q value inconsistent with photon residual" );
    }
    Q_value = true_Q;
    *file_Q = true_Q;
  }
  else
  {
    // calculate the mass of the residual from the Q value
    double true_mRes;
    if( rest_masses->mProd == rest_masses->mProj )
    {
      // elastic or inelastic collision
      true_mRes = rest_masses->mTarg - Q_value;
    }
    else
    {
      // knock-on reaction
      true_mRes = ( rest_masses->mTarg - Q_value ) + 
        ( rest_masses->mProj - rest_masses->mProd );
    }
    // check the input mass of residual
    if( ( rest_masses->mRes > 0.0 ) && ( std::abs( true_mRes - rest_masses->mRes ) > etol*true_mRes ) )
    {
      Msg::Warning( "Rel::relativistic_masses::setup_masses",
		    "mass of residual inconsistent with Q value" );
    }
    to_save->mRes = true_mRes;
  }

  // total mass of all particles
  M_total = rest_masses->mTarg + rest_masses->mRes + 
       rest_masses->mProj + rest_masses->mProd;
}
// ------------- Rel::relativistic_masses::get_threshold --------------------
// Calculates the threshold kinetic energy
void Rel::relativistic_masses::get_threshold( )
{
  if( Q_value >= 0.0 )
  {
    threshold = 0.0;
  }
  else
  {
    threshold = - M_total * Q_value / ( 2.0*rest_masses->mTarg );
  }
}

// ********* class Rel::relativistic_map **********************
// ------------- Rel::relativistic_map::set_boost --------------------
// Sets up the boost to the lab frame
void Rel::relativistic_map::set_boost( double T_in_lab )
{
  Tin_lab = T_in_lab;
  pin_lab = relativistic_F::p_from_T( Tin_lab, relMasses->rest_masses->mProj );
  double initial_rest_energy = relMasses->rest_masses->mTarg +
    relMasses->rest_masses->mProj;
  Minkowski = initial_rest_energy * initial_rest_energy +
    2.0 * T_in_lab * relMasses->rest_masses->mTarg;
  Minkowski = sqrt( Minkowski );
  pin_cm = pin_lab*relMasses->rest_masses->mTarg / Minkowski;
  sinh_chi = pin_cm/relMasses->rest_masses->mTarg;
  cosh_chi = sqrt( 1.0 + sinh_chi*sinh_chi );
  get_p_cm_out( );
}
// ------------- Rel::relativistic_map::boost --------------------
// Boosts the ejected particle from the center-of-mass to the lab frame
void Rel::relativistic_map::boost( double mu_cm, double *T_lab, double *p_lab_parallel )
{
  double E_cm = relMasses->rest_masses->mProd + T_cm_out;
  double p_cm_parallel = p_cm_out*mu_cm;
  double p_cm_perpendicular = p_cm_out*sqrt( 1.0 - mu_cm*mu_cm );
  *p_lab_parallel = sinh_chi*E_cm + cosh_chi*p_cm_parallel;
  double p_lab_sq = ( *p_lab_parallel )*( *p_lab_parallel ) +
    p_cm_perpendicular * p_cm_perpendicular;
  *T_lab = relativistic_F::T_from_p( sqrt( p_lab_sq ),
				     relMasses->rest_masses->mProd );
}
// ------------- Rel::relativistic_map::get_p_cm_out --------------------
// Calculates the center-of-mass energy and momentum for discrete 2-body reactions
void Rel::relativistic_map::get_p_cm_out( )
{
  double temp = relMasses->M_total*relMasses->Q_value + 2*relMasses->rest_masses->mTarg* Tin_lab;
  if( temp < 0.0 ) temp = 0.0;
  p_cm_out = sqrt( temp*( temp + 4*relMasses->rest_masses->mRes*relMasses->rest_masses->mProd ) ) /
    ( 2.0 * Minkowski );
  T_cm_out = relativistic_F::T_from_p( p_cm_out,
				       relMasses->rest_masses->mProd );
}
// ------------- Rel::relativistic_map::get_T_lab_out --------------------
// Calculates the lab-frame kinetic energy for discrete 2-body reactions
double Rel::relativistic_map::get_T_lab_out( double mu_cm )
{
  // boost to lab frame
  double p_lab_parallel;
  double Tout_lab;
  boost( mu_cm, &Tout_lab, &p_lab_parallel );
  return Tout_lab;
}
// ------------- Rel::relativistic_map::get_mu_cm --------------------
// Gets the value of mu_cm given the laboratory kinetic energy of outgoing double.
double Rel::relativistic_map::get_mu_cm( double T_lab )
{
  // This routine uses the fact that the equation for mu_cm is a quadratic
  double p_lab_sq = T_lab*( 2*relMasses->rest_masses->mProd + T_lab );
  double A = p_cm_out*sinh_chi;
  A *= A;  // square it
  double E_prod = relMasses->rest_masses->mProd + T_cm_out;
  double B = 2*E_prod*p_cm_out*sinh_chi*cosh_chi;
  double C = E_prod*sinh_chi;
  C = C*C + p_cm_out*p_cm_out - p_lab_sq;
  double root_1;
  double root_2;
  int num_roots = math_F::quadratic( A, B, C, &root_1, &root_2 );
  if( num_roots != 2 )
  {
    Msg::FatalError( "Rel::relativistic_map::get_mu_cm",
		     "wrong number of roots" );
  }

  static double mu_tol = Global.Value( "looser_tol" );
  if( root_2 > 1.0 )
  {
    if( root_2 > 1.0 + mu_tol )
    {
      Msg::Warning ( "Rel::relativistic_map::get_mu_cm",
		     Msg::pastenum( "mu too big: ", root_2 ) );
    }
    root_2 = 1.0;
  }
  if( root_2 < -1.0 )
  {
    if( ( root_2 < -1.0 - mu_tol ) && ( Tin_lab - relMasses->threshold < mu_tol*T_lab ) )
    {
      Msg::Warning ( "Rel::relativistic_map::get_mu_cm",
		     Msg::pastenum( "mu too small: ", root_2 ) );
    }
    root_2 = -1.0;
  }
  return root_2;
}
// ------------- Rel::relativistic_map::get_E_mu_lab --------------------
// Gets the laboratory kinetic energy and cosine of outgoing particle
void Rel::relativistic_map::get_E_mu_lab( double mu_cm, double *Tout_lab, double *mu_lab )
{
  double p_cm_perpendicular = p_cm_out*sqrt( 1.0 - mu_cm*mu_cm );
  double p_lab_parallel;
  boost( mu_cm, Tout_lab, &p_lab_parallel );
  *mu_lab = p_lab_parallel/sqrt( p_lab_parallel*p_lab_parallel +
    p_cm_perpendicular*p_cm_perpendicular );
}
// ------------- Rel::relativistic_map::Newton_min_Eout --------------------
// For exothermic reactions, returns the Newtonian incident energy for minimal outgoing lab kinetic energy
double Rel::relativistic_map::Newton_min_Eout( )
{
  // For mu = -1 find E_in such that dV_0/dE_in = dvcm_out/dE_in
  double mass_in = relMasses->rest_masses->mTarg + relMasses->rest_masses->mProj;
  double gamma = relMasses->rest_masses->mProj * relMasses->rest_masses->mProd /( mass_in * mass_in );
  double beta_eject = relMasses->rest_masses->mRes / ( relMasses->rest_masses->mProd + relMasses->rest_masses->mRes );
  double alpha = beta_eject * relMasses->rest_masses->mTarg / mass_in;

  double E_min_in = 0.0;

  if( relMasses->Q_value * ( alpha - gamma ) > 0.0 )
  {
    E_min_in = - gamma * beta_eject * relMasses->Q_value / ( alpha * ( gamma - alpha ) );
  }
  else if( relMasses->Q_value * ( alpha - gamma ) < 0.0 )
  {
    E_min_in = beta_eject * relMasses->Q_value / ( gamma - alpha );
  }
  else //  if( relMasses->Q_value * ( alpha - gamma ) == 0.0 )
  {
    E_min_in = relMasses->threshold;
  }

  // don't go below the threshold
  if( E_min_in < relMasses->threshold )
  {
    E_min_in = relMasses->threshold;
  }

  return E_min_in;
}
// ------------- Rel::relativistic_map::Newton_zero_Eout --------------------
// Returns the incident kinetic energy for zero outgoing lab-frame kinetic energy, Newtonian
double Rel::relativistic_map::Newton_zero_Eout( )
{
  double mass_in = relMasses->rest_masses->mTarg + relMasses->rest_masses->mProj;
  double gamma = relMasses->rest_masses->mProj * relMasses->rest_masses->mProd /( mass_in * mass_in );
  double beta_eject = relMasses->rest_masses->mRes / ( relMasses->rest_masses->mProd + relMasses->rest_masses->mRes );
  double alpha = beta_eject * relMasses->rest_masses->mTarg / mass_in;

  double E_zero_in;

  // special for gamma as the residual
  if( relMasses->rest_masses->mRes == 0.0 )
  {
    E_zero_in = relMasses->Q_value * relMasses->Q_value / ( 2 * mass_in * gamma );
  }
  else if( gamma == alpha )
  {
    E_zero_in = 0.0;
  }
  else
  {
    E_zero_in = -beta_eject*relMasses->Q_value/( alpha - gamma );
  }
  return E_zero_in;
}

// ************************ clsss Rel::relativistic_param ******************
// ------------- Rel::relativistic_param::zero_Eout --------------------
// Returns the incident kinetic energy for zero outgoing lab-frame kinetic energy, relativistic
double Rel::relativistic_param::zero_Eout( )
{
  // set up the void* parameters for the root finder
  void *params = static_cast< void * >( this );

  // get the incident kinetic energy for the Newtonian approximation (It's low)
  double T_in = relMap.Newton_zero_Eout( );
  if( T_in <= 0.0 )
  {
    return 0.0;
  }

  mu_cm = -1.0;
  double p_out = relativistic_F::p_out_lab( T_in, params );
  Ddvec::dd_entry low_Tin( T_in, p_out );
  // a high bound
  T_in *= 1.1;
  p_out = relativistic_F::p_out_lab( T_in, params );
  Ddvec::dd_entry high_Tin( T_in, p_out );

  // for the root finder
  static double E_tol = Global.Value( "looser_tol" );
  double this_Tin = math_F::zeroin( relativistic_F::p_out_lab, 0.0,
				    low_Tin, high_Tin, params, E_tol );
  return this_Tin;
}
// ------------- Rel::relativistic_param::find_hit --------------------
// Finds the incident kinetic energy for given T_lab and mucm.
// Returns the solutions.
double Rel::relativistic_param::find_hit( double T_lab, double mucm, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  mu_cm = mucm;
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  Ddvec::dd_entry use_pair_0 = pair_0;
  Ddvec::dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "looser_tol" );
  double this_Tin = math_F::zeroin( relativistic_F::T_out_lab, T_lab,
				    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}
// ------------- Rel::relativistic_param::find_bottom --------------------
// For mu < 0 find the incident kinetic energy which minimizes Tout
// On exit pair_0 and pair_1 are tighter bounds.
Ddvec::dd_entry Rel::relativistic_param::find_bottom( double mu,
						 const Ddvec::dd_entry &pair_in_0,
					  const Ddvec::dd_entry &pair_in_1, double tol )
{
  if( ( mu >= 0.0 ) || ( relMap.relMasses->Q_value == 0.0 ) )
  {
    return Ddvec::dd_entry( 0.0, 0.0 );
  }

  // start at the extremes
  Ddvec::dd_entry pair_0( pair_in_0 );
  Ddvec::dd_entry pair_1( pair_in_1 );
  // parameters for the function relativistic_F::T_out_lab
  mu_cm = mu;
  void *params = static_cast< void *>( this );
  // This is also the result if the initial bounds already meet the tolerance.
  double Ein = 0.5 * ( pair_0.x + pair_1.x );
  double Tout = relativistic_F::T_out_lab( Ein, params );
  
  while( std::abs( pair_1.x - pair_0.x ) > tol * pair_1.x )
  {
    Ein = math_F::parabola_bottom( relativistic_F::T_out_lab, pair_0,
					  pair_1, params );
    // don't go above or below the initial bounds
    double Ein_mid = 0.5 * ( pair_1.x + pair_0.x );
    if( Ein < pair_in_0.x )
    {
      Ein = 0.5 * ( Ein_mid + pair_in_0.x );
    }
    else if( Ein > pair_in_1.x )
    {
      Ein = 0.5 * ( Ein_mid + pair_in_1.x );
    }
    Tout = relativistic_F::T_out_lab( Ein, params );

    // set up the next iteration
    if( Ein < pair_0.x )
    {
      pair_1.x = pair_0.x;
      pair_1.y = pair_0.y;
      pair_0.x = Ein;
      pair_0.y = Tout;
    }
    else if( Ein <  pair_1.x )
    {
      if( pair_0.y < pair_1.y )
      {
	pair_1.x = Ein;
	pair_1.y = Tout;
      }
      else
      {
	pair_0.x = Ein;
	pair_0.y = Tout;
      }
    }
    else // Ein >=  pair_1.x
    {
      pair_0.x = pair_1.x;
      pair_0.y = pair_1.y;
      pair_1.x = Ein;
      pair_1.y = Tout;
    }
  }
  Ddvec::dd_entry ans( Ein, Tout );
  return ans;
}

// ********* class Rel::two_step_relativistic_masses **********************
// ------------- Rel::two_step_relativistic_masses::set_masses ----------
// Sets up the map from center-of-mass to laboratory coordinates
void Rel::two_step_relativistic_masses::set_masses( Maps::particleInfo *step1_particles,
   double file_Q1, Maps::particleInfo *step2_particles, double file_Q2 )
{
  static double etol = Global.Value( "looser_tol" );

  // information for step 1
  step1Masses.rest_masses = step1_particles;
  step1Masses.Q_value = file_Q1;
  
  // Check the mass of the excited product
  double true_mProd = step1_particles->mProj + step1_particles->mTarg -
    step1_particles->mRes - file_Q1;

  if( ( step1_particles->mProd > 0.0 ) &&
      ( std::abs( step1_particles->mProd - true_mProd ) > etol*true_mProd ) )
  {
    Msg::Warning( "Rel::two_step_relativistic_masses::set_masses",
	     "mass of product from step 1 inconsistent" );
  }
  step1_particles->mProd = true_mProd;

  // get the total mass for step 1
  step1Masses.M_total = step1_particles->mTarg + step1_particles->mProj +
    step1_particles->mProd + step1_particles->mRes;

  // set up step 2
  step2_particles->mProj = 0.0;
  step2_particles->mTarg = step1_particles->mProd;

  // If file_Q2 is given in the data, check it.
  double true_Q = step2_particles->mTarg -
    ( step2_particles->mProd +
      step2_particles->mRes );

  if( file_Q2 > 0.0 )
  {
    if( std::abs( true_Q - file_Q2 ) > etol*true_Q )
    {
      Msg::Warning( "Rel::two_step_relativistic_masses::set_masses",
		    "Q value of step 2 inconsistent " );
    }
  }
  
  step2Masses.rest_masses = step2_particles;
  step2Masses.Q_value = true_Q;

  // get the total mass for step 2
  step2Masses.M_total = step2_particles->mTarg +
    step2_particles->mProd + step2_particles->mRes;
}
// ------------- Rel::two_step_relativistic_masses::get_p_cm2 ----------
// Returns the outgoing center-of-mass momentum for step 2
double Rel::two_step_relativistic_masses::get_p_cm2( )
{
  double S_2 = step2Masses.rest_masses->mTarg;
  double temp = step2Masses.M_total * step2Masses.Q_value;
  double p_cm2_out = std::sqrt( temp*
			     ( temp + 4*step2Masses.rest_masses->mRes *
			       step2Masses.rest_masses->mProd ) ) / ( 2.0*S_2);

  return p_cm2_out;
}

// ********* class Rel::two_step_relativistic_map **********************
// ------------- Rel::two_step_relativistic_map::set_boost_Theta ----------
//! Sets up the boost and rotation for step 1
void Rel::two_step_relativistic_map::set_boost_Theta( double T_in, double mucm_1 )
{
  double Tout_lab1;
  double mu_lab;

  step1Map.relMasses = &twoStepMasses->step1Masses;
  step1Map.set_boost( T_in );
  step1Map.get_E_mu_lab( mucm_1, &Tout_lab1, &mu_lab );

  // the step-1 boost
  double ratio = Tout_lab1 / twoStepMasses->step1Masses.rest_masses->mProd;
  cosh_chi1 = 1.0 + ratio;
  sinh_chi1 = std::sqrt( ratio * ( 2.0 + ratio ) );

  // save the rotation
  mu_lab1 = mu_lab;
}
// ------------- Rel::two_step_relativistic_map::do_boost_Theta ----------
// Does the boost to the lab frame, without the rotation
void Rel::two_step_relativistic_map::do_boost_Theta( const Rel::Ep_vector &Ep_cm,
						Rel::Ep_vector *Ep_lab )
{
  Ep_lab->Ep[ 0 ] = cosh_chi1 * Ep_cm.Ep[ 0 ] + sinh_chi1 * Ep_cm.Ep[ 1 ];
  Ep_lab->Ep[ 1 ] = sinh_chi1 * Ep_cm.Ep[ 0 ] + cosh_chi1 * Ep_cm.Ep[ 1 ];
  Ep_lab->Ep[ 2 ] = Ep_cm.Ep[ 2 ];
  Ep_lab->Ep[ 3 ] = Ep_cm.Ep[ 3 ];
}
// ------------- Rel::two_step_relativistic_map::rotate_Theta ----------
// Does the rotation to coordinants parallel to the incident particle
void Rel::two_step_relativistic_map::rotate_Theta( const Rel::Ep_vector &Ep_lab1,
					      Rel::Ep_vector *Ep_lab )
{
  double sin_phi = std::sqrt( 1.0 - mu_lab1 * mu_lab1 );
  
  Ep_lab->Ep[ 0 ] = Ep_lab1.Ep[ 0 ];
  Ep_lab->Ep[ 1 ] = mu_lab1 * Ep_lab1.Ep[ 1 ] - sin_phi * Ep_lab1.Ep[ 2 ];
  Ep_lab->Ep[ 2 ] = sin_phi * Ep_lab1.Ep[ 1 ] + mu_lab1 * Ep_lab1.Ep[ 2 ];
  Ep_lab->Ep[ 3 ] = Ep_lab1.Ep[ 3 ];
}
// ------------- Rel::two_step_relativistic_map::two_step_start_E_lab ----------
// Initial steps for getting the laboratory kinetic energy of the outgoing particle for 2-step reactions
void Rel::two_step_relativistic_map::two_step_start_E_lab( double T_in, double mucm_1,
						      double mucm_2, double w,
						      Rel::Ep_vector *Ep_lab )
{
  // for step 1
  set_boost_Theta( T_in, mucm_1 );
  
  // for step 2
  double p_cm2_out = twoStepMasses->get_p_cm2( );
  double T_cm2_out = relativistic_F::T_from_p( p_cm2_out,
			twoStepMasses->step2Masses.rest_masses->mProd );

  // do the boost
  // At O' in coordinates aligned with the outgoing particle of step 1 we have
  Rel::Ep_vector Ep_cm;
  Ep_cm.Ep[ 0 ] = twoStepMasses->step2Masses.rest_masses->mProd + T_cm2_out;
  Ep_cm.Ep[ 1 ] = p_cm2_out * mucm_2;
  double sin_theta2 = std::sqrt( 1.0 - mucm_2*mucm_2 );
  double cos_w = std::cos( w );
  double sin_w = std::sin( w );
  Ep_cm.Ep[ 2 ] = p_cm2_out * sin_theta2 * cos_w;
  Ep_cm.Ep[ 3 ] = p_cm2_out * sin_theta2 * sin_w;

  // boost to O, parallel to outgoing particle of step 1
  do_boost_Theta( Ep_cm, Ep_lab );
}
// ------------- Rel::two_step_relativistic_map::two_step_get_E_lab ----------
// Gets the laboratory kinetic energy of the outgoing particle for 2-step reactions
double Rel::two_step_relativistic_map::two_step_get_E_lab( double T_in, double mucm1,
						      double mucm2 )
{
  Rel::Ep_vector Ep_lab;

  // get the energy and momentum in the rotated lab frame
  two_step_start_E_lab( T_in, mucm1, mucm2, 0.0, &Ep_lab );
  double p_lab = Ep_lab.get_p( );

  // The result
  double Tout_lab = relativistic_F::T_from_p( p_lab,
			 twoStepMasses->step2Masses.rest_masses->mProd );
  return Tout_lab;
}
// ------------- Rel::two_step_relativistic_map::two_step_get_mu_lab ----------
// Gets the laboratory kinetic energy and cosine of outgoing particle for 2-step reactions
double Rel::two_step_relativistic_map::two_step_get_mu_lab( double T_in, double mucm_1,
							  double mucm_2, double w )
{
  Rel::Ep_vector Ep_lab;

  // get the energy and momentum in the rotated lab frame
  two_step_start_E_lab( T_in, mucm_1, mucm_2, w, &Ep_lab );

  // rotate to the direction of the incident particle
  Rel::Ep_vector Ep_lab_out;
  rotate_Theta( Ep_lab, &Ep_lab_out );
  double p_lab = Ep_lab_out.get_p( );

  // The results
  double mu_lab = Ep_lab_out.Ep[ 1 ] / p_lab;
  return mu_lab;
}
// ------------- Rel::two_step_relativistic_map::two_step_get_E_mu_lab ----------
// Gets the laboratory kinetic energy and cosine of outgoing particle for 2-step reactions
void Rel::two_step_relativistic_map::two_step_get_E_mu_lab( double T_in, double mucm_1,
						       double mucm_2, double w,
			      double *Tout_lab, double *mu_lab )
{
  Rel::Ep_vector Ep_lab;

  // get the energy and momentum in the rotated lab frame
  two_step_start_E_lab( T_in, mucm_1, mucm_2, w, &Ep_lab );

  // rotate to the direction of the incident particle
  Rel::Ep_vector Ep_lab_out;
  rotate_Theta( Ep_lab, &Ep_lab_out );
  double p_lab = Ep_lab_out.get_p( );

  // The results
  *Tout_lab = relativistic_F::T_from_p( p_lab,
		      twoStepMasses->step2Masses.rest_masses->mProd );
  *mu_lab = Ep_lab_out.Ep[ 1 ] / p_lab;
}

// ********* class Rel::relativistic_2_step_param **********************
// ------------- Rel::relativistic_2_step_param::get_mucm1 ----------
// Returns the direction cosine for the first step for given incident kinetic energy and desired outgoing kinetic energy
double Rel::relativistic_2_step_param::get_mucm1( double T_in, double mucm2, double Tout_lab )
{
  // set up the void* parameters for the root finder
  Tin_lab = T_in;
  mu_cm_2 = mucm2;
  void *params = static_cast< void * >( this );

  // the initial pairs
  Ddvec::dd_entry pair_0;
  Ddvec::dd_entry pair_1;
  pair_0.x = -1.0;
  pair_0.y = relativistic_F::mucm1_T_out_lab( -1.0, params );
  pair_1.x = 1.0;
  pair_1.y = relativistic_F::mucm1_T_out_lab( 1.0, params );
  
  static double E_tol = Global.Value( "looser_tol" );
  double mucm1 = math_F::zeroin( relativistic_F::mucm1_T_out_lab, Tout_lab,
				    pair_0, pair_1, params, E_tol );
  return mucm1;
}
// ------------- Rel::relativistic_2_step_param::get_mucm2 ----------
// Returns the direction cosine for the second step for desired outgoing kinetic energy
double Rel::relativistic_2_step_param::get_mucm2( double Tout_lab )
{
  // set up the void* parameters for the root finder
  void *params = static_cast< void * >( this );

  // the initial pairs
  Ddvec::dd_entry pair_0;
  Ddvec::dd_entry pair_1;
  pair_0.x = -1.0;
  pair_0.y = relativistic_F::mucm2_T_out_lab( -1.0, params );
  pair_1.x = 1.0;
  pair_1.y = relativistic_F::mucm2_T_out_lab( 1.0, params );
  
  static double E_tol = Global.Value( "looser_tol" );
  double mucm2 = math_F::zeroin( relativistic_F::mucm2_T_out_lab, Tout_lab,
				    pair_0, pair_1, params, E_tol );
  return mucm2;
}
// ------------- Rel::relativistic_2_step_param::find_hit ----------
// Finds the incident energy for given T_out_lab
double Rel::relativistic_2_step_param::find_hit( double T_out_lab,
                    const Ddvec::dd_entry &pair_0, const Ddvec::dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  Ddvec::dd_entry use_pair_0 = pair_0;
  Ddvec::dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "looser_tol" );
  double this_Tin = math_F::zeroin( relativistic_F::two_step_T_out_lab, T_out_lab,
                                    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}
// ------------- Rel::relativistic_2_step_param::find_bottom --------------------
// Find the incident kinetic energy which minimizes Tout
// On exit pair_0 and pair_1 are tighter bounds.
Ddvec::dd_entry Rel::relativistic_2_step_param::find_bottom( const Ddvec::dd_entry &pair_in_0,
  const Ddvec::dd_entry &pair_in_1, double tol )
{
  // start at the extremes
  Ddvec::dd_entry pair_0( pair_in_0 );
  Ddvec::dd_entry pair_1( pair_in_1 );
  
  // parameters for the function relativistic_F::T_out_lab
  void *params = static_cast< void *>( this );
  // This is also the result if the initial bounds already meet the tolerance.
  double Ein = 0.5 * ( pair_0.x + pair_1.x );
  double Tout = relativistic_F::two_step_T_out_lab( Ein, params );
  while( std::abs( pair_1.x - pair_0.x ) > tol * pair_1.x )
  {
    Ein = math_F::parabola_bottom( relativistic_F::two_step_T_out_lab, pair_0,
					  pair_1, params );
    // don't go below the threshold
    if( Ein < relTwoStepMap.twoStepMasses->step1Masses.threshold )
    {
      Ein = relTwoStepMap.twoStepMasses->step1Masses.threshold;
    }
    Tout = relativistic_F::two_step_T_out_lab( Ein, params );

    // set up the next iteration
    if( Ein <  pair_0.x )
    {
      pair_1.x = pair_0.x;
      pair_1.y = pair_0.y;
      pair_0.x = Ein;
      pair_0.y = Tout;
    }
    else if( Ein <  pair_1.x )
    {
      if( pair_0.y < pair_1.y )
      {
	pair_1.x = Ein;
	pair_1.y = Tout;
      }
      else
      {
	pair_0.x = Ein;
	pair_0.y = Tout;
      }
    }
    else // Ein >=  pair_1.x
    {
      pair_0.x = pair_1.x;
      pair_0.y = pair_1.y;
      pair_1.x = Ein;
      pair_1.y = Tout;
    }
  }
  Ddvec::dd_entry ans( Ein, Tout );
  return ans;
}
// ------------- Rel::relativistic_2_step_param::find_lowest_bottom --------------------
// Find the incident kinetic energy which minimizes Eout for mucm1 = mucm2 = -1
Ddvec::dd_entry Rel::relativistic_2_step_param::find_lowest_bottom( )
{
  static double etol = Global.Value( "looser_tol" );

  Ddvec::dd_entry ans;
  // pairs for the root finder, zeroin
  Ddvec::dd_entry pair_0; // ( Ein, Tout ) for Ein below the minimum
  Ddvec::dd_entry pair_1; // ( Ein, Tout ) for Ein above the minimum

  // parameters for the function relativistic_F::T_out_lab
  void *params = static_cast< void *>( this );

  // This coding assumes that the outgoing kinetic energy is never zero.

  // use the root finder
  pair_0.x = relTwoStepMap.twoStepMasses->step1Masses.threshold;
  pair_0.y = relativistic_F::two_step_T_out_lab( pair_0.x, params );
  pair_1.x = relTwoStepMap.step1Map.Newton_min_Eout( );
  mu_cm_1 = -1.0;
  mu_cm_2 = -1.0;
  pair_1.y = relativistic_F::two_step_T_out_lab( pair_1.x, params );
  ans = find_bottom( pair_0, pair_1, etol );

  return ans;
}

// ***************** functions ********************
// ------------- relativistic_F::p_from_T --------------------
// Gets the momentum from the kinetic energy
double relativistic_F::p_from_T( double T, double E0 )
{
  double p = sqrt( T*( 2*E0 + T ) );
  return p;
}
// ------------- relativistic_F::T_from_p ------------------
// Gets the momentum from the kinetic energy
double relativistic_F::T_from_p( double p, double E0 )
{
  double p_sq = p*p;
  double T = p_sq/( E0 + sqrt( E0*E0 + p_sq ) );
  return T;
}
// ------------- relativistic_F::T_out_lab --------------------
// Returns the lab-frame kinetic energy of the outgoing particle
double relativistic_F::T_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a Rel::relativistic_param *
  Rel::relativistic_param *map =
    static_cast< Rel::relativistic_param * >( params );
  // Sets up the boost to the lab frame
  map->relMap.set_boost( T_lab_in );
  double T = map->relMap.get_T_lab_out( map->mu_cm );
  return T;
}

// ------------- relativistic_F::p_out_lab --------------------
// Returns the lab-frame  parallel momentum of the outgoing particle
double relativistic_F::p_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a Rel::relativistic_param *
  Rel::relativistic_param *map =
    static_cast< Rel::relativistic_param * >( params );
  // Sets up the boost to the lab frame
  map->relMap.set_boost( T_lab_in );
  double p_parallel;  // parallel momentum in lab frame
  double T_lab;  // lab frame kinetic energy
  map->relMap.boost( map->mu_cm, &T_lab, &p_parallel );
  return p_parallel;
}

// ------------- relativistic_F::two_step_T_out_lab --------------------
// For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of T_in
double relativistic_F::two_step_T_out_lab( double T_in, void *params )
{
  // The parameters are really a Rel::relativistic_2_step_param *
  Rel::relativistic_2_step_param *map =
    static_cast< Rel::relativistic_2_step_param * >( params );

  double T_out_lab = map->relTwoStepMap.two_step_get_E_lab( T_in,
				map->mu_cm_1, map->mu_cm_2 );
  return T_out_lab;
}

// ------------- relativistic_F::mucm1_T_out_lab --------------------
// For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of mucm_1
double relativistic_F::mucm1_T_out_lab( double mucm_1, void *params )
{
  // The parameters are really a Rel::relativistic_2_step_param *
  Rel::relativistic_2_step_param *map =
    static_cast< Rel::relativistic_2_step_param * >( params );

  double T_out_lab = map->relTwoStepMap.two_step_get_E_lab( map->Tin_lab,
	       mucm_1, map->mu_cm_2 );
  return T_out_lab;
}

// ------------- relativistic_F::mucm2_T_out_lab --------------------
// For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of mucm_2
double relativistic_F::mucm2_T_out_lab( double mucm_2, void *params )
{
  // The parameters are really a Rel::relativistic_2_step_param *
  Rel::relativistic_2_step_param *map =
    static_cast< Rel::relativistic_2_step_param * >( params );

  double T_out_lab = map->relTwoStepMap.
    two_step_get_E_lab( map->Tin_lab,
		    map->mu_cm_1, mucm_2 );
  return T_out_lab;
}
