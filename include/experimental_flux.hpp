#ifndef __experimental_flux_hpp__
#define __experimental_flux_hpp__

#include <string>
#include <vector>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_integration.h>

#include "constants.hpp"
#include "spectral_flux.hpp"

// Constant numbers for precision etc.
//const double radint_abs_prec = 1.0e-1, radint_rel_prec = 1.0e-6;
const double ergint_abs_prec = 0.0, ergint_rel_prec = 1.0e-6;
const int gagg_method = 5, int_space_size = 1e5;

struct exp_setup { int n_bins; double bin_lo; double bin_delta; double r_max; double b_field; double length; std::string eff_exposure_file; };
struct exp_flux_params { double mass; double length; double r_max; double erg; double rad; OneDInterpolator* eff_exp; SolarModel* s; gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; OneDInterpolator* eff_exp; OneDInterpolator* spectral_flux; };


double rho_integrand_gagg(double rho, void * params) {
  // Retrieve parameters and other integration variables.
  struct exp_flux_params * p1 = (struct exp_flux_params *)params;
  double erg = (p1->erg);
  double r = (p1->rad);
  SolarModel *s = p1->s;

  double cylinder = rho*rho - r*r;
  cylinder = rho/sqrt(cylinder);

  return cylinder*(s->Gamma_P_Primakoff(erg, rho));
}

double rad_integrand_gagg(double rad, void * params) {
  struct exp_flux_params * p2 = (struct exp_flux_params *)params;
  p2->rad = rad;
  SolarModel *s = p2->s;
  double r_max = std::min(1.0, s->r_hi);

  gsl_function f1;
  f1.function = &rho_integrand_gagg;
  f1.params = p2;

  double result, error;
  gsl_integration_qag (&f1, rad, r_max, 0.01*ergint_abs_prec, 0.01*ergint_rel_prec, int_space_size, method1, p2->w1, &result, &error);

  result = rad*result;
  return result;
}

double erg_integrand_gagg(double erg, void * params)
{
  const double eVm = 1.0e7*gev2cm;
  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar/(60.0*60.0*24.0*365.0)));
  const double conversion_prob_factor = pow(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm),2);

  struct exp_flux_params * p3 = (struct exp_flux_params *)params;
  p3->erg = erg;
  SolarModel *s = p3->s;
  double r_min = s->r_lo, r_max = std::min(p3->r_max, s->r_hi);
  double m = p3->mass;
  double length = (p3->length)/eVm;

  double argument = 0.25*1.0e-3*length*m*m/erg;
  double sincsq = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = p3->eff_exp->interpolate(erg);

  gsl_function f2;
  f2.function = &rad_integrand_gagg;
  f2.params = p3;

  double spectral_flux, spectral_flux_error;
  gsl_integration_qag (&f2, r_min, r_max, 0.1*ergint_abs_prec, 0.1*ergint_rel_prec, int_space_size, method1, p3->w2, &spectral_flux, &spectral_flux_error);

  return 0.5*gsl_pow_2(erg/pi)*exposure*(1.0e-4*factor*spectral_flux/365.0)*conversion_prob_factor*sincsq;
}

double erg_integrand_gagg_from_file(double erg, void * params) {
  const double eVm = gev2cm*1.0E7;
  struct exp_flux_params_file * p = (struct exp_flux_params_file *)params;
  double m = p->mass;
  double length = (p->length)/eVm;

  double argument = 0.25*1.0e-3*length*m*m/erg;
  double temp = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = p->eff_exp->interpolate(erg);
  // N.B. Here we assume axion is massless in stellar interior:
  double exp_flux = p->spectral_flux->interpolate(erg);

  return temp*exposure*exp_flux;
};

//std::vector<double> axion_photon_counts (double mass, double gagg, std::vector<double> bins, bool save_output, std::string output_path);
std::vector<double> axion_photon_counts (double mass, double gagg, exp_setup *setup, std::string spectral_flux_file, bool save_output, std::string output_path);
std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s, bool save_output, std::string output_path);
//std::vector<double> axion_electron_counts (double mass, double gaee, std::vector<double> bins, std::vector<double> peak_positions, bool save_output, std::string output_path);
//std::vector<double> axion_electron_counts (double mass, double gaee, exp_setup setup, std::vector<double> peak_positions, bool save_output, std::string output_path);



#endif // defined __experimental_flux_hpp__
