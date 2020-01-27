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
const double ergint_abs_prec = 0.0, ergint_rel_prec = 1.0e-2;
const int gagg_method = 5;

struct exp_setup { int n_bins; double bin_lo; double bin_delta; double r_max; double b_field; double length; std::string eff_exposure_file; };
struct exp_flux_params { double mass; double length; double r_max; double erg; double rad; OneDInterpolator* eff_exp; SolarModel* s; gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; OneDInterpolator* eff_exp; OneDInterpolator* spectral_flux; };

double erg_integrand_gagg_from_file(double erg, void * params);
double erg_integrand_gagg(double erg, void * params);

//std::vector<double> axion_photon_counts (double mass, double gagg, std::vector<double> bins, bool save_output, std::string output_path);
std::vector<double> axion_photon_counts (double mass, double gagg, exp_setup *setup, std::string spectral_flux_file, bool save_output, std::string output_path);
std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s, bool save_output, std::string output_path);
//std::vector<double> axion_electron_counts (double mass, double gaee, std::vector<double> bins, std::vector<double> peak_positions, bool save_output, std::string output_path);
//std::vector<double> axion_electron_counts (double mass, double gaee, exp_setup setup, std::vector<double> peak_positions, bool save_output, std::string output_path);



#endif // defined __experimental_flux_hpp__
