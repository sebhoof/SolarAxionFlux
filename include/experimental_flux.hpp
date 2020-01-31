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
//const double ergint_abs_prec = 0.0, ergint_rel_prec = 1.0e-3;
const int ergint_from_file_abs_prec = 0.0;
const double ergint_from_file_rel_prec = 1.0e-3, ergint_from_file_method = 5;

struct exp_setup { int n_bins; double bin_lo; double bin_delta; double r_max; double b_field; double length; std::string eff_exposure_file; };
struct erg_integration_params { double mass; double length; double r_max; OneDInterpolator* eff_exp; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; OneDInterpolator* eff_exp; OneDInterpolator* spectral_flux; };

double erg_integrand_from_file(double erg, void * params);
double erg_integrand(double erg, void * params);

//std::vector<double> axion_photon_counts (double mass, double gagg, std::vector<double> bins, bool save_output, std::string output_path);
std::vector<double> axion_photon_counts (double mass, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s);
std::vector<double> axion_electron_counts (double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_electron_counts_full (double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s);



#endif // defined __experimental_flux_hpp__
