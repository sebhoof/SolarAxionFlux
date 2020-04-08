#ifndef __experimental_flux_hpp__
#define __experimental_flux_hpp__

#include <string>
#include <vector>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_integration.h>

#include "constants.hpp"
#include "spectral_flux.hpp"

// Conversion probability into (massless) axions, (gagg*B*L/2)^2, in a reference magnetic field with B = 9.0 T, L = 9.26 m, and gagg = 10^-10 GeV^-1.
const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
const double conversion_prob_factor = gsl_pow_2(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm));

// Conversion probability correction for massive axions.
double conversion_prob_correction(double mass, double erg, double length);

// Constant numbers for precision etc.
const int ergint_from_file_abs_prec = 0.0;
const double ergint_from_file_rel_prec = 1.0e-3, ergint_from_file_method = 5;

// Data structs to pass variables to functions and integrators.
struct exp_setup { int n_bins; double bin_lo; double bin_delta; double erg_resolution; double r_max; double b_field; double length; double (*eff_exposure)(double); };
//struct erg_integration_params { double mass; double length; double r_max; double (*eff_exposure)(double); SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct erg_integration_params { double mass; double length; double r_max; double (*eff_exposure)(double); SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_cquad_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; double (*eff_exposure)(double); OneDInterpolator* spectral_flux; };
struct convolution_params { double sigma; double erg0; OneDInterpolator* spectral_flux; };

// Wrapper functions for integrating the axion spectra.
double erg_integrand_from_file(double erg, void * params);
double erg_integrand(double erg, void * params);

//double convolution_kernel(double erg, void * params);

// Effective exposures and setups (in seconds x cm) for various experiments.
// CAST 2007 results [hep-ex/0702006; 2004 CCD measurement]
double eff_exposure_cast_2007 (double erg);
exp_setup cast_2007_setup = { 20, 0.8, 0.3, 0, 0.231738, 9.0, 9.26, &eff_exposure_cast_2007 };

// Functions to calculate the spectrum with finite energy resolution convolution kernel.
std::vector<double> convolved_spectrum_from_file(std::vector<double> ergs, double support[2], double resolution, std::string filename);

// Functions to calculate the counts in all bins of a helioscope experiment, given an experimental configuration.
std::vector<double> axion_photon_counts_from_file(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_photon_counts_full(double mass, double gagg, exp_setup *setup, SolarModel *s);
std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_electron_counts_full(double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s);
std::vector<std::vector<double>> axion_photon_counts_from_file(exp_setup *setup, std::vector<double> masses, std::string spectral_flux_file_gagg, std::string spectral_flux_file_gaee = "", std::string saveas = "");

#endif // defined __experimental_flux_hpp__
