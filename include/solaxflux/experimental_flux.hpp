// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __experimental_flux_hpp__
#define __experimental_flux_hpp__

#include <string>
#include <vector>
#include <algorithm>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_integration.h>

#include "constants.hpp"
#include "utils.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"

// Data for various existing and future helioscope experiments
enum experiment { CAST2007, CAST2017_A, CAST2017_B, CAST2017_C, CAST2017_D, CAST2017_E, CAST2017_F, CAST2017_G, CAST2017_H, CAST2017_I, CAST2017_J, CAST2017_K, CAST2017_L, IAXO, BABYIAXO, IAXOPLUS };
const std::map<std::string,experiment> experiment_name = { {"CAST2007",CAST2007}, {"CAST2017_A",CAST2017_A}, {"CAST2017_B",CAST2017_B}, {"CAST2017_C",CAST2017_C}, {"CAST2017_D",CAST2017_D}, {"CAST2017_E",CAST2017_E},
                                                           {"CAST2017_F",CAST2017_F}, {"CAST2017_G",CAST2017_G}, {"CAST2017_H",CAST2017_H}, {"CAST2017_I",CAST2017_I}, {"CAST2017_J",CAST2017_J}, {"CAST2017_K",CAST2017_K},
                                                           {"CAST2017_L",CAST2017_L}, {"IAXO", IAXO}, {"babyIAXO", BABYIAXO}, {"IAXOplus", IAXOPLUS} };

struct exp_setup { int n_bins; double bin_lo; double bin_delta; double erg_resolution; double r_max; double b_field; double length; std::string dataset; };
// CAST 2007 results [hep-ex/0702006; 2004 CCD measurement]
exp_setup cast_2007_setup = { 20, 0.8, 0.3, 0, 0.231738, 9.0, 9.26, "CAST2007" };
const std::vector<int> obs_cast2007 { 1, 3, 1, 1, 1, 2, 1, 2, 0, 2, 0, 1, 0, 2, 2, 0, 2, 1, 2, 2 };
const std::vector<double> bkg_cast2007 { 2.286801272, 1.559182673, 2.390746817, 1.559182673, 2.598637835, 1.039455092, 0.727618599, 1.559182673, 1.247346181, 1.455237199, 1.871019235, 0.831564073, 1.663128217, 1.247346181,
                                         1.143400636, 1.663128217, 1.247346181, 1.247346181, 2.286801272, 1.247346181};
// CAST 2017 results [1705.02290; all detectors]
exp_setup cast_2017_setup = { 10, 2.0, 0.5, 0, 1.0, 9.0, 9.26, "" };

// Possible IAXO setups [arXiv:1904.09155]
exp_setup iaxo_setup = { 14, 1.0, 0.5, 0.1, 1.0, 2.5, 20.0, "IAXO" };
exp_setup baby_iaxo_setup = { 0, 0, 0, 0, 1.0, 2.0, 10.0, "babyIAXO" };
exp_setup iaxo_plus_setup = { 0, 0, 0, 0, 1.0, 3.5, 22.0, "IAXOplus" };

// Conversion probability into (massless) axions, (gagg*B*L/2)^2, in a reference magnetic field with B = 9.0 T, L = 9.26 m, and gagg = 10^-10 GeV^-1.
const double conversion_prob_factor = gsl_pow_2(0.5*(1.0e-3*g_agg)*(9.0/eV2T)*(9.26/eVm));
// Conversion probability correction for massive axions.
double conversion_prob_correction(double mass, double erg, double length);

// Constant numbers for precision etc.
const int ergint_from_file_abs_prec = 0.0;
const double ergint_from_file_rel_prec = 1.0e-3, ergint_from_file_method = 5;

// Data structs to pass variables to functions and integrators.
//struct erg_integration_params { double mass; double length; double r_max; std::string dataset; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct erg_integration_params { double mass; double length; double r_max; std::string dataset; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_cquad_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; std::string dataset; OneDInterpolator* spectral_flux; double support [2]; double sigma; };
struct simple_convolution_params { double sigma; double erg0; OneDInterpolator* spectral_flux; };
struct convolution_params { double erg0; exp_flux_params_file* p; };

// Wrapper functions for integrating the axion spectra.
double erg_integrand_from_file(double erg, void * params);
double erg_integrand(double erg, void * params);

// double convolution_kernel(double erg, void * params);

// Effective exposures and setups (in seconds x cm) for various experiments.
double eff_exposure(double erg, std::string dataset);
double eff_exposure(double erg, experiment dataset);

// Functions to calculate the spectrum with finite energy resolution convolution kernel.
std::vector<double> convolved_spectrum_from_file(std::vector<double> ergs, double support[2], double resolution, std::string filename);

// Functions to calculate the counts in all bins of a helioscope experiment, given an experimental configuration.
std::vector<double> axion_photon_counts_from_file(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_photon_counts_full(double mass, double gagg, exp_setup *setup, SolarModel *s);
std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_electron_counts_full(double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s);
std::vector<std::vector<double>> axion_reference_counts_from_file(exp_setup *setup, std::vector<double> masses, std::string spectral_flux_file_gagg, std::string spectral_flux_file_gaee = "", std::string saveas = "", bool save_convolved_spectra=false);

std::vector<double> counts_prediciton_from_file(double mass, double gagg, std::string reference_counts_file, double gaee = 0);

#endif // defined __experimental_flux_hpp__
