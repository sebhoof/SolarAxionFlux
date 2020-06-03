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

// Conversion probability into (massless) axions, (gagg*B*L/2)^2, in a reference magnetic field with B = 9.0 T, L = 9.26 m, and gagg = 10^-10 GeV^-1.
const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
const double conversion_prob_factor = gsl_pow_2(0.5*(1.0e-3*g_agg)*(9.0/eV2T)*(9.26/eVm));

// Conversion probability correction for massive axions.
double conversion_prob_correction(double mass, double erg, double length);

// Constant numbers for precision etc.
const int ergint_from_file_abs_prec = 0.0;
const double ergint_from_file_rel_prec = 1.0e-3, ergint_from_file_method = 5;

// Data structs to pass variables to functions and integrators.
struct exp_setup { int n_bins; double bin_lo; double bin_delta; double erg_resolution; double r_max; double b_field; double length; std::string dataset; };
//struct erg_integration_params { double mass; double length; double r_max; std::string dataset; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; gsl_integration_workspace* w2; };
struct erg_integration_params { double mass; double length; double r_max; std::string dataset; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_cquad_workspace* w1; gsl_integration_workspace* w2; };
struct exp_flux_params_file { double mass; double length; std::string dataset; OneDInterpolator* spectral_flux; };
struct convolution_params { double sigma; double erg0; OneDInterpolator* spectral_flux; };

// Wrapper functions for integrating the axion spectra.
double erg_integrand_from_file(double erg, void * params);
double erg_integrand(double erg, void * params);

//double convolution_kernel(double erg, void * params);

// Effective exposures and setups (in seconds x cm) for various experiments.
double eff_exposure(double erg, std::string dataset);
// CAST 2007 results [hep-ex/0702006; 2004 CCD measurement]
exp_setup cast_2007_setup = { 20, 0.8, 0.3, 0, 0.231738, 9.0, 9.26, "CAST2007" };
const std::vector<int> obs_cast2007 {1, 3, 1, 1, 1, 2, 1, 2, 0, 2, 0, 1, 0, 2, 2, 0, 2, 1, 2, 2};
const std::vector<double> bkg_cast2007 {2.286801272, 1.559182673, 2.390746817, 1.559182673, 2.598637835, 1.039455092, 0.727618599, 1.559182673, 1.247346181, 1.455237199, 1.871019235, 0.831564073, 1.663128217, 1.247346181, 1.143400636, 1.663128217,
                                        1.247346181, 1.247346181, 2.286801272, 1.247346181};
// CAST 2017 results [1705.02290; all detectors]
exp_setup cast_2017_setup = { 10, 2.0, 0.5, 0, 1.0, 9.0, 9.26, "" };
const std::vector<std::vector<int>> obs_cast2017 { {0, 3, 3, 0, 0, 1, 3, 3, 3, 3}, {5, 5, 5, 3, 3, 0, 5, 2, 2, 1}, {3, 3, 1, 2, 2, 2, 4, 5, 4, 3}, {1, 5, 5, 2, 1, 2, 2, 5, 4, 0}, {2, 3, 2, 2, 2, 1, 0, 2, 1, 1}, {3, 5, 1, 4, 1, 2, 0, 3, 2, 2},
                                                   {3, 4, 4, 5, 1, 2, 3, 2, 3, 2}, {2, 1, 0, 1, 3, 2, 2, 3, 0, 1}, {1, 2, 2, 1, 3, 0, 0, 1, 4, 0}, {2, 1, 3, 1, 1, 0, 1, 1, 5, 5}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 2, 1, 0, 0, 0, 0, 0, 0, 0} };

const std::vector<std::vector<double>> bkg_cast2017 { {0.926256, 1.96148, 1.79803, 1.30766, 1.30766, 1.96148, 2.61531, 2.77877, 2.94223, 2.07045}, {3.68151, 4.86486, 4.99634, 3.55003, 2.49817, 3.28707, 2.89262, 3.68151, 3.48429, 3.41855},
                                                      {2.54573, 3.18216, 4.45502, 2.86394, 2.29116, 2.29116, 3.30945, 3.75495, 3.62766, 3.56402}, {2.72482, 5.5794, 3.95748, 2.40044, 2.27069, 2.33556, 3.37359, 3.43847, 3.24384, 3.11408},
                                                      {1.44613, 2.30066, 2.43213, 1.70906, 1.97199, 1.24893, 1.24893, 2.23493, 2.16919, 2.23493}, {1.30963, 2.94666, 2.35733, 2.55377, 2.02992, 1.50607, 2.16088, 2.75022, 2.29185, 2.29185},
                                                      {2.33334, 2.74167, 2.21667, 2.80001, 2.21667, 1.75001, 2.62501, 2.21667, 2.80001, 2.33334}, {1.74724, 2.37125, 2.68326, 1.62243, 2.05924, 1.74724, 1.49763, 1.74724, 1.18563, 2.24645},
                                                      {1.72998, 3.45995, 1.79405, 1.72998, 1.9222, 1.72998, 2.69107, 2.24256, 1.98627, 2.11442}, {1.89627, 2.25182, 2.96292, 1.4222, 1.65924, 1.65924, 1.95553, 2.1333, 1.71849, 2.07404},
                                                      {0.0150685, 0.0568493, 0.060274, 0.0150685, 0.0150685, 0.00753425, 0.0267123, 0.0150685, 0.0267123, 0.0116438},
                                                      {0.0409574, 0.226904, 0.243287, 0.0532447, 0.0188404, 0.0344043, 0.0417766, 0.0409574, 0.0409574, 0.0286702} };

// Functions to calculate the spectrum with finite energy resolution convolution kernel.
std::vector<double> convolved_spectrum_from_file(std::vector<double> ergs, double support[2], double resolution, std::string filename);

// Functions to calculate the counts in all bins of a helioscope experiment, given an experimental configuration.
std::vector<double> axion_photon_counts_from_file(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_photon_counts_full(double mass, double gagg, exp_setup *setup, SolarModel *s);
std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file);
std::vector<double> axion_electron_counts_full(double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s);
std::vector<std::vector<double>> axion_reference_counts_from_file(exp_setup *setup, std::vector<double> masses, std::string spectral_flux_file_gagg, std::string spectral_flux_file_gaee = "", std::string saveas = "");

std::vector<double> counts_prediciton_from_file(double mass, double gagg, std::string reference_counts_file, double gaee = 0);

#endif // defined __experimental_flux_hpp__
