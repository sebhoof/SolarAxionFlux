// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>

#include "utils.hpp"
#include "solar_model.hpp"

const double distance_factor = pow(radius_sol/(1.0e-2*keV2cm),3) / ( pow(1.0e2*distance_sol,2) * (1.0e6*hbar) ); // = Rsol^3 [in keV^-3] / (2 pi^2 d^2 [in cm^2] * 1 keV*s

/////////////////////////////////////////////////////
//  Integration routines for the solar axion flux  //
/////////////////////////////////////////////////////

// Variables and wrapper functions for solar model integration routines
// Variables to define the behaviour of the GSL integrators.
// Integration over the full Sun (1D)
const int int_method_1d = 5, int_space_size_1d = 1e6;
const double int_abs_prec_1d = 0.0, int_rel_prec_1d = 1.0e-3;
// 1D-Integration for files
const int int_method_file = 5, int_space_size_file = 1e6;
const double int_abs_prec_file = 0.0, int_rel_prec_file = 1.0e-4;
// Integration over the central Solar disc (2D)
const int int_method_2d = 5, int_space_size_2d = 1e6, int_space_size_2d_cquad = 1e6;
const double int_abs_prec_2d = 0.0, int_rel_prec_2d = 1.0e-3;
// Parameter structs for GSL integrators.
struct solar_model_integration_parameters_1d { double erg; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_function* f; gsl_integration_workspace* w; };
struct solar_model_integration_parameters_2d { double erg; double rad; double r_1; double r_2; SolarModel* s; double (SolarModel::*integrand)(double, double);
                                               gsl_function* f1; gsl_integration_cquad_workspace* w1; gsl_function* f2; gsl_integration_cquad_workspace* w2; };
struct solar_model_integration_params_custom { double erg; SolarModel* sol; Isotope isotope; };
// Function wrappers for GSL integration over various Solar geometries.
// Integration over the full Sun (1D)
double rho_integrand_1d(double rho, void * params);
double erg_integrand_1d(double erg, void * params);
// Integration over the central Solar disc (2D)
double rho_integrand_2d(double rho, void * params);
double rad_integrand_2d(double rad, void * params);

// General functions for various integration routines
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="", Isotope isotope={});
std::vector<std::vector<double> > calculate_total_flux_solar_disc_at_fixed_radii(double erg_lo, double erg_hi, std::vector<double> radii, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="");
std::vector<std::vector<double> > calculate_spectral_flux_solar_disc_at_fixed_radii(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="", Isotope isotope={});
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="", Isotope isotope={});

// Convenience functions for integrating specific processes
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > calculate_spectral_flux_Primakoff(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas = "");
std::vector<double> calculate_spectral_flux_plasmon(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > calculate_spectral_flux_plasmon(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_plasmon(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionphoton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > calculate_spectral_flux_axionphoton(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionphoton(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas = "");
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > calculate_spectral_flux_axionelectron(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas= "");
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas = "");

// General functions to allow for custom integration routines of non-SolarModel-type functions
std::vector<double> calculate_spectral_flux_custom(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas = "", Isotope isotope={});

// Convienience functions for integrating specific custom processes
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_ee(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_opacity_element(std::vector<double> ergs, SolarModel &s, std::string element, std::string saveas = "");

// Function to perform simple integrations from a text file
double integrated_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file, bool includes_electron_interactions = true);

#endif // defined __spectral_flux_hpp__
