// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <sstream>
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

// Variables that define the behaviour of the GSL integrators and wrapper functions for solar model integration routines
// 1D-Integration for files
const int int_method_file = 5, int_space_size_file = 1e6;
const double int_abs_prec_file = 0.0, int_rel_prec_file = 1.0e-4;
struct solar_model_integration_params_custom { double erg; SolarModel* sol; Isotope isotope; };

// Integration over the full Sun (1D), see Eq. (2.42) in [arXiv:2101.08789]
const int int_method_1d = 5, int_space_size_1d = 1e6;
const double int_abs_prec_1d = 0.0, int_rel_prec_1d = 1.0e-3;
struct solar_model_integration_parameters_1d { double erg; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_function* f; gsl_integration_workspace* w; };
double r_integrand_1d(double r, void * params);
double erg_integrand_1d(double erg, void * params);

// Integration over the central Solar disc (2D), see (2.45) in [arXiv:2101.08789]
const int int_method_2d = 5, int_space_size_2d = 1e6;
const double int_abs_prec_2d = 0.0, int_rel_prec_2d = 1.0e-3;
struct solar_model_integration_parameters_2d { double erg; double rho; double rho_0; double rho_1; SolarModel* s; double (SolarModel::*integrand)(double, double);
  gsl_function* f1; gsl_integration_workspace* w1; gsl_function* f2; gsl_integration_cquad_workspace* w2; };
double r_integrand_2d(double r, void * params);
double z_integrand_2d(double r, void * params);
double rho_integrand_2d(double rho, void * params);
double erg_integrand_2d(double erg, void * params);

// General functions for various integration routines; see Eq. (2.42) and (2.45) in [arXiv:2101.08789]
std::vector<std::vector<double> > calculate_d2Phi_a_domega_drho(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho(std::vector<double> ergs, double rho_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas = "", Isotope isotope = {});
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_between_rhos(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas = "", bool use_ring_geometry=false, Isotope isotope = {});
std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho(std::vector<double> ergs, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas = "", Isotope isotope = {});
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(double erg_lo, double erg_hi, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas = "");

// Convenience functions for integrating specific processes
std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_plasmon(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_plasmon(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_plasmon(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_axionphoton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionphoton(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionphoton(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas = "");
std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas= "");

// General functions to allow for custom integration routines of non-SolarModel-type functions
std::vector<double> calculate_spectral_flux_custom(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas = "", Isotope isotope = {});

// Convienience functions for integrating specific custom processes
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_opacity_element(std::vector<double> ergs, SolarModel &s, std::string element, std::string saveas = "");

// Function to perform simple integrations from a text file
double integrated_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file, bool includes_electron_interactions = true);

#endif // defined __spectral_flux_hpp__
