#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <vector>

#include <gsl/gsl_integration.h>

#include "utils.hpp"
//for spectral flux
const double ref_erg_value = 2.0;
const int int_method_1 = 5, int_method_2 = 2, int_space_size = 1e8;
const double int_abs_prec = 0.0, int_rel_prec = 1.0e-2;
const double abs_prec2 = 0.0, rel_prec2 = 1.0e-3;
struct integration_params { double erg; SolarModel* sol; int iz; };
//struct solar_disc_integration_params { double erg; double rad; double r_max; SolarModel* s; double (*integrand)(double, double); gsl_integration_workspace* w1; };
struct solar_disc_integration_params { double erg; double rad; double r_max; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; };

double rho_integrand(double rho, void * params);
double rad_integrand(double rad, void * params);

//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, int iz, double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas);
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, int iz, double r_max, SolarModel &s, double (*integrand)(double, double));
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas);
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (*integrand)(double, double));
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, int iz, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas);
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, int iz, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double));
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas);
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double));
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, int iz, SolarModel &s, double (*integrand)(double, void*), std::string saveas);
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, int iz, SolarModel &s, double (*integrand)(double, void*));
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas);
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*));

std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas);
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas);
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz);
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz, std::string saveas);
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas);
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas);
//for total flux
struct integration_params2 {SolarModel* sol; double (*integrand)(double, void*); int iz; };
double spectral_flux_integrand(double erg, void * params );
double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s,int iz);
#endif // defined __spectral_flux_hpp__
