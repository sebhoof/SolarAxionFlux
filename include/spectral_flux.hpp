#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <vector>

#include <gsl/gsl_integration.h>

#include "utils.hpp"
//for spectral flux
const int method1 = 5;
const int method2 = 2;
const double abs_prec1 = 0.0, rel_prec1 = 1.0e-3;
const double abs_prec2 = 0.0, rel_prec2 = 1.0e-3;
struct integration_params { double erg; SolarModel* sol; int iz; };
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz);
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*));
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz);
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz,std::string saveas);
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s,std::string saveas);
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s);
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s,std::string saveas);
//for total flux
struct integration_params2 {SolarModel* sol; double (*integrand)(double, void*); int iz; };
double spectral_flux_integrand(double erg, void * params );
double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s,int iz);
#endif // defined __spectral_flux_hpp__
