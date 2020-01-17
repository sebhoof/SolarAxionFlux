#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <vector>

#include <gsl/gsl_integration.h>

#include "utils.hpp"

const int method1 = 5;
const double abs_prec1 = 1.0e-6, rel_prec1 = 1.0e-6;
struct integration_params { double erg; SolarModel* sol; int iz; };
struct integration_params2 {SolarModel* sol; double* integrand(double, void*); int iz; };
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz);
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*));
void calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s,std::string saveas);
void calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas);
void calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s,std::string saveas);
void calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz);
void calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz,std::string saveas);
void calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s,std::string saveas);
void calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s,std::string saveas);
void calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s);
void calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s,std::string saveas);

#endif // defined __spectral_flux_hpp__
