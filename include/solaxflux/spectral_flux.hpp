#ifndef __spectral_flux_hpp__
#define __spectral_flux_hpp__

#include <iostream>
#include <vector>
#include <random>

#include <gsl/gsl_integration.h>

#include "utils.hpp"
#include "solar_model.hpp"

enum SpectrumModes {table, analytical, solar_model};

class AxionMCGenerator {
  public:
    AxionMCGenerator();
    // Compute custom spectrum in energy range of interest
    AxionMCGenerator(SolarModel s, SolarModelMemberFn process, double omega_min, double omega_max, double omega_delta = 0.001, double r_max = 1.0);
    AxionMCGenerator(std::string spectrum_file);
    ~AxionMCGenerator();
    void init(std::string inv_cdf_file);
    void save_inv_cdf_to_file(std::string inv_cdf_file);
    double evaluate_inv_cdf(double x);
    std::vector<double> draw_axion_energies(int n);
    double get_norm();
  private:
    void init_inv_cdf_interpolator();
    bool mc_generator_ready = false;
    double integrated_norm;
    std::vector<double> inv_cdf_data_x;
    std::vector<double> inv_cdf_data_erg;
    gsl_interp_accel* inv_cdf_acc;
    gsl_spline* inv_cdf;
};

class AxionSpectrum {
  public:
    AxionSpectrum();
    AxionSpectrum(std::string file); // Allowed arguments depend on file can be 2, 3, or 4 columns
    AxionSpectrum(SolarModel sol, SolarModelMemberFn process);
    AxionSpectrum(double norm, double default_g1, double a, double b); // In this mode, only one coupling allowed
    ~AxionSpectrum();
    double axion_flux(double erg, double g1); // Issue warnings for wrong use of all these functions.
    double axion_flux(double erg, double g1, double g2);
    double axion_flux(double erg, double r, double g1, double g2);
    std::vector<double> axion_flux(std::vector<double> ergs, double g1);
    std::vector<double> axion_flux(std::vector<double> ergs, double g1, double g2);
    std::vector<double> axion_flux(std::vector<double> ergs, double r, double g1, double g2);
    std::vector<std::vector<double>> axion_flux(std::vector<double> ergs, std::vector<double> radii, double g1, double g2, bool saveas="");
  private:
    SpectrumModes mode;
    AxionMCGenerator mc;
    double default_g1 = 1.0e6*g_agg;
    double default_g2 = g_aee;
    double default_r = 1;
    double r_lo = -1;
    double r_hi = 0;
};

// Integrate the spectral flux
double spectral_flux_integrand(double erg, void * params);
// Parameter structs for GSL integrators.
struct integration_params { double erg; SolarModel* sol; Isotope isotope; };
//struct solar_disc_integration_params { double erg; double rad; double r_max; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_workspace* w1; };
struct solar_disc_integration_params { double erg; double rad; double r_max; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_cquad_workspace* w1; };
struct integration_params2 { SolarModel* sol; double (*integrand)(double, void*); Isotope isotope; };

// Various overloaded routines to calculate the Solar spectral axion flux.
// TODO: Simplify the structure of these with default values, etc.
// TODO: Possible solution is sth like std::vector<double> calculate_spectral_flux_process(std::vector<double> ergs, std::string process_name, Isotope isotope, std::string saveas);

//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, Isotope isotope, double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas);
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, Isotope isotope, double r_max, SolarModel &s, double (*integrand)(double, double));
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas);
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (*integrand)(double, double));
// -> solar_model.hpp/cpp
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, Isotope isotope, SolarModel &s, double (*integrand)(double, void*), std::string saveas = "");
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas = "");
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas = "");
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas = "");
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_isotope(std::vector<double> ergs, SolarModel &s, Isotope isotope, std::string saveas = "");
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s,std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas= "");
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas = "");
double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s, Isotope isotope);

// For simple integrated flux
double integrated_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file, bool includes_electron_interactions = true);

#endif // defined __spectral_flux_hpp__
