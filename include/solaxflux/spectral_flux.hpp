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

enum SpectrumModes { table, analytical, solar_model, undefined };

class AxionSpectrum {
  public:
    AxionSpectrum();
    void init_table_mode(std::string file, double g1, double g2=g_aee);
    AxionSpectrum(std::string file, double g1, double g2=g_aee); // Allowed arguments depend on file can be 2, 3, or 4 columns
    //AxionSpectrum(SolarModel* sol, SolarModelMemberFn process1=&SolarModel::Gamma_P_Primakoff, SolarModelMemberFn process2=&SolarModel::Gamma_P_all_electron);
    void init_solar_model_mode(SolarModel* sol, SolarModelMemberFn process2=&SolarModel::Gamma_P_all_electron);
    AxionSpectrum(SolarModel* sol, SolarModelMemberFn process2=&SolarModel::Gamma_P_all_electron); // process1 is always &SolarModel::Gamma_P_Primakoff
    void init_analytical_mode(double norm, double g1, double a, double b);
    AxionSpectrum(double norm, double g1, double a, double b); // In this mode, only one coupling allowed
    ~AxionSpectrum();
    void switch_mode(SpectrumModes new_mode);
    SpectrumModes get_class_mode();
    std::tuple<int, double, double> get_table_parameters();
    std::vector<double> get_analytical_parameters();
    //double axion_flux(double erg);
    double axion_flux(double erg, double g1); // Issue warnings for wrong use of all these functions.
    double axion_flux(double erg, double g1, double g2);
    double axion_flux(double erg, double r, double g1, double g2);
    std::vector<double> axion_flux(std::vector<double> ergs, double g1);
    std::vector<double> axion_flux(std::vector<double> ergs, double g1, double g2);
    std::vector<double> axion_flux(std::vector<double> ergs, double r, double g1, double g2);
    std::vector<std::vector<double>> axion_flux(std::vector<double> ergs, std::vector<double> radii, double g1, double g2);
  private:
    SpectrumModes mode = undefined;
    // Variables for 'table' mode
    int table_submode = -1;
    int pts;
    void init_interp(gsl_interp_accel*& acc, gsl_spline*& interp, const double* x, const double* y);
    void init_numbered_interp(const int index, const double* x, const double* y);
    std::vector<gsl_interp_accel*> acc;
    std::vector<gsl_spline*> spline;
    std::vector<std::pair<gsl_interp_accel*,gsl_interp_accel*>> acc_2d;
    std::vector<gsl_spline2d*> spline_2d;
    double default_g1 = 1.0e6*g_agg;
    double default_g2 = g_aee;
    // Variables for 'solar_model' mode
    SolarModel* s;
    SolarModelMemberFn function1 = &SolarModel::Gamma_P_Primakoff;
    SolarModelMemberFn function2;
    bool solar_model_okay = false;
    // Variables for 'analytical' mode
    std::vector<double> analytical_parameters;
    // General parameters
    //double default_r = 1;
};

class AxionMCGenerator1D {
  public:
    AxionMCGenerator1D();
    AxionMCGenerator1D(std::vector<double> ergs, std::vector<double> flux);
    AxionMCGenerator1D(std::string file, bool is_already_inv_cdf_file=false);
    void change_parameters(double erg_min, double erg_max, double erg_delta);
    AxionMCGenerator1D(SolarModel* sol, double g1=1.0e6*g_agg, double g2=g_aee, double r=1);
    AxionMCGenerator1D(double a, double b);
    AxionMCGenerator1D(AxionSpectrum* spectrum, double g1=1.0e6*g_agg, double g2=g_aee, double r=1);
    ~AxionMCGenerator1D();
    void save_inv_cdf_to_file(std::string inv_cdf_file);
    double evaluate_inv_cdf(double x);
    std::vector<double> draw_axion_energies(int n);
    std::vector<double> draw_axion_energies(int n, double g1, double g2);
    double get_norm();
    std::vector<double> get_omega_params();
  private:
    bool simple_mc_generator_ready = false;
    bool analytical_mc_generator_ready = false;
    void init_inv_cdf_interpolator();
    void init_from_spectral_data(std::vector<double> ergs, std::vector<double> flux);
    std::pair<double,double> analytical_parameters;
    void init_with_local_spectrum(double g1, double g2, double r);
    bool full_mc_generator_ready = false;
    double omega_min = 0.1;
    double omega_max = 10.0;
    double omega_delta = 0.001;
    std::vector<double> generate_ergs();
    //double default_g1 = 1.0e6*g_agg;
    //double default_g2 = g_aee;
    double default_r = 1;
    double integrated_norm;
    AxionSpectrum sp;
    std::vector<double> inv_cdf_data_x;
    std::vector<double> inv_cdf_data_erg;
    gsl_interp_accel* inv_cdf_acc;
    gsl_spline* inv_cdf;
};

class AxionMCGenerator2D {
  public:
    AxionMCGenerator2D();
    AxionMCGenerator2D(std::vector<std::vector<double> > data, bool is_already_inv_cdf_data=false);
    AxionMCGenerator2D(std::string file, bool is_already_inv_cdf_file=false);
    AxionMCGenerator2D(SolarModel &sol, std::vector<double> ergs, std::vector<double> rads, double gaee=0, std::string save_fluxes_prefix="");
    ~AxionMCGenerator2D();
    void save_inv_cdf_to_file(std::string inv_cdf_file);
    double evaluate_inv_cdf_rad(double x);
    double evaluate_inv_cdf_erg_given_rad(double x, double rad);
    std::vector<double> draw_axion_radii(int n);
    std::vector<double> draw_axion_energies_given_radii(std::vector<double> radii);
    std::vector<double> draw_axion_energies(int n);
  private:
    bool mc_generator_ready = false;
    std::vector<double> radii;
    OneDInterpolator inv_cdf_rad;
    std::vector<OneDInterpolator> inv_cdf_ergs;
    //TwoDInterpolator inv_cdf_grid;
    void init_inv_cdf_interpolator();
    void init_from_spectral_data(std::vector<std::vector<double>> data);
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
