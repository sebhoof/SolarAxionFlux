#ifndef __solar_model_hpp__
#define __solar_model_hpp__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"
#include "utils.hpp"

// SolarModel class: Provides a container to store a (tabulated) Solar model and functions to return its properties.
class SolarModel {
  public:
    SolarModel();
    SolarModel(std::string file, opacitycode set_opcode, const bool set_raffelt_approx = false);
    ~SolarModel();
    SolarModel& operator=(SolarModel&&);
    // Delete copy constructor and assignment operator to avoid shallow copies
    SolarModel(const SolarModel&) = delete;
    SolarModel operator=(const SolarModel&) = delete;
    const opacitycode opcode;

    // Solar properties
    int lookup_isotope_index(Isotope isotope);
    // Routine to return the screening parameter kappa^2 (kappa^-1 = Debye-Hueckel radius).
    double kappa_squared(double r);
    // Routine to return the temperature in keV.
    double temperature_in_keV(double r);
    // Routine to return density
    double density(double r);
    // Routines to return ion number density times charge^2 (for each isotope and for all)
    double z2_n_iz(double r, int isotope_index);
    double z2_n_iz(double r, Isotope isotope);
    double z2_n(double r);
    // Routine to return ion density for each isotope
    double n_iz(double r, int isotope_index);
    double n_iz(double r, Isotope isotope);
    double n_element(double r, std::string el_name);
    double H_mass_fraction(double r);
    double He_mass_fraction(double r);
    double metallicity(double r);
    // alpha is the expected contribution of all metals to z2_n per nucleon density: z2_n = (X + Y + \alpha Z)*\rho / m_u
    double alpha(double r);
    // Routine to return electron density
    double n_electron(double r);
    // Routine to return the plasma frequency squared.
    double omega_pl_squared(double r);

    // Rates for various axion production channels
    double Gamma_P_ff(double omega, double r, int isotope_index);
    double Gamma_P_ff(double omega, double r, Isotope isotope);
    double Gamma_P_ff(double omega, double r);
    double Gamma_P_ee(double omega, double r);
    double Gamma_P_Compton(double omega, double r);
    double Gamma_P_opacity(double omega, double r, std::string element);
    double Gamma_P_opacity(double omega, double r);
    double Gamma_P_Primakoff(double omega, double r);
    double Gamma_P_all_electron(double omega, double r);

    // Calculate the solar axion spectrum for axion-photon and axion-electron interactions
    std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, double r_max=1.0);
    std::vector<double> calculate_spectral_flux_all_electron(std::vector<double> ergs, double r_max=1.0);

    // Interpolation routines for the opacity data
    double op_grid_interp_erg(double u, int ite, int jne, std::string element);
    double tops_grid_interp_erg(double erg, float t, float rho);
    double opas_grid_interp_erg(double erg, double r);
    double opacity_table_interpolator_op(double omega, double r, std::string element);
    double opacity_table_interpolator_op2(double omega, double r, std::string element);
    double opacity_table_interpolator_tops(double omega, double r);
    double opacity_table_interpolator_opas(double omega, double r);
    double opacity_element(double omega, double r, std::string element);
    // N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
    double opacity_element(double omega, double r, Isotope isotope);
    double opacity(double omega, double r);
    // Metadata and information
    double get_r_lo();
    double get_r_hi();
    double get_gagg_ref_value_in_inverse_GeV();
    double get_gaee_ref_value();
    std::string get_solaxlib_name_and_version();
    bool is_initialised();
  private:
    void init_interp(gsl_interp_accel*& acc, gsl_spline*& interp, const double* x, const double* y);
    void init_numbered_interp(const int index, const double* x, const double* y);
    bool initialisation_status = false;
    // Use the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf) eq. 16 a, default is false
    bool raffelt_approx;
    // Solar data
    std::string solarmodel_name;
    // Min. and max. radius of the solar model file (distance r from the centre of the Sun in units of the solar radius)
    double r_lo, r_hi;
    ASCIItableReader data;
    ASCIItableReader data_alpha;
    std::map<Isotope,int> isotope_index_map;
    int num_tracked_isotopes;
    int num_interp_pts;
    std::vector<Isotope> tracked_isotopes;
    std::vector<gsl_interp_accel*> accel;
    std::vector<gsl_spline*> linear_interp;
    std::map< std::string, std::map<std::pair<int,int>, gsl_interp_accel*> > opacity_acc_op;
    std::map< std::string, std::map<std::pair<int,int>, gsl_spline*> > opacity_lin_interp_op;
    std::map<std::pair<float,float>, gsl_interp_accel*> opacity_acc_tops;
    std::map<std::pair<float,float>, gsl_spline*> opacity_lin_interp_tops;
    std::map<double, gsl_interp_accel*> opacity_acc_opas;
    std::map<double, gsl_spline*> opacity_lin_interp_opas;
    std::vector<std::vector<float>> tops_grid;
    std::vector<float> tops_temperatures;
    std::vector<float> tops_densities;
    // TODO. Maybe convert n_isotope and z2_n_isotope into maps like n_element with std::map<Isotope, ...> etc.
    std::vector<gsl_interp_accel*> n_isotope_acc;
    std::vector<gsl_spline*> n_isotope_lin_interp;
    std::vector<gsl_interp_accel*> z2_n_isotope_acc;
    std::vector<gsl_spline*> z2_n_isotope_lin_interp;
    std::map<std::string, gsl_interp_accel*> n_element_acc;
    std::map<std::string, gsl_spline*> n_element_lin_interp;
};

// Variables and wrapper functions for solar model integration routines
// Variables to define the behaviour of the GSL integrators.
// Integration over the full Sun (1D)
const int int_method_1d = 5, int_space_size_1d = 1e8;
const double int_abs_prec_1d = 0.0, int_rel_prec_1d = 1.0e-3;
// Integration over the central Solar disc (2D)
const int int_method_2d = 5, int_space_size_2d = 1e6, int_space_size_2d_cquad = 1e6;
const double int_abs_prec_2d = 0.0, int_rel_prec_2d = 1.0e-2;
// Parameter structs for GSL integrators.
struct solar_model_integration_parameters { double erg; double rad; double r_max; SolarModel* s; double (SolarModel::*integrand)(double, double); gsl_integration_cquad_workspace* w1; };
// Function wrappers for GSL integration over various Solar geometries.
// Integration over the full Sun (1D)
double rho_integrand_1d(double rho, void * params);
// Integration over the central Solar disc (2D)
double rho_integrand_1d(double rho, void * params);
double rad_integrand_2d(double rad, void * params);

std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, Isotope isotope, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="");
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="");
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas="", Isotope isotope={});

#endif // defined __solar_model_hpp__
