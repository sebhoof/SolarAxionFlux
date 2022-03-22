// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __solar_model_hpp__
#define __solar_model_hpp__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"
#include "utils.hpp"

// SolarModel class: Provides a container to store a (tabulated) Solar model and functions to return its properties.
class SolarModel {
  public:
    // Constructors, destructors, operators...
    SolarModel();
    SolarModel(std::string path_to_model_file, opacitycode opcode_tag = OP, const bool set_raffelt_approx = false);
    SolarModel(std::string path_to_model_file, std::string opcode_name, const bool set_raffelt_approx = false);
    ~SolarModel();
    SolarModel& operator=(SolarModel&&);
    // Delete copy constructor and assignment operator to avoid shallow copies
    SolarModel(const SolarModel&) = delete;
    SolarModel operator=(const SolarModel&) = delete;

    // Solar properties
    // Interpolate numerated data
    double interp_index(int i, double r);
    // Fast lookup function for the standard isotope index
    int lookup_isotope_index(Isotope isotope);
    // Solar plasma temperature in keV
    double temperature_in_keV(double r);
    // Solar plasma density in g cm^-3
    double density(double r);
    // Routines to return ion number density (in cm^-3) times (normalised) charge^2
    double z2_n_iz(double r, int isotope_index);
    double z2_n_iz(double r, Isotope isotope);
    double z2_n(double r); // sum over all isotopes
    // Routine to return ion density (in cm^-3) for each isotope
    double n_iz(double r, int isotope_index);
    double n_iz(double r, Isotope isotope);
    double n_element(double r, std::string el_name);
    double mass_fraction(double r, std::string element);
    // Metallicity Z
    double metallicity(double r);
    // alpha is the expected contribution of all metals to z2_n per nucleon density: z2_n = (X + Y + alpha*Z)*density / m_u
    double alpha(double r);
    // The electron density of the plasma (in cm^-3)
    double n_electron(double r);
    // Screening parameter kappa^2 (in keV^2; kappa^-1 = Debye-Hueckel radius)
    double kappa_squared(double r);
    // Plasma frequency squared (in keV^2)
    double omega_pl_squared(double r);
    double r_from_omega_pl(double omega_pl);

    // Solar B-field
    double bfield(double r);

    // Production rates for the various axion production channels
    double Gamma_ff(double omega, double r, int isotope_index);
    double Gamma_ff(double omega, double r, std::vector<int> isotope_indices);
    double Gamma_ff(double omega, double r, Isotope isotope);
    double Gamma_ff(double omega, double r); // sum over all isotopes
    double Gamma_ee(double omega, double r);
    double Gamma_Compton(double omega, double r);
    double Gamma_opacity(double omega, double r, std::string element);
    double Gamma_opacity(double omega, double r, Isotope isotope); // overloaded for convenience; opacity only depends on chemical element, not on isotope
    double Gamma_opacity(double omega, double r); // sum over all elements
    double Gamma_all_electron(double omega, double r); // sum over all axion-electron interactions
    double Gamma_Primakoff(double omega, double r); // The usual Primakoff rate, but incl. corrections for lower energies < 1 keV
    double Gamma_LP(double omega, double r);
    double Gamma_LP_Rosseland(double omega, double r); // using Rosseland opacities
    double Gamma_TP(double omega, double r); // only non-resonant part (m_a = 0)
    double Gamma_TP_Rosseland(double omega, double r); // using Rosseland opacities; only non-resonant part (m_a = 0)
    double Gamma_plasmon(double omega, double r); // all plasmon interactions
    double Gamma_all_photon(double omega, double r); // sum over all axion-photon interactions
    // General nuclear transition and most improtant iron 57 
    double Gamma_nuclear(double omega, double r, Nucleartransition trans);
    double Gamma_Fe57(double omega, double r);
    // Interpolation routines for the opacity data
    double op_grid_interp_erg(double u, int ite, int jne, std::string element);
    double tops_grid_interp_erg(double erg, float t, float rho);
    double opas_grid_interp_erg(double erg, double r);
    double opacity_table_interpolator_op(double omega, double r, std::string element);
    double opacity_table_interpolator_tops(double omega, double r);
    double opacity_table_interpolator_opas(double omega, double r);
    double opacity_element(double omega, double r, std::string element);

    // Interpolation routines for the ionisation tables from the Opacity Project
    double ionisationsqr_grid(int ite, int jne, std::string element);
    double ionisationsqr_element(double r, std::string element);

    // N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
    double opacity_element(double omega, double r, Isotope isotope);
    double opacity(double omega, double r);
    std::vector<double> log10_rosseland_opacity(std::vector<double> radii);
    double interpolate_rosseland_opacity(double r);

    // Electron degeneracy-related functions
    double calc_electron_chemical_potential(double r);
    double electron_chemical_potential(double r);
    std::vector<std::vector<double> > calc_electron_degeneracy_factor(std::vector<double> ergs, std::vector<double> all_radii);
    std::vector<double> calc_averaged_electron_degeneracy_factor(std::vector<double> radii);
    double avg_degeneracy_factor(double r);

    // B-field correction
    void set_bfields(double b_rad, double b_tach, double b_outer); // Set B-fields in tesla
    std::vector<double> get_bfields();
    // Set the opacity correction of the Opacity Project values according to opacity*(1 + delta), with
    // delta = a + b * log10(T(0)/T(r)) / log10(T(0)/T(r_CZ)), where r_CZ = location of convective zone
    void set_opacity_correction(double a, double b);
    std::vector<double> get_opacity_correction();
    double apply_opacity_correction_factor(double r);

    // Metadata and information from the Solar model
    void save_solar_model_data(std::string output_file_root, std::vector<double> ergs, int n_radii=1000);
    double get_r_lo();
    double get_r_hi();
    std::vector<double> get_supported_radii(std::vector<double> radii);
    std::vector<double> get_all_radii();
    double get_gagg_ref_value_in_inverse_GeV();
    double get_gaee_ref_value();
    std::string get_solaxlib_name_and_version();
    std::string get_solar_model_name();
    std::string get_opacitycode_name();
    bool is_initialised();
    
    double interpolated_integrand(double omega, double r);
    gsl_interp_accel* acc_interp_integrand;
    gsl_spline* interp_integrand;
    
  private:
    // INFO
    const opacitycode opcode;
    bool initialisation_status = false;
    // Use the approximation by Raffelt, PRD 33 (1986) 4, Eq. (16a); default is false
    bool raffelt_approx;
    // Solar model file name (derived from path)
    std::string solar_model_name;
    // PROPERTIES
    // Min., max. radius of the solar model file (distance r from the centre of the Sun; in units of the Solar radius)
    double r_lo, r_hi;
    // Opacity corrections; default value = 0
    double opacity_correction_a = 0.0;
    double opacity_correction_b = 0.0;
    // B-field reference values
    double bfield_rad_T = 3.0e3;
    double bfield_tach_T = 50.0;
    double bfield_outer_T = 4.0;
    // DATA AND INTERPOLATION
    ASCIItableReader data;
    ASCIItableReader data_rosseland_opacity;
    std::map<Isotope,int> isotope_index_map;
    int num_tracked_isotopes;
    int num_interp_pts;
    std::vector<Isotope> tracked_isotopes;
    std::vector<gsl_interp_accel*> accel;
    std::vector<gsl_spline*> linear_interp;
    std::vector<gsl_interp_accel*> n_isotope_acc;
    std::vector<gsl_spline*> n_isotope_lin_interp;
    std::vector<gsl_interp_accel*> z2_n_isotope_acc;
    std::vector<gsl_spline*> z2_n_isotope_lin_interp;
    std::map<std::string, gsl_interp_accel*> n_element_acc;
    std::map<std::string, gsl_spline*> n_element_lin_interp;
    std::map< std::string, std::map<std::pair<int,int>, gsl_interp_accel*> > opacity_acc_op;
    std::map< std::string, std::map<std::pair<int,int>, gsl_spline*> > opacity_lin_interp_op;
    std::map<std::pair<float,float>, gsl_interp_accel*> opacity_acc_tops;
    std::map<std::pair<float,float>, gsl_spline*> opacity_lin_interp_tops;
    std::map<double, gsl_interp_accel*> opacity_acc_opas;
    std::map<double, gsl_spline*> opacity_lin_interp_opas;
    std::vector<std::vector<float>> tops_grid;
    std::vector<float> tops_temperatures;
    std::vector<float> tops_densities;
    std::map< std::pair<int,int>,  std::map<std::string, double> > element_ionisationsqr;
    // private routines to initialise internal interpolators.
    void init_interp(gsl_interp_accel*& acc, gsl_spline*& interp, const double* x, const double* y);
    void init_numbered_interp(const int index, const double* x, const double* y);
};

// Typedef of SolarModel member function as 'SolarModelMemberFn'
typedef double (SolarModel::*SolarModelMemberFn)(double,double);
const std::map<std::string, SolarModelMemberFn> map_interaction_name_to_function {
  {"Primakoff", &SolarModel::Gamma_Primakoff}, {"Compton", &SolarModel::Gamma_Compton}, {"ee", &SolarModel::Gamma_ee},
  {"ff", &SolarModel::Gamma_ff}, {"opacity", &SolarModel::Gamma_opacity}, {"all_electron", &SolarModel::Gamma_all_electron}
};
SolarModelMemberFn get_SolarModel_function_pointer(std::string interaction_name);

std::string standard_header(SolarModel *s);

#endif // defined __solar_model_hpp__
