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
#include <gsl/gsl_spline.h>
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
    // Screening parameter kappa^2 (kappa^-1 = Debye-Hueckel radius)
    double kappa_squared(double r);
    // Plasma frequency squared (in keV^2)
    double omega_pl_squared(double r);

    // Production rates for the various axion production channels
    double Gamma_P_ff(double omega, double r, int isotope_index);
    double Gamma_P_ff(double omega, double r, Isotope isotope);
    double Gamma_P_ff(double omega, double r); // sum over all isotopes
    double Gamma_P_ee(double omega, double r);
    double Gamma_P_Compton(double omega, double r);
    double Gamma_P_opacity(double omega, double r, std::string element);
    double Gamma_P_opacity(double omega, double r, Isotope isotope); // overloaded for convenience; opacity only depends on chemical element, not on isotope
    double Gamma_P_opacity(double omega, double r); // sum over all elements
    double Gamma_P_Primakoff(double omega, double r);
    double Gamma_P_all_electron(double omega, double r); // sum over all axion-electron interactions

    // Interpolation routines for the opacity data
    double op_grid_interp_erg(double u, int ite, int jne, std::string element);
    double tops_grid_interp_erg(double erg, float t, float rho);
    double opas_grid_interp_erg(double erg, double r);
    double opacity_table_interpolator_op(double omega, double r, std::string element);
    double opacity_table_interpolator_op2(double omega, double r, std::string element);
    double opacity_table_interpolator_tops(double omega, double r);
    double opacity_table_interpolator_opas(double omega, double r);
    double opacity_element(double omega, double r, std::string element);

    // Interpolation routines for the ionisation tables from the Opacity Project
    double ionisationsqr_grid(int ite, int jne, std::string element);
    double ionisationsqr_element(double r, std::string element);

    // N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
    double opacity_element(double omega, double r, Isotope isotope);
    double opacity(double omega, double r);
    // Set the opacity correction of the Opacity Project values according to opacity*(1 + delta), with
    // delta = a + b * log10(T(0)/T(r)) / log10(T(0)/T(r_CZ)), where r_CZ = location of convective zone
    void set_opacity_correction(double a, double b);
    std::vector<double> get_opacity_correction();
    double apply_opacity_correction_factor(double r);

    // Metadata and information from the Solar model
    double get_r_lo();
    double get_r_hi();
    std::vector<double> get_supported_radii(std::vector<double> radii);
    double get_gagg_ref_value_in_inverse_GeV();
    double get_gaee_ref_value();
    std::string get_solaxlib_name_and_version();
    std::string get_solar_model_name();
    std::string get_opacitycode_name();
    bool is_initialised();

  private:
    const opacitycode opcode;
    bool initialisation_status = false;
    // Use the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf) eq. 16 a, default is false
    bool raffelt_approx;
    // Solar model file name (derived from path)
    std::string solar_model_name;
    // Min., max. radius of the solar model file (distance r from the centre of the Sun; in units of the Solar radius)
    double r_lo, r_hi;
    // Opacity corrections; default value = 0
    double opacity_correction_a = 0.0;
    double opacity_correction_b = 0.0;
    ASCIItableReader data;
    ASCIItableReader data_alpha;
    std::map<Isotope,int> isotope_index_map;
    int num_tracked_isotopes;
    int num_interp_pts;
    std::vector<Isotope> tracked_isotopes;
    std::vector<gsl_interp_accel*> accel;
    std::vector<gsl_spline*> linear_interp;
    // TODO. Maybe convert n_isotope and z2_n_isotope into maps like n_element with std::map<Isotope, ...> etc.
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
const std::map<std::string, SolarModelMemberFn> map_interaction_name_to_function { {"Primakoff",&SolarModel::Gamma_P_Primakoff}, {"Compton",&SolarModel::Gamma_P_Compton}, {"ee",&SolarModel::Gamma_P_ee}, {"ff",&SolarModel::Gamma_P_ff},
                                                                                  {"opacity",&SolarModel::Gamma_P_opacity}, {"all_electron",&SolarModel::Gamma_P_all_electron} };
SolarModelMemberFn get_SolarModel_function_pointer(std::string interaction_name);

#endif // defined __solar_model_hpp__
