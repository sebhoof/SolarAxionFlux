// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "constants.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

// The name of the module and other info
void module_info();
// A simple unit test for the Python wrapper
void test_module();
void py11_save_Primakoff_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root, std::string process = "Primakoff");
void py11_save_spectral_flux_for_varied_opacities(double erg_min, double erg_max, int n_ergs, double a, double b, std::string solar_model_file, std::string output_file_root);
void py11_save_reference_counts(std::vector<double> masses, std::string dataset, std::string ref_spectrum_file_gagg, std::string ref_spectrum_file_gaee, std::string output_file_name);
std::vector<double> py11_interpolate_saved_reference_counts(double mass, double gagg, std::string reference_counts_file, double gaee);
void py11_calculate_inverse_cdfs_from_solar_model(std::string solar_model_file, std::vector<double> radii, std::vector<double> energies, double gaee, std::string save_output_prefix);
std::vector<std::vector<double> > py11_draw_mc_samples_from_file(std::string mc_file_prefix, int n);
void py11_save_solar_model(std::string solar_model_file, std::string out_file);

#endif // defined __spectral_flux_hpp__
