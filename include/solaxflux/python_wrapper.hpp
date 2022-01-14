// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "constants.hpp"
#include "utils.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"
#include "tests.hpp"

// The name of the module and other info
void module_info();
// A simple unit test for the Python wrapper
void test_module();
// Save the relevant quantities for axion physics as separate text files
void py11_save_solar_model(std::vector<double> ergs, std::string solar_model_file, std::string output_file_root, int n_radii = 1000);
// Calculate the spectral flux for different energies and radii
void py11_save_spectral_flux_for_different_radii(std::vector<double> ergs, std::vector<double> radii, std::string solar_model_file, std::string output_file_root, std::string process = "Primakoff", std::string op_code = "OP");
// Calculate the spectral flux for different energies (radius = 1) with corrected opacities and magnetic fields
void py11_save_varied_spectral_flux(std::vector<double> ergs, std::string solar_model_file, std::string output_file_root, double a = 0.0, double b = 0.0, std::vector<double> c = {0.0, 0.0, 0.0});
// Calculate reference counts for a named helioscope experiment/dataset
std::vector<std::vector<double> > py11_calculate_reference_counts(std::vector<double> masses, std::string dataset, std::string ref_spectrum_file_gagg, std::string ref_spectrum_file_gaee, std::string output_file_name);

#endif // defined __spectral_flux_hpp__
