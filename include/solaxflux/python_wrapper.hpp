#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <pybind11/pybind11.h>

#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

void test_module();
void py11_save_Primakoff_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root, std::string process = "Primakoff");

#endif // defined __spectral_flux_hpp__
