#ifndef __python_wrapper_hpp__
#define __python_wrapper_hpp__

#include <pybind11/pybind11.h>

#include "constants.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

void test_module();
void py11_save_Primakoff_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root, std::string process = "Primakoff");
void py11_save_spectral_flux_for_varied_opacities(double erg_min, double erg_max, int n_ergs, double a, double b, std::string solar_model_file, std::string output_file_root);

#endif // defined __spectral_flux_hpp__
