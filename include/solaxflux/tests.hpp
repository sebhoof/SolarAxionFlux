// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __tests_hpp__
#define __tests_hpp__

#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

// A simple selection of unit tests for the library.
void run_unit_test() {
  auto t_start = std::chrono::high_resolution_clock::now();
  std::cout << "\n     ### This is the " LIBRARY_NAME " library ###\n" << std::endl;
  std::cout << "# Testing the Solar Model routines (this may take 15 mins or longer)..." << std::endl;

  auto t1s = std::chrono::high_resolution_clock::now();
  std::string solar_model_name = SOLAXFLUX_DIR "/data/solar_models/SolarModel_B16-AGSS09.dat";
  SolarModel s;
  try {
    s = SolarModel("WRONG_FILENAME");
  } catch(XFileNotFound& err) {
    std::cout << "# Oops, the wrong filename was used... This will throw an error like the one below:" << std::endl;
    std::cout << err.what() << std::endl;
    std::cout << "# The error was caught and handled by using the correct filename. The test will continue now..." << std::endl;
    s = SolarModel(solar_model_name, OP, true);
  }
  auto t1e = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Setting up the Solar model '" << solar_model_name << "' took " << std::chrono::duration_cast<std::chrono::seconds>(t1e-t1s).count() << " seconds." << std::endl;

  const int n_erg_values = 500;
  const int n_erg_values_LP = 1000;
  std:: vector<double> test_ergs;
  for (int k=0; k<n_erg_values; k++) { test_ergs.push_back(0.1+k*11.9/n_erg_values); } //0.1
  std:: vector<double> test_ergs_LP;
  for (int k=0; k<n_erg_values_LP; k++) {
      test_ergs_LP.push_back((0.001*gsl_pow_int(1.006,k)));
  }
  const int n_rad_values = 6;
  std:: vector<double> test_rads;
  for (int k=0; k<n_rad_values; k++) { test_rads.push_back(k*1.0/(n_rad_values-1)); }

  std::cout << "\n# Test isotope- and element-specific functions..." << std::endl;
  std::cout << "Ratio of 3He and total He (3He + 4He) number densities at 0.5 Rsol: " << s.n_iz(0.5, {"He", 3}) / s.n_element(0.5, "He") << " (should be approx. 0.000251)." << std::endl;

  auto t2s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, "results/primakoff.dat");
  auto t2e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full Primakoff spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::milliseconds>(t11e-t11s).count()/1000.0 << " seconds." << std::endl;

  auto t3s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating non-resonant transversal plasmon spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_P_TP, "results/TP.dat");
  auto t3e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full TP spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t12e-t12s).count() << " seconds." << std::endl;
    
  auto t4s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating non-resonant transversal plasmon spectrum (Rosseland)..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_P_TP_Rosseland, "results/TP_Rosseland.dat");
  auto t4e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full TP spectrum (Rosseland) (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t12e-t12s).count() << " seconds." << std::endl;
    
  auto t5s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating resonant longitudinal plasmon spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs_LP, s, &SolarModel::Gamma_P_LP, "results/LP.dat");
  auto t5e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the LP spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t13e-t13s).count() << " seconds." << std::endl;

  auto t6s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating resonant longitudinal plasmon spectrum (Rosseland)..." << std::endl;
  calculate_spectral_flux(test_ergs_LP, s, &SolarModel::Gamma_P_LP_Rosseland, "results/LP_Rosseland.dat");
  auto t6e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the LP spectrum (Rosseland) (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t131e-t131s).count() << " seconds." << std::endl;

  auto t7s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating Primakoff spectrum for " << n_rad_values << " different radii..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, test_rads, s, "results/primakoff_different_radii.dat");
  auto t7e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the Primakoff spectrum (" << n_erg_values << " energy values) for " << n_rad_values << " different radii took " << std::chrono::duration_cast<std::chrono::seconds>(t14e-t14s).count() << " seconds." << std::endl;

  auto t8s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating integrated Primakoff flux between [0,50] keV for " << n_rad_values << " different radii..." << std::endl;
  calculate_total_flux_solar_disc_at_fixed_radii(0.0, 50.0, test_rads, s, &SolarModel::Gamma_P_Primakoff, "results/primakoff_integrated_fluxes.dat");
  auto t8e = std::chrono::high_resolution_clock::now();
  std::cout << "# alculating integrated Primakoff flux for " << n_rad_values << " different radii took " << std::chrono::duration_cast<std::chrono::seconds>(t15e-t15s).count() << " seconds." << std::endl;

  auto t9s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(test_ergs, s, "results/compton.dat");
  auto t9e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the Compton spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::milliseconds>(t16e-t16s).count()/1000.0 << " seconds." << std::endl;

  auto t10s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating total ff spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s, "results/all_ff.dat");
  auto t10e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the ff spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t2e-t2s).count() << " seconds." << std::endl;

  auto t11s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Computing opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux_opacity(test_ergs, s, "results/metals.dat");
  auto t11e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the opacity spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t3e-t3s).count() << " seconds." << std::endl;

  auto t12s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Computing full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s, "results/all_gaee.dat");
  auto t12e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t4e-t4s).count() << " seconds." << std::endl;

  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Finished testing! Total runtime: " << std::chrono::duration_cast<std::chrono::minutes>(t_end-t_start).count() << " mins." << std::endl;
}

#endif // defined __tests_hpp__
