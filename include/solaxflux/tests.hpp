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
  std:: vector<double> test_ergs;
  for (int k=0; k<n_erg_values; k++) { test_ergs.push_back(0.1+k*11.9/n_erg_values); }
  const int n_rad_values = 6;
  std:: vector<double> test_rads;
  for (int k=0; k<n_rad_values; k++) { test_rads.push_back(k*1.0/(n_rad_values-1)); }

  std::cout << "\n# Test isotope- and element-specific functions..." << std::endl;
  std::cout << "Ratio of 3He and total He (3He + 4He) number densities at 0.5 Rsol: " << s.n_iz(0.5, {"He", 3}) / s.n_element(0.5, "He") << " (should be approx. 0.000251)." << std::endl;

  std::cout << "\n# Calculating Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, "results/primakoff.dat");

  std::cout << "\n# Calculating non-resonant transversal plasmon spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_P_TP, "results/TP.dat");

  auto t11s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating Primakoff spectrum for " << n_rad_values << " different radii..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, test_rads, s, "results/primakoff_different_radii.dat");
  auto t11e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the Primakoff spectrum (" << n_erg_values << " energy values) for " << n_rad_values << " different radii took " << std::chrono::duration_cast<std::chrono::seconds>(t11e-t11s).count() << " seconds." << std::endl;

  auto t12s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating integrated Primakoff flux between [0,50] keV for " << n_rad_values << " different radii..." << std::endl;
  calculate_total_flux_solar_disc_at_fixed_radii(0.0, 50.0, test_rads, s, &SolarModel::Gamma_P_Primakoff, "results/primakoff_integrated_fluxes.dat");
  auto t12e = std::chrono::high_resolution_clock::now();
  std::cout << "# alculating integrated Primakoff flux for " << n_rad_values << " different radii took " << std::chrono::duration_cast<std::chrono::seconds>(t12e-t12s).count() << " seconds." << std::endl;

  std::cout << "\n# Calculating Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(test_ergs, s, "results/compton.dat");

  auto t2s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating total ff spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s, "results/all_ff.dat");
  auto t2e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the ff spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t2e-t2s).count() << " seconds." << std::endl;

  auto t3s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Computing opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux_opacity(test_ergs, s, "results/metals.dat");
  auto t3e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the opacity spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t3e-t3s).count() << " seconds." << std::endl;

  auto t4s = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Computing full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s, "results/all_gaee.dat");
  auto t4e = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_erg_values << " energy values) took " << std::chrono::duration_cast<std::chrono::seconds>(t4e-t4s).count() << " seconds." << std::endl;

  // BEGIN EXPERIMENTAL
  //std::cout << "\n# Generating MC object from Primakoff spectrum file, drawing 10 random axions energies, and saving inverse cdf data..." << std::endl;
  //AxionMCGenerator1D mc_engine ("results/primakoff.dat");
  //std::vector<double> mc_energies = mc_engine.draw_axion_energies(10);
  //std::cout << "MC'ed axion energies (in keV):" << std::endl;
  //std::cout << "[ "; for (auto erg = mc_energies.begin(); erg != mc_energies.end(); ++erg) { std::cout << *erg << " "; }; std::cout << "]" << std::endl;
  //mc_engine.save_inv_cdf_to_file("results/primakoff_inv_cdf.dat");
  // END EXPERIMENTAL
  //std::cout << "\n# Computing counts in CAST2007 experiment from axion-photon interactions (full calculation)..." << std::endl;
  //axion_photon_counts_full(1.0e-3, 1.0e-10, &cast_2007_setup, &s);
  //auto t6 = std::chrono::high_resolution_clock::now();
  //std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::minutes>(t6-t5).count() << " minutes." << std::endl;
  //double lowerg = 0.1;
  //double higherg = 10.0;
  //std::cout << "# Calculating full flux between " << lowerg << " keV and " << higherg << " keV." << std::endl;
  //std::cout << calculate_flux(lowerg,higherg,s,{}) << std::endl;

  // BEGIN EXPERIMENTAL (TO REVIEW!)
  // std::cout << "\n# Computing reference counts in CAST2007 experiment from all interactions (from files)..." << std::endl;
  // auto t5s = std::chrono::high_resolution_clock::now();
  // std::vector<double> axion_masses = { 1.0e-6, 2.0e-6, 5.0e-6, 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 };
  // axion_reference_counts_from_file(&cast_2007_setup, axion_masses, "results/primakoff.dat", "results/all_gaee.dat", "results/reference_counts_cast2007.dat");
  // auto t5e = std::chrono::high_resolution_clock::now();
  // std::cout << "# Calculating " << 2*axion_masses.size()*cast_2007_setup.n_bins << " reference counts took " << std::chrono::duration_cast<std::chrono::seconds>(t5e-t5s).count() << " seconds." << std::endl;
  // //std::cout << "# Compute counts in CAST2007 experiment from axion-electron interactions..." << std::endl;
  // //axion_electron_counts_full(1.0e-3, 1.0e-13, 1.0e-10, &setup, &s);
  // //auto t8 = std::chrono::high_resolution_clock::now();
  // //std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::minutes>(t8-t7).count() << " minutes." << std::endl;
  // //std::cout << "\n# Compute counts in CAST2007 experiment from axion-electron interactions (from spectrum file)..." << std::endl;
  // //axion_electron_counts(1.0e-3, 1.0e-13, 1.0e-10, &cast_2007_setup, "results/all_gaee.dat");
  // std::cout << "\n# Computing predicted counts in CAST2007 experiment (m_a = 0.004 eV, g_agamma = 10^-10/GeV, g_ae = 0 or 10^-12)..." << std::endl;
  // std::vector<double> counts = counts_prediciton_from_file(4.0e-3, 1.0e-10, "results/reference_counts_cast2007.dat", 0);
  // //auto t7 = std::chrono::high_resolution_clock::now();
  // std::cout << "Counts (g_agamma): " << counts[0]; for (auto c = counts.begin()+1; c != counts.end(); ++c) { std::cout << " | " << *c; }; std::cout << std::endl;
  // //std::cout << "# Initial setup and one calculation took " << std::chrono::duration_cast<std::chrono::seconds>(t7-t6).count() << " seconds." << std::endl;
  // counts = counts_prediciton_from_file(4.0e-3, 1.0e-10, "results/reference_counts_cast2007.dat", 1.0e-12);
  // //auto t8 = std::chrono::high_resolution_clock::now();
  // std::cout << "Counts (all): " << counts[0]; for (auto c = counts.begin()+1; c != counts.end(); ++c) { std::cout << " | " << *c; }; std::cout << std::endl;
  // //std::cout << "# Setup step does not have to be repeated; one additional calculation took " << std::chrono::duration_cast<std::chrono::seconds>(t8-t7).count() << " seconds." << std::endl;
  // END EXPERIMENTAL (TO REVIEW!)

  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Finished testing! Total runtime: " << std::chrono::duration_cast<std::chrono::minutes>(t_end-t_start).count() << " mins." << std::endl;
}

#endif // defined __tests_hpp__
