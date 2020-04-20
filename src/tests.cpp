#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

int main() {
  auto t_start = std::chrono::high_resolution_clock::now();
  std::cout << "\n     ### This is the " << LIBRARY_NAME << " library ###\n" << std::endl;
  std::cout << "# Testing the Solar Model routines..." << std::endl;

  std::string solar_model_name = "data/SolarModel_B16-AGSS09.dat";
  SolarModel s (solar_model_name,OP,true);
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Setting up the Solar model '" << solar_model_name << "' took "
            << std::chrono::duration_cast<std::chrono::seconds>(t1-t_start).count() << " seconds." << std::endl;
  const int n_test_values = 1000;
  std:: vector<double> test_ergs;
  for (int k=0; k<n_test_values; k++) { test_ergs.push_back(0.1+11.9/n_test_values*(k)); };
  ASCIItableReader javis_data("results/2013_redondo_all.dat");
  std::vector<double> javis_ergs = javis_data[0];

  std::cout << "\n# Test isotope- and element-specific functions..." << std::endl;
  std::cout << "Ratio of 3He and total He (3He + 4He) number densities at 0.5 Rsol: " << s.n_iz(0.5, {"He", 3}) / s.n_element(0.5, "He") << " (should be approx. 0.000251)." << std::endl;
  std::cout << "\n# Calculating Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, "results/primakoff.dat");
  std::cout << "\n# Generating MC object from Primakoff spectrum file, drawing 10 random axions energies, and saving inverse cdf data..." << std::endl;
  AxionMCGenerator mc_engine ("results/primakoff.dat");
  std::vector<double> mc_energies = mc_engine.draw_axion_energies(10);
  std::cout << "MC'ed axion energies (in keV):" << std::endl;
  std::cout << "[ "; for (auto erg = mc_energies.begin(); erg != mc_energies.end(); ++erg) { std::cout << *erg << " "; }; std::cout << "]" << std::endl;
  mc_engine.save_inv_cdf_to_file("results/primakoff_inv_cdf.dat");
  std::cout << "\n# Calculating Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(test_ergs, s, "results/compton.dat");
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Calculating FF spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s, "results/all_ff.dat");
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the FF spectrum took "
            << std::chrono::duration_cast<std::chrono::seconds>(t3-t2).count() << " seconds." << std::endl;
  std::cout << "\n# Computing opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux_opacity(test_ergs, s, "results/metals.dat");
  auto t4 = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Computing full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s, "results/all_gaee.dat");
  auto t5 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_test_values << " energy values) took "
            << std::chrono::duration_cast<std::chrono::seconds>(t5-t4).count() << " seconds." << std::endl;
  //std::cout << "\n# Computing counts in CAST2007 experiment from axion-photon interactions (full calculation)..." << std::endl;
  //axion_photon_counts_full(1.0e-3, 1.0e-10, &cast_2007_setup, &s);
  //auto t6 = std::chrono::high_resolution_clock::now();
  //std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::minutes>(t6-t5).count() << " minutes." << std::endl;
  //double lowerg = 0.1;
  //double higherg = 10.0;
  //std::cout << "# Calculating full flux between " << lowerg << " keV and " << higherg << " keV." << std::endl;
  //std::cout << calculate_flux(lowerg,higherg,s,{}) << std::endl;
  std::cout << "\n# Computing reference counts in CAST2007 experiment from all interactions (from files)..." << std::endl;
  std::vector<double> axion_masses = {1.0e-6, 2.0e-6, 5.0e-6, 1.0e-5, 2.0e-5, 5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0};
  axion_reference_counts_from_file(&cast_2007_setup, axion_masses, "results/primakoff.dat", "results/all_gaee.dat", "results/reference_counts_cast2007.dat");
  auto t6 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the reference counts took " << std::chrono::duration_cast<std::chrono::seconds>(t6-t5).count() << " seconds." << std::endl;
  //std::cout << "# Compute counts in CAST2007 experiment from axion-electron interactions..." << std::endl;
  //axion_electron_counts_full(1.0e-3, 1.0e-13, 1.0e-10, &setup, &s);
  //auto t8 = std::chrono::high_resolution_clock::now();
  //std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::minutes>(t8-t7).count() << " minutes." << std::endl;
  //std::cout << "\n# Compute counts in CAST2007 experiment from axion-electron interactions (from spectrum file)..." << std::endl;
  //axion_electron_counts(1.0e-3, 1.0e-13, 1.0e-10, &cast_2007_setup, "results/all_gaee.dat");
  std::cout << "\n# Computing predicted counts in CAST2007 experiment (m_a = 0.004 eV, g_agamma = 10^-10/GeV, g_ae = 0 or 10^-12)..." << std::endl;
  std::vector<double> counts = counts_prediciton_from_file(4.0e-3, 1.0e-10, "results/reference_counts_cast2007.dat", 0);
  //auto t7 = std::chrono::high_resolution_clock::now();
  std::cout << "Counts (g_agamma): " << counts[0]; for (auto c = counts.begin()+1; c != counts.end(); ++c) { std::cout << " | " << *c; }; std::cout << std::endl;
  //std::cout << "# Initial setup and one calculation took " << std::chrono::duration_cast<std::chrono::seconds>(t7-t6).count() << " seconds." << std::endl;
  counts = counts_prediciton_from_file(4.0e-3, 1.0e-10, "results/reference_counts_cast2007.dat", 1.0e-12);
  //auto t8 = std::chrono::high_resolution_clock::now();
  std::cout << "Counts (all): " << counts[0]; for (auto c = counts.begin()+1; c != counts.end(); ++c) { std::cout << " | " << *c; }; std::cout << std::endl;
  //std::cout << "# Setup step does not have to be repeated; one additional calculation took " << std::chrono::duration_cast<std::chrono::seconds>(t8-t7).count() << " seconds." << std::endl;

  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << "\n# Finished testing! Total runtime: " << std::chrono::duration_cast<std::chrono::minutes>(t_end-t_start).count() << " mins." << std::endl;
  return 0;
}
