#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"

int main() {
  const std::string sol_agss09met = "data/SolarModel_AGSS09met.dat";
  const std::string sol_agss09met_old = "data/SolarModel_AGSS09met_old.dat";
  const std::string sol_agss09ph = "data/SolarModel_AGSS09met.dat";
  const std::string sol_gs98 = "data/SolarModel_GS98.dat";

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "# Computing axion spectra for all Solar models..." << std::endl;

  SolarModel agss09met (sol_agss09met, OP, true);
  SolarModel agss09met_old (sol_agss09met_old, OP, true);
  SolarModel agss09ph (sol_agss09ph, OP, true);
  SolarModel gs98 (sol_gs98, OP, true);

  auto t2 = std::chrono::high_resolution_clock::now();
  auto t12 = std::chrono::duration_cast<std::chrono::minutes>(t2-t1).count();
  std::cout << "# Setting up all Solar models took " << t12 << " minutes." << std::endl;

  const int n_test_values = 1000;
  std::vector<double> test_ergs;
  for (int k = 0; k < n_test_values; k++) { test_ergs.push_back(0.1+11.9/n_test_values*(k)); };

  std::cout << "# Computing Primakoff spectra for all models..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, agss09met, "gagg_agss09met");
  calculate_spectral_flux_Primakoff(test_ergs, agss09met_old, "gagg_agss09met_old");
  calculate_spectral_flux_Primakoff(test_ergs, agss09ph, "gagg_agss09ph");
  calculate_spectral_flux_Primakoff(test_ergs, gs98, "gagg_gs98");
  auto t3 = std::chrono::high_resolution_clock::now();
  auto t23 = std::chrono::duration_cast<std::chrono::seconds>(t3-t2).count();
  std::cout << "# Computing Primakoff spectra took " << t23 << " seconds." << std::endl;

  std::cout << "# Computing axion-electron spectra for all models..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, agss09met, "gaee_agss09met");
  calculate_spectral_flux_axionelectron(test_ergs, agss09met_old, "gaee_agss09met_old");
  calculate_spectral_flux_axionelectron(test_ergs, agss09ph, "gaee_agss09ph");
  calculate_spectral_flux_axionelectron(test_ergs, gs98, "gaee_gs98");
  auto t4 = std::chrono::high_resolution_clock::now();
  auto t34 = std::chrono::duration_cast<std::chrono::minutes>(t4-t3).count();
  std::cout << "# Computing axion-electron spectra took " << t34 << " minutes." << std::endl;

  std::cout << "# Finished all calculations!" << std::endl;
  return 0;
}
