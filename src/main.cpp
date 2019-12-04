#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"

int main() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Testing the Solar Model routines..." << std::endl;

  std::string solar_model_name = "data/SolarModel_AGSS09ph.dat";
  SolarModel s (solar_model_name);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Setting up the Solar model '" << solar_model_name << "' took "
            << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds." << std::endl;

  std::vector<double> ergs;
  for (double erg = 0.3; erg < 9.2; erg += 0.1) { ergs.push_back(erg); };

  double test_compton = s.Gamma_P_Compton (1.5, 0.1);
  std::cout << "Test Gamma_Compton value: " << test_compton << std::endl;

  double test_ee = s.Gamma_P_ee (1.5, 0.1);
  std::cout << "Test Gamma_bremsstrahlung value: " << test_ee << std::endl;

  double test_element = s.Gamma_P_element (1.5, 0.1, 1);
  std::cout << "Test element value: " << test_element << std::endl;

  double test_ff = s.Gamma_P_ff(1.5, 0.1, 1);
  std::cout << "Test ff full value: " << test_ff << std::endl;

  double test_pri = s.Gamma_P_Primakoff(1.5, 0.1);
  std::cout << "Test Primakoff value: " << test_pri << std::endl;

  std::cout << "Compute Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(ergs, s);

  std::cout << "Compute Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(ergs, s);

  std::cout << "Compute FF spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(ergs, s);

  auto t4 = std::chrono::high_resolution_clock::now();
  std::cout << "Compute full axion-electron spectrum..." << std::endl;
  ASCIItableReader javis_data("results/2013_redondo_all.dat");
  std::vector<double> javis_masses = javis_data[0];
  calculate_spectral_flux_axionelectron(javis_masses, s);
  auto t5 = std::chrono::high_resolution_clock::now();
  std::cout << "Calculating the full axion-electron spectrum for 23,577 mass values took "
            << std::chrono::duration_cast<std::chrono::seconds>(t5-t4).count()/60.0 << " minutes." << std::endl;

  std::cout << "Finished testing!" << std::endl;
  return 0;
}
