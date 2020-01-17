#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"

int main() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "# Testing the Solar Model routines..." << std::endl;

  std::string solar_model_name = "data/SolarModel_Javi_red.dat";
  SolarModel s (solar_model_name,OP,true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "# Setting up the Solar model '" << solar_model_name << "' took "
            << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds." << std::endl;
  int n_test_values = 2000;
  std:: vector<double> test_ergs;
  for (int k=0; k<n_test_values;k++ ) {test_ergs.push_back(0.1+11.9/n_test_values*(k));}
  ASCIItableReader javis_data("results/2013_redondo_all.dat");
  std::vector<double> javis_ergs = javis_data[0];
    
  std::cout << "# Compute Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, "primakoff");
  std::cout << "# Compute Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(test_ergs, s,"compton");
  std::cout << "# Compute FF spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s,"all_ff");
  std::cout << "# Compute opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux_opacity(test_ergs, s,"metals");
  auto t4 = std::chrono::high_resolution_clock::now();
  std::cout << "# Compute full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s,"all_gaee");
  auto t5 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_test_values << " energy values) took "
            << std::chrono::duration_cast<std::chrono::seconds>(t5-t4).count() << " seconds." << std::endl;
  std::cout << "# Finished testing!" << std::endl;
  return 0;
}

