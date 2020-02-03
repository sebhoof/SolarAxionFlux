#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"

int main() {
  const std::vector<std::string> model_files = {"data/SolarModel_GS98.dat", "data/SolarModel_AGS05.dat", "data/SolarModel_AGSS09met.dat",
    "data/SolarModel_AGSS09met_old.dat", "data/SolarModel_AGSS09ph.dat", "data/SolarModel_BP00.dat", "data/SolarModel_BP04.dat",
    "data/SolarModel_BS05-OP.dat", "data/SolarModel_BS05-AGSOP.dat"};
  const std::vector<std::string> model_names = {"gs98", "ags05", "agss09met", "agss09met_old", "agss09ph", "bp00", "bp04", "bs05op", "bs05agsop"};
  const int num_models = model_files.size();
  const int n_test_values = 1000;
  std::vector<double> test_ergs;
  for (int k = 0; k < n_test_values; k++) { test_ergs.push_back(0.1+11.9/n_test_values*(k)); };

  auto t0 = std::chrono::high_resolution_clock::now();
  std::cout << "# Computing axion spectra for all Solar models..." << std::endl;

  for (int i = 0; i < num_models; i++) {
    auto t1 = std::chrono::high_resolution_clock::now();
    SolarModel sol (model_files[i], OP, true);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto t12 = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    std::cout << "# Setting up model " << model_files[i] << " took " << t12 << " seconds." << std::endl;

    calculate_spectral_flux_Primakoff(test_ergs, sol, "gagg_"+model_names[i]);
    auto t3 = std::chrono::high_resolution_clock::now();
    auto t23 = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    std::cout << "# Computing the Primakoff spectrum for model " << model_files[i] << " took " << t23 << " seconds." << std::endl;

    calculate_spectral_flux_axionelectron(test_ergs, sol, "gaee_"+model_names[i]);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto t34 = std::chrono::duration_cast<std::chrono::seconds>(t4-t3).count();
    std::cout << "# Computing axion-electron spectra took " << t34 << " seconds." << std::endl;
  };

  auto t5 = std::chrono::high_resolution_clock::now();
  auto t05 = std::chrono::duration_cast<std::chrono::minutes>(t5-t0).count();
  std::cout << "# Finished all calculations after " << t05 << " minutes!" << std::endl;

  return 0;
}
