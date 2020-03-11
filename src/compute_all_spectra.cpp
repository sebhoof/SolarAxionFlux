#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"

int main() {
  // Paths of all Solar model files
  const std::vector<std::string> model_files = {"data/SolarModel_GS98.dat", "data/SolarModel_AGS05.dat", "data/SolarModel_AGSS09.dat", "data/SolarModel_AGSS09ph.dat",
        "data/SolarModel_BP98.dat", "data/SolarModel_BP00.dat", "data/SolarModel_BP04.dat", "data/SolarModel_BS05-OP.dat", "data/SolarModel_BS05-AGSOP.dat", "data/SolarModel_B16-GS98.dat",
        "data/SolarModel_B16-AGSS09.dat"};
  // List all opacity codes
  const std::vector<opacitycode> opacity_codes = {OP, OPAS, LEDCOP, ATOMIC};

  // Set up energy values of the spectra, etc.
  const std::vector<std::string> model_names = {"gs98", "ags05", "agss09", "agss09ph", "bp98", "bp00", "bp04", "bs05op", "bs05agsop", "b16gs98", "b16agss09"};
  const int num_models = model_files.size();
  const int num_opacity_codes = opacity_codes.size();
  const int n_test_values = 20000;
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

    calculate_spectral_flux_Primakoff(test_ergs, sol, "results/gagg_"+model_names[i]+".dat");
    auto t3 = std::chrono::high_resolution_clock::now();
    auto t23 = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    std::cout << "# Computing the Primakoff spectrum for model " << model_files[i] << " took " << t23 << " seconds." << std::endl;

    calculate_spectral_flux_axionelectron(test_ergs, sol, "results/gaee_"+model_names[i]+".dat");
    auto t4 = std::chrono::high_resolution_clock::now();
    auto t34 = std::chrono::duration_cast<std::chrono::seconds>(t4-t3).count();
    std::cout << "# Computing axion-electron spectra took " << t34 << " seconds." << std::endl;
  };

  auto t5 = std::chrono::high_resolution_clock::now();
  auto t05 = std::chrono::duration_cast<std::chrono::minutes>(t5-t0).count();
  std::cout << "# Finished all Solar model calculations after " << t05 << " minutes!" << std::endl;


  std::cout << "# Computing axion spectra for all opacity codes..." << std::endl;

  for (int i = 0; i < num_opacity_codes; i++) {
    auto t7 = std::chrono::high_resolution_clock::now();
    SolarModel sol ("data/SolarModel_AGSS09.dat", opacity_codes[i], true);
    auto t8 = std::chrono::high_resolution_clock::now();
    auto t78 = std::chrono::duration_cast<std::chrono::seconds>(t8-t7).count();
    std::cout << "# Setting up model AGSS09 for opacity code " << opacitycode_names[i] << " took " << t78 << " seconds." << std::endl;

    calculate_spectral_flux_axionelectron(test_ergs, sol, "results/gaee_"+opacitycode_names[i]+".dat");
    auto t9 = std::chrono::high_resolution_clock::now();
    auto t89 = std::chrono::duration_cast<std::chrono::seconds>(t9-t8).count();
    std::cout << "# Computing axion-electron spectra took " << t89 << " seconds." << std::endl;
  };

  auto t10 = std::chrono::high_resolution_clock::now();
  auto t510 = std::chrono::duration_cast<std::chrono::minutes>(t10-t5).count();
  std::cout << "# Finished all opacity code calculations after " << t510 << " minutes!" << std::endl;

  return 0;
}
