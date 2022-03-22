// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "tests.hpp"
/*
int main() {
  run_unit_test();
  return 0;
}
 */
#include <iostream>
#include <string>
#include <vector>
#include <random>

#include "solar_model.hpp"
#include "spectral_flux.hpp"
#include "utils.hpp"

int main () {
  std::vector<double> all_radii, all_ergs;
  for (int i = 0; i < 101; ++i) { all_radii.push_back(0.01*i); }
  for (int i = 0; i < 200; ++i) { all_ergs.push_back(0.01*(i+1)); }
  std::cout << "Setting up solar model..." << std::endl;
  SolarModel s (SOLAXFLUX_DIR "/data/solar_models/SolarModel_B16-AGSS09.dat");
  std::cout << "Finished." << std::endl;
//  print_current_time();
//  calculate_d2Phi_a_domega_drho(all_ergs, all_radii, s, &SolarModel::Gamma_Primakoff, "matrix_P.dat");
//  print_current_time();
//  integrate_d2Phi_a_domega_drho_between_rhos(all_ergs, all_radii, s, &SolarModel::Gamma_Primakoff, "matrix_P_rings.dat");
//  print_current_time();
  integrate_d2Phi_a_domega_drho_between_rhos(all_ergs, all_radii, s, &SolarModel::Gamma_Primakoff, "matrix_P_true_rings.dat", true);
  print_current_time();
  // calculate_d2Phi_a_domega_drho(all_ergs, all_radii, s, &SolarModel::Gamma_all_electron, "matrix_gaee.dat");
  // print_current_time();
  //integrate_d2Phi_a_domega_drho_between_rhos(all_ergs, all_radii, s, &SolarModel::Gamma_all_electron, "matrix_gaee_rings.dat");
  // print_current_time();
  integrate_d2Phi_a_domega_drho_between_rhos(all_ergs, all_radii, s, &SolarModel::Gamma_all_electron, "matrix_gaee_true_rings.dat", true);
  print_current_time();

  return 0;
}
