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

using namespace std::chrono;

high_resolution_clock::time_point time_now() { return high_resolution_clock::now(); };

// A simple selection of unit tests for the library.
void run_unit_test() {
  auto t_start = time_now();
  std::cout << "\n     ### This is the " LIBRARY_NAME " library ###\n" << std::endl;
  std::cout << "# Testing the Solar Model routines (this should take about 5-10 mins)...\n" << std::endl;

  auto t1s = time_now();
  std::string solar_model_name = SOLAXFLUX_DIR "/data/solar_models/SolarModel_B16-AGSS09.dat";
  std::string output_path = SOLAXFLUX_DIR "/results/";
  SolarModel s;
  try {
    s = SolarModel("WRONG_FILENAME");
  } catch(XFileNotFound& err) {
    std::cout << "# Oops, the wrong filename was used. This will throw an exception:" << std::endl;
    std::cout << err.what() << std::endl;
    std::cout << "# The exception was caught and handled by using the correct filename. The test will continue now..." << std::endl;
    s = SolarModel(solar_model_name, OP, true);
  }
  auto t1e = time_now();
  std::cout << "\n# Setting up the Solar model '" << solar_model_name << "' took " << duration_cast<seconds>(t1e-t1s).count() << " seconds." << std::endl;
  const int n_erg_values = 500;
  const int n_erg_values_LP = 1000;
  const int n_erg_values_Fe57 = 1000;
  std:: vector<double> test_ergs;
  for (int k=0; k<n_erg_values; k++) { test_ergs.push_back(0.1+k*11.9/n_erg_values); }
  std:: vector<double> test_ergs_LP;
  for (int k=0; k<n_erg_values_LP; k++) { test_ergs_LP.push_back((0.001*gsl_pow_int(1.006,k))); }
  std::vector<double> test_ergs_Fe;
  double width = 0.1;
  double step = width / n_erg_values_Fe57;
  for (int k=0; k<n_erg_values_Fe57; k++) { test_ergs_Fe.push_back(14.4 + (k - n_erg_values_Fe57/2) * step); }
  const int n_rad_values = 6;
  std:: vector<double> test_rads;
  for (int k=0; k<n_rad_values; k++) { test_rads.push_back(k*1.0/(n_rad_values-1)); }

  std::cout << "\n# Test isotope- and element-specific functions..." << std::endl;
  std::cout << "Ratio of 3He and total He (3He + 4He) number densities at 0.5 Rsol: " << s.n_iz(0.5, {"He", 3}) / s.n_element(0.5, "He") << " (should be approx. 0.000251)." << std::endl;
  const double wpl = 0.2;
  std::cout << "Radius where plasma frequency is " << wpl << " keV: " << s.r_from_omega_pl(wpl) << " Rsol (should be approx. 0.144 Rsol)." << std::endl;
  double res_width = -gsl_expm1(-sqrt(s.omega_pl_squared(0))/ s.temperature_in_keV(0))*s.opacity(sqrt(s.omega_pl_squared(0)),0);
  std::cout << "Width of LP resonance in the core: " << res_width << " keV (should be approx. 0.00189 keV)." << std::endl;

  auto t2s = time_now();
  std::cout << "\n# Calculating Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, output_path + "primakoff.dat");
  auto t2e = time_now();
  std::cout << "# Calculating the full Primakoff spectrum (" << n_erg_values << " energy values) took " << duration_cast<milliseconds>(t2e-t2s).count()/1000.0 << " seconds." << std::endl;

  auto t3s = time_now();
  std::cout << "\n# Calculating non-resonant transversal plasmon spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_TP, output_path + "TP.dat");
  auto t3e = time_now();
  std::cout << "# Calculating the full TP spectrum (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t3e-t3s).count() << " seconds." << std::endl;

  auto t4s = time_now();
  std::cout << "\n# Calculating non-resonant transversal plasmon spectrum (Rosseland)..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_TP_Rosseland, output_path + "TP_Rosseland.dat");
  auto t4e = time_now();
  std::cout << "# Calculating the full TP spectrum (Rosseland) (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t4e-t4s).count() << " seconds." << std::endl;

  auto t5s = time_now();
  std::cout << "\n# Calculating resonant longitudinal plasmon spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs_LP, s, &SolarModel::Gamma_LP, output_path + "LP.dat");
  auto t5e = time_now();
  std::cout << "# Calculating the LP spectrum (" << n_erg_values_LP << " energy values) took " << duration_cast<seconds>(t5e-t5s).count() << " seconds." << std::endl;

  auto t6s = time_now();
  std::cout << "\n# Calculating resonant longitudinal plasmon spectrum (Rosseland)..." << std::endl;
  calculate_spectral_flux(test_ergs_LP, s, &SolarModel::Gamma_LP_Rosseland, output_path + "LP_Rosseland.dat");
  auto t6e = time_now();
  std::cout << "# Calculating the LP spectrum (Rosseland) (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t6e-t6s).count() << " seconds." << std::endl;

  auto t7s = time_now();
  std::cout << "\n# Calculating Primakoff spectrum for " << n_rad_values << " different radii..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, test_rads, s, output_path + "primakoff_different_radii.dat");
  auto t7e = time_now();
  std::cout << "# Calculating the Primakoff spectrum (" << n_erg_values << " energy values) for " << n_rad_values << " different radii took " << duration_cast<seconds>(t7e-t7s).count() << " seconds." << std::endl;

  auto t8s = time_now();
  std::cout << "\n# Calculating integrated Primakoff flux between [0,50] keV for " << n_rad_values << " different radii..." << std::endl;
  calculate_total_flux_solar_disc_at_fixed_radii(0.0, 50.0, test_rads, s, &SolarModel::Gamma_Primakoff, output_path + "primakoff_integrated_fluxes.dat");
  auto t8e = time_now();
  std::cout << "# Calculating integrated Primakoff flux for " << n_rad_values << " different radii took " << duration_cast<seconds>(t8e-t8s).count() << " seconds." << std::endl;

  auto t9s = time_now();
  std::cout << "\n# Calculating Compton spectrum..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_Compton, output_path + "compton.dat");
  auto t9e = time_now();
  std::cout << "# Calculating the Compton spectrum (" << n_erg_values << " energy values) took " << duration_cast<milliseconds>(t9e-t9s).count()/1000.0 << " seconds." << std::endl;

  auto t10s = time_now();
  std::cout << "\n# Calculating total ff spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s, output_path + "all_ff.dat");
  auto t10e = time_now();
  std::cout << "# Calculating the ff spectrum (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t10e-t10s).count() << " seconds." << std::endl;

  auto t11s = time_now();
  std::cout << "\n# Computing opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux(test_ergs, s, &SolarModel::Gamma_opacity, output_path + "metals.dat");
  auto t11e = time_now();
  std::cout << "# Calculating the opacity spectrum (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t11e-t11s).count() << " seconds." << std::endl;

  auto t12s = time_now();
  std::cout << "\n# Computing full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s, output_path + "all_gaee.dat");
  auto t12e = time_now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_erg_values << " energy values) took " << duration_cast<seconds>(t12e-t12s).count() << " seconds." << std::endl;
   
  auto t13s = time_now();
  std::cout << "\n# Computing spectrum from Fe57 nuclear transition" << std::endl;
  std::vector<double> specFeFlux = calculate_spectral_flux(test_ergs_Fe, s, &SolarModel::Gamma_Fe57, output_path + "Fe57.dat");
  auto t13e = time_now();
  std::cout << "# Calculating the spectrum from Fe57 nuclear transition for (" << n_erg_values_Fe57 << " energy values) took " << duration_cast<seconds>(t13e-t13s).count() << " seconds." << std::endl;
  double integratedflux = 0.0;
  for (std::vector<double>::iterator it=specFeFlux.begin() ; it != specFeFlux.end(); ++it) { integratedflux += *it; }
  integratedflux *= step;
  std::cout << "The integrated Fe57 flux: " << integratedflux << " g_eff^2 cm^-2 s^-1 (should be 5.0565e+23 g_eff^2 cm^-2 s^-1)." << std::endl;

  auto t_end = time_now();
  std::cout << "\n# Finished testing! Total runtime: " << duration_cast<minutes>(t_end-t_start).count() << " mins." << std::endl;
}

#endif // defined __tests_hpp__
