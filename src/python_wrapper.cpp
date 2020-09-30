// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "python_wrapper.hpp"

void module_info() {
  std::cout << "This is the " << LIBRARY_NAME << " Python interface." << std::endl;
}

void test_module() {
  std::cout << "Test successful!" << std::endl;
}

void py11_save_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root, std::string process) {
  SolarModel s (solar_model_file, OP);

  std::vector<double> ergs;
  double erg_stepsize = (erg_max - erg_min)/double(n_ergs);
  for (int i = 0; i < n_ergs+1; i++) { ergs.push_back(erg_min + i*erg_stepsize); };

  if ( (rad_min < s.get_r_lo()) || (rad_max > s.get_r_hi()) ) {
    std::cout << "WARNING! Changing min. and/or max. radius to min./max. radius available in the selected Solar model." << std::endl;
  };
  // Check min/max and avoid Python roundoff errors
  rad_min = std::min(std::max(rad_min, 1.000001*s.get_r_lo()), 0.999999*s.get_r_hi());
  rad_max = std::max(std::min(rad_max, 0.999999*s.get_r_hi()), 1.000001*s.get_r_lo());
  std::cout << "Energies: " <<  n_ergs+1 << " values in [" << erg_min << ", " << erg_max << "] keV, " << n_rads << " radii [" << rad_min << ", " << rad_max << ") R_sol."<< std::endl;
  double rad_stepsize = (rad_max - rad_min)/double(n_rads);

  // Do the calculation for each radius
  for (int j = 0; j < n_rads; j++) {
    double r = rad_min + j*rad_stepsize;
    std::string output_file = output_file_root+"_"+std::to_string(j)+".dat";
    //std::ofstream output;
    //output.open(output_file);
    //output << "# Solar axion spectrum using solar model '" << solar_model_file << "' and r_max = " << r << "." << std::endl;
    //output.close();
    if (process == "Primakoff") {
      calculate_spectral_flux_Primakoff(ergs, s, r, output_file);
    } else if (process == "electron") {
      calculate_spectral_flux_axionelectron(ergs, s, r, output_file);
    } else {
      terminate_with_error("ERROR! the process '"+process+"' is not a valid option. Choose 'Primakoff' or 'electron'.");
    };
  };
}

void py11_save_spectral_flux_for_varied_opacities(double erg_min, double erg_max, int n_ergs, double a, double b, std::string solar_model_file, std::string output_file_root) {
  std::cout << "Setting up Solar model from file " << solar_model_file << " and parameters (" << a << ", " << b << ") ..." << std::endl;
  SolarModel s (solar_model_file, OP);
  s.set_opacity_correction(a, b);

  std::vector<double> ergs;
  double erg_stepsize = (erg_max - erg_min)/double(n_ergs);
  for (int i = 0; i < n_ergs+1; i++) { ergs.push_back(erg_min + i*erg_stepsize); };

  std::string output_file = output_file_root+"_Primakoff.dat";
  calculate_spectral_flux_Primakoff(ergs, s, output_file);
  output_file = output_file_root+"_electron.dat";
  calculate_spectral_flux_axionelectron(ergs, s, output_file);
}

void py11_save_reference_counts(std::vector<double> masses, std::string dataset, std::string ref_spectrum_file_gagg, std::string ref_spectrum_file_gaee, std::string output_file_name) {
  exp_setup setup;
  if (dataset == "CAST2007") {
    setup = cast_2007_setup;
  } else if (dataset == "IAXO") {
    setup = iaxo_setup;
  } else {
    setup = cast_2017_setup;
    setup.dataset = dataset;
  };
  axion_reference_counts_from_file(&setup, masses, ref_spectrum_file_gagg, ref_spectrum_file_gaee, output_file_name, true);
}

std::vector<double> py11_interpolate_saved_reference_counts(double mass, double gagg, std::string reference_counts_file, double gaee) {
  return counts_prediciton_from_file(mass, gagg, reference_counts_file, gaee);
}

void py11_calculate_inverse_cdfs_from_solar_model (std::string solar_model_file, std::vector<double> radii, std::vector<double> energies, double gaee, std::string save_output_prefix) {
  SolarModel s (solar_model_file, OP);
  std::cout << s.temperature_in_keV(0.9) << std::endl;
  std::cout << s.temperature_in_keV(0.95) << std::endl;
  std::cout << s.temperature_in_keV(0.985) << std::endl;
  std::cout << s.kappa_squared(0.985) << std::endl;
  std::cout << s.omega_pl_squared(0.985) << std::endl;
  AxionMCGenerator2D mc (s, energies, radii, gaee, save_output_prefix);
  mc.save_inv_cdf_to_file(save_output_prefix+"_inverseCDFdata");
}

std::vector<std::vector<double> > py11_draw_mc_samples_from_file (std::string mc_file_prefix, int n) {
  AxionMCGenerator2D mc (mc_file_prefix, true);
  std::vector<double> radii = mc.draw_axion_radii(n);
  std::vector<double> energies = mc.draw_axion_energies_given_radii(radii);
  std::vector<std::vector<double> > result = {radii, energies};
  return result;
}

void py11_save_solar_model(std::string solar_model_file, std::string out_file) {
  SolarModel s (solar_model_file, OP);
  double r_lo = s.get_r_lo();
  double r_hi = s.get_r_hi();
  double delta_r = (r_hi - r_lo) / 100.0;
  std::vector<double> radii, temperature, omega_pl, kappa_squared;
  for (int i = 0; i <= 100; i++) {
    double r = r_lo + i*delta_r;
    radii.push_back(r);
    temperature.push_back( s.temperature_in_keV(r) );
    omega_pl.push_back( sqrt(s.omega_pl_squared(r)) );
    kappa_squared.push_back( sqrt(s.kappa_squared(r)) );
  };

  std::vector<std::vector<double>> buffer = { radii, temperature, omega_pl, kappa_squared };
  save_to_file(out_file, buffer);
}


PYBIND11_MODULE(pyaxionflux, m) {
    m.doc() = "Solar Axion Flux library functions";
    m.def("test_module", &test_module, "A simple routine to test the Python interface.");
    m.def("calculate_spectrum", &py11_save_spectral_flux_for_different_radii, "Integrates 'Primakoff' or 'electron' flux from Solar model file with signature (erg_min, erg_max, n_ergs, rad_min, rad_max, n_rads, solar_model_file, output_file_root, process = 'Primakoff').");
    m.def("calculate_spectra_with_varied_opacities", &py11_save_spectral_flux_for_varied_opacities, "Integrates fluxes from Solar model file with signature (erg_min, erg_max, n_ergs, a, b, solar_model_file, output_file_root).");
    m.def("calculate_reference_counts", &py11_save_reference_counts, "Calculate reference counts for experiment with signature (masses, dataset, ref_spectrum_file_gagg, ref_spectrum_file_gaee, output_file_name).");
    m.def("interpolate_reference_counts", &py11_interpolate_saved_reference_counts, "Interpolate reference counts from file with signature (mass, gagg, reference_counts_file, gaee=0).");
    m.def("calculate_inverse_cdfs_from_solar_model", &py11_calculate_inverse_cdfs_from_solar_model, "Calculate (solar_model_file, radii, energies, gaee, save_output_prefix).");
    m.def("draw_mc_samples_from_file", &py11_draw_mc_samples_from_file, ".");
    m.def("save_solar_model", &py11_save_solar_model, ".");
}
