// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "python_wrapper.hpp"

using namespace pybind11::literals;

PYBIND11_MODULE(pyaxionflux, m) {
  m.doc() = "Python wrapper for " LIBRARY_NAME " library functions";

  m.def("module_info", &module_info, "Basic information about the library.");
  m.def("test_module", &test_module, "A few simple unit tests of the library.");
  m.def("save_solar_model", &py11_save_solar_model, "Export the relevant information from a solar model.", "ergs"_a, "solar_model_file"_a, "output_file_root"_a, "n_radii"_a=1000);
  m.def("calculate_spectra", &py11_save_spectral_flux_for_different_radii, "Integrates 'Primakoff' and/or 'ABC' flux from solar model file for different radii.",  "ergs"_a, "radii"_a, "solar_model_file"_a, "output_file_root"_a, "process"_a="Primakoff", "op_code"_a="OP");
  m.def("calculate_spectra_with_varied_opacities", &py11_save_spectral_flux_for_varied_opacities, "Integrates fluxes from solar model file and varies the opacities.", "ergs"_a, "solar_model_file"_a, "output_file_root"_a, "a"_a=0, "b"_a=0);
  m.def("calculate_reference_counts", &py11_calculate_reference_counts, "Calculate reference counts for each bin of a known experiment.", "masses"_a, "dataset"_a, "spectrum_file_P"_a, "spectrum_file_ABC"_a="", "output_file_name"_a="");
  // EXPERIMENTAL
  m.def("interpolate_reference_counts", &py11_interpolate_saved_reference_counts, "EXPERIMENTAL. Interpolate reference counts from file with signature (mass, gagg, reference_counts_file, gaee=0).");
  m.def("calculate_inverse_cdfs_from_solar_model", &py11_calculate_inverse_cdfs_from_solar_model, "EXPERIMENTAL. Calculate the inverse CDF for the axion spectrum from a solar model with signature (solar_model_file, radii, energies, gaee, save_output_prefix).");
  m.def("draw_mc_samples_from_file", &py11_draw_mc_samples_from_file, "EXPERIMENTAL. Draw MC samples from a tabulated spectrum.");
}

void module_info() {
  std::cout << "This is the " << LIBRARY_NAME << " Python interface. The source code was installed in " << SOLAXFLUX_DIR << "." << std::endl;
}

void test_module() {
  run_unit_test();
  std::cout << "INFO. "<< LIBRARY_NAME << " unit test was successful!" << std::endl;
}

void py11_save_solar_model(std::vector<double> ergs, std::string solar_model_file, std::string output_file_root, int n_radii) {
  SolarModel s (solar_model_file, OP);
  double r_lo = s.get_r_lo();
  double r_hi = s.get_r_hi();
  double delta_r = (r_hi - r_lo) / double(n_radii);
  std::vector<double> radii_1, temperature, n_e, n_bar, omega_pl, kappa_squared;
  std::vector<double> radii_2, opacities;
  for (int i = 0; i <= n_radii; i++) {
    double r = r_lo + i*delta_r;
    radii_1.push_back(r);
    temperature.push_back( s.temperature_in_keV(r) );
    n_e.push_back( s.n_electron(r) );
    n_bar.push_back( s.z2_n(r) );
    omega_pl.push_back( sqrt(s.omega_pl_squared(r)) );
    kappa_squared.push_back( sqrt(s.kappa_squared(r)) );
    for (auto erg = ergs.begin(); erg != ergs.end(); ++erg) {
      radii_2.push_back(r);
      opacities.push_back( s.opacity(*erg, r) );
    }
  }

  std::string comment_1 = "Axion quantities from solar model "+s.get_solar_model_name()+" calculated by " LIBRARY_NAME\
                          ".\nColumns: Radius r [R_sol] | Temperature T [keV] | Electron density n_e [cm^-3] | n_bar [cm^-3] | Plasma frequency omega_pl [keV] | Screening scale kappa_s [keV]";
  std::vector<std::vector<double>> buffer_1 = { radii_1, temperature, n_e, n_bar, omega_pl, kappa_squared };
  save_to_file(output_file_root+"_model.dat", buffer_1, comment_1);
  std::string comment_2 = "Opacites for solar model "+s.get_solar_model_name()+" and opacity code "+s.get_opacitycode_name()+" calculated by" LIBRARY_NAME\
                          ".\nColumns: Radius r [R_sol] | Axion energy [keV] | Opacity [keV]";
  std::vector<std::vector<double>> buffer_2 = { radii_2,  opacities };
  save_to_file(output_file_root+"_opacities.dat", buffer_2, comment_2);
}

void py11_save_spectral_flux_for_different_radii(std::vector<double> ergs, std::vector<double> radii, std::string solar_model_file, std::string output_file_root, std::string process, std::string op_code) {
  SolarModel s (solar_model_file, op_code);

  int n_radii = radii.size();
  int n_ergs = ergs.size();
  double rad_min = *std::min_element(radii.begin(), radii.end());
  double rad_max = *std::max_element(radii.begin(), radii.end());
  double erg_min = *std::min_element(ergs.begin(), ergs.end());
  double erg_max = *std::max_element(ergs.begin(), ergs.end());

  if ( (rad_min < s.get_r_lo()) || (rad_max > s.get_r_hi()) ) { std::cout << "WARNING! Changing min. and/or max. radius to min./max. radius available in the selected Solar model." << std::endl; }
  // Check min/max and avoid Python roundoff errors
  rad_min = std::min(std::max(rad_min, 1.000001*s.get_r_lo()), 0.999999*s.get_r_hi());
  rad_max = std::max(std::min(rad_max, 0.999999*s.get_r_hi()), 1.000001*s.get_r_lo());
  std::cout << "INFO. Performing calculation for ";
  if (n_ergs > 1) { std::cout << n_ergs << " energies in [" << erg_min << ", " << erg_max << "] keV and "; } else { std::cout << "one energy value (" << erg_min << " keV) and "; }
  if (n_radii > 1) { std::cout << n_radii << " radii in [" << rad_min << ", " << rad_max << "] R_sol." << std::endl; } else { std::cout << "one radius value (" << rad_min << " R_sol)." << std::endl; }

  // Do the calculation for each radius
  for (int i = 0; i < n_radii; i++) {
    std::string output_file = output_file_root;
    if (n_radii > 1) { output_file += "_"+std::to_string(i); }
    if (process == "Primakoff") {
      calculate_spectral_flux_Primakoff(ergs, s, radii[i], output_file+"_P.dat");
    } else if (process == "ABC") {
      calculate_spectral_flux_axionelectron(ergs, s, radii[i], output_file+"_ABC.dat");
    } else if (process == "both") {
      calculate_spectral_flux_Primakoff(ergs, s, radii[i], output_file+"_P.dat");
      calculate_spectral_flux_axionelectron(ergs, s, radii[i], output_file+"_ABC.dat");
    } else {
      std::string err_msg = "The process '"+process+"' is not a valid option. Choose 'Primakoff', 'ABC', or 'both'.";
      throw XUnsupportedOption(err_msg);
    }
  }
}

void py11_save_spectral_flux_for_varied_opacities(std::vector<double> ergs, std::string solar_model_file, std::string output_file_root, double a, double b) {
  SolarModel s (solar_model_file, OP);
  s.set_opacity_correction(a, b);
  std::vector<double> check = s.get_opacity_correction();
  std::cout << "INFO. Setting up Solar model from file " << solar_model_file << " and opacity correction parameters (" << check[0] << ", " << check[1] << ")." << std::endl;

  std::string output_file = output_file_root+"_Primakoff.dat";
  calculate_spectral_flux_Primakoff(ergs, s, output_file);
  output_file = output_file_root+"_electron.dat";
  calculate_spectral_flux_axionelectron(ergs, s, output_file);
}

std::vector<std::vector<double> > py11_calculate_reference_counts(std::vector<double> masses, std::string dataset, std::string spectrum_file_P, std::string spectrum_file_ABC, std::string output_file_name) {
  exp_setup setup;
  if (dataset == "CAST2007") {
    setup = cast_2007_setup;
  } else if (dataset.compare(0, 9, "CAST2017_") == 0) {
    setup = cast_2017_setup;
    setup.dataset = dataset;
  } else {
    std::string err_msg = "Unkown option for 'dataset' variable. Use ";
    for (auto it = experiment_name.begin(); it != --experiment_name.end(); ++it) { err_msg += it->first + ", "; }; err_msg += "or " + (--experiment_name.end())->first + ".";
    throw XUnsupportedOption(err_msg);
  }

  return axion_reference_counts_from_file(&setup, masses, spectrum_file_P, spectrum_file_ABC, output_file_name, true);
}

// EXPERIMENTAL
std::vector<double> py11_interpolate_saved_reference_counts(double mass, double gagg, std::string reference_counts_file, double gaee) {
  return counts_prediciton_from_file(mass, gagg, reference_counts_file, gaee);
}

void py11_calculate_inverse_cdfs_from_solar_model(std::string solar_model_file, std::vector<double> radii, std::vector<double> energies, double gaee, std::string save_output_prefix) {
  SolarModel s (solar_model_file, OP);
  std::cout << s.temperature_in_keV(0.9) << std::endl;
  std::cout << s.temperature_in_keV(0.95) << std::endl;
  std::cout << s.temperature_in_keV(0.985) << std::endl;
  std::cout << s.kappa_squared(0.985) << std::endl;
  std::cout << s.omega_pl_squared(0.985) << std::endl;
  AxionMCGenerator2D mc (s, energies, radii, gaee, save_output_prefix);
  mc.save_inv_cdf_to_file(save_output_prefix+"_inverseCDFdata");
}

std::vector<std::vector<double> > py11_draw_mc_samples_from_file(std::string mc_file_prefix, int n) {
  AxionMCGenerator2D mc (mc_file_prefix, true);
  std::vector<double> radii = mc.draw_axion_radii(n);
  std::vector<double> energies = mc.draw_axion_energies_given_radii(radii);
  std::vector<std::vector<double> > result = { radii, energies };
  return result;
}
