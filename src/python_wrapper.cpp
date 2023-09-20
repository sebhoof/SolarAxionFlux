// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "python_wrapper.hpp"

using namespace pybind11::literals;

PYBIND11_MODULE(pyaxionflux, m) {
  m.doc() = "Python wrapper for " LIBRARY_NAME " library functions";

  try { pybind11::module_::import("numpy"); } catch (...) { return; }

  m.def("module_info", &module_info, "Basic information about the library.");
  m.def("test_module", &test_module, "A few simple unit tests of the library.");
  pybind11::class_<SolarModel>(m, "SolarModel", "A simplified reduced implementation of the C++ SolarModel class in Python.")
    .def(pybind11::init([](std::string file) { return new SolarModel(file); }), "Class constructor using only the path to the solar model file.", "solar_model_file"_a)
    .def("temperature", pybind11::vectorize(&SolarModel::temperature_in_keV), "Solar model temperature (in keV)", "radius"_a)
    .def("kappa_squared", pybind11::vectorize(&SolarModel::kappa_squared), "Screening scale squared (in keV^2)", "radius"_a)
    .def("omega_pl_squared", pybind11::vectorize(&SolarModel::omega_pl_squared), "Plasma frequency squared (in keV^2)", "radius"_a)
    .def("n_e", pybind11::vectorize(&SolarModel::n_electron), "Electron density (in keV^3)", "radius"_a)
    .def("z2_n", pybind11::vectorize(&SolarModel::z2_n), "Charge-square-weighted ion density (in keV^3)", "radius"_a)
    .def("degeneracy_factor", pybind11::vectorize(&SolarModel::avg_degeneracy_factor), "Electron degeneracy factor", "radius"_a)
    .def("primakoff_rate", pybind11::vectorize(&SolarModel::Gamma_Primakoff), "Primakoff production rate", "omega"_a, "radius"_a)
    .def("abc_rate", pybind11::vectorize(&SolarModel::Gamma_all_electron), "Production rate for ABC processes", "omega"_a, "radius"_a)
    .def("save_solar_model_data", &SolarModel::save_solar_model_data, "Save all solar model data relevant for axion computations.", "output_file_root"_a, "ergs"_a, "n_radii"_a=1000)
  ;
  m.def("calculate_spectra", &py11_calc_spectral_flux_up_to_rmax, "Integrates 'Primakoff' and/or 'ABC' flux from solar model file up to radius rmax.",  "ergs"_a, "rmax"_a, "solar_model"_a, "output_file_root"_a="", "process"_a="Primakoff");
  m.def("save_spectra", &py11_save_spectral_flux_for_different_radii, "Integrates 'Primakoff' and/or 'ABC' flux from solar model file for different radii and saves the results as a text file.",  "ergs"_a, "radii"_a, "solar_model_file"_a, "output_file_root"_a, "process"_a="Primakoff", "op_code"_a="OP");
  const std::vector<double> v1 = { 1.0, 20.0 };
  m.def("calculate_fluxes_on_solar_disc", &py11_calc_integrated_flux_up_to_different_radii, "Integrated flux within different radii on the solar disc.", "radii"_a, "s"_a, "output_file_root"_a="", "erg_limits"_a=v1, "process"_a="Primakoff");
  const std::vector<double> v2 = { 3.0e3, 50.0, 4.0 };
  m.def("calculate_varied_spectra", &py11_save_varied_spectral_flux, "Integrates fluxes from solar model file and varies the opacities.", "ergs"_a, "solar_model_file"_a, "output_file_root"_a, "a"_a=0, "b"_a=0, "c"_a=v2);
  m.def("calculate_reference_counts", &py11_calculate_reference_counts, "Calculate reference counts for each bin of a known experiment.", "masses"_a, "dataset"_a, "spectrum_file_P"_a, "spectrum_file_ABC"_a="", "output_file_name"_a="");
}

void module_info() {
  std::cout << "This is the " << LIBRARY_NAME << " Python interface. The source code was installed in " << SOLAXFLUX_DIR << "." << std::endl;
}

void test_module() {
  run_unit_test();
  std::cout << "INFO. "<< LIBRARY_NAME << " unit test was successful!" << std::endl;
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

  if (process == "Primakoff") {
    integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(ergs, radii, s, output_file_root+"_P.dat");
  } else if (process == "ABC") {
    integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(ergs, radii, s, output_file_root+"_ABC.dat");
  } else if (process == "plasmon") {
    integrate_d2Phi_a_domega_drho_up_to_rho_plasmon(ergs, radii, s, output_file_root+"_plasmon.dat");
  } else if (process == "all") {
      integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(ergs, radii, s, output_file_root+"_P.dat");
      integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(ergs, radii, s, output_file_root+"_ABC.dat");
      integrate_d2Phi_a_domega_drho_up_to_rho_plasmon(ergs, radii, s, output_file_root+"_plasmon.dat");
  } else {
    std::string err_msg = "The process '"+process+"' is not a valid option. Choose 'ABC', 'plasmon', 'Primakoff', or 'all'.";
    throw XUnsupportedOption(err_msg);
  }
}

std::vector<std::vector<double> > py11_calc_spectral_flux_up_to_rmax(std::vector<double> ergs, double rmax, SolarModel *s, std::string output_file_root, std::string process) {
  int n_ergs = ergs.size();
  double erg_min = *std::min_element(ergs.begin(), ergs.end());
  double erg_max = *std::max_element(ergs.begin(), ergs.end());

  if ( (rmax < s->get_r_lo()) || (rmax > s->get_r_hi()) ) { std::cout << "WARNING! Changing min. and/or max. radius to min./max. radius available in the selected Solar model." << std::endl; }
  // Check min/max and avoid Python roundoff errors
  rmax = std::max(std::min(rmax, 0.999999*s->get_r_hi()), 1.000001*s->get_r_lo());

  std::cout << "INFO. Performing calculation for ";
  if (n_ergs > 1) { std::cout << n_ergs << " energies in [" << erg_min << ", " << erg_max << "] keV and "; } else { std::cout << "one energy value (" << erg_min << " keV) and "; }
  std::cout << "up to rmax = " << rmax << " R_sol." << std::endl;

  std::vector<std::vector<double> > res;
  if (process == "Primakoff") {
    res = integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(ergs, rmax, *s, output_file_root+"_P.dat");
  } else if (process == "ABC") {
    res = integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(ergs, rmax, *s, output_file_root+"_ABC.dat");
  } else {
    std::string err_msg = "The process '"+process+"' is not a valid option. Choose 'ABC' or 'Primakoff'.";
    throw XUnsupportedOption(err_msg);
  }
  // Remove the zero entries and rvals column
  // TODO. This should be removed from the library itself.
  res.erase(res.begin());
  int n_half = res[0].size()/2;
  res[0].erase(res[0].begin(), res[0].begin()+n_half);
  res[1].erase(res[1].begin(), res[1].begin()+n_half);
  return res;
}

std::vector<std::vector<double> > py11_calc_integrated_flux_up_to_different_radii(std::vector<double> radii, SolarModel *s, std::string output_file_root, std::vector<double> erg_limits, std::string process) {
  int n_radii = radii.size();
  double rad_min = *std::min_element(radii.begin(), radii.end());
  double rad_max = *std::max_element(radii.begin(), radii.end());

  if ( (rad_min < s->get_r_lo()) || (rad_max > s->get_r_hi()) ) { std::cout << "WARNING! Changing min. and/or max. radius to min./max. radius available in the selected Solar model." << std::endl; }
  // Check min/max and avoid Python roundoff errors
  rad_min = std::min(std::max(rad_min, 1.000001*s->get_r_lo()), 0.999999*s->get_r_hi());
  rad_max = std::max(std::min(rad_max, 0.999999*s->get_r_hi()), 1.000001*s->get_r_lo());
  std::cout << "INFO. Performing calculation for energies in [" << erg_limits[0] << ", " << erg_limits[1] << "] keV and ";
  if (n_radii > 1) { std::cout << n_radii << " radii in [" << rad_min << ", " << rad_max << "] R_sol." << std::endl; } else { std::cout << "one radius value (" << rad_min << " R_sol)." << std::endl; }

  std::vector<std::vector<double> > result;
  std::string saveas_p = (output_file_root != "") ? output_file_root+"_P.dat" : "";
  std::string saveas_e = (output_file_root != "") ? output_file_root+"_ABC.dat" : "";
  std::string saveas_pl = (output_file_root != "") ? output_file_root+"_plasmon.dat" : "";
  if (process == "Primakoff") {
    result = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_Primakoff, saveas_p);
  } else if (process == "ABC") {
    result = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_all_electron, saveas_e);
  } else if (process == "plasmon") {
    result = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_plasmon, saveas_pl);
  } else if (process == "all") {
    result = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_Primakoff, saveas_p);
    std::vector<std::vector<double> > temp1 = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_all_electron, saveas_e);
    result.push_back(temp1[1]);
    std::vector<std::vector<double> > temp2 = integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(erg_limits[0], erg_limits[1], radii, *s, &SolarModel::Gamma_plasmon, saveas_pl);
    result.push_back(temp2[1]);
  } else {
    std::string err_msg = "The process '"+process+"' is not a valid option. Choose 'ABC', 'plasmon', 'Primakoff', or 'all'.";
    throw XUnsupportedOption(err_msg);
  }

  return result;
}

void py11_save_varied_spectral_flux(std::vector<double> ergs, std::string solar_model_file, std::string output_file_root, double a, double b, std::vector<double> c) {
  SolarModel s (solar_model_file, OP);
  s.set_opacity_correction(a, b);
  s.set_bfields(c[0], c[1], c[2]);

  std::vector<double> check_1 = s.get_opacity_correction();
  std::vector<double> check_2 = s.get_bfields();
  std::cout << "INFO. Setting up Solar model from file " << solar_model_file << " and opacity correction parameters (" << check_1[0] << ", " << check_1[1] << ") "\
               "and B-fields (" << check_2[0] << ", " << check_2[1] << ", " << check_2[2] << ")." << std::endl;

  std::string output_file = output_file_root+"_Primakoff.dat";
  fully_integrate_d2Phi_a_domega_drho_in_rho_Primakoff(ergs, s, output_file);
  output_file = output_file_root+"_ABC.dat";
  fully_integrate_d2Phi_a_domega_drho_in_rho_axionelectron(ergs, s, output_file);
  if (c[0]+c[1]+c[2] > 0) {
    output_file = output_file_root+"_plasmon.dat";
    fully_integrate_d2Phi_a_domega_drho_in_rho_plasmon(ergs, s, output_file);
  }
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
