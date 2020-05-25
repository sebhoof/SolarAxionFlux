#include "python_wrapper.hpp"

void test_module() {
  std::cout << "This is the " << LIBRARY_NAME << " Python interface." << std::endl;
  std::cout << "Test successful!" << std::endl;
}

void py11_save_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root, std::string process) {
  // Setup steps
  std::cout << "Setting up Solar model..." << std::endl;
  SolarModel s (solar_model_file,OP);
  std::cout << "Done!" << std::endl;

  std:: vector<double> ergs;
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
  SolarModel s (solar_model_file,OP);
  s.set_opacity_correction(a, b);
  std::cout << "Done!" << std::endl;

  std:: vector<double> ergs;
  double erg_stepsize = (erg_max - erg_min)/double(n_ergs);
  for (int i = 0; i < n_ergs+1; i++) { ergs.push_back(erg_min + i*erg_stepsize); };

  std::string output_file = output_file_root+"_Primakoff.dat";
  calculate_spectral_flux_Primakoff(ergs, s, 1.0, output_file);
  output_file = output_file_root+"_electron.dat";
  calculate_spectral_flux_axionelectron(ergs, s, 1.0, output_file);
}

PYBIND11_MODULE(pyaxionflux, m) {
    m.doc() = "Solar Axion Flux library functions";

    m.def("test_module", &test_module, "A simple routine to test the Python interface.");
    m.def("calculate_spectrum", &py11_save_spectral_flux_for_different_radii, "Integrated Primakoff flux from file with signature (erg_min, erg_max, n_ergs, rad_min, rad_max, n_rads, solar_model_file, output_file_root, process = 'Primakoff').");
    m.def("calculate_spectra_with_varied_opacities", &py11_save_spectral_flux_for_varied_opacities, "Integrated Primakoff flux from file with signature (erg_min, erg_max, n_ergs, a, b, solar_model_file, output_file_root).");
}
