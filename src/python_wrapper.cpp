#include "python_wrapper.hpp"

void py11_save_Primakoff_spectral_flux_for_different_radii(double erg_min, double erg_max, int n_ergs, double rad_min, double rad_max, int n_rads, std::string solar_model_file, std::string output_file_root) {
  // Setup steps
  double erg_stepsize = (erg_max - erg_min)/double(n_ergs-1);
  double rad_stepsize = (rad_max - rad_min)/double(n_rads-1);
  std::cout << "Setting up Solar model..." << std::endl;
  SolarModel s (solar_model_file,OP);
  std::cout << "Done!" << std::endl;
  std:: vector<double> ergs;
  for (int i = 0; i < n_ergs; i++) { ergs.push_back(erg_min + i*erg_stepsize); };

  // Do the calculation for each radius
  for (int j = 0; j < n_rads; j++) {
    double r = rad_min + j*rad_stepsize;
    std::string output_file = output_file_root+"_"+std::to_string(j)+".dat";
    calculate_spectral_flux_Primakoff(ergs, s, r, output_file);
  };
}

PYBIND11_MODULE(pyaxionflux, m) {
    m.doc() = "Solar Axion Flux library functions";

    m.def("Primakoff_flux", &integrated_Primakoff_flux_from_file, "Integrated Primakoff flux from file with signature (erg_min, erg_max, filename).");
    m.def("integrand_Primakoff_flux", &py11_save_Primakoff_spectral_flux_for_different_radii, "Integrated Primakoff flux from file with signature (erg_min, erg_max, n_ergs, rad_min, rad_max, n_rads, solar_model_file, output_file_root).");
}
