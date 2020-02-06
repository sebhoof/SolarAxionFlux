#include <iostream>
#include <vector>
#include <chrono>

#include "utils.hpp"
#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

int main() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "# Testing the Solar Model routines..." << std::endl;

  std::string solar_model_name = "data/SolarModel_AGSS09met.dat";
  SolarModel s (solar_model_name,OP,true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "# Setting up the Solar model '" << solar_model_name << "' took "
            << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count() << " seconds." << std::endl;
  const int n_test_values = 1000;
  std:: vector<double> test_ergs;
  for (int k=0; k<n_test_values; k++ ) { test_ergs.push_back(0.1+11.9/n_test_values*(k)); };
  ASCIItableReader javis_data("results/2013_redondo_all.dat");
  std::vector<double> javis_ergs = javis_data[0];

  std::cout << "# Calculating Primakoff spectrum..." << std::endl;
  calculate_spectral_flux_Primakoff(test_ergs, s, "primakoff");
  std::cout << "# Calculating Compton spectrum..." << std::endl;
  calculate_spectral_flux_Compton(test_ergs, s, "compton");
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating FF spectrum..." << std::endl;
  calculate_spectral_flux_all_ff(test_ergs, s, "all_ff");
  auto t4 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the FF spectrum took "
            << std::chrono::duration_cast<std::chrono::seconds>(t4-t3).count() << " seconds." << std::endl;
  std::cout << "# Compute opacity contribution (only metals in OP case)..." << std::endl;
  calculate_spectral_flux_opacity(test_ergs, s, "metals");
  auto t5 = std::chrono::high_resolution_clock::now();
  std::cout << "# Compute full axion-electron spectrum..." << std::endl;
  calculate_spectral_flux_axionelectron(test_ergs, s, "all_gaee");
  auto t6 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the full axion-electron spectrum (" << n_test_values << " energy values) took "
            << std::chrono::duration_cast<std::chrono::seconds>(t6-t5).count() << " seconds." << std::endl;
  std::cout << "# Compute counts in CAST2007 experiment from axion-photon interactions..." << std::endl;
  axion_photon_counts_full(1.0e-3, 1.0e-10, &cast_2007_setup, &s);
  auto t7 = std::chrono::high_resolution_clock::now();
  std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::seconds>(t7-t6).count() << " seconds." << std::endl;
  //std::cout << "# Compute counts in CAST2007 experiment from axion-electron interactions..." << std::endl;
  //axion_electron_counts_full(1.0e-3, 1.0e-13, 1.0e-10, &setup, &s);
  //auto t8 = std::chrono::high_resolution_clock::now();
  //std::cout << "# Calculating the counts took " << std::chrono::duration_cast<std::chrono::minutes>(t8-t7).count() << " minutes." << std::endl;
  std::cout << "# Compute counts in CAST2007 experiment from axion-electron interactions (from spectrum file)..." << std::endl;
  axion_electron_counts(1.0e-3, 1.0e-13, 1.0e-10, &cast_2007_setup, "results/all_gaee.dat");

  double lowerg = 0.1;
  double higherg = 10.0;
//  std::cout << "# Calculating full flux between " << lowerg << " keV and " << higherg << " keV." << std::endl;
//  std::cout << calculate_flux(lowerg,higherg,s,0);
  std::cout << "# Finished testing!" << std::endl;
  return 0;
}
//compare 4 different op codes
/*
int main() {
//compare different opacities in each 1 kev bin
    std::string solar_model_name = "data/SolarModel_AGSS09.dat";
    std::cout << "# compare different opacities in each 1 kev bin" << std::endl;
    std:: cout << "# building solar models..." << std::endl;
    SolarModel sOPAS (solar_model_name,OPAS,true);
    SolarModel sOP (solar_model_name,OP,true);
    SolarModel sATOMIC (solar_model_name,ATOMIC,true);
    SolarModel sLEDCOP (solar_model_name,LEDCOP,true);
//computing fluxes
    std:: cout << "# computing fluxes..." << std::endl;
    std::vector<double> erg_boundaries = {0.1};
    for (int k=1; k<13;k++ ) {erg_boundaries.push_back(double(k));}
    std::ofstream output;
    output.open("results/compare_opacities.dat");
    output << "fluxes in each bin for four opacity codes" << std::endl;
    output << "OP OPAS ATOMIC LEDCOP" << std::endl;
    for (int k=0; k<12; k++) {
        output << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sOP,0) << " "
            << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sOPAS,0) << " "
            << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sATOMIC,0) << " "
            << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sLEDCOP,0) << std::endl;
    }
    output.close();
}
 */
// compare five different solar models
/*
int main() {
 //compare different solarmodels in each 1 kev bin
     std::string solar_model_name = "data/SolarModel_AGSS09.dat";
     std::cout << "# compare different opacities in each 1 kev bin" << std::endl;
     std:: cout << "# building solar models..." << std::endl;
     SolarModel sAGSS09 ("data/SolarModel_AGSS09.dat",OP,true);
     SolarModel sAGSS09ph ("data/SolarModel_AGSS09ph.dat",OP,true);
     SolarModel sGS98 ("data/SolarModel_GS98.dat",OP,true);
    SolarModel sAGS05 ("data/SolarModel_AGS05.dat",OP,true);
    SolarModel sJavi ("data/SolarModel_Javi_red.dat",OP,true);
 //computing fluxes
     std:: cout << "# computing fluxes..." << std::endl;
     std::vector<double> erg_boundaries = {0.1};
     for (int k=1; k<13;k++ ) {erg_boundaries.push_back(double(k));}
     std::ofstream output;
     output.open("results/compare_SolarModels.dat");
     output << "fluxes in each bin for five solar models" << std::endl;
     output << "AGSS09 AGSS09ph GS98 AGS05 JaviRed" << std::endl;
     for (int k=0; k<12; k++) {
         output << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sAGSS09,0) << " "
             << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sAGSS09ph,0) << " "
             << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sGS98,0) << " "
            << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sAGS05,0) << " "
             << calculate_flux(erg_boundaries[k],erg_boundaries[k+1],sJavi,0) << std::endl;
     }
     output.close();
}
*/
