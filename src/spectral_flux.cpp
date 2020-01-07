#include "spectral_flux.hpp"

double integrand_Primakoff(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_Primakoff(erg, r));
}

double integrand_Compton(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_Compton(erg, r));
}

double integrand_weightedCompton(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);
  double u = erg/(sol->temperature_in_keV(r));

  return 0.5*gsl_pow_2(r*erg/pi)*0.5*(1.0 - 1.0/gsl_expm1(u))*(sol->Gamma_P_Compton(erg, r));
}

double integrand_element(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  int iz = (p->iz);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_element(erg, r, iz));
}

// Includes FF flux and ee contribution
double integrand_all_ff(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  double sum = 0.0;
  //sum += sol->Gamma_P_ff(erg, r);
  for (int iz = 0; iz < 2; iz++) { sum += sol->Gamma_P_ff(erg, r, iz); };
  return 0.5*gsl_pow_2(r*erg/pi)*(sum + sol->Gamma_P_ee(erg, r));
}

double integrand_all_axionelectron(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);
  double u = erg/(sol->temperature_in_keV(r));

  double element_contrib = 0.0;
  // Add H, He contributions from ff approximation
  //  element_contrib += sol->Gamma_P_ff(erg, r);
  for (int iz = 0; iz < 2; iz++) { element_contrib += sol->Gamma_P_ff(erg, r, iz); };
  // Add opacity terms all non-H or He elements (metals)
  for (int iz = 2; iz < n_op_elements; iz++) { element_contrib += sol->Gamma_P_element(erg, r, iz); };
  return 0.5*gsl_pow_2(r*erg/pi)*(element_contrib + sol->Gamma_P_Compton(erg, r) + sol->Gamma_P_ee(erg, r));
}

// same function with safe method
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz,std::string saveas) {
  // = Rsol [in keV^-1] / (2 pi^2 d^2 [in m^2]) * 1; express 1 as 1/(year*keV)
  // Better: per cm^2 per s per keV?
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar/(60.0*60.0*24.0*365.0)));
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1.
  // Integrate function *1/norm and rescale result *norm at the end

  double result, error;
  gsl_function f;
  f.function = integrand;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6);
  std::ofstream output;
  if (saveas != ""){
      output.open("results/"+saveas+".dat");}
  for (int i=0; i<ergs.size(); ++i) {
    integration_params p = {ergs[i], &s, iz};
    f.params = &p;
    gsl_integration_qag (&f, s.r_lo, s.r_hi, abs_prec1, rel_prec1, 1e6, method1, w, &result, &error);
    if (saveas!= ""){ output << ergs[i] << " " << factor*result << std::endl;}
    
    //std::cout << ergs[i] << " " << factor*result << std::endl;
  };
  if (saveas!= ""){ output.close();}
  gsl_integration_workspace_free (w);
}
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz) {
 
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar/(60.0*60.0*24.0*365.0)));
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1.
  // Integrate function *1/norm and rescale result *norm at the end

  double result, error;
  gsl_function f;
  f.function = integrand;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6);

  for (int i=0; i<ergs.size(); ++i) {
    integration_params p = {ergs[i], &s, iz};
    f.params = &p;
    gsl_integration_qag (&f, s.r_lo, s.r_hi, abs_prec1, rel_prec1, 1e6, method1, w, &result, &error);
  };
  gsl_integration_workspace_free (w);
}
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*)) { calculate_spectral_flux(ergs, s, integrand, 0); }
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*),std::string saveas) { calculate_spectral_flux(ergs, s, integrand, 0,saveas); }

void calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_Primakoff); }
void calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_Primakoff,saveas); }
void calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_Compton); }
void calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_Compton,saveas); }
void calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_weightedCompton); }
void calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_weightedCompton,saveas); }
void calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz) { calculate_spectral_flux(ergs, s, &integrand_element, iz); }
void calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz,std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_element, iz,saveas); }
void calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_all_ff); }
void calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s,std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_all_ff,saveas); }
void calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_all_axionelectron); }
void calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s,std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_all_axionelectron,saveas); }
