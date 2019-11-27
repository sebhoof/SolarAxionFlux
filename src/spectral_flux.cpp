#include "spectral_flux.hpp"

double integrand_Primakoff(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return erg*erg*(sol->Gamma_P_Primakoff(erg, r));
}

double integrand_Compton(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return erg*erg*(sol->Gamma_P_Compton(erg, r));
}

double integrand_element(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  int iz = (p->iz);
  SolarModel* sol = (p->sol);

  return erg*erg*(sol->Gamma_P_element(erg, r, iz));
}


void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz) {
  // = Rsol [in keV^-1] / (2 pi^2 d^2 [in m^2]) * 1; express 1 as 1/(year*keV)
  const double factor = pow(radius_sol/(1.0e7*gev2cm),3)/(2.0*pi*pi*pow(distance_sol,2)) / (1.0e6*hbar/(60.0*60.0*24.0*365.0));
  double flux;

  gsl_function f;
  f.function = integrand;

  for (int i=0; i<ergs.size(); ++i) {
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6);
    integration_params p = {ergs[i], &s, iz};
    f.params = &p;
    gsl_integration_qag (&f, s.r_lo, s.r_hi, abs_prec1, rel_prec1, 1e6, method1, w, &result, &error);
    flux = factor*result;
    std::cout << ergs[i] << " " << flux << std::endl;
  };
}

void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*)) { calculate_spectral_flux(ergs, s, integrand, 0); }

void calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_Primakoff); }

void calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_Compton); }

void calculate_spectral_flux_element(std::vector<double> ergs, SolarModel &s, int iz) { calculate_spectral_flux(ergs, s, &integrand_element, iz); }
