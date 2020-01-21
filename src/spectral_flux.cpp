#include "spectral_flux.hpp"

double integrand_Primakoff(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_Primakoff(erg, r));
};

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
  if (erg == 0) {return 0;}
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
  if (sol->raffelt_approx == false) {
        for (int iz = 0; iz < 2; iz++) { sum += sol->Gamma_P_ff(erg, r, iz); };
  } else {
        sum += sol->Gamma_P_ff(erg, r);
  }
  return 0.5*gsl_pow_2(r*erg/pi)*(sum + sol->Gamma_P_ee(erg, r));
}

double integrand_all_axionelectron(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);
  if (sol->opcode == OP){
      double element_contrib = 0.0;
      double u = erg/(sol->temperature_in_keV(r));
      // Add H, He contributions from ff approximation
      //  element_contrib += sol->Gamma_P_ff(erg, r);
      if (sol->raffelt_approx == false) {
            for (int iz = 0; iz < 2; iz++) { element_contrib += sol->Gamma_P_ff(erg, r, iz); };
      } else {
            element_contrib += sol->Gamma_P_ff(erg, r);
      }
      // Add opacity terms all non-H or He elements (metals)
      for (int iz = 2; iz < n_op_elements; iz++) { element_contrib += sol->Gamma_P_element(erg, r, iz); };
      return 0.5*gsl_pow_2(r*erg/pi)*(element_contrib + sol->Gamma_P_Compton(erg, r) + sol->Gamma_P_ee(erg, r));
  }
  if ((sol->opcode == LEDCOP) || (sol->opcode == ATOMIC)){
      double u = erg/(sol->temperature_in_keV(r));
      double reducedCompton = 0.5*(1.0 - 1.0/gsl_expm1(u)) * sol->Gamma_P_Compton(erg, r);
      return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_opacity (erg,r)+reducedCompton + sol->Gamma_P_ee(erg, r));
  }
    if (sol->opcode == OPAS) {
        double u = erg/(sol->temperature_in_keV(r));
        double reducedCompton = 0.5*(1.0 - 1.0/gsl_expm1(u)) * sol->Gamma_P_Compton(erg, r);
        return 0.5*gsl_pow_2(r*erg/pi)*sol->Gamma_P_opacity (erg,r);
    }
  return 0;
}
double integrand_opacity(double r, void * params) {
  // Retrieve parameters and other integration variables.
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);
  if (sol->opcode == OP){
      double element_contrib = 0.0;
      // Add opacity terms all non-H or He elements (metals)
      for (int iz = 2; iz < n_op_elements; iz++) { element_contrib += sol->Gamma_P_element(erg, r, iz); };
      return 0.5*gsl_pow_2(r*erg/pi)*element_contrib ;
  }
  if ((sol->opcode == OPAS) || (sol->opcode == LEDCOP) || (sol->opcode == ATOMIC) ){
      return 0.5*gsl_pow_2(r*erg/pi)* sol->Gamma_P_opacity (erg,r) ;
  }
  return 0;
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

/*
// Javis integrator to check how much different it makes
void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz,std::string saveas) {
  // = Rsol [in keV^-1] / (2 pi^2 d^2 [in m^2]) * 1; express 1 as 1/(year*keV)
  // Better: per cm^2 per s per keV?
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar/(60.0*60.0*24.0*365.0)));
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1.
  // Integrate function *1/norm and rescale result *norm at the end
  std::ofstream output;
  if (saveas != ""){
      output.open("results/"+saveas+".dat");}
  for (int i=0; i<ergs.size(); ++i) {
      double result = 0;
      integration_params p = {ergs[i], &s, iz};
      for (int k = 0; k < s.data.getnrow()-1; k++){
          double rlow = s.data["radius"][k];
          double rhigh = s.data["radius"][k+1];
          double integrandlow, integrandhigh;
          double ulow = ergs[i]/(s.data["temperature"][k]*(1.0E-3*K2eV));
          double uhigh = ergs[i]/(s.data["temperature"][k+1]*(1.0E-3*K2eV));
          if ((0.12 < ulow) && (ulow < 18)) {
              integrandlow = integrand(rlow,&p);
          } else {
              integrandlow = 0;
          }
          if ((0.12 < uhigh) && (uhigh < 18)) {
              integrandhigh = integrand(rhigh,&p);
          } else {
              integrandhigh = 0;
          }
          result += 0.5*(integrandlow + integrandhigh) * (rhigh-rlow);
      }
      if (saveas!= ""){ output << ergs[i] << " " << factor*result << std::endl;}
  };
  if (saveas!= ""){ output.close();}
}
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz) {calculate_spectral_flux_javi(ergs, s, integrand,iz,"");}
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*)) { calculate_spectral_flux_javi(ergs, s, integrand, 0); }
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*),std::string saveas) { calculate_spectral_flux_javi(ergs, s, integrand, 0,saveas); }
 */

 // will become integrator for spectral flux to compute total flux in some energy range
double spectral_flux_integrand(double erg, void * params ){
    struct integration_params2 * p2 = (struct integration_params2 *)params;
    int iz = (p2->iz);
    SolarModel* s = (p2->sol);
    const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar/(60.0*60.0*24.0*365.0)));
    const double normfactor = 1.0e20;
    double result, error;
    gsl_function f;
    f.function = p2->integrand;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e6);
    integration_params p = {erg, s, iz};
    f.params = &p;
    gsl_integration_qag (&f, s->r_lo, s->r_hi, abs_prec1, rel_prec1, 1e6, method1, w, &result, &error);
    gsl_integration_workspace_free (w);
    return factor*result/normfactor;
}

double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s,int iz){
    const double normfactor = 1.0e20;
    double result, error;
    gsl_function f;
    f.function = spectral_flux_integrand;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e8);
    integration_params2 p2 = {&s,&integrand_all_axionelectron,iz};
    f.params = &p2;
    gsl_integration_qag (&f, lowerlimit, upperlimit, abs_prec2, rel_prec2, 1e8, method2, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result*normfactor;
}
void calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), int iz) {calculate_spectral_flux(ergs, s, integrand,iz,"");}
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
void calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s) { calculate_spectral_flux(ergs, s, &integrand_opacity); }
void calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s,std::string saveas) { calculate_spectral_flux(ergs, s, &integrand_opacity,saveas); }
