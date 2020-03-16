#include "spectral_flux.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Various integrands for the different contributions/combinations of contributions to the solar axion flux. //
// All in units of axions / (cm^2 s keV).                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Primakoff contribution [ref].
double integrand_Primakoff(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_Primakoff(erg, r));
}

// Compton contribution [ref]
double integrand_Compton(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_Compton(erg, r));
}

// Weighted Compton contribution [ref]
double integrand_weightedCompton(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  if (erg == 0) {return 0;}
  SolarModel* sol = (p->sol);
  double u = erg/(sol->temperature_in_keV(r));

  return 0.5*gsl_pow_2(r*erg/pi)*0.5*(1.0 - 1.0/gsl_expm1(u))*(sol->Gamma_P_Compton(erg, r));
}

double integrand_opacity_element(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  std::string el_name = (p->isotope).name();
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_opacity(erg, r, el_name));
}

double integrand_opacity(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double result = 0.0;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  if (sol->opcode == OP) {
    double element_contrib = 0.0;
    // Add opacity terms all non-H or He elements (metals)
    for (int k = 2; k < num_op_elements; k++) { element_contrib += sol->Gamma_P_opacity(erg, r, op_element_names[k]); };
    result = 0.5*gsl_pow_2(r*erg/pi)*element_contrib;
  }
  if ((sol->opcode == OPAS) || (sol->opcode == LEDCOP) || (sol->opcode == ATOMIC)) {
      result = 0.5*gsl_pow_2(r*erg/pi) * sol->Gamma_P_opacity(erg, r);
  };

  return result;
}

// Includes FF flux and ee contribution [ref].
double integrand_all_ff(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_ff(erg, r) + sol->Gamma_P_ee(erg, r));
}

double integrand_all_axionelectron(double r, void * params) {
  struct integration_params * p = (struct integration_params *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);
  if (sol->opcode == OP) {
      double element_contrib = 0.0;
      if (sol->raffelt_approx == false) {
        element_contrib += sol->Gamma_P_opacity(erg, r, "H") + sol->Gamma_P_opacity(erg, r, "He");
      } else {
        element_contrib += sol->Gamma_P_ff(erg, r);
      };
      for (int i = 2; i < num_op_elements; i++) { element_contrib += sol->Gamma_P_opacity(erg, r, op_element_names[i]); };
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

std::vector<double> calculate_spectral_flux(std::vector<double> ergs, Isotope isotope, SolarModel &s, double (*integrand)(double, void*), std::string saveas) {
  // Constant factor for consistent units, i.e. integrated flux will be in units of cm^-2 s^-1 keV^-1.
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / ( pow(1.0e2*distance_sol,2) * (1.0e6*hbar) );
  // = Rsol^3 [in keV^-3] / (2 pi^2 d^2 [in cm^2] * 1 [1 corresponds to s x keV))
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1. Integrate function *1/norm and rescale result *norm at the end.
  std::vector<double> result;

  gsl_function f;
  f.function = integrand;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);

  std::ofstream output;
  if (saveas != "") { output.open(saveas); };

  for (int i=0; i<ergs.size(); ++i) {
    double integral, error;
    integration_params p = {ergs[i], &s, isotope};
    f.params = &p;
    gsl_integration_qag (&f, s.r_lo, s.r_hi, int_abs_prec, int_rel_prec, int_space_size, int_method_1, w, &integral, &error);
    result.push_back(factor*integral);
    if (saveas!= ""){ output << ergs[i] << " " << factor*integral << std::endl; };
  };

  if (saveas!= "") { output.close(); };
  gsl_integration_workspace_free (w);

  return result;
}

double rho_integrand(double rho, void * params) {
  // Retrieve parameters and other integration variables.
  struct solar_disc_integration_params * p1 = (struct solar_disc_integration_params *)params;
  double erg = (p1->erg);
  double rad = (p1->rad);
  SolarModel *s = p1->s;

  // Normalising factor ~ max values, which occur for rho = r_lo
  // double norm_factor1 = (s->*(p1->integrand))(ref_erg_value, s->r_lo);
  const double norm_factor1 = 1.0;
  double cylinder = rho*rho - rad*rad;
  cylinder = rho/sqrt(cylinder);

  //std::cout << "# DEBUG INFO: rho = " << rho << ", erg = " << erg << ", res = " << cylinder * ( (s->*(p1->integrand))(erg, rho) ) / norm_factor1 << std::endl;

  //return cylinder*(p1->integrand(erg, rho));
  return cylinder * ( (s->*(p1->integrand))(erg, rho) ) / norm_factor1;
}

double rad_integrand(double rad, void * params) {
  struct solar_disc_integration_params * p2 = (struct solar_disc_integration_params *)params;
  p2->rad = rad;
  SolarModel *s = p2->s;
  double r_max = std::min(p2->r_max, s->r_hi);

  gsl_function f1;
  f1.function = &rho_integrand;
  f1.params = p2;

  double result, error;
  //gsl_integration_qag(&f1, rad, r_max, 0.01*int_abs_prec, 0.01*int_rel_prec, int_space_size, int_method_1, p2->w1, &result, &error);
  //gsl_integration_qags(&f1, rad, r_max, 0.01*int_abs_prec, 0.01*int_rel_prec, int_space_size, p2->w1, &result, &error);
  gsl_integration_cquad(&f1, rad, r_max, 0.01*int_abs_prec, 0.01*int_rel_prec, p2->w1, &result, &error, NULL);
  //std::cout << "# DEBUG INFO: rad = " << rad << ", integral 1 = " << result << std::endl;

  result = rad*result;
  return result;
}

std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, Isotope isotope, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas) {
  // Constant factor for consistent units, i.e. integrated flux will be in units of cm^-2 s^-1 keV^-1.
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / ( pow(1.0e2*distance_sol,2) * (1.0e6*hbar) );
  // = Rsol^3 [in keV^-3] / (2 pi^2 d^2 [in cm^2] * 1 [1 corresponds to s x keV))
  std::vector<double> result;

  //gsl_integration_workspace * w1 = gsl_integration_workspace_alloc(int_space_size);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_cquad);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc(int_space_size);

  double r_min = s.r_lo;
  r_max = std::min(r_max, s.r_hi);
  //double norm_factor1 = (s.*integrand)(ref_erg_value, r_min);
  const double norm_factor1 = 1.0;

  //solar_disc_integration_params p2 { 0.0, 0.0, r_max, &s, integrand, w1 };
  //double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff;
  //solar_disc_integration_params p2 { 0.0, 0.0, r_max, &s, func_ptr, w1 };
  solar_disc_integration_params p2 { 0.0, 0.0, r_max, &s, integrand, w1 };
  gsl_function f2;
  f2.function = &rad_integrand;

  std::ofstream output;
  if (saveas != "") { output.open(saveas); };

  //std::cout << "# DEBUG INFO: r in [" << r_min << ", " << r_max << "] ..." << std::endl;

  for (int i=0; i<ergs.size(); ++i) {
    double integral, error;
    p2.erg = ergs[i];
    f2.params = &p2;
    // was 0.1*factor
    gsl_integration_qag (&f2, r_min, r_max, int_abs_prec, int_rel_prec, int_space_size, int_method_1, w2, &integral, &error);
    result.push_back(factor*norm_factor1*integral);
    if (saveas!= ""){ output << ergs[i] << " " << factor*norm_factor1*integral << std::endl; };
  };

  if (saveas!= "") { output.close(); };
  //gsl_integration_workspace_free (w1);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_workspace_free (w2);

  return result;
};

/*
// Javis integrator to check how much different it makes
void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), Isotope isotope, std::string saveas) {
  std::ofstream output;
  if (saveas != ""){
      output.open("results/"+saveas+".dat");}
  for (int i=0; i<ergs.size(); ++i) {
      double result = 0;
      integration_params p = {ergs[i], &s, isotope};
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
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), Isotope isotope) {calculate_spectral_flux_javi(ergs, s, integrand,iz,"");}
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*)) { calculate_spectral_flux_javi(ergs, s, integrand, 0); }
 void calculate_spectral_flux_javi(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*),std::string saveas) { calculate_spectral_flux_javi(ergs, s, integrand, 0,saveas); }
 */

 // Generic integrator to compute the spectral flux in some energy range.
double spectral_flux_integrand(double erg, void * params) {
  // Constant factor for consistent units, i.e. integrated flux will be in units of cm^-2 s^-1 keV^-1.
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / ( pow(1.0e2*distance_sol,2) * (1.0e6*hbar) );
  // = Rsol^3 [in keV^-3] / (2 pi^2 d^2 [in cm^2] * 1 [1 corresponds to s x keV))
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1. Integrate function *1/norm and rescale result *norm at the end.
  struct integration_params2 * p2 = (struct integration_params2 *)params;
  Isotope isotope = (p2->isotope);
  SolarModel* s = (p2->sol);
  const double normfactor = 1.0;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  double result, error;
  gsl_function f;
  f.function = p2->integrand;
  integration_params p = {erg, s, isotope};
  f.params = &p;
  gsl_integration_qag (&f, s->r_lo, s->r_hi, int_abs_prec, int_rel_prec, int_space_size, int_method_1, w, &result, &error);
  gsl_integration_workspace_free (w);
  return factor*result/normfactor;
}

double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s, Isotope isotope) {
    const double normfactor = 1.0e20;
    double result, error;
    gsl_function f;
    f.function = spectral_flux_integrand;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1e8);
    integration_params2 p2 = {&s, &integrand_all_axionelectron, isotope};
    f.params = &p2;
    gsl_integration_qag (&f, lowerlimit, upperlimit, abs_prec2, rel_prec2, 1e8, int_method_2, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result*normfactor;
}

double flux_from_file_integrand(double erg, void * params) {
  OneDInterpolator * interp = (OneDInterpolator *)params;
  //std::cout << "DEBUG INFO. flux_from_file_integrand(" << erg << " keV) = " << interp->interpolate(erg) << " ." << std::endl;
  return interp->interpolate(erg);
}

double integrated_Primakoff_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file) {
  double result, error;

  OneDInterpolator spectral_flux (spectral_flux_file);
  if ( (erg_min < spectral_flux.lower()) || (erg_max > spectral_flux.upper()) ) {
    terminate_with_error("ERROR! The integration boundaries given to 'integrated_flux_from_file' are incompatible with the min/max available energy in the file "+spectral_flux_file+".");
  };

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  gsl_function f;
  f.function = &flux_from_file_integrand;
  f.params = &spectral_flux;

  gsl_integration_qag(&f, erg_min, erg_max, abs_prec2, rel_prec2, int_space_size, int_method_1, w, &result, &error);

  gsl_integration_workspace_free (w);

  return result;
}


////////////////////////////////////////////////////////////////////
// Overloaded versions of the functions above for convenient use. //
////////////////////////////////////////////////////////////////////

//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas) { return calculate_spectral_flux_solar_disc(ergs, r_max, 0, s, integrand, saveas); }
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,Isotope isotope, double r_max, SolarModel &s, double (*integrand)(double, double)) { return calculate_spectral_flux_solar_disc(ergs, r_max, isotope, s, integrand, ""); }
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,double r_max, SolarModel &s, double (*integrand)(double, double)) { return calculate_spectral_flux_solar_disc(ergs, r_max, 0, s, integrand); }
std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas) { std::string NONE = ""; return calculate_spectral_flux_solar_disc(ergs, NONE, r_max, s, integrand, saveas); }
std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas) { std::string NONE = ""; return calculate_spectral_flux(ergs, NONE, s, integrand, saveas); }
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_Primakoff, saveas); }
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas) { double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff; return calculate_spectral_flux_solar_disc(ergs, r_max, s, integrand, saveas); }
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_Compton,saveas); }
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_weightedCompton,saveas); }
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, std::string element, SolarModel &s) { return calculate_spectral_flux(ergs, element, s, &integrand_opacity_element); }
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, std::string element, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, element, s, &integrand_opacity_element, saveas); }
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_all_ff,saveas); }
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s,std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_all_axionelectron,saveas); }
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_opacity, saveas); }
