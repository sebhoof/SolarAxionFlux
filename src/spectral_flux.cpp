// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "spectral_flux.hpp"

/////////////////////////////////////////////////////
//  Integration routines for the solar axion flux  //
/////////////////////////////////////////////////////

// N.B. The solar axion fluxes will be calculated in units of [axions cm^-2 s^-1 keV^-1]

std::string standard_header(SolarModel *s) {
  double g_ag = s->get_gagg_ref_value_in_inverse_GeV();
  double g_ae = s->get_gaee_ref_value();
  std::string timedate = current_time_string();
  std::stringstream res;
  res << timedate << " | Calculation performed with " LIBRARY_NAME "; g_agamma = " << std::scientific << g_ag << " GeV^-1, g_ae = " << g_ae << ".\n";
  return res.str();
};

// Below are all standard integration routines
double rho_integrand_1d(double rho, void * params) {
  struct solar_model_integration_parameters_1d * p1 = (struct solar_model_integration_parameters_1d *)params;
  double erg = (p1->erg);
  SolarModel *s = p1->s;

  return 2.0 * gsl_pow_2(0.5*erg*rho/pi) * (s->*(p1->integrand))(erg, rho); // N.B. Factor of 2 from integration over theta angle.
}

double erg_integrand_1d(double erg, void * params) {
  double result, error;
  struct solar_model_integration_parameters_1d * p2 = (struct solar_model_integration_parameters_1d *)params;
  p2->erg = erg;

  if ((p2->integrand == &SolarModel::Gamma_LP) || (p2->integrand == &SolarModel::Gamma_LP_Rosseland)) {
      std::vector<double> radii;
      double res = p2->s->r_from_omega_pl(erg);
      double low = p2->s->get_r_lo();
      double high = std::min(p2->s->get_r_hi(), 0.98);  // 0.99 maximum set by hand to avoid missing opacity data
      if ((res > low) && (res < high)) {radii ={low, res , high};}
      else {radii ={low, high};}
      gsl_integration_qagp(p2->f, &radii[0], radii.size(), int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, p2->w, &result, &error);}

  else {
  gsl_integration_qag(p2->f, p2->s->get_r_lo(), p2->s->get_r_hi(), int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, int_method_1d, p2->w, &result, &error);
  //static std::vector<double> radii = { p2->s->get_r_lo(), p2->s->r_from_omega_pl(erg), p2->s->get_r_hi() };
  //gsl_integration_qagp(p2->f, &radii[0], radii.size(), int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, p2->w, &result, &error);
  //gsl_integration_qags(p2->f, p2->s->get_r_lo(), 0.9, int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, p2->w, &result, &error);
  }

  return result;
}

double rho_integrand_2d(double rho, void * params) {
  struct solar_model_integration_parameters_2d * p3 = (struct solar_model_integration_parameters_2d *)params;
  double erg = (p3->erg);
  double rad = (p3->rad);

  double cylinder = rho*rho - rad*rad;
  cylinder = rho/sqrt(cylinder);

  return 2.0 * cylinder * gsl_pow_2(0.5*erg/pi)*( (p3->s->*(p3->integrand))(erg, rho) ); // N.B. Factor of 2 from Z_2 symmtery in z-axis.
}

double rad_integrand_2d(double rad, void * params) {
  struct solar_model_integration_parameters_2d * p2 = (struct solar_model_integration_parameters_2d *)params;
  p2->rad = rad;

  auto t1 = std::chrono::high_resolution_clock::now();
  double result, error;
  size_t n_evals;
  //gsl_integration_qag(&f2, rad, p2->s->get_r_hi(), 0.01*int_abs_prec_2d, 0.01*int_rel_prec_2d, int_space_size_2d, int_method_2d, p2->w1, &result, &error);
  //gsl_integration_qags(&f2, rad, p2->s->get_r_hi(), 0.1*int_abs_prec_2d, 0.1*int_rel_prec_2d, int_space_size_2d, p2->w1, &result, &error);
  gsl_integration_cquad(p2->f2, rad, 0.999999999*p2->s->get_r_hi(), 0.1*int_abs_prec_2d, 0.1*int_rel_prec_2d, p2->w2, &result, &error, &n_evals);
  auto t2 = std::chrono::high_resolution_clock::now();

  result *= rad;
  return result;
}

double erg_integrand_2d(double erg, void * params) {
  struct solar_model_integration_parameters_2d * p1 = (struct solar_model_integration_parameters_2d *)params;
  p1->erg = erg;

  double result, error;
  size_t n_evals;

  gsl_integration_cquad(p1->f1, p1->r_1, p1->r_2, int_abs_prec_2d, int_rel_prec_2d, p1->w1, &result, &error, &n_evals);
  return result;
}

std::vector<std::vector<double> > calculate_d2Phi_a_domega_drho(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas) {
  std::vector<double> all_ergs, all_radii, results;

  gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);

  std::vector<double> valid_rhos = s.get_supported_radii(rhos);
  double rho_min = valid_rhos.front();
  double rho_max = valid_rhos.back();
  int n_rho_vals = valid_rhos.size();

  gsl_function f2;
  f2.function = &rho_integrand_2d;
  solar_model_integration_parameters_2d p { 0.0, 0.0, 0.0, 0.0, &s, integrand, NULL, NULL, &f2, w2 };
  f2.params = &p;
  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
      all_radii.push_back(rho_min);
      all_ergs.push_back(*erg);
      results.push_back(0);
  }
  for (auto rho = valid_rhos.begin(); rho != valid_rhos.end(); rho++) {
    p.r_2 = *rho;
    for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
      p.erg = *erg;
      all_radii.push_back(p.r_2);
      all_ergs.push_back(p.erg);
      results.push_back( distance_factor*rad_integrand_2d(p.r_2, &p) );
    }
  }

  gsl_integration_cquad_workspace_free(w2);

  std::vector<std::vector<double> > buffer = { all_radii, all_ergs, results };
  std::string comment = standard_header(&s);
  comment += "Differential flux on the solar disc by " LIBRARY_NAME ".\nColumns: Radius on solar disc [R_sol] | Energy [keV] | Differential axion flux [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}


std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho(std::vector<double> ergs, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  std::vector<double> integrals;
  std::ofstream output;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_1d);
  gsl_function f;
  f.function = rho_integrand_1d;
  solar_model_integration_parameters_1d p { 0.0, &s, integrand, &f, w };
  f.params = &p;

  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) { integrals.push_back( distance_factor*erg_integrand_1d(*erg, &p) ); }

  gsl_integration_workspace_free(w);

  std::vector<std::vector<double> > buffer = { ergs, integrals };
  std::string comment = standard_header(&s);
  comment += "Spectral flux over full solar volume.\nColumns: energy values [keV] | axion flux [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_and_for_omega_interval(double erg_lo, double erg_hi, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas) {
  std::vector<double> results, errors;
  std::vector<double> valid_rhos = s.get_supported_radii(rhos);
  double r_min = valid_rhos.front();
  double r_max = valid_rhos.back();
  int n_r_vals = valid_rhos.size();
  std::vector<double> relevant_peaks = get_relevant_peaks(erg_lo, erg_hi);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_2d);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);
  gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);

  gsl_function f;
  f.function = &erg_integrand_2d;
  gsl_function f1;
  f1.function = &rad_integrand_2d;
  gsl_function f2;
  f2.function = &rho_integrand_2d;
  solar_model_integration_parameters_2d p { 0.0, 0.0, 0.0, 0.0, &s, integrand, &f1, w1, &f2, w2 };
  f.params = &p;
  f1.params = &p;
  f2.params = &p;

  for (int i = 0; i < n_r_vals; ++i) {
    double integral, error;

    if (i > 0) {
      p.r_1 = valid_rhos[i-1];
    } else {
      p.r_1 = s.get_r_lo();
    }

    p.r_2 = valid_rhos[i];
    if (p.r_2 > p.r_1) {
      //gsl_integration_qagiu(&f, 0.0, 10.0*int_abs_prec_2d, 10.0*int_rel_prec_2d, int_space_size_2d, w, &integral, &error); // Alternative integration from 0 -> infinity; too slow.
      gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), 10.0*int_abs_prec_2d, 10.0*int_rel_prec_2d, int_space_size_2d, w, &integral, &error);
    } else {
      integral = 0;
      error = 0;
    }

    if (i > 0) {
      results.push_back( results.back() + distance_factor*integral );
      errors.push_back( sqrt(gsl_pow_2(errors.back()) + gsl_pow_2(distance_factor*error)) );
    } else {
      results.push_back(distance_factor*integral);
      errors.push_back(distance_factor*error);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_cquad_workspace_free(w2);

  std::vector<std::vector<double>> buffer = { valid_rhos, results, errors };
  std::string comment = standard_header(&s);
  comment += "Total spectral flux for a given radius.\nColumns: Radius on solar disc [R_sol], Axion flux [cm^-2 s^-1] |  Axion flux error estimate [cm^-2 s^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  std::vector<double> all_ergs, all_radii, results;

  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);
  gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);

  std::vector<double> valid_rhos = s.get_supported_radii(rhos);
  double r_min = valid_rhos.front();
  double r_max = valid_rhos.back();
  int n_r_vals = valid_rhos.size();

  gsl_function f1;
  f1.function = &rad_integrand_2d;
  gsl_function f2;
  f2.function = &rho_integrand_2d;
  solar_model_integration_parameters_2d p { 0.0, 0.0, 0.0, 0.0, &s, integrand, &f1, w1, &f2, w2 };
  f1.params = &p;
  f2.params = &p;
  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
      all_radii.push_back(r_min);
      all_ergs.push_back(*erg);
      results.push_back(0);
  }
  for (int i = 1; i < n_r_vals; ++i) {
    p.r_1 = r_min;
    p.r_2 = valid_rhos[i];
    for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
      all_radii.push_back(p.r_2);
      all_ergs.push_back(*erg);
      results.push_back( distance_factor*erg_integrand_2d(*erg, &p) );
    }
  }

  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_cquad_workspace_free(w2);

  std::vector<std::vector<double> > buffer = { all_radii, all_ergs, results };
  std::string comment = standard_header(&s);
  comment += "Spectral flux over full solar disc at fixed radius.\nColumns: Radius on solar disc [R_sol] | Energy [keV] | Axion flux [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho(std::vector<double> ergs, double rho_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  // Check if rho_max >= biggest available radius and switch to faster 1D integration if that's the case
  if (rho_max >= s.get_r_hi()) {
    std::vector<std::vector<double> > tmp1 = fully_integrate_d2Phi_a_domega_drho_in_rho(ergs, s, integrand, saveas, isotope);
    std::vector<double> tmp2 (tmp1[0].size(), s.get_r_hi());
    std::vector<std::vector<double> > all_results = { tmp2, tmp1[0], tmp1[1] };
    return all_results;
  } else {
    std::vector<double> rhos = { rho_max };
    std::vector<std::vector<double> > all_results = integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rhos, s, integrand, saveas, isotope);
    return all_results;
  }
}

std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return fully_integrate_d2Phi_a_domega_drho_in_rho(ergs, s, &SolarModel::Gamma_Primakoff, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rhos, s, &SolarModel::Gamma_Primakoff, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_Primakoff(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rho_max, s, &SolarModel::Gamma_Primakoff, saveas);
}

std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_plasmon(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return fully_integrate_d2Phi_a_domega_drho_in_rho(ergs, s, &SolarModel::Gamma_plasmon, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_plasmon(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rhos, s, &SolarModel::Gamma_plasmon, saveas);
}

std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_plasmon(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rho_max, s, &SolarModel::Gamma_plasmon, saveas);
}

std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_axionphoton(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return fully_integrate_d2Phi_a_domega_drho_in_rho(ergs, s, &SolarModel::Gamma_all_photon, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionphoton(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rhos, s, &SolarModel::Gamma_all_photon, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionphoton(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rho_max, s, &SolarModel::Gamma_all_photon, saveas);
}

std::vector<std::vector<double> > fully_integrate_d2Phi_a_domega_drho_in_rho_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return fully_integrate_d2Phi_a_domega_drho_in_rho(ergs, s, &SolarModel::Gamma_all_electron, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(std::vector<double> ergs, std::vector<double> rhos, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rhos, s, &SolarModel::Gamma_all_electron, saveas);
}

std::vector<std::vector<double> > integrate_d2Phi_a_domega_drho_up_to_rho_axionelectron(std::vector<double> ergs, double rho_max, SolarModel &s, std::string saveas) {
  return integrate_d2Phi_a_domega_drho_up_to_rho(ergs, rho_max, s, &SolarModel::Gamma_all_electron, saveas);
}


// Here, we also define some simple, custom integration routines similar to the ones defined above
// Weighted Compton contribution
double integrand_weightedCompton(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  if (erg == 0) {return 0;}
  SolarModel* sol = (p->sol);
  double u = erg/(sol->temperature_in_keV(r));

  return 0.5*gsl_pow_2(r*erg/pi)*0.5*(1.0 - 1.0/gsl_expm1(u))*(sol->Gamma_Compton(erg, r));
}

// Includes FF flux and ee contribution as in arXiv:1310.0823
double integrand_all_ff(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_ff(erg, r) + sol->Gamma_ee(erg, r));
}

// Calculate the flux from opacity for one element only
double integrand_opacity_element(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  std::string el_name = (p->isotope).get_element_name();
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_opacity(erg, r, el_name));
}

std::vector<double> calculate_spectral_flux_custom(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas, Isotope isotope) {
  std::vector<double> results, errors;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_1d);
  gsl_function f;
  f.function = integrand;
  solar_model_integration_params_custom p = { 0.0, &s, isotope };
  f.params = &p;

  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
    double integral, error;
    p.erg = *erg;
    gsl_integration_qag(&f, s.get_r_lo(), s.get_r_hi(), int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, int_method_1d, w, &integral, &error);
    results.push_back(distance_factor*integral);
    errors.push_back(distance_factor*error);
  }

  gsl_integration_workspace_free(w);

  std::vector<std::vector<double>> buffer = { ergs, results, errors };
  std::string comment = standard_header(&s);
  comment += "Spectral flux over full solar volume.\nColumns: energy values [keV] | axion flux [cm^-2 s^-1 keV^-1] | axion flux error estimate [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return results;
}

std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux_custom(ergs, s, &integrand_weightedCompton, saveas);
}

std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux_custom(ergs, s, &integrand_all_ff,saveas);
}

std::vector<double> calculate_spectral_flux_opacity_element(std::vector<double> ergs, SolarModel &s, std::string element, std::string saveas) {
  Isotope isotope (element);
  return calculate_spectral_flux_custom(ergs, s, &integrand_opacity_element, saveas, isotope);
}

// Additional integration routines for integrating the content of a file
double flux_integrand_from_file(double erg, void * params) {
  OneDInterpolator * interp = (OneDInterpolator *)params;
  return interp->interpolate(erg);
}

double integrated_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file, bool includes_electron_interactions) {
  double result, error;

  OneDInterpolator spectral_flux_interpolator (spectral_flux_file);
  if ( (erg_min < spectral_flux_interpolator.lower()) || (erg_max > spectral_flux_interpolator.upper()) ) {
    erg_min = std::max(erg_min, spectral_flux_interpolator.lower());
    erg_max = std::min(erg_max, spectral_flux_interpolator.upper());
    std::cout << "WARNING! Integration boundaries for 'integrated_flux_from_file' are incompatible with the min/max available energy in the file "+spectral_flux_file
              << ". Setting integration region to overlap: [" << erg_min << ", " << erg_max << "]." << std::endl;
  }

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_file);
  gsl_function f;
  f.function = &flux_integrand_from_file;
  f.params = &spectral_flux_interpolator;

  if (includes_electron_interactions) {
    std::vector<double> relevant_peaks = get_relevant_peaks(erg_min, erg_max);
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), 10.0*int_abs_prec_file, 10.0*int_rel_prec_file, int_space_size_file, w, &result, &error);
  } else {
    gsl_integration_qag(&f, erg_min, erg_max, int_abs_prec_file, int_rel_prec_file, int_space_size_file, int_method_file, w, &result, &error);
  }

  gsl_integration_workspace_free(w);

  return result;
}
