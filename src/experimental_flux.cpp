#include "experimental_flux.hpp"

double erg_integrand_from_file(double erg, void * params) {
  struct exp_flux_params_file * p = (struct exp_flux_params_file *)params;
  double mass = p->mass;
  double length = (p->length)/eVm;
  //double norm_factor = (p->spectral_flux->interpolate(ref_erg_value))*(24.0*60.0*60.0*(p->eff_exp->interpolate(ref_erg_value)))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*m*m/ref_erg_value)/pi));

  double argument = 0.25*1.0e-3*length*mass*mass/erg;
  double sincsq = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = 24.0*60.0*60.0*(p->eff_exp->interpolate(erg));
  // N.B. Here we assume axion is massless in stellar interior:
  double exp_flux = p->spectral_flux->interpolate(erg);

  //return exposure*exp_flux*sincsq/norm_factor;
  return exposure*exp_flux*sincsq;
}

std::vector<double> axion_photon_counts(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;
  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  //const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double conversion_prob_factor = gsl_pow_2(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm));
  //const double prefactor = distance_factor*conversion_prob_factor;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  double length = (setup->length)/eVm;
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  OneDInterpolator spectral_flux (spectral_flux_file);
  //double norm_factor = (spectral_flux.interpolate(ref_erg_value))*(24.0*60.0*60.0*eff_exp.interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*mass*mass/ref_erg_value)/pi));

  double gagg_result, gagg_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p { mass, setup->length, &eff_exp, &spectral_flux };
  gsl_function f;
  f.function = &erg_integrand_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    // TODO: Should set abs prec. threshold to ~ 0.001 counts? Would need correct units for energy integrand.
    //       Massive downside: only valid for given gagg... cannot simply rescale results. Computational cost not worth it?!
    gsl_integration_qag(&f, erg_lo, erg_hi, ergint_from_file_abs_prec, ergint_from_file_rel_prec, int_space_size, ergint_from_file_method, w, &gagg_result, &gagg_error);
    //double res = gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*prefactor*norm_factor*gagg_result;
    double counts = gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gagg_result;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };
  gsl_integration_workspace_free (w);

  return result;
}

double erg_integrand(double erg, void * params)
{
  struct erg_integration_params * p3 = (struct erg_integration_params *)params;
  SolarModel *s = p3->s;
  double r_min = s->r_lo, r_max = std::min(p3->r_max, s->r_hi);
  double m = p3->mass;
  double length = (p3->length)/eVm;

  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*24.0*60.0*60.0*(p3->eff_exp->interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*m*m/ref_erg_value)/pi));

  double argument = 0.25*1.0e-3*length*m*m/erg;
  double sincsq = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = 24.0*60.0*60.0*(p3->eff_exp->interpolate(erg));

  struct solar_disc_integration_params  p2 { erg, 0.0, r_max, s, p3->integrand, p3->w1 };

  gsl_function f2;
  f2.function = &rad_integrand;
  f2.params = &p2;

  double spectral_flux, spectral_flux_error;
  gsl_integration_qag (&f2, r_min, r_max, 0.1*int_abs_prec , 0.1*int_rel_prec , int_space_size, int_method_1, p3->w2, &spectral_flux, &spectral_flux_error);

  //std::cout << "erg = " << erg << ", integral 2 = " << 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3 << std::endl;

  return 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3;
}

std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;

  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double conversion_prob_factor = gsl_pow_2(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm));
  const double factor = distance_factor*conversion_prob_factor;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  double length = setup->length;
  double norm_factor1 = s->Gamma_P_Primakoff(ref_erg_value, s->r_lo);
  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*24.0*60.0*60.0*(eff_exp.interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*mass*mass/ref_erg_value)/pi));

  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff;

  erg_integration_params p3 = { mass, length, setup->r_max, &eff_exp, s, integrand, w1, w2 };
  gsl_function f3;
  f3.function = &erg_integrand;
  f3.params = &p3;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    double gagg_result, gagg_error;
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag (&f3, erg_lo, erg_hi, int_abs_prec, int_rel_prec, int_space_size, int_method_1, w3, &gagg_result, &gagg_error);
    double counts = factor*norm_factor1*norm_factor3*gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*gagg_result;
    //std::cout << "integral 3 = " << gagg_result << std::endl;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };

  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}

std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;

  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  //const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double conversion_prob_factor = gsl_pow_2(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm));
  //const double prefactor = distance_factor*conversion_prob_factor;
  const double all_peaks [32] = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  double length = (setup->length)/eVm;
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  OneDInterpolator spectral_flux (spectral_flux_file);
  //double norm_factor = (spectral_flux.interpolate(ref_erg_value))*(24.0*60.0*60.0*eff_exp.interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*mass*mass/ref_erg_value)/pi));

  double gaee_result, gaee_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p { mass, setup->length, &eff_exp, &spectral_flux };
  gsl_function f;
  f.function = &erg_integrand_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    std::vector<double> relevant_peaks;
    relevant_peaks.push_back(erg_lo);
    for (int i = 0; i < 32; i++) { if ( (erg_lo < all_peaks[i]) && (all_peaks[i] < erg_hi) ) { relevant_peaks.push_back(all_peaks[i]); }; };
    relevant_peaks.push_back(erg_hi);
    // TODO: Should set abs prec. threshold to ~ 0.001 counts? Would need correct units for energy integrand.
    //       Massive downside: only valid for given gagg... cannot simply rescale results. Computational cost not worth it?!
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), ergint_from_file_abs_prec, ergint_from_file_rel_prec, int_space_size, w, &gaee_result, &gaee_error);
    //double res = gsl_pow_2((gaee/1.0e-13)*(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*prefactor*norm_factor*gaee_result;
    double counts = gsl_pow_2((gaee/1.0e-13)*(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gaee_result;
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };
  gsl_integration_workspace_free (w);

  return result;
}

std::vector<double> axion_electron_counts_full (double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;

  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double conversion_prob_factor = gsl_pow_2(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm));
  const double factor = distance_factor*conversion_prob_factor;
  const double all_peaks [32] = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  double length = setup->length;

  static double norm_factor1 = s->Gamma_P_all_electron(ref_erg_value, s->r_lo);
  static double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*(24.0*60.0*60.0*eff_exp.interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*mass*mass/ref_erg_value)/pi));

  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_all_electron;

  erg_integration_params p3 = { mass, length, setup->r_max, &eff_exp, s, integrand, w1, w2 };
  gsl_function f3;
  f3.function = &erg_integrand;
  f3.params = &p3;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    double gaee_result, gaee_error;
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    std::vector<double> relevant_peaks;
    relevant_peaks.push_back(erg_lo);
    for (int i = 0; i < 32; i++) { if ( (erg_lo < all_peaks[i]) && (all_peaks[i] < erg_hi) ) { relevant_peaks.push_back(all_peaks[i]); }; };
    relevant_peaks.push_back(erg_hi);
    //gsl_integration_qag (&f3, erg_lo, erg_hi, int_abs_prec, int_rel_prec, int_space_size, gagg_method, w3, &gagg_result, &gagg_error);
    gsl_integration_qagp(&f3, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec, int_rel_prec, int_space_size, w3, &gaee_result, &gaee_error);
    double counts = factor*norm_factor1*norm_factor3*gsl_pow_2((gagg/1.0e-10)*(gaee/1.0e-13)*(setup->b_field/9.0)*(setup->length/9.26))*gaee_result;
    //std::cout << "integral 3 = " << gagg_result << std::endl;
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };

  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}
