#include "experimental_flux.hpp"

double erg_integrand_gagg_from_file(double erg, void * params) {
  const double eVm = gev2cm*1.0e7;
  struct exp_flux_params_file * p = (struct exp_flux_params_file *)params;
  double m = p->mass;
  double length = (p->length)/eVm;

  double argument = 0.25*1.0e-3*length*m*m/erg;
  double temp = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = p->eff_exp->interpolate(erg);
  // N.B. Here we assume axion is massless in stellar interior:
  double exp_flux = p->spectral_flux->interpolate(erg);

  return temp*exposure*exp_flux;
}

std::vector<double> axion_photon_counts(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file, bool save_output, std::string output_path) {
  std::vector<double> result;
  // prefactor_gagg = (keV/eV)^6 * (1 cm^2/eVcm^2) * (1 day/eVs) * (10^10 cm/eVcm) * (10^-19 eV^-1)^4 * ((9.26 m/eVm) * (9.0 T/(T/eV^2) ))^2 / (128 pi^3)
  const double prefactor_gagg = 29302.30262;
  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  OneDInterpolator spectral_flux (spectral_flux_file);

  double gagg_result, gagg_error;
  gsl_integration_workspace * v = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p { mass, setup->length, &eff_exp, &spectral_flux };
  gsl_function f;
  f.function = &erg_integrand_gagg_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    // TODO: Should set abs prec. threshold to ~ 0.001 counts? Would need correct units for energy integrand.
    //       Massive downside: only valid for given gagg... cannot simply rescale results. Computational cost not worth it?!
    gsl_integration_qag(&f, erg_lo, erg_hi, int_abs_prec, int_rel_prec, int_space_size, gagg_method, v, &gagg_result, &gagg_error);
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(prefactor_gagg*gagg_result));
    result.push_back(gsl_pow_2(gagg*1.0e19)*prefactor_gagg*gagg_result);
  };
  gsl_integration_workspace_free (v);

  return result;
}

double erg_integrand_gagg(double erg, void * params)
{
  const double eVm = 1.0e7*gev2cm;

  struct exp_flux_params * p3 = (struct exp_flux_params *)params;
  p3->erg = erg;
  SolarModel *s = p3->s;
  double r_min = s->r_lo, r_max = std::min(p3->r_max, s->r_hi);
  double m = p3->mass;
  double length = (p3->length)/eVm;
  double (SolarModel::*func_ptr)(double, double) = &SolarModel::Gamma_P_Primakoff;

  static double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*24.0*60.0*60.0*(p3->eff_exp->interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*m*m/ref_erg_value)/pi));

  double argument = 0.25*1.0e-3*length*m*m/erg;
  double sincsq = gsl_pow_2(gsl_sf_sinc(argument/pi));
  double exposure = 24.0*60.0*60.0*(p3->eff_exp->interpolate(erg));

  struct solar_disc_integration_params  p2 { erg, 0.0, r_max, s, func_ptr, p3->w1 };

  gsl_function f2;
  f2.function = &rad_integrand;
  f2.params = &p2;

  double spectral_flux, spectral_flux_error;
  gsl_integration_qag (&f2, r_min, r_max, 0.1*int_abs_prec , 0.1*int_rel_prec , int_space_size, gagg_method, p3->w2, &spectral_flux, &spectral_flux_error);

  //std::cout << "erg = " << erg << ", integral 2 = " << 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3 << std::endl;

  return 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3;
}

std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s, bool save_output, std::string output_path) {
  std::vector<double> result;

  const double eVm = 1.0e7*gev2cm;
  const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3;
  const double distance_factor = 1.0e-4*pow(radius_sol/(1.0e-2*keV2cm),3) / (pow(distance_sol,2) * (1.0e6*hbar));
  const double conversion_prob_factor = pow(0.5*1.0e-19*(9.0/eV2T)*(9.26/eVm),2);
  const double factor = distance_factor*conversion_prob_factor;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  OneDInterpolator eff_exp (setup->eff_exposure_file);
  double length = setup->length;
  static double norm_factor1 = s->Gamma_P_Primakoff(ref_erg_value, s->r_lo);
  static double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*24.0*60.0*60.0*(eff_exp.interpolate(ref_erg_value))*gsl_pow_2(gsl_sf_sinc((0.25*1.0e-3*length*mass*mass/ref_erg_value)/pi));

  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  exp_flux_params p3 = { mass, length, setup->r_max, 0.0, 0.0, &eff_exp, s, w1, w2 };
  gsl_function f3;
  f3.function = &erg_integrand_gagg;
  f3.params = &p3;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    double gagg_result, gagg_error;
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag (&f3, erg_lo, erg_hi, int_abs_prec, int_rel_prec, int_space_size, gagg_method, w3, &gagg_result, &gagg_error);
    double counts = factor*norm_factor1*norm_factor3*gsl_pow_2((gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*gagg_result;
    //std::cout << "integral 3 = " << gagg_result << std::endl;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };

  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}
