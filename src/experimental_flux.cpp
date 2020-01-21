#include "experimental_flux.hpp"

std::vector<double> axion_photon_counts (double mass, double gagg, exp_setup *setup, std::string spectral_flux_file, bool save_output, std::string output_path) {
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
    gsl_integration_qag (&f, erg_lo, erg_hi, ergint_abs_prec, ergint_rel_prec, int_space_size, gagg_method, v, &gagg_result, &gagg_error);
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(prefactor_gagg*gagg_result));
    result.push_back(gsl_pow_2(gagg*1.0e19)*prefactor_gagg*gagg_result);
  };
  gsl_integration_workspace_free (v);

  return result;
};

std::vector<double> axion_photon_counts_full (double mass, double gagg, exp_setup *setup, SolarModel *s, bool save_output, std::string output_path) {
  std::vector<double> result;
  // prefactor_gagg = (keV/eV)^6 * (1 cm^2/eVcm^2) * (1 day/eVs) * (10^10 cm/eVcm) * (10^-19 eV^-1)^4 * ((9.26 m/eVm) * (9.0 T/(T/eV^2) ))^2 / (128 pi^3)
  //const double prefactor_gagg = 29302.30262;
  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  OneDInterpolator eff_exp (setup->eff_exposure_file);

  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  exp_flux_params p3 =  { mass, setup->length, setup->r_max, 0.0, 0.0, &eff_exp, s, w1, w2 };
  gsl_function f3;
  f3.function = &erg_integrand_gagg;
  f3.params = &p3;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    double gagg_result, gagg_error;
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag (&f3, erg_lo, erg_hi, ergint_abs_prec, ergint_rel_prec, int_space_size, gagg_method, w3, &gagg_result, &gagg_error);
    double counts = gsl_pow_2((gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*gagg_result;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };

  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
};
