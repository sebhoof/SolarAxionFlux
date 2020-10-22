// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "experimental_flux.hpp"


///////////////////////////////////////////////////////
//  Experimental axion-photon conversion and spectra //
///////////////////////////////////////////////////////

// Conversion probability correction for massive axions (mass in eV)
double conversion_prob_correction(double mass, double erg, double length) {
  if (mass > 0) {
    double argument = 0.25*1.0e-3*(length/eVm)*mass*mass/erg;
    return gsl_pow_2(gsl_sf_sinc(argument/pi));
  }
  return 1.0;
}

// Effective exposures (in seconds x cm^2) for various experiments.
double eff_exposure_cast2007(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2007_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_a(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_A_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_b(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_B_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_c(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_C_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_d(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_D_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_e(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_E_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_f(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_F_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_g(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_G_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_h(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_H_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_i(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_I_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_j(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_J_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_k(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_K_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

double eff_exposure_cast2017_l(double erg) {
  static OneDInterpolator eff_exp ("data/exposures/CAST2017_L_EffectiveExposure.dat");
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

// Baseline IAXO exposure from [arXiv:1904.09155]
double eff_exposure_iaxo(double erg) {
  const double eff = 0.7*0.8;
  const double time = 3.0 * 0.5 * 365.0*24.0*60.0*60.0;
  const double area = 2.3e4;
  const double eff_exp = eff*time*area;
  return eff_exp;
}

// babyIAXO exposure from [arXiv:1904.09155]
double eff_exposure_babyiaxo(double erg) {
  const double eff = 0.35*0.7;
  const double time = 1.5 * 0.5 * 365.0*24.0*60.0*60.0;
  const double area = 0.77e4;
  const double eff_exp = eff*time*area;
  return eff_exp;
}

// IAXO+ exposure from [arXiv:1904.09155]
double eff_exposure_iaxoplus(double erg) {
  const double eff = 0.7*0.8;
  const double time = 5.0 * 0.5 * 365.0*24.0*60.0*60.0;
  const double area = 3.9e4;
  const double eff_exp = eff*time*area;
  return eff_exp;
}

double eff_exposure(double erg, std::string dataset) {
  experiment exp = experiment_name.at(dataset);
  return eff_exposure(erg, exp);
}

double eff_exposure(double erg, experiment dataset) {
  double res = 0.0;
  switch (dataset) {
    case(CAST2007):
      res = eff_exposure_cast2007(erg); break;
    case(CAST2017_A):
      res = eff_exposure_cast2017_a(erg); break;
    case(CAST2017_B):
      res = eff_exposure_cast2017_b(erg); break;
    case(CAST2017_C):
      res = eff_exposure_cast2017_c(erg); break;
    case(CAST2017_D):
      res = eff_exposure_cast2017_d(erg); break;
    case(CAST2017_E):
      res = eff_exposure_cast2017_e(erg); break;
    case(CAST2017_F):
      res = eff_exposure_cast2017_f(erg); break;
    case(CAST2017_G):
      res = eff_exposure_cast2017_g(erg); break;
    case(CAST2017_H):
      res = eff_exposure_cast2017_h(erg); break;
    case(CAST2017_I):
      res = eff_exposure_cast2017_i(erg); break;
    case(CAST2017_J):
      res = eff_exposure_cast2017_j(erg); break;
    case(CAST2017_K):
      res = eff_exposure_cast2017_k(erg); break;
    case(CAST2017_L):
      res = eff_exposure_cast2017_l(erg); break;
    case(IAXO):
      res = eff_exposure_iaxo(erg); break;
    case(BABYIAXO):
      res = eff_exposure_babyiaxo(erg); break;
    case(IAXOPLUS):
      res = eff_exposure_iaxoplus(erg); break;
    default:
      terminate_with_error("ERROR! Data set not known!");
  }
  return res;
}


//////////////////////////////////////////////////////////////////////////////////////
//  Modified and extended integration routines (possibly include energy dispersion) //
//////////////////////////////////////////////////////////////////////////////////////

// Various wrapper and helper functions

// Additional integration routines for integrating the content of a file (with effective exposure; no convolution)
double exp_flux_integrand_from_file(double erg, void * params) {
  struct exp_flux_from_file_integration_parameters * p = (struct exp_flux_from_file_integration_parameters *)params;

  double sincsq = conversion_prob_correction(p->mass, erg, p->length);
  double exposure = eff_exposure(erg, p->dataset);
  // N.B. Here we assume axion is massless in stellar interior:
  double exp_flux = p->spectral_flux->interpolate(erg);

  return exposure*exp_flux*sincsq;
}

double simple_convolution_kernel(double erg, void * params) {
  struct simple_convolution_params * p = (struct simple_convolution_params *)params;
  double two_sigma2 = 2.0*gsl_pow_2(p->sigma);
  return (p->spectral_flux->interpolate(erg))*exp(-gsl_pow_2(p->erg0 - erg)/two_sigma2)/sqrt(two_sigma2*pi);
}

double convolution_kernel(double erg, void * params) {
  struct convolution_params * p = (struct convolution_params *)params;
  double sincsq = conversion_prob_correction(p->p->mass, erg, p->p->length);
  double exposure = eff_exposure(erg, p->p->dataset);
  double exp_flux = p->p->spectral_flux->interpolate(erg);
  double two_sigma2 = 2.0*gsl_pow_2(p->p->sigma);
  return exposure*exp_flux*sincsq * exp(-gsl_pow_2(p->erg0 - erg)/two_sigma2)/sqrt(two_sigma2*pi);
}

// Additional integration routines for integrating the content of a file (with effective exposure and convolution)
double convolved_exp_flux_integrand_from_file(double erg, void * params) {
  double result, error;
  struct exp_flux_from_file_integration_parameters * p = (struct exp_flux_from_file_integration_parameters *)params;
  struct convolution_params q = { erg, p };
  std::vector<double> relevant_peaks = get_relevant_peaks(p->support[0], p->support[1]);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size_file);

  gsl_function f;
  f.function = &convolution_kernel;
  f.params = &q;

  gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w, &result, &error);

  gsl_integration_workspace_free(w);

  return result;
}

// TODO OUTDATED! Basically the alternative to 'from file integration'?
double erg_integrand(double erg, void * params) {
  struct erg_integration_params * p3 = (struct erg_integration_params *)params;
  SolarModel *s = p3->s;
  double r_min = s->get_r_lo(), r_max = std::min(p3->r_max, s->get_r_hi());

  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*(eff_exposure(ref_erg_value, p3->dataset))*conversion_prob_correction(p3->mass, ref_erg_value, p3->length);

  double sincsq = conversion_prob_correction(p3->mass, erg, p3->length);
  double exposure = eff_exposure(erg, p3->dataset);

  gsl_function f2;
  f2.function = &rad_integrand_2d;
  struct solar_model_integration_parameters_2d p2 { erg, 0.0, 0.0, 0.0, s, p3->integrand, &f2, p3->w1, &f2, p3->w1 };

  f2.params = &p2;

  double spectral_flux, spectral_flux_error;
  gsl_integration_qag (&f2, r_min, r_max, 0.1*int_abs_prec_file , 0.1*int_rel_prec_file , int_space_size_file, int_method_file, p3->w2, &spectral_flux, &spectral_flux_error);

  return 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3;
}


// TODO: Here actual integration routines from file
std::vector<double> convolved_spectrum_from_file(std::vector<double> ergs, double support[2], double resolution, std::string filename) {
  std::vector<double> result;
  OneDInterpolator spectrum (filename);

  double flux, flux_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size_file);

  struct simple_convolution_params p = { resolution, 0, &spectrum };

  gsl_function f;
  f.function = &simple_convolution_kernel;
  f.params = &p;

  std::vector<double> relevant_peaks = get_relevant_peaks(support[0], support[1]);
  for (auto erg0 = ergs.begin(); erg0 != ergs.end(); ++erg0) {
    p.erg0 = *erg0;
    std::cout << *erg0 << " " << flux << std::endl;
    //gsl_integration_qags(&f, support[0], support[1], int_abs_prec, int_rel_prec, int_space_size_file, w, &flux, &flux_error);
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w, &flux, &flux_error);
    //gsl_integration_qag (&f, support[0], support[1], int_abs_prec, int_rel_prec, int_space_size_file, int_method_1, w, &flux, &flux_error);
    result.push_back(flux);
  }
  gsl_integration_workspace_free(w);

  return result;
}

// Functions to calculate the counts in all bins of a helioscope experiment
std::vector<double> axion_photon_counts_from_file(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;
  OneDInterpolator spectral_flux;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  double support [2] = { bin_lo, bin_hi };

  // Do convolution if necessary.
  double erg_resolution = setup->erg_resolution;
  if (erg_resolution > 0) {
    ASCIItableReader temp (spectral_flux_file);
    std::vector<double> ergs = temp[0];
    std::vector<double> flux = convolved_spectrum_from_file(ergs, support, erg_resolution, spectral_flux_file);
    spectral_flux = OneDInterpolator(ergs, flux);
    save_to_file(spectral_flux_file+"_convolved", {ergs, flux}, "");
  } else {
    spectral_flux = OneDInterpolator(spectral_flux_file);
  }

  double gagg_result, gagg_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size_file);
  exp_flux_from_file_integration_parameters p { mass, setup->length, setup->dataset, &spectral_flux, {support[0], support[1]}, setup->erg_resolution };
  gsl_function f;
  f.function = &exp_flux_integrand_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag(&f, erg_lo, erg_hi, int_abs_prec_file, int_rel_prec_file, int_space_size_file, int_method_file, w, &gagg_result, &gagg_error);
    double counts = gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gagg_result;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  }
  gsl_integration_workspace_free (w);

  return result;
}

// Integration routines for full experimental counts.
// CAVE: Can be very slow!
// TODO check for overlap/outdated?
std::vector<double> axion_photon_counts_full(double mass, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;

  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double factor = distance_factor*conversion_prob_factor;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double norm_factor1 = s->Gamma_P_Primakoff(ref_erg_value, s->get_r_lo());
  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*eff_exposure(ref_erg_value, setup->dataset)*conversion_prob_correction(mass, ref_erg_value, setup->length);

  //gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size_file);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size_file);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size_file);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff;

  erg_integration_params p3 = { mass, setup->length, setup->r_max, setup->dataset, s, integrand, w1, w2 };
  gsl_function f3;
  f3.function = &erg_integrand;
  f3.params = &p3;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    double gagg_result, gagg_error;
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag (&f3, erg_lo, erg_hi, int_abs_prec_file, int_rel_prec_file, int_space_size_file, int_method_file, w3, &gagg_result, &gagg_error);
    double counts = factor*norm_factor1*norm_factor3*gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*gagg_result;
    //std::cout << "integral 3 = " << gagg_result << std::endl;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  }

  //gsl_integration_workspace_free (w1);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}

std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;
  const double all_peaks [32] = { 0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                  5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785 };
  static OneDInterpolator spectral_flux (spectral_flux_file);

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);

  double gaee_result, gaee_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size_file);
  exp_flux_from_file_integration_parameters p { mass, setup->length, setup->dataset, &spectral_flux, setup->erg_resolution };
  gsl_function f;
  f.function = &exp_flux_integrand_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    std::vector<double> relevant_peaks;
    relevant_peaks.push_back(erg_lo);
    for (int i = 0; i < 32; i++) { if ( (erg_lo < all_peaks[i]) && (all_peaks[i] < erg_hi) ) { relevant_peaks.push_back(all_peaks[i]); } }
    relevant_peaks.push_back(erg_hi);
    // TODO: Should set abs prec. threshold to ~ 0.001 counts? Would need correct units for energy integrand.
    //       Massive downside: only valid for given gagg... cannot simply rescale results. Computational cost not worth it?!
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w, &gaee_result, &gaee_error);
    double counts = gsl_pow_2((gaee/1.0e-13)*(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gaee_result;
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  }
  gsl_integration_workspace_free (w);

  return result;
}

// Integration routines for full experimental counts.
// CAVE: Can be very slow!
// TODO check for overlap/outdated?
std::vector<double> axion_electron_counts_full(double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;
  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double factor = distance_factor*conversion_prob_factor;
  const double all_peaks [32] = { 0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                  5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785 };

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);

  double norm_factor1 = s->Gamma_P_all_electron(ref_erg_value, s->get_r_lo());
  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*eff_exposure(ref_erg_value,setup->dataset)*conversion_prob_correction(mass, ref_erg_value, setup->length);

  //gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size_file);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size_file);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size_file);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_all_electron;

  erg_integration_params p3 = { mass, setup->length, setup->r_max, setup->dataset, s, integrand, w1, w2 };
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
    for (int i = 0; i < 32; i++) { if ( (erg_lo < all_peaks[i]) && (all_peaks[i] < erg_hi) ) { relevant_peaks.push_back(all_peaks[i]); } }
    relevant_peaks.push_back(erg_hi);
    //gsl_integration_qag (&f3, erg_lo, erg_hi, int_abs_prec, int_rel_prec, int_space_size_file, gagg_method, w3, &gagg_result, &gagg_error);
    gsl_integration_qagp(&f3, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w3, &gaee_result, &gaee_error);
    double counts = factor*norm_factor1*norm_factor3*gsl_pow_2((gagg/1.0e-10)*(gaee/1.0e-13)*(setup->b_field/9.0)*(setup->length/9.26))*gaee_result;
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  }

  //gsl_integration_workspace_free (w1);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
//  Routines to compute 'reference count files' for the various experiments //
//////////////////////////////////////////////////////////////////////////////

// Return relative counts at reference values of the coupling.
std::vector<std::vector<double>> axion_reference_counts_from_file(exp_setup *setup, std::vector<double> masses, std::string spectral_flux_file_gagg, std::string spectral_flux_file_gaee, std::string saveas, bool save_convolved_spectra) {
  std::vector<std::vector<double>> result;
  std::vector<std::vector<double>> convolved_spectra_gagg, convolved_spectra_gaee;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  double overall_factor = gsl_pow_2((setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor;
  double erg_resolution = setup->erg_resolution;

  /*
  OneDInterpolator spectral_flux_gagg, spectral_flux_gaee;
  // TODO: Do convolution if necessary -> Actually, this should be done only later!
  if (erg_resolution > 0) {

    ASCIItableReader temp (spectral_flux_file_gagg);
    std::vector<double> ergs = temp[0];
    std::vector<double> flux_gagg = convolved_spectrum_from_file(ergs, support, erg_resolution, spectral_flux_file_gagg);
    spectral_flux_gagg = OneDInterpolator(ergs, flux_gagg);
    save_to_file(spectral_flux_file_gagg+"_conv", {ergs, flux_gagg}, "");

    if (spectral_flux_file_gaee != "") {
      temp = ASCIItableReader (spectral_flux_file_gaee);
      ergs = temp[0];
      std::vector<double> flux_gaee = convolved_spectrum_from_file(ergs, support, erg_resolution, spectral_flux_file_gaee);
      spectral_flux_gaee = OneDInterpolator(ergs, flux_gaee);
      save_to_file(spectral_flux_file_gaee+"_conv", {ergs, flux_gaee}, "");
    }

  } else {
    spectral_flux_gagg = OneDInterpolator(spectral_flux_file_gagg);
    if (spectral_flux_file_gaee != "") { spectral_flux_gaee = OneDInterpolator(spectral_flux_file_gaee); }
  }
  */

  double gagg_result, gagg_error, gaee_result, gaee_error;
  OneDInterpolator spectral_flux_gagg (spectral_flux_file_gagg);
  OneDInterpolator spectral_flux_gaee;
  if (spectral_flux_file_gaee != "") { spectral_flux_gaee = OneDInterpolator(spectral_flux_file_gaee); }

  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size_file);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size_file);
  exp_flux_from_file_integration_parameters p1 { 0, setup->length, setup->dataset, &spectral_flux_gagg, {bin_lo, bin_hi}, setup->erg_resolution };
  exp_flux_from_file_integration_parameters p2 { 0, setup->length, setup->dataset, &spectral_flux_gaee, {bin_lo, bin_hi}, setup->erg_resolution };
  gsl_function f1, f2;
  f1.params = &p1;
  f2.params = &p2;
  if (erg_resolution > 0) {
    f1.function = &convolved_exp_flux_integrand_from_file;
    f2.function = &convolved_exp_flux_integrand_from_file;
  } else {
    f1.function = &exp_flux_integrand_from_file;
    f2.function = &exp_flux_integrand_from_file;
  }

  std::vector<std::vector<double>> relevant_peaks;
  if (spectral_flux_file_gaee != "") {
    for (int bin = 0; bin < n_bins; ++bin) {
      double erg_lo = bin_lo + bin*bin_delta;
      relevant_peaks.push_back(get_relevant_peaks(erg_lo, erg_lo+bin_delta));
    }
  }

  std::vector<double> gagg_ergs, gaee_ergs;
  if ((erg_resolution > 0) && save_convolved_spectra) {
    ASCIItableReader temp1 (spectral_flux_file_gagg);
    gagg_ergs = temp1[0];
    if (spectral_flux_file_gaee != "") {
      ASCIItableReader temp2 (spectral_flux_file_gaee);
      gaee_ergs = temp1[0];
    }
  }

  std::vector<double> expanded_masses, bin_centres, results_gagg, results_gaee;
  std::vector<double> convolved_spectra_masses_gagg, convolved_spectra_masses_gaee, convolved_spectra_energies_gagg, convolved_spectra_energies_gaee, convolved_spectra_results_gagg, convolved_spectra_results_gaee;
  for (auto mass = masses.begin(); mass != masses.end(); mass++) {
    p1.mass = *mass;
    p2.mass = *mass;
    if ((erg_resolution > 0) && save_convolved_spectra) {
      for (auto erg = gagg_ergs.begin(); erg != gagg_ergs.end(); ++erg) {
        convolved_spectra_masses_gagg.push_back(*mass);
        convolved_spectra_energies_gagg.push_back(*erg);
        convolved_spectra_results_gagg.push_back(GSL_FN_EVAL(&f1,*erg));
      }
      if (spectral_flux_file_gaee != "") {
        for (auto erg = gaee_ergs.begin(); erg != gaee_ergs.end(); ++erg) {
          convolved_spectra_masses_gaee.push_back(*mass);
          convolved_spectra_energies_gaee.push_back(*erg);
          convolved_spectra_results_gaee.push_back(GSL_FN_EVAL(&f2,*erg));
        }
      }
    }
    for (int bin = 0; bin < n_bins; ++bin) {
      expanded_masses.push_back(*mass);
      double erg_lo = bin_lo + bin*bin_delta;
      double erg_hi = erg_lo + bin_delta;
      bin_centres.push_back(0.5*(erg_lo + erg_hi));
      // TODO: Improve results with QAWO adaptive integration for oscillatory functions?! factor out sin^2()?
      gsl_integration_qag(&f1, erg_lo, erg_hi, int_abs_prec_file, int_rel_prec_file, int_space_size_file, int_method_file, w1, &gagg_result, &gagg_error);
      results_gagg.push_back(overall_factor*gagg_result);
      if (spectral_flux_file_gaee != "") {
        gsl_integration_qagp(&f2, &relevant_peaks[bin][0], relevant_peaks[bin].size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w2, &gaee_result, &gaee_error);
        results_gaee.push_back(overall_factor*gaee_result);
      }
    }
  }

  std::string header = "Reference counts for g_agamma = 10^-10 1/GeV and g_ae = 10^-13\nColumns: Axion mass [eV] | Energy bin centre [keV] | Counts from Primakoff";
  result.push_back(expanded_masses);
  result.push_back(bin_centres);
  result.push_back(results_gagg);
  if (spectral_flux_file_gaee != "") { result.push_back(results_gaee); header += " | Counts from axion-electron"; }
  save_to_file(saveas, result, header);

  header = "Reference flux (convolved) for g_agamma = 10^-10 1/GeV\nColumns: Axion mass [eV] | Energy [keV] | Primakoff flux [s^-1 cm^-1 keV^-1]";
  convolved_spectra_gagg.push_back(convolved_spectra_masses_gagg);
  convolved_spectra_gagg.push_back(convolved_spectra_energies_gagg);
  convolved_spectra_gagg.push_back(convolved_spectra_results_gagg);
  save_to_file(spectral_flux_file_gagg+"_convolved", convolved_spectra_gagg, header);

  if (spectral_flux_file_gaee != "") {
    convolved_spectra_gaee.push_back(convolved_spectra_masses_gaee);
    convolved_spectra_gaee.push_back(convolved_spectra_energies_gaee);
    convolved_spectra_gaee.push_back(convolved_spectra_results_gaee);
    save_to_file(spectral_flux_file_gaee+"_convolved", convolved_spectra_gaee, header);
  }

  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);

  return result;
}

std::vector<double> counts_prediciton_from_file(double mass, double gagg, std::string reference_counts_file, double gaee) {
  std::vector<double> result;

  // Only setup the interpolators once while the function is being used...
  static std::string reference_counts_file_bak = "";
  static int n_cols;
  static int n_bins;
  static int n_masses;
  static double min_m = 1.0e-4;
  static double lgm0 = -4.0; // Axions with m = 10^-4 eV are effectivly massless
  static std::vector<double> log_masses;
  static std::vector<std::vector<double>> ref_counts_gagg;
  static std::vector<OneDInterpolator> interp_ref_counts_gagg;
  static std::vector<std::vector<double>> ref_counts_gaee;
  static std::vector<OneDInterpolator> interp_ref_counts_gaee;

  // ... unless the user decides to change the reference counts file.
  bool new_reference_counts_file = (reference_counts_file != reference_counts_file_bak);
  if (new_reference_counts_file) {
    reference_counts_file_bak = reference_counts_file;
    ASCIItableReader data (reference_counts_file);
    n_cols = data.getncol();
    int n_rows = data.getnrow();

    // Make sure that the table is processed correctly with any formatting
    // Extract unique mass values from the file
    std::vector<double> m = data[0];
    sort(m.begin(), m.end());
    m.erase(unique(m.begin(), m.end()), m.end());
    n_masses = m.size();
    log_masses = std::vector<double> (n_masses);
    min_m = m[0];
    if ((min_m > 0) && (min_m < 1.0e-4)) { lgm0 = log10(min_m); }
    if (n_masses > 1) { if ((min_m = 0) && (m[1] < 1.0e-4)) { lgm0 = log10(m[1]) - 100.0; } }
    for (int k=0; k<n_masses; ++k) { log_masses[k] = safe_log10(m[k],lgm0); }
    // Extract unique energy bin values from the file
    std::vector<double> bins = data[1];
    sort(bins.begin(), bins.end());
    bins.erase(unique(bins.begin(), bins.end()), bins.end());
    n_bins = bins.size();

    // Resize the vectors according to the inputs of the file.
    interp_ref_counts_gagg = std::vector<OneDInterpolator> (n_bins);
    ref_counts_gagg.resize(n_bins, std::vector<double> (n_masses));
    if (n_cols > 3) {
      interp_ref_counts_gaee = std::vector<OneDInterpolator> (n_bins);
      ref_counts_gaee.resize(n_bins, std::vector<double> (n_masses));
    }

    // Assign the data from the file to the appropriate places.
    for (int i=0; i<n_rows; ++i) {
      auto itj = std::find(bins.begin(), bins.end(), data[1][i]);
      int j = std::distance(bins.begin(), itj);
      auto itk = std::find(m.begin(), m.end(), data[0][i]);
      int k = std::distance(m.begin(), itk);
      ref_counts_gagg[j][k] = data[2][i];
      if (n_cols > 3) { ref_counts_gaee[j][k] = data[3][i]; }
    }

    if (n_masses > 1) {
      for (int j=0; j<n_bins; ++j) {
        //OneDInterpolator temp_gagg (log_masses, ref_counts_gagg[j]);
        //interp_ref_counts_gagg[j] = std::move(temp_gagg);
        interp_ref_counts_gagg[j] = OneDInterpolator(log_masses, ref_counts_gagg[j]);
        if (n_cols > 3) {
          OneDInterpolator temp_gaee (log_masses, ref_counts_gaee[j]);
          interp_ref_counts_gaee[j] = std::move(temp_gaee);
        }
      }
    }
  }

  // Reference values are gagg = 10^-10/GeV and gaee = 10^-13.
  double gagg_rel_sq = gagg*gagg/1.0e-20;
  double gaee_rel_sq = gaee*gaee/1.0e-26;
  double lgm = safe_log10(mass,lgm0);
  if (n_masses > 1) {
    for (int i = 0; i < n_bins; ++i) {
      double temp = gagg_rel_sq*interp_ref_counts_gagg[i].interpolate(lgm);
      if (n_cols > 3) { temp += gaee_rel_sq*interp_ref_counts_gaee[i].interpolate(lgm); }
      result.push_back(gagg_rel_sq*temp);
    }
  } else {
    terminate_with_error_if(mass != min_m, "ERROR! Your reference counts file is only valid for an axion mass of m = "+std::to_string(min_m)+" eV!");
    for (int i = 0; i < n_bins; ++i) {
      double temp = gagg_rel_sq*ref_counts_gagg[i][0];
      if (n_cols > 3) { temp += gaee_rel_sq*ref_counts_gaee[i][0]; }
      result.push_back(gagg_rel_sq*temp);
    }
  }

  return result;
}
