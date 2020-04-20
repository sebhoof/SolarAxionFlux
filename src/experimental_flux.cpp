#include "experimental_flux.hpp"

////////////////////////////////////////
//  Experimental setup and coherence  //
////////////////////////////////////////

// Conversion probability correction for massive axions.
double conversion_prob_correction(double mass, double erg, double length) {
  double argument = 0.25*1.0e-3*(length/eVm)*mass*mass/erg;
  return gsl_pow_2(gsl_sf_sinc(argument/pi));
}

// Effective exposures (in seconds x cm) for various experiments.
double eff_exposure_cast_2007(double erg) {
  static OneDInterpolator eff_exp ("data/CAST2007_EffectiveExposure.dat");
  // Effective exposure file is in cm x days.
  return 24.0*60.0*60.0*eff_exp.interpolate(erg);
}

/////////////////////////
//  Wrapper functions  //
/////////////////////////

double erg_integrand_from_file(double erg, void * params) {
  struct exp_flux_params_file * p = (struct exp_flux_params_file *)params;

  double sincsq = conversion_prob_correction(p->mass, erg, p->length);
  double exposure = p->eff_exposure(erg);
  // N.B. Here we assume axion is massless in stellar interior:
  double exp_flux = p->spectral_flux->interpolate(erg);

  return exposure*exp_flux*sincsq;
}

double erg_integrand(double erg, void * params) {
  struct erg_integration_params * p3 = (struct erg_integration_params *)params;
  SolarModel *s = p3->s;
  double r_min = s->get_r_lo(), r_max = std::min(p3->r_max, s->get_r_hi());

  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*(p3->eff_exposure(ref_erg_value))*conversion_prob_correction(p3->mass, ref_erg_value, p3->length);

  double sincsq = conversion_prob_correction(p3->mass, erg, p3->length);
  double exposure = p3->eff_exposure(erg);

  struct solar_model_integration_parameters p2 { erg, 0.0, r_max, s, p3->integrand, p3->w1 };

  gsl_function f2;
  f2.function = &rad_integrand_2d;
  f2.params = &p2;

  double spectral_flux, spectral_flux_error;
  gsl_integration_qag (&f2, r_min, r_max, 0.1*int_abs_prec , 0.1*int_rel_prec , int_space_size, int_method_1, p3->w2, &spectral_flux, &spectral_flux_error);

  return 0.5*gsl_pow_2(erg/pi)*exposure*spectral_flux*sincsq/norm_factor3;
}

double convolution_kernel(double erg, void * params) {
  struct convolution_params * p = (struct convolution_params *)params;
  double two_sigma2 = 2.0*gsl_pow_2(p->sigma);
  return (p->spectral_flux->interpolate(erg))*exp(-gsl_pow_2(p->erg0 - erg)/two_sigma2)/sqrt(two_sigma2*pi);
}

std::vector<double> convolved_spectrum_from_file(std::vector<double> ergs, double support[2], double resolution, std::string filename) {
  std::vector<double> result;
  OneDInterpolator spectrum (filename);

  double flux, flux_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);

  struct convolution_params p = {resolution, 0, &spectrum};

  gsl_function f;
  f.function = &convolution_kernel;
  f.params = &p;

  for (auto erg0 = ergs.begin(); erg0 != ergs.end(); ++erg0) {
    p.erg0 = *erg0;
    gsl_integration_qag (&f, support[0], support[1], int_abs_prec, int_rel_prec , int_space_size, int_method_1, w, &flux, &flux_error);
    result.push_back(flux);
  };
  gsl_integration_workspace_free (w);

  return result;
}

////////////////////////////////////////////////////////////////////////////////
//  Functions to calculate the counts in all bins of a helioscope experiment  //
////////////////////////////////////////////////////////////////////////////////

std::vector<double> axion_photon_counts_from_file(double mass, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;
  OneDInterpolator spectral_flux;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);

  // Do convolution if necessary.
  double erg_resolution = setup->erg_resolution;
  if (erg_resolution > 0) {
    ASCIItableReader temp (spectral_flux_file);
    std::vector<double> ergs = temp[0];
    double support [2] = {bin_lo, bin_hi};
    std::vector<double> flux = convolved_spectrum_from_file(ergs, support, erg_resolution, spectral_flux_file);
    spectral_flux = OneDInterpolator(ergs, flux);
    save_to_file(spectral_flux_file+"_convolved", {ergs, flux}, "");
  } else {
    spectral_flux = OneDInterpolator(spectral_flux_file);
  };

  double gagg_result, gagg_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p { mass, setup->length, setup->eff_exposure, &spectral_flux };
  gsl_function f;
  f.function = &erg_integrand_from_file;
  f.params = &p;

  double erg_lo, erg_hi = bin_lo;
  for (int bin = 0; bin < n_bins; ++bin) {
    erg_lo = erg_hi;
    erg_hi += bin_delta;
    gsl_integration_qag(&f, erg_lo, erg_hi, ergint_from_file_abs_prec, ergint_from_file_rel_prec, int_space_size, ergint_from_file_method, w, &gagg_result, &gagg_error);
    double counts = gsl_pow_2(gsl_pow_2(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gagg_result;
    printf("gagg | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };
  gsl_integration_workspace_free (w);

  return result;
}

std::vector<double> axion_photon_counts_full(double mass, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;

  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double factor = distance_factor*conversion_prob_factor;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double norm_factor1 = s->Gamma_P_Primakoff(ref_erg_value, s->get_r_lo());
  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*(setup->eff_exposure(ref_erg_value))*conversion_prob_correction(mass, ref_erg_value, setup->length);

  //gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_cquad);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff;

  erg_integration_params p3 = { mass, setup->length, setup->r_max, setup->eff_exposure, s, integrand, w1, w2 };
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

  //gsl_integration_workspace_free (w1);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}

std::vector<double> axion_electron_counts(double mass, double gaee, double gagg, exp_setup *setup, std::string spectral_flux_file) {
  std::vector<double> result;
  const double all_peaks [32] = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};
  static OneDInterpolator spectral_flux (spectral_flux_file);

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);

  double gaee_result, gaee_error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p { mass, setup->length, setup->eff_exposure, &spectral_flux };
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
    double counts = gsl_pow_2((gaee/1.0e-13)*(gagg/1.0e-10)*(setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor*gaee_result;
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };
  gsl_integration_workspace_free (w);

  return result;
}

std::vector<double> axion_electron_counts_full(double mass, double gaee, double gagg, exp_setup *setup, SolarModel *s) {
  std::vector<double> result;
  const double distance_factor = 1.0e-4*gsl_pow_3(radius_sol/(1.0e-2*keV2cm)) / (gsl_pow_2(distance_sol) * (1.0e6*hbar));
  const double factor = distance_factor*conversion_prob_factor;
  const double all_peaks [32] = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);

  double norm_factor1 = s->Gamma_P_all_electron(ref_erg_value, s->get_r_lo());
  double norm_factor3 = 0.5*gsl_pow_2(ref_erg_value/pi)*(setup->eff_exposure(ref_erg_value))*conversion_prob_correction(mass, ref_erg_value, setup->length);

  //gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_cquad);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (int_space_size);

  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_all_electron;

  erg_integration_params p3 = { mass, setup->length, setup->r_max, setup->eff_exposure, s, integrand, w1, w2 };
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
    printf("gaee | % 6.4f [%3.2f, %3.2f] % 4.3e\n", log10(mass), erg_lo, erg_hi, log10(counts));
    result.push_back(counts);
  };

  //gsl_integration_workspace_free (w1);
  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_workspace_free (w2);
  gsl_integration_workspace_free (w3);

  return result;
}


// Return relative counts at reference values of the coupling.
std::vector<std::vector<double>> axion_reference_counts_from_file(exp_setup *setup, std::vector<double> masses, std::string spectral_flux_file_gagg, std::string spectral_flux_file_gaee, std::string saveas) {
  std::vector<std::vector<double>> result;

  int n_bins = setup->n_bins;
  double bin_lo = setup->bin_lo;
  double bin_delta = setup->bin_delta;
  double bin_hi = bin_lo + bin_delta*double(n_bins);
  double overall_factor = gsl_pow_2((setup->b_field/9.0)*(setup->length/9.26))*conversion_prob_factor;
  double support [2] = {bin_lo, bin_hi};

  OneDInterpolator spectral_flux_gagg, spectral_flux_gaee;
  // Do convolution if necessary.
  double erg_resolution = setup->erg_resolution;
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
    };
  } else {
    spectral_flux_gagg = OneDInterpolator(spectral_flux_file_gagg);
    if (spectral_flux_file_gaee != "") { spectral_flux_gaee = OneDInterpolator(spectral_flux_file_gaee); };
  };

  double gagg_result, gagg_error, gaee_result, gaee_error;
  gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (int_space_size);
  gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (int_space_size);
  exp_flux_params_file p1 { 0, setup->length, setup->eff_exposure, &spectral_flux_gagg };
  exp_flux_params_file p2 { 0, setup->length, setup->eff_exposure, &spectral_flux_gaee };
  gsl_function f1, f2;
  f1.function = &erg_integrand_from_file;
  f1.params = &p1;
  f2.function = &erg_integrand_from_file;
  f2.params = &p2;

  std::vector<std::vector<double>> relevant_peaks;
  if (spectral_flux_file_gaee != "") {
    for (int bin = 0; bin < n_bins; ++bin) {
      double erg_lo = bin_lo + bin*bin_delta;
      relevant_peaks.push_back(get_relevant_peaks(erg_lo, erg_lo+bin_delta));
    };
  };

  std::vector<double> expanded_masses, bin_centres, results_gagg, results_gaee;
  for (auto mass = masses.begin(); mass != masses.end(); mass++) {
    p1.mass = *mass;
    p2.mass = *mass;
    for (int bin = 0; bin < n_bins; ++bin) {
      expanded_masses.push_back(*mass);
      double erg_lo = bin_lo + bin*bin_delta;
      double erg_hi = erg_lo + bin_delta;
      bin_centres.push_back(0.5*(erg_lo + erg_hi));
      gsl_integration_qag(&f1, erg_lo, erg_hi, ergint_from_file_abs_prec, ergint_from_file_rel_prec, int_space_size, ergint_from_file_method, w1, &gagg_result, &gagg_error);
      results_gagg.push_back(overall_factor*gagg_result);
      if (spectral_flux_file_gaee != "") {
        gsl_integration_qagp(&f2, &relevant_peaks[bin][0], relevant_peaks[bin].size(), ergint_from_file_abs_prec, ergint_from_file_rel_prec, int_space_size, w2, &gaee_result, &gaee_result);
        results_gaee.push_back(overall_factor*gaee_result);
      };
    };
  };

  result.push_back(expanded_masses);
  result.push_back(bin_centres);
  result.push_back(results_gagg);
  if (spectral_flux_file_gaee != "") { result.push_back(results_gaee); };

  save_to_file(saveas, result, "Reference counts for g_agamma = 10^-11 1/GeV and g_ae = 10^-13\nColumns: Axion mass [eV] | Energy bin centre [keV] | Counts from Primakoff | Counts from axion-electron");
  gsl_integration_workspace_free (w1);
  gsl_integration_workspace_free (w2);

  return result;
}

double safe_log10(double x, double lgx0) {
  double result = lgx0;
  if (x > 0) { result = log10(x); };
  return result;
}

std::vector<double> counts_prediciton_from_file(double mass, double gagg, std::string reference_counts_file, double gaee) {
  std::vector<double> result;

  // Only setup the interpolators once while the function is being used...
  static std::string reference_counts_file_bak = "";
  static int n_cols;
  static int n_bins;
  static double lgm0;
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

    // Make sure that the table is processed correctly with any formatting.
    // Extract unique mass values from the file.
    std::vector<double> m = data[0];
    sort(m.begin(), m.end());
    m.erase(unique(m.begin(), m.end()), m.end());
    int n_masses = m.size();
    // TODO: Need checks/safeguards for n_masses = 1 and n_masses = 1 with m[0] = 0.
    log_masses = std::vector<double> (n_masses);
    lgm0 = -10; // Axions with m = 10^-10 eV are as good as massless.
    if ((m[0] > 0) && (m[0] < 1.0e-10)) { lgm0 = log10(m[0]); };
    for (int k=0; k<n_masses; ++k) { log_masses[k] = safe_log10(m[k],lgm0); };
    // Extract unique energy bin values from the file.
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
    };

    // Assign the data from the file to the appropriate places.
    for (int i=0; i<n_rows; ++i) {
      auto itj = std::find(bins.begin(), bins.end(), data[1][i]);
      int j = std::distance(bins.begin(), itj);
      auto itk = std::find(m.begin(), m.end(), data[0][i]);
      int k = std::distance(m.begin(), itk);
      ref_counts_gagg[j][k] = data[2][i];
      if (n_cols > 3) { ref_counts_gaee[j][k] = data[3][i]; };
    };

    for (int j=0; j<n_bins; ++j) {
      //OneDInterpolator temp_gagg (log_masses, ref_counts_gagg[j]);
      //interp_ref_counts_gagg[j] = std::move(temp_gagg);
      interp_ref_counts_gagg[j] = OneDInterpolator(log_masses, ref_counts_gagg[j]);
      if (n_cols > 3) {
        OneDInterpolator temp_gaee (log_masses, ref_counts_gaee[j]);
        interp_ref_counts_gaee[j] = std::move(temp_gaee);
      };
    };
  };

  // Reference values are gagg = 10^-10/GeV and gaee = 10^-13.
  double gagg_rel_sq = gagg*gagg/1.0e-20;
  double gaee_rel_sq = gaee*gaee/1.0e-26;
  double lgm = safe_log10(mass,lgm0);
  for (int i = 0; i < n_bins; ++i) {
    double temp = gagg_rel_sq*interp_ref_counts_gagg[i].interpolate(lgm);
    if (n_cols > 3) { temp += gaee_rel_sq*interp_ref_counts_gaee[i].interpolate(lgm); };
    result.push_back(gagg_rel_sq*temp);
  };

  return result;
}
