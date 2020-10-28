// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "spectral_flux.hpp"

/////////////////////////////////////////////////////
//  Integration routines for the solar axion flux  //
/////////////////////////////////////////////////////

// N.B. The solar axion fluxes will be calculated in units of [axions cm^-2 s^-1 keV^-1]

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

  gsl_integration_qag(p2->f, p2->s->get_r_lo(), p2->s->get_r_hi(), int_abs_prec_1d, int_rel_prec_1d, int_space_size_1d, int_method_1d, p2->w, &result, &error);

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

std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  std::vector<double> results;
  std::ofstream output;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_1d);
  gsl_function f;
  f.function = rho_integrand_1d;
  solar_model_integration_parameters_1d p { 0.0, &s, integrand, &f, w };
  f.params = &p;

  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) { results.push_back( distance_factor*erg_integrand_1d(*erg, &p) ); }

  gsl_integration_workspace_free(w);

  std::vector<std::vector<double> > buffer = { ergs, results };
  std::string comment = "Spectral flux over full solar volume by " LIBRARY_NAME ".\nColumns: energy values [keV] | axion flux [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return results;
}

std::vector<std::vector<double> > calculate_total_flux_solar_disc_at_fixed_radii(double erg_lo, double erg_hi, std::vector<double> radii, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas) {
  std::vector<double> results, errors;
  std::vector<double> valid_radii = s.get_supported_radii(radii);
  double r_min = valid_radii.front();
  double r_max = valid_radii.back();
  int n_r_vals = valid_radii.size();
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
      p.r_1 = valid_radii[i-1];
    } else {
      p.r_1 = s.get_r_lo();
    }

    p.r_2 = valid_radii[i];
    if (p.r_2 > p.r_1) {
      //gsl_integration_qagiu(&f, 0.0, 10.0*int_abs_prec_2d, 10.0*int_rel_prec_2d, int_space_size_2d, w, &integral, &error); // Alternative integration from 0 -> infinity; too slow.
      gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(),10.0*int_abs_prec_2d, 10.0*int_rel_prec_2d, int_space_size_2d, w, &integral, &error);
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

  std::vector<std::vector<double>> buffer = { valid_radii, results, errors };
  std::string comment = "Total spectral flux for a given radius by " LIBRARY_NAME ".\nColumns: Radius on solar disc [R_sol], Axion flux [cm^-2 s^-1 keV^-1] |  Axion flux error estimate [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}

std::vector<std::vector<double> > calculate_spectral_flux_solar_disc_at_fixed_radii(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  std::vector<double> all_ergs, all_radii, results;

  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);
  gsl_integration_cquad_workspace * w2 = gsl_integration_cquad_workspace_alloc(int_space_size_2d_cquad);

  std::vector<double> valid_radii = s.get_supported_radii(radii);
  double r_min = valid_radii.front();
  double r_max = valid_radii.back();
  int n_r_vals = valid_radii.size();

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
    p.r_2 = valid_radii[i];
    for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
      all_radii.push_back(p.r_2);
      all_ergs.push_back(*erg);
      results.push_back( distance_factor*erg_integrand_2d(*erg, &p) );
    }
  }

  gsl_integration_cquad_workspace_free(w1);
  gsl_integration_cquad_workspace_free(w2);

  std::vector<std::vector<double> > buffer = { all_radii, all_ergs, results };
  std::string comment = "Spectral flux over full solar disc at fixed radius by " LIBRARY_NAME ".\nColumns: Radius on solar disc [R_sol] | Energy [keV] | Axion flux [cm^-2 s^-1 keV^-1]";
  save_to_file(saveas, buffer, comment);

  return buffer;
}

std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs, double r_max, SolarModel &s, double (SolarModel::*integrand)(double, double), std::string saveas, Isotope isotope) {
  // Check if r_max >= max. available radius and switch to faster 1D integration if that's the case
  if (r_max >= s.get_r_hi()) {
    return calculate_spectral_flux(ergs, s, integrand, saveas, isotope);
  } else {
    std::vector<double> radii = { r_max };
    std::vector<std::vector<double> > result = calculate_spectral_flux_solar_disc_at_fixed_radii(ergs, radii, s, integrand, saveas, isotope);
    return result[0];
  }
}

std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux(ergs, s, &SolarModel::Gamma_P_Primakoff, saveas);
}

std::vector<std::vector<double> > calculate_spectral_flux_Primakoff(std::vector<double> ergs, std::vector<double> radii, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux_solar_disc_at_fixed_radii(ergs, radii, s, &SolarModel::Gamma_P_Primakoff, saveas);
}

std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas) {
  return calculate_spectral_flux_solar_disc(ergs, r_max, s, &SolarModel::Gamma_P_Primakoff, saveas);
}

std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas) {
  return calculate_spectral_flux(ergs, s, &SolarModel::Gamma_P_Compton, saveas);
}

std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux(ergs, s, &SolarModel::Gamma_P_all_electron, saveas);
}
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas) {
  return calculate_spectral_flux_solar_disc(ergs, r_max, s, &SolarModel::Gamma_P_all_electron, saveas);
}

std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas) {
  return calculate_spectral_flux(ergs, s, &SolarModel::Gamma_P_opacity, saveas);
}


// Here, we also define some simple, custom integration routines similar to the ones defined above
// Weighted Compton contribution
double integrand_weightedCompton(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  if (erg == 0) {return 0;}
  SolarModel* sol = (p->sol);
  double u = erg/(sol->temperature_in_keV(r));

  return 0.5*gsl_pow_2(r*erg/pi)*0.5*(1.0 - 1.0/gsl_expm1(u))*(sol->Gamma_P_Compton(erg, r));
}

// Includes FF flux and ee contribution as in arXiv:1310.0823
double integrand_all_ff(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_ff(erg, r) + sol->Gamma_P_ee(erg, r));
}

// Calculate the flux from opacity for one element only
double integrand_opacity_element(double r, void * params) {
  struct solar_model_integration_params_custom * p = (struct solar_model_integration_params_custom *)params;
  double erg = (p->erg);
  std::string el_name = (p->isotope).get_element_name();
  SolarModel* sol = (p->sol);

  return 0.5*gsl_pow_2(r*erg/pi)*(sol->Gamma_P_opacity(erg, r, el_name));
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
  std::string comment = "Spectral flux over full solar volume by " LIBRARY_NAME ".\nColumns: energy values [keV] | axion flux [axions / cm^2 s keV] | axion flux error estimate [axions / cm^2 s keV]";
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
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec_file, int_rel_prec_file, int_space_size_file, w, &result, &error);
  } else {
    gsl_integration_qag(&f, erg_min, erg_max, int_abs_prec_file, int_rel_prec_file, int_space_size_file, int_method_file, w, &result, &error);
  }

  gsl_integration_workspace_free(w);

  return result;
}

//////////////////////////////////////////////////////////
//  More advanced methods and classes for the spectrum  //
//////////////////////////////////////////////////////////

AxionSpectrum::AxionSpectrum() {}

void AxionSpectrum::init_numbered_interp(const int index, const double* x, const double* y) {
  acc[index] = gsl_interp_accel_alloc();
  spline[index] = gsl_spline_alloc(gsl_interp_linear, pts);
  gsl_spline_init(spline[index], x, y, pts);
}

void AxionSpectrum::init_table_mode(std::string file, double g1, double g2) {
  mode = table;
  ASCIItableReader data (file);
  pts = data.getnrow();
  int n_cols = data.getncol();
  table_submode = n_cols-1;
  default_g1 = g1;
  default_g2 = g2;
  terminate_with_error_if((n_cols<2)||(n_cols>4), "ERROR! Your input file '"+file+"' appears to have less than 2 or more than 4 columns!");
  //data.resize(n_cols);
  // TODO: Interpolation of the specturm is better for log10(flux), but need a safe and numerically sound way to deal with 0 w/o causing problems in integration.
  //for (int i=0; i<n_cols-1; ++i) { data[i] = std::move(data[i]); }
  //data[n_cols-1] = std::move(safe_log10(data[n_cols-1]));
  if (table_submode < 3) {
    acc.resize(n_cols-1);
    spline.resize(n_cols-1);
    init_numbered_interp(0, &data[0][0], &data[1][0]);
    if (table_submode == 2) { init_numbered_interp(1, &data[0][0], &data[2][0]); }
  } else {
    int n_grids = 1 + (g2>0);
    acc_2d.resize(n_grids);
    spline_2d.resize(n_grids);
    // Make sure that entries are properly sorted
    std::vector<double> x_vec = data[0];
    std::sort(x_vec.begin(), x_vec.end());
    x_vec.erase(std::unique(x_vec.begin(), x_vec.end()), x_vec.end());
    int nx = x_vec.size();
    std::vector<double> y_vec = data[1];
    std::sort(y_vec.begin(), y_vec.end());
    y_vec.erase(std::unique(y_vec.begin(), y_vec.end()), y_vec.end());
    int ny = y_vec.size();
    terminate_with_error_if(nx*ny != pts, "ERROR! Number of data points ("+std::to_string(pts)+") in '"+file+"' inconsistent with number of unique 'x' and 'y' values ("+std::to_string(nx)+" and "+std::to_string(ny)+")! Check the formatting.");

    const double* x = &x_vec[0];
    const double* y = &y_vec[0];
    // Allocate memory for "z" values array in gsl format
    std::vector<double*> z (n_grids);
    for (int i=0; i<n_grids; ++i) {
      z[i] = (double*) malloc(nx * ny * sizeof(double));
      spline_2d[i] = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
      acc_2d[i].first = gsl_interp_accel_alloc();
      acc_2d[i].second = gsl_interp_accel_alloc();
    }

    // Determine first and last "x" and "y" values and grid step size.
    double x_lo = x_vec.front();
    double x_up = x_vec.back();
    double y_lo = y_vec.front();
    double y_up = y_vec.back();
    double x_delta = (x_up-x_lo) / (nx-1);
    double y_delta = (y_up-y_lo) / (ny-1);

    // Intialise grids
    for (int j=0; j<pts; ++j) {
      // Determine appropriate indices for the grid points.
      double temp = (data[0][j]-x_lo) / x_delta;
      int ind_x = (int) (temp+0.5);
      temp = (data[1][j]-y_lo) / y_delta;
      int ind_y = (int) (temp+0.5);
      for (int i=0; i<n_grids; ++i) { gsl_spline2d_set(spline_2d[i], z[i], ind_x, ind_y, data[i+2][j]); }
    }
    for (int i=0; i<n_grids; ++i) { gsl_spline2d_init(spline_2d[i], x, y, z[i], nx, ny); }
  }
}

AxionSpectrum::AxionSpectrum(std::string file, double g1, double g2) { init_table_mode(file, g1, g2); }

void AxionSpectrum::init_solar_model_mode(SolarModel* sol, SolarModelMemberFn process2) {
  mode = solar_model; // Do not set default_g1, default_g2 here! This info in the solar model.
  s = sol;
  function2 = process2;
  solar_model_okay = s->is_initialised();
  terminate_with_error_if(not(solar_model_okay), "ERROR! The instance of SolarModel used to initialse an instance of AxionSpectrum has not been initialised properly.");
}

AxionSpectrum::AxionSpectrum(SolarModel* sol, SolarModelMemberFn process2) { init_solar_model_mode(sol, process2); }

void AxionSpectrum::init_analytical_mode(double norm, double g1, double a, double b) {
  mode = analytical; // Do not set default_g1 here! This info in the parameters.
  analytical_parameters = { norm, g1, a, b };
}

AxionSpectrum::AxionSpectrum(double norm, double g1, double a, double b) { init_analytical_mode(norm, g1, a, b); }

void AxionSpectrum::switch_mode(SpectrumModes new_mode) {
  if (new_mode != mode) {
    bool check = false;
    switch (new_mode) {
      case table :
        check = (table_submode >= 1) && (table_submode <= 3);
        break;
      case solar_model :
        check = solar_model_okay;
        break;
      case analytical :
        check = (analytical_parameters.size() == 4);
        break;
      default :
        std::cout << "WARNING! You tried to change to an unknown mode." << std::endl;
    }
    if (check) {
      mode = new_mode;
    } else {
      std::cout << "WARNING! Changing mode of AxionSpectrum is not permitted; the class instance has not been initialised for the new mode." << std::endl;
    }
  }
}

AxionSpectrum::~AxionSpectrum() {
  for (auto s : spline) { gsl_spline_free(s); }
  for (auto s : spline_2d) { gsl_spline2d_free(s); }
  for (auto a : acc) { gsl_interp_accel_free(a); }
  for (auto a : acc_2d) { gsl_interp_accel_free(a.first); gsl_interp_accel_free(a.second); }
}

std::tuple<int, double, double> AxionSpectrum::get_table_parameters() { return std::make_tuple(table_submode, default_g1, default_g2); }

std::vector<double> AxionSpectrum::get_analytical_parameters() { return analytical_parameters; }

SpectrumModes AxionSpectrum::get_class_mode() {
  std::cout << "INFO. The current mode is '";
  switch (mode) {
    case table : std::cout << "table"; break;
    case analytical : std::cout << "analytical"; break;
    case solar_model : std::cout << "solar_model"; break;
    case undefined : std::cout << "undefined"; break;
    default : std::cout << "n/a'. This is a bug, please report it"; break;
  }
  std::cout << "'." << std::endl;
  return mode;
}

std::vector<std::vector<double>> AxionSpectrum::axion_flux(std::vector<double> ergs, std::vector<double> radii, double g1, double g2) {
  std::vector<std::vector<double>> result;
  static bool g2_warning_issued = false;
  static bool r_warning_issued = false;
  for (auto r = radii.begin(); r != radii.end(); r++) {
    //std::vector<double> g1_flux
    //std::vector<double> g2_flux;
    std::vector<double> total_flux;
    switch(mode) {
      case table :
        if ((table_submode == 1) && (g2 > 0) && not(g2_warning_issued)) { std::cout << "WARNING! You don't have a second spectrum loaded! The value of g2 > 0 has no effect." << std::endl; g2_warning_issued = true; }
        if ((table_submode < 3) && (*r < 1) && not(r_warning_issued)) { std::cout << "WARNING! You didn't provide data for radii! The value of r < 1 has no effect." << std::endl; r_warning_issued = true; }
        for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
          double ref_g1_flux = 0;
          double ref_g2_flux = 0;
          if (table_submode < 4) {
            ref_g1_flux = gsl_pow_2(g1/default_g1) * gsl_spline_eval(spline[0], *erg, acc[0]);
            if (table_submode > 1) { ref_g2_flux = gsl_pow_2(g2/default_g2) * gsl_spline_eval(spline[1], *erg, acc[1]); }
          } else {
            ref_g1_flux = gsl_pow_2(g1/default_g1) * gsl_spline2d_eval(spline_2d[0], *r, *erg, acc_2d[0].first, acc_2d[0].second);
            if (default_g2 > 0) { ref_g2_flux = gsl_pow_2(g2/default_g2) * gsl_spline2d_eval(spline_2d[1], *r, *erg, acc_2d[1].first, acc_2d[1].second); }
          }
          //g1_flux.push_back(ref_g1_flux);
          //g2_flux.push_back(ref_g2_flux);
          total_flux.push_back(ref_g1_flux+ref_g2_flux);
        }
        break;
      case solar_model :
      {
        static double ref_g1 = s->get_gagg_ref_value_in_inverse_GeV();
        static double ref_g2 = s->get_gaee_ref_value();
        //std::vector<double> ref_g1_fluxes = s->calculate_spectral_flux_any(ergs, function1, *r);
        std::vector<double> ref_g1_fluxes = calculate_spectral_flux_solar_disc(ergs, *r, *s, function1);
        //std::vector<double> ref_g2_fluxes = s->calculate_spectral_flux_any(ergs, function2, *r);
        std::vector<double> ref_g2_fluxes = calculate_spectral_flux_solar_disc(ergs, *r, *s, function2);
        for (int i=0; i<ergs.size(); ++i) {
          double ref_g1_flux = gsl_pow_2(g1/ref_g1) * ref_g1_fluxes[i];
          double ref_g2_flux = gsl_pow_2(g2/ref_g2) * ref_g2_fluxes[i];
          //g1_flux.push_back(ref_g1_flux);
          //g2_flux.push_back(ref_g2_flux);
          total_flux.push_back(ref_g1_flux+ref_g2_flux);
        }
        break;
      }
      case analytical :
        for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
          double ref_g1_flux = analytical_parameters[0] * gsl_pow_2(g1/analytical_parameters[1]) * pow(*erg,analytical_parameters[2]) * exp(-analytical_parameters[3]*(*erg));
          //g1_flux.push_back(ref_g1_flux);
          total_flux.push_back(ref_g1_flux);
        }
        break;
      case undefined :
        terminate_with_error("ERROR! The AxionSpectrum class has not been initialised properly.");
        break;
      default :
        terminate_with_error("FATAL ERROR! The AxionSpectrum mode was somehow set to an unknown value; this is a bug, please report it.");
    }
    result.push_back(total_flux);
  }
  return result;
}

double AxionSpectrum::axion_flux(double erg, double g1) { return axion_flux(erg, g1, 0); }

double AxionSpectrum::axion_flux(double erg, double g1, double g2) { return axion_flux(erg, 1, g1, 0); }

double AxionSpectrum::axion_flux(double erg, double r, double g1, double g2) {
  std::vector<double> ergs = { erg };
  std::vector<double> result = axion_flux(ergs, r, g1, g2);
  return result[0];
}

std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double g1) { return axion_flux(ergs, g1, 0); }

std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double g1, double g2) { return axion_flux(ergs, 1, g1, 0); }

std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double r, double g1, double g2) {
  std::vector<double> radii = { r };
  std::vector<std::vector<double>> result = axion_flux(ergs, radii, g1, g2);
  return result[0];
}

// Monte Carlo-related functions
AxionMCGenerator1D::AxionMCGenerator1D() {
  inv_cdf_acc = gsl_interp_accel_alloc();
  inv_cdf = gsl_spline_alloc(gsl_interp_linear, 2);
}

void AxionMCGenerator1D::init_from_spectral_data(std::vector<double> ergs, std::vector<double> flux) {
  int pts = ergs.size();
  inv_cdf_data_erg = std::vector<double> (pts);
  inv_cdf_data_x = std::vector<double> (pts);

  double norm = 0.0;
  inv_cdf_data_erg[0] = ergs[0];
  inv_cdf_data_x[0] = norm;
  for(int i=1; i<pts; ++i) {
    // Trapozoidal rule integration.
    norm += 0.5 * (ergs[i] - ergs[i-1]) * (flux[i] + flux[i-1]);
    inv_cdf_data_erg[i] = ergs[i];
    inv_cdf_data_x[i] = norm;
  }

  integrated_norm = norm;
  for (int i=0; i<pts; ++i) { inv_cdf_data_x[i] = inv_cdf_data_x[i]/integrated_norm; }

  init_inv_cdf_interpolator();
}

AxionMCGenerator1D::AxionMCGenerator1D(std::vector<double> ergs, std::vector<double> flux) { init_from_spectral_data(ergs, flux); }

AxionMCGenerator1D::AxionMCGenerator1D(std::string file, bool is_already_inv_cdf_file) {
  if (is_already_inv_cdf_file) {
    ASCIItableReader inv_cdf_data (file);
    inv_cdf_data_x = inv_cdf_data[0];
    inv_cdf_data_erg = inv_cdf_data[1];

    terminate_with_error_if(inv_cdf_data_x.end()[-2] > 1.0, "ERROR! Sanity check for MC generator failed! The second to last entry of your inverse CDF is greater than 1.");
    integrated_norm = 1.0;

    init_inv_cdf_interpolator();
  } else {
    ASCIItableReader spectrum_data (file);
    init_from_spectral_data(spectrum_data[0], spectrum_data[1]);
  }
}

void AxionMCGenerator1D::change_parameters(double erg_min, double erg_max, double erg_delta) {
  omega_min = std::min(erg_min,erg_max);
  omega_max = std::max(erg_min,erg_max);
  omega_delta = erg_delta;
}

std::vector<double> AxionMCGenerator1D::generate_ergs() {
  std::vector<double> result;
  int n_omega_vals = int((omega_max-omega_min)/omega_delta);
  for (int i=0; i<n_omega_vals; ++i) { result.push_back(omega_min + i*omega_delta); }
  return result;
}

void AxionMCGenerator1D::init_with_local_spectrum(double g1, double g2, double r) {

  inv_cdf_data_erg = generate_ergs();
  int n_omega_vals = inv_cdf_data_erg.size();
  inv_cdf_data_x = std::vector<double> (n_omega_vals);
  default_r = r;

  //std::vector<double> flux1 = sol->calculate_spectral_flux_Primakoff(inv_cdf_data_erg, r);
  //std::vector<double> flux2 = sol->calculate_spectral_flux_all_electron(inv_cdf_data_erg, r);

  std::vector<double> flux1 = sp.axion_flux(inv_cdf_data_erg, g1, g2, r);
  std::vector<double> flux2 = sp.axion_flux(inv_cdf_data_erg, g1, g2, r);

  double norm = 0.0;
  inv_cdf_data_x[0] = norm;
  for(int i=1; i<n_omega_vals; ++i) {
    // Trapozoidal rule integration.
    norm += 0.5 * omega_delta * (flux1[i] + flux1[i-1]);
    norm += 0.5 * omega_delta * (flux2[i] + flux2[i-1]);
    inv_cdf_data_x[i] = norm;
  }

  integrated_norm = norm;
  for (int i=0; i<n_omega_vals; ++i) { inv_cdf_data_x[i] = inv_cdf_data_x[i]/integrated_norm; }

  init_inv_cdf_interpolator();
  full_mc_generator_ready = true;
}

AxionMCGenerator1D::AxionMCGenerator1D(SolarModel* sol, double g1, double g2, double r) {
  sp = AxionSpectrum(sol);
  init_with_local_spectrum(g1, g2, r);
}

AxionMCGenerator1D::AxionMCGenerator1D(double a, double b) {
  analytical_parameters = { a, b };
  analytical_mc_generator_ready = true;
}

AxionMCGenerator1D::AxionMCGenerator1D(AxionSpectrum* spectrum, double g1, double g2, double r) {
  SpectrumModes mode = spectrum->get_class_mode();
  switch(mode) {
    case table :
    {
      auto table_params = spectrum->get_table_parameters();
      int table_submode = std::get<0>(table_params);
      double g1 = std::get<1>(table_params);
      double g2 = std::get<2>(table_params);
      inv_cdf_data_erg = generate_ergs();
      std::vector<double> flux = spectrum->axion_flux(inv_cdf_data_erg, r, g1, g2);
      init_from_spectral_data(inv_cdf_data_erg, flux);
      if (table_submode > 2) {
        sp = *spectrum;
        full_mc_generator_ready = true;
      }
      break;
    }
    case analytical :
    {
      std::vector<double> p = spectrum->get_analytical_parameters();
      analytical_parameters = { p[2], p[3] };
      analytical_mc_generator_ready = true;
      break;
    }
    case solar_model :
    {
      sp = *spectrum;
      init_with_local_spectrum(g1, g2, r);
      break;
    }
    default : std::cout << "WARNING! The mode of AxionMCGenerator1D, used in the construction of AxionSpectrum, is inappropriate.";
  }
}

void AxionMCGenerator1D::init_inv_cdf_interpolator() {
  int pts = inv_cdf_data_x.size();
  terminate_with_error_if(pts < 2, "ERROR! You tried to initialise the MC generator with less than two data points.");

  const double* prob = &inv_cdf_data_x[0];
  const double* erg = &inv_cdf_data_erg[0];

  inv_cdf_acc = gsl_interp_accel_alloc();
  inv_cdf = gsl_spline_alloc(gsl_interp_linear, pts);
  gsl_spline_init(inv_cdf, prob, erg, pts);

  simple_mc_generator_ready = true;
}

void AxionMCGenerator1D::save_inv_cdf_to_file(std::string inv_cdf_file) {
  const std::string comment = "Inverse CDF obtained by " LIBRARY_NAME ".\nColumns: Random variable | Inverse CDF of energy distribution";
  std::vector<std::vector<double>> buffer;
  buffer.push_back(inv_cdf_data_x);
  buffer.push_back(inv_cdf_data_erg);
  save_to_file(inv_cdf_file, buffer, comment);
}

double AxionMCGenerator1D::evaluate_inv_cdf(double x) {
  double result;
  if (analytical_mc_generator_ready) {
    static double ap1 = 1.0 + analytical_parameters.first;
    static double invb = 1.0 / analytical_parameters.second;
    result = gsl_cdf_gamma_Pinv(x, ap1, invb);
  } else if (simple_mc_generator_ready) {
    result = gsl_spline_eval(inv_cdf, x, inv_cdf_acc);
  } else {
    terminate_with_error("ERROR! The (simple version of the) MC generator has not been initialised properly!");
  }
  return result;
}

std::vector<double> AxionMCGenerator1D::draw_axion_energies(int n) {
  std::vector<double> result;

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng (seed);
  std::uniform_real_distribution<double> unif(0, 1);

  //std::cout << "Test U[0,1]: " << unif(rng) << std::endl;
  for (int i=0; i<n; ++i) { result.push_back(evaluate_inv_cdf( unif(rng) )); }

  return result;
}

std::vector<double> AxionMCGenerator1D::draw_axion_energies(int n, double g1, double g2) {
  terminate_with_error_if(not(full_mc_generator_ready), "ERROR! The (full version of the) MC generator has not been initialised properly!");

  std::vector<double> ergs = generate_ergs();
  std::vector<double> flux = sp.axion_flux(ergs, default_r, g1, g2);
  AxionMCGenerator1D mc (ergs,flux);

  return mc.draw_axion_energies(n);
}

double AxionMCGenerator1D::get_norm() {
  terminate_with_error_if(not(simple_mc_generator_ready), "ERROR! The (simple version of the) MC generator has not been initialised properly!");
  return integrated_norm;
}

AxionMCGenerator1D::~AxionMCGenerator1D() {
  gsl_spline_free(inv_cdf);
  gsl_interp_accel_free(inv_cdf_acc);
}

AxionMCGenerator2D::AxionMCGenerator2D() { }

AxionMCGenerator2D::AxionMCGenerator2D(std::string file, bool is_already_inv_cdf_file) {
  if (is_already_inv_cdf_file) {
    //inv_cdf_grid = TwoDInterpolator(file);
    //init_inv_cdf_interpolator();
  } else {
    terminate_with_error("ERROR! 2D interpolation + integration of files is currently not supported.");
    //ASCIItableReader spectrum_data (file);
    //init_from_spectral_data(spectrum_data);
  }
}

AxionMCGenerator2D::AxionMCGenerator2D(std::vector<std::vector<double> > data, bool is_already_inv_cdf_file) {
  if (is_already_inv_cdf_file) {
    //inv_cdf_grid = TwoDInterpolator(data);
    //init_inv_cdf_interpolator();
  } else {
    terminate_with_error("ERROR! 2D interpolation + integration of spectral data is currently not supported.");
    //init_from_spectral_data(data);
  }
}

AxionMCGenerator2D::AxionMCGenerator2D(SolarModel &sol, std::vector<double> ergs, std::vector<double> rads, double gaee, std::string save_fluxes_prefix) {
  std::vector<std::vector<double> > total_fluxes, fluxes, total_gaee_fluxes, gaee_fluxes;
  std::vector<double> cdf_rads;
  std::string filename = "";
  double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff;
  if (save_fluxes_prefix != "") { filename = save_fluxes_prefix+"_total_Primakoff.dat"; }
  total_fluxes = calculate_total_flux_solar_disc_at_fixed_radii(0.0, 50.0, rads, sol, integrand, filename);
  int n_rad_entries = total_fluxes[0].size();
  if (save_fluxes_prefix != "") { filename = save_fluxes_prefix+"_Primakoff.dat"; }
  fluxes = calculate_spectral_flux_solar_disc_at_fixed_radii(ergs, rads, sol, integrand, filename);
  int n_entries = fluxes[0].size();
  if (gaee > 0) {
    integrand = &SolarModel::Gamma_P_all_electron;
    if (save_fluxes_prefix != "") { filename = save_fluxes_prefix+"_total_all_electron.dat"; }
    total_gaee_fluxes = calculate_total_flux_solar_disc_at_fixed_radii(0.0, 50.0, rads, sol, integrand, filename);
    if (save_fluxes_prefix != "") { filename = save_fluxes_prefix+"_all_electron.dat"; }
    gaee_fluxes = calculate_spectral_flux_solar_disc_at_fixed_radii(ergs, rads, sol, integrand, filename);
    for (int i = 0; i < n_rad_entries; ++i) { total_fluxes[1][i] += gsl_pow_2(gaee/1.0e-13)*total_gaee_fluxes[1][i]; }
    for (int i = 0; i < n_entries; ++i) { fluxes[2][i] += gsl_pow_2(gaee/1.0e-13)*gaee_fluxes[2][i]; }
  }

  double norm = 0;
  int k = 0;
  std::vector<double> norms (n_rad_entries, 0.0);
  inv_cdf_ergs.resize(n_rad_entries);
  for (int i = 0; i < n_rad_entries; i++) {
    if (i > 0) { norm += 0.5*(total_fluxes[0][i] - total_fluxes[0][i-1])*(total_fluxes[1][i] + total_fluxes[1][i-1]); }
    //double omega_min = sqrt(sol.omega_pl_squared(total_fluxes[0][i]));
    //double omega_max = sol.temperature_in_keV(total_fluxes[0][i]);
    k += 1;
    if (i < n_rad_entries-1) {
      while (fluxes[0][k] < total_fluxes[0][i+1]) {
        norms[i] += 0.5*(fluxes[1][k] - fluxes[1][k-1])*(fluxes[2][k] + fluxes[2][k-1]);
        k += 1;
      }
    } else {
      while (k < fluxes[0].size()) {
        norms[i] += 0.5*(fluxes[1][k] - fluxes[1][k-1])*(fluxes[2][k] + fluxes[2][k-1]);
        k += 1;
      }
    }
    std::cout << "Norm: " << norms[i] << std::endl;
  }

  double integral1 = 0;
  k = 0;
  for (int i = 0; i < n_rad_entries; i++) {
    if (i > 0) {
      integral1 += 0.5*(total_fluxes[0][i] - total_fluxes[0][i-1])*(total_fluxes[1][i] + total_fluxes[1][i-1]);
      cdf_rads.push_back(integral1/norm);
    } else {
      cdf_rads.push_back(0);
    }
    std::cout << total_fluxes[0][i] << " | " << cdf_rads[i] << std::endl;
    double integral2 = 0;
    std::vector<double> temp_ergs;
    std::vector<double> temp_cdf_ergs;
     temp_ergs.push_back(fluxes[1][k]);
     temp_cdf_ergs.push_back(0);
     k += 1;
    if (i < n_rad_entries-1) {
      while (fluxes[0][k] < total_fluxes[0][i+1]) {
        integral2 += 0.5*(fluxes[1][k] - fluxes[1][k-1])*(fluxes[2][k] + fluxes[2][k-1]);
        temp_ergs.push_back(fluxes[1][k]);
        temp_cdf_ergs.push_back(integral2/norms[i]);
        std::cout << total_fluxes[0][i] << " vs " << fluxes[0][k] << " | " << fluxes[1][k] << " | " << integral2/norms[i] << std::endl;
        k += 1;
      }
    } else {
      while (k < n_entries) {
        integral2 += 0.5*(fluxes[1][k] - fluxes[1][k-1])*(fluxes[2][k] + fluxes[2][k-1]);
        temp_ergs.push_back(fluxes[1][k]);
        temp_cdf_ergs.push_back(integral2/norms[i]);
        std::cout << total_fluxes[0][i] << " vs " << fluxes[0][k] << " | " << fluxes[1][k] << " | " << integral2/norms[i] << std::endl;
        k += 1;
      }
    }

    std::vector<std::vector<double> > buffer1 = { temp_cdf_ergs, temp_ergs };
    OneDInterpolator interp (buffer1);
    inv_cdf_ergs[i] = std::move(interp);
  }

  std::vector<std::vector<double> > buffer2 = { cdf_rads, total_fluxes[0] };
  inv_cdf_rad = OneDInterpolator(buffer2);

  //std::vector<std::vector<double> > buffer2 = { fluxes[0], cdf_ergs, fluxes[1] };
  //inv_cdf_grid = TwoDInterpolator(buffer2);

  radii = total_fluxes[0];

  mc_generator_ready = true;
}

//AxionMCGenerator2D::init_from_spectral_data(std::vector<std::vector<double>> data, std::vector<std::vector<double>> data);

void AxionMCGenerator2D::init_inv_cdf_interpolator() {
  //energies = inv_cdf_grid.get_unique_y_vals();
  //radii = inv_cdf_grid.get_unique_y_vals();
  mc_generator_ready = true;
}

void AxionMCGenerator2D::save_inv_cdf_to_file(std::string inv_cdf_file) {
  std::string comment = "Inverse CDF obtained by " LIBRARY_NAME ".\nColumns: Radius on solar disc [Rsol] | Random variable | Inverse CDF of energy distribution, given radius.";
  std::vector<std::vector<double> > buffer (3);
  for (int i=0; i < radii.size(); i++) {
    std::cout << "Loop Printig data... " << i << " / " << radii.size() << std::endl;
    std::vector<std::vector<double> > data = inv_cdf_ergs[i].get_data();
    int n_ergs = data[0].size();
    std::vector<double> rads (n_ergs, radii[i]);
    buffer[0].insert(buffer[0].end(), rads.begin(), rads.end());
    buffer[1].insert(buffer[1].end(), data[0].begin(), data[0].end());
    buffer[2].insert(buffer[2].end(), data[1].begin(), data[1].end());
  }
  //save_to_file(inv_cdf_file, inv_cdf_grid.get_data(), comment);
  save_to_file(inv_cdf_file+".dat", buffer, comment);
  comment = "Inverse CDF obtained by " LIBRARY_NAME ".\nColumns: Random variable | Inverse CDF of radius distribution";
  save_to_file(inv_cdf_file+"_aux.dat", inv_cdf_rad.get_data(), comment);
}

double AxionMCGenerator2D::evaluate_inv_cdf_rad(double x) {
  terminate_with_error_if(not(mc_generator_ready), "ERROR! The 2D MC generator has not been initialised properly!");
  return inv_cdf_rad.interpolate(x);
}

double AxionMCGenerator2D::evaluate_inv_cdf_erg_given_rad(double x, double rad) {
  terminate_with_error_if(not(mc_generator_ready), "ERROR! The 2D MC generator has not been initialised properly!");
  double result;
  auto it = lower_bound(radii.begin(), radii.end(), rad);
  size_t index = std::distance(radii.begin(), it);
  if ( (it == radii.begin()) || (it == radii.end()) ) {
    result = inv_cdf_ergs[index].interpolate(x);
  } else {
    double r1 = *(it-1);
    double r2 = *it;
    double e1 = inv_cdf_ergs[index-1].interpolate(x);
    double e2 = inv_cdf_ergs[index].interpolate(x);
    result = e1 + (rad - r1)*(e2 - e1)/(r2 - r1);
  }
  return result;
  //return inv_cdf_grid.interpolate(rad, x);
}

std::vector<double> AxionMCGenerator2D::draw_axion_radii(int n) {
  std::vector<double> result;

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng (seed);
  std::uniform_real_distribution<double> unif(0, 1);

  for (int i=0; i<n; ++i) { result.push_back(inv_cdf_rad.interpolate( unif(rng) )); }

  return result;
}

std::vector<double> AxionMCGenerator2D::draw_axion_energies_given_radii(std::vector<double> radii) {
  std::vector<double> result;

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng (seed);
  std::uniform_real_distribution<double> unif(0, 1);

  for (auto r = radii.begin(); r != radii.end(); ++r) {
    result.push_back(evaluate_inv_cdf_erg_given_rad( unif(rng), *r ));
  }

  return result;
}

std::vector<double> AxionMCGenerator2D::draw_axion_energies(int n) {
  std::vector<double> radii = draw_axion_radii(n);
  return draw_axion_energies_given_radii(radii);
}

AxionMCGenerator2D::~AxionMCGenerator2D() { }
