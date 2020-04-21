#include "spectral_flux.hpp"

AxionSpectrum::AxionSpectrum() {};

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
  terminate_with_error_if((n_cols<2)||(n_cols>4), "ERROR! Your infput file '"+file+"' appears to have less than 2 or more than 4 columns!");
  //data.resize(n_cols);
  // TODO: Interpolation of the specturm is better for log10(flux), but need a safe and numerically sound way to deal with 0 w/o causing problems in integration.
  //for (int i=0; i<n_cols-1; ++i) { data[i] = std::move(data[i]); };
  //data[n_cols-1] = std::move(safe_log10(data[n_cols-1]));
  if (table_submode < 3) {
    acc.resize(n_cols-1);
    spline.resize(n_cols-1);
    init_numbered_interp(0, &data[0][0], &data[1][0]);
    if (table_submode == 2) { init_numbered_interp(1, &data[0][0], &data[2][0]); };
  } else {
    int n_grids = 1 + (g2>0);
    acc_2d.resize(n_grids);
    spline_2d.resize(n_grids);
    // Make sure that entries are properly sorted
    std::vector<double> x_vec = data[0];
    sort(x_vec.begin(), x_vec.end());
    x_vec.erase(unique(x_vec.begin(), x_vec.end()), x_vec.end());
    int nx = x_vec.size();
    std::vector<double> y_vec = data[1];
    sort(y_vec.begin(), y_vec.end());
    y_vec.erase(unique(y_vec.begin(), y_vec.end()), y_vec.end());
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
    };

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
      for (int i=0; i<n_grids; ++i) {
        gsl_spline2d_set(spline_2d[i], z[i], ind_x, ind_y, data[i+2][j]);
      };
    };
    for (int i=0; i<n_grids; ++i) { gsl_spline2d_init(spline_2d[i], x, y, z[i], nx, ny); };
  };
}

AxionSpectrum::AxionSpectrum(std::string file, double g1, double g2) { init_table_mode(file, g1, g2); };

void AxionSpectrum::init_solar_model_mode(SolarModel* sol, SolarModelMemberFn process2) {
  mode = solar_model; // Do not set default_g1, default_g2 here! This info in the solar model.
  s = sol;
  function2 = process2;
  solar_model_okay = s->is_initialised();
  terminate_with_error_if(not(solar_model_okay), "ERROR! The instance of SolarModel used to initialse an instance of AxionSpectrum has not been initialised properly.");
};

AxionSpectrum::AxionSpectrum(SolarModel* sol, SolarModelMemberFn process2) { init_solar_model_mode(sol, process2); };

void AxionSpectrum::init_analytical_mode(double norm, double g1, double a, double b) {
  mode = analytical; // Do not set default_g1 here! This info in the parameters.
  analytical_parameters = {norm, g1, a, b};
};

AxionSpectrum::AxionSpectrum(double norm, double g1, double a, double b) { init_analytical_mode(norm, g1, a, b); };

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
    };
    if (check) {
      mode = new_mode;
    } else {
      std::cout << "WARNING! Changing mode of AxionSpectrum is not permitted; the class instance has not been initialised for the new mode." << std::endl;
    };
  };
}

AxionSpectrum::~AxionSpectrum() {}; // TODO: Correct destructor

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
        if ((table_submode == 1) && (g2 > 0) && not(g2_warning_issued)) { std::cout << "WARNING! You don't have a second spectrum loaded! The value of g2 > 0 has no effect." << std::endl; g2_warning_issued = true; };
        if ((table_submode < 3) && (*r < 1) && not(r_warning_issued)) { std::cout << "WARNING! You didn't provide data for radii! The value of r < 1 has no effect." << std::endl; r_warning_issued = true; };
        for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
          double ref_g1_flux = 0;
          double ref_g2_flux = 0;
          if (table_submode < 4) {
            ref_g1_flux = gsl_pow_2(g1/default_g1) * gsl_spline_eval(spline[0], *erg, acc[0]);
            if (table_submode > 1) { ref_g2_flux = gsl_pow_2(g2/default_g2) * gsl_spline_eval(spline[1], *erg, acc[1]); };
          } else {
            ref_g1_flux = gsl_pow_2(g1/default_g1) * gsl_spline2d_eval(spline_2d[0], *r, *erg, acc_2d[0].first, acc_2d[0].second);
            if (default_g2 > 0) { ref_g2_flux = gsl_pow_2(g2/default_g2) * gsl_spline2d_eval(spline_2d[1], *r, *erg, acc_2d[1].first, acc_2d[1].second); };
          };
          //g1_flux.push_back(ref_g1_flux);
          //g2_flux.push_back(ref_g2_flux);
          total_flux.push_back(ref_g1_flux+ref_g2_flux);
        };
        break;
      case solar_model :
      {
        static double ref_g1 = s->get_gagg_ref_value_in_inverse_GeV();
        static double ref_g2 = s->get_gaee_ref_value();
        std::vector<double> ref_g1_fluxes = s->calculate_spectral_flux_any(ergs, function1, *r);
        std::vector<double> ref_g2_fluxes = s->calculate_spectral_flux_any(ergs, function2, *r);
        for (int i=0; i<ergs.size(); ++i) {
          double ref_g1_flux = gsl_pow_2(g1/ref_g1) * ref_g1_fluxes[i];
          double ref_g2_flux = gsl_pow_2(g2/ref_g2) * ref_g2_fluxes[i];
          //g1_flux.push_back(ref_g1_flux);
          //g2_flux.push_back(ref_g2_flux);
          total_flux.push_back(ref_g1_flux+ref_g2_flux);
        };
        break;
      }
      case analytical :
        for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
          double ref_g1_flux = analytical_parameters[0] * gsl_pow_2(g1/analytical_parameters[1]) * pow(*erg,analytical_parameters[2]) * exp(-analytical_parameters[3]*(*erg));
          //g1_flux.push_back(ref_g1_flux);
          total_flux.push_back(ref_g1_flux);
        };
        break;
      case undefined :
        terminate_with_error("ERROR! The AxionSpectrum class has not been initialised properly.");
        break;
      default :
        terminate_with_error("FATAL ERROR! The AxionSpectrum mode was somehow set to an unknown value; this is a bug, please report it.");
    };
    result.push_back(total_flux);
  };
  return result;
}

double AxionSpectrum::axion_flux(double erg, double g1) { return axion_flux(erg, g1, 0); }
double AxionSpectrum::axion_flux(double erg, double g1, double g2) { return axion_flux(erg, 1, g1, 0); }
double AxionSpectrum::axion_flux(double erg, double r, double g1, double g2) {
  std::vector<double> ergs = {erg};
  std::vector<double> result = axion_flux(ergs, r, g1, g2);
  return result[0];
}
std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double g1) { return axion_flux(ergs, g1, 0); };
std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double g1, double g2) { return axion_flux(ergs, 1, g1, 0); };
std::vector<double> AxionSpectrum::axion_flux(std::vector<double> ergs, double r, double g1, double g2) {
  std::vector<double> radii = {r};
  std::vector<std::vector<double>> result = axion_flux(ergs, radii, g1, g2);
  return result[0];
}


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
  //double result = 0.0;
  double erg = (p->erg);
  SolarModel* sol = (p->sol);

  // TODO: Check if this is conistent now!
  //if (sol->opcode == OP) {
  //  double element_contrib = 0.0;
  //  // Add opacity terms all non-H or He elements (metals)
  //  for (int k = 2; k < num_op_elements; k++) { element_contrib += sol->Gamma_P_opacity(erg, r, op_element_names[k]); };
  //  result = 0.5*gsl_pow_2(r*erg/pi)*element_contrib;
  //}
  //if ((sol->opcode == OPAS) || (sol->opcode == LEDCOP) || (sol->opcode == ATOMIC)) {
  //    result = 0.5*gsl_pow_2(r*erg/pi) * sol->Gamma_P_opacity(erg, r);
  //};

  //return result;
  return 0.5*gsl_pow_2(r*erg/pi) * sol->Gamma_P_opacity(erg, r);
}

// Includes FF flux and ee contribution as in arxiv[1310.0823].
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

  return 0.5*gsl_pow_2(r*erg/pi) * sol->Gamma_P_all_electron(erg,r);
}

std::vector<double> calculate_spectral_flux(std::vector<double> ergs, Isotope isotope, SolarModel &s, double (*integrand)(double, void*), std::string saveas) {
  // Constant factor for consistent units, i.e. integrated flux will be in units of cm^-2 s^-1 keV^-1.
  const double factor = pow(radius_sol/(1.0e-2*keV2cm),3) / ( pow(1.0e2*distance_sol,2) * (1.0e6*hbar) );
  // = Rsol^3 [in keV^-3] / (2 pi^2 d^2 [in cm^2] * 1 [1 corresponds to s x keV))
  // TODO: Define double norm = f(2.0) and add it to the integration_params with default norm = 1. Integrate function *1/norm and rescale result *norm at the end.
  std::vector<double> results, errors;

  gsl_function f;
  f.function = integrand;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);

  std::ofstream output;
  //if (saveas != "") {
  //  output.open(saveas);
  //  output << "# Spectral flux over full solar volume by " << LIBRARY_NAME << ".\n# Columns: energy values [keV], axion flux [axions/cm^2 s keV], axion flux error estimate [axions/cm^2 s keV]" << std::endl;
  //};

  for (auto erg = ergs.begin(); erg != ergs.end(); erg++) {
    double integral, error;
    integration_params p = {*erg, &s, isotope};
    f.params = &p;
    gsl_integration_qag (&f, s.get_r_lo(), s.get_r_hi(), int_abs_prec, int_rel_prec, int_space_size, int_method_1, w, &integral, &error);
    results.push_back(factor*integral);
    errors.push_back(factor*error);
    //if (saveas != ""){ output << *erg << " " << factor*integral << factor*error << std::endl; };
  };

  //if (saveas!= "") { output.close(); };
  gsl_integration_workspace_free (w);

  std::vector<std::vector<double>> buffer = {ergs, results, errors};
  std::string comment = "Spectral flux over full solar volume by "+LIBRARY_NAME+".\nColumns: energy values [keV], axion flux [axions / cm^2 s keV], axion flux error estimate [axions / cm^2 s keV]";
  if (saveas != ""){ save_to_file(saveas, buffer, comment); };

  return results;
}

std::vector<double> calculate_spectral_flux(std::vector<double> ergs, SolarModel &s, double (*integrand)(double, void*), std::string saveas) { std::string NONE = ""; return calculate_spectral_flux(ergs, NONE, s, integrand, saveas); }

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
  gsl_integration_qag (&f, s->get_r_lo(), s->get_r_hi(), 0.1*int_abs_prec, 0.1*int_rel_prec, int_space_size, int_method_1, w, &result, &error);
  gsl_integration_workspace_free (w);
  return factor*result/normfactor;
}

double calculate_flux(double lowerlimit, double upperlimit, SolarModel &s, Isotope isotope) {
    const double normfactor = 1.0e20;
    double result, error;
    gsl_function f;
    f.function = spectral_flux_integrand;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
    integration_params2 p2 = {&s, &integrand_all_axionelectron, isotope};
    //integration_params2 p2 = {&s, &integrand_Primakoff, isotope};
    f.params = &p2;
    gsl_integration_qag (&f, lowerlimit, upperlimit, int_abs_prec, int_rel_prec, int_space_size, int_method_2, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result*normfactor;
}

double flux_from_file_integrand(double erg, void * params) {
  OneDInterpolator * interp = (OneDInterpolator *)params;
  //std::cout << "DEBUG INFO. flux_from_file_integrand(" << erg << " keV) = " << interp->interpolate(erg) << " ." << std::endl;
  return interp->interpolate(erg);
}

double integrated_flux_from_file(double erg_min, double erg_max, std::string spectral_flux_file, bool includes_electron_interactions) {
  // Peak positions for axion electron interactions
  const std::vector<double> all_peaks = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163, 5.64394,
                                         5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};
  double result, error;

  OneDInterpolator spectral_flux (spectral_flux_file);
  if ( (erg_min < spectral_flux.lower()) || (erg_max > spectral_flux.upper()) ) {
    terminate_with_error("ERROR! The integration boundaries given to 'integrated_flux_from_file' are incompatible with the min/max available energy in the file "+spectral_flux_file+".");
  };

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (int_space_size);
  gsl_function f;
  f.function = &flux_from_file_integrand;
  f.params = &spectral_flux;

  if (includes_electron_interactions) {
    std::vector<double> relevant_peaks;
    relevant_peaks.push_back(erg_min);
    for (auto peak_erg = all_peaks.begin(); peak_erg != all_peaks.end(); peak_erg++) { if ( (erg_min < *peak_erg) && (*peak_erg < erg_max) ) { relevant_peaks.push_back(*peak_erg); }; };
    relevant_peaks.push_back(erg_max);
    gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), int_abs_prec, int_rel_prec, int_space_size, w, &result, &error);
  } else {
    gsl_integration_qag(&f, erg_min, erg_max, abs_prec2, rel_prec2, int_space_size, int_method_1, w, &result, &error);
  };

  gsl_integration_workspace_free (w);

  return result;
}


////////////////////////////////////////////////////////////////////
// Overloaded versions of the functions above for convenient use. //
////////////////////////////////////////////////////////////////////

//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,double r_max, SolarModel &s, double (*integrand)(double, double), std::string saveas) { return calculate_spectral_flux_solar_disc(ergs, r_max, 0, s, integrand, saveas); }
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,Isotope isotope, double r_max, SolarModel &s, double (*integrand)(double, double)) { return calculate_spectral_flux_solar_disc(ergs, r_max, isotope, s, integrand, ""); }
//std::vector<double> calculate_spectral_flux_solar_disc(std::vector<double> ergs,double r_max, SolarModel &s, double (*integrand)(double, double)) { return calculate_spectral_flux_solar_disc(ergs, r_max, 0, s, integrand); }
// -> solar_model.cpp/hpp
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_Primakoff, saveas); }
std::vector<double> calculate_spectral_flux_Primakoff(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas) { double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_Primakoff; return calculate_spectral_flux_solar_disc(ergs, r_max, s, integrand, saveas); }
std::vector<double> calculate_spectral_flux_Compton(std::vector<double> ergs, SolarModel &s,std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_Compton, saveas); }
std::vector<double> calculate_spectral_flux_weightedCompton(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_weightedCompton, saveas); }
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, std::string element, SolarModel &s) { return calculate_spectral_flux(ergs, element, s, &integrand_opacity_element); }
std::vector<double> calculate_spectral_flux_element(std::vector<double> ergs, std::string element, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, element, s, &integrand_opacity_element, saveas); }
std::vector<double> calculate_spectral_flux_all_ff(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_all_ff,saveas); }
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s,std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_all_axionelectron, saveas); }
std::vector<double> calculate_spectral_flux_axionelectron(std::vector<double> ergs, SolarModel &s, double r_max, std::string saveas) { double (SolarModel::*integrand)(double, double) = &SolarModel::Gamma_P_all_electron; return calculate_spectral_flux_solar_disc(ergs, r_max, s, integrand, saveas); }
std::vector<double> calculate_spectral_flux_opacity(std::vector<double> ergs, SolarModel &s, std::string saveas) { return calculate_spectral_flux(ergs, s, &integrand_opacity, saveas); }

////////////////////////////////////
// Monte Carlo-related functions. //
////////////////////////////////////

AxionMCGenerator::AxionMCGenerator() {};

AxionMCGenerator::AxionMCGenerator(std::string spectrum_file) {
  ASCIItableReader spectrum_data (spectrum_file);
  int pts = spectrum_data.getnrow();
  inv_cdf_data_erg = std::vector<double> (pts);
  inv_cdf_data_x = std::vector<double> (pts);

  double norm = 0.0;
  inv_cdf_data_erg[0] = spectrum_data[0][0];
  inv_cdf_data_x[0] = norm;
  for(int i=1; i<pts; ++i) {
    // Trapozoidal rule integration.
    norm += 0.5 * (spectrum_data[0][i] - spectrum_data[0][i-1]) * (spectrum_data[1][i] + spectrum_data[1][i-1]);
    inv_cdf_data_erg[i] = spectrum_data[0][i];
    inv_cdf_data_x[i] = norm;
  };

  integrated_norm = norm;
  for (int i=0; i<pts; ++i) { inv_cdf_data_x[i] = inv_cdf_data_x[i]/integrated_norm; };

  init_inv_cdf_interpolator();
}

AxionMCGenerator::~AxionMCGenerator() {
  gsl_spline_free(inv_cdf);
  gsl_interp_accel_free(inv_cdf_acc);
}

void AxionMCGenerator::init(std::string inv_cdf_file) {
  ASCIItableReader inv_cdf_data (inv_cdf_file);
  inv_cdf_data_x = inv_cdf_data[0];
  inv_cdf_data_erg = inv_cdf_data[1];

  terminate_with_error_if(inv_cdf_data_x.end()[-2] > 1.0, "ERROR! Sanity check for MC generator failed! The second to last entry of your inverse CDF is greater than 1.");
  integrated_norm = 1.0;

  init_inv_cdf_interpolator();
}

void AxionMCGenerator::init_inv_cdf_interpolator() {
  int pts = inv_cdf_data_x.size();
  terminate_with_error_if(not(pts >= 2), "ERROR! You tried to initialise the MC generator with less than two data points.");

  const double* prob = &inv_cdf_data_x[0];
  const double* erg = &inv_cdf_data_erg[0];

  inv_cdf_acc = gsl_interp_accel_alloc();
  inv_cdf = gsl_spline_alloc(gsl_interp_linear, pts);
  gsl_spline_init(inv_cdf, prob, erg, pts);

  mc_generator_ready = true;
}

void AxionMCGenerator::save_inv_cdf_to_file(std::string inv_cdf_file) {
  const std::string comment = "Inverse CDF obtained by "+LIBRARY_NAME;
  std::vector<std::vector<double>> buffer;
  buffer.push_back(inv_cdf_data_x);
  buffer.push_back(inv_cdf_data_erg);
  save_to_file(inv_cdf_file, buffer, comment);
}

double AxionMCGenerator::evaluate_inv_cdf(double x) {
  terminate_with_error_if(not(mc_generator_ready), "ERROR! Your MC generator has not been initialised properly!");
  return gsl_spline_eval(inv_cdf, x, inv_cdf_acc);
}

std::vector<double> AxionMCGenerator::draw_axion_energies(int n) {
  std::vector<double> result;

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 rng (seed);
  std::uniform_real_distribution<double> unif(0, 1);

  //std::cout << "Test U[0,1]: " << unif(rng) << std::endl;
  for (int i=0; i<n; ++i) { result.push_back(evaluate_inv_cdf( unif(rng) )); };

  return result;
}

double AxionMCGenerator::get_norm() {
  terminate_with_error_if(not(mc_generator_ready), "ERROR! Your MC generator has not been initialised properly!");
  return integrated_norm;
}

AxionMCGenerator::AxionMCGenerator(SolarModel s, SolarModelMemberFn process, double omega_min, double omega_max, double omega_delta, double r_max) {
  int n_omega_vals = int((omega_max-omega_min)/omega_delta);
  inv_cdf_data_erg = std::vector<double> (n_omega_vals);
  for (int i=0; i<n_omega_vals; ++i) { inv_cdf_data_erg[i] = omega_min + i*omega_delta; };

  if (r_max < 1.0) {
    inv_cdf_data_x = calculate_spectral_flux_solar_disc(inv_cdf_data_erg, r_max, s, process, "");
  } else {
    // TODO: Use non-disc integration here AND/OR(!) check in calculate_spectral_flux_solar_disc if r_max > s.get_r_hi() and switch there!
    inv_cdf_data_x = calculate_spectral_flux_solar_disc(inv_cdf_data_erg, 1.0, s, process, "");
  };

  double norm = 0.0;
  inv_cdf_data_x[0] = norm;
  for(int i=1; i<n_omega_vals; ++i) {
    // Trapozoidal rule integration.
    norm += 0.5 * omega_delta * (inv_cdf_data_x[i] + inv_cdf_data_x[i-1]);
    inv_cdf_data_x[i] = norm;
  };

  integrated_norm = norm;
  for (int i=0; i<n_omega_vals; ++i) { inv_cdf_data_x[i] = inv_cdf_data_x[i]/integrated_norm; };

  init_inv_cdf_interpolator();
}

/* REQUIRES e.g. the BOOST library
// TODO: Replace this with one's own routines?
// Inverse CDF for the Primakoff spectrum approximation
// d Phi_a (omega) / d omega ~ omega^a * exp(-b * omega/keV)
double inv_cdf_Primakoff_analytic_approx(double x, double a, double b) {
  double ap1 = 1.0 + a;
  return boost::gamma_q_inv(ap1, (1.0 - x)/pow(b,ap1) ) / b;
}
*/
