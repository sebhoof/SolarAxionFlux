#include "utils.hpp"

void terminate_with_error(std::string err_string) {
  std::cerr << err_string << std::endl;
  exit(EXIT_FAILURE);
};

// ASCII Table Reader by Christoph Weniger
int ASCIItableReader::read(std::string filename) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (in.fail()) {
    std::cout << "ERROR. Failed loading: " << filename << std::endl;
    exit(-1);
  }
  std::string line;
  while(std::getline(in, line)) {
    if (line[0] == '#') continue;
    std::stringstream ss(line);

    size_t i = 0;
    double tmp;
    while(ss >> tmp) {
      if ( i+1 > data.size() ) data.resize(i+1);
      data[i].push_back(tmp);
      i++;
    }
  }
  in.close();
  return 0;
}

void ASCIItableReader::setcolnames(std::vector<std::string> names) {
  if ( (int) names.size() == data.size() ) {
    size_t i = 0;
    for (auto it = names.begin(); it != names.end(); it++) {
      colnames[*it] = i;
      i++;
    }
  } else {
    std::cout << "Warning in ASCIItableReader: Column number incompatible." << std::endl;
  }
}

// Constructors
SolarModel::SolarModel() {};
SolarModel::SolarModel(std::string file)
{
  data = ASCIItableReader(file);
  int pts = data.getnrow();
  // Terminate if number of columns is wrong; i.e. the wrong solar model file format.
  if (data.getncol() != 35) { terminate_with_error("ERROR! Solar model file '"+file+"' not compatible with this code!"); };
  data.setcolnames("mass", "radius", "temperature", "rho", "Pressure", "Luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_C13", "X_N14", "X_N15", "X_O16", "X_O17", "X_O18", "X_Ne", "X_Na", "X_Mg", "X_Al", "X_Si", "X_P", "X_S", "X_Cl", "X_Ar",
                   "X_K", "X_Ca", "X_Sc", "X_Ti", "X_V", "X_Cr", "X_Mn", "X_Fe", "X_Co", "X_Ni");

  // Extract the radius from the files (in units of the solar radius).
  r_lo = data["radius"][0];
  r_hi = data["radius"][pts-1];

  // Extract the temperature from the files (has to be converted into keV) & calculate the screening scale kappa_s_squared.
  // Initialise necessary variables for the screening scale calculation.
  std::vector<double> temperature;
  //std::vector<double> density;
  std::vector<double> kappa_s_sq;
  std::vector<double> w_pl_sq;
  std::vector<double> n_e;
  std::vector<std::vector<double>> rho_iz (15);
  std::vector<std::vector<double>> z2_n_iz (15);
  // Multiplicative factor: (4 pi alpha_EM / atomic_mass_unit) x (1 g/cm^3) in units of keV^3
  const double factor = 4.0*pi*alpha_EM*gsl_pow_3(1.0E+6*gev2cm)/((1.0E+9*eV2g)*atomic_mass_unit);
  // Atomic weight of species i (exact weight if isotope is known OR estimate from average solar abundance from data if available OR estimate from natural terrestrial abundance).
  const double A_vals [29] = {1.007825, 4.002603, 3.016029, 12.000000, 13.003355, 14.003074, 15.000109, 15.994915, 16.999132, 17.999160,
                              20.1312812, 22.989769, 24.3055, 26.9815385, 28.085, 30.973762, 32.0675, 35.4515, 36.275403, 39.0983, 40.078, 44.955908, 47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934};
  // Ionisation of species i assuming full ionisation.
  const double Z_vals [29] = {1.0,      2.0,      2.0,       6.0,       6.0,       7.0,       7.0,       8.0,       8.0,       8.0,
                              10.0,       11.0,      12.0,    13.0,       14.0,   15.0,      16.0,    17.0,    18.0,      19.0,    20.0,   21.0,      22.0,   23.0,    24.0,    25.0,      26.0,   27.0,      28.0};

  #ifdef DEBUG_MODE
    std::cout << "DEBUGGING INFO for solar models:\nradius/Rsol T/K kappa_s^2/keV^2 omega_pl^2/keV^2" << std::endl;
  #endif

  // Linearly extrapolate the data in the solar model file to r = 0 if necessary.
  if (r_lo > 0)
  {
    double r0 = data["radius"][0], r1 = data["radius"][1];

    for (int i = 0; i < 29; i++)
    {
      double intercept = (r0*data[i][1]-r1*data[i][0])/(r0-r1);
      //data[i].insert(data[i].begin(), intercept);
      data.prepend_data(intercept, i);
    };
    r_lo = 0.0;
    pts += 1;
  };
  const double* radius = &data["radius"][0];
  // Calculate the necessary quantities -- T(r), kappa_s^2(r) and omega_pl^2(r) -- and store them internally.
  for (int i = 0; i < pts; i++)
  {
    double sum = 0.0;
    double ne = 0.0;
    temperature.push_back((1.0E-3*K2eV)*data["temperature"][i]);
    for (int j = 0; j < 29; j++)
    {
      double temp = data[j+6][i]*Z_vals[j]/A_vals[j];
      ne += temp;
      sum += temp*(1.0 + Z_vals[j]);
    };
    kappa_s_sq.push_back(factor*sum*data["rho"][i]/temperature[i]);
    w_pl_sq.push_back(factor*ne*data["rho"][i]/(1.0E+6*m_electron));
    double rhorel = data["rho"][i]/((1.0E+9*eV2g)*atomic_mass_unit);
    n_e.push_back(ne*rhorel);
    for (int j = 0; j < 15; j++)
    {
      double z2_n = 0.0;
      double dens = 0.0;
      for (auto it = op_elements[j].begin(); it != op_elements[j].end(); ++it)
      {
        z2_n += Z_vals[*it-6]*Z_vals[*it-6]*data[*it][i]/A_vals[*it-6];
        dens += data[*it][i];
      };
      z2_n_iz[j].push_back(z2_n*rhorel);
      rho_iz[j].push_back(atomic_mass_unit*dens*rhorel);
    };
    #ifdef DEBUG_MODE
      printf("%5.4f %1.6e %1.6e %1.6e\n", data["radius"][i], temperature[i], kss, wpls);
    #endif
  };
  // Set up the interpolating functions for temperature and screening scale.
  accel[0] = gsl_interp_accel_alloc ();
  linear_interp[0] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* temp_vals = &temperature[0];
  gsl_spline_init (linear_interp[0], radius, temp_vals, pts);
  accel[1] = gsl_interp_accel_alloc ();
  linear_interp[1] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* kappa_squared_vals = &kappa_s_sq[0];
  gsl_spline_init (linear_interp[1], radius, kappa_squared_vals, pts);
  accel[2] = gsl_interp_accel_alloc ();
  linear_interp[2] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* omega_pl_squared_vals = &w_pl_sq[0];
  gsl_spline_init (linear_interp[2], radius, omega_pl_squared_vals, pts);
  accel[3] = gsl_interp_accel_alloc ();
  linear_interp[3] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* n_e_vals = &n_e[0];
  gsl_spline_init (linear_interp[3], radius, n_e_vals, pts);
  //accel[4] = gsl_interp_accel_alloc ();
  //linear_interp[4] = gsl_spline_alloc (gsl_interp_linear, pts);
  //const double* density_vals = &density[0];
  //gsl_spline_init (linear_interp[4], radius, density_vals, pts);

  //logger() << LogTags::info << "Initialisation of solar model from file '"+file+"' complete!" << std::endl;
  //logger() << LogTags::debug << "Entries in model file: " << pts << " for solar radius in [" << data["radius"][0] << ", " << data["radius"][pts-1] << "]." << EOM;

  // Do opacities
  for (int iz = 0; iz < 15; iz++) {
    z2_n_iz_acc.push_back( gsl_interp_accel_alloc() );
    z2_n_iz_lin_interp.push_back( gsl_spline_alloc(gsl_interp_linear, pts) );
    const double* z2_n_iz_vals = &z2_n_iz[iz][0];
    gsl_spline_init (z2_n_iz_lin_interp[iz], radius, z2_n_iz_vals, pts);
    rho_iz_acc.push_back( gsl_interp_accel_alloc() );
    rho_iz_lin_interp.push_back( gsl_spline_alloc(gsl_interp_linear, pts) );
    const double* rho_iz_vals = &rho_iz[iz][0];
    gsl_spline_init (rho_iz_lin_interp[iz], radius, rho_iz_vals, pts);

    std::map<std::pair<int,int>, gsl_interp_accel*> temp1;
    std::map<std::pair<int,int>, gsl_spline*> temp2;

    for (int j = 0; j < 197; j++){
      std::string op_filename = "data/opacity_tables/opacity_table_"+std::to_string(iz+1)+"_"+std::to_string(op_grid[j][0])+"_"+std::to_string(op_grid[j][1])+".dat";
      ASCIItableReader op_data = ASCIItableReader(op_filename);

      // Determine the number of interpolated mass values.
      int op_pts = op_data[0].size();
      auto pr = std::make_pair(op_grid[j][0], op_grid[j][1]);
      //std::cout << pr.first << " " << pr.second << std::endl;
      temp1[pr] = gsl_interp_accel_alloc();
      temp2[pr] = gsl_spline_alloc (gsl_interp_linear, op_pts);
      const double* omega = &op_data[0][0];
      const double* opacity = &op_data[1][0];
      gsl_spline_init (temp2[pr], omega, opacity, op_pts);
    };
    opacity_acc.push_back(temp1);
    opacity_lin_interp.push_back(temp2);
  };
}

// Move assignment operator
SolarModel& SolarModel::operator=(SolarModel &&model)
{
  if (this != &model)
  {
    std::swap(data,model.data);
    std::swap(accel,model.accel);
    std::swap(linear_interp, model.linear_interp);
  }
  return *this;
}

// Class destructor
SolarModel::~SolarModel()
{
  for (auto interp : linear_interp)
    gsl_spline_free (interp);
  for (auto acc : accel)
    gsl_interp_accel_free (acc);
}

// Routine to return the temperature (in keV) of the zone around the distance r from the centre of the Sun.
double SolarModel::temperature_in_keV(double r) { return gsl_spline_eval(linear_interp[0], r, accel[0]); }

// Routine to return the screening paramter kappa^2 in units of keV^2 (kappa^-1 = Debye-Hueckel radius).
double SolarModel::kappa_squared(double r)
{
  // Interpolated value, directly from the Solar model.
  return gsl_spline_eval(linear_interp[1], r, accel[1]);
}

// Routine to return the number density times Z^2 of ion iz in the zone around the distance r from the centre of the Sun.
double SolarModel::z2_n_iz(double r, int iz) { return gsl_spline_eval(z2_n_iz_lin_interp[iz], r, z2_n_iz_acc[iz]); }

// Routine to return the mass density of ion iz in the zone around the distance r from the centre of the Sun.
double SolarModel::rho_iz(double r, int iz) { return gsl_spline_eval(rho_iz_lin_interp[iz], r, rho_iz_acc[iz]); }

// Routine to return the electron density in the zone around the distance r from the centre of the Sun.
double SolarModel::n_e(double r) { return gsl_spline_eval(linear_interp[3], r, accel[3]); }

// Routine to return the plasma freqeuency squared (in keV^2) of the zone around the distance r from the centre of the Sun.
double SolarModel::omega_pl_squared(double r) { return gsl_spline_eval(linear_interp[2], r, accel[2]); }

/* //Old integrand function
double integrand1(double t, void * params) {
  // Retrieve parameters and other integration variables.
  double y = *(double *) params;

  return gsl_pow_3(t)/gsl_pow_2(t*t + y*y);
};


double integrand2(double x, void * params) {
  // Retrieve parameters and other integration variables.
  struct integrand_params * p = (struct integrand_params *)params;
  double u = (p->u);
  double y = (p->y);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1E6);
  double result, error;

  gsl_function f;
  f.function = &integrand1;
  f.params = &y;

  //gsl_set_error_handler_off();
  gsl_integration_qag (&f, sqrt(x*x + u) - x, sqrt(x*x + u) - x, 0.1*abs_prec, 0.1*rel_prec, 1e6, method, w, &result, &error);
  //printf ("GSL status: %s\n", gsl_strerror (status));
  //gsl_integration_qags(&F, rad, rmax, 1e-1*abs_prec, 1e-1*rel_prec, 1E6, w, &result, &error);
  gsl_integration_workspace_free (w);

  return x * exp(-x*x) * result;
}

double aux_function(double u, double y) {
  integrand_params p = {u, y};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1E6);
  double result, error;

  gsl_function f;
  f.function = &integrand2;
  f.params = &p;

  //gsl_set_error_handler_off();
  //gsl_integration_qag (&f, sqrt(x*x + u) - x, sqrt(x*x + u) - x, 0.1*abs_prec, 0.1*rel_prec, 1e6, method, w, &result, &error);
  gsl_integration_qagiu (&f, 0, abs_prec, rel_prec, 1e6, w, &result, &error);
  //printf ("GSL status: %s\n", gsl_strerror (status));
  //gsl_integration_qags(&F, rad, rmax, 1e-1*abs_prec, 1e-1*rel_prec, 1E6, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result;
}
*/

double integrand(double x, void * params) {
  // Retrieve parameters and other integration variables.
  struct integrand_params * p = (struct integrand_params *)params;
  double u = (p->u);
  double y = (p->y);

  double y2 = y*y;
  double a2 = sqrt(x*x + u) + x;
  a2 = a2*a2;
  double b2 = sqrt(x*x + u) - x;
  b2 = b2*b2;
  double analytical_integral = 0.5 * ( 1.0/(1.0 + y2/b2) - 1.0/(1.0 + y2/a2) + log((a2 + y2) / (b2 + y2)) );

  return x * exp(-x*x) * analytical_integral;
}


double aux_function(double u, double y) {
  integrand_params p = {u, y};

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1E6);
  double result, error;

  gsl_function f;
  f.function = &integrand;
  f.params = &p;

  //gsl_set_error_handler_off();
  //gsl_integration_qag (&f, sqrt(x*x + u) - x, sqrt(x*x + u) - x, 0.1*abs_prec, 0.1*rel_prec, 1e6, method, w, &result, &error);
  gsl_integration_qagiu (&f, 0, abs_prec, rel_prec, 1e6, w, &result, &error);
  //printf ("GSL status: %s\n", gsl_strerror (status));
  //gsl_integration_qags(&F, rad, rmax, 1e-1*abs_prec, 1e-1*rel_prec, 1E6, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result;
}


// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ff_full(double omega, double r, int iz) {
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * pow(alpha_EM*1.0e-13,2) * pow(1.0e6*gev2cm,6);
  double u = omega/temperature_in_keV(r);
  double y_red = sqrt(kappa_squared(r)/(2.0*1.0e6*m_electron*temperature_in_keV(r)));
  return prefactor1 * n_e(r)*z2_n_iz(r,iz)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(1.0e6*m_electron,3.5));
}

// Calculate the e-e bremsstrahlung contribution; from Eq. (2.18) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ee(double omega, double r) {
  // N.B. "y" and "prefactor2" are different from the "y_red" and "prefactor1" above.
  const double prefactor2 = (4.0*sqrt(pi)/3.0) * pow(alpha_EM*1.0e-13,2) * pow(1.0e6*gev2cm,6);
  double u = omega/temperature_in_keV(r);
  double y = sqrt(kappa_squared(r)/(1.0e6*m_electron*temperature_in_keV(r)));
  return prefactor2 * n_e(r)*n_e(r)*exp(-u)*aux_function(u,y) / (omega*sqrt(temperature_in_keV(r))*pow(1.0e6*m_electron,3.5));
}

// Calculate the Compton contribution; from Eq. (2.19) of [arXiv:1310.0823]
double SolarModel::Gamma_P_Compton (double omega, double r) {
  const double prefactor3 = (alpha_EM/3.0) * pow(1.0e-13/(1.0e6*m_electron),2) * pow(1.0e6*gev2cm,3);
  double u = omega/temperature_in_keV(r);
  double v = omega/(1.0e6*m_electron);
  return prefactor3 * v*v*n_e(r)/gsl_expm1(u);
}

// Read off interpolated elements
double SolarModel::op_grid_interp_erg (double u, int ite, int jne, int iz) {
  auto key = std::make_pair(ite,jne);
  return gsl_spline_eval(opacity_lin_interp[iz].at(key), log(u), opacity_acc[iz].at(key));
};
  //double result = interp1d(np.log(s[:,0]),s[:,1],bounds_error=False,fill_value=0,kind='linear')(np.log(u))
//  re
//}

double SolarModel::opacity_table_interpolator (double omega, double r, int iz) {
  // Need tempeature in Kelvin
  double temperature = temperature_in_keV(r)/(1.0e-3*K2eV);
  double ne = n_e(r);
  double ite = 40.0*log10(temperature);
  double jne = 4.0*log10(ne);
  int ite2 = int(ceil(20.0*log10(temperature))*2);
  int ite1 = ite2 - 2;
  // Need omega in Kelvin
  double u1 = omega/(1.0e-3*K2eV*pow(10,double(ite1)/40.0));
  double u2 = omega/(1.0e-3*K2eV*pow(10,double(ite2)/40.0));
  int jne2 = int(ceil(log10(ne)*2)*2);
  int jne1 = jne2 - 2;
  double t1 = (ite-double(ite1))/2.0;
  double t2 = (jne-double(jne1))/2.0;
  std::cout << "Test: " << u1 << " " << u2 << " | " << ite1 << " " << ite2 << " | " << jne1 << " " << jne2 << std::endl;
  double result = (1.0-t1)*(1.0-t2)*op_grid_interp_erg(u1,ite1,jne1,iz)
                + (1.0-t1)*t2*op_grid_interp_erg(u1,ite1,jne2,iz)
                + t1*(1.0-t2)*op_grid_interp_erg(u2,ite2,jne1,iz)
                + t1*t2*op_grid_interp_erg(u2,ite2,jne2,iz);
  return result;
};

double SolarModel::opacity (double omega, double r, int iz) {
  const double prefactor4 = a_Bohr*a_Bohr*(1.0e6*gev2cm);
  double u = omega/temperature_in_keV(r);
  return rho_iz(r, iz)*opacity_table_interpolator(omega, r, iz)*gsl_expm1(-u);
};

double SolarModel::Gamma_P_element (double omega, double r, int iz) {
  const double prefactor5 = 0.5*1.0e-13*1.0e-13*4.0*pi/alpha_EM;
  double u = omega/temperature_in_keV(r);
  double result = 0.0;
  if (u < 17.5) {
    double v = omega/(1.0e6*m_electron);
    result = prefactor5*v*v*opacity(omega,r,iz)/gsl_expm1(u);
  };

  return result;
}

double SolarModel::Gamma_P_Primakoff (double erg, double r) {
  // N.B. gagg = 10^-16 keV^-1 = 10^-19 eV^-1
  const double prefactor6 = 1.0e-16*1.0e-16/(32.0*pi*alpha_EM);

  // Get kappa_s^2, omega_plasma^2 and the temperature.
  double ks_sq = kappa_squared(r);
  double w_pl_sq = omega_pl_squared(r);
  double T_in_keV = temperature_in_keV(r);

  // Calculate the flux.
  double x = 4.0*(erg*erg)/ks_sq;
  double y = w_pl_sq/(erg*erg);
  double energy_factor = sqrt(1.0 - y)/gsl_expm1(erg/T_in_keV);
  double rate = (ks_sq*T_in_keV)*((1.0 + 1.0/x)*gsl_log1p(x) - 1.0);

  return  prefactor6*energy_factor*rate;
};
