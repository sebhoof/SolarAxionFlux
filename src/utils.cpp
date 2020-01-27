#include "utils.hpp"

void terminate_with_error(std::string err_string) {
  std::cerr << err_string << std::endl;
  exit(EXIT_FAILURE);
};
void my_handler (const char * reason, const char * file, int line, int gsl_errno) {
    if (gsl_errno == GSL_EDOM) { }
    else {
        std::cout << reason << " in " << file << " line: " << line << std::endl;
        std::cout << "error number: " << gsl_errno << std::endl;
        abort();
    }
}

/// Check if a file exists
bool file_exists(const std::string& filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

// Initialiser for the OneDInterpolator class.
void OneDInterpolator::init(std::string file, std::string type)
{
  // Check if file exists.
  if (not(file_exists(file)))
  {
    //DarkBit_error().raise(LOCAL_INFO, "ERROR! File '"+file+"' not found!");
    terminate_with_error("ERROR! File '"+file+"' not found!");
  } else {
    //logger() << LogTags::debug << "Reading data from file '"+file+"' and interpolating it with '"+type+"' method." << EOM;
  };
  // Read numerical values from data file.
  ASCIItableReader tab (file);
  tab.setcolnames("x", "y");
  // Initialise gsl interpolation routine.
  int pts = tab["x"].size();
  const double* x = &tab["x"][0];
  const double* y = &tab["y"][0];
  acc = gsl_interp_accel_alloc ();
  if (type == "cspline")
  {
    spline = gsl_spline_alloc (gsl_interp_cspline, pts);
  }
  else if (type == "linear")
  {
    spline = gsl_spline_alloc (gsl_interp_linear, pts);
  }
  else
  {
    //DarkBit_error().raise(LOCAL_INFO, "ERROR! Interpolation type '"+type+"' not known to class OneDInterpolator.\n       Available types: 'linear' and 'cspline'.");
    terminate_with_error("ERROR! Interpolation type '"+type+"' not known to class OneDInterpolator.\n       Available types: 'linear' and 'cspline'.");
  };
  gsl_spline_init (spline, x, y, pts);
  // Get first and last value of the "x" component.
  lo = tab["x"].front();
  up = tab["x"].back();
};

// Overloaded class creators for the OneDInterpolator class using the init function above.
OneDInterpolator::OneDInterpolator(std::string file, std::string type) { init(file, type); };
OneDInterpolator::OneDInterpolator(std::string file) { init(file, "linear"); };
OneDInterpolator::OneDInterpolator() {};

// Move assignment operator
OneDInterpolator& OneDInterpolator::operator=(OneDInterpolator&& interp)
{
  if(this != &interp)
  {
    std::swap(acc,interp.acc);
    std::swap(spline,interp.spline);
    std::swap(lo,interp.lo);
    std::swap(up,interp.up);
  }
  return *this;
}

// Destructor
OneDInterpolator::~OneDInterpolator()
{
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}

// Routine to access interpolated values.
double OneDInterpolator::interpolate(double x) { return gsl_spline_eval(spline, x, acc); };

// Routines to return upper and lower boundaries of interpolating function
double OneDInterpolator::lower() { return lo; };
double OneDInterpolator::upper() { return up; };

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
SolarModel::SolarModel() : opcode(OP) {};    // default constructor (not functional)
SolarModel::SolarModel(std::string file,opacitycode set_opcode) : SolarModel(file,set_opcode,true){};  //if raff approx is not set, then true
SolarModel::SolarModel(std::string file, opacitycode set_opcode, bool set_raffelt_approx) : opcode(set_opcode) {
  if ((set_opcode != OP) && (file != "data/SolarModel_AGSS09.dat")){
      std::cout << "Warning: The chosen opacity code is only compatible with the solar model AGSS09." << std::endl;
      std::cout << "         Results will be inconsistent." << std::endl;
  }
  raffelt_approx = set_raffelt_approx;
  data = ASCIItableReader(file);
  int pts = data.getnrow();
  // Terminate if number of columns is wrong; i.e. the wrong solar model file format.
  int n_cols = data.getncol();
  if (n_cols == 35)
  {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_C13", "X_N14", "X_N15", "X_O16", "X_O17", "X_O18", "X_Ne", "X_Na", "X_Mg", "X_Al", "X_Si", "X_P", "X_S", "X_Cl", "X_Ar",
                     "X_K", "X_Ca", "X_Sc", "X_Ti", "X_V", "X_Cr", "X_Mn", "X_Fe", "X_Co", "X_Ni");
  }
  else if (n_cols == 12)
  {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_N14", "X_O16");
  }
  else
  {
    terminate_with_error("ERROR! Solar model file '"+file+"' not compatible with this code!");
  };

  // Extract the radius from the files (in units of the solar radius).
  r_lo = data["radius"][0];
  r_hi = data["radius"][pts-1];

  // Extract the temperature (has to be converted into keV), density, electron density, and ion density * charge^2  from the files.
  // Initialise necessary variables for the screening scale calculation.
  std::vector<double> temperature;
  std::vector<double> n_e;
  std::vector<double> n_e_Raff;
  std::vector<double> z2_n_Raff;
  std::vector<double> density;
  std::vector<std::vector<double>> n_iz (n_op_elements);
  std::vector<std::vector<double>> z2_n_iz (n_op_elements);
  // Multiplicative factor: (4 pi alpha_EM / atomic_mass_unit) x (1 g/cm^3) in units of keV^3
  const double factor = 4.0*pi*alpha_EM*gsl_pow_3(keV2cm)/((1.0E+9*eV2g)*atomic_mass_unit);
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

    for (int i = 0; i < n_cols; i++)
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
//    temperature
    temperature.push_back((1.0E-3*K2eV)*data["temperature"][i]);     //temperature
//    density
    density.push_back(data["rho"][i]);

    
//    electron density
    double rhorel = data["rho"][i]/((1.0E+9*eV2g)*atomic_mass_unit);
//    electron density from Raffelt
    double ne_Raff = 0.5 * (1.0 + data["X_H1"][i]) * data["rho"][i] /(atomic_mass_unit*eV2g*1.0E+9);
    n_e_Raff.push_back(ne_Raff);
//    electron density form pressure (currently not used)
    double radiation_pressure = 4.0/3.0*5.678e-15*pow(data["temperature"][i],4.0);                    //radiation pressure
    double ion_number_dens = 0.0;
    for (int l = 0; l<29;l++) {ion_number_dens += data[6+l][i]*rhorel/A_vals[l];}
    double ne_press = (data["Pressure"][i]-radiation_pressure)/data["temperature"][i]/1.381e-16-ion_number_dens;
//    electron density from summing over all elements (full ionisation)
    double ne = 0.0;
    for (int j = 0; j < 29; j++) { ne += data[j+6][i]*Z_vals[j]/A_vals[j];}
    n_e.push_back(ne*rhorel);
//    ion density weighted by charge^2
//    ion density weighted by charge^2 from Raffelt
    z2_n_Raff.push_back(data["rho"][i] /(atomic_mass_unit*eV2g*1.0E+9));
//    ion density weighted by charge^2 from summing over all elements (full ionisation)
    for (int k = 0; k < n_op_elements; k++)
    {
      double z2_n = 0.0;
      double n = 0.0;
      for (auto it = op_elements[k].begin(); it != op_elements[k].end(); ++it)
      {
        int j = *it;
        n += data[j+6][i]/A_vals[j];
        z2_n += Z_vals[j]*Z_vals[j]*data[j+6][i]/A_vals[j];
      };
      z2_n_iz[k].push_back(z2_n*rhorel);
      n_iz[k].push_back(n*rhorel);
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
  const double* n_e_vals = &n_e[0];
  gsl_spline_init (linear_interp[1], radius, n_e_vals, pts);
  accel[2] = gsl_interp_accel_alloc ();
  linear_interp[2] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* n_e_vals_Raff = &n_e_Raff[0];
  gsl_spline_init (linear_interp[2], radius, n_e_vals_Raff, pts);
  accel[3] = gsl_interp_accel_alloc ();
  linear_interp[3] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* z2_n_vals_Raff = &z2_n_Raff[0];
  gsl_spline_init (linear_interp[3], radius, z2_n_vals_Raff, pts);
  accel[4] = gsl_interp_accel_alloc ();
  linear_interp[4] = gsl_spline_alloc (gsl_interp_linear, pts);
  const double* density_vals = &density[0];
  gsl_spline_init (linear_interp[4], radius, density_vals, pts);

  //Quantities depnding on specfific element
  for (int iz = 0; iz < n_op_elements; iz++) {
    z2_n_iz_acc.push_back( gsl_interp_accel_alloc() );
    z2_n_iz_lin_interp.push_back( gsl_spline_alloc(gsl_interp_linear, pts) );
    const double* z2_n_iz_vals = &z2_n_iz[iz][0];
    gsl_spline_init (z2_n_iz_lin_interp[iz], radius, z2_n_iz_vals, pts);
    n_iz_acc.push_back( gsl_interp_accel_alloc() );
    n_iz_lin_interp.push_back( gsl_spline_alloc(gsl_interp_linear, pts) );
    const double* n_iz_vals = &n_iz[iz][0];
    gsl_spline_init (n_iz_lin_interp[iz], radius, n_iz_vals, pts);
    //OP opacities
    if (opcode == OP) {
        std::map<std::pair<int,int>, gsl_interp_accel*> temp1;
        std::map<std::pair<int,int>, gsl_spline*> temp2;
        for (int j = 0; j < op_grid_size; j++){
          std::string op_filename = "data/opacity_tables/OP/opacity_table_"+std::to_string(iz+1)+"_"+std::to_string(op_grid[j][0])+"_"+std::to_string(op_grid[j][1])+".dat";
          ASCIItableReader op_data = ASCIItableReader(op_filename);

          // Determine the number of interpolated mass values.
          int op_pts = op_data[0].size();
          auto pr = std::make_pair(op_grid[j][0], op_grid[j][1]);
          temp1[pr] = gsl_interp_accel_alloc();
          temp2[pr] = gsl_spline_alloc (gsl_interp_linear, op_pts);
          const double* omega = &op_data[0][0];
          const double* s = &op_data[1][0];
          //std::vector<double> op_data_exp;
          //for (int c=0; c < op_pts; c++){op_data_exp.push_back(exp(op_data[0][c]));}
          //const double* omega_exp = &op_data_exp[0];
          gsl_spline_init (temp2[pr], omega, s, op_pts);
        };
        opacity_acc_op.push_back(temp1);
        opacity_lin_interp_op.push_back(temp2);
    }
  };
  //LEDCOP & ATOMIC (both TOPS) opacitites
  std::string name;
  if (opcode == LEDCOP){
      tops_grid = ledcop_grid;
      tops_temperatures = ledcop_temperatures;
      tops_densities = ledcop_densities;
      name = "LEDCOP";
  }
  if (opcode == ATOMIC){
        tops_grid = atomic_grid;
        tops_temperatures = atomic_temperatures;
        tops_densities = atomic_densities;
        name = "ATOMIC";
  }
  if ((opcode == LEDCOP) || (opcode == ATOMIC)) {
      for (int j = 0; j < tops_grid.size(); j++){
          std::stringstream Tstream;
          std::stringstream rhostream;
          Tstream << std::fixed << std::setprecision(3) << tops_grid[j][0];
          rhostream << std::fixed << std::setprecision(3) << tops_grid[j][1];
          std::string tops_filename = "data/opacity_tables/"+name+"/T"+Tstream.str()+"Rho"+rhostream.str()+".dat";
          ASCIItableReader tops_data = ASCIItableReader(tops_filename);
          // Determine the number of interpolated energy values.
          int tops_pts = tops_data[0].size();
          auto pr = std::make_pair(tops_grid[j][0], tops_grid[j][1]);
          opacity_acc_tops[pr] = gsl_interp_accel_alloc();
          opacity_lin_interp_tops[pr] = gsl_spline_alloc (gsl_interp_linear, tops_pts);
          const double* omega = &tops_data[0][0];
          const double* s = &tops_data[1][0];
          gsl_spline_init (opacity_lin_interp_tops[pr], omega, s, tops_pts);
        };
    }
    if (opcode == OPAS) {
        for (int j = 0 ; j < opas_radii.size(); j++) {
            std::stringstream Rstream;
            Rstream << std::fixed << std::setprecision(2) << opas_radii[j];
            std::string opas_filename = "data/opacity_tables/OPAS/R" +Rstream.str() +".dat";
            ASCIItableReader opas_data = ASCIItableReader(opas_filename);
            // Determine the number of interpolated energy values.
            int opas_pts = opas_data[0].size();
            double rad = opas_radii[j];
            opacity_acc_opas[rad] = gsl_interp_accel_alloc();
            opacity_lin_interp_opas[rad] = gsl_spline_alloc (gsl_interp_linear, opas_pts);
            const double* omega = &opas_data[0][0];
            const double* s = &opas_data[1][0];
            gsl_spline_init (opacity_lin_interp_opas[rad], omega, s, opas_pts);
        }
    }
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
  for (auto interp : linear_interp) { gsl_spline_free (interp); };
  for (auto interp : n_iz_lin_interp) { gsl_spline_free (interp); };
  for (auto interp : z2_n_iz_lin_interp) { gsl_spline_free (interp); };
  for (auto acc : accel) { gsl_interp_accel_free (acc); };
  for (auto acc : n_iz_acc) { gsl_interp_accel_free (acc); };
  for (auto acc : z2_n_iz_acc) { gsl_interp_accel_free (acc); };
}

// Routine to return the temperature (in keV) of the zone around the distance r from the centre of the Sun.
double SolarModel::temperature_in_keV(double r) { return gsl_spline_eval(linear_interp[0], r, accel[0]); }
double SolarModel::density(double r) { return gsl_spline_eval(linear_interp[4], r, accel[4]); }

// Routine to return the screening paramter kappa^2 in units of keV^2 (kappa^-1 = Debye-Hueckel radius).
double SolarModel::kappa_squared(double r)
{
  return 4.0*pi*alpha_EM/temperature_in_keV(r)*(z2_n(r)+n_e(r))*gsl_pow_3(keV2cm);
}
// Routine to return the number density times Z^2 of ion iz in the zone around the distance r from the centre of the Sun with and without the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf).
double SolarModel::z2_n_iz(double r, int iz) { return gsl_spline_eval(z2_n_iz_lin_interp[iz], r, z2_n_iz_acc[iz]); }
double SolarModel::z2_n(double r){
    if (raffelt_approx == false) {
        double sum = 0.0;
        for (int iz = 0; iz < n_op_elements; iz++) {sum+=SolarModel::z2_n_iz(r,iz);}
        return sum;
    } else {
        return gsl_spline_eval(linear_interp[3], r, accel[3]);
    }
}
// Routine to return the number density of ion iz in the zone around the distance r from the centre of the Sun.
double SolarModel::n_iz(double r, int iz) { return gsl_spline_eval(n_iz_lin_interp[iz], r, n_iz_acc[iz]); }

// Routine to return the electron density in the zone around the distance r from the centre of the Sun with and without the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf).
double SolarModel::n_e(double r) {
    if (raffelt_approx == false) {
        return gsl_spline_eval(linear_interp[1], r, accel[1]);
    } else {
        return gsl_spline_eval(linear_interp[2], r, accel[2]);
    }
}
// Routine to return the plasma freqeuency squared (in keV^2) of the zone around the distance r from the centre of the Sun.
double SolarModel::omega_pl_squared(double r) { return 4.0*pi*alpha_EM/m_electron*n_e(r)*gsl_pow_3(keV2cm); }
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

// TODO: For compatability, introduce a function that does not depend on iz; this then sums all contrib for any opacity code.

// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ff(double omega, double r, int iz) {
  if (omega == 0) {return 0;}
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * pow(alpha_EM*g_aee,2) * gsl_pow_6(keV2cm);
  double u = omega/temperature_in_keV(r);
  double y_red = sqrt(kappa_squared(r)/(2.0*m_electron*temperature_in_keV(r)));
  return prefactor1 * n_e(r)*z2_n_iz(r,iz)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
}
double SolarModel::Gamma_P_ff(double omega, double r) {
  if (omega == 0) {return 0;}
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * pow(alpha_EM*g_aee,2) * gsl_pow_6(keV2cm);
  double u = omega/temperature_in_keV(r);
  double y_red = sqrt(kappa_squared(r)/(2.0*m_electron*temperature_in_keV(r)));
  return prefactor1 * n_e(r)*z2_n(r)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
}


// Calculate the e-e bremsstrahlung contribution; from Eq. (2.18) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ee(double omega, double r) {
  // N.B. "y" and "prefactor2" are different from the "y_red" and "prefactor1" above.
    if (omega == 0) {return 0;}
  const double prefactor2 = (4.0*sqrt(pi)/3.0) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);
  double u = omega/temperature_in_keV(r);
  double y = sqrt(kappa_squared(r)/(m_electron*temperature_in_keV(r)));
  return prefactor2 * n_e(r)*n_e(r)*exp(-u)*aux_function(u,y) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
}

// Calculate the Compton contribution; from Eq. (2.19) of [arXiv:1310.0823]
double SolarModel::Gamma_P_Compton (double omega, double r) {
  if (omega == 0) {return 0;}
  const double prefactor3 = (alpha_EM/3.0) * pow(g_aee/(m_electron),2) * pow(keV2cm,3);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor3 * v*v*n_e(r)/gsl_expm1(u);
}

// Read off interpolated elements
double SolarModel::op_grid_interp_erg (double u, int ite, int jne, int iz) {
    auto key = std::make_pair(ite,jne);
    if (opacity_lin_interp_op[iz].find(key) == opacity_lin_interp_op[iz].end()) {
        if (unavailable_OP.find(key) == unavailable_OP.end()){
            std::cout << "OP data for " << op_elements_simple[iz] << ", ite=" << ite << " and jne=" << jne << " not found!"  << std::endl;
        }
        return 0;
    }
    double result =  gsl_spline_eval(opacity_lin_interp_op[iz].at(key), log(u), opacity_acc_op[iz].at(key));
    if (gsl_isnan(result) == true) {return 0;}
    return result;
};
double SolarModel::tops_grid_interp_erg (double erg, float T, float rho) {
  auto key = std::make_pair(T,rho);
  if (opacity_lin_interp_tops.find(key) == opacity_lin_interp_tops.end()) {
      std::cout << "Grid point {" << T << ", " << rho << "} not found" << std::endl;
      return 0;
  }
  double result = gsl_spline_eval(opacity_lin_interp_tops.at(key), erg, opacity_acc_tops.at(key));
  if (gsl_isnan(result) == true) {return 0;}
  return result;
};
double SolarModel::opas_grid_interp_erg (double erg, double r) {
    if (opacity_lin_interp_opas.find(r) == opacity_lin_interp_opas.end()) {
        std::cout << "OPAS data for R=" <<r << " not found!"  << std::endl;
        return 0;
    }
    double result = gsl_spline_eval(opacity_lin_interp_opas.at(r), erg, opacity_acc_opas.at(r));
    if (gsl_isnan(result) == true) {return 0;}
    return result;
};
  //double result = interp1d(np.log(s[:,0]),s[:,1],bounds_error=False,fill_value=0,kind='linear')(np.log(u))
//  re
//}

double SolarModel::opacity_table_interpolator_op2 (double omega, double r, int iz) {
  // Need temperature in Kelvin
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
  double result = (1.0-t1)*(1.0-t2)*op_grid_interp_erg(u1,ite1,jne1,iz) + (1.0-t1)*t2*op_grid_interp_erg(u1,ite1,jne2,iz)
                + t1*(1.0-t2)*op_grid_interp_erg(u2,ite2,jne1,iz) + t1*t2*op_grid_interp_erg(u2,ite2,jne2,iz);

  if (result < 0) {
    std::cout << "ERROR! Negative opacity!" << std::endl;
    std::cout << "Test 1: " << u1 << " " << u2 << " | " << ite1 << " " << ite2 << " | " << jne1 << " " << jne2 << std::endl;
    std::cout << "Test 2: " << t1 << " " << t2 << " | " << op_grid_interp_erg(u1,ite1,jne1,iz) << " " << op_grid_interp_erg(u1,ite2,jne2,iz) << " | " << result << std::endl;
  };
  return result;
};
double SolarModel::opacity_table_interpolator_op (double omega, double r, int iz) {
  // Need temperature in Kelvin
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
  double result = pow(pow(op_grid_interp_erg(u1,ite1,jne1,iz),1.0-t2)*pow(op_grid_interp_erg(u1,ite1,jne2,iz),t2),1.0-t1)*  pow(pow(op_grid_interp_erg(u2,ite2,jne1,iz),1.0-t2)*pow(op_grid_interp_erg(u2,ite2,jne2,iz),t2),t1);
  if (result < 0) {
    std::cout << "ERROR! Negative opacity!" << std::endl;
    std::cout << "Test 1: " << u1 << " " << u2 << " | " << ite1 << " " << ite2 << " | " << jne1 << " " << jne2 << std::endl;
    std::cout << "Test 2: " << t1 << " " << t2 << " | " << op_grid_interp_erg(u1,ite1,jne1,iz) << " " << op_grid_interp_erg(u1,ite2,jne2,iz) << " | " << result << std::endl;
  };
  return result;
};
double SolarModel::opacity_table_interpolator_tops (double omega, double r) {
  double temperature = temperature_in_keV(r);
  double rho = density(r);
  int lenT = tops_temperatures.size();
  int lenrho = tops_densities.size();
  if ((rho <= tops_densities[0]) || (temperature<=tops_temperatures[0]))  { return 0; }
  float Tlow,Tup;
  for (int k = 0; k < lenT; k++) {
        if (tops_temperatures[k] < temperature) {
            Tlow = tops_temperatures[k];
            Tup = tops_temperatures[k+1];
        }
  }
  float rholow, rhoup;
  for (int k = 0; k < lenrho; k++) {
        if (tops_densities[k] < rho) {
            rholow = tops_densities[k];
            rhoup = tops_densities[k+1];
        }
  }
  double t1 = (temperature-double(Tlow))/double(Tup-Tlow);
  double t2 = (rho-double(rholow))/double(rhoup-rholow);
  double result = pow(pow(tops_grid_interp_erg(omega,Tlow,rholow),1.0-t2)*pow(tops_grid_interp_erg(omega,Tlow,rhoup),t2),1.0-t1)*   pow(pow(tops_grid_interp_erg(omega,Tup,rholow),1.0-t2)*pow(tops_grid_interp_erg(omega,Tup,rhoup),t2),t1);
  if (result < 0) {
    std::cout << "ERROR! Negative opacity!" << std::endl;
  };
  return result;
};
double SolarModel::opacity_table_interpolator_opas (double omega, double r) {
    if (r > opas_radii.back()) {return 0;}
    int lenR = opas_radii.size();
    double Rlow, Rup;
    for (int k = 0; k < lenR; k++) {
          if (opas_radii[k] < r) {
              Rlow = opas_radii[k];
              Rup = opas_radii[k+1];
          }
    }
    double Tup = temperature_in_keV(Rlow);
    double Tlow = temperature_in_keV(Rup);
    double temperature = temperature_in_keV(r);
    double t1 = (temperature-Tlow)/(Tup-Tlow);
    double result = pow(opas_grid_interp_erg(omega*1000.0,Rup),1.0-t1)*pow(opas_grid_interp_erg(omega*1000.0,Rlow),t1);
    if (result < 0) {
      std::cout << "ERROR! Negative opacity!" << std::endl;
    };
    return result;
}

double SolarModel::opacity_element (double omega, double r, int iz) {
    gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);
    if (opcode != OP) {
        std::cout << "Warning: Chosen opacity code does not provide opacities for indivdual elements" << std::endl;
        return 0;
    }
    const double prefactor4 = a_Bohr*a_Bohr*(keV2cm);
    double u = omega/temperature_in_keV(r);
    double result = prefactor4*n_iz(r, iz)*opacity_table_interpolator_op(omega, r, iz)*(-gsl_expm1(-u));
    gsl_set_error_handler (old_handler);
    return result;
};
double SolarModel::opacity (double omega, double r){
    gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);
    double result = 0;
    if (opcode == OP) {
        double sum = 0;
        for (int iz = 0; iz < n_op_elements; iz++) { sum += opacity_element(omega,r,iz); };
        result =  sum;
    }
    if ((opcode == LEDCOP) || (opcode == ATOMIC)){
        result = opacity_table_interpolator_tops(omega, r)*density(r)*keV2cm;
    }
    if (opcode == OPAS) {
        result = opacity_table_interpolator_opas(omega, r)*density(r)*keV2cm;
    }
    gsl_set_error_handler (old_handler);
    return result;
}
double SolarModel::Gamma_P_element (double omega, double r, int iz) {
  const double prefactor5 = 0.5*g_aee*g_aee/(4.0*pi*alpha_EM);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor5*v*v*opacity_element(omega,r,iz)/gsl_expm1(u);
}
double SolarModel::Gamma_P_opacity (double omega, double r) {
  const double prefactor5 = 0.5*g_aee*g_aee/(4.0*pi*alpha_EM);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor5*v*v*opacity(omega,r)/gsl_expm1(u);
}
double SolarModel::Gamma_P_Primakoff (double erg, double r) {
  // N.B. gagg = 10^-16 keV^-1 = 10^-19 eV^-1
  if (erg == 0) {return 0;}
  const double prefactor6 = g_agg*g_agg/(32.0*pi);

  // Get kappa_s^2, omega_plasma^2 and the temperature.
  double ks_sq = kappa_squared(r);
  double w_pl_sq = omega_pl_squared(r);
  double T_in_keV = temperature_in_keV(r);

  // Calculate the flux.
  double x = 4.0*(erg*erg)/ks_sq;
  if (w_pl_sq/(erg*erg) > 1.0) {return 0;}
  double phase_factor = 2.0*sqrt(1.0 - w_pl_sq/(erg*erg))/gsl_expm1(erg/T_in_keV);
  double rate = (ks_sq*T_in_keV)*((1.0 + 1.0/x)*gsl_log1p(x) - 1.0);

  return  prefactor6*phase_factor*rate;
};
