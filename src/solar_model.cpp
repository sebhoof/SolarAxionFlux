// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#include "solar_model.hpp"

// Constructors
SolarModel::SolarModel() : opcode(OP) {} // N.B. We don't need dummy memory allocation for GSL since destructor checks if the vectors containing them are empty

SolarModel::SolarModel(std::string path_to_model_file, opacitycode opcode_tag, bool set_raffelt_approx) : opcode(opcode_tag) {
  std::string path_to_data, model_file_name;
  locate_data_folder(path_to_model_file, path_to_data, model_file_name);
  if ((opcode_tag != OP) && (model_file_name != "SolarModel_AGSS09.dat")) {
    std::cout << "WARNING. The chosen opacity code is only compatible with the solar model AGSS09." << std::endl;
    std::cout << "         Results will be inconsistent." << std::endl;
  }
  // Set whether to use approximations from https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf equation 16 a or alternatively sum over all elements assuming full ionisation
  raffelt_approx = set_raffelt_approx;
  solar_model_name = model_file_name;
  data = ASCIItableReader(path_to_model_file);
  int pts = data.getnrow();
  // Terminate if number of columns is wrong; i.e. the wrong solar model file format.
  int n_cols = data.getncol();
  num_tracked_isotopes = n_cols-6;
  if (n_cols == 35) {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_C13", "X_N14", "X_N15", "X_O16", "X_O17", "X_O18", "X_Ne", "X_Na", "X_Mg", "X_Al", "X_Si", "X_P", "X_S", "X_Cl", "X_Ar",
                     "X_K", "X_Ca", "X_Sc", "X_Ti", "X_V", "X_Cr", "X_Mn", "X_Fe", "X_Co", "X_Ni");
    tracked_isotopes = { {"H",1}, {"He",4}, {"He",3}, {"C",12}, {"C",13}, {"N",14}, {"N",15}, {"O",16}, {"O",17}, {"O",18}, {"Ne",0}, {"Na",0}, {"Mg",0}, {"Al",0}, {"Si",0}, {"P",0}, {"S",0}, {"Cl",0}, {"Ar",0}, {"K",0}, {"Ca",0}, {"Sc",0},
                         {"Ti",0}, {"V",0}, {"Cr",0}, {"Mn",0}, {"Fe",0}, {"Co",0}, {"Ni",0} };
  } else if (n_cols == 12) {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_N14", "X_O16");
    tracked_isotopes = { {"H",1}, {"He",4}, {"He",3}, {"C",12}, {"N",14}, {"O",16} };
  } else {
    //terminate_with_error("ERROR! Solar model file '"+path_to_model_file+"' not compatible with this code (wrong number of columns)!");
    std::string err_msg = "Solar model file '"+path_to_model_file+"' not compatible with this code (wrong number of columns)!";
    throw XSanityCheck(err_msg);
  }

  // Initialise isotope-index map.
  for (int j = 0; j < num_tracked_isotopes; j++) { isotope_index_map[tracked_isotopes[j]] = j; }

  // Extract the radius from the files (in units of the solar radius).
  r_lo = data["radius"][0];
  r_hi = data["radius"][pts-1];

  // Extract the temperature (has to be converted into keV), density, electron density, and ion density * charge^2  from the files.
  // Initialise necessary variables for the screening scale calculation.
  std::vector<double> temperature, n_e, n_e_Raff, density, n_total, z2_n_total, alpha;
  std::vector<std::vector<double>> n_isotope (num_tracked_isotopes);
  std::vector<std::vector<double>> z2_n_isotope (num_tracked_isotopes);
  std::vector<std::vector<double>> n_op_element (num_op_elements);
  // Multiplicative factor: (4 pi alpha_EM / atomic_mass_unit) x (1 g/cm^3) in units of keV^3
  const double factor = 4.0*pi*alpha_EM*gsl_pow_3(keV2cm)/((1.0E+9*eV2g)*atomic_mass_unit);
  // Linearly extrapolate the data in the solar model file to r = 0 if necessary.
  if (r_lo > 0) {
    double r0 = data["radius"][0], r1 = data["radius"][1];

    for (int i = 0; i < n_cols; i++) {
      double intercept = (r0*data[i][1]-r1*data[i][0])/(r0-r1);
      //data[i].insert(data[i].begin(), intercept);
      data.prepend_data(intercept, i);
    }
    r_lo = 0;
    pts += 1;
  }
  num_interp_pts = pts;

  const double* radius = &data["radius"][0];
  // Calculate the necessary quantities and store them internally
  for (int i = 0; i < pts; i++) {
    // Temperature: already given in the Solar model; convert to keV
    temperature.push_back((1.0E-3*K2eV)*data["temperature"][i]);

    // Density: already given in the Solar model
    density.push_back(data["rho"][i]);

    // Electron density: currently 2 options
    // (1) Raffelt's approximation (assuming only fully ionised hydrogen and helium, Y = 1-X (raffelt_approx = true)
    // (2) Assuming full ionisation and summing over all elements (raffelt_approx = false)
    double rhorel = data["rho"][i]/((1.0E+9*eV2g)*atomic_mass_unit); // Density in atoms per cm^3
    double ne_Raff = 0.5 * (1.0 + data["X_H1"][i]) * data["rho"][i] /(atomic_mass_unit*eV2g*1.0E+9); // Electron density for Raffelt's approx
    n_e_Raff.push_back(ne_Raff);
    double ne = 0.0;
    // Electron density from summing over all elements (full ionisation)
    for (int j = 0; j < num_tracked_isotopes; j++) {
      // Isotope corresponds to column (index in 'tracked_isotopes' + 6) in data!
      Isotope isotope = tracked_isotopes[j];
      ne += data[j+6][i] * isotope.z_val() / atomic_weight(isotope);
    }
    n_e.push_back(ne*rhorel);

    // Ion number density and ion number density weighted by charge^2 from summing over all elements (full ionisation) and individual element densities
    double n_total_val = 0.0, z2_n_total_val = 0.0;
    for (int j = 0; j < num_tracked_isotopes; j++) {
      Isotope isotope = tracked_isotopes[j];
      double n = data[j+6][i] * rhorel / atomic_weight(isotope);
      n_total_val += n;
      n_isotope[j].push_back(n);
      double z2_n = gsl_pow_2(isotope.z_val()) * n;
      z2_n_total_val += z2_n;
      z2_n_isotope[j].push_back(z2_n);
    }
    n_total.push_back(n_total_val);
    z2_n_total.push_back(z2_n_total_val);

    // Calculate number density for an element from the OP Code
    for (int k = 0; k < num_op_elements; k++) {
      double temp = 0.0;
      std::string element = op_element_names[k];
      for (int j = 0; j < num_tracked_isotopes; j++) {
        Isotope isotope = tracked_isotopes[j];
        if (isotope.get_element_name() == element) { temp += data[j+6][i] * rhorel / atomic_weight(isotope); }
      }
      n_op_element[k].push_back(temp);
    }
  }

  // Set up the interpolating functions quantities so far
  accel.resize(7);
  linear_interp.resize(7);
  init_numbered_interp(0, radius, &temperature[0]); // Temperature
  init_numbered_interp(1, radius, &n_e[0]); // n_e from summing over elements (full ionisation)
  init_numbered_interp(2, radius, &n_e_Raff[0]); // n_e from Raffelt
  init_numbered_interp(3, radius, &density[0]); // Density
  init_numbered_interp(4, radius, &n_total[0]); // Number density of all ions (currently not used)
  init_numbered_interp(5, radius, &z2_n_total[0]); // Z^2 x number density (full ionisation)

  // Quantities depending on specfific isotope or element
  n_isotope_acc.resize(num_tracked_isotopes);
  n_isotope_lin_interp.resize(num_tracked_isotopes);
  z2_n_isotope_acc.resize(num_tracked_isotopes);
  z2_n_isotope_lin_interp.resize(num_tracked_isotopes);
  for (int j = 0; j < num_tracked_isotopes; j++) {
    init_interp(n_isotope_acc[j], n_isotope_lin_interp[j], radius, &n_isotope[j][0]); // Ion density for each isotope
    init_interp(z2_n_isotope_acc[j], z2_n_isotope_lin_interp[j], radius, &z2_n_isotope[j][0]); // Ion density weighted by charge^2 for each isotope (full ionisation)
  }
  // N.B. The maps don't need to be resized etc. since the access operator '[]' will silently create the map if it is missing.
  for (int k = 0; k < num_op_elements; k++) {
      std::string element = op_element_names[k];
      init_interp(n_element_acc[element], n_element_lin_interp[element], radius, &n_op_element[k][0]); // Ion density for each element (= summed isotopes with same charge)
  }

  // Read squared ionisation from ionisation tables
  for (int j = 0; j < op_grid_size; j++){
      std::string op_filename = path_to_data+"ionisation_tables/ionisation_table_"+std::to_string(op_grid[j][0])+"_"+std::to_string(op_grid[j][1])+".dat";
      ASCIItableReader ion_data = ASCIItableReader(op_filename);
      ion_data.setcolnames("atomic number", "ionisation", "ionisationsqr");
      std::map<std::string, double> temp_ionisationsqr;
      for (int k = 0; k < num_op_elements; k++) {
          std::string element = op_element_names[k];
          temp_ionisationsqr[element] = ion_data["ionisationsqr"][k];
      }
      auto pr = std::make_pair(op_grid[j][0], op_grid[j][1]);
      element_ionisationsqr[pr] = temp_ionisationsqr;
  }

  // Opacity tables setup for interpolating functions (only for chosen opacity code)
  // Do we use OP opacities?
  if (opcode == OP) {
    for (int k = 0; k < num_op_elements; k++) {
      std::string element = op_element_names[k];
      // Initialise grid values
      std::map<std::pair<int,int>, gsl_interp_accel*> temp_acc;
      std::map<std::pair<int,int>, gsl_spline*> temp_interp;
      for (int j = 0; j < op_grid_size; j++) {
        std::string op_filename = path_to_data+"opacity_tables/OP/opacity_table_"+element+"_"+std::to_string(op_grid[j][0])+"_"+std::to_string(op_grid[j][1])+".dat";
        ASCIItableReader op_data = ASCIItableReader(op_filename);

        // Determine the number of interpolated mass values.
        int op_pts = op_data[0].size();
        auto pr = std::make_pair(op_grid[j][0], op_grid[j][1]);
        temp_acc[pr] = gsl_interp_accel_alloc();
        temp_interp[pr] = gsl_spline_alloc (gsl_interp_linear, op_pts);
        const double* omega = &op_data[0][0];
        const double* s = &op_data[1][0];
        gsl_spline_init (temp_interp[pr], omega, s, op_pts);
      }
      opacity_acc_op[element] = temp_acc;
      opacity_lin_interp_op[element] = temp_interp;
    }
  }
  //  Do we use LEDCOP & ATOMIC (both TOPS) opacitites?
  //std::string name;
  if (opcode == LEDCOP){
      tops_grid = ledcop_grid;
      tops_temperatures = ledcop_temperatures;
      tops_densities = ledcop_densities;
      //name = "LEDCOP";
  }
  if (opcode == ATOMIC){
        tops_grid = atomic_grid;
        tops_temperatures = atomic_temperatures;
        tops_densities = atomic_densities;
        //name = "ATOMIC";
  }
  if ((opcode == LEDCOP) || (opcode == ATOMIC)) {
    for (int j = 0; j < tops_grid.size(); j++){
      std::stringstream Tstream;
      std::stringstream rhostream;
      Tstream << std::fixed << std::setprecision(3) << tops_grid[j][0];
      rhostream << std::fixed << std::setprecision(3) << tops_grid[j][1];
      //std::string tops_filename = path_to_data+"opacity_tables/"+name+"/T"+Tstream.str()+"Rho"+rhostream.str()+".dat";
      std::string tops_filename = path_to_data+"opacity_tables/"+get_opacitycode_name()+"/T"+Tstream.str()+"Rho"+rhostream.str()+".dat";
      ASCIItableReader tops_data = ASCIItableReader(tops_filename);
      // Determine the number of interpolated energy values.
      int tops_pts = tops_data[0].size();
      auto pr = std::make_pair(tops_grid[j][0], tops_grid[j][1]);
      opacity_acc_tops[pr] = gsl_interp_accel_alloc();
      opacity_lin_interp_tops[pr] = gsl_spline_alloc(gsl_interp_linear, tops_pts);
      const double* omega = &tops_data[0][0];
      const double* s = &tops_data[1][0];
      gsl_spline_init (opacity_lin_interp_tops[pr], omega, s, tops_pts);
    }
  }

  //  Do we use OPAS opacities?
  if (opcode == OPAS) {
    for (int j = 0 ; j < opas_radii.size(); j++) {
      std::stringstream Rstream;
      Rstream << std::fixed << std::setprecision(2) << opas_radii[j];
      std::string opas_filename = path_to_data+"opacity_tables/OPAS/R"+Rstream.str()+".dat";
      ASCIItableReader opas_data = ASCIItableReader(opas_filename);
      // Determine the number of interpolated energy values.
      int opas_pts = opas_data[0].size();
      double rad = opas_radii[j];
      opacity_acc_opas[rad] = gsl_interp_accel_alloc();
      opacity_lin_interp_opas[rad] = gsl_spline_alloc(gsl_interp_linear, opas_pts);
      const double* omega = &opas_data[0][0];
      const double* s = &opas_data[1][0];
      gsl_spline_init (opacity_lin_interp_opas[rad], omega, s, opas_pts);
    }
  }

  // Calculate alpha values from tables
  if (alpha_available.find(solar_model_name) != alpha_available.end()) {
     data_alpha = ASCIItableReader(path_to_data+"alpha_tables/alpha"+model_file_name.substr(model_file_name.find("_")));
  } else {
     data_alpha = ASCIItableReader(path_to_data+"alpha_tables/alpha_B16-AGSS09.dat");
  }
  int n_cols_alpha = data_alpha.getncol();
  data_alpha.setcolnames("radius", "alpha");
  int pts_alpha = data_alpha.getnrow();
  const double* radius_alpha = &data_alpha["radius"][0];
  for (int i = 0; i < pts_alpha; i++) {alpha.push_back(data_alpha["alpha"][i]);}
  accel[6] = gsl_interp_accel_alloc ();
  linear_interp[6] = gsl_spline_alloc (gsl_interp_linear, pts_alpha);
  const double* alpha_vals = &alpha[0];
  gsl_spline_init (linear_interp[6], radius_alpha, alpha_vals, pts_alpha);

  // All done! Set Solar model class to be correctly initialised
  initialisation_status = true;
}

SolarModel::SolarModel(std::string path_to_model_file, std::string opcode_name, bool set_raffelt_approx) : SolarModel::SolarModel(path_to_model_file, opacitycode_tag.at(opcode_name), set_raffelt_approx) {}

// Class destructor
SolarModel::~SolarModel() {
  for (auto interp : linear_interp) { gsl_spline_free(interp); }
  for (auto interp : n_isotope_lin_interp) { gsl_spline_free(interp); }
  for (auto interp : z2_n_isotope_lin_interp) { gsl_spline_free(interp); }
  for (auto map : n_element_lin_interp) { gsl_spline_free(map.second); }
  for (auto map_1 : opacity_lin_interp_op) {
    for (auto map_2 : map_1.second) { gsl_spline_free(map_2.second); }
  }
  for (auto map : opacity_lin_interp_tops) { gsl_spline_free(map.second); }
  for (auto map : opacity_lin_interp_opas) { gsl_spline_free(map.second); }
  for (auto acc : accel) { gsl_interp_accel_free(acc); }
  for (auto acc : n_isotope_acc) { gsl_interp_accel_free(acc); }
  for (auto acc : z2_n_isotope_acc) { gsl_interp_accel_free(acc); }
  for (auto map : n_element_acc) { gsl_interp_accel_free(map.second); }
  for (auto map_1 : opacity_acc_op) {
    for (auto map_2 : map_1.second) { gsl_interp_accel_free(map_2.second); }
  }
  for (auto map : opacity_acc_tops) { gsl_interp_accel_free(map.second); }
  for (auto map : opacity_acc_opas) { gsl_interp_accel_free(map.second); }
}

// Move assignment operator
SolarModel& SolarModel::operator=(SolarModel &&src) {
  if (this != &src) {
    std::swap(initialisation_status,src.initialisation_status);
    std::swap(raffelt_approx,src.raffelt_approx);
    std::swap(solar_model_name,src.solar_model_name);
    std::swap(data,src.data);
    std::swap(data_alpha,src.data_alpha);
    std::swap(tracked_isotopes,src.tracked_isotopes);
    std::swap(num_interp_pts,src.num_interp_pts);
    std::swap(r_lo,src.r_lo);
    std::swap(r_hi,src.r_hi);
    std::swap(num_tracked_isotopes,src.num_tracked_isotopes);
    std::swap(isotope_index_map,src.isotope_index_map);
    std::swap(tracked_isotopes,src.tracked_isotopes);
    std::swap(accel,src.accel);
    std::swap(linear_interp,src.linear_interp);
    std::swap(opacity_acc_op,src.opacity_acc_op);
    std::swap(opacity_lin_interp_op,src.opacity_lin_interp_op);
    std::swap(tops_grid,src.tops_grid);
    std::swap(tops_temperatures,src.tops_temperatures);
    std::swap(tops_densities,src.tops_densities);
    std::swap(opacity_acc_tops,src.opacity_acc_tops);
    std::swap(opacity_lin_interp_tops,src.opacity_lin_interp_tops);
    std::swap(opacity_acc_opas,src.opacity_acc_opas);
    std::swap(opacity_lin_interp_opas,src.opacity_lin_interp_opas);
    std::swap(n_isotope_acc,src.n_isotope_acc);
    std::swap(n_isotope_lin_interp,src.n_isotope_lin_interp);
    std::swap(z2_n_isotope_acc,src.z2_n_isotope_acc);
    std::swap(z2_n_isotope_lin_interp,src.z2_n_isotope_lin_interp);
    std::swap(n_element_acc,src.n_element_acc);
    std::swap(n_element_lin_interp,src.n_element_lin_interp);
    std::swap(element_ionisationsqr,src.element_ionisationsqr);
  }
  return *this;
}

// Initialise the interpolators
void SolarModel::init_interp(gsl_interp_accel*& acc, gsl_spline*& interp, const double* x, const double* y) {
  acc = gsl_interp_accel_alloc();
  interp = gsl_spline_alloc(gsl_interp_linear, num_interp_pts);
  gsl_spline_init(interp, x, y, num_interp_pts);
}

void SolarModel::init_numbered_interp(const int index, const double* x, const double* y) { init_interp(accel[index], linear_interp[index], x, y); }

// Isotope index lookup
int SolarModel::lookup_isotope_index(Isotope isotope) { return isotope_index_map.at(isotope); }

// Routine to return the various Solar quantities as a function of radius (see hpp file)
double SolarModel::temperature_in_keV(double r) { return gsl_spline_eval(linear_interp[0], r, accel[0]); }
double SolarModel::density(double r) { return gsl_spline_eval(linear_interp[3], r, accel[3]); }
double SolarModel::kappa_squared(double r) { return 4.0*pi*alpha_EM/temperature_in_keV(r)*(z2_n(r)+n_electron(r))*gsl_pow_3(keV2cm); }
double SolarModel::alpha(double r) {
  if (alpha_available.find(solar_model_name) != alpha_available.end()) {
    return gsl_spline_eval(linear_interp[6], r, accel[6]);
  } else {
    double result = 0;
    for (int k = 2; k < num_op_elements; k++) {
      std::string element = op_element_names[k];
      result += mass_fraction(r, element) / metallicity(r) * ionisationsqr_element(r, element) / atomic_weight({element,0});
    }
    return result;
  }
}
double SolarModel::n_element(double r, std::string element) { return gsl_spline_eval(n_element_lin_interp.at(element), r, n_element_acc.at(element)); }
double SolarModel::mass_fraction(double r, std::string element){ return n_element(r,element)*atomic_weight({element,0})*(1.0E+9*eV2g)*atomic_mass_unit/density(r); }
double SolarModel::z2_n_iz(double r, int isotope_index) { return gsl_spline_eval(z2_n_isotope_lin_interp[isotope_index], r, z2_n_isotope_acc[isotope_index]); }
// N.B. Convenience function below (may be slow for many calls!)
double SolarModel::z2_n_iz(double r, Isotope isotope) { int isotope_index = lookup_isotope_index(isotope); return z2_n_iz(r, isotope_index); }
double SolarModel::z2_n(double r) {
  if (alpha_available.find(solar_model_name) != alpha_available.end()) {
    return (mass_fraction(r,"H") + mass_fraction(r,"He") + alpha(r)*metallicity(r)) * density(r)/((1.0E+9*eV2g)*atomic_mass_unit);
  } else {
    return gsl_spline_eval(linear_interp[5], r, accel[5]);  // full ionisation
  }
}
double SolarModel::n_iz(double r, int isotope_index) { return gsl_spline_eval(n_isotope_lin_interp[isotope_index], r, n_isotope_acc[isotope_index]); }
// N.B. Convenience function below (may be slow for many calls!)
double SolarModel::n_iz(double r, Isotope isotope) { int isotope_index = lookup_isotope_index(isotope); return n_iz(r, isotope_index); }
double SolarModel::n_electron(double r) {
  if (raffelt_approx == false) {
    return gsl_spline_eval(linear_interp[1], r, accel[1]);
  } else {
    return gsl_spline_eval(linear_interp[2], r, accel[2]);
  }
}
double SolarModel::metallicity(double r){ return 1.0 - mass_fraction(r, "H") - mass_fraction(r, "He"); }
// Routine to return the plasma freqeuency squared (in keV^2) of the zone around the distance r from the centre of the Sun.
double SolarModel::omega_pl_squared(double r) { return 4.0*pi*alpha_EM/m_electron*n_electron(r)*gsl_pow_3(keV2cm); }

// Opacity correction factor
double SolarModel::apply_opacity_correction_factor(double r) {
  static double temp_0 = temperature_in_keV(r_lo);
  const double r_cz_theo = radius_cz;
  const double temp_cz = temperature_in_keV(r_cz_theo);
  double result = 1.0;
  if (r < r_cz_theo) {
    result *= ( 1.0 + opacity_correction_a + opacity_correction_b * log10(temp_0/temperature_in_keV(r)) / log10(temp_0/temp_cz) );
    result = std::max(result, 0.0);
  }
  return result;
}

// Opacity for individual isotope (only possible for OP)
double SolarModel::opacity_element(double omega, double r, std::string element) {
  const double prefactor4 = a_Bohr*a_Bohr*(keV2cm);

  gsl_error_handler_t* old_handler = gsl_set_error_handler(&my_gsl_handler); // GSL error handler modified to avoid boundary error (fill value = 0)
  //terminate_with_error_if(opcode != OP, "ERROR! Chosen opacity code does not provide opacities for indivdual elements.");
  if (opcode != OP) {
    std::string err_msg = "The chosen opacity code ("+get_opacitycode_name()+") does not provide opacities for indivdual elements.";
    throw XUnsupportedOption(err_msg);
  }

  double u = omega/temperature_in_keV(r);
  double result = prefactor4*n_element(r, element)*opacity_table_interpolator_op(omega, r, element)*(-gsl_expm1(-u));
  gsl_set_error_handler(old_handler); // reset the GSL error handler

  return result*apply_opacity_correction_factor(r);
}
// N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
double SolarModel::opacity_element(double omega, double r, Isotope isotope) { return opacity_element(omega, r, isotope.get_element_name()); }

// Opacity for total solar mixture
double SolarModel::opacity(double omega, double r) {
  double result = 0.0;

  if (opcode == OP) {
    for (int k = 0; k < num_op_elements; k++) { result += opacity_element(omega, r, op_element_names[k]); }
  } else if ((opcode == LEDCOP) || (opcode == ATOMIC)) {
    result = opacity_table_interpolator_tops(omega, r)*density(r)*keV2cm;
  } else if (opcode == OPAS) {
    result = opacity_table_interpolator_opas(omega, r)*density(r)*keV2cm;
  }

  return result*apply_opacity_correction_factor(r);
}



double SolarModel::bfield(double r) {
  const double b_rad = 3.0e3/(1.0e6*eV2T); // in keV^2
  const double b_tach = 40.0/(1.0e6*eV2T); // in keV^2
  const double b_outer = 3.5/(1.0e6*eV2T); // in keV^2
  const double lambda = 10.0*radius_cz + 1.0;
  const double lambda_factor = (1.0 + lambda)*pow(1.0 + 1.0/lambda, lambda);

  double result = 0.0;
  if (r < radius_cz+size_tach) {
    double x = gsl_pow_2(r/radius_cz);
    if (x < 1.0) { result += b_rad*(1.0+bfield_correction_rad)*lambda_factor*x*pow(1.0-x, lambda); }
    double y = gsl_pow_2((r - radius_cz)/size_tach);
    if (y < 1.0) { result += b_tach*(1.0+bfield_correction_tach)*(1.0-y); }
  } else {
    double z = gsl_pow_2((r - radius_outer)/size_outer);
    if (z < 1.0) { result += b_outer*(1.0+bfield_correction_outer)*(1.0 - z); }
  }

  return result;
}

// Auxiliary function required for FF and ee contribution
const double abs_prec_aux_fun = 0;
const double rel_prec_aux_fun = 1.0e-4;
//const int int_method_aux_fun = 5;
const int int_space_size_aux_fun = 1e6;
struct integrand_params_aux_fun { double u; double y; };

double integrand(double x, void * params) {
  // Retrieve parameters and other integration variables.
  struct integrand_params_aux_fun * p = (struct integrand_params_aux_fun *)params;
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
  integrand_params_aux_fun p = { u, y };
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_aux_fun);
  double result, error;
  gsl_function f;
  f.function = &integrand;
  f.params = &p;
  //gsl_set_error_handler_off();
  //gsl_integration_qag (&f, sqrt(x*x + u) - x, sqrt(x*x + u) - x, 0.1*abs_prec, 0.1*rel_prec, 1e6, method, w, &result, &error);
  gsl_integration_qagiu(&f, 0, abs_prec_aux_fun, rel_prec_aux_fun, int_space_size_aux_fun, w, &result, &error);
  //printf ("GSL status: %s\n", gsl_strerror (status));
  //gsl_integration_qags(&F, rad, rmax, 1e-1*abs_prec, 1e-1*rel_prec, 1E6, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

struct integrand_params_rosseland { SolarModel* s; double r;};
double rosseland_integrand(double omega, void * params){
    const double prefactor = 15.0/(4.0*gsl_pow_4(pi));
    struct integrand_params_rosseland * p = (struct integrand_params_rosseland *)params;
    double r = (p->r);
    // double opac = p->s->opacity(u * p->s->temperature_in_keV(r), r);
    double opac = p->s->opacity(omega, r);
    double u = omega/(p->s->temperature_in_keV(r));
    //double opac = 1.0;
    if (opac == 0) {
        return 0;
    }
    //double R = 15.0 / (4.0 * gsl_pow_4(pi)) * gsl_pow_4(u)*exp(u) / gsl_pow_2(gsl_expm1(u));
    //double R = 15.0 / (4.0 * gsl_pow_4(pi)) * gsl_pow_4(u)*exp(u) / gsl_pow_2(exp(u) - 1.0);
    double R = prefactor * gsl_sf_exp(u) *gsl_pow_2(u/gsl_sf_exprel(u));
    // double R = 15.0 / (4.0 * gsl_pow_4(pi)) * gsl_pow_3(u)/((1.0 - gsl_sf_exp(-u))*gsl_sf_exprel(u));
    // double R = 15.0 / (4.0 * pi*pi*pi*pi) * u*u*u*u *exp(u) / ((exp(u)-1)*(exp(u)-1));
    return R / opac ;
}
double SolarModel::rosseland_opacity(double r) {
    integrand_params_rosseland p = { this, r };
    double result, error;
    size_t neval;
    gsl_function f;
    f.function = &rosseland_integrand;
    f.params = &p;
    //gsl_integration_workspace * w = gsl_integration_workspace_alloc(int_space_size_aux_fun);
    gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc(int_space_size_aux_fun);
    std::vector<double> relevant_peaks = get_relevant_peaks(0.1, 19.0);
    //gsl_integration_qagiu(&f, 0, abs_prec_aux_fun, rel_prec_aux_fun, int_space_size_aux_fun, w, &result, &error);
    //gsl_integration_qag(&f, 3.0, 5.0, 0, 1.0e-4, int_space_size_aux_fun, 5, w, &result, &error);
    //gsl_integration_qagp(&f, &relevant_peaks[0], relevant_peaks.size(), 0.0, 1.0e-3, int_space_size_aux_fun, w, &result, &error);
    //gsl_integration_qags(&f, 0.1, 19.0, 0, 1.0e-4, int_space_size_aux_fun, w, &result, &error);
    gsl_integration_cquad(&f, 0.1, 19.0, 0.0, 1.0e-4, v, &result, &error, &neval);
    //gsl_integration_workspace_free (w);
    gsl_integration_cquad_workspace_free(v);
    //gsl_integration_qng(&f, 4, 5 , abs_prec_aux_fun, 1.0e-2, &result, &error,&neval);
    return temperature_in_keV(r)/result;
}


// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823] (assuming full ionisation) for one isotope
double SolarModel::Gamma_P_ff(double omega, double r, int isotope_index) {
  if (omega == 0) { return 0; }
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);
  double u = omega/temperature_in_keV(r);
  double y_red = sqrt(kappa_squared(r)/(2.0*m_electron*temperature_in_keV(r)));
  return prefactor1 * n_electron(r)*z2_n_iz(r,isotope_index)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
}

// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823] (assuming full ionisation) for on isotope
// N.B. Convenience function below (may be slow for many calls!)  (assuming full ionisation)
double SolarModel::Gamma_P_ff(double omega, double r, Isotope isotope) { int isotope_index = lookup_isotope_index(isotope); return Gamma_P_ff(omega, r, isotope_index); }

// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823] (assuming full ionisation)
double SolarModel::Gamma_P_ff(double omega, double r) {
  double result = 0;
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);

  if (omega > 0) {
    if (raffelt_approx == false) {
      // Assume full ionisation and only take H and He
      static int iso_ind_1 = lookup_isotope_index({"H",1}), iso_ind_2 = lookup_isotope_index({"He",3}), iso_ind_3 = lookup_isotope_index({"He",4});
      result = Gamma_P_ff(omega, r, iso_ind_1) + Gamma_P_ff(omega, r, iso_ind_2) + Gamma_P_ff(omega, r, iso_ind_3);
    } else {
      // Contributions from all elements, no full ionisation assumed when alpha available
      double u = omega/temperature_in_keV(r);
      double y_red = sqrt(kappa_squared(r)/(2.0*m_electron*temperature_in_keV(r)));
      result = prefactor1 * n_electron(r)*z2_n(r)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
    }
  }

  return result;
}

// Calculate the e-e bremsstrahlung contribution; from Eq. (2.18) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ee(double omega, double r) {
  // N.B. "y" and "prefactor2" are different from the "y_red" and "prefactor1" above.
  const double prefactor2 = (4.0*sqrt(pi)/3.0) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);
  if (omega > 0) {
    double u = omega/temperature_in_keV(r);
    double y = sqrt(kappa_squared(r)/(m_electron*temperature_in_keV(r)));
    return prefactor2 * gsl_pow_2(n_electron(r)) * exp(-u) * aux_function(u,y) / (omega * sqrt(temperature_in_keV(r)) * pow(m_electron,3.5));
  } else {
    return 0;
  }
}

// Calculate the Compton contribution; from Eq. (2.19) of [arXiv:1310.0823]
double SolarModel::Gamma_P_Compton(double omega, double r) {
  const double prefactor3 = (alpha_EM/3.0) * pow(g_aee/(m_electron),2) * pow(keV2cm,3);
  if (omega > 0) {
    double u = omega/temperature_in_keV(r);
    double v = omega/m_electron;
    return prefactor3 * v*v*n_electron(r)/gsl_expm1(u);
  } else {
    return 0;
  }
}

// opacity contribution from one isotope; first term of Eq. (2.21) of [arXiv:1310.0823]
double SolarModel::Gamma_P_opacity(double omega, double r, std::string element) {
  const double prefactor5 = 0.5*g_aee*g_aee/(4.0*pi*alpha_EM);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor5*v*v*opacity_element(omega,r,element)/gsl_expm1(u);
}

double SolarModel::Gamma_P_opacity(double omega, double r, Isotope isotope) {
  std::string element = isotope.get_element_name();
  return Gamma_P_opacity(omega, r, element);
}

// full opacity contribution; first term of Eq. (2.21) of [arXiv:1310.0823]
double SolarModel::Gamma_P_opacity(double omega, double r) {
  const double prefactor5 = 0.5*g_aee*g_aee/(4.0*pi*alpha_EM);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor5*v*v*opacity(omega,r)/gsl_expm1(u);
}

double SolarModel::Gamma_P_all_electron(double erg, double r) {
  double result = 0;
  if (opcode == OP) {
    double element_contrib = 0.0;
    static int iso_ind_1 = lookup_isotope_index({"H",1}), iso_ind_2 = lookup_isotope_index({"He",3}), iso_ind_3 = lookup_isotope_index({"He",4});
    element_contrib += Gamma_P_ff(erg, r, iso_ind_1) + Gamma_P_ff(erg, r, iso_ind_2) + Gamma_P_ff(erg, r, iso_ind_3);
    for (int k = 2; k < num_op_elements; k++) { element_contrib += Gamma_P_opacity(erg, r, op_element_names[k]); }
    result = element_contrib + Gamma_P_Compton(erg, r) + Gamma_P_ee(erg, r);
  } else if ((opcode == LEDCOP) || (opcode == ATOMIC)) {
    double u = erg/temperature_in_keV(r);
    double reducedCompton = 0.5*(1.0 - 1.0/gsl_expm1(u)) * Gamma_P_Compton(erg, r);
    result = Gamma_P_opacity(erg, r) + reducedCompton + Gamma_P_ee(erg, r);
  } else if (opcode == OPAS) {
    result = Gamma_P_opacity(erg, r);
  } else {
    //terminate_with_error("ERROR! Unkown option for 'opcode' argument. Use OP, LEDCOP, ATOMIC, or OPAS.");
    std::string err_msg = "Unkown option for 'opcode' argument. Use ";
    for (auto it = opacitycode_name.begin(); it != --opacitycode_name.end(); ++it) { err_msg += it->second + ", "; }; err_msg += "or " + (--opacitycode_name.end())->second + ".";
    throw XUnsupportedOption(err_msg);
  }
  return result;
}

double SolarModel::Gamma_P_Primakoff(double erg, double r) {
  const double prefactor6 = g_agg*g_agg/(32.0*pi);
  if (erg > 0) {
    double w_pl_sq = omega_pl_squared(r);
    if (w_pl_sq/(erg*erg) > 1.0) { return 0; }
    double ks_sq = kappa_squared(r);
    double T_in_keV = temperature_in_keV(r);
    double x = 4.0*(erg*erg)/ks_sq;
    double phase_factor = 2.0*sqrt(1.0 - w_pl_sq/(erg*erg))/gsl_expm1(erg/T_in_keV);
    double rate = (ks_sq*T_in_keV)*((1.0 + 1.0/x)*gsl_log1p(x) - 1.0);
    return prefactor6*phase_factor*rate;
  } else {
    return 0;
  }
}

double SolarModel::Gamma_P_LP(double erg, double r) {
  const double prefactor8 = g_agg*g_agg;
  return 0.0;
}
//simplified version
/*
double SolarModel::Gamma_P_TP(double omega, double r) {
  const double prefactor9 = g_agg*g_agg ;
  const double geom_factor = 1.0;  // factor accounting for observers position (1.0 = angular average)
  double u = omega/temperature_in_keV(r);
  double average_b_field_sq = gsl_pow_2(bfield(r)) * 1.0e-12 /3.0;  // in keV^4
  double result = geom_factor * prefactor9 * average_b_field_sq * opacity(omega,r) * exp(-u) * gsl_pow_2(omega) / gsl_pow_2(omega_pl_squared(r));
  return result;
}
*/

//Full version
double SolarModel::Gamma_P_TP(double omega, double r) {
  double u = omega/temperature_in_keV(r);
  //static OneDInterpolator ross_op (SOLAXFLUX_DIR "/data/opacity_tables/arXiv_1601_01930v2_fig10_opacity.txt");
  //double gamma = -gsl_expm1(-u)*ross_op.interpolate(r);
  //double gamma = - rosseland_opacity(r) * gsl_expm1(-u);
  double gamma = -gsl_expm1(-u)*opacity(omega, r);
  double DeltaPsq = gsl_pow_2(0.5*omega_pl_squared(r)/omega);
  double average_b_field_sq = gsl_pow_2(bfield(r))/(3.0*2.0*pi);
  double DeltaTsq = g_agg*g_agg * average_b_field_sq /4.0;
  const double geom_factor = 1.0;  // factor accounting for observers position (1.0 = angular average)
  const double photon_polarization = 2.0;
  double result = geom_factor * photon_polarization * gamma * DeltaTsq / ( (DeltaPsq+gsl_pow_2(0.5*gamma)) * gsl_expm1(u) ) ;
  return result;
}

double SolarModel::Gamma_P_all_photon(double erg, double r) {
  return Gamma_P_Primakoff(erg, r) + Gamma_P_TP(erg, r);
}

// Read off interpolated elements for op, tops and opas
double SolarModel::op_grid_interp_erg(double u, int ite, int jne, std::string element) {
  double result = 0;
  auto grid_position = std::make_pair(ite,jne);

  if (unavailable_OP.find(grid_position) == unavailable_OP.end()) {
    if (opacity_lin_interp_op.find(element) == opacity_lin_interp_op.end()) {
      std::string err_msg = "OP data for element "+element+" does not exist.";
      throw XUnsupportedOption(err_msg);
      //terminate_with_error("ERROR! OP data for element "+element+" does not exist.");
    } else if (opacity_lin_interp_op.at(element).find(grid_position) == opacity_lin_interp_op.at(element).end()) {
      std::cout << "WARNING. OP data for " << element << " at position ite = " << ite << " and jne = " << jne << " does not exist."  << std::endl;
    } else {
      gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_gsl_handler);   // error handler modified to avoid boundary error (fill value = 0)
      result = gsl_spline_eval(opacity_lin_interp_op.at(element).at(grid_position), log(u), opacity_acc_op.at(element).at(grid_position));
      gsl_set_error_handler (old_handler);    //  reset the error handler
      if (gsl_isnan(result)) { result = 0; }
    }
  }
  return result;
}

double SolarModel::tops_grid_interp_erg(double erg, float t, float rho) {
  double result = 0;
  auto key = std::make_pair(t,rho);
  if (opacity_lin_interp_tops.find(key) == opacity_lin_interp_tops.end()) {
    std::cout << "WARNING. Grid point {" << t << ", " << rho << "} not found!" << std::endl;
  } else {
    gsl_error_handler_t* old_handler = gsl_set_error_handler(&my_gsl_handler); // GSL error handler modified to avoid boundary error (fill value = 0)
    result = gsl_spline_eval(opacity_lin_interp_tops.at(key), erg, opacity_acc_tops.at(key));
    gsl_set_error_handler(old_handler); // reset the GSL error handler
    if (gsl_isnan(result)) { return 0; }
  }
  return result;
}

double SolarModel::opas_grid_interp_erg(double erg, double r) {
  if (opacity_lin_interp_opas.find(r) == opacity_lin_interp_opas.end()) {
    std::cout << "WARNING. OPAS data for R = " << r << " not found!"  << std::endl;
    return 0;
  }

  gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_gsl_handler); // GSL error handler modified to avoid boundary error (fill value = 0)
  double result = gsl_spline_eval(opacity_lin_interp_opas.at(r), erg, opacity_acc_opas.at(r));
  gsl_set_error_handler (old_handler); // reset the GSL error handler
  if (gsl_isnan(result) == true) { return 0; }
  return result;
}

//  double linear interpolation on solar grid (currently not used)
// double SolarModel::opacity_table_interpolator_op2(double omega, double r, std::string element) {
//   // Need temperature in Kelvin
//   double temperature = temperature_in_keV(r)/(1.0e-3*K2eV);
//   double ne = n_electron(r);
//   double ite = 40.0*log10(temperature);
//   double jne = 4.0*log10(ne);
//   int ite2 = int(ceil(20.0*log10(temperature))*2);
//   int ite1 = ite2 - 2;
//   // Need omega in Kelvin
//   double u1 = omega/(1.0e-3*K2eV*pow(10,double(ite1)/40.0));
//   double u2 = omega/(1.0e-3*K2eV*pow(10,double(ite2)/40.0));
//   int jne2 = int(ceil(log10(ne)*2)*2);
//   int jne1 = jne2 - 2;
//   double t1 = (ite-double(ite1))/2.0;
//   double t2 = (jne-double(jne1))/2.0;
//   double result = (1.0-t1)*(1.0-t2)*op_grid_interp_erg(u1,ite1,jne1,element) + (1.0-t1)*t2*op_grid_interp_erg(u1,ite1,jne2,element)
//                 + t1*(1.0-t2)*op_grid_interp_erg(u2,ite2,jne1,element) + t1*t2*op_grid_interp_erg(u2,ite2,jne2,element);
//
//   if (result < 0) {
//     std::cout << "ERROR! Negative opacity!" << std::endl;
//     std::cout << "Test 1: " << u1 << " " << u2 << " | " << ite1 << " " << ite2 << " | " << jne1 << " " << jne2 << std::endl;
//     std::cout << "Test 2: " << t1 << " " << t2 << " | " << op_grid_interp_erg(u1,ite1,jne1,element) << " " << op_grid_interp_erg(u1,ite2,jne2,element) << " | " << result << std::endl;
//   }
//   return result;
// }

//  double logarithmic interpolation on solar grid (used for all codes)
double SolarModel::opacity_table_interpolator_op(double omega, double r, std::string element) {
  // Need temperature in Kelvin
  double temperature = temperature_in_keV(r)/(1.0e-3*K2eV);
  double ne = n_electron(r);
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
  double result = pow(pow(op_grid_interp_erg(u1,ite1,jne1,element),1.0-t2)*pow(op_grid_interp_erg(u1,ite1,jne2,element),t2),1.0-t1)*  pow(pow(op_grid_interp_erg(u2,ite2,jne1,element),1.0-t2)*pow(op_grid_interp_erg(u2,ite2,jne2,element),t2),t1);
  if (result < 0) {
    std::string err_msg = "Negative opacity from SolarModel::opacity_table_interpolator_op.";
    throw XSanityCheck(err_msg);
    //std::cout << "ERROR! Negative opacity!" << std::endl;
  }
  return result;
}

//  double logarithmic interpolation for all ionisations
double SolarModel::ionisationsqr_element(double r, std::string element) {
  // Need temperature in Kelvin
  double temperature = temperature_in_keV(r)/(1.0e-3*K2eV);
  double ne = n_electron(r);
  double ite = 40.0*log10(temperature);
  double jne = 4.0*log10(ne);
  int ite2 = int(ceil(20.0*log10(temperature))*2);
  int ite1 = ite2 - 2;
  int jne2 = int(ceil(log10(ne)*2)*2);
  int jne1 = jne2 - 2;
  double t1 = (ite-double(ite1))/2.0;
  double t2 = (jne-double(jne1))/2.0;
  double result = pow(pow(ionisationsqr_grid(ite1,jne1,element),1.0-t2)*pow(ionisationsqr_grid(ite1,jne2,element),t2),1.0-t1)*  pow(pow(ionisationsqr_grid(ite2,jne1,element),1.0-t2)*pow(ionisationsqr_grid(ite2,jne2,element),t2),t1);
  if (result < 0) {
    std::string err_msg = "Negative ionisation.";
    throw XSanityCheck(err_msg);
    //std::cout << "ERROR! Negative ionisation!" << std::endl;
  }
  return result;
}


double SolarModel::opacity_table_interpolator_tops(double omega, double r) {
  double temperature = temperature_in_keV(r);
  double rho = density(r);
  int lenT = tops_temperatures.size();
  int lenrho = tops_densities.size();
  tops_temperatures[0];
  tops_densities[0];
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
  double result = pow(pow(tops_grid_interp_erg(omega,Tlow,rholow),1.0-t2) * pow(tops_grid_interp_erg(omega,Tlow,rhoup),t2),1.0-t1) * pow(pow(tops_grid_interp_erg(omega,Tup,rholow),1.0-t2) * pow(tops_grid_interp_erg(omega,Tup,rhoup),t2),t1);
  if (result < 0) { std::cout << "ERROR! Negative opacity!" << std::endl; }
  return result;
}

double SolarModel::opacity_table_interpolator_opas(double omega, double r) {
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
  if (result >= 0) {
    return result;
  } else {
    std::cout << "ERROR! Negative opacity (" << result << ") at r = " << r << ", omega = " << omega << "! Will set opacity to 0 and continue but this is a bug, so please report it." << std::endl;
    return 0;
  }
}


// Read off ionisation states
double SolarModel::ionisationsqr_grid(int ite, int jne, std::string element) {
  double result = 0.0;
  auto grid_position = std::make_pair(ite,jne);
  if (unavailable_OP.find(grid_position) == unavailable_OP.end()) {
    if (element_ionisationsqr.find(grid_position) == element_ionisationsqr.end()) {
        std::cout << "WARNING. OP Ionisation data for " << element << " at position ite = " << ite << " and jne = " << jne << " does not exist."  << std::endl;
    } else if (element_ionisationsqr.at(grid_position).find(element) == element_ionisationsqr.at(grid_position).end()) {
        terminate_with_error("ERROR! OP  Ionisation data for element "+element+" does not exist.");
    } else  {
      result = element_ionisationsqr.at(grid_position).at(element);
      if (gsl_isnan(result) == true) { return 0; }
    }
  }
  return result;
}

// TODO: Utilise this + the new typedef
SolarModelMemberFn get_SolarModel_function_pointer(std::string interaction_name) {
  SolarModelMemberFn integrand;
  auto iter = map_interaction_name_to_function.find(interaction_name);
  if (iter != map_interaction_name_to_function.end()) {
    integrand = map_interaction_name_to_function.at(interaction_name);
  } else {
    std::string avail_keys = "";
    for (auto& x: map_interaction_name_to_function) { avail_keys += x.first + " "; }
    std::cout << "NON-FATAL ERROR! The interaction '"+interaction_name+"' is not available. Enter one of the valid options below." << std::endl;
    std::cout << avail_keys << std::endl;
    std::string new_interaction_name;
    std::cout << "Enter a valid interaction name here: "; std::cin >> new_interaction_name;
    integrand = get_SolarModel_function_pointer(new_interaction_name);
  }
  return integrand;
}

// Metadata and information
double SolarModel::get_r_lo() { return r_lo; }
double SolarModel::get_r_hi() { return r_hi; }

std::vector<double> SolarModel::get_supported_radii(std::vector<double> radii) {
  std::vector<double> supported_radii;
  auto it = std::min_element(radii.begin(), radii.end());
  double rad_min = *it;
  it = std::max_element(radii.begin(), radii.end());
  double rad_max = *it;
  if ((r_lo > rad_min) || (r_hi < rad_max)) { std::cout << "WARNING. Radii do not agree with min/max radius values available in Solar model! Unsupported radii will be ignored." << std::endl; }
  supported_radii.push_back(r_lo);
  for (auto r = radii.begin(); r !=radii.end(); ++r) { if ((r_lo < *r) && (*r < r_hi)) { supported_radii.push_back(*r); } }
  //supported_radii.push_back(0.95*r_hi);
  return supported_radii;
}

double SolarModel::get_gagg_ref_value_in_inverse_GeV() { return 1.0e6*g_agg; }
double SolarModel::get_gaee_ref_value() { return g_aee; }
std::string SolarModel::get_solaxlib_name_and_version() { return LIBRARY_NAME; }
std::string SolarModel::get_solar_model_name() { return solar_model_name; }
std::string SolarModel::get_opacitycode_name() { return opacitycode_name.at(opcode); }
void SolarModel::set_bfield_correction(double c_rad, double c_tach, double c_outer) {
  bfield_correction_rad = c_rad;
  bfield_correction_tach = c_tach;
  bfield_correction_outer = c_outer;
}
std::vector<double> SolarModel::get_bfield_correction() {
  std::vector<double> result = { bfield_correction_rad, bfield_correction_tach, bfield_correction_outer };
  return result;
}
void SolarModel::set_opacity_correction(double a, double b) { opacity_correction_a = a; opacity_correction_b = b; }
std::vector<double> SolarModel::get_opacity_correction() {
  std::vector<double> result = { opacity_correction_a, opacity_correction_b };
  return result;
}
bool SolarModel::is_initialised() { return initialisation_status; }
