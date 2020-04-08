// ------------
//  DISCLAIMER
// ------------
// The ASCIItableReader and SolarModel classes are based on classes with the same name in GAMBIT, which was published under a BSD license (reproduced below).
// Copyright (c) 2017, The GAMBIT Collaboration All rights reserved.
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//  (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//  (3) Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "utils.hpp"

void terminate_with_error(std::string err_string) {
  std::cerr << err_string << std::endl;
  exit(EXIT_FAILURE);
}

void terminate_with_error_if(bool condition, std::string err_string) {
  if (condition) { terminate_with_error(err_string); };
}

void my_handler(const char * reason, const char * file, int line, int gsl_errno) {
    if (gsl_errno == GSL_EDOM) { }
    else {
        std::cout << reason << " in " << file << " line: " << line << std::endl;
        std::cout << "error number: " << gsl_errno << std::endl;
        abort();
    }
}

// Check if a file exists
bool file_exists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void save_to_file(std::string path, std::vector<std::vector<double>> buffer, std::string comment, bool overwrite) {
  //std::cout << "Saving results to " << path << "..." << std::endl;

  if (file_exists(path)) {
    if (overwrite) {
      std::cout << "WARNING! File " << path << " exists and will be overwritten..." << std::endl;
    } else {
      std::cout << "File " << path << " exists! Now saving to " << path << "_new" << std::endl;
      path += "_new";
    };
  };

  std::ofstream output;
  output.open(path);
  if (comment != "") {
    size_t pos = 0;
    const std::string newline = "\n";
    const std::string comment_newline = "\n# ";
    while ((pos = comment.find(newline, pos)) != std::string::npos) {
     comment.replace(pos, newline.length(), comment_newline);
     pos += comment_newline.length();
   };
   output << "# " << comment << std::endl; };
  int n_cols = buffer.size();
  if (n_cols > 0) {
    int n_rows = buffer[0].size();
    if (n_rows > 0) {
      for (int i=0; i<n_rows; ++i) {
        output << buffer[0][i];
        for (int j=1; j<n_cols; ++j) {
          output << " " << buffer[j][i];
        };
        output << std::endl;
      };
    } else {
      std::cout << "WARNING! The data you are trying to write to " << path << " is empty! Created an empty file." << std::endl;
    };
  };
  output.close();
}

OneDInterpolator::OneDInterpolator() {}

void OneDInterpolator::init(std::vector<double> x, std::vector<double> y, std::string type) {
  // Initialise gsl interpolation routine.
  int pts = x.size();
  const double* x_ptr = &x[0];
  const double* y_ptr = &y[0];
  acc = gsl_interp_accel_alloc ();
  if (type == "cspline") {
    spline = gsl_spline_alloc (gsl_interp_cspline, pts);
  } else if (type == "linear") {
    spline = gsl_spline_alloc (gsl_interp_linear, pts);
  } else {
    terminate_with_error("ERROR! Interpolation type '"+type+"' not known to class OneDInterpolator.\n       Available types: 'linear' and 'cspline'.");
  };
  gsl_spline_init (spline, x_ptr, y_ptr, pts);
  // Get first and last value of the "x" component.
  lo = x.front();
  up = x.back();
}

// Initialiser for the OneDInterpolator class.
OneDInterpolator::OneDInterpolator(std::string file, std::string type) {
  // Check if file exists.
  terminate_with_error_if(not(file_exists(file)), "ERROR! File '"+file+"' not found!");
  // Read numerical values from data file.
  ASCIItableReader tab (file);
  //tab.setcolnames("x", "y"); -> Would raise warning if we try to interpolate data with errors

  init(tab[0], tab[1], type);
}

OneDInterpolator::OneDInterpolator(std::vector<double> x, std::vector<double> y, std::string type) { init(x, y, type); }

// Move assignment operator
OneDInterpolator& OneDInterpolator::operator=(OneDInterpolator&& interp) {
  if(this != &interp) {
    std::swap(acc,interp.acc);
    std::swap(spline,interp.spline);
    std::swap(lo,interp.lo);
    std::swap(up,interp.up);
  }
  return *this;
}

// Destructor
OneDInterpolator::~OneDInterpolator() {
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}

// Routine to access interpolated values.
double OneDInterpolator::interpolate(double x) { return gsl_spline_eval(spline, x, acc); }

// Routines to return upper and lower boundaries of interpolating function.
double OneDInterpolator::lower() { return lo; }
double OneDInterpolator::upper() { return up; }

// ASCII Table Reader originally by Christoph Weniger
int ASCIItableReader::read(std::string filename) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (in.fail()) {
    std::cout << "ERROR! Failed loading: " << filename << std::endl;
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

void Isotope::init(std::string s, int a) {
  element_name = s;
  isotope_a_value = a;
  // Map of element name -> Z value. CAVE: We run into problems if the const map is defined outside of this function.
  const std::map<std::string, int> z_value_map { {"H", 1}, {"He", 2}, {"C", 6}, {"N", 7}, {"O", 8}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17},
                                                 {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28} };
  if (z_value_map.find(s) == z_value_map.end()) {
    std::string avail_keys = "";
    for (auto& p: z_value_map) { avail_keys += p.first + " "; };
    terminate_with_error("ERROR! Element named '"+s+"' not found! Available keys are:\n"+avail_keys);
  } else {
    element_z_value = z_value_map.at(s);
  };
}
Isotope::Isotope(std::string s, int a) { init(s,a); };
Isotope::Isotope(std::pair<std::string,int> p) { init(p.first,p.second); };
// This is for convenience in order to define elements as an Isotope; for now: empty contructor
// TODO: do init(s,-1) here or allow strings like "He_3" etc. Also, new feature: if el_a_value < -1, trigger adding up all values for the same el_name
Isotope::Isotope(std::string s) { };
bool Isotope::operator< (const Isotope& other) const { return (other.name() < element_name) || ( (other.name() == element_name) && (other.a_val() < isotope_a_value) ); }
bool Isotope::operator== (const Isotope& other) const { return ((other.name() == element_name) && (other.a_val() == isotope_a_value)); }
std::string Isotope::name() const { return element_name; };
std::string Isotope::index_name() const { return "X_"+element_name+std::to_string(isotope_a_value); };
int Isotope::a_val() const { return isotope_a_value; };
int Isotope::z_val() const { return element_z_value; };
bool Isotope::same_z(Isotope *isotope) { return element_name == isotope->name(); };

double atomic_weight(Isotope isotope) { return isotope_avg_weight.at(isotope); }


// Constructors
SolarModel::SolarModel() : opcode(OP) {};    // default constructor (not functional)
SolarModel::SolarModel(std::string file, opacitycode set_opcode, bool set_raffelt_approx) : opcode(set_opcode) {
  if ((set_opcode != OP) && (file != "data/SolarModel_AGSS09.dat")) {
      std::cout << "WARNING. The chosen opacity code is only compatible with the solar model AGSS09." << std::endl;
      std::cout << "         Results will be inconsistent." << std::endl;
  }
//  set whether to use approximations from https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf equation 16 a or alternatively sum over all elements assuming full ionisation
  raffelt_approx = set_raffelt_approx;
  solarmodel_name = file;
  data = ASCIItableReader(file);
  int pts = data.getnrow();
  // Terminate if number of columns is wrong; i.e. the wrong solar model file format.
  int n_cols = data.getncol();
  num_tracked_isotopes = n_cols-6;
  if (n_cols == 35) {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_C13", "X_N14", "X_N15", "X_O16", "X_O17", "X_O18", "X_Ne", "X_Na", "X_Mg", "X_Al", "X_Si", "X_P", "X_S", "X_Cl", "X_Ar",
                     "X_K", "X_Ca", "X_Sc", "X_Ti", "X_V", "X_Cr", "X_Mn", "X_Fe", "X_Co", "X_Ni");
    tracked_isotopes = {{"H",1}, {"He",4}, {"He",3}, {"C",12}, {"C",13}, {"N",14}, {"N",15}, {"O",16}, {"O",17}, {"O",18}, {"Ne",0}, {"Na",0}, {"Mg",0}, {"Al",0}, {"Si",0}, {"P",0}, {"S",0}, {"Cl",0}, {"Ar",0}, {"K",0}, {"Ca",0}, {"Sc",0},
                     {"Ti",0}, {"V",0}, {"Cr",0}, {"Mn",0}, {"Fe",0}, {"Co",0}, {"Ni",0}};
  } else if (n_cols == 12) {
    data.setcolnames("mass", "radius", "temperature", "rho", "pressure", "luminosity", "X_H1", "X_He4", "X_He3", "X_C12", "X_N14", "X_O16");
    tracked_isotopes = {{"H",1}, {"He",4}, {"He",3}, {"C",12}, {"N",14}, {"O",16}};
  } else {
    terminate_with_error("ERROR! Solar model file '"+file+"' not compatible with this code!");
  };

  // Initialise isotope-index map.
  for (int j = 0; j < num_tracked_isotopes; j++) { isotope_index_map[tracked_isotopes[j]] = j; };

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
    };
    r_lo = 0.0;
    pts += 1;
  };
  num_interp_pts = pts;

  const double* radius = &data["radius"][0];
  //   Calculate the necessary quantities and store them internally.
  for (int i = 0; i < pts; i++) {
    //  temperature
    temperature.push_back((1.0E-3*K2eV)*data["temperature"][i]);
    //  density
    density.push_back(data["rho"][i]);
    //  electron density
    // there are three ways of computing the electron density,
    //      1) raffelt approximation (assuming only fully ionised hydrogen and helium (Y = 1 - X)) (raffelt_approx = true)
    //      2) assuming full ionisation and summing over all elements (raffelt_approx = false)
    //      3) deduce from the pressure column of the model. first subtract radiation and ion oressure. (currently not implemented)
    double rhorel = data["rho"][i]/((1.0E+9*eV2g)*atomic_mass_unit);
    //  electron density from Raffelt
    double ne_Raff = 0.5 * (1.0 + data["X_H1"][i]) * data["rho"][i] /(atomic_mass_unit*eV2g*1.0E+9);
    n_e_Raff.push_back(ne_Raff);
    double ne = 0.0;
    //  electron density from summing over all elements (full ionisation)
    for (int j = 0; j < num_tracked_isotopes; j++) {
      // Isotope corresponds to column (index in 'tracked_isotopes' + 6) in data!
      Isotope isotope = tracked_isotopes[j];
      ne += data[j+6][i] * isotope.z_val() / atomic_weight(isotope);
    };
    n_e.push_back(ne*rhorel);
    //  ion density weighted by charge^2 from summing over all elements (full ionisation) and individual element densities
    double n_total_val = 0.0, z2_n_total_val = 0.0;
    for (int j = 0; j < num_tracked_isotopes; j++) {
      Isotope isotope = tracked_isotopes[j];
      double n = data[j+6][i] * rhorel / atomic_weight(isotope);
      n_total_val += n;
      n_isotope[j].push_back(n);
      double z2_n = gsl_pow_2(isotope.z_val()) * n;
      z2_n_total_val += z2_n;
      z2_n_isotope[j].push_back(z2_n);
    };

    n_total.push_back(n_total_val);
    z2_n_total.push_back(z2_n_total_val);

    // Calculate n_element for OP Code
    for (int k = 0; k < num_op_elements; k++) {
      double temp = 0.0;
      std::string element = op_element_names[k];
      for (int j = 0; j < num_tracked_isotopes; j++) {
        Isotope isotope = tracked_isotopes[j];
        if (isotope.name() == element) { temp += data[j+6][i] * rhorel / atomic_weight(isotope);
        };
      };
      n_op_element[k].push_back(temp);
    };
  };

  // Set up the interpolating functions quantities independent of iz.
  init_numbered_interp(0, radius, &temperature[0]); // Temperature
  init_numbered_interp(1, radius, &n_e[0]); // n_e from summing over elements (full ionisation)
  init_numbered_interp(2, radius, &n_e_Raff[0]); // n_e from Raffelt
  init_numbered_interp(3, radius, &density[0]); // Density
  init_numbered_interp(4, radius, &n_total[0]); // Number density of all ions (currently not used)
  init_numbered_interp(5, radius, &z2_n_total[0]); // Z^2 x number density (full ionisation)

  //  Quantities depending on specfific isotope or element
  n_isotope_acc.resize(num_tracked_isotopes);
  n_isotope_lin_interp.resize(num_tracked_isotopes);
  z2_n_isotope_acc.resize(num_tracked_isotopes);
  z2_n_isotope_lin_interp.resize(num_tracked_isotopes);
  for (int j = 0; j < num_tracked_isotopes; j++) {
    init_interp(n_isotope_acc[j], n_isotope_lin_interp[j], radius, &n_isotope[j][0]); // Ion density for each isotope
    init_interp(z2_n_isotope_acc[j], z2_n_isotope_lin_interp[j], radius, &z2_n_isotope[j][0]); // Ion density weighted by charge^2 for each isotope (full ionisation)
  };
  // N.B. The maps don't need to be resized or similar, since the access operator '[]' will silently create the map if it is missing.
  for (int k = 0; k < num_op_elements; k++) {
      std::string element = op_element_names[k];
      init_interp(n_element_acc[element], n_element_lin_interp[element], radius, &n_op_element[k][0]); // Ion density for each element (= summed isotopes with same Z)
  };
  //  OPACITY TABLES set up interpolating functions (only for chosen opacity code)
  //  OP opacities
  if (opcode == OP) {
    for (int k = 0; k < num_op_elements; k++) {
      std::string element = op_element_names[k];
      // Initialise grid values
      std::map<std::pair<int,int>, gsl_interp_accel*> temp_acc;
      std::map<std::pair<int,int>, gsl_spline*> temp_interp;
      for (int j = 0; j < op_grid_size; j++) {
        std::string op_filename = "data/opacity_tables/OP/opacity_table_"+element+"_"+std::to_string(op_grid[j][0])+"_"+std::to_string(op_grid[j][1])+".dat";
        ASCIItableReader op_data = ASCIItableReader(op_filename);

        // Determine the number of interpolated mass values.
        int op_pts = op_data[0].size();
        auto pr = std::make_pair(op_grid[j][0], op_grid[j][1]);
        temp_acc[pr] = gsl_interp_accel_alloc();
        temp_interp[pr] = gsl_spline_alloc (gsl_interp_linear, op_pts);
        const double* omega = &op_data[0][0];
        const double* s = &op_data[1][0];
        gsl_spline_init (temp_interp[pr], omega, s, op_pts);
      };
      opacity_acc_op[element] = temp_acc;
      opacity_lin_interp_op[element] = temp_interp;
    };
  };

  //  LEDCOP & ATOMIC (both TOPS) opacitites
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

  // OPAS
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
  //alpha
  if (alpha_available.find(solarmodel_name) != alpha_available.end()) {
     data_alpha =  ASCIItableReader("data/alpha"+file.substr(file.find("_")));
  }
  else {
     data_alpha =  ASCIItableReader("data/alpha_B16-AGSS09.dat");
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
}

// Move assignment operator
SolarModel& SolarModel::operator=(SolarModel &&model) {
  if (this != &model) {
    std::swap(data,model.data);
    std::swap(accel,model.accel);
    std::swap(linear_interp, model.linear_interp);
  };
  return *this;
}

// Class destructor
SolarModel::~SolarModel() {
  for (auto interp : linear_interp) { gsl_spline_free (interp); };
  for (auto interp : n_isotope_lin_interp) { gsl_spline_free (interp); };
  for (auto interp : z2_n_isotope_lin_interp) { gsl_spline_free (interp); };
  for (auto map : n_element_lin_interp) { gsl_spline_free (map.second); };
  for (auto map_1 : opacity_lin_interp_op) {
    for (auto map_2 : map_1.second) { gsl_spline_free (map_2.second); };
  };
  for (auto map : opacity_lin_interp_tops) { gsl_spline_free (map.second); };
  for (auto map : opacity_lin_interp_opas) { gsl_spline_free (map.second); };
  for (auto acc : accel) { gsl_interp_accel_free (acc); };
  for (auto acc : n_isotope_acc) { gsl_interp_accel_free (acc); };
  for (auto acc : z2_n_isotope_acc) { gsl_interp_accel_free (acc); };
  for (auto map : n_element_acc) { gsl_interp_accel_free (map.second); };
  for (auto map_1 : opacity_acc_op) {
    for (auto map_2 : map_1.second) { gsl_interp_accel_free (map_2.second); };
  };
  for (auto map : opacity_acc_tops) { gsl_interp_accel_free (map.second); };
  for (auto map : opacity_acc_opas) { gsl_interp_accel_free (map.second); };
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

// Routine to return the temperature (in keV) of the zone around the distance r from the centre of the Sun.
double SolarModel::temperature_in_keV(double r) { return gsl_spline_eval(linear_interp[0], r, accel[0]); }
// Same for density
double SolarModel::density(double r) { return gsl_spline_eval(linear_interp[3], r, accel[3]); }
// Routine to return the screening paramter kappa^2 in units of keV^2 (kappa^-1 = Debye-Hueckel radius).
double SolarModel::kappa_squared(double r) { return 4.0*pi*alpha_EM/temperature_in_keV(r)*(z2_n(r)+n_electron(r))*gsl_pow_3(keV2cm); }
// alpha is the expected contribution of all metals to z2_n per nucleon density: z2_n = (X + Y + \alpha Z)*\rho / m_u
double SolarModel::alpha(double r) {
    if (alpha_available.find(solarmodel_name) != alpha_available.end()) {
        return gsl_spline_eval(linear_interp[6], r, accel[6]);
    } else { return 4.0; }
}
double SolarModel::n_element(double r, std::string element) { return gsl_spline_eval(n_element_lin_interp.at(element), r, n_element_acc.at(element)); }
double SolarModel::H_mass_fraction(double r) { return n_element(r,"H")*atomic_weight({"H",1})*(1.0E+9*eV2g)*atomic_mass_unit/density(r); }
double SolarModel::He_mass_fraction(double r){ return n_element(r,"He")*atomic_weight({"He",4})*(1.0E+9*eV2g)*atomic_mass_unit/density(r); }
// Routine to return the number density times charge^2 of ion iz in the zone around the distance r from the centre of the Sun (full ionisation)
double SolarModel::z2_n_iz(double r, int isotope_index) { return gsl_spline_eval(z2_n_isotope_lin_interp[isotope_index], r, z2_n_isotope_acc[isotope_index]); }
// N.B. Convenience function below (may be slow for many calls!)
double SolarModel::z2_n_iz(double r, Isotope isotope) { int isotope_index = lookup_isotope_index(isotope); return z2_n_iz(r, isotope_index); }
// returns total z2_n without assuming full ionisation for some solar models where alpha is available
double SolarModel::z2_n(double r) {
    if (alpha_available.find(solarmodel_name) != alpha_available.end()) {
        return (H_mass_fraction(r) + He_mass_fraction(r)+ alpha(r) * metallicity(r)) * density(r)/((1.0E+9*eV2g)*atomic_mass_unit);
    } else {
        return  gsl_spline_eval(linear_interp[5], r, accel[5]);  // full ionisation
    }
}

// Routine to return the number density of ion iz in the zone around the distance r from the centre of the Sun.
double SolarModel::n_iz(double r, int isotope_index) { return gsl_spline_eval(n_isotope_lin_interp[isotope_index], r, n_isotope_acc[isotope_index]); }
// N.B. Convenience function below (may be slow for many calls!)
double SolarModel::n_iz(double r, Isotope isotope) { int isotope_index = lookup_isotope_index(isotope); return n_iz(r, isotope_index); }

// Routine to return the electron density in the zone around the distance r
double SolarModel::n_electron(double r) {
    if (raffelt_approx == false) { return gsl_spline_eval(linear_interp[1], r, accel[1]);}
    else { return gsl_spline_eval(linear_interp[2], r, accel[2]);}
}
double SolarModel::metallicity(double r){ return 1.0 - H_mass_fraction(r) - He_mass_fraction(r);}
// Routine to return the plasma freqeuency squared (in keV^2) of the zone around the distance r from the centre of the Sun.
double SolarModel::omega_pl_squared(double r) { return 4.0*pi*alpha_EM/m_electron*n_electron(r)*gsl_pow_3(keV2cm); }

//  auxiliary function required for FF and ee contribution
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

// Calculate the free-free contribution; from Eq. (2.17) of [arXiv:1310.0823] (assuming full ionisation)  for one isotope
double SolarModel::Gamma_P_ff(double omega, double r, int isotope_index) {
  //std::cout << "DEBUG INFO. FF at (" << omega << ", " << r << ", " << isotope_index << ") calculation..." << std::endl;
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
  double result = 0.0;
  const double prefactor1 = (8.0*sqrt(pi)/(3.0*sqrt(2.0))) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);

  if (omega == 0) { return 0; }

  if (raffelt_approx == false) { //(assuming full ionisation)  and only take H and He
    static int iso_ind_1 = lookup_isotope_index({"H",1}), iso_ind_2 = lookup_isotope_index({"He",3}), iso_ind_3 = lookup_isotope_index({"He",4});
    result = Gamma_P_ff(omega, r, iso_ind_1) + Gamma_P_ff(omega, r, iso_ind_2) + Gamma_P_ff(omega, r, iso_ind_3);
  } else { //contributions from all elements, no full ionisation assumed when alpha available
    double u = omega/temperature_in_keV(r);
    double y_red = sqrt(kappa_squared(r)/(2.0*m_electron*temperature_in_keV(r)));
    result = prefactor1 * n_electron(r)*z2_n(r)*exp(-u)*aux_function(u,y_red) / (omega*sqrt(temperature_in_keV(r))*pow(m_electron,3.5));
  };
  return result;
}

// Calculate the e-e bremsstrahlung contribution; from Eq. (2.18) of [arXiv:1310.0823]
double SolarModel::Gamma_P_ee(double omega, double r) {
  // N.B. "y" and "prefactor2" are different from the "y_red" and "prefactor1" above.
  //std::cout << "DEBUG INFO. ee at (" << omega << ", " << r << ") calculation..." << std::endl;
  if (omega == 0) {return 0;}
  const double prefactor2 = (4.0*sqrt(pi)/3.0) * gsl_pow_2(alpha_EM*g_aee) * gsl_pow_6(keV2cm);
  double u = omega/temperature_in_keV(r);
  double y = sqrt(kappa_squared(r)/(m_electron*temperature_in_keV(r)));
  return prefactor2 * gsl_pow_2(n_electron(r)) * exp(-u) * aux_function(u,y) / (omega * sqrt(temperature_in_keV(r)) * pow(m_electron,3.5));
}

// Calculate the Compton contribution; from Eq. (2.19) of [arXiv:1310.0823]
double SolarModel::Gamma_P_Compton(double omega, double r) {
  //std::cout << "DEBUG INFO. Compton at (" << omega << ", " << r << ") calculation..." << std::endl;
  if (omega == 0) {return 0;}
  const double prefactor3 = (alpha_EM/3.0) * pow(g_aee/(m_electron),2) * pow(keV2cm,3);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor3 * v*v*n_electron(r)/gsl_expm1(u);
}

// Read off interpolated elements for op, tops and opas
double SolarModel::op_grid_interp_erg(double u, int ite, int jne, std::string element) {
  double result = 0.0;
  auto grid_position = std::make_pair(ite,jne);

  if (unavailable_OP.find(grid_position) == unavailable_OP.end()) {
    if (opacity_lin_interp_op.find(element) == opacity_lin_interp_op.end()) {
      terminate_with_error("ERROR! OP data for element "+element+" does not exist.");
    } else if (opacity_lin_interp_op.at(element).find(grid_position) == opacity_lin_interp_op.at(element).end()) {
      std::cout << "WARNING. OP data for " << element << " at position ite = " << ite << " and jne = " << jne << " does not exist."  << std::endl;
    } else {
      gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);   // error handler modified to avoid boundary error (fill value = 0)
      result = gsl_spline_eval(opacity_lin_interp_op.at(element).at(grid_position), log(u), opacity_acc_op.at(element).at(grid_position));
      gsl_set_error_handler (old_handler);    //  reset the error handler
      if (gsl_isnan(result) == true) { return 0; };
    };
  };

  return result;
};

double SolarModel::tops_grid_interp_erg(double erg, float t, float rho) {
  auto key = std::make_pair(t,rho);
  if (opacity_lin_interp_tops.find(key) == opacity_lin_interp_tops.end()) {
      std::cout << "Grid point {" << t << ", " << rho << "} not found" << std::endl;
      return 0;
  }
  gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);   // error handler modified to avoid boundary error (fill value = 0)
  double result = gsl_spline_eval(opacity_lin_interp_tops.at(key), erg, opacity_acc_tops.at(key));
  gsl_set_error_handler (old_handler);    //  reset the error handler
  if (gsl_isnan(result) == true) { return 0; };
  return result;
};

double SolarModel::opas_grid_interp_erg(double erg, double r) {
    if (opacity_lin_interp_opas.find(r) == opacity_lin_interp_opas.end()) {
        std::cout << "OPAS data for R = " << r << " not found!"  << std::endl;
        return 0;
    }
    gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);   // error handler modified to avoid boundary error (fill value = 0)
    double result = gsl_spline_eval(opacity_lin_interp_opas.at(r), erg, opacity_acc_opas.at(r));
    gsl_set_error_handler (old_handler);    //  reset the error handler
    if (gsl_isnan(result) == true) { return 0; }
    return result;
};

//  double linear interpolation on solar grid (currently not used)
double SolarModel::opacity_table_interpolator_op2(double omega, double r, std::string element) {
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
  double result = (1.0-t1)*(1.0-t2)*op_grid_interp_erg(u1,ite1,jne1,element) + (1.0-t1)*t2*op_grid_interp_erg(u1,ite1,jne2,element)
                + t1*(1.0-t2)*op_grid_interp_erg(u2,ite2,jne1,element) + t1*t2*op_grid_interp_erg(u2,ite2,jne2,element);

  if (result < 0) {
    std::cout << "ERROR! Negative opacity!" << std::endl;
    std::cout << "Test 1: " << u1 << " " << u2 << " | " << ite1 << " " << ite2 << " | " << jne1 << " " << jne2 << std::endl;
    std::cout << "Test 2: " << t1 << " " << t2 << " | " << op_grid_interp_erg(u1,ite1,jne1,element) << " " << op_grid_interp_erg(u1,ite2,jne2,element) << " | " << result << std::endl;
  };
  return result;
};

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
    std::cout << "ERROR! Negative opacity!" << std::endl;
  };
  return result;
};

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
  double result = pow(pow(tops_grid_interp_erg(omega,Tlow,rholow),1.0-t2)*pow(tops_grid_interp_erg(omega,Tlow,rhoup),t2),1.0-t1)*   pow(pow(tops_grid_interp_erg(omega,Tup,rholow),1.0-t2)*pow(tops_grid_interp_erg(omega,Tup,rhoup),t2),t1);
  if (result < 0) {
    std::cout << "ERROR! Negative opacity!" << std::endl;
  };
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
    if (result < 0) {
      std::cout << "ERROR! Negative opacity!" << std::endl;
    };
    return result;
}

//  opacity for individual isotope (only possible for OP)
double SolarModel::opacity_element(double omega, double r, std::string element) {
    gsl_error_handler_t* old_handler = gsl_set_error_handler (&my_handler);    //other ahndler to avoid boundary errors
    if (opcode != OP) {
        std::cout << "WARNING. Chosen opacity code does not provide opacities for indivdual elements" << std::endl;
        return 0;
    }
    const double prefactor4 = a_Bohr*a_Bohr*(keV2cm);
    double u = omega/temperature_in_keV(r);
    double result = prefactor4*n_element(r, element)*opacity_table_interpolator_op(omega, r, element)*(-gsl_expm1(-u));
    gsl_set_error_handler (old_handler);
    return result;
};

// N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
// TODO However, maybe we care about opacity from isotope, i.e. not time n_element but times n_isotope...
double SolarModel::opacity_element(double omega, double r, Isotope isotope) { return opacity_element(omega, r, isotope.name()); }

//  opacity for total solar mixture
double SolarModel::opacity(double omega, double r) {
    double result = 0.0;
    if (opcode == OP) {
        for (int k = 0; k < num_op_elements; k++) { result += opacity_element(omega,r,op_element_names[k]); };
    }
    if ((opcode == LEDCOP) || (opcode == ATOMIC)){
        result = opacity_table_interpolator_tops(omega, r)*density(r)*keV2cm;
    }
    if (opcode == OPAS) {
        result = opacity_table_interpolator_opas(omega, r)*density(r)*keV2cm;
    }
    return result;
}
// opacity contribution from one isotope; first term of Eq. (2.21) of [arXiv:1310.0823]
double SolarModel::Gamma_P_opacity(double omega, double r, std::string element) {
  //std::cout << "DEBUG INFO. Opacity at (" << omega << ", " << r << ", " << element << ") calculation..." << std::endl;
  const double prefactor5 = 0.5*g_aee*g_aee/(4.0*pi*alpha_EM);
  double u = omega/temperature_in_keV(r);
  double v = omega/m_electron;
  return prefactor5*v*v*opacity_element(omega,r,element)/gsl_expm1(u);
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
    for (int k = 2; k < num_op_elements; k++) { element_contrib += Gamma_P_opacity(erg, r, op_element_names[k]); };
    result = element_contrib + Gamma_P_Compton(erg, r) + Gamma_P_ee(erg, r);
  } else if ((opcode == LEDCOP) || (opcode == ATOMIC)) {
    double u = erg/temperature_in_keV(r);
    double reducedCompton = 0.5*(1.0 - 1.0/gsl_expm1(u)) * Gamma_P_Compton(erg, r);
    result = Gamma_P_opacity(erg, r) + reducedCompton + Gamma_P_ee(erg, r);
  } else if (opcode == OPAS) {
    result = Gamma_P_opacity (erg, r);
  } else {
      terminate_with_error("ERROR! Unkown option for 'opcode' attribute. Use OP, LEDCOP, ATOMIC, or OPAS.");
  };
    return result;
}

double SolarModel::Gamma_P_Primakoff(double erg, double r) {
  if (erg == 0) {return 0;}
  const double prefactor6 = g_agg*g_agg/(32.0*pi);
  double ks_sq = kappa_squared(r);
  double w_pl_sq = omega_pl_squared(r);
  double T_in_keV = temperature_in_keV(r);
  double x = 4.0*(erg*erg)/ks_sq;
  if (w_pl_sq/(erg*erg) > 1.0) {return 0;}
  double phase_factor = 2.0*sqrt(1.0 - w_pl_sq/(erg*erg))/gsl_expm1(erg/T_in_keV);
  double rate = (ks_sq*T_in_keV)*((1.0 + 1.0/x)*gsl_log1p(x) - 1.0);
  return  prefactor6*phase_factor*rate;
}

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

std::vector<double> get_relevant_peaks(double erg_lo, double erg_hi) {
  const std::vector<double> all_peaks = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163,
                                         5.64394, 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};
  std::vector<double> result;
  result.push_back(erg_lo);
  for (auto peak_erg = all_peaks.begin(); peak_erg != all_peaks.end(); peak_erg++) { if ( (erg_lo < *peak_erg) && (*peak_erg < erg_hi) ) { result.push_back(*peak_erg); }; };
  result.push_back(erg_hi);

  return result;
}
