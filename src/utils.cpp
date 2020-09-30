// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

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

void locate_data_folder(std::string path_to_model_file, std::string &path_to_data, std::string &model_file_name) {
  auto pos = path_to_model_file.find_last_of("/");
  if (pos!= std::string::npos) {
    path_to_data = path_to_model_file.substr(0, pos);
    model_file_name = path_to_model_file.substr(pos+1, path_to_model_file.length());
  } else {
    path_to_data = ".";
    model_file_name = path_to_model_file;
  };
  path_to_data += "/../"; // Since we expect model file to be in data/solar_models.
};

void save_to_file(std::string path, std::vector<std::vector<double>> buffer, std::string comment, bool overwrite) {
  if (path != "") {
    if (file_exists(path)) {
      if (overwrite) {
        //std::cout << "WARNING! File " << path << " exists and will be overwritten..." << std::endl;
      } else {
        std::cout << "File " << path << " exists! Now saving to " << path << "_new" << std::endl; path += "_new";
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
      output << "# " << comment << std::endl;
    };
    output << std::scientific << std::setprecision(8);
    int n_cols = buffer.size();
    if (n_cols > 0) {
      int n_rows = buffer[0].size();
      if (n_rows > 0) {
        for (int i=0; i<n_rows; ++i) {
          output << buffer[0][i];
          for (int j=1; j<n_cols; ++j) { output << " " << buffer[j][i]; };
          output << std::endl;
        };
      } else {
        std::cout << "WARNING! The data you are trying to write to " << path << " is empty! Created an empty file." << std::endl;
      };
    };
    output.close();
  };
}

OneDInterpolator::OneDInterpolator() {
  // NOTE. Allocate memory in the default constructor for the move-assign operator (via std::swap) works with the destructor.
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, 2);
}

void OneDInterpolator::init(const std::vector<double> &x, const std::vector<double> &y, std::string type) {
  // Initialise gsl interpolation routine.
  int pts = x.size();
  const double* x_ptr = &x[0];
  const double* y_ptr = &y[0];
  acc = gsl_interp_accel_alloc();
  if (type == "cspline") {
    //spline_type = gsl_interp_cspline;
    spline = gsl_spline_alloc (gsl_interp_cspline, pts);
  } else if (type == "linear") {
    //spline_type = gsl_interp_linear;
    spline = gsl_spline_alloc(gsl_interp_linear, pts);
  } else {
    terminate_with_error("ERROR! Interpolation type '"+type+"' not known to class OneDInterpolator.\n       Available types: 'linear' and 'cspline'.");
  };
  gsl_spline_init(spline, x_ptr, y_ptr, pts);
  lo = x.front();
  up = x.back();
}

void OneDInterpolator::init(std::string type) { init(data[0], data[1], type); }

// Initialiser for the OneDInterpolator class.
OneDInterpolator::OneDInterpolator(std::string file, std::string type) {
  // Check if file exists.
  terminate_with_error_if(not(file_exists(file)), "ERROR! File '"+file+"' for interpolation not found!");
  // Read numerical values from data file.
  ASCIItableReader tab (file);
  //tab.setcolnames("x", "y"); -> Would raise warning if we try to interpolate data with errors

  data = tab.get_data();
  init(type);
  //init(tab[0], tab[1], type);
}

OneDInterpolator::OneDInterpolator(const std::vector<double> &x, const std::vector<double> &y, std::string type) { init(x, y, type); }
OneDInterpolator::OneDInterpolator(std::vector<std::vector<double> > table, std::string type) { data = table; init(type); }

// Move constructor
//OneDInterpolator::OneDInterpolator(OneDInterpolator&& src) {
//    std::cout << "Entering custom move constructor..." << std::endl;
//    lo = std::move(src.lo);
//    up = std::move(src.up);
//    acc = std::move(src.acc);
//    spline = std::move(src.spline);
//}

// Move constructor
//OneDInterpolator::OneDInterpolator(OneDInterpolator&& src) {
//  std::swap(acc,src.acc);
//  std::swap(spline,src.spline);
//  std::swap(lo,src.lo);
//  std::swap(up,src.up);
//}

// Move assignment operator
OneDInterpolator& OneDInterpolator::operator=(OneDInterpolator&& src) {
  if(this != &src) {
    std::swap(acc,src.acc);
    std::swap(spline,src.spline);
    std::swap(data,src.data);
    std::swap(lo,src.lo);
    std::swap(up,src.up);
  }
  return *this;
}

// Destructor
OneDInterpolator::~OneDInterpolator() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

// Routine to access interpolated values.
double OneDInterpolator::interpolate(double x) { return gsl_spline_eval(spline, x, acc); }

std::vector<double> OneDInterpolator::interpolate(std::vector<double> x) {
  std::vector<double> result;
  for (auto it = x.begin(); it != x.end(); it++) { result.push_back(interpolate(*it)); };
  return result;
}

// Routines to return upper and lower boundaries of interpolating function.
double OneDInterpolator::lower() { return lo; }
double OneDInterpolator::upper() { return up; }


TwoDInterpolator::TwoDInterpolator() {
  // NOTE. Allocate memory in the default constructor for the move-assign operator (via std::swap) works with the destructor.
  x_acc = gsl_interp_accel_alloc();
  y_acc = gsl_interp_accel_alloc();
  spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, 2, 2);
}

// Move assignment operator
TwoDInterpolator& TwoDInterpolator::operator=(TwoDInterpolator&& interp)
{
  if(this != &interp)
  {
    std::swap(x_acc,interp.x_acc);
    std::swap(y_acc,interp.y_acc);
    std::swap(spline,interp.spline);
    std::swap(x_lo,interp.x_lo);
    std::swap(x_up,interp.x_up);
    std::swap(y_lo,interp.y_lo);
    std::swap(y_up,interp.y_up);
  }
  return *this;
}

TwoDInterpolator::TwoDInterpolator(std::string file, std::string type) {
  terminate_with_error_if(not(file_exists(file)), "ERROR! File '"+file+"' for interpolation not found!");
  ASCIItableReader table (file);
  data = table.get_data();
  init(type);
}

TwoDInterpolator::TwoDInterpolator(std::vector<std::vector<double>> table, std::string type) { data = table; init(type); }

// Initialiser for the AxionInterpolator class.
void TwoDInterpolator::init(std::string type) {
  std::cout << "Entering TwoDInterpolator init..." << std::endl;
  unique_x_vals = data[0];
  std::sort(unique_x_vals.begin(), unique_x_vals.end());
  unique_x_vals.erase(std::unique(unique_x_vals.begin(), unique_x_vals.end()), unique_x_vals.end());
  int nx = unique_x_vals.size();
  unique_y_vals = data[1];
  std::sort(unique_y_vals.begin(), unique_y_vals.end());
  unique_y_vals.erase(std::unique(unique_y_vals.begin(), unique_y_vals.end()), unique_y_vals.end());
  int ny = unique_y_vals.size();
  int n_grid_pts = data[2].size();

  terminate_with_error_if(nx*ny != n_grid_pts, "ERROR! The number of grid points ("+std::to_string(n_grid_pts)+") for TwoDInterpolator does not equal the number of unique 'x' and 'y' values ("+std::to_string(nx)+" and "+std::to_string(ny)+")!.");

  const double* x = &unique_x_vals[0];
  const double* y = &unique_y_vals[0];
  // Allocate memory for "z" values array in gsl format
  double* z = (double*) malloc(nx * ny * sizeof(double));

  if (type == "bicubic") {
    spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nx, ny);
  } else if (type == "bilinear") {
    spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
  } else {
    terminate_with_error("ERROR! Interpolation type '"+type+"' not known to class TwoDInterpolator.\n       Available types: 'bilinear' and 'bicubic'.");
  };

  x_acc = gsl_interp_accel_alloc();
  y_acc = gsl_interp_accel_alloc();

  // Determine first and last "x" and "y" values and grid step size.
  x_lo = unique_x_vals.front();
  x_up = unique_x_vals.back();
  y_lo = unique_y_vals.front();
  y_up = unique_y_vals.back();
  double x_delta = (x_up-x_lo) / (nx-1);
  double y_delta = (y_up-y_lo) / (ny-1);

  // Intialise grid.
  for (int i = 0; i < n_grid_pts; i++)
  {
    // Determine appropriate indices for the grid points.
    double temp = (data[0][i]-x_lo) / x_delta;
    int ind_x = (int) (temp+0.5);
    temp = (data[1][i]-y_lo) / y_delta;
    int ind_y = (int) (temp+0.5);

    gsl_spline2d_set(spline, z, ind_x, ind_y, data[2][i]);
  };
    gsl_spline2d_init (spline, x, y, z, nx, ny);
    std::cout << "TwoDInterpolator init done!" << std::endl;
};

TwoDInterpolator::~TwoDInterpolator() {
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(x_acc);
  gsl_interp_accel_free(y_acc);
}

// Routine to access interpolated values.
double TwoDInterpolator::interpolate(double x, double y) { return gsl_spline2d_eval(spline, x, y, x_acc, y_acc); };

// Routine to check if a point is inside the interpolating box.
bool TwoDInterpolator::is_inside_box(double x, double y) { return ((x >= x_lo) && (x <= x_up) && (y >= y_lo) && (y <= y_up)); };


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

std::vector<double> get_relevant_peaks(double erg_lo, double erg_hi) {
  const std::vector<double> all_peaks = {0.653029, 0.779074, 0.920547, 0.956836, 1.02042, 1.05343, 1.3497, 1.40807, 1.46949, 1.59487, 1.62314, 1.65075, 1.72461, 1.76286, 1.86037, 2.00007, 2.45281, 2.61233, 3.12669, 3.30616, 3.88237, 4.08163,
                                         5.64394, 5.76064, 6.14217, 6.19863, 6.58874, 6.63942, 6.66482, 7.68441, 7.74104, 7.76785};
  std::vector<double> result;
  result.push_back(erg_lo);
  for (auto peak_erg = all_peaks.begin(); peak_erg != all_peaks.end(); peak_erg++) { if ( (erg_lo < *peak_erg) && (*peak_erg < erg_hi) ) { result.push_back(*peak_erg); }; };
  result.push_back(erg_hi);

  return result;
}
