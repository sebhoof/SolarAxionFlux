// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __utils_hpp__
#define __utils_hpp__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>

#include <sys/stat.h> // Needed to check if file exists before we can expect C++14 std

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"

/////////////////////////////////////////////////////
//  General functions  (I/O, error handling, ...)  //
/////////////////////////////////////////////////////
class XFileNotFound : public std::runtime_error::runtime_error {
  public:
    XFileNotFound(): std::runtime_error::runtime_error("ERROR! File not found.") {}
    XFileNotFound(std::string file_name): std::runtime_error::runtime_error("ERROR! The file '"+file_name+"' was not found.") {}
};

class XSanityCheck : public std::runtime_error::runtime_error {
  public:
    XSanityCheck(): std::runtime_error::runtime_error("ERROR! The code failed a sanity check at runtime.") {}
    XSanityCheck(std::string err_msg): std::runtime_error::runtime_error("ERROR! "+err_msg) {}
};

class XUnsupportedOption : public std::runtime_error::runtime_error {
  public:
    XUnsupportedOption(): std::runtime_error::runtime_error("ERROR! Illegal option passed to the code.") {}
    XUnsupportedOption(std::string err_msg): std::runtime_error::runtime_error("ERROR! "+err_msg) {}
};

void my_global_exception_handler();
// Globally enforce the custom termination behaviour.
const auto terminator { std::set_terminate(my_global_exception_handler) };
void my_gsl_handler(const char * reason, const char * file, int line, int gsl_errno);

void terminate_with_error(std::string err_string);
void terminate_with_error_if(bool condition, std::string err_string);
bool file_exists(const std::string& filename);
void locate_data_folder(std::string path_to_model_file, std::string &path_to_data, std::string &model_file_name);
void save_to_file(std::string path, std::vector<std::vector<double>> buffer, std::string comment = "", bool overwrite = true);

// ASCIItableReader class from GAMBIT (main author: Christoph Weniger)
class ASCIItableReader {
  public:
    ASCIItableReader(std::string filename) { read(filename); }
    ASCIItableReader() {}
    ~ASCIItableReader() {}

    int read(std::string filename);
    void setcolnames(std::vector<std::string> names);

    template <typename... Args>
    void setcolnames(std::string name, Args... args) {
      std::vector<std::string> vec;
      vec.push_back(name);
      setcolnames(vec, args...);
    }
    template <typename... Args>
    void setcolnames(std::vector<std::string> vec, std::string name, Args... args) {
      vec.push_back(name);
      setcolnames(vec, args...);
    }

    void prepend_data(double datum, int i) { data[i].insert(data[i].begin(), datum); }
    const std::vector<double> & operator[] (int i) { return data[i]; }
    const std::vector<double> & operator[] (std::string name) { return data[colnames[name]]; }
    std::vector<std::vector<double> > get_data() { return data; }
    int getncol() { return data.size(); }
    int getnrow() { return data[0].size(); }

  private:
    std::vector<std::vector<double>> data;
    std::map<std::string, int> colnames;
};

///////////////////////////////
//  Interpolation functions  //
///////////////////////////////

// OneDInterpolator class: Provides a general one-dimensional interpolation container based on the gsl library.
// Can be declared static for efficiency & easy one-time initialisation of interpolating functions.
class OneDInterpolator {
  public:
    // Overloaded class creators for the OneDInterpolator class using the init function below.
    OneDInterpolator(std::string file, std::string type = "linear");
    OneDInterpolator(std::vector<std::vector<double> > table, std::string type = "linear");
    OneDInterpolator(const std::vector<double> &x, const std::vector<double> &y, std::string type = "linear");
    OneDInterpolator();
    OneDInterpolator (OneDInterpolator&&) = default;
    OneDInterpolator& operator=(OneDInterpolator&&);
    // Destructor
    ~OneDInterpolator();
    // Delete copy constructor and assignment operator to avoid shallow copies
    OneDInterpolator(const OneDInterpolator&) = delete;
    OneDInterpolator operator=(const OneDInterpolator&) = delete;
    // Routines to access interpolated values.
    double interpolate(double x);
    std::vector<double> interpolate(std::vector<double> x);
    // Routine to access upper and lower boundaries of available data.
    double lower();
    double upper();
    std::vector<std::vector<double> > get_data() { return data; }
  private:
    std::vector<std::vector<double> > data;
    // Initialiser for the OneDInterpolator class.
    void init(const std::vector<double> &x, const std::vector<double> &y, std::string type);
    void init(std::string type);
    // The gsl objects for the interpolating functions.
    gsl_interp_accel *acc;
    gsl_spline *spline;
    // Upper and lower boundaries available for the interpolating function.
    double lo, up;
};

// TwoDInterpolator class: Provides a general one-dimensional interpolation container based on the gsl library.
class TwoDInterpolator {
  public:
    // Overloaded class creators for the AxionInterpolator class using the init function below.
    TwoDInterpolator();
    TwoDInterpolator(std::string file, std::string type="bilinear");
    TwoDInterpolator(std::vector<std::vector<double> > table, std::string type="bilinear");
    TwoDInterpolator& operator=(TwoDInterpolator&&);
    // Destructor.
    ~TwoDInterpolator();
    // Delete copy constructor and assignment operator to avoid shallow copies.
    TwoDInterpolator(const TwoDInterpolator&) = delete;
    TwoDInterpolator operator=(const TwoDInterpolator&) = delete;
    // Routine to access interpolated values.
    double interpolate(double x, double y);
    // Routine to check if a point is inside the interpolating box.
    bool is_inside_box(double x, double y);
    std::vector<double> get_unique_x_vals() { return unique_x_vals; }
    std::vector<double> get_unique_y_vals() { return unique_y_vals; }
    std::vector<std::vector<double> > get_data() { return data; }
  private:
    // Initialiser for the TwoDInterpolator class.
    void init(std::string type);
    // The gsl objects for the interpolating functions that need to be available to the class routines.
    gsl_interp_accel *x_acc;
    gsl_interp_accel *y_acc;
    gsl_spline2d *spline;
    std::vector<std::vector<double> > data;
    std::vector<double> unique_x_vals;
    std::vector<double> unique_y_vals;
    // Upper and lower "x" and "y" values available to the interpolating function.
    double x_lo, y_lo, x_up, y_up;
};

////////////////////////////////////////////////////
//  General physics helper functions and classes  //
////////////////////////////////////////////////////

// Isotope class
class Isotope {
   public:
      Isotope();
      Isotope(std::string s, int a);
      Isotope(std::pair<std::string,int> p);
      Isotope(std::string s);
      ~Isotope() {}

      bool operator< (const Isotope& other) const;
      bool operator== (const Isotope& other) const;

      void init(std::string s, int a);
      std::string get_element_name() const;
      std::string index_name() const;
      int a_val() const;
      int z_val() const;
      bool same_z(Isotope *a);
   private:
      std::string element_name;
      int element_z_value;
      int isotope_a_value;
};

// Average weights of isotopes (depending of availability: exact, in the Sun, on Earth)
const std::map<Isotope, double> isotope_avg_weight { {{"H",1}, 1.007825}, {{"H",0}, 1.008}, {{"He",4}, 4.002603}, {{"He",3}, 3.016029},
    {{"He",0}, 4.002602}, {{"C",12}, 12.000000}, {{"C",13}, 13.003355}, {{"C",0}, 12.011}, {{"N",14}, 14.003074}, {{"N",15}, 15.000109}, {{"N",0}, 14.0067}, {{"O",16}, 15.994915},
    {{"O",17}, 16.999132}, {{"O",18}, 17.999160}, {{"O",0}, 15.999}, {{"Ne",0}, 20.1312812}, {{"Na",0}, 22.989769}, {{"Mg",0}, 24.3055}, {{"Al",0}, 26.9815385}, {{"Si",0}, 28.085}, {{"P",0}, 30.973762}, {{"S",0}, 32.0675},
                {{"Cl",0}, 35.4515}, {{"Ar",0}, 36.275403}, {{"K",0}, 39.0983}, {{"Ca",0}, 40.078}, {{"Sc",0}, 44.955908}, {{"Ti",0}, 47.867}, {{"V",0}, 50.9415}, {{"Cr",0}, 51.9961}, {{"Mn",0}, 54.938044},
                                                      {{"Fe",0}, 55.845}, {{"Co",0}, 58.933194}, {{"Ni",0}, 58.6934} };
double atomic_weight(Isotope isotope);

// Map opacity code names to labels; catalogue of available opacity data
enum opacitycode { OP, OPAS, LEDCOP, ATOMIC };
const std::map<opacitycode,std::string> opacitycode_name = { {OP,"OP"}, {OPAS,"OPAS"}, {LEDCOP,"LEDCOP"}, {ATOMIC,"ATOMIC"} };
const std::set<std::string> alpha_available = { "data/solar_models/SolarModel_AGS05.dat", "data/solar_models/SolarModel_AGSS09.dat", "data/solar_models/SolarModel_AGSS09ph.dat", "data/solar_models/SolarModel_B16-AGSS09.dat",
                                                "data/solar_models/SolarModel_B16-GS98.dat","data/solar_models/SolarModel_GS98.dat" };
const int op_grid_size = 213;
const int op_grid [op_grid_size][2] = { {150,54}, {150,56}, {152,54}, {152,56}, {156,68}, {158,68}, {158,70}, {160,64}, {160,66}, {160,68}, {160,70}, {162,64}, {162,66}, {162,68}, {162,70}, {164,66}, {164,68}, {164,70}, {164,72}, {166,66}, {166,68},
{166,70}, {166,72}, {166,74}, {168,68}, {168,70}, {168,72}, {168,74}, {168,76}, {170,70}, {170,72}, {170,74}, {172,72}, {172,74}, {170,76}, {172,76}, {174,74}, {174,76}, {174,78}, {176,74}, {176,76}, {176,78}, {176,80}, {178,76}, {178,78}, {178,80},
{180,76}, {180,78}, {180,80}, {182,78}, {182,80}, {184,78}, {184,80}, {184,82}, {186,78}, {186,80}, {186,82}, {188,80}, {188,82}, {190,80}, {190,82}, {190,84}, {192,80}, {192,82}, {192,84}, {194,80}, {194,82}, {194,84}, {196,80}, {196,82}, {196,84},
{198,82}, {198,84}, {198,86}, {200,82}, {200,84}, {200,86}, {202,82}, {202,84}, {202,86}, {204,82}, {204,84}, {204,86}, {206,82}, {206,84}, {206,86}, {208,84}, {208,86}, {208,88}, {210,84}, {210,86}, {210,88}, {212,84}, {212,86}, {212,88}, {214,84},
{214,86}, {214,88}, {216,84}, {216,86}, {216,88}, {218,86}, {218,88}, {220,86}, {220,88}, {220,90}, {222,86}, {222,88}, {222,90}, {224,86}, {224,88}, {224,90}, {226,86}, {226,88}, {226,90}, {228,86}, {228,88}, {228,90}, {230,86}, {230,88}, {230,90},
{232,88}, {232,90}, {232,92}, {234,88}, {234,90}, {234,92}, {236,88}, {236,90}, {236,92}, {238,88}, {238,90}, {238,92}, {240,88}, {240,90}, {240,92}, {242,88}, {242,90}, {242,92}, {244,90}, {244,92}, {246,90}, {246,92}, {246,94}, {248,90}, {248,92},
{248,94}, {250,90}, {250,92}, {250,94}, {252,90}, {252,92}, {252,94}, {254,90}, {254,92}, {254,94}, {256,90}, {256,92}, {256,94}, {256,96}, {258,92}, {258,94}, {258,96}, {260,92}, {260,94}, {260,96}, {260,98}, {262,92}, {262,94}, {262,96}, {262,98},
{264,94}, {264,96}, {264,98}, {266,94}, {266,96}, {266,98}, {266,100}, {268,96}, {268,98}, {268,100}, {270,96}, {270,98}, {270,100}, {270,102}, {272,96}, {272,98}, {272,100}, {272,102}, {274,98}, {274,100}, {274,102}, {276,98}, {276,100}, {276,102},
{278,100}, {278,102}, {278,104}, {280,100}, {280,102}, {280,104}, {282,100}, {282,102}, {282,104}, {284,100}, {284,102}, {284,104}, {286,102}, {286,104}, {286,106}, {288,102}, {288,104}, {288,106} };
const std::set<std::pair<int,int>> unavailable_OP = { {152,68}, {152,70}, {154,68}, {154,70}, {156,70} };
const int num_op_elements = 17;
const std::string op_element_names [num_op_elements] = { "H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "Ca", "Cr", "Mn", "Fe", "Ni" };
const std::vector<std::vector<int>> op_elements = { {0}, {1,2}, {3,4}, {5,6}, {7,8,9}, {10}, {11}, {12}, {13}, {14}, {16}, {18}, {20}, {24}, {25}, {26}, {28} };
const std::vector<std::vector<float> > ledcop_grid =  { {0.5,10.175}, {0.5,13.684}, {0.5,18.308}, {0.6,10.175}, {0.6,13.684}, {0.6,18.308}, {0.6,24.268}, {0.6,31.802}, {0.6,41.156}, {0.8,13.684}, {0.8,18.308}, {0.8,24.268}, {0.8,31.802},
                                                        {0.8,41.156}, {0.8,52.611}, {0.8,66.544}, {1.0,31.802}, {1.0,41.156}, {1.0,52.611}, {1.0,66.544}, {1.0,83.466}, {1.0,103.442}, {1.0,124.995}, {1.25,52.611}, {1.25,66.544}, {1.25,83.466},
                                                        {1.25,103.442}, {1.25,124.995}, {1.25,143.111}, {1.25,150.5}, {1.5,103.442}, {1.5,124.995}, {1.5,143.111}, {1.5,150.5} };
const std::vector<float> ledcop_temperatures = { 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5 };
const std::vector<float> ledcop_densities = { 10.175, 13.684, 18.308, 24.268, 31.802, 41.156, 52.611, 66.544, 83.466, 103.442, 124.995, 143.111, 150.5 };
const std::vector<std::vector<float> > atomic_grid = { {0.55,10.175}, {0.55,13.684},  {0.55,18.308}, {0.6,10.175}, {0.6,13.684},  {0.6,18.308}, {0.6,24.268}, {0.65,13.684}, {0.65,18.308}, {0.65,24.268}, {0.7,18.308}, {0.7,24.268}, {0.7,31.802},
                                                       {0.7,41.156}, {0.8,18.308}, {0.8,24.268}, {0.8,31.802}, {0.8,41.156}, {0.8,52.611}, {0.9,31.802}, {0.9,41.156}, {0.9,52.611}, {0.9,66.544}, {1.0,41.156}, {1.0,52.611}, {1.0,66.544},
                                                       {1.0,83.466}, {1.0,103.442}, {1.125,52.611}, {1.125,66.544}, {1.125,83.466}, {1.125,103.442}, {1.125,124.995}, { 1.25,83.466}, {1.25,103.442}, {1.25,124.995}, {1.25,143.111}, {1.25,150.5},
                                                       {1.375,103.442}, {1.375,124.995}, {1.375,143.111}, {1.375,150.5} };
const std::vector<float> atomic_temperatures = { 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.125, 1.25, 1.375 };
const std::vector<float> atomic_densities = { 10.175, 13.684, 18.308, 24.268, 31.802, 41.156, 52.611, 66.544, 83.466, 103.442, 124.995, 143.111, 150.5 };
const std::vector<double> opas_radii = { 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71 };

// Get the (approximate) locations of the peaks in the axion-electron spectrum to allow for more accurate integration
std::vector<double> get_relevant_peaks(double erg_lo, double erg_hi);


///////////////////////////////////
//  Other convenience functions  //
///////////////////////////////////

double safe_log10(double x, double lgx0);

#endif // defined __utils_hpp__
