#ifndef __utils_hpp__
#define __utils_hpp__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <sstream>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"

void terminate_with_error(std::string err_string);
void my_handler (const char * reason, const char * file, int line, int gsl_errno);
const double abs_prec = 0, rel_prec = 1.0e-4;
const int method = 5;
struct integrand_params {double u; double y;};

// Note: For future versions, use a list of 'Isotope's that are available in a SolarModel
// el_a_value=0 is a default for 'only tracked in bulk';
// Columns is Solar model can stay as they are.
// -> Variable iz becomes type 'Isotope'
// For the OP code, just use the op list with names(!) as below and automatically sum over all Isotopes with iz[0] == IsotopeName!
// If a_val() > 0, use this, else get default bulk value from global table (maybe create own namespace).
class Isotope {
   public:
      Isotope() {};
      Isotope(std::string s, int a) { element_name = s; isotope_a_value = a; };
      Isotope(std::pair<std::string,int> p) { element_name = p.first; isotope_a_value = p.second; };
      // This is for convenience in order to define elements as an Isotope
      // TODO: if el_a_value < -1, trigger adding up all values for el_name
      Isotope(std::string s) { element_name = s; isotope_a_value = -1; };

      bool operator< (const Isotope& other) const;
      bool operator== (const Isotope& other) const;

      std::string name() const;
      std::string index_name() const;
      int a_val() const;
      int z_val() const;
      bool same_z(Isotope *a);
   private:
      std::string element_name;
      int isotope_a_value;
};

// Map of isotope name -> Z value
const std::map<std::string, int> element_z_value ({ {"H", 1}, {"He", 2}, {"C", 6}, {"N", 7}, {"O", 8}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17},
                                               {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28} });

const std::map<Isotope, double> isotope_avg_weight ({ {{"H",1}, 1.007825}, {{"He",4}, 4.002603}, {{"He",3}, 3.016029}, {{"C",12}, 12.000000}, {{"C",13}, 13.003355}, {{"N",14}, 14.003074}, {{"N",15}, 15.000109}, {{"O",16}, 15.994915},
                                                      {{"O",17}, 16.999132}, {{"O",18}, 17.999160}, {{"Ne",0}, 20.1312812}, {{"Na",0}, 22.989769}, {{"Mg",0}, 24.3055}, {{"Al",0}, 26.9815385}, {{"Si",0}, 28.085}, {{"P",0}, 30.973762}, {{"S",0}, 32.0675},
                                                      {{"Cl",0}, 35.4515}, {{"Ar",0}, 36.275403}, {{"K",0}, 39.0983}, {{"Ca",0}, 40.078}, {{"Sc",0}, 44.955908}, {{"Ti",0}, 47.867}, {{"V",0}, 50.9415}, {{"Cr",0}, 51.9961}, {{"Mn",0}, 54.938044},
                                                      {{"Fe",0}, 55.845}, {{"Co",0}, 58.933194}, {{"Ni",0}, 58.6934} });

double atomic_weight(Isotope isotope);

enum opacitycode {OP, OPAS, LEDCOP, ATOMIC};
const std::string opacitycode_names [4] = {"OP", "OPAS", "LEDCOP", "ATOMIC"};

const int op_grid_size = 212;
const int op_grid [op_grid_size][2] = {{150,54}, {150,56}, {152,54}, {152,56}, {156,68}, {158,68}, {158,70}, {160,64}, {160,66}, {160,68}, {160,70}, {162,64}, {162,66}, {162,68}, {162,70}, {164,66}, {164,68}, {164,70}, {164,72}, {166,66}, {166,68}, {166,70}, {166,72}, {166,74}, {168,68}, {168,70}, {168,72}, {168,74}, {170,70}, {170,72}, {170,74}, {172,72}, {172,74}, {170,76}, {172,76}, {174,74}, {174,76}, {174,78}, {176,74}, {176,76}, {176,78}, {176,80}, {178,76}, {178,78}, {178,80}, {180,76}, {180,78}, {180,80}, {182,78}, {182,80}, {184,78}, {184,80}, {184,82}, {186,78}, {186,80}, {186,82}, {188,80}, {188,82}, {190,80}, {190,82}, {190,84}, {192,80}, {192,82}, {192,84}, {194,80}, {194,82}, {194,84}, {196,80}, {196,82}, {196,84}, {198,82}, {198,84}, {198,86}, {200,82}, {200,84}, {200,86}, {202,82}, {202,84}, {202,86}, {204,82}, {204,84}, {204,86}, {206,82}, {206,84}, {206,86}, {208,84}, {208,86}, {208,88}, {210,84}, {210,86}, {210,88}, {212,84}, {212,86}, {212,88}, {214,84}, {214,86}, {214,88}, {216,84}, {216,86}, {216,88}, {218,86}, {218,88}, {220,86}, {220,88}, {220,90}, {222,86}, {222,88}, {222,90}, {224,86}, {224,88}, {224,90}, {226,86}, {226,88}, {226,90}, {228,86}, {228,88}, {228,90}, {230,86}, {230,88}, {230,90}, {232,88}, {232,90}, {232,92}, {234,88}, {234,90}, {234,92}, {236,88}, {236,90}, {236,92}, {238,88}, {238,90}, {238,92}, {240,88}, {240,90}, {240,92}, {242,88}, {242,90}, {242,92}, {244,90}, {244,92}, {246,90}, {246,92}, {246,94}, {248,90}, {248,92}, {248,94}, {250,90}, {250,92}, {250,94}, {252,90}, {252,92}, {252,94}, {254,90}, {254,92}, {254,94}, {256,90}, {256,92}, {256,94}, {256,96}, {258,92}, {258,94}, {258,96}, {260,92}, {260,94}, {260,96}, {260,98}, {262,92}, {262,94}, {262,96}, {262,98}, {264,94}, {264,96}, {264,98}, {266,94}, {266,96}, {266,98}, {266,100}, {268,96}, {268,98}, {268,100}, {270,96}, {270,98}, {270,100}, {270,102}, {272,96}, {272,98}, {272,100}, {272,102}, {274,98}, {274,100}, {274,102}, {276,98}, {276,100}, {276,102}, {278,100}, {278,102}, {278,104}, {280,100}, {280,102}, {280,104}, {282,100}, {282,102}, {282,104}, {284,100}, {284,102}, {284,104}, {286,102}, {286,104}, {286,106}, {288,102}, {288,104}, {288,106}};
const std::set<std::pair<int,int>> unavailable_OP = {{152,68}, {152,70}, {154,68}, {154,70}, {156,70}};
const int num_op_elements = 17;
const std::string op_element_names [num_op_elements] = {"H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "Ca", "Cr", "Mn", "Fe", "Ni"};
const std::vector<std::vector<int>> op_elements = {{0}, {1,2}, {3,4}, {5,6}, {7,8,9}, {10}, {11}, {12}, {13}, {14}, {16}, {18}, {20}, {24}, {25}, {26}, {28}};
const std::vector<std::vector<float>> ledcop_grid =  {{0.5,10.175}, {0.5,13.684}, {0.5,18.308}, {0.6,10.175}, {0.6,13.684}, {0.6,18.308}, {0.6,24.268}, {0.6,31.802}, {0.6,41.156}, {0.8,13.684}, {0.8,18.308}, {0.8,24.268}, {0.8,31.802}, {0.8,41.156}, {0.8,52.611}, {0.8,66.544}, {1.0,31.802}, {1.0,41.156}, {1.0,52.611}, {1.0,66.544}, {1.0,83.466}, {1.0,103.442}, {1.0,124.995}, {1.25,52.611}, {1.25,66.544}, {1.25,83.466}, {1.25,103.442}, {1.25,124.995}, {1.25,143.111}, {1.25,150.5}, {1.5,103.442}, {1.5,124.995}, {1.5,143.111}, {1.5,150.5}};
const std::vector<float> ledcop_temperatures = {0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5};
const std::vector<float> ledcop_densities = {10.175, 13.684, 18.308, 24.268, 31.802, 41.156, 52.611, 66.544, 83.466, 103.442, 124.995, 143.111, 150.5};
const std::vector<std::vector<float>> atomic_grid = {{0.55,10.175}, {0.55,13.684},  {0.55,18.308}, {0.6,10.175}, {0.6,13.684},  {0.6,18.308}, {0.6,24.268}, {0.65,13.684}, {0.65,18.308}, {0.65,24.268}, {0.7,18.308}, {0.7,24.268}, {0.7,31.802}, {0.7,41.156}, {0.8,18.308}, {0.8,24.268}, {0.8,31.802}, {0.8,41.156}, {0.8,52.611}, {0.9,31.802}, {0.9,41.156}, {0.9,52.611}, {0.9,66.544}, {1.0,41.156}, {1.0,52.611}, {1.0,66.544}, {1.0,83.466}, {1.0,103.442}, {1.125,52.611}, {1.125,66.544}, {1.125,83.466}, {1.125,103.442}, {1.125,124.995}, {1.25,83.466}, {1.25,103.442}, {1.25,124.995}, {1.25,143.111}, {1.25,150.5}, {1.375,103.442}, {1.375,124.995}, {1.375,143.111}, {1.375,150.5}};
const std::vector<float> atomic_temperatures = {0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.125, 1.25, 1.375};
const std::vector<float> atomic_densities = {10.175, 13.684, 18.308, 24.268, 31.802, 41.156, 52.611, 66.544, 83.466, 103.442, 124.995, 143.111, 150.5};
const std::vector<double> opas_radii = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71};

bool file_exists(const std::string& filename);
void save_to_file(std::string path, std::vector<std::vector<double>> data, bool overwrite);
void save_to_file(std::string path, std::vector<std::vector<double>> data);

// OneDInterpolator class: Provides a general 1-D interpolation container based on the gsl library.
// Can be declared static for efficiency & easy one-time initialisation of interpolating functions.
class OneDInterpolator
{
  public:
    // Overloaded class creators for the OneDInterpolator class using the init function below.
    OneDInterpolator(std::string file, std::string type);
    OneDInterpolator(std::string file);
    OneDInterpolator();
    OneDInterpolator& operator=(OneDInterpolator&&);
    // Destructor
    ~OneDInterpolator();
    // Delete copy constructor and assignment operator to avoid shallow copies
    OneDInterpolator(const OneDInterpolator&) = delete;
    OneDInterpolator operator=(const OneDInterpolator&) = delete;
    // Routine to access interpolated values.
    double interpolate(double x);
    // Routine to access upper and lower boundaries of available data.
    double lower();
    double upper();
  private:
    // Initialiser for the OneDInterpolator class.
    void init(std::string file, std::string type);
    // The gsl objects for the interpolating functions.
    gsl_interp_accel *acc;
    gsl_spline *spline;
    // Upper and lower boundaries available for the interpolating function.
    double lo;
    double up;
};

class ASCIItableReader {
  public:
    ASCIItableReader(std::string filename) { read(filename); };
    ASCIItableReader() {};  // Dummy initializer
    ~ASCIItableReader() {}

    int read(std::string filename);
    void setcolnames(std::vector<std::string> names);

    template <typename... Args>
    void setcolnames(std::string name, Args... args)
    {
      std::vector<std::string> vec;
      vec.push_back(name);
      setcolnames(vec, args...);
    }
    template <typename... Args>
    void setcolnames(std::vector<std::string> vec, std::string name, Args... args)
    {
      vec.push_back(name);
      setcolnames(vec, args...);
    }

    void prepend_data(double datum, int i) { data[i].insert(data[i].begin(), datum); };
    const std::vector<double> & operator[] (int i) { return data[i]; };
    const std::vector<double> & operator[] (std::string name) { return data[colnames[name]]; };
    int getncol() { return data.size(); }
    int getnrow() { return data[0].size(); }

  private:
    std::vector<std::vector<double> > data;
    std::map<std::string, int> colnames;
};


// SolarModel class: Provides a container to store a (tabulated) Solar model and functions to return its properties.
class SolarModel
{
  public:
    SolarModel();
    SolarModel(std::string file, opacitycode set_opcode, const bool set_raffelt_approx);
    SolarModel(std::string file, opacitycode set_opcode);
    ~SolarModel();
    SolarModel& operator=(SolarModel&&);
    // Delete copy constructor and assignment operator to avoid shallow copies
    SolarModel(const SolarModel&) = delete;
    SolarModel operator=(const SolarModel&) = delete;
    // Min. and max. radius of the solar model file (distance r from the centre of the Sun in units of the solar radius)
    bool raffelt_approx;  // whether to use the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf) egs. 16 a-c, default is true
    const opacitycode opcode;
    double r_lo, r_hi;
    int lookup_isotope_index(Isotope isotope);
    // Routine to return the screening parameter kappa^2 (kappa^-1 = Debye-Hueckel radius).
    double kappa_squared(double r);
    // Routine to return the temperature in keV.
    double temperature_in_keV(double r);
    // Routine to return density
    double density(double r);
    // Routines to return ion number density times charge^2 (for each isotope and for all)
    double z2_n_iz(double r, int isotope_index);
    double z2_n_iz(double r, Isotope isotope);
    double z2_n(double r);
    // Routine to return ion density for each isotope
    double n_iz(double r, int isotope_index);
    double n_iz(double r, Isotope isotope);
    // Routine
    double n_element(double r, std::string el_name);
    // Routine to return electron density
    double n_electron(double r);
    // Routine to return the plasma frequency squared.
    double omega_pl_squared(double r);
    // Rates for various axion prodcution processes
    double Gamma_P_ff(double omega, double r, int isotope_index);
    double Gamma_P_ff(double omega, double r, Isotope isotope);
    double Gamma_P_ff(double omega, double r);
    double Gamma_P_ee(double omega, double r);
    double Gamma_P_Compton(double omega, double r);
    double Gamma_P_opacity(double omega, double r, std::string element);
    double Gamma_P_opacity(double omega, double r);
    double Gamma_P_Primakoff(double omega, double r);
    double Gamma_P_all_electron(double omega, double r);
    // interpolating the opacity data
    double op_grid_interp_erg(double u, int ite, int jne, std::string element);
    double tops_grid_interp_erg(double erg, float t, float rho);
    double opas_grid_interp_erg(double erg, double r);
    double opacity_table_interpolator_op(double omega, double r, std::string element);
    double opacity_table_interpolator_op2(double omega, double r, std::string element);
    double opacity_table_interpolator_tops(double omega, double r);
    double opacity_table_interpolator_opas(double omega, double r);
    double opacity_element(double omega, double r, std::string element);
    // N.B. Opacity only depends on chemical properties; below just overloaded for convenience;
    double opacity_element(double omega, double r, Isotope isotope);
    double opacity(double omega, double r);
  private:
    // solar data
    ASCIItableReader data;
    std::map<Isotope,int> isotope_index_map;
    int num_tracked_isotopes;
    std::vector<Isotope> tracked_isotopes;
    gsl_interp_accel *accel[7];
    gsl_spline *linear_interp[7];
    //std::vector< std::map<std::pair<int,int>, gsl_interp_accel*> > opacity_acc_op;
    //std::vector< std::map<std::pair<int,int>, gsl_spline*> > opacity_lin_interp_op;
    std::map< std::string, std::map<std::pair<int,int>, gsl_interp_accel*> > opacity_acc_op;
    std::map< std::string, std::map<std::pair<int,int>, gsl_spline*> > opacity_lin_interp_op;
    std::map<std::pair<float,float>, gsl_interp_accel*> opacity_acc_tops;
    std::map<std::pair<float,float>, gsl_spline*> opacity_lin_interp_tops;
    std::map<double, gsl_interp_accel*> opacity_acc_opas;
    std::map<double, gsl_spline*> opacity_lin_interp_opas;
    std::vector<std::vector<float>> tops_grid;
    std::vector<float> tops_temperatures;
    std::vector<float> tops_densities;
    std::vector<gsl_interp_accel*> n_isotope_acc;
    std::vector<gsl_spline*> n_isotope_lin_interp;
    std::vector<gsl_interp_accel*> z2_n_isotope_acc;
    std::vector<gsl_spline*> z2_n_isotope_lin_interp;
    std::map<std::string, gsl_interp_accel*> n_element_acc;
    std::map<std::string, gsl_spline*> n_element_lin_interp;
};

#endif // defined __utils_hpp__
