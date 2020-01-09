#ifndef __utils_hpp__
#define __utils_hpp__

#include <iostream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include "constants.hpp"

void terminate_with_error(std::string err_string);

const double abs_prec = 1.0e-4, rel_prec = 1.0e-6;
const int method = 5;
struct integrand_params {double u; double y;};

// Map of element name -> Z value
const std::map<std::string, int> el_z_value ({{"H", 1}, {"He", 2}, {"C", 6}, {"N", 7}, {"O", 8}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17},
                                              {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}});;

// Note: For future versions, use a list of 'Element's that are available in a SolarModel
// el_a_value=0 is a default for 'only tracked in bulk';
// Columns is Solar model can stay as they are.
// -> Variable iz becomes type 'Element'
// For the OP code, just use the op list with names(!) as below and automatically sum over all elements with iz[0] == ElementName!
// If a_val() > 0, use this, else get default bulk value from global table (maybe create own namespace).
class Element {
   public:
      Element() {};
      Element(std::string s, int a) {el_name = s; el_a_value = a; };
      Element(std::pair<std::string,int> p) {el_name = p.first; el_a_value = p.second; };
      Element(std::string s) {el_name = s; el_a_value = 0; };

      std::string name() { return el_name; };
      int a_val() { return el_a_value; };
      int z_val() { return el_z_value.at(el_name); };
   private:
      std::string el_name;
      int el_a_value;
};
enum opacitycode {OP, OPAS, LEDCOP, ATOMIC};
const int op_grid [197][2] = {{150,54}, {150,56}, {152,54}, {152,56}, {160,64}, {160,66}, {162,64}, {162,66}, {164,66}, {164,68}, {166,66}, {166,68}, {166,70}, {168,68}, {168,70}, {168,72}, {168,74}, {170,70}, {170,72}, {170,74}, {172,72}, {172,74}, {172,76}, {174,74}, {174,76}, {174,78}, {176,74}, {176,76}, {176,78}, {178,76}, {178,78}, {180,76}, {180,78}, {180,80}, {182,78}, {182,80}, {184,78}, {184,80}, {186,78}, {186,80}, {186,82}, {188,80}, {188,82}, {190,80}, {190,82}, {190,84}, {192,80}, {192,82}, {192,84}, {194,80}, {194,82}, {194,84}, {196,80}, {196,82}, {196,84}, {198,82}, {198,84}, {198,86}, {200,82}, {200,84}, {200,86}, {202,82}, {202,84}, {202,86}, {204,82}, {204,84}, {204,86}, {206,82}, {206,84}, {206,86}, {208,84}, {208,86}, {208,88}, {210,84}, {210,86}, {210,88}, {212,84}, {212,86}, {212,88}, {214,84}, {214,86}, {214,88}, {216,84}, {216,86}, {216,88}, {218,86}, {218,88}, {220,86}, {220,88}, {220,90}, {222,86}, {222,88}, {222,90}, {224,86}, {224,88}, {224,90}, {226,86}, {226,88}, {226,90}, {228,86}, {228,88}, {228,90}, {230,86}, {230,88}, {230,90}, {232,88}, {232,90}, {232,92}, {234,88}, {234,90}, {234,92}, {236,88}, {236,90}, {236,92}, {238,88}, {238,90}, {238,92}, {240,88}, {240,90}, {240,92}, {242,88}, {242,90}, {242,92}, {244,90}, {244,92}, {246,90}, {246,92}, {246,94}, {248,90}, {248,92}, {248,94}, {250,90}, {250,92}, {250,94}, {252,90}, {252,92}, {252,94}, {254,90}, {254,92}, {254,94}, {256,90}, {256,92}, {256,94}, {256,96}, {258,92}, {258,94}, {258,96}, {260,92}, {260,94}, {260,96}, {260,98}, {262,92}, {262,94}, {262,96}, {262,98}, {264,94}, {264,96}, {264,98}, {266,94}, {266,96}, {266,98}, {266,100}, {268,96}, {268,98}, {268,100}, {270,96}, {270,98}, {270,100}, {270,102}, {272,96}, {272,98}, {272,100}, {272,102}, {274,98}, {274,100}, {274,102}, {276,98}, {276,100}, {276,102}, {278,100}, {278,102}, {278,104}, {280,100}, {280,102}, {280,104}, {282,100}, {282,102}, {282,104}, {284,100}, {284,102}, {284,104}, {286,102}, {286,104}, {286,106}, {288,102}, {288,104}, {288,106}};

const float ledcop_grid [26][2] = {{0.5,10.175}, {0.5,13.684}, {0.6,10.175}, {0.6,13.684}, {0.6,18.308}, {0.6,24.268}, {0.6,31.802}, {0.8,18.308}, {0.8,24.268}, {0.8,31.802}, {0.8,41.156}, {0.8,52.611}, {1.0,41.156}, {1.0,52.611}, {1.0,66.544}, {1.0,83.466}, {1.0,103.442}, {1.25,66.544}, {1.25,83.466}, {1.25,103.442}, {1.25,124.995}, {1.25,143.111}, {1.25,150.5}, {1.5,124.995}, {1.5,143.111}, {1.5,150.5}};
const float ledcop_temperatures [13] = {0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5};
const float ledcop_densities [13] = {10.175, 13.684, 18.308, 24.268, 31.802, 41.156, 52.611, 66.544, 83.466, 103.442, 124.995, 143.111, 150.5};

const int n_op_elements = 17;
const std::string op_elements_simple [n_op_elements] = {"H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "Ca", "Cr", "Mn", "Fe", "Ni"};
const std::vector<std::vector<int>> op_elements = {{0}, {1,2}, {3,4}, {5,6}, {7,8,9}, {10}, {11}, {12}, {13}, {14}, {16}, {18}, {20}, {24}, {25}, {26}, {28}};

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
    SolarModel(std::string file, opacitycode set_opcode,const bool set_raffelt_approx);
    SolarModel(std::string file, opacitycode set_opcode);
    ~SolarModel();
    SolarModel& operator=(SolarModel&&);
    // Delete copy constructor and assignment operator to avoid shallow copies
    SolarModel(const SolarModel&) = delete;
    SolarModel operator=(const SolarModel&) = delete;
    // Min. and max. radius of the solar model file (distance r from the centre of the Sun in units of the solar radius)
    bool raffelt_approx;  // whether to use the approximation by Raffelt (https://wwwth.mpp.mpg.de/members/raffelt/mypapers/198601.pdf), default is false
    opacitycode opcode;
    double r_lo, r_hi;
    // Routine to return the screening parameter kappa^2 (kappa^-1 = Debye-Hueckel radius).
    double kappa_squared(double r);
    // Routine to return the temperature in keV.
    double temperature_in_keV(double r);
    double density(double r);
    // ...
    double z2_n_iz(double r, int iz);
    double z2_n(double r);
    // ...
    double n_iz(double r, int iz);
    // Raffelt approximation
    double n_e(double r);
    // Routine to return the plasma frequency squared.
    double omega_pl_squared(double r);

    // ...
    double Gamma_P_ff(double omega, double r, int iz);
    double Gamma_P_ff(double omega, double r);
    double Gamma_P_ee(double omega, double r);
    double Gamma_P_Compton (double omega, double r);
    double op_grid_interp_erg (double u, int ite, int jne, int iz);
    double ledcop_grid_interp_erg (double erg, float T, float rho);
    double opacity_table_interpolator (double omega, double r, int iz);
    double opacity_table_interpolator (double omega, double r);
    double opacity_table_interpolator2 (double omega, double r, int iz);
    double opacity_element (double omega, double r, int iz);
    double opacity (double omega, double r);
    double Gamma_P_element (double omega, double r, int iz);
    double Gamma_P_opacity (double omega, double r);
    double Gamma_P_Primakoff (double omega, double r);
  private:
    ASCIItableReader data;
    gsl_interp_accel *accel[5];
    gsl_spline *linear_interp[5];
    std::vector< std::map<std::pair<int,int>, gsl_interp_accel*> > opacity_acc_op;
    std::vector< std::map<std::pair<int,int>, gsl_spline*> > opacity_lin_interp_op;
    std::map<std::pair<float,float>, gsl_interp_accel*> opacity_acc_ledcop;
    std::map<std::pair<float,float>, gsl_spline*> opacity_lin_interp_ledcop;
    std::vector<gsl_interp_accel*> n_iz_acc;
    std::vector<gsl_spline*> n_iz_lin_interp;
    std::vector<gsl_interp_accel*> z2_n_iz_acc;
    std::vector<gsl_spline*> z2_n_iz_lin_interp;
    
};

#endif // defined __utils_hpp__
