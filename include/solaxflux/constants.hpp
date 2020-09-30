// Copyright 2020 Sebastian Hoof & Lennert J. Thormaehlen
// See the LICENSE file for the license conditions and a disclaimer

#ifndef __constants_hpp__
#define __constants_hpp__

#include <string>

#include <gsl/gsl_math.h>

// Library name
const std::string LIBRARY_NAME = "SolarAxionFlux v0.6b";

// Mathematical constants
const double pi = M_PI;

// Physical constants
const double hbar = 6.582119514e-25; // in [GeV s]
const double radius_sol = 6.957e8; // [m]; nominal value of the Solar radius, as recommend by https://www.iau.org/static/resolutions/IAU2015_English.pdf
const double distance_sol = 1.49597870700e11; // [m]; taken to be 1 au (with the value of 1 au as recommnded by https://www.iau.org/static/resolutions/IAU2012_English.pdf)
const double alpha_EM = 7.2973525664e-3; // fine structure constant
const double atomic_mass_unit = 0.931494028; // [GeV]
const double gev2cm = 197.327053e-16; // [cm / GeV^-1]
const double eV2g = 1.782661907e-33; // [eV / g]
const double s2cm = 2.99792458e10; // [cm / s]
const double keV2cm = 197.327053e-10; // [cm / keV^-1]
const double eVm = gev2cm*1.0e7; // [m / eV^-1]
const double K2eV = 8.6173303e-5; // [eV / K]
const double m_electron = 0.5109989461e3; // electron mass [keV]
const double a_Bohr = 5.29177210903e-9; // Bohr radius [cm]
const double g_aee = 1.0E-13; // default axion-electron coupling
const double g_agg = 1.0E-16; // default axion-photon coupling [keV^-1]
const double eV2T = sqrt(4.0*pi)*1.4440271*1.0e-3; // [T / eV^2]

#endif // defined __constants_hpp__
