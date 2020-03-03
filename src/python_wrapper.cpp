#include <pybind11/pybind11.h>

#include "spectral_flux.hpp"
#include "experimental_flux.hpp"

PYBIND11_MODULE(pyaxionflux, m) {
    m.doc() = "Solar Axion Flux library functions";

    m.def("conversion_prob_correction", &conversion_prob_correction, "The axion-photon conversion probability mass correction.");
}
