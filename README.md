# SolarAxionFlux

Code to calculate the expected flux from axion-photon and axion-electron conversion inside the Sun.

## Results

## Installation

The code is written in C++ and contains Python files for data processing tasks. It only depends on the GSL library. While we provide results for existing Solar models and opacity tables, these need to be provided by the user in a specific format to use the full computational chain.

Example for how to install the code on Mac OS:
* Use e.g. `brew install gsl` to install the GSL library
* Clone this repo via `git clone https://github.com/sebhoof/SolarAxionFlux [foldername]`
* Do `mkdir build; cd build/` in the project directory
* In most cases `cmake ..` and then `make` should build everything
* Do `./bin/test_axionflux` in the project folder to run a test
