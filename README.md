# Solar Axion Flux

<em><font size="4">SolarAxionFlux library to calculate the expected flux from axion-photon and axion-electron conversion inside the Sun.</font></em>

## Results

Plots, etc.

## Installation

The code is written in C++ and contains Python files for data processing tasks. It only depends on the GSL library. While we provide results for existing Solar models and opacity tables, these need to be provided by the user in a specific format to use the full computational chain.

A quick guide on how to install and test the code.
* For Mac OS use e.g. [Homebrew](https://brew.sh) to install the GSL library via `brew install gsl`. For Linux use `sudo apt-get install libgsl-dev` instead. If you have no admin privileges, you need to [install the GSL library from source](https://www.gnu.org/software/gsl/).
* Clone this repo via `git clone https://github.com/sebhoof/SolarAxionFlux [foldername]`, where `[foldername]` can be replaced with a folder name of your choice.
* Follow this up by `cd [foldername]; mkdir build; cd build/`.
* In most cases `cmake ..` and then `make` should build everything.
* The `test_axionflux` executable in the `bin` directory runs a simple test program.
