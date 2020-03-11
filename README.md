# Solar Axion Flux

<em><font size="4">A C++ library to calculate the expected flux from axion-photon and axion-electron interactions inside the Sun.</font></em>

_**This code is currently still under development!**_

## Example results

Here, we show a number of plots comparing different Solar models and opacity codes. More details and references will be given in a future publication.

## Installation

The code is written in C++ and contains Python files for data processing tasks. It only depends on the GSL library and, optionally, on Python 3.x to build Python wrappers (currently experimental and unsupported). We re-distribute (in adjusted form) Solar models and opacity tables, which should be acknoledged appropriately using the references stated below.

A quick guide on how to install and test the code.
* For Mac OS use e.g. [Homebrew](https://brew.sh) to install the GSL library via `brew install gsl`. For Linux use `sudo apt-get install libgsl-dev` instead. If you do not have admin privileges on either system, you need to [install the GSL library from source](https://www.gnu.org/software/gsl/).
* Clone this repo via `git clone https://github.com/sebhoof/SolarAxionFlux [foldername]`, where `[foldername]` can be replaced by a folder name of your choice.
* Follow this up by `cd [foldername]; mkdir build; cd build/`.
* In most cases `cmake ..` and then `make` should build everything. If this fails, try to specify the compiler that you want to use via `cmake -D CMAKE_CXX_COMPILER=[compiler executable name or path] ..`.
* The `test_axionflux` executable in the `bin` directory runs a simple test program. It needs to be exectued from the `[foldername]` directory to find the Solar models located in `data/`, i.e. `./bin/test_axionflux`.

## References

For using the Solar model files and opacity code results, please cite the following works

### The code

We currently do not have a DOI or reference to go with this code. Please contact the developers if you wish to acknowledge this code in your work.

### Solar models

### Opacity codes
