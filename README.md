# Solar Axion Flux

<em><font size="4">A C++ library to calculate the expected flux from axion-photon and axion-electron interactions inside the Sun.</font></em>

_**This code is currently under development by [Sebastian Hoof](mailto:hoof@uni-goettingen.de) and [Lennert Thormaehlen](mailto:l.thormaehlen@thphys.uni-heidelberg.de).**_

This code has been published under the BSD 3-clause license. Consult the [LICENSE](LICENSE) file for details.

## Example results

Here, we show a number of plots comparing different Solar models and opacity codes. More details and references will be given in a future publication.

## Installation

The code is written in C++ and contains Python files for data processing tasks. It only depends on the GSL library and, optionally, on Python v3.x to build Python wrappers with pybind11 (requires cython).

A complete guide on how to install and test the code (requires CMAKE v3.12 or higher) is given below:
* For Mac OS: use e.g. [Homebrew](https://brew.sh) to install CMAKE via `brew install cmake`.
* For Mac OS: use e.g. [Homebrew](https://brew.sh) to install the GSL library via `brew install gsl`. For Linux: use `sudo apt-get install libgsl-dev` instead. If you do not have admin privileges on either system, you need to [install the GSL library from source](https://www.gnu.org/software/gsl/). We recommend using GSL library v2.4 or higher.
* Clone this repo via `git clone https://github.com/sebhoof/SolarAxionFlux [foldername]`, where `[foldername]` can be replaced by a folder name of your choice.
* Use the latest `master` branch (no need to do anything) or checkout the latest tagged version as a new branch (`git checkout v0.5b -b [some_branch_name]`).
* Set up a directory via `cd [foldername]`, `mkdir build`, and `cd build/`.
* In most cases `cmake ..` and then `make` should build everything. If this fails, consult the [Troubleshooting](#troubleshooting) section.

## Testing
The `test_library` executable in the `bin/` directory runs a simple test program. If you installed the Python frontend, you can also run this test from your Python 3 terminal or notebook via `from lib import pyaxionflux as afl` (assuming you are in `[foldername]` or added `[foldername]` to your `$PYTHONPATH` variable) and `afl.test_module()`. In either case, the output can be found in the `[foldername]/results/` folder.

## References

We re-distribute (in adjusted form) Solar models and opacity tables, which should be acknowledged appropriately using the references stated below. Please check for updates as this list is being updated during the development stage.

### The code

We currently do not have a DOI or reference to go with this code. Please contact us if you wish to acknowledge this code in your work.

### Solar models
* BP98 [arXiv:astro-ph/9805135](https://arxiv.org/astro-ph/abs/astro-ph/9805135)
* BP00 [arXiv:astro-ph/0010346](https://arxiv.org/astro-ph/abs/astro-ph/0010346)
* BP04 [arXiv:astro-ph/0402114](https://arxiv.org/astro-ph/abs/astro-ph/0402114)
* BS05-OP, BS05-AGSOP [arXiv:astro-ph/0412440](https://arxiv.org/astro-ph/abs/astro-ph/0412440)
* AGS05 [arXiv:0909.2668](https://arxiv.org/astro-ph/abs/0909.2668)
* GS98, AGSS09(met), AGSS09ph [arXiv:0909.2668](https://arxiv.org/astro-ph/abs/0909.2668), [arXiv:0910.3690](https://arxiv.org/astro-ph/abs/0910.3690)
* B16-GS98, B16-AGSS09 [arXiv:1611.09867](https://arxiv.org/astro-ph/abs/1611.09867)

### Opacities
* LEDCOP [APS Conf. Series **75** (1995)](https://ui.adsabs.harvard.edu/abs/1995ASPC...78...51M)
* OP [arXiv:astro-ph/0410744](https://arxiv.org/astro-ph/abs/astro-ph/0410744), [arXiv:astro-ph/0411010](https://arxiv.org/astro-ph/abs/astro-ph/0411010)
* OPAS [ApJ **754** 1 (1012)](https://doi.org/10.1088/0004-637X/745/1/10), [ApJ Suppl. Series **220** 1 (1015)](https://doi.org/10.1088/0067-0049/220/1/2)
* ATOMIC [arXiv:1601.01005](https://arxiv.org/astro-ph/abs/1601.01005)

## Troubleshooting
* "I get some compiler related error." Try specifying the compiler that you want to use via `cmake -D CMAKE_CXX_COMPILER=[compiler executable name or path] ..`
* "CMAKE can't find the GSL library." You can give CMAKE a hint of where to find the desired version of the GSL library via `GSL_ROOT=[path to GSL folder] cmake ..`
* "I get some pybind11 error during CMAKE." Usually, there are multiple Python versions installed (by the system, package managers, virtual environments, ...), which are likely to be the root of the problem. Check that CMAKE recognises the desired Python version that you want to use. You man need to install `cython` via `pip`. The Python frontend can be disable via `cmake -D PYTHON_SUPPORT=OFF ..`.
* "I managed to build the Python library, but I can't run the test." Check the name of the library in the `lib/` folder (e.g. `pyaxionflux.cpython-38-darwin.so`) which indicates the Python version used to build it (in this case, Python 3.8), which may help identifying the appropriate executable on your system.
