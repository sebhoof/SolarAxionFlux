cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(SolarAxionFlux VERSION 1.1 DESCRIPTION "Calculate the solar axion flux from axion-photon and axion-electron interactions" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

# Must use GNUInstallDirs to install libraries into correct locations on all platforms.
include(GNUInstallDirs)

# Set preprocessor macros
add_compile_definitions(SOLAXFLUX_DIR="${PROJECT_SOURCE_DIR}")
add_compile_definitions(LIBRARY_NAME="SolarAxionFlux v${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

# Set some compiler flags for non-release-type builds (okay for a small library)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Check for GNU Scientific Library (GSL) and add it
find_package(GSL 2.3 REQUIRED)
include_directories("${GSL_INCLUDE_DIRS}")
set(LIBRARIES ${LIBRARIES} ${GSL_LIBRARIES})

# Add an option to disable Python support
option(PYTHON_SUPPORT "Adds a Python module library using the PYBIND11 headers." ON)

# Check for Python 3 and PYBIND11 (if desired)
set(INTERNAL_PYBIND11_DIR "${PROJECT_SOURCE_DIR}/external/pybind11")
if(PYTHON_SUPPORT)
  find_package(Python3 COMPONENTS Interpreter Development QUIET)
  if(Python3_FOUND)
    # Include FetchContent to clone PYBIND11 from git repo
    include(FetchContent)
    FetchContent_Declare(pybind11
                         GIT_REPOSITORY https://github.com/pybind/pybind11.git
                         GIT_TAG v2.10.4
                         GIT_PROGRESS ON
                         GIT_SHALLOW ON
                         GIT_SUBMODULES ""
                         GIT_CONFIG advice.detachedHead=False
                         SOURCE_DIR ${INTERNAL_PYBIND11_DIR})

    if(EXISTS ${INTERNAL_PYBIND11_DIR})
      message("-- Python ${Python3_VERSION} and pybind11 found!")
      set(PYBIND11_INSTALLED 1)
      add_subdirectory(${INTERNAL_PYBIND11_DIR})
    else()
      message("-- Python ${Python3_VERSION} found! Will now download pybind11 header files...")
      FetchContent_Populate(pybind11)
      add_subdirectory(${INTERNAL_PYBIND11_DIR})
    endif()
  else()
    message("-- Optional program Python 3 not found. Skipping pybind11 integration...")
  endif()
else()
  message("-- PYTHON_SUPPORT has been set to OFF. Skipping pybind11 integration...")
endif()

add_subdirectory(src)

# Generate thisLibrary scripts
include(cmake/thisSolarAxionFluxLib.cmake)
