cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(SolarAxionFlux VERSION 0.1 DESCRIPTION "Calculate the solar axion flux" LANGUAGES CXX)
if (CMAKE_VERSION VERSION_LESS "3.1")
   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
   endif ()
else ()
   set (CMAKE_CXX_STANDARD 11)
endif ()


# Must use GNUInstallDirs to install libraries into correct locations on all platforms.
include(GNUInstallDirs)

# Check for GNU Scientific Library (GSL) and add it
include(cmake/FindGSL.cmake)
if(GSL_FOUND)
   if (NOT GSL_INCLUDE_DIRS STREQUAL "")
      include_directories("${GSL_INCLUDE_DIRS}")
   endif()
else()
   message("GSL libraries not found.")
endif()
set(LIBRARIES ${LIBRARIES} ${GSL_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory(src)
