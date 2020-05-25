set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
set(LIB_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/lib)
set(HEADER_FILES_DIR ${CMAKE_SOURCE_DIR}/include/solaxflux)
set(HEADER_FILES ${HEADER_FILES_DIR}/spectral_flux.hpp ${HEADER_FILES_DIR}/experimental_flux.hpp ${HEADER_FILES_DIR}/solar_model.hpp ${HEADER_FILES_DIR}/utils.hpp ${HEADER_FILES_DIR}/constants.hpp)

add_library(axionflux SHARED experimental_flux.cpp spectral_flux.cpp solar_model.cpp utils.cpp ${HEADER_FILES})
set_target_properties(axionflux PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})
# set_target_properties(axionflux PROPERTIES PUBLIC_HEADER ${HEADER_FILES})
target_include_directories(axionflux PRIVATE ${CMAKE_SOURCE_DIR}/include/solaxflux ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(axionflux ${LIBRARIES})
set_target_properties(axionflux PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIB_OUTPUT_DIR})

# Disable global install steps for now...
# install(TARGETS axionflux EXPORT MyLibraryConfig
#    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})

# # This makes the project importable from the install directory
# # Put config file in per-project dir (name MUST match), can also
# # just go into 'cmake'.
# install(EXPORT MyLibraryConfig DESTINATION share/MyLibrary/cmake)

# # This makes the project importable from the build directory
# export(TARGETS axionflux FILE MyLibraryConfig.cmake)

# Create executable for the unit tests
add_executable(test_library tests.cpp ${HEADER_FILES})
target_include_directories(test_library PRIVATE ${CMAKE_SOURCE_DIR}/include/solaxflux ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_library PUBLIC axionflux)
set_target_properties(test_library PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Create the PYBIND11 library
if(PYBIND11_INSTALLED)
  pybind11_add_module(pyaxionflux python_wrapper.cpp)
  target_include_directories(pyaxionflux PRIVATE ${CMAKE_SOURCE_DIR}/include/solaxflux ${CMAKE_SOURCE_DIR}/include/pybind11 ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(pyaxionflux PUBLIC axionflux)
  set_target_properties(pyaxionflux PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIB_OUTPUT_DIR})
endif()
