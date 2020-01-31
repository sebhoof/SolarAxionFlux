set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(HEADER_FILES ${CMAKE_SOURCE_DIR}/include/spectral_flux.hpp ${CMAKE_SOURCE_DIR}/include/experimental_flux.hpp ${CMAKE_SOURCE_DIR}/include/utils.hpp ${CMAKE_SOURCE_DIR}/include/constants.hpp)

add_library(axionflux SHARED spectral_flux.cpp experimental_flux.cpp utils.cpp ${HEADER_FILES})
target_include_directories(axionflux PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:../include>)
target_link_libraries(axionflux ${LIBRARIES})

install(TARGETS axionflux EXPORT MyLibraryConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT MyLibraryConfig DESTINATION share/MyLibrary/cmake)

# This makes the project importable from the build directory
export(TARGETS axionflux FILE MyLibraryConfig.cmake)

# This is for the unit tests
add_executable(test_axionflux main.cpp ${HEADER_FILES})
target_include_directories(test_axionflux PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(test_axionflux axionflux)
set_target_properties(test_axionflux
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

add_executable(compute_all_spectra compute_all_spectra.cpp ${HEADER_FILES})
target_include_directories(compute_all_spectra PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(compute_all_spectra axionflux)
set_target_properties(compute_all_spectra
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)
