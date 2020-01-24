#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "axionflux" for configuration ""
set_property(TARGET axionflux APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(axionflux PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libaxionflux.so"
  IMPORTED_SONAME_NOCONFIG "@rpath/libaxionflux.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS axionflux )
list(APPEND _IMPORT_CHECK_FILES_FOR_axionflux "${_IMPORT_PREFIX}/lib/libaxionflux.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
