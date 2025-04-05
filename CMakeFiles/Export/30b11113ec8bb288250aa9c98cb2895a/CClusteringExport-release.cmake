#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "CClustering" for configuration "Release"
set_property(TARGET CClustering APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(CClustering PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/CClustering/libCClustering.so"
  IMPORTED_SONAME_RELEASE "libCClustering.so"
  )

list(APPEND _cmake_import_check_targets CClustering )
list(APPEND _cmake_import_check_files_for_CClustering "${_IMPORT_PREFIX}/lib/CClustering/libCClustering.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
