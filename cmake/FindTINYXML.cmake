# Search for tinyxml2 via find_package(tinyxml2)
include(CMakeFindDependencyMacro)
find_dependency(tinyxml2 NO_MODULE)

# If this does not work, search for tinyxml2 via pkg-config
if(NOT tinyxml2_FOUND OR NOT TARGET tinyxml2::tinyxml2)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(tinyxml2 QUIET IMPORTED_TARGET tinyxml2)
    if(TARGET PkgConfig::tinyxml2 AND NOT TARGET tinyxml2::tinyxml2)
      # Define tinyxml2::tinyxml2 so in the code we just use tinyxml2::tinyxml2
      # as if the package was found via find_package(tinyxml2)
      add_library(tinyxml2::tinyxml2 INTERFACE IMPORTED)
      # Equivalent to target_link_libraries INTERFACE, but compatible with CMake 3.10
      set_target_properties(tinyxml2::tinyxml2 PROPERTIES INTERFACE_LINK_LIBRARIES PkgConfig::tinyxml2)
      set(tinyxml2_FOUND TRUE)
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TINYXML DEFAULT_MSG tinyxml2_FOUND)
