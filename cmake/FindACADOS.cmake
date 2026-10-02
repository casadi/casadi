# Search for acados via its acadosConfig.cmake, installed in <prefix>/cmake
if(NOT acados_DIR)
  foreach(ACADOS_PREFIX $ENV{ACADOS_INSTALL_DIR} $ENV{ACADOS_SOURCE_DIR})
    if(EXISTS "${ACADOS_PREFIX}/cmake/acadosConfig.cmake")
      set(acados_DIR "${ACADOS_PREFIX}/cmake")
      break()
    endif()
  endforeach()
endif()
find_package(acados CONFIG QUIET)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACADOS DEFAULT_MSG acados_FOUND)
