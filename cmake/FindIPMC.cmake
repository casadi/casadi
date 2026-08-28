# Locate ipmc solver
#
# Sets:
#   IPMC_FOUND
#   IPMC_INCLUDE_DIR     directory to add to the include path
#   IPMC_LIBRARY         absolute path to the ipmc shared object
#   ipmc::ipmc        imported target

set(_ipmc_hints)
foreach(_h "${IPMC_ROOT}" "${IPMC_DIR}" "$ENV{IPMC_ROOT}" "$ENV{IPMC_DIR}")
  if(_h)
    list(APPEND _ipmc_hints "${_h}" "${_h}/build")
  endif()
endforeach()

if(_ipmc_hints)
  set(_ipmc_extra NO_DEFAULT_PATH)
else()
  set(_ipmc_extra)
endif()

# The public header is <root>/ipmc/ipmc.h, so the include directory is the
# root of the ipmc checkout itself.
find_path(IPMC_INCLUDE_DIR
  NAMES ipmc/ipmc.h
  HINTS ${_ipmc_hints}
  PATH_SUFFIXES include
  ${_ipmc_extra}
)

find_library(IPMC_LIBRARY
  NAMES ipmc
  HINTS ${_ipmc_hints}
  PATH_SUFFIXES ipmc lib lib64
  ${_ipmc_extra}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPMC
  REQUIRED_VARS IPMC_LIBRARY IPMC_INCLUDE_DIR
)

if(IPMC_FOUND AND NOT TARGET ipmc::ipmc)
  add_library(ipmc::ipmc UNKNOWN IMPORTED)
  set_target_properties(ipmc::ipmc PROPERTIES
    IMPORTED_LOCATION "${IPMC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${IPMC_INCLUDE_DIR}"
  )
  # blasfeo is attached by the caller, not here: this module is run from the
  # proprietary-package loop in the top-level CMakeLists, which is well before
  # the blasfeo target exists.  See the WITH_IPMC block there.
endif()

mark_as_advanced(IPMC_INCLUDE_DIR IPMC_LIBRARY)
