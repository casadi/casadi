# Get package info using pkg-config
find_package(PkgConfig)
pkg_search_module(IPOPT ipopt)

# Callback support
if(IPOPT_INCLUDEDIR)
  if(EXISTS ${IPOPT_INCLUDEDIR}/IpIpoptData.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpOrigIpoptNLP.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpTNLPAdapter.hpp  AND EXISTS ${IPOPT_INCLUDEDIR}/IpDenseVector.hpp AND EXISTS ${IPOPT_INCLUDEDIR}/IpExpansionMatrix.hpp)
    set(IPOPT_FOUND_EXTRA TRUE)
  else ()
    message(STATUS "Detected an IPOPT configuration without development headers. Build will proceed, but without callback functionality. To enable it, see https://sourceforge.net/apps/trac/casadi/wiki/enableIpoptCallback")
  endif ()
endif(IPOPT_INCLUDEDIR)

# Find sIPOPT library and add this to the list of libraries
if(IPOPT_LIBRARIES)
  find_library(SIPOPT_LIBRARY sipopt HINTS ${IPOPT_LIBDIR})
  if(SIPOPT_LIBRARY)
#       set(IPOPT_LIBRARIES "${SIPOPT_LIBRARY} ${IPOPT_LIBRARIES}")
    set(WITH_SIPOPT TRUE)
  else ()
    message(STATUS "Detected an IPOPT configuration without support for parametric sensitivities (sIPOPT). Build will proceed, but without sIPOPT functionality.")
  endif ()
endif(IPOPT_LIBRARIES)

# Set standard flags
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS)
