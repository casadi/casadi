# - Try to find FICO Xpress Optimizer
#
# Searches for the Xpress C API (libxprs / xprs.dll, xprs.h).
# A typical Xpress install sets the XPRESSDIR environment variable
# (the directory containing include/ and lib/).  XPRESS_HOME is also accepted.
#
# Sets:
#  XPRESS_FOUND         - True if Xpress was found
#  XPRESS_INCLUDE_DIRS  - Path to xprs.h
#  XPRESS_LIBRARIES     - The xprs library to link
#  Target xpress::xpress

find_path(XPRESS_INCLUDE_DIR
  NAMES xprs.h
  HINTS
    ENV XPRESSDIR
    ENV XPRESS_HOME
    ENV XPRESS_DIR
  PATH_SUFFIXES include
  PATHS
    "/opt/xpressmp"
    "/opt/fico/xpressmp"
    "C:/xpressmp"
    "C:/Program Files/FICO/Xpress/xpressmp"
)

find_library(XPRESS_LIBRARY
  NAMES xprs
  HINTS
    ENV XPRESSDIR
    ENV XPRESS_HOME
    ENV XPRESS_DIR
  PATH_SUFFIXES lib
  PATHS
    "/opt/xpressmp"
    "/opt/fico/xpressmp"
    "C:/xpressmp"
    "C:/Program Files/FICO/Xpress/xpressmp"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XPRESS DEFAULT_MSG
  XPRESS_LIBRARY XPRESS_INCLUDE_DIR)

if(XPRESS_FOUND)
  set(XPRESS_INCLUDE_DIRS "${XPRESS_INCLUDE_DIR}")
  set(XPRESS_LIBRARIES "${XPRESS_LIBRARY}")

  if(NOT TARGET xpress::xpress)
    add_library(xpress::xpress INTERFACE IMPORTED)
    set_target_properties(xpress::xpress PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${XPRESS_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${XPRESS_LIBRARY}"
    )
  endif()
endif()

mark_as_advanced(XPRESS_INCLUDE_DIR XPRESS_LIBRARY)
