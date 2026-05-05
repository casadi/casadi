# - Try to find MOSEK Optimizer
#
# Searches for the MOSEK C API (libmosek64 / mosek64.dll, mosek.h).
# A typical MOSEK install sets the MOSEKDIR or MOSEK_HOME environment
# variable, pointing to the directory that contains include/, h/, and bin/
# (or platform/<...>/h, platform/<...>/bin).
#
# Sets:
#  MOSEK_FOUND         - True if MOSEK was found
#  MOSEK_INCLUDE_DIRS  - Path to mosek.h
#  MOSEK_LIBRARIES     - The mosek library to link
#  Target mosek::mosek

find_path(MOSEK_INCLUDE_DIR
  NAMES mosek.h
  HINTS
    ENV MOSEKDIR
    ENV MOSEK_HOME
    ENV MOSEK_DIR
    ENV MOSEK
  PATH_SUFFIXES h include
  PATHS
    "/opt/mosek"
    "/opt/mosek/10/tools/platform/linux64x86"
    "/opt/mosek/10.1/tools/platform/linux64x86"
    "/opt/mosek/10.2/tools/platform/linux64x86"
    "/opt/mosek/11/tools/platform/linux64x86"
    "$ENV{HOME}/mosek"
    "C:/Program Files/Mosek"
)

find_library(MOSEK_LIBRARY
  NAMES mosek64 mosek
  HINTS
    ENV MOSEKDIR
    ENV MOSEK_HOME
    ENV MOSEK_DIR
    ENV MOSEK
  PATH_SUFFIXES bin lib
  PATHS
    "/opt/mosek"
    "/opt/mosek/10/tools/platform/linux64x86"
    "/opt/mosek/10.1/tools/platform/linux64x86"
    "/opt/mosek/10.2/tools/platform/linux64x86"
    "/opt/mosek/11/tools/platform/linux64x86"
    "$ENV{HOME}/mosek"
    "C:/Program Files/Mosek"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG
  MOSEK_LIBRARY MOSEK_INCLUDE_DIR)

if(MOSEK_FOUND)
  set(MOSEK_INCLUDE_DIRS "${MOSEK_INCLUDE_DIR}")
  set(MOSEK_LIBRARIES "${MOSEK_LIBRARY}")

  if(NOT TARGET mosek::mosek)
    add_library(mosek::mosek INTERFACE IMPORTED)
    set_target_properties(mosek::mosek PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${MOSEK_LIBRARY}"
    )
  endif()
endif()

mark_as_advanced(MOSEK_INCLUDE_DIR MOSEK_LIBRARY)
