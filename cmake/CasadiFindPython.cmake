# Locate Python in a way that works across the CMake versions we support.
#
# FindPythonInterp / FindPythonLibs are deprecated.  On CMake >= 3.18 we use
# the modern FindPython3 with the header-only Development.Module component, so
# manylinux / abi3 module builds don't require libpython to be present.  On
# older CMake we fall back to the legacy modules.  Either way the legacy
# variables the rest of the build relies on are populated:
#   PYTHONINTERP_FOUND, PYTHON_EXECUTABLE, PYTHON_VERSION_STRING,
#   PYTHON_INCLUDE_PATH, PYTHON_INCLUDE_DIRS, PYTHON_LIBRARY, PYTHON_LIBRARIES
#
# Usage:
#   casadi_find_python()                       # interpreter only
#   casadi_find_python(DEVELOPMENT)            # + headers (and lib on Windows)
#   casadi_find_python(DEVELOPMENT REQUIRED)   # error out if not found
#   casadi_find_python(QUIET)                  # suppress status messages

include_guard(GLOBAL)

macro(casadi_find_python)
  set(_cfp_dev OFF)
  set(_cfp_flags "")
  foreach(_cfp_arg ${ARGN})
    if(_cfp_arg STREQUAL "DEVELOPMENT")
      set(_cfp_dev ON)
    elseif(_cfp_arg STREQUAL "REQUIRED" OR _cfp_arg STREQUAL "QUIET")
      list(APPEND _cfp_flags ${_cfp_arg})
    endif()
  endforeach()

  if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
    if(_cfp_dev)
      find_package(Python3 COMPONENTS Interpreter Development.Module ${_cfp_flags})
    else()
      find_package(Python3 COMPONENTS Interpreter ${_cfp_flags})
    endif()
    # Map the modern results onto the legacy variable names.
    set(PYTHONINTERP_FOUND ${Python3_Interpreter_FOUND})
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
    set(PYTHON_VERSION_STRING ${Python3_VERSION})
    if(_cfp_dev)
      set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS})
      set(PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
      set(PYTHON_LIBRARY ${Python3_LIBRARIES})
      set(PYTHON_LIBRARIES ${Python3_LIBRARIES})
    endif()
  else()
    find_package(PythonInterp "3" ${_cfp_flags})
    if(_cfp_dev)
      find_package(PythonLibs "3" ${_cfp_flags})
    endif()
  endif()
endmacro()
