cmake_minimum_required(VERSION 3.10.2)

# Workaround for OS X
if(APPLE)
  if(WITH_PYTHON3)
  set(PYTHON_CONFIG_NAME python3-config)
  else()
  set(PYTHON_CONFIG_NAME python-config)
  endif()
  find_program(PYTHON_CONFIG_EXECUTABLE
               NAMES ${PYTHON_CONFIG_NAME} DOC "python-config executable")
  if(PYTHON_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${PYTHON_CONFIG_EXECUTABLE} --prefix
                    OUTPUT_VARIABLE PYTHON_PREFIX_STRING
                    RESULT_VARIABLE PYTHON_PREFIX_FAILED
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT PYTHON_PREFIX_FAILED)
      file(GLOB LOCAL_PYTHON_NAME RELATIVE ${PYTHON_PREFIX_STRING}/lib "${PYTHON_PREFIX_STRING}/lib/python*.*")
      find_library(PYTHON_LIBRARY NAMES ${LOCAL_PYTHON_NAME}
        PATHS ${PYTHON_PREFIX_STRING}/lib NO_DEFAULT_PATH)
      file(GLOB LOCAL_PYTHON_NAME "${PYTHON_PREFIX_STRING}/include/python*")
      find_path(PYTHON_INCLUDE_DIR
        NAMES Python.h
        PATHS ${LOCAL_PYTHON_NAME} NO_DEFAULT_PATH)
    endif()
    message(STATUS ${PYTHON_INCLUDE_DIR})
  endif()
endif()

function(CasadiGetPythonVersionMajor VERSION_MAJOR_RESULT INCLUDE_PATH)
  set(MATCH_STRING "#define[ \t]+PY_MAJOR_VERSION[ \t]+([0-9]+)")
  file(STRINGS "${INCLUDE_PATH}/patchlevel.h" PY_MAJOR_VERSION REGEX ${MATCH_STRING})
  string(REGEX REPLACE "${MATCH_STRING}" "\\1" PYTHON_VERSION_MAJOR "${PY_MAJOR_VERSION}")
  set(${VERSION_MAJOR_RESULT} ${PYTHON_VERSION_MAJOR} PARENT_SCOPE)
endfunction()

if(SWIG_EXPORT)
  if(WITH_PYTHON3)
    set(PYTHON_VERSION_MAJOR 3)
  else()
    set(PYTHON_VERSION_MAJOR 2)
  endif()
else()
  if(PYTHON_INCLUDE_DIR)
    set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIR})
  else()
    # Find packages
    if(WITH_PYTHON3)
    set(MINPYVERSION "3")
    else()
    set(MINPYVERSION "")
    endif()

    # We don't really need the libs, but the headers.
    find_package(PythonInterp ${MINPYVERSION})
    find_package(PythonLibs ${MINPYVERSION})
    if(NOT PYTHON_INCLUDE_PATH)
      message(SEND_ERROR "Python include headers not found, but required." )
    endif()
  endif()
  
  CasadiGetPythonVersionMajor(PYTHON_VERSION_MAJOR "${PYTHON_INCLUDE_PATH}")
  include_directories(${PYTHON_INCLUDE_PATH})
  
  
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)

if(WITH_PYTHON_INTERRUPTS)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()
if(WITH_COPYSIGN_UNDEF)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DCASADI_WITH_COPYSIGN_UNDEF")
endif()

set(PYTHONFLAG "")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DPy_USING_UNICODE")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-noproxydel")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-fastunpack")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-modernargs")

# -stubs emits PEP-484 .pyi files alongside the .py shadow.  Only the
# patched casadi fork of SWIG supports this; vanilla SWIG rejects the
# flag.  Feature-detect so the build still works with stock SWIG
# (without pyright-friendly type hints).  PEP-484 stub syntax is
# inherently Python 3 — never emit .pyi on a Python 2 build.  Under
# SWIG_IMPORT (sdist consumer) SWIG may be absent entirely; skip the
# feature-detect and let the sdist carry the pre-generated .pyi.
set(_emit_stubs OFF)
if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: consuming pre-generated SWIG output; no local stub generation")
elseif(NOT "${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  message(STATUS "Python 2 target; skipping .pyi stub generation (PEP-484 syntax is Python 3 only)")
else()
  execute_process(
    COMMAND ${SWIG_EXECUTABLE} -python -help
    OUTPUT_VARIABLE _swig_python_help
    ERROR_VARIABLE  _swig_python_help_err
    RESULT_VARIABLE _swig_help_rc
  )
  string(FIND "${_swig_python_help}${_swig_python_help_err}" "-stubs" _swig_has_stubs)
  if(NOT _swig_has_stubs EQUAL -1)
    set(_emit_stubs ON)
    set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-stubs")
    message(STATUS "SWIG supports -stubs and Python 3 target; generating PEP-484 type stubs (.pyi)")
  else()
    message(STATUS "SWIG lacks -stubs support; skipping .pyi stub generation. "
                   "For typed Python bindings, build with the casadi-patched SWIG "
                   "fork (see ~/programs/swig or equivalent).")
  endif()
endif()

# Emit the tiny stub-forwarder artifacts in the build dir at configure
# time.  They're one-line files (forwarder + py.typed marker) and
# embedding them here avoids needing to track generated-looking files
# in source control.
set(_stubs_gen_dir ${CMAKE_CURRENT_BINARY_DIR}/stubs)
file(WRITE ${_stubs_gen_dir}/__init__.pyi "from .casadi import *\n")
file(WRITE ${_stubs_gen_dir}/py.typed     "partial\n")
if(WITH_PYTHON_GIL_RELEASE)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-threads")
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DCASADI_WITH_PYTHON_GIL_RELEASE")
endif()
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
set(PYTHONFLAG "3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-py3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON3")
endif()

if(WITH_PYTHON_GIL_RELEASE)
  set(CMAKE_SWIG_NAME_SUFFIX "_gil_release")
else()
  set(CMAKE_SWIG_NAME_SUFFIX "")
endif()

# Main header files
if(NOT SKIP_CONFIG_H_GENERATION)
  configure_file(../swig_config.h.cmake "${CMAKE_CURRENT_BINARY_DIR}/swig_config.h" ESCAPE_QUOTES)
endif()

if(SWIG_EXPORT)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra)
  set(CMAKE_SWIG_WRAP_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/source/dummy)
  
  configure_file(setup.py.in "${PROJECT_SOURCE_DIR}/setup.py" @ONLY)
  file(COPY requirements.txt DESTINATION ${PROJECT_SOURCE_DIR})
  file(COPY MANIFEST.in DESTINATION ${PROJECT_SOURCE_DIR})

  # target*/extra/ is the self-contained Python package surface that a
  # downstream SWIG_IMPORT build consumes.  SWIG writes casadi.py (and
  # casadi.pyi under -stubs) here; we need to also drop in __init__.py
  # and, for Python 3, the stub forwarder + py.typed marker so that
  # the SWIG_IMPORT install step has everything it needs without
  # having to reach back into the source tree.
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
       DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra)
  if(_emit_stubs)
    file(COPY ${_stubs_gen_dir}/__init__.pyi
              ${_stubs_gen_dir}/py.typed
         DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra)
  endif()
endif()

if(SWIG_IMPORT)
  set(PYTHON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/source/casadiPYTHON_wrap${CMAKE_SWIG_NAME_SUFFIX}.cxx)
  set_source_files_properties(${PYTHON_FILE}  PROPERTIES  CPLUSPLUS ON)
else()
  # Generate SWIG wrapper
  set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_module_initialize(casadi python)
  swig_add_source_to_module(casadi FALSE PYTHON_FILE ../casadi.i)
  get_directory_property(swig_extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${PYTHON_FILE}")
endif()

if(SWIG_EXPORT)
  add_custom_target(python_source DEPENDS ${PYTHON_FILE})
endif()

# Make target
add_custom_target(python DEPENDS _casadi)

add_library(_casadi MODULE ${PYTHON_FILE})

set_target_properties(_casadi PROPERTIES PREFIX "")
if(WIN32 AND NOT CYGWIN)
  set_target_properties(_casadi PROPERTIES SUFFIX ".pyd")
  set(CASADI_PYTHON_LIBRARY_SUFFIX ".pyd")
else()
  set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
endif()

if(APPLE AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
  set_target_properties(_casadi PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-register ${MAYBE_WERROR}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR} -Wno-maybe-uninitialized")
  if(WIN32)
    # Fix for error: '_hypot' has not been declared in 'std'
    # cfr https://github.com/tudat-team/tudatpy/commit/ea9ef90a462825a65741f89ac1776413868f52cb
    set(COMPILE_FLAGS "${COMPILE_FLAGS} -D_hypot=hypot")
    # Fix for error: casadiPYTHON_wrap.cxx:(.text+0xd25cc): undefined reference to `__imp_Py_InitModule4'
    # cfr https://github.com/tpaviot/pythonocc-core/commit/9aef6c8e86c5b77c33742acf36e45db0b2f9ca23#diff-1e7de1ae2d059d21e1dd75d5812d5a34b0222cef273b7c3a2af62eb747f9d20a
    if(${CMAKE_SIZEOF_VOID_P} MATCHES "8")
      set(COMPILE_FLAGS "${COMPILE_FLAGS} -DMS_WIN64")
    endif()
  endif()
else()
  set(COMPILE_FLAGS)
endif()

if(MSVC)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} /bigobj")
endif()

message("COMPILE_FLAGS: ${COMPILE_FLAGS}")

set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")

if(WIN32) # dlls don't allow missing symbols
  if(NOT PYTHON_LIBRARY)
    message(SEND_ERROR "Python libraries not found, but required." )
  endif()
  target_link_libraries(_casadi ${PYTHON_LIBRARY})
endif()
target_link_libraries(_casadi casadi)

# Custom installation command for Python
add_custom_target(install_python
COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install_python
  -D CMAKE_INSTALL_PREFIX="${PYTHON_PREFIX}"
  -P cmake_install.cmake
)
add_dependencies(install_python _casadi)

# Install C++ wrapper library
install(TARGETS _casadi
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

if (SWIG_IMPORT)
# Install Python proxy classes + the generated PEP-484 stub (alongside .py).
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra/casadi.pyi
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
  OPTIONAL
)
else()
install(FILES
  ${PROJECT_BINARY_DIR}/swig/python/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)
install(FILES
  ${PROJECT_BINARY_DIR}/swig/python/casadi.pyi
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
  OPTIONAL
)
endif()

# Ship __init__.pyi forwarder + py.typed marker whenever the main
# stub (casadi.pyi) will be installed.  Two cases:
#   - SWIG_IMPORT (pip install from sdist): the sdist producer already
#     made the call — if target*/extra/casadi.pyi exists there, we ship
#     the trio.  The local SWIG may be absent or vanilla; don't gate
#     on `_emit_stubs`.
#   - Normal / SWIG_EXPORT build: use `_emit_stubs`, which tracks
#     whether we just generated casadi.pyi ourselves.
# Python 2 never gets stubs either way (PEP-484 is Py3 syntax).
set(_ship_stub_trio OFF)
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  if(SWIG_IMPORT)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra/casadi.pyi")
      set(_ship_stub_trio ON)
    endif()
  elseif(_emit_stubs)
    set(_ship_stub_trio ON)
  endif()
endif()
if(_ship_stub_trio)
  install(FILES
    ${_stubs_gen_dir}/__init__.pyi
    ${_stubs_gen_dir}/py.typed
    DESTINATION "${PYTHON_PREFIX}/casadi"
    COMPONENT install_python
  )
endif()

if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
  # Install Python tools
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
    DESTINATION "${PYTHON_PREFIX}/casadi"
    COMPONENT install_python
    USE_SOURCE_PERMISSIONS
    PATTERN .pyc EXCLUDE
    PATTERN .svn EXCLUDE
    PATTERN structure.py EXCLUDE
  )
else()
  # Install Python tools
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
    DESTINATION "${PYTHON_PREFIX}/casadi"
    COMPONENT install_python
    USE_SOURCE_PERMISSIONS
    PATTERN .pyc EXCLUDE
    PATTERN .svn EXCLUDE
    PATTERN structure3.py EXCLUDE
  )
endif()

# Install Python package initialization
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

# Install pip metadata files to ensure that casadi installed via CMake is listed by pip list
# See https://packaging.python.org/specifications/recording-installed-packages/
# and https://packaging.python.org/en/latest/specifications/core-metadata/#core-metadata
option(CASADI_PYTHON_PIP_METADATA_INSTALL "Use CMake to install Python pip metadata. Set to off if some other tool already installs it." OFF)
mark_as_advanced(CASADI_PYTHON_PIP_METADATA_INSTALL)
set(CASADI_PYTHON_PIP_METADATA_INSTALLER "cmake" CACHE STRING "Specify the string to identify the pip Installer. Default: cmake, change this if you are using another tool.")
mark_as_advanced(CASADI_PYTHON_PIP_METADATA_INSTALLER)
if (CASADI_PYTHON_PIP_METADATA_INSTALL)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/METADATA "")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Metadata-Version: 2.1\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Name: casadi\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Version: ${CASADI_VERSION}\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/INSTALLER "${CASADI_PYTHON_PIP_METADATA_INSTALLER}\n")
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/METADATA" "${CMAKE_CURRENT_BINARY_DIR}/INSTALLER"
    DESTINATION ${PYTHON_PREFIX}/casadi-${CASADI_VERSION}.dist-info
    COMPONENT install_python)
endif()

# Example of how to extend CasADi with additional features
if (WITH_EXTENDING_CASADI)
  set_source_files_properties(../extending_casadi/extending_casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_add_module(extending_casadi python ../extending_casadi/extending_casadi.i)
  swig_link_libraries(extending_casadi extending_casadi)
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  endif()
  set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_extending_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

  add_custom_target(extending_casadi_python DEPENDS _extending_casadi extending_casadi)

  install(TARGETS _extending_casadi
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${PROJECT_BINARY_DIR}/swig/python/extending_casadi.py
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../extending_casadi/__init__.py
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

endif()
