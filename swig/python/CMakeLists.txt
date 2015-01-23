cmake_minimum_required(VERSION 2.8.6)

# Python front-end
add_custom_target(python DEPENDS _casadi)

# Workaround for OS X
if(APPLE)
  find_program(PYTHON_CONFIG_EXECUTABLE
               NAMES python-config DOC "python-config executable")
  if(PYTHON_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${PYTHON_CONFIG_EXECUTABLE} --prefix
                    OUTPUT_VARIABLE PYTHON_PREFIX_STRING
                    RESULT_VARIABLE PYTHON_PREFIX_FAILED
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT PYTHON_PREFIX_FAILED)
      find_library(PYTHON_LIBRARY NAMES python2.7
        PATHS ${PYTHON_PREFIX_STRING}/lib NO_DEFAULT_PATH)
      find_path(PYTHON_INCLUDE_DIR
        NAMES Python.h
        PATHS ${PYTHON_PREFIX_STRING}/include/python2.7 NO_DEFAULT_PATH)
    endif()
    message(STATUS ${PYTHON_INCLUDE_DIR})
  endif()
endif()

# Find packages
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(NUMPY REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${NUMPY_INCLUDED_DIRS})

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)

if(WITH_PYTHON_INTERRUPTS)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()

set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DPy_USING_UNICODE")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-noproxydel")

# Generate SWIG wrapper
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)

swig_add_module(casadi python ../casadi.i)
swig_link_libraries(casadi ${PYTHON_LIBRARIES} casadi)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR}")
endif()

# the name used for the shared module for the python interface is determined
get_target_property(CASADI_PYTHON_LIBRARY_SUFFIX ${SWIG_MODULE_casadi_REAL_NAME} SUFFIX)
if(NOT ${CASADI_PYTHON_LIBRARY_SUFFIX})
  if( WIN32 )
    set(CASADI_PYTHON_LIBRARY_SUFFIX ".pyd")
  else()
    set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
  endif()
endif()
set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

# Workaround bug in cmake 2.8.3 and lower. See #704.
add_custom_command(
  TARGET python
  COMMAND ${CMAKE_COMMAND}
  ARGS -E touch ${PROJECT_BINARY_DIR}/_casadi_python_dummy_target
)

if(NOT DEFINED PYTHON_PREFIX)
  if(PYTHONINTERP_FOUND)
    #determining the installation path for python libraries
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/find_installation_path.py"
      OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
      RESULT_VARIABLE PYTHON_INSTALLATION_PATH_RESULT
      ERROR_VARIABLE PYTHON_INSTALLATION_PATH_ERROR)
  endif()
  if(NOT "${PYTHON_INSTALLATION_PATH_RESULT}" STREQUAL "0")
message(FATAL_ERROR "Something went wrong detecting your python setup...
${PROJECT_BINARY_DIR} \$ ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/find_installation_path.py
Returned the following error messages:
${PYTHON_INSTALLATION_PATH_ERROR}
You may try to fix this problem, or supply PYTHON_PREFIX manually.

e.g.
\$ cmake -DPYTHON_PREFIX=\"/usr/lib/python2.7/dist-packages\" ..

The idea is to supply a directory that python knows about (included in PYTHONPATH).

")
  else()
    # Replace backslashes with normal slashes
    if(MSYS)
      string(REGEX REPLACE "\\\\" "/" PYTHON_INSTALLATION_PATH ${PYTHON_INSTALLATION_PATH})
    endif()

    set(PYTHON_PREFIX  ${PYTHON_INSTALLATION_PATH} CACHE PATH "Python installation prefix")
  endif()
endif()

# Construct a dynamic loader for python 
# NOTE(@jaeandersson): Needed?
file(WRITE "${PROJECT_BINARY_DIR}/python/casadi/casadi_loader.py" "failed_modules={}\nloaded_modules=[]\n")
file(APPEND "${PROJECT_BINARY_DIR}/python/casadi/casadi_loader.py"
"try:\n"
"  from casadi import *\n"
"  loaded_modules.append('casadi')\n"
"except Exception as e:\n"
"  failed_modules['casadi'] = str(e)\n"
)

# NOTE(@jaeandersson): What is this?
install(DIRECTORY ${PROJECT_BINARY_DIR}/python/casadi
  DESTINATION "${PYTHON_PREFIX}"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
  PATTERN .py EXCLUDE
)

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

# Install Python proxy classes
install(FILES ${PROJECT_BINARY_DIR}/swig/python/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

# Install Python tools
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
)

# Install Python package initialization
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

# For DLL-based systems, rpath is not available,
# so copy the dlls to the same directory as the .pyd files
file(GLOB dlls "${PROJECT_BINARY_DIR}/lib/*.dll")
install(FILES ${dlls}
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)
