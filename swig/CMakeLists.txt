cmake_minimum_required(VERSION 2.6)

find_package(SWIG REQUIRED)

# NOTE: you are not allowd to generate header files with configure_file in this module.
# See https://github.com/casadi/casadi/issues/646

# Enforce cmake2.4+ style behaviour for mixed filename/full path linker directives
cmake_policy(SET CMP0003 NEW)

# cmake_policy propagation rules
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

set(CASADI_SWIG_FLAGS "")

# Include the use swig file
include(${CMAKE_MODULE_PATH}/UseSWIGMOD.cmake)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../doc/doc.i)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_DOXDOC")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=casadi")

set_source_files_properties(casadi.i  PROPERTIES  CPLUSPLUS ON)

set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "${CASADI_SWIG_FLAGS}")
set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-fvirtual")

if(WITH_WERROR)
  # This will make the swig command return an error so that "make" will fail,
  # but swig will emit the output files anyways so subsequent calls to "make"
  # will succeed. At least buildbot will report failures.
  #
  # I'm pretty sure this is also OK also in windows, see:
  # http://www.swig.org/Doc3.0/Warnings.html#Warnings_nn7
  set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "-Werror")
  set(MAYBE_WERROR "-Werror")
else()
  set(MAYBE_WERROR "")
endif()

if(WITH_PYTHON)

# Python front-end
#add_subdirectory(python)

# Find packages
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH})

find_package(NUMPY)
if(NUMPY_FOUND)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${NUMPY_INCLUDED_DIRS})
  set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DWITH_NUMPY")
endif()

if(WITH_PYTHON_INTERRUPTS)
  set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DWITH_PYTHON_INTERRUPTS")
endif()

set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DPy_USING_UNICODE")
set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-noproxydel")

swig_add_module(casadi python casadi.i)
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

set(CASADI_PYTHON_LIBRARY_casadi ${SWIG_MODULE_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)
# adding some extra commands to the target _casadi which is declared automatically when constructing a swig-python interface

add_custom_target(_python_modules DEPENDS _casadi)

# Workaround bug in cmake 2.8.3 and lower. See #704.
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E touch ${PROJECT_BINARY_DIR}/_python_modules
)

# workaround for Visual Studio! The default way is not general enough and should be changed!
if(MSVC)
  set(PYTHONLIBPATH_casadi ${LIBRARY_OUTPUT_PATH}/Release/${CASADI_PYTHON_LIBRARY_casadi})
else()
  set(PYTHONLIBPATH_casadi ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY_casadi})
endif()

# add_custom_command(
#   TARGET _casadi
#   COMMAND ${CMAKE_COMMAND}
#   ARGS -E copy ${PYTHONLIBPATH_casadi} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY_casadi}
# )

if(NOT DEFINED PYTHON_PREFIX)
  if(PYTHONINTERP_FOUND)
    #determining the installation path for python libraries
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/python/find_installation_path.py"
      OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
      RESULT_VARIABLE PYTHON_INSTALLATION_PATH_RESULT
      ERROR_VARIABLE PYTHON_INSTALLATION_PATH_ERROR)
  endif()
  if(NOT "${PYTHON_INSTALLATION_PATH_RESULT}" STREQUAL "0")
message(FATAL_ERROR "Something went wrong detecting your python setup...
${PROJECT_BINARY_DIR} \$ ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/find_installation_path.py
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
file(WRITE "${PROJECT_BINARY_DIR}/python/casadi/casadi_loader.py" "failed_modules={}\nloaded_modules=[]\n")
file(APPEND "${PROJECT_BINARY_DIR}/python/casadi/casadi_loader.py"
"try:\n"
"  from casadi import *\n"
"  loaded_modules.append('casadi')\n"
"except Exception as e:\n"
"  failed_modules['casadi'] = str(e)\n"
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/python/casadi
  DESTINATION "${PYTHON_PREFIX}"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
  PATTERN .py EXCLUDE
)

#  For DLL-based systems, rpath is not available,
# so copy the dlls to the same directory as the .pyd files
file(GLOB dlls "${PROJECT_BINARY_DIR}/lib/*.dll")
install(FILES ${dlls}
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

install(TARGETS _casadi
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/tools
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

install(FILES ${PROJECT_BINARY_DIR}/swig/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

#creation of a custom installation command
add_custom_target(install_python
COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install_python
  -D CMAKE_INSTALL_PREFIX="${PYTHON_PREFIX}"
  -P cmake_install.cmake
)

add_dependencies(install_python _casadi)

endif(WITH_PYTHON)

# MATLAB front-end
if(WITH_MATLAB)
  #add_subdirectory(matlab)

  find_package(MATLAB REQUIRED)

  include_directories(${MATLAB_INCLUDE_DIR})
  add_definitions(${MATLAB_EXTRA_CXXFLAGS})
  add_custom_target(_matlab_modules DEPENDS casadi_matlab)

  swig_add_module(casadi_matlab matlab casadi.i)

endif(WITH_MATLAB)

# JSON/Haskell front-end
if(WITH_JSON)
#add_subdirectory(json)

add_custom_target(_json_modules DEPENDS casadi_json)

# fooXML_wrap.xml are a big files which Joris wants to delete immedately after using.
# These files are dependencies of json outputs, but they don't exist according to cmake.
# So we silence these warnings. THIS IS A HACK AND WE NEED A BETTER SOLUTION
if(POLICY CMP0046)
  cmake_policy(SET CMP0046 OLD)
endif()

swig_add_module_xml(casadi_xml casadi.i)
add_custom_target(casadi_json
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/json/generate.py ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml casadi
  COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml
)
add_dependencies(casadi_json casadi_xml
                 ${CMAKE_CURRENT_SOURCE_DIR}/json/generate.py
                 ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml)

#install(FILES ${PROJECT_BINARY_DIR}/swig/casadi.json
#  DESTINATION share/casadi
#  COMPONENT componentjson
#)

#set(CPACK_MONOLITHIC_INSTALL FALSE)
set(CPACK_TGZ_COMPONENT_INSTALL TRUE)

set(CPACK_COMPONENT_COMPONENTJSON_DISPLAY_NAME "Json meta-data")
set(CPACK_COMPONENT_COMPONENTJSON_DESCRIPTION
  "Installs meta-data about CasADi in json form.")

set(CPACK_COMPONENTS_GROUPING "ALL_COMPONENTS_IN_ONE")
set(CPACK_COMPONENTS_ALL Unspecified componentjson)
# This will not actually work since we are in a EXCLUDE_FROM_ALL
# http://cmake.3232098.n2.nabble.com/Yet-another-install-vs-exclude-from-all-thread-td4967491.html

# Pass meta information
configure_file(
  ${CMAKE_SOURCE_DIR}/swig/linklist.txt.cmake
  ${PROJECT_BINARY_DIR}/swig/linklist.txt
  ESCAPE_QUOTES
)

include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_CONTACT "${PACKAGE_MAINTAINER} <contact@casadi.org>")
set(CPACK_PACKAGE_NAME "libcasadi")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_DESCRIPTION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt")
set(CPACK_DEBIAN_PACKAGE_LICENSE "LGPL")

set(THREE_PART_VERSION_REGEX "[0-9]+\\.[0-9]+\\.[0-9]+")

# Breaks up a string in the form n1.n2.n3 into three parts and stores
# them in major, minor, and patch.  version should be a value, not a
# variable, while major, minor and patch should be variables.
macro(THREE_PART_VERSION_TO_VARS version major minor patch)
  if(${version} MATCHES ${THREE_PART_VERSION_REGEX})
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+\\+?" "\\1" ${major} "${version}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+\\+?" "\\1" ${minor} "${version}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)\\+?" "\\1" ${patch} "${version}")
  else()
    message("macro(THREE_PART_VERSION_TO_VARS ${version} ${major} ${minor} ${patch}")
    message(FATAL_ERROR "Problem parsing version string, I can't parse it properly.")
  endif()
endmacro()

THREE_PART_VERSION_TO_VARS(${PACKAGE_VERSION} CPACK_PACKAGE_VERSION_MAJOR CPACK_PACKAGE_VERSION_MINOR CPACK_PACKAGE_VERSION_PATCH)

if(ENABLE_STATIC)
  set(CPACK_PACKAGE_FILE_NAME "libcasadi-static-${PACKAGE_VERSION_FULL}-${CMAKE_SYSTEM_PROCESSOR}")
else()
  set(CPACK_PACKAGE_FILE_NAME "libcasadi-shared-${PACKAGE_VERSION_FULL}-${CMAKE_SYSTEM_PROCESSOR}")
endif()

# if version end on a plus, it is a late build
if(${PACKAGE_VERSION} MATCHES "\\+$")
set(CPACK_PACKAGE_VERSION_PATCH "99")
endif()

include(CPack)

endif(WITH_JSON)
