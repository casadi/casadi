cmake_minimum_required(VERSION 2.6)

# NOTE: you are not allowd to generate header files with configure_file in this module.
# See https://github.com/casadi/casadi/issues/646

# Enforce cmake2.4+ style behaviour for mixed filename/full path linker directives
if(COMMAND cmake_policy)
cmake_policy(SET CMP0003 NEW)
endif()

# cmake_policy propagation rules
if(COMMAND cmake_policy)
if(POLICY CMP0011)
cmake_policy(SET CMP0011 OLD)
endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

# Find packages
find_package(PythonLibs)
find_package(PythonInterp)
find_package(NUMPY)
if(WITH_MATLAB)
  find_package(MATLAB)
endif()
set(CASADI_SWIG_FLAGS "")

# Include the use swig file
include(${CMAKE_MODULE_PATH}/UseSWIGMOD.cmake)

if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
  if(WITH_IPOPT_CALLBACK)
    set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT_CALLBACK")
  endif()
endif()

if(LAPACK_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_LAPACK")
endif()

# enabling openmp support if requesed
if(OPENMP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_OPENMP")
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../doc/doc.i)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_DOXDOC")
endif()

# Find python library and add include path for python headers
if(PYTHONLIBS_FOUND)
  include_directories(${PYTHON_INCLUDE_PATH})
endif()

if(NUMPY_FOUND)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${NUMPY_INCLUDED_DIRS})
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_NUMPY")
endif()

if(WITH_PYTHON_INTERRUPTS)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-I${CMAKE_CURRENT_SOURCE_DIR}/swig_lib")
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DCASADI_MODULE")

set(CASADI_SWIG_MAIN_MODULE casadi_core)

set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DPy_USING_UNICODE")
set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-noproxydel")

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=casadi")
foreach(my_module ${CASADI_SWIG_MODULES})
set_source_files_properties(${my_module}.i  PROPERTIES  CPLUSPLUS ON)
endforeach()
set_source_files_properties(${CASADI_SWIG_MAIN_MODULE}.i  PROPERTIES SWIG_FLAGS "-DSWIG_MAIN_MODULE")

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
endif()

# PYTHON interface
if(WITH_WERROR)
  set(MAYBE_WERROR "-Werror")
else()
  set(MAYBE_WERROR "")
endif()

foreach(my_module ${CASADI_SWIG_MODULES})
swig_add_module(${my_module} python ${my_module}.i)
swig_link_libraries(${my_module} ${PYTHON_LIBRARIES} ${my_module})
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set_target_properties(_${my_module} PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(_${my_module} PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR}")
endif()
endforeach()

# the name used for the shared module for the python interface is determined
get_target_property(CASADI_PYTHON_LIBRARY_SUFFIX ${SWIG_MODULE_${CASADI_SWIG_MAIN_MODULE}_REAL_NAME} SUFFIX)

if(NOT ${CASADI_PYTHON_LIBRARY_SUFFIX})
  if( WIN32 )
    set(CASADI_PYTHON_LIBRARY_SUFFIX ".pyd")
  else()
    set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
  endif()
endif()

foreach(my_module ${CASADI_SWIG_MODULES})
  set(CASADI_PYTHON_LIBRARY_${my_module} ${SWIG_MODULE_${my_module}_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})
endforeach()

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)
# adding some extra commands to the target _casadi which is declared automatically when constructing a swig-python interface

set(modules_python)
foreach(my_module ${CASADI_SWIG_MODULES})
  list(APPEND modules_python _${my_module})
endforeach()
add_custom_target(_python_modules DEPENDS ${modules_python}) # the python interface is created with the command "make python"

if(WITH_MATLAB)
  set(modules_matlab) 
  include_directories(${MATLAB_INCLUDE_DIR})
  add_definitions(${MATLAB_EXTRA_CXXFLAGS})
  foreach(my_module ${CASADI_SWIG_MAIN_MODULE})
    list(APPEND modules_matlab ${my_module}_matlab)
  endforeach()
  add_custom_target(_matlab_modules DEPENDS ${modules_matlab}) # the matlab interface is created with the command "make matlab"
endif()

set(modules_json)
foreach(my_module ${CASADI_SWIG_MODULES})
  list(APPEND modules_json ${my_module}_json)
endforeach()
add_custom_target(_json_modules DEPENDS ${modules_json}) # the json interface is created with the command "make json"

# Workaround bug in cmake 2.8.3 and lower. See #704.
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E touch ${PROJECT_BINARY_DIR}/_python_modules
)

# workaround for Visual Studio! The default way is not general enough and should be changed!
if(MSVC)
  foreach(my_module ${CASADI_SWIG_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/Release/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
else()
  foreach(my_module ${CASADI_SWIG_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
endif()

foreach(my_module ${CASADI_SWIG_MODULES})
add_custom_command(
  TARGET _${my_module}
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PYTHONLIBPATH_${my_module}} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY_${my_module}}
)
endforeach()

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
foreach(my_module ${CASADI_SWIG_MODULES})
file(APPEND "${PROJECT_BINARY_DIR}/python/casadi/casadi_loader.py"
"try:\n"
"  from ${my_module} import *\n"
"  loaded_modules.append('${my_module}')\n"
"except Exception as e:\n"
"  failed_modules['${my_module}'] = str(e)\n"
)
endforeach()

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

foreach(my_module ${modules_python})
  install(TARGETS ${my_module}
    DESTINATION "${PYTHON_PREFIX}/casadi"
    COMPONENT install_python
  )
endforeach()

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
foreach(my_module ${CASADI_SWIG_MODULES})
  install(FILES ${PROJECT_BINARY_DIR}/swig/${my_module}.py
    DESTINATION "${PYTHON_PREFIX}/casadi"
    COMPONENT install_python
  )
endforeach()

#creation of a custom installation command
add_custom_target(install_python
COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install_python
  -D CMAKE_INSTALL_PREFIX="${PYTHON_PREFIX}"
  -P cmake_install.cmake
)
foreach(my_module ${CASADI_SWIG_MODULES})
add_dependencies(install_python _${my_module})
endforeach()

# MATLAB interface
if(WITH_MATLAB)
  foreach(my_module ${CASADI_SWIG_MAIN_MODULE})
    swig_add_module(${my_module}_matlab matlab ${my_module}.i)
  endforeach()
endif()
  
# json interface

# fooXML_wrap.xml are a big files which Joris wants to delete immedately after using.
# These files are dependencies of json outputs, but they don't exist according to cmake.
# So we silence these warnings. THIS IS A HACK AND WE NEED A BETTER SOLUTION
if(COMMAND cmake_policy)
  if(CMAKE_MAJOR_VERSION GREATER 2)
    cmake_policy(SET CMP0046 OLD)
  endif()
endif()

foreach(my_module ${CASADI_SWIG_MODULES})
  swig_add_module_xml(${my_module}_xml ${my_module}.i)
endforeach()

foreach(my_module ${CASADI_SWIG_MODULES})
add_custom_target(${my_module}_json
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/json/generate.py ${PROJECT_BINARY_DIR}/swig/${my_module}XML_wrap.xml ${my_module}
  COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/swig/${my_module}XML_wrap.xml
)
add_dependencies(${my_module}_json ${my_module}_xml
                 ${CMAKE_CURRENT_SOURCE_DIR}/json/generate.py
                 ${PROJECT_BINARY_DIR}/swig/${my_module}XML_wrap.xml)
endforeach()

#foreach(my_module ${CASADI_SWIG_MODULES})
#install(FILES ${PROJECT_BINARY_DIR}/swig/${my_module}.json
#  DESTINATION share/casadi
#  COMPONENT componentjson
#)
#endforeach()


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
