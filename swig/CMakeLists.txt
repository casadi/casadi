cmake_minimum_required (VERSION 2.6)

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

find_package(PythonLibs)
find_package(NUMPY)



set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules)

set(CASADI_SWIG_FLAGS "")
set(SWIG_LINK_LIBRARIES_LIST "")

# Include the use swig file
include (../cmake_modules/UseSWIGMOD.cmake)

if(NOT WITH_DEPRECATED)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITHOUT_PRE_1_9_X")
endif()

set(CASADI_SWIG_INTERFACES)

if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT")
  if(WITH_IPOPT_CALLBACK)
    set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT_CALLBACK")
  endif()
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} ipopt)
endif()

if(WITH_SUNDIALS)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SUNDIALS")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} sundials)
endif()

if(QPOASES_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_QPOASES")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} qpoases)
endif()

if(DSDP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_DSDP")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} dsdp)
endif()

if(WITH_CSPARSE)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CSPARSE")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} csparse)
endif()

if(KNITRO_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_KNITRO")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} knitro)
endif()

if(CPLEX_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CPLEX")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} cplex)
endif()

# if(GSL_FOUND)
#   set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_GSL")
#   set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_gsl_interface ${GSL_LIBRARIES})
# endif()

if(OOQP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_OOQP")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} ooqp)
endif()

if(SQIC_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SQIC")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} sqic)
endif()

if(SLICOT_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SLICOT")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} slicot)
endif()

if(WORHP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_WORHP")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} worhp)
endif()

if(SNOPT_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SNOPT")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} snopt)
endif()

if(LAPACK_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_LAPACK")
  set(CASADI_SWIG_INTERFACES ${CASADI_SWIG_INTERFACES} lapack)
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
  include_directories (${PYTHON_INCLUDE_PATH})
endif()

if(NUMPY_FOUND)
  include_directories (${CMAKE_CURRENT_SOURCE_DIR})
  include_directories (${NUMPY_INCLUDED_DIRS})
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_NUMPY")
endif()

if(WITH_PYTHON_INTERRUPTS)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-I${CMAKE_CURRENT_SOURCE_DIR}/swig_lib")
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DCASADI_MODULE")
  
# set the global swig flags to empty
# set(CMAKE_SWIG_FLAGS "-outdir ${PROJECT_BINARY_DIR}")
# SEPARATE_ARGUMENTS(CMAKE_SWIG_FLAGS)

# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
# set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# SEPARATE_ARGUMENTS(CASADI_SWIG_FLAGS)

set(CASADI_SWIG_MAIN_MODULE casadi_main)
set(CASADI_SWIG_SEPARATE_MODULES casadi_main casadi_primitive casadi_primitive_tools casadi_noncore)
foreach(my_interface ${CASADI_SWIG_INTERFACES})
set(CASADI_SWIG_SEPARATE_MODULES ${CASADI_SWIG_SEPARATE_MODULES} casadi_interface_${my_interface})
endforeach()

#set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-noproxydel") # This should be active during python, and nonactive during octave. How?

set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DPy_USING_UNICODE")
set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-noproxydel")

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=${CASADI_SWIG_MAIN_MODULE}")
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
set_source_files_properties (${my_module}.i  PROPERTIES  CPLUSPLUS ON)
endforeach()
set_source_files_properties (${CASADI_SWIG_MAIN_MODULE}.i  PROPERTIES SWIG_FLAGS "-DSWIG_MAIN_MODULE")

set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "${CASADI_SWIG_FLAGS}")

# PYTHON interface
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
SWIG_ADD_MODULE (${my_module} python ${my_module}.i)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(_${my_module} PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-delete-non-virtual-dtor -Wno-self-assign")
endif()
endforeach()

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  SWIG_LINK_LIBRARIES (${my_module} ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control casadi_control casadi_integration casadi_convex_programming  casadi_nonlinear_programming)
endforeach()
foreach(my_interface ${CASADI_SWIG_INTERFACES})
  SWIG_LINK_LIBRARIES (casadi_interface_${my_interface} casadi_${my_interface}_interface)
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

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  set(CASADI_PYTHON_LIBRARY_${my_module} ${SWIG_MODULE_${my_module}_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})
endforeach()

# a python library is built in the build directory inside swig/python 
make_directory(${PROJECT_BINARY_DIR}/python/casadi)
# adding some extra commands to the target _casadi which is declared automatically when constructing a swig-python interface

set(modules_python)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  list(APPEND modules_python _${my_module})
endforeach()
add_custom_target(_python_modules DEPENDS ${modules_python}) # the python interface is created with the command "make python"
set(modules_octave)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  list(APPEND modules_octave ${my_module}_octave)
endforeach()
add_custom_target(_octave_modules DEPENDS ${modules_octave}) # the octave interface is created with the command "make octave"
set(modules_haskell)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  list(APPEND modules_haskell ${my_module}_tree.hs)
endforeach()
add_custom_target(_haskell_modules DEPENDS ${modules_haskell}) # the haskell interface is created with the command "make octave"

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_BINARY_DIR}/swig/${my_module}.py ${PROJECT_BINARY_DIR}/python/casadi/${my_module}.py
)
endforeach()
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/casadi.py ${PROJECT_BINARY_DIR}/python/casadi/casadi.py
)

# Workaround bug in cmake 2.8.3 and lower. See #704.
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E touch ${PROJECT_BINARY_DIR}/_python_modules
)
add_custom_command(
  TARGET _octave_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E touch ${PROJECT_BINARY_DIR}/_octave_modules
)

# workaround for Visual Studio! The default way is not general enough and should be changed!
if(MSVC)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/Release/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
else(MSVC)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
endif(MSVC)

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
add_custom_command(
  TARGET _${my_module}
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PYTHONLIBPATH_${my_module}} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY_${my_module}}
)
endforeach()
  
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)

  if(NOT DEFINED PYTHON_PREFIX )
    #determining the installation path for python libraries
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/python/find_installation_path.py"
      OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
      RESULT_VARIABLE PYTHON_INSTALLATION_PATH_RESULT
      ERROR_VARIABLE PYTHON_INSTALLATION_PATH_ERROR
    )
    
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

  #since in the main CMakeLists.txt we used "add_subdirectory(swig EXCLUDE_FROM_ALL)" the following install command won't add anything to "make install"
  install(DIRECTORY ${PROJECT_BINARY_DIR}/python/casadi 
    DESTINATION ${PYTHON_PREFIX}
    COMPONENT install_python
    USE_SOURCE_PERMISSIONS
    PATTERN .pyc EXCLUDE
    PATTERN .svn EXCLUDE
    PATTERN .py EXCLUDE
  )
  
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/tools
    DESTINATION ${PYTHON_PREFIX}/casadi
    COMPONENT install_python
    USE_SOURCE_PERMISSIONS
    PATTERN .pyc EXCLUDE
    PATTERN .svn EXCLUDE
  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
    DESTINATION ${PYTHON_PREFIX}/casadi
    COMPONENT install_python
  )
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/casadi.py
    DESTINATION ${PYTHON_PREFIX}/casadi
    COMPONENT install_python
  )
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
    install(FILES ${PROJECT_BINARY_DIR}/swig/${my_module}.py 
      DESTINATION ${PYTHON_PREFIX}/casadi
      COMPONENT install_python
    )
  endforeach()
endif()

  #creation of a custom installation command
add_custom_target(install_python
COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install_python
  -P cmake_install.cmake
)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
add_dependencies(install_python DEPENDS _${my_module})
endforeach()


# OCTAVE interface
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  SWIG_ADD_MODULE (${my_module}_octave octave ${my_module}.i octrun.swg)
  SWIG_LINK_LIBRARIES (${my_module}_octave ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control casadi_control  casadi_integration casadi_convex_programming casadi_nonlinear_programming)
  set_target_properties(${my_module}_octave PROPERTIES PREFIX "")
  set_target_properties(${my_module}_octave PROPERTIES SUFFIX ".oct")
  set_target_properties(${my_module}_octave PROPERTIES OUTPUT_NAME "${my_module}")
endforeach()
foreach(my_interface ${CASADI_SWIG_INTERFACES})
  SWIG_LINK_LIBRARIES (casadi_interface_${my_interface}_octave casadi_${my_interface}_interface)
endforeach()

# haskell interface
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  SWIG_ADD_MODULE_XML (${my_module}_haskell ${my_module}.i)
endforeach()

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
ADD_CUSTOM_COMMAND(
   OUTPUT ${my_module}_tree.hs 
   COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/haskell/generate.py ${PROJECT_BINARY_DIR}/swig/${my_module}XML_wrap.xml ${my_module}
   DEPENDS ${PROJECT_BINARY_DIR}/swig/${my_module}XML_wrap.xml ${CMAKE_CURRENT_SOURCE_DIR}/haskell/generate.py
)
endforeach()

foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
install(FILES ${PROJECT_BINARY_DIR}/swig/${my_module}_tree.hs ${PROJECT_BINARY_DIR}/swig/${my_module}_classes.hs ${PROJECT_BINARY_DIR}/swig/${my_module}_swiginclude.hpp
  DESTINATION share/casadi
  COMPONENT componenthaskell
)
endforeach()


#SET(CPACK_MONOLITHIC_INSTALL FALSE)
SET(CPACK_TGZ_COMPONENT_INSTALL TRUE)

set(CPACK_COMPONENT_COMPONENTHASKELL_DISPLAY_NAME "Haskell headers")
set(CPACK_COMPONENT_COMPONENTHASKELL_DESCRIPTION 
   "Necessary headers to install a Haskell interface to CasADi")
   
set(CPACK_COMPONENTS_GROUPING "ALL_COMPONENTS_IN_ONE")
set(CPACK_COMPONENTS_ALL Unspecified componenthaskell)
# This will not actually work since we are in a EXCLUDE_FROM_ALL
# http://cmake.3232098.n2.nabble.com/Yet-another-install-vs-exclude-from-all-thread-td4967491.html

set(LINK_LIBRARIES_ALL ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control casadi_control casadi_integration casadi_convex_programming casadi_nonlinear_programming)

# Pass meta information
CONFIGURE_FILE(
  ${CMAKE_SOURCE_DIR}/swig/linklist.txt.cmake
  ${PROJECT_BINARY_DIR}/swig/linklist.txt
  ESCAPE_QUOTES
)

# begin Patch swig to get concat functionality
find_package(OCTAVE)
if(OCTAVE_FOUND)
include_directories (${OCTAVE_INCLUDE_DIR})
EXECUTE_PROCESS(COMMAND ${SWIG_EXECUTABLE} -swiglib OUTPUT_VARIABLE ORIG_SWIGDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

ADD_CUSTOM_COMMAND(
   OUTPUT octrun.swg
   COMMAND ${CMAKE_COMMAND} -E remove octrunORIG.swg octrun.swg
   COMMAND ${SWIG_EXECUTABLE} -octave -o octrunORIG.swg -co octrun.swg
   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/octave/swig_lib .
   COMMAND ${OCTAVE_EXECUTABLE} -q octrun_patch.m
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/octave/swig_lib/octrun_patch.m ${ORIG_SWIGDIR}/octave/octrun.swg
)
else(OCTAVE_FOUND)
#Don't do anything
ADD_CUSTOM_COMMAND(
   OUTPUT octrun.swg
   COMMENT "No octave patch will be applied (octave executable not found)"
)
endif(OCTAVE_FOUND)
# end patch

ADD_CUSTOM_COMMAND(
 TARGET ${CASADI_SWIG_MAIN_MODULE}_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/octave/@swig_ref ${LIBRARY_OUTPUT_PATH}/@swig_ref
)
ADD_CUSTOM_COMMAND(
 TARGET ${CASADI_SWIG_MAIN_MODULE}_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E remove ${LIBRARY_OUTPUT_PATH}/casadi.oct
)
ADD_CUSTOM_COMMAND(
 TARGET ${CASADI_SWIG_MAIN_MODULE}_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/octave/casadi_multi.m ${LIBRARY_OUTPUT_PATH}/casadi.m
)
ADD_CUSTOM_COMMAND(
 TARGET ${CASADI_SWIG_MAIN_MODULE}_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/octave/casadi_helpers.m ${LIBRARY_OUTPUT_PATH}/casadi_helpers.m
)


include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_CONTACT "${PACKAGE_MAINTAINER} <contact@casadi.org>")
set(CPACK_PACKAGE_NAME "libcasadi")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_DESCRIPTION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt")
SET(CPACK_DEBIAN_PACKAGE_LICENSE "LGPL")

set(THREE_PART_VERSION_REGEX "[0-9]+\\.[0-9]+\\.[0-9]+")

# Breaks up a string in the form n1.n2.n3 into three parts and stores
# them in major, minor, and patch.  version should be a value, not a
# variable, while major, minor and patch should be variables.
MACRO(THREE_PART_VERSION_TO_VARS version major minor patch)
  IF(${version} MATCHES ${THREE_PART_VERSION_REGEX})
    STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+\\+?" "\\1" ${major} "${version}")
    STRING(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+\\+?" "\\1" ${minor} "${version}")
    STRING(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)\\+?" "\\1" ${patch} "${version}")
  ELSE(${version} MATCHES ${THREE_PART_VERSION_REGEX})
    MESSAGE("MACRO(THREE_PART_VERSION_TO_VARS ${version} ${major} ${minor} ${patch}")
    MESSAGE(FATAL_ERROR "Problem parsing version string, I can't parse it properly.")
  ENDIF(${version} MATCHES ${THREE_PART_VERSION_REGEX})
ENDMACRO(THREE_PART_VERSION_TO_VARS)

THREE_PART_VERSION_TO_VARS(${PACKAGE_VERSION} CPACK_PACKAGE_VERSION_MAJOR CPACK_PACKAGE_VERSION_MINOR CPACK_PACKAGE_VERSION_PATCH)

set(PACKAGE_VERSION_FULL "${git_describe}")
if ((PACKAGE_VERSION_FULL STREQUAL "") OR (NOT PACKAGE_VERSION MATCHES "\\+"))
  set(PACKAGE_VERSION_FULL "${PACKAGE_VERSION}")
endif()

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


