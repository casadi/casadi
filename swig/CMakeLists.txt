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

if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT")
  if(WITH_IPOPT_CALLBACK)
    set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT_CALLBACK")
  endif()
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_ipopt_interface ${IPOPT_LIBRARIES})
endif()

if(WITH_SUNDIALS)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SUNDIALS")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_sundials_interface ${SUNDIALS_LIBRARIES})
endif()

if(ACADO_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_ACADO")
  if(ACADO_SHARED_LIBRARIES)
    set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_acado_interface casadi_qpoases_interface ${ACADO_SHARED_LIBRARIES})
  else(ACADO_SHARED_LIBRARIES)
    if(ACADO_LIBRARIES)
      set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_acado_interface casadi_qpoases_interface ${ACADO_LIBRARIES})
    endif(ACADO_LIBRARIES)
  endif(ACADO_SHARED_LIBRARIES)
endif()

if(QPOASES_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_QPOASES")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_qpoases_interface ${QPOASES_LIBRARIES})
endif()

if(DSDP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_DSDP")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_dsdp_interface ${DSDP_LIBRARIES})
endif()

if(WITH_CSPARSE)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CSPARSE")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_csparse_interface ${CSPARSE_LIBRARIES})
endif()

if(KNITRO_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_KNITRO")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_knitro_interface ${KNITRO_LIBRARIES})
endif()

if(CPLEX_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CPLEX")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_cplex_interface ${CPLEX_LIBRARIES})
endif()

# if(GSL_FOUND)
#   set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_GSL")
#   set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_gsl_interface ${GSL_LIBRARIES})
# endif()

if(OOQP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_OOQP")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_ooqp_interface ${OOQP_LIBRARIES})
endif()

if(WORHP_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_WORHP")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_worhp_interface ${WORHP_LIBRARIES})
endif()

if(LAPACK_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_LAPACK")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} casadi_lapack_interface ${LAPACK_LIBRARIES})
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
  
# set the global swig flags to empty
# set(CMAKE_SWIG_FLAGS "-outdir ${PROJECT_BINARY_DIR}")
# SEPARATE_ARGUMENTS(CMAKE_SWIG_FLAGS)

# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
# set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# SEPARATE_ARGUMENTS(CASADI_SWIG_FLAGS)

set(CASADI_SWIG_MAIN_MODULE casadi_main)
set(CASADI_SWIG_SEPARATE_MODULES casadi_main casadi_primitive casadi_primitive_tools casadi_noncore)

#set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-noproxydel") # This should be active during python, and nonactive during octave. How?

set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-DPy_USING_UNICODE")
set(USE_SWIG_FLAGS_python ${USE_SWIG_FLAGS_python} "-noproxydel")

if(WITH_SWIG_SPLIT)
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SWIG_SPLIT")
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DSWIG_TYPE_TABLE=${CASADI_SWIG_MAIN_MODULE}")
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
set_source_files_properties (${my_module}.i  PROPERTIES  CPLUSPLUS ON)
endforeach()
set_source_files_properties (${CASADI_SWIG_MAIN_MODULE}.i  PROPERTIES SWIG_FLAGS "-DSWIG_MAIN_MODULE")
else(WITH_SWIG_SPLIT)
set_source_files_properties (casadi.i PROPERTIES
  SWIG_FLAGS "-DSWIG_MAIN_MODULE"
  CPLUSPLUS ON
)
endif(WITH_SWIG_SPLIT)

set(USE_SWIG_FLAGS ${USE_SWIG_FLAGS} "${CASADI_SWIG_FLAGS}")

# PYTHON interface
if(WITH_SWIG_SPLIT)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
SWIG_ADD_MODULE (${my_module} python ${my_module}.i)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(_${my_module} PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-delete-non-virtual-dtor -Wno-self-assign")
endif()
endforeach()
else(WITH_SWIG_SPLIT)
SWIG_ADD_MODULE (casadi python casadi.i) #the target corresponding this module is called "_casadi"
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-delete-non-virtual-dtor -Wno-self-assign")
endif()
endif(WITH_SWIG_SPLIT)

if(WITH_SWIG_SPLIT)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
    SWIG_LINK_LIBRARIES (${my_module} ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control ${TINYXML_LIBRARIES} casadi_integration casadi_nonlinear_programming casadi ${CASADI_DEPENDENCIES})
  endforeach()
  else()
    SWIG_LINK_LIBRARIES (casadi ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control ${TINYXML_LIBRARIES} casadi_integration casadi_nonlinear_programming casadi ${CASADI_DEPENDENCIES})
endif()

# the name used for the shared module for the python interface is determined
if(WITH_SWIG_SPLIT)
  get_target_property(CASADI_PYTHON_LIBRARY_SUFFIX ${SWIG_MODULE_${CASADI_SWIG_MAIN_MODULE}_REAL_NAME} SUFFIX)
else()
  get_target_property(CASADI_PYTHON_LIBRARY_SUFFIX ${SWIG_MODULE_casadi_REAL_NAME} SUFFIX)
endif()

if(NOT ${CASADI_PYTHON_LIBRARY_SUFFIX})
  if( WIN32 )
    set(CASADI_PYTHON_LIBRARY_SUFFIX ".pyd")
  else()
    set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
  endif()
endif()

if(WITH_SWIG_SPLIT)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
    set(CASADI_PYTHON_LIBRARY_${my_module} ${SWIG_MODULE_${my_module}_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})
  endforeach()
else()
  set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})
endif()

# a python library is built in the build directory inside swig/python 
make_directory(${PROJECT_BINARY_DIR}/python/casadi)
# adding some extra commands to the target _casadi which is declared automatically when constructing a swig-python interface

if(WITH_SWIG_SPLIT)
  add_custom_target(_python_modules DEPENDS _casadi_main _casadi_primitive _casadi_primitive_tools _casadi_noncore) # the python interface is created with the command "make python"
  add_custom_target(_octave_modules DEPENDS casadi_main_octave casadi_primitive_octave casadi_primitive_tools_octave casadi_noncore_octave) # the octave interface is created with the command "make octave"
else(WITH_SWIG_SPLIT)
  add_custom_target(_python_modules DEPENDS _casadi ) # the python interface is created with the command "make python"
  add_custom_target(_octave_modules DEPENDS casadi_octave ) # the octave interface is created with the command "make octave"
endif(WITH_SWIG_SPLIT)
  
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py ${PROJECT_BINARY_DIR}/python/casadi/__init__.py
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/python/tools ${PROJECT_BINARY_DIR}/python/casadi/tools
)

if(WITH_SWIG_SPLIT)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_BINARY_DIR}/${my_module}.py ${PROJECT_BINARY_DIR}/python/casadi/${my_module}.py
)
endforeach()
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/casadi.py ${PROJECT_BINARY_DIR}/python/casadi/casadi.py
)
else(WITH_SWIG_SPLIT)
add_custom_command(
  TARGET _python_modules
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_BINARY_DIR}/swig/casadi.py ${PROJECT_BINARY_DIR}/python/casadi/casadi.py
) 
endif(WITH_SWIG_SPLIT)

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
  if(WITH_SWIG_SPLIT)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/Release/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
  else(WITH_SWIG_SPLIT)
  set(PYTHONLIBPATH ${LIBRARY_OUTPUT_PATH}/Release/${CASADI_PYTHON_LIBRARY})
  endif(WITH_SWIG_SPLIT)
else(MSVC)
  if(WITH_SWIG_SPLIT)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  set(PYTHONLIBPATH_${my_module} ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY_${my_module}})
  endforeach()
  else(WITH_SWIG_SPLIT)
  set(PYTHONLIBPATH ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY})
  endif(WITH_SWIG_SPLIT)
endif(MSVC)

if(WITH_SWIG_SPLIT)
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  add_custom_command(
    TARGET _${my_module}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy ${PYTHONLIBPATH_${my_module}} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY_${my_module}}
  )
  endforeach()
else(WITH_SWIG_SPLIT)
add_custom_command(
  TARGET _casadi
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PYTHONLIBPATH} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY}
)
endif(WITH_SWIG_SPLIT)
  
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)

  #determining the installation path for python libraries
  EXEC_PROGRAM( "${PYTHON_EXECUTABLE}"
    ARGS "${CMAKE_CURRENT_SOURCE_DIR}/python/find_installation_path.py"
    OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
  )

  # Replace backslashes with normal slashes
  if(MSYS)
    string(REGEX REPLACE "\\\\" "/" PYTHON_INSTALLATION_PATH ${PYTHON_INSTALLATION_PATH})
  endif()

  if(NOT DEFINED PYTHON_PREFIX )
    set(PYTHON_PREFIX  ${PYTHON_INSTALLATION_PATH} CACHE INTERNAL "Python installation prefix")
  endif()

  #since in the main CMakeLists.txt we used "add_subdirectory(swig EXCLUDE_FROM_ALL)" the following install command won't add anything to "make install"
  install(DIRECTORY ${PROJECT_BINARY_DIR}/python/casadi 
    DESTINATION ${PYTHON_PREFIX}
    COMPONENT install_python
    USE_SOURCE_PERMISSIONS
    PATTERN .svn EXCLUDE
  )
  #creation of a custom installation command
  if(WITH_SWIG_SPLIT)
    add_custom_target(install_python
    COMMAND ${CMAKE_COMMAND}
      -D COMPONENT=install_python
      -P cmake_install.cmake
    )
  foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  add_dependencies(install_python DEPENDS _${my_module})
  endforeach()

  else()
  add_custom_target(install_python
    DEPENDS _casadi
    COMMAND ${CMAKE_COMMAND}
      -D COMPONENT=install_python
      -P cmake_install.cmake
  )
  endif()
endif()

# OCTAVE interface
if(WITH_SWIG_SPLIT)
foreach(my_module ${CASADI_SWIG_SEPARATE_MODULES})
  SWIG_ADD_MODULE (${my_module}_octave octave ${my_module}.i octrun.swg)
  SWIG_LINK_LIBRARIES (${my_module}_octave ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control ${TINYXML_LIBRARIES} casadi_integration casadi_nonlinear_programming casadi ${CASADI_DEPENDENCIES})
  set_target_properties(${my_module}_octave PROPERTIES PREFIX "")
  set_target_properties(${my_module}_octave PROPERTIES SUFFIX ".oct")
  set_target_properties(${my_module}_octave PROPERTIES OUTPUT_NAME "${my_module}")
endforeach()
else()
  SWIG_ADD_MODULE (casadi_octave octave casadi.i octrun.swg)
  SWIG_LINK_LIBRARIES (casadi_octave ${SWIG_LINK_LIBRARIES_LIST} casadi_optimal_control ${TINYXML_LIBRARIES} casadi_integration casadi_nonlinear_programming casadi ${CASADI_DEPENDENCIES})
  set_target_properties(casadi_octave PROPERTIES PREFIX "")
  set_target_properties(casadi_octave PROPERTIES SUFFIX ".oct")
  set_target_properties(casadi_octave PROPERTIES OUTPUT_NAME "casadi_interface")
endif()

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

if(WITH_SWIG_SPLIT)
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
else(WITH_SWIG_SPLIT)
ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/octave/@swig_ref ${LIBRARY_OUTPUT_PATH}/@swig_ref
)
ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E remove ${LIBRARY_OUTPUT_PATH}/casadi.oct
)
ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/octave/casadi_single.m ${LIBRARY_OUTPUT_PATH}/casadi.m
)
ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/octave/casadi_helpers.m ${LIBRARY_OUTPUT_PATH}/casadi_helpers.m
)
endif(WITH_SWIG_SPLIT)


