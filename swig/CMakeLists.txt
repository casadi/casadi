cmake_minimum_required (VERSION 2.0)

#if(COMMAND cmake_policy)
#cmake_policy(SET CMP0003 NEW)
#endif(COMMAND cmake_policy)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules)
project(casadi_swig C CXX)

find_package (SWIG REQUIRED QUIET)
find_package(MPI QUIET)
if(MPI_FOUND)
   include_directories (${MPI_INCLUDE_PATH})
   set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
endif(MPI_FOUND)

set(CASADI_SWIG_FLAGS -includeall)

find_package (IPOPT QUIET)
if(IPOPT_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_IPOPT")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} ipopt_interface ${IPOPT_LIBRARIES})
endif(IPOPT_FOUND)

find_package (SUNDIALS QUIET)
if(SUNDIALS_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SUNDIALS")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} sundials_interface ${SUNDIALS_LIBRARIES})
endif(SUNDIALS_FOUND)

find_package (ACADO QUIET)
if(ACADO_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_ACADO")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} acado_interface ${ACADO_LIBRARIES})
endif(ACADO_FOUND)

if(WITH_SUPERLU)
find_package (BLAS QUIET)
if(BLAS_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_SUPERLU")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} superlu_interface superlu)
endif(BLAS_FOUND)
endif(WITH_SUPERLU)

set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CSPARSE")
set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} csparse_interface csparse)

find_package (LAPACK QUIET)
if(LAPACK_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_LAPACK")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} lapack_interface ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

find_package (KNITRO QUIET)
if(KNITRO_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_KNITRO")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} knitro_interface ${KNITRO_LIBRARIES})
endif(KNITRO_FOUND)

find_package (CPLEX QUIET)
if(CPLEX_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_CPLEX")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} cplex_interface ${CPLEX_LIBRARIES})
endif(CPLEX_FOUND)

find_package (LIFTOPT QUIET)
if(LIFTOPT_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_LIFTOPT")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} liftopt_interface ${LIFTOPT_LIBRARIES})
endif(LIFTOPT_FOUND)

if(WITH_GSL)
find_package (GSL)
find_package (BLAS)
if(GSL_FOUND)
if(BLAS_FOUND)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_GSL")
  set(SWIG_LINK_LIBRARIES_LIST ${SWIG_LINK_LIBRARIES_LIST} gsl_interface ${GSL_LIBRARIES})
endif(BLAS_FOUND)
endif(GSL_FOUND)
endif(WITH_GSL)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../doc/doc.i)
  set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_DOXDOC")
endif()

# Find SWIG and include the use swig file
include (${SWIG_USE_FILE})

# Find python library and add include path for python headers
find_package (PythonLibs REQUIRED)

include_directories (${PYTHON_INCLUDE_PATH})

find_package (NUMPY)
if(NUMPY_FOUND)
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${NUMPY_INCLUDED_DIRS})
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_NUMPY")
endif(NUMPY_FOUND)

if(WITH_PYTHON_INTERRUPTS)
set(CASADI_SWIG_FLAGS ${CASADI_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif(WITH_PYTHON_INTERRUPTS)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories (${CMAKE_CURRENT_SOURCE_DIR}/../build/lib)

# set the global swig flags to empty
# set(CMAKE_SWIG_FLAGS "-outdir ${PROJECT_BINARY_DIR}")
# SEPARATE_ARGUMENTS(CMAKE_SWIG_FLAGS)

# SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
# SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# SEPARATE_ARGUMENTS(CASADI_SWIG_FLAGS)
set_source_files_properties (casadi.i PROPERTIES
  SWIG_FLAGS "-DPy_USING_UNICODE ${CASADI_SWIG_FLAGS}"
  CPLUSPLUS ON
)

# PYTHON interface
SWIG_ADD_MODULE (casadi python casadi.i) #the target corresponding this module is called "_casadi"
# add_dependencies(_casadi casadi) #workaround to add dependencies of casadi to _casadi (SWIG_ADD_MODULE does not add dependencies from c++ source files)
IF(NUMPY_FOUND)
cmake_policy(SET CMP0003 NEW)
SWIG_LINK_LIBRARIES (casadi ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} optimal_control tinyxml casadi )
ELSE(NUMPY_FOUND)
SWIG_LINK_LIBRARIES (casadi ${PYTHON_LIBRARIES}  ${SWIG_LINK_LIBRARIES_LIST} optimal_control tinyxml casadi )
ENDIF(NUMPY_FOUND)
# the name used for the shared module for the python interface is determined
get_target_property(CASADI_PYTHON_LIBRARY_SUFFIX ${SWIG_MODULE_casadi_REAL_NAME} SUFFIX)
if(NOT ${CASADI_PYTHON_LIBRARY_SUFFIX})
  set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
endif(NOT ${CASADI_PYTHON_LIBRARY_SUFFIX})
set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

# a python library is built in the build directory inside swig/python 
make_directory(${PROJECT_BINARY_DIR}/python/casadi)
# adding some extra commands to the target _casadi which is declared automatically when constructing a swig-python interface
add_custom_command(
  TARGET _casadi
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py ${PROJECT_BINARY_DIR}/python/casadi/__init__.py
)
add_custom_command(
  TARGET _casadi
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${PROJECT_BINARY_DIR}/casadi.py ${PROJECT_BINARY_DIR}/python/casadi/casadi.py
)
add_custom_command(
  TARGET _casadi
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy ${LIBRARY_OUTPUT_PATH}/${CASADI_PYTHON_LIBRARY} ${PROJECT_BINARY_DIR}/python/casadi/${CASADI_PYTHON_LIBRARY}
)
find_package(PythonInterp)
IF(PYTHONINTERP_FOUND)
  #determining the installation path for python libraries
  EXEC_PROGRAM( "${PYTHON_EXECUTABLE}"
    ARGS "${CMAKE_CURRENT_SOURCE_DIR}/python/find_installation_path.py"
    OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
  )
  #since in the main CMakeLists.txt we used "add_subdirectory(swig EXCLUDE_FROM_ALL)" the following install command won't add anything to "make install"
  install(DIRECTORY ${PROJECT_BINARY_DIR}/python/casadi 
    DESTINATION ${PYTHON_INSTALLATION_PATH}
    COMPONENT install_python
  )
  #creation of a custom installation command
  add_custom_target(install_python
    DEPENDS _casadi
    COMMAND ${CMAKE_COMMAND}
      -D COMPONENT=install_python
      -P cmake_install.cmake
  )
ENDIF(PYTHONINTERP_FOUND)

# OCTAVE interface
SWIG_ADD_MODULE (casadi_octave octave casadi.i)
SWIG_LINK_LIBRARIES (casadi_octave ${SWIG_LINK_LIBRARIES_LIST} optimal_control tinyxml casadi )
set_target_properties(casadi_octave PROPERTIES PREFIX "")
set_target_properties(casadi_octave PROPERTIES SUFFIX ".oct")
set_target_properties(casadi_octave PROPERTIES OUTPUT_NAME "casadi_interface")

ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/xcasadi.m ${LIBRARY_OUTPUT_PATH}/casadi.m
)
 
ADD_CUSTOM_COMMAND(
 TARGET casadi_octave
 POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E remove ${LIBRARY_OUTPUT_PATH}/casadi.oct
)


