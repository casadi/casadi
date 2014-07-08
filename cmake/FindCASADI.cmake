find_path(CASADI_INCLUDE_DIR 
  casadi/core/sx/sx_element.hpp
  HINTS $ENV{CASADI_PREFIX}/include
)

if(CASADI_INCLUDE_DIR)
  set(CASADI_INCLUDE_DIR ${CASADI_INCLUDE_DIR} ${CASADI_INCLUDE_DIR}/casadi)
  set(CASADI_FOUND_INCLUDE TRUE)
  message(STATUS "Found CasADi include dir: ${CASADI_INCLUDE_DIR}")
else()
  message(STATUS "Could not find CasADi include dir")
endif()

set(CASADI_LIBS_LIST
  casadi_cplex_interface
  casadi_nlpsolver_ipopt
  casadi_lapack_interface
  casadi_sundials_interface
  casadi_csparse_interface
  casadi_knitro_interface
  casadi_optimal_control
  casadi_integration
  casadi_nonlinear_programming
  casadi_csparse
  casadi_tinyxml
  casadi)

foreach(LIB in ${CASADI_LIBS_LIST})
  find_library(CASADI_LIB_${LIB}
    NAMES ${LIB}
    HINTS ${CASADI_INCLUDE_DIR}/build/lib $ENV{CASADI_PREFIX}/lib)
  if(CASADI_LIB_${LIB})
    #message(STATUS "Found ${LIB}: ${CASADI_LIB_${LIB}}")
    set(CASADI_LIBRARIES ${CASADI_LIBRARIES} ${CASADI_LIB_${LIB}})
  else()
    #message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(CASADI_LIBRARIES)
  message(STATUS "Found CasADi libs: ${CASADI_LIBRARIES}")
else()
  message(STATUS "Could not find CasADi libs")
endif()

if(CASADI_FOUND_INCLUDE AND CASADI_LIBRARIES)
  set(CASADI_FOUND TRUE)
endif()
