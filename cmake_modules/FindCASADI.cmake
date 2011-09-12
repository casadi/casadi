FIND_PATH(CASADI_INCLUDE_DIR 
  casadi/casadi/sx/sx.hpp
)

IF (CASADI_INCLUDE_DIR)
  SET(CASADI_FOUND_INCLUDE TRUE)
  MESSAGE(STATUS "Found CasADi include dir: ${CASADI_INCLUDE_DIR}")
ELSE (CASADI_INCLUDE_DIR)
  MESSAGE(STATUS "Could not find CasADi include dir")
ENDIF (CASADI_INCLUDE_DIR)

SET(CASADI_LIBS_LIST
  cplex_interface
  ipopt_interface
  lapack_interface
  sundials_interface
  superlu_interface
  csparse_interface
  knitro_interface
  optimal_control
  superlu
  csparse
  tinyxml
  casadi
)

FOREACH(LIB in ${CASADI_LIBS_LIST})
  FIND_LIBRARY(CASADI_LIB_${LIB}
    NAMES ${LIB}
    HINTS ${CASADI_INCLUDE_DIR}/build/lib)
  IF(CASADI_LIB_${LIB})
    #MESSAGE(STATUS "Found ${LIB}: ${CASADI_LIB_${LIB}}")
    SET(CASADI_LIBRARIES ${CASADI_LIBRARIES} ${CASADI_LIB_${LIB}})
  ELSE(CASADI_LIB_${LIB})
    #MESSAGE(STATUS "Could not find lib${LIB}")
  ENDIF(CASADI_LIB_${LIB})
ENDFOREACH(LIB)

IF (CASADI_LIBRARIES)
   MESSAGE(STATUS "Found CasADi libs: ${CASADI_LIBRARIES}")
ELSE (CASADI_LIBRARIES)
   MESSAGE(STATUS "Could not find CasADi libs")
ENDIF (CASADI_LIBRARIES)

IF(CASADI_FOUND_INCLUDE AND CASADI_LIBRARIES)
  SET(CASADI_FOUND TRUE)
ENDIF(CASADI_FOUND_INCLUDE AND CASADI_LIBRARIES)





