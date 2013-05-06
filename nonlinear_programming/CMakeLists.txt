cmake_minimum_required(VERSION 2.6)
SET(NONLINEAR_PROGRAMMING_SRCS
  symbolic_nlp.hpp         symbolic_nlp.cpp
  sqp_method.hpp          sqp_method.cpp          sqp_internal.hpp          sqp_internal.cpp
  scpgen.hpp              scpgen.cpp              scpgen_internal.hpp       scpgen_internal.cpp
  nlp_qp_solver.hpp       nlp_qp_solver.cpp       nlp_qp_internal.hpp       nlp_qp_internal.cpp
  nlp_implicit_solver.hpp nlp_implicit_solver.cpp nlp_implicit_internal.hpp nlp_implicit_internal.cpp
  newton_implicit_solver.hpp newton_implicit_solver.cpp newton_implicit_internal.hpp newton_implicit_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_nonlinear_programming STATIC ${NONLINEAR_PROGRAMMING_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_nonlinear_programming SHARED ${NONLINEAR_PROGRAMMING_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_nonlinear_programming
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/nonlinear_programming
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

