cmake_minimum_required(VERSION 2.6)
set(NONLINEAR_PROGRAMMING_SRCS
  symbolic_nlp.hpp         symbolic_nlp.cpp
  sqp_method.hpp          sqp_method.cpp          sqp_internal.hpp          sqp_internal.cpp
  stabilized_sqp_method.hpp          stabilized_sqp_method.cpp          stabilized_sqp_internal.hpp          stabilized_sqp_internal.cpp
  scpgen.hpp              scpgen.cpp              scpgen_internal.hpp       scpgen_internal.cpp
  nlp_qp_solver.hpp       nlp_qp_solver.cpp       nlp_qp_internal.hpp       nlp_qp_internal.cpp
  nlp_implicit_solver.hpp nlp_implicit_solver.cpp nlp_implicit_internal.hpp nlp_implicit_internal.cpp
  newton_implicit_solver.hpp newton_implicit_solver.cpp newton_implicit_internal.hpp newton_implicit_internal.cpp
  simple_homotopy_nlp_solver.hpp simple_homotopy_nlp_solver.cpp simple_homotopy_nlp_internal.hpp simple_homotopy_nlp_internal.cpp
)

casadi_library(casadi_nonlinear_programming ${NONLINEAR_PROGRAMMING_SRCS})

install(DIRECTORY ./
  DESTINATION include/casadi/nonlinear_programming
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
