cmake_minimum_required(VERSION 2.6)

casadi_library(casadi_nonlinear_programming
  symbolic_nlp.hpp         symbolic_nlp.cpp
  nlp_qp_solver.hpp       nlp_qp_solver.cpp       nlp_qp_internal.hpp       nlp_qp_internal.cpp
  nlp_implicit_solver.hpp nlp_implicit_solver.cpp nlp_implicit_internal.hpp nlp_implicit_internal.cpp
  newton_implicit_solver.hpp newton_implicit_solver.cpp newton_implicit_internal.hpp newton_implicit_internal.cpp
  simple_homotopy_nlp_solver.hpp simple_homotopy_nlp_solver.cpp simple_homotopy_nlp_internal.hpp simple_homotopy_nlp_internal.cpp
)
target_link_libraries(casadi_nonlinear_programming casadi)

# SQPMethod -  A basic SQP method
casadi_library(casadi_nlpsolver_sqpmethod
  sqp_method.hpp          sqp_method.cpp          sqp_internal.hpp          sqp_internal.cpp)
target_link_libraries(casadi_nlpsolver_sqpmethod casadi_nonlinear_programming)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_library(casadi_nlpsolver_scpgen
  scpgen.hpp              scpgen.cpp              scpgen_internal.hpp       scpgen_internal.cpp)
target_link_libraries(casadi_nlpsolver_scpgen casadi_nonlinear_programming)

# StabilizedSQP -  An implementation of Stabilized SQP
casadi_library(casadi_nlpsolver_stabilizedsqp
  stabilized_sqp_method.hpp          stabilized_sqp_method.cpp          stabilized_sqp_internal.hpp          stabilized_sqp_internal.cpp)
target_link_libraries(casadi_nlpsolver_stabilizedsqp casadi_nonlinear_programming)

install(DIRECTORY ./
  DESTINATION include/casadi/nonlinear_programming
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

