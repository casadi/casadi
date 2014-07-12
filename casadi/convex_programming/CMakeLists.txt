cmake_minimum_required(VERSION 2.6)

casadi_library(casadi_lpsolver_qp
  qp_lp_solver.cpp     qp_lp_solver.hpp     qp_lp_internal.cpp     qp_lp_internal.hpp)

casadi_library(casadi_qpsolver_qcqp
  qcqp_qp_solver.cpp   qcqp_qp_solver.hpp   qcqp_qp_internal.cpp   qcqp_qp_internal.hpp)

casadi_library(casadi_socpsolver_sdp
  sdp_socp_solver.cpp  sdp_socp_solver.hpp  sdp_socp_internal.cpp  sdp_socp_internal.hpp)

casadi_library(casadi_stabilizedqpsolver_qp
  qp_stabilizer.cpp    qp_stabilizer.hpp    qp_stabilizer_internal.cpp      qp_stabilizer_internal.hpp)

casadi_library(casadi_linearsolver_symbolicqr
  symbolic_qr.hpp         symbolic_qr.cpp         symbolic_qr_internal.hpp         symbolic_qr_internal.cpp)

if(WITH_CSPARSE)
  casadi_library(casadi_qcqpsolver_socp
    socp_qcqp_solver.cpp socp_qcqp_solver.hpp socp_qcqp_internal.cpp socp_qcqp_internal.hpp)
  target_link_libraries(casadi_qcqpsolver_socp casadi_linearsolver_csparsecholesky)

  casadi_library(casadi_sdqpsolver_sdp
    sdp_sdqp_solver.cpp  sdp_sdqp_solver.hpp  sdp_sdqp_internal.cpp  sdp_sdqp_internal.hpp)
  target_link_libraries(casadi_sdqpsolver_sdp casadi_linearsolver_csparsecholesky)
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/convex_programming
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
