cmake_minimum_required(VERSION 2.6)
set(CONVEX_PROGRAMMING_SRCS
  qp_lp_solver.cpp     qp_lp_solver.hpp     qp_lp_internal.cpp     qp_lp_internal.hpp
  qcqp_qp_solver.cpp   qcqp_qp_solver.hpp   qcqp_qp_internal.cpp   qcqp_qp_internal.hpp
  sdp_socp_solver.cpp  sdp_socp_solver.hpp  sdp_socp_internal.cpp  sdp_socp_internal.hpp
  qp_stabilizer.cpp    qp_stabilizer.hpp    qp_stabilizer_internal.cpp      qp_stabilizer_internal.hpp)


if(WITH_CSPARSE)
list(APPEND CONVEX_PROGRAMMING_SRCS
  socp_qcqp_solver.cpp socp_qcqp_solver.hpp socp_qcqp_internal.cpp socp_qcqp_internal.hpp
  sdp_sdqp_solver.cpp  sdp_sdqp_solver.hpp  sdp_sdqp_internal.cpp  sdp_sdqp_internal.hpp
)
endif()

casadi_library(casadi_convex_programming ${CONVEX_PROGRAMMING_SRCS})
if(WITH_CSPARSE)
  target_link_libraries(casadi_convex_programming casadi_csparse_interface)
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/convex_programming
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
