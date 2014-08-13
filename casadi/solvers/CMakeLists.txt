cmake_minimum_required(VERSION 2.6)

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  fixed_step_integrator.hpp
  fixed_step_integrator.cpp
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_meta.cpp)

# Collocation integrator
casadi_plugin(Integrator collocation
  implicit_fixed_step_integrator.hpp
  implicit_fixed_step_integrator.cpp
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_meta.cpp
  )
target_link_libraries(casadi_integrator_collocation casadi_integrator_rk)

# Old (all-at-once) collocation integrator
casadi_plugin(Integrator oldcollocation
  old_collocation_integrator.hpp
  old_collocation_integrator.cpp
  old_collocation_integrator_meta.cpp)

casadi_plugin(ImplicitFunction newton
  newton.hpp newton.cpp newton_meta.cpp)

casadi_plugin(HomotopyNlpSolver simple
  simple_homotopy_nlp.hpp simple_homotopy_nlp.cpp simple_homotopy_nlp_meta.cpp)

# SQPMethod -  A basic SQP method
casadi_plugin(NlpSolver sqpmethod
  sqpmethod.hpp sqpmethod.cpp sqpmethod_meta.cpp)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_plugin(NlpSolver scpgen
  scpgen.hpp scpgen.cpp scpgen_meta.cpp)

# StabilizedSQP -  An implementation of Stabilized SQP
casadi_plugin(NlpSolver stabilizedsqp
  stabilized_sqp.hpp stabilized_sqp.cpp stabilized_sqp_meta.cpp)

# Reformulations
casadi_plugin(QpSolver nlp
  qp_to_nlp.hpp qp_to_nlp.cpp qp_to_nlp_meta.cpp)
casadi_plugin(ImplicitFunction nlp
  implicit_to_nlp.hpp implicit_to_nlp.cpp implicit_to_nlp_meta.cpp)
casadi_plugin(LpSolver qp
  lp_to_qp.cpp lp_to_qp.hpp lp_to_qp_meta.cpp)
casadi_plugin(QpSolver qcqp
  qp_to_qcqp.cpp qp_to_qcqp.hpp qp_to_qcqp_meta.cpp)
casadi_plugin(SocpSolver sdp
  socp_to_sdp.cpp socp_to_sdp.hpp socp_to_sdp_meta.cpp)
casadi_plugin(StabilizedQpSolver qp
  stabilized_qp_to_qp.cpp stabilized_qp_to_qp.hpp stabilized_qp_to_qp_meta.cpp)
casadi_plugin(LinearSolver symbolicqr
  symbolic_qr.hpp symbolic_qr.cpp symbolic_qr_meta.cpp
)
if(WITH_CSPARSE)
  casadi_plugin(QcqpSolver socp
    qcqp_to_socp.cpp qcqp_to_socp.hpp qcqp_to_socp_meta.cpp)
  target_link_libraries(casadi_qcqpsolver_socp casadi_linearsolver_csparsecholesky)
  casadi_plugin(SdqpSolver sdp
    sdqp_to_sdp.cpp sdqp_to_sdp.hpp sdqp_to_sdp_meta.cpp)
  target_link_libraries(casadi_sdqpsolver_sdp
    casadi_linearsolver_csparsecholesky)
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/solvers
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
