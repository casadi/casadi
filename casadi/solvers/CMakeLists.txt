cmake_minimum_required(VERSION 2.6)

# Base classes for fixed-step integrators
casadi_library(casadi_integrators
  fixed_step_integrator.hpp
  fixed_step_integrator.cpp
  implicit_fixed_step_integrator.hpp
  implicit_fixed_step_integrator.cpp)

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_meta.cpp)
target_link_libraries(casadi_integrator_rk casadi_integrators)

# Collocation integrator
casadi_plugin(Integrator collocation
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_meta.cpp)
target_link_libraries(casadi_integrator_collocation casadi_integrators)

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

casadi_plugin(DpleSolver simple
  simple_indef_dple_internal.hpp
  simple_indef_dple_internal.cpp
  simple_indef_dple_internal_meta.cpp)

casadi_plugin(DpleSolver condensing
  condensing_indef_dple_internal.hpp
  condensing_indef_dple_internal.cpp
  condensing_indef_dple_internal_meta.cpp)
  
casadi_plugin(DpleSolver lifting
  lifting_indef_dple_internal.hpp
  lifting_indef_dple_internal.cpp
  lifting_indef_dple_internal_meta.cpp)  

casadi_plugin(DleSolver simple
  simple_indef_dle_internal.hpp
  simple_indef_dle_internal.cpp
  simple_indef_dle_internal_meta.cpp)
  
casadi_plugin(LrDleSolver smith
  smith_lr_dle_internal.hpp
  smith_lr_dle_internal.cpp
  smith_lr_dle_internal_meta.cpp)

casadi_plugin(LrDpleSolver lifting
  lifting_lr_dple_internal.hpp
  lifting_lr_dple_internal.cpp
  lifting_lr_dple_internal_meta.cpp)  
  
casadi_plugin(LrDleSolver fixed_smith
  fixed_smith_lr_dle_internal.hpp
  fixed_smith_lr_dle_internal.cpp
  fixed_smith_lr_dle_internal_meta.cpp)

casadi_plugin(DleSolver fixed_smith
  fixed_smith_dle_internal.hpp
  fixed_smith_dle_internal.cpp
  fixed_smith_dle_internal_meta.cpp)
  
casadi_plugin(LrDpleSolver dple
  dple_to_lr_dple.hpp
  dple_to_lr_dple.cpp
  dple_to_lr_dple_meta.cpp) 
  
casadi_plugin(LrDleSolver dle
  dle_to_lr_dle.hpp
  dle_to_lr_dle.cpp
  dle_to_lr_dle_meta.cpp) 
  
casadi_plugin(DleSolver lrdle
  lr_dle_to_dle.hpp
  lr_dle_to_dle.cpp
  lr_dle_to_dle_meta.cpp)
  
casadi_plugin(DleSolver dple
  dple_to_dle.hpp
  dple_to_dle.cpp
  dple_to_dle_meta.cpp)

casadi_plugin(DpleSolver lrdple
  lr_dple_to_dple.hpp
  lr_dple_to_dple.cpp
  lr_dple_to_dple_meta.cpp)
  
casadi_plugin(CleSolver simple
  simple_indef_cle_internal.hpp
  simple_indef_cle_internal.cpp
  simple_indef_cle_internal_meta.cpp)  

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
  casadi_plugin(SdqpSolver sdp
    sdqp_to_sdp.cpp sdqp_to_sdp.hpp sdqp_to_sdp_meta.cpp)
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/solvers
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
