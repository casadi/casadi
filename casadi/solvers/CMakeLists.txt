cmake_minimum_required(VERSION 2.6)

# Explicit Runge-Kutta integrator
casadi_plugin(Integrator rk
  fixed_step_integrator_internal.hpp
  fixed_step_integrator_internal.cpp
  rk_integrator_internal.hpp
  rk_integrator_internal.cpp
  rk_integrator_internal_meta.cpp)

# Collocation integrator
casadi_plugin(Integrator collocation
  implicit_fixed_step_integrator_internal.hpp
  implicit_fixed_step_integrator_internal.cpp
  collocation_integrator_internal.hpp
  collocation_integrator_internal.cpp
  collocation_integrator_internal_meta.cpp
  )
target_link_libraries(casadi_integrator_collocation casadi_integrator_rk)

# Old (all-at-once) collocation integrator
casadi_plugin(Integrator oldcollocation
  old_collocation_integrator_internal.hpp
  old_collocation_integrator_internal.cpp
  old_collocation_integrator_internal_meta.cpp)

casadi_plugin(ImplicitFunction newton
  newton_implicit_internal.hpp   newton_implicit_internal.cpp
  newton_implicit_internal_meta.cpp)

casadi_plugin(HomotopyNlpSolver simple
   simple_homotopy_nlp_internal.hpp   simple_homotopy_nlp_internal.cpp
   simple_homotopy_nlp_internal_meta.cpp
   )

# SQPMethod -  A basic SQP method
casadi_plugin(NlpSolver sqpmethod   sqp_internal.hpp   sqp_internal.cpp
   sqp_internal_meta.cpp
)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_plugin(NlpSolver scpgen   scpgen_internal.hpp   scpgen_internal.cpp
  scpgen_internal_meta.cpp
)

# StabilizedSQP -  An implementation of Stabilized SQP
casadi_plugin(NlpSolver stabilizedsqp   stabilized_sqp_internal.hpp   stabilized_sqp_internal.cpp
  stabilized_sqp_internal_meta.cpp
)

# Reformulations
casadi_plugin(QpSolver nlp nlp_qp_internal.hpp   nlp_qp_internal.cpp
nlp_qp_internal_meta.cpp
)
casadi_plugin(ImplicitFunction nlp   nlp_implicit_internal.hpp   nlp_implicit_internal.cpp
nlp_implicit_internal_meta.cpp
)
casadi_plugin(LpSolver qp   qp_lp_internal.cpp   qp_lp_internal.hpp
qp_lp_internal_meta.cpp
)
casadi_plugin(QpSolver qcqp   qcqp_qp_internal.cpp   qcqp_qp_internal.hpp
qcqp_qp_internal_meta.cpp
)
casadi_plugin(SocpSolver sdp   sdp_socp_internal.cpp   sdp_socp_internal.hpp
sdp_socp_internal_meta.cpp
)
casadi_plugin(StabilizedQpSolver qp   qp_stabilizer_internal.cpp   qp_stabilizer_internal.hpp
qp_stabilizer_internal_meta.cpp
)
casadi_plugin(LinearSolver symbolicqr   symbolic_qr_internal.hpp   symbolic_qr_internal.cpp
symbolic_qr_internal_meta.cpp
)
if(WITH_CSPARSE)
  casadi_plugin(QcqpSolver socp   socp_qcqp_internal.cpp   socp_qcqp_internal.hpp
    socp_qcqp_internal_meta.cpp
  )
  target_link_libraries(casadi_qcqpsolver_socp   casadi_linearsolver_csparsecholesky
  )
  casadi_plugin(SdqpSolver sdp   sdp_sdqp_internal.cpp    sdp_sdqp_internal.hpp
    sdp_sdqp_internal_meta.cpp
  )
  target_link_libraries(casadi_sdqpsolver_sdp   casadi_linearsolver_csparsecholesky
  )
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/solvers
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
