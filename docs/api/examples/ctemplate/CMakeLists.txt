if(WITH_LAPACK AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_integrator_cvodes
  casadi_integrator_idas
  casadi_linsol_lapacklu
  casadi_linsol_lapackqr
  casadi_linsol_csparse
)
endif()
