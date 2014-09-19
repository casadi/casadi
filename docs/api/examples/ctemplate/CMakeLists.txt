if(LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_integrator_cvodes
  casadi_integrator_idas
  casadi_linearsolver_lapacklu
  casadi_linearsolver_lapackqr
  casadi_linearsolver_csparse
)
endif()
