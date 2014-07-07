if(LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_integrator_cvodes
  casadi_integrator_idas
  casadi_lapack_interface
  casadi_csparse_interface
)
endif()
