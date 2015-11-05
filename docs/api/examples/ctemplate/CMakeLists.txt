if(LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_ivpsol_cvodes
  casadi_ivpsol_idas
  casadi_linearsolver_lapacklu
  casadi_linearsolver_lapackqr
  casadi_linearsolver_csparse
)
endif()
