if(LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_sundials_interface
  casadi_lapack_interface
  casadi_csparse_interface
)
endif()
