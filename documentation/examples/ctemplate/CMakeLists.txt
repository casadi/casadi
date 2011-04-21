find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(LIFTOPT QUIET)

IF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  sundials_interface
  superlu_interface 
  lapack_interface
  csparse_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  superlu
  csparse
  ${BLAS_LIBRARIES}
)
ENDIF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)
