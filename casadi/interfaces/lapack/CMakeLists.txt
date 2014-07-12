cmake_minimum_required(VERSION 2.6)
include_directories(../)

# Linear solver using dense LU factorization in LAPACK
casadi_library(casadi_linearsolver_lapacklu
  lapack_lu_dense.hpp
  lapack_lu_dense.cpp)
target_link_libraries(casadi_linearsolver_lapacklu ${LAPACK_LIBRARIES})

# Linear solver using dense QR factorization in LAPACK
casadi_library(casadi_linearsolver_lapackqr
  lapack_qr_dense.hpp
  lapack_qr_dense.cpp)
target_link_libraries(casadi_linearsolver_lapackqr ${LAPACK_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/lapack
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
