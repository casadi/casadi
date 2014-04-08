cmake_minimum_required(VERSION 2.6)

include_directories(../)
set(LAPACK_INTERFACE_SRCS
  lapack_lu_dense.hpp
  lapack_lu_dense.cpp
  lapack_qr_dense.hpp
  lapack_qr_dense.cpp)

casadi_library(casadi_lapack_interface ${LAPACK_INTERFACE_SRCS})
target_link_libraries(casadi_lapack_interface ${LAPACK_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/lapack
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
