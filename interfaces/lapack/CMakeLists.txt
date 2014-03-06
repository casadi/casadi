cmake_minimum_required(VERSION 2.6)

include_directories(../)
set(LAPACK_INTERFACE_SRCS
  lapack_lu_dense.hpp
  lapack_lu_dense.cpp
  lapack_qr_dense.hpp
  lapack_qr_dense.cpp
)

if(ENABLE_STATIC)
add_library(casadi_lapack_interface STATIC ${LAPACK_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_lapack_interface SHARED ${LAPACK_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_lapack_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/lapack
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

