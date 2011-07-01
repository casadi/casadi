cmake_minimum_required(VERSION 2.0)

include_directories(../)
set(LAPACK_INTERFACE_SRCS
  lapack_lu_dense.hpp
  lapack_lu_dense.cpp
  lapack_qr_dense.hpp
  lapack_qr_dense.cpp
)

if(ENABLE_STATIC)
add_library(lapack_interface STATIC ${LAPACK_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(lapack_interface SHARED ${LAPACK_INTERFACE_SRCS})
endif(ENABLE_SHARED)
