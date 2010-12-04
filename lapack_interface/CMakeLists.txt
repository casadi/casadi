cmake_minimum_required(VERSION 2.0)

include_directories(../)
set(LAPACK_INTERFACE_SRCS
  lapack_lu.hpp
  lapack_lu.cpp
)

add_library(lapack_interface STATIC ${LAPACK_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(lapack_interface_shared SHARED ${LAPACK_INTERFACE_SRCS})
endif(ENABLE_SHARED)
