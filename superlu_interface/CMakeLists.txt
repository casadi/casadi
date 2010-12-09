cmake_minimum_required(VERSION 2.0)

include_directories(../)
set(SUPERLU_INTERFACE_SRCS
  superlu.hpp
  superlu.cpp
  superlu_internal.hpp
  superlu_internal.cpp
)

add_library(superlu_interface STATIC ${SUPERLU_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(superlu_interface_shared SHARED ${SUPERLU_INTERFACE_SRCS})
endif(ENABLE_SHARED)
