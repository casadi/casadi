cmake_minimum_required(VERSION 2.0)

set(CSPARSE_INTERFACE_SRCS
  csparse.hpp
  csparse.cpp
  csparse_internal.hpp
  csparse_internal.cpp
)

add_library(csparse_interface STATIC ${CSPARSE_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(csparse_interface_shared SHARED ${CSPARSE_INTERFACE_SRCS})
endif(ENABLE_SHARED)
