cmake_minimum_required(VERSION 2.0)

include_directories(../)
set(SUPERLU_INTERFACE_SRCS
  superlu.hpp
  superlu.cpp
  superlu_internal.hpp
  superlu_internal.cpp
)

if(ENABLE_STATIC)
add_library(superlu_interface STATIC ${SUPERLU_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(superlu_interface SHARED ${SUPERLU_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS superlu_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
