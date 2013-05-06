cmake_minimum_required(VERSION 2.6)
include_directories(${CSPARSE_INCLUDE_DIR})

set(CSPARSE_INTERFACE_SRCS
  csparse.hpp
  csparse.cpp
  csparse_internal.hpp
  csparse_internal.cpp
)

if(ENABLE_STATIC)
  add_library(casadi_csparse_interface STATIC ${CSPARSE_INTERFACE_SRCS})
endif()

if(ENABLE_SHARED)
  add_library(casadi_csparse_interface SHARED ${CSPARSE_INTERFACE_SRCS})
endif()

install(TARGETS casadi_csparse_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/csparse
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

