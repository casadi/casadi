cmake_minimum_required(VERSION 2.0)
SET(NONLINEAR_PROGRAMMING_SRCS
  sqp_method.hpp
  sqp_method.cpp
  sqp_internal.hpp
  sqp_internal.cpp
)

if(ENABLE_STATIC)
add_library(nonlinear_programming STATIC ${NONLINEAR_PROGRAMMING_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(nonlinear_programming SHARED ${NONLINEAR_PROGRAMMING_SRCS})
endif(ENABLE_SHARED)
install(TARGETS nonlinear_programming
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/nonlinear_programming
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

