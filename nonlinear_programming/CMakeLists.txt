cmake_minimum_required(VERSION 2.0)
SET(NONLINEAR_PROGRAMMING_SRCS
  symbolic_nlp.hpp     symbolic_nlp.cpp
  sqp_method.hpp	sqp_method.cpp		sqp_internal.hpp 		sqp_internal.cpp
  ip_method.hpp	        ip_method.cpp		ip_internal.hpp 		ip_internal.cpp
  lifted_sqp.hpp	lifted_sqp.cpp		lifted_sqp_internal.hpp 	lifted_sqp_internal.cpp
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

