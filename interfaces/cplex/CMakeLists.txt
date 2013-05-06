cmake_minimum_required(VERSION 2.6)
include_directories(${CPLEX_INCLUDE_DIR})

SET(CPLEX_INTERFACE_SRCS
  cplex_solver.hpp
  cplex_solver.cpp
  cplex_internal.hpp
  cplex_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_cplex_interface STATIC ${CPLEX_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_cplex_interface SHARED ${CPLEX_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_cplex_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./
   DESTINATION include/casadi/interfaces/cplex
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)


