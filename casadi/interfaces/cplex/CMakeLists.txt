cmake_minimum_required(VERSION 2.6)
include_directories(${CPLEX_INCLUDE_DIR})

set(CPLEX_INTERFACE_SRCS
  cplex_solver.hpp
  cplex_solver.cpp
  cplex_internal.hpp
  cplex_internal.cpp)

if(ENABLE_STATIC)
add_library(casadi_cplex_interface STATIC ${CPLEX_INTERFACE_SRCS})
endif()
if(ENABLE_SHARED)
add_library(casadi_cplex_interface SHARED ${CPLEX_INTERFACE_SRCS})
endif()
target_link_libraries(casadi_cplex_interface casadi ${CPLEX_LIBRARIES})

set(CASADI_MODULES ${CASADI_MODULES} casadi_cplex_interface PARENT_SCOPE)

install(TARGETS casadi_cplex_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/cplex
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


