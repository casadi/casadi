cmake_minimum_required(VERSION 2.0)
include_directories(${LIFTOPT_INCLUDE_DIRS})

SET(LIFTOPT_INTERFACE_SRCS
  liftopt_solver.hpp
  liftopt_solver.cpp
  liftopt_internal.hpp
  liftopt_internal.cpp
  casadi_lifter.hpp
  casadi_lifter.cpp
)

if(ENABLE_STATIC)
add_library(casadi_liftopt_interface STATIC ${LIFTOPT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_liftopt_interface SHARED ${LIFTOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_liftopt_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/liftopt
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

