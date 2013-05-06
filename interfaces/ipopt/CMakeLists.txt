cmake_minimum_required(VERSION 2.6)
include_directories(${IPOPT_INCLUDE_DIRS})
link_directories(${IPOPT_LIBRARY_DIRS})

SET(IPOPT_INTERFACE_SRCS
  ipopt_solver.hpp
  ipopt_solver.cpp
  ipopt_internal.hpp
  ipopt_internal.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp
)

if(WITH_IPOPT_CALLBACK)
  add_definitions(-DWITH_IPOPT_CALLBACK)
endif()

if(WITH_SIPOPT)
  add_definitions(-DWITH_SIPOPT)
endif()

if(ENABLE_STATIC)
add_library(casadi_ipopt_interface STATIC ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_ipopt_interface SHARED ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_ipopt_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/ipopt
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

