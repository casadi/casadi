cmake_minimum_required(VERSION 2.6)
include_directories(${IPOPT_INCLUDE_DIRS})
link_directories(${IPOPT_LIBRARY_DIRS})

set(NLPSOLVER_IPOPT_SRCS
  ipopt_internal.hpp
  ipopt_internal.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp)
  
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
else()
  add_definitions(-DHAVE_CSTDDEF)
endif()

if(WITH_IPOPT_CALLBACK)
  add_definitions(-DWITH_IPOPT_CALLBACK)
endif()

if(IPOPT_FOUND)
  add_definitions(${IPOPT_CFLAGS_OTHER})
endif()

if(WITH_SIPOPT)
  add_definitions(-DWITH_SIPOPT)
endif()

casadi_library(casadi_nlpsolver_ipopt ${NLPSOLVER_IPOPT_SRCS})
target_link_libraries(casadi_nlpsolver_ipopt ${IPOPT_LIBRARIES})

set_property(TARGET casadi_nlpsolver_ipopt PROPERTY PKG_CONFIG_REQUIRES ipopt)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/ipopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
