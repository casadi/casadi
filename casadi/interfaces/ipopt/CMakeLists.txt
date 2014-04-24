cmake_minimum_required(VERSION 2.6)
include_directories(${IPOPT_INCLUDE_DIRS})
link_directories(${IPOPT_LIBRARY_DIRS})

set(IPOPT_INTERFACE_SRCS
  ipopt_solver.hpp
  ipopt_solver.cpp
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

casadi_library(casadi_ipopt_interface ${IPOPT_INTERFACE_SRCS})
target_link_libraries(casadi_ipopt_interface ${IPOPT_LIBRARIES})

set_property(TARGET casadi_ipopt_interface PROPERTY PKG_CONFIG_REQUIRES ipopt)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/ipopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
