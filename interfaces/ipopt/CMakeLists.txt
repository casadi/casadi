cmake_minimum_required(VERSION 2.0)
include_directories(${IPOPT_INCLUDE_DIRS})
link_directories(${IPOPT_LIBRARY_DIRS})

SET(IPOPT_INTERFACE_SRCS
  ipopt_solver.hpp
  ipopt_solver.cpp
  ipopt_internal.hpp
  ipopt_internal.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp
  ipopt_qp_internal.hpp
  ipopt_qp_internal.cpp
  ipopt_qp_solver.hpp
  ipopt_qp_solver.cpp
)

IF(IPOPT_FOUND_EXTRA)
  add_definitions(-DWITH_IPOPT_CALLBACK)
ENDIF(IPOPT_FOUND_EXTRA)


if(ENABLE_STATIC)
add_library(ipopt_interface STATIC ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(ipopt_interface SHARED ${IPOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS ipopt_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/ipopt
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

