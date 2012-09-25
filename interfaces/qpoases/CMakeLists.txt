cmake_minimum_required(VERSION 2.0)

find_package(QPOASES REQUIRED QUIET)
include_directories(${QPOASES_INCLUDE_DIR})

set(QPOASES_INTERFACE_SRCS
  qpoases_internal.hpp
  qpoases_internal.cpp
  qpoases_solver.hpp
  qpoases_solver.cpp
)

if(ENABLE_STATIC)
add_library(qpoases_interface STATIC ${QPOASES_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(qpoases_interface SHARED ${QPOASES_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS qpoases_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/qpoases
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
