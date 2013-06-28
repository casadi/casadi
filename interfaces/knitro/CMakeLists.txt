cmake_minimum_required(VERSION 2.6)
include_directories(${KNITRO_INCLUDE_DIR})

SET(KNITRO_INTERFACE_SRCS
  knitro_solver.hpp
  knitro_solver.cpp
  knitro_internal.hpp
  knitro_internal.cpp
)


if(ENABLE_STATIC)
add_library(casadi_knitro_interface STATIC ${KNITRO_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_knitro_interface SHARED ${KNITRO_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_knitro_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/knitro
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

