cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SQIC_INCLUDE_DIR})

set(SQIC_INTERFACE_SRCS
  sqic_internal.hpp
  sqic_internal.cpp
  sqic_solver.hpp
  sqic_solver.cpp
)

set(CASADI_DEPENDENCIES ${CMAKE_DL_LIBS})

if(ENABLE_STATIC)
add_library(casadi_sqic_interface STATIC ${SQIC_INTERFACE_SRCS})
add_library(sqic_fortran_c_bridge STATIC wsqic.f90)
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_sqic_interface SHARED ${SQIC_INTERFACE_SRCS})
add_library(sqic_fortran_c_bridge SHARED wsqic.f90)
endif(ENABLE_SHARED)
install(TARGETS casadi_sqic_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(TARGETS sqic_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/sqic
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

