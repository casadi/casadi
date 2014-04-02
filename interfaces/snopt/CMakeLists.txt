cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SNOPT_INCLUDE_DIR})

set(SNOPT_INTERFACE_SRCS
  snopt_internal.hpp
  snopt_internal.cpp
  snopt_solver.hpp
  snopt_solver.cpp
)

set(CASADI_DEPENDENCIES ${CMAKE_DL_LIBS})

if(ENABLE_STATIC)
add_library(casadi_snopt_interface STATIC ${SNOPT_INTERFACE_SRCS})
add_library(snopt_fortran_c_bridge STATIC wsnopt.f90)
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_snopt_interface SHARED ${SQIC_INTERFACE_SRCS})
add_library(snopt_fortran_c_bridge SHARED wsnopt.f90)
endif(ENABLE_SHARED)
install(TARGETS casadi_snopt_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(TARGETS snopt_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./
   DESTINATION include/casadi/interfaces/snopt
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
