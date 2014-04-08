cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SNOPT_INCLUDE_DIR})

set(SNOPT_INTERFACE_SRCS
  snopt_internal.hpp
  snopt_internal.cpp
  snopt_solver.hpp
  snopt_solver.cpp)

casadi_library(casadi_snopt_interface ${SNOPT_INTERFACE_SRCS})

if(ENABLE_STATIC)
add_library(snopt_fortran_c_bridge STATIC wsnopt.f90)
endif()
if(ENABLE_SHARED)
add_library(snopt_fortran_c_bridge SHARED wsnopt.f90)
endif()
target_link_libraries(casadi_snopt_interface snopt_fortran_c_bridge ${CMAKE_DL_LIBS})
target_link_libraries(snopt_fortran_c_bridge ${SNOPT_LIBRARIES})

install(TARGETS snopt_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/snopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

