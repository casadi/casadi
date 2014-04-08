cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SQIC_INCLUDE_DIR})

embed_resource("sqic" "wsqic.f90" resource_file)

message("including: ${resource_file}")

set(SQIC_INTERFACE_SRCS
  sqic_internal.hpp
  sqic_internal.cpp
  sqic_solver.hpp
  sqic_solver.cpp
  stabilized_sqic_internal.hpp
  stabilized_sqic_internal.cpp
  stabilized_sqic_solver.hpp
  stabilized_sqic_solver.cpp
  ${resource_file})

casadi_library(casadi_sqic_interface ${SQIC_INTERFACE_SRCS})
if(ENABLE_STATIC)
add_library(sqic_fortran_c_bridge STATIC wsqic.f90)
endif()
if(ENABLE_SHARED)
add_library(sqic_fortran_c_bridge SHARED wsqic.f90)
endif()
target_link_libraries(casadi_sqic_interface sqic_fortran_c_bridge ${CMAKE_DL_LIBS})
target_link_libraries(sqic_fortran_c_bridge ${SQIC_LIBRARIES})

install(TARGETS sqic_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sqic
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
