cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SQIC_INCLUDE_DIR})

embed_resource("sqic" "wsqic.f90" resource_file)

message("including: ${resource_file}")

casadi_library(casadi_qpsolver_sqic
  sqic_internal.hpp
  sqic_internal.cpp
  ${resource_file})

casadi_library(casadi_stabilizedqpsolver_sqic
  stabilized_sqic_internal.hpp
  stabilized_sqic_internal.cpp)
target_link_libraries(casadi_stabilizedqpsolver_sqic casadi_qpsolver_sqic)

if(ENABLE_STATIC)
add_library(sqic_fortran_c_bridge STATIC wsqic.f90)
endif()
if(ENABLE_SHARED)
add_library(sqic_fortran_c_bridge SHARED wsqic.f90)
endif()
target_link_libraries(casadi_qpsolver_sqic sqic_fortran_c_bridge ${CMAKE_DL_LIBS})
target_link_libraries(casadi_stabilizedqpsolver_sqic sqic_fortran_c_bridge ${CMAKE_DL_LIBS})
target_link_libraries(sqic_fortran_c_bridge ${SQIC_LIBRARIES})

install(TARGETS sqic_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sqic
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
