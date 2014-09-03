cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SQIC_INCLUDE_DIR})

embed_resource("sqic" "wsqic.f90" resource_file)

message("including: ${resource_file}")

casadi_plugin(QpSolver sqic
  sqic_interface.hpp
  sqic_interface.cpp
  sqic_interface_meta.cpp
  ${resource_file})

casadi_plugin(StabilizedQpSolver sqic
  stabilized_sqic_interface.hpp
  stabilized_sqic_interface.cpp
  stabilized_sqic_interface_meta.cpp
  ${resource_file})

if(ENABLE_STATIC)
add_library(sqic_fortran_c_bridge STATIC wsqic.f90)
endif()
if(ENABLE_SHARED)
add_library(sqic_fortran_c_bridge SHARED wsqic.f90)
endif()
casadi_plugin_link_libraries(QpSolver sqic sqic_fortran_c_bridge ${CMAKE_DL_LIBS})
casadi_plugin_link_libraries(StabilizedQpSolver sqic sqic_fortran_c_bridge ${CMAKE_DL_LIBS})
target_link_libraries(sqic_fortran_c_bridge ${SQIC_LIBRARIES})

install(TARGETS sqic_fortran_c_bridge
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sqic
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
