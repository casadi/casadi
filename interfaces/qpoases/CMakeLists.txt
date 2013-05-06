cmake_minimum_required(VERSION 2.6)
include_directories(${QPOASES_INCLUDE_DIR})

set(QPOASES_INTERFACE_SRCS
  qpoases_internal.hpp
  qpoases_internal.cpp
  qpoases_solver.hpp
  qpoases_solver.cpp
)

if(ENABLE_STATIC)
add_library(casadi_qpoases_interface STATIC ${QPOASES_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_qpoases_interface SHARED ${QPOASES_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_qpoases_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/qpoases
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpoases_interface PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
