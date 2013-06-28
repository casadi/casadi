cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${WORHP_INCLUDE_DIR})

set(WORHP_INTERFACE_SRCS
  worhp_internal.hpp
  worhp_internal.cpp
  worhp_solver.hpp
  worhp_solver.cpp
)

if(ENABLE_STATIC)
add_library(casadi_worhp_interface STATIC ${WORHP_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_worhp_interface SHARED ${WORHP_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_worhp_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/worhp
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_worhp_interface PROPERTIES COMPILE_FLAGS "-Wno-unused-function")
endif()

