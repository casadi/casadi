cmake_minimum_required(VERSION 2.6)
SET(CONTROL_SRCS
  dple_solver.hpp
  dple_solver.cpp
  dple_internal.hpp
  dple_internal.cpp
  simple_indef_dple_solver.hpp
  simple_indef_dple_solver.cpp
  simple_indef_dple_internal.hpp
  simple_indef_dple_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_control STATIC ${CONTROL_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_control SHARED ${CONTROL_SRCS})
endif(ENABLE_SHARED)
target_link_libraries(casadi_control casadi)
set(CASADI_MODULES ${CASADI_MODULES} casadi_control PARENT_SCOPE)
install(TARGETS casadi_control
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/control
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
