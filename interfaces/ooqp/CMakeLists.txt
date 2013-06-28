cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${OOQP_INCLUDE_DIRS})

set(OOQP_INTERFACE_SRCS
  ooqp_internal.hpp
  ooqp_internal.cpp
  ooqp_solver.hpp
  ooqp_solver.cpp
)

if(ENABLE_STATIC)
add_library(casadi_ooqp_interface STATIC ${OOQP_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_ooqp_interface SHARED ${OOQP_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_ooqp_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/ooqp
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

