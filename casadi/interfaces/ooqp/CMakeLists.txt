cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${OOQP_INCLUDE_DIRS})

casadi_plugin(QpSolver ooqp
  ooqp_interface.hpp
  ooqp_interface.cpp
  ooqp_interface_meta.cpp)

casadi_plugin_link_libraries(QpSolver ooqp ${OOQP_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/ooqp
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)


