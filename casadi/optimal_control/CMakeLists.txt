cmake_minimum_required(VERSION 2.6)

casadi_library(casadi_optimal_control
  variable.hpp                   variable.cpp
  symbolic_ocp.hpp               symbolic_ocp.cpp
  xml_node.hpp                   xml_node.cpp
  direct_single_shooting.hpp     direct_single_shooting.cpp     direct_single_shooting_internal.hpp     direct_single_shooting_internal.cpp
  direct_multiple_shooting.hpp   direct_multiple_shooting.cpp   direct_multiple_shooting_internal.hpp   direct_multiple_shooting_internal.cpp
  direct_collocation.hpp         direct_collocation.cpp         direct_collocation_internal.hpp         direct_collocation_internal.cpp
)
set_property(GLOBAL APPEND PROPERTY CASADI_SWIG_MODULES casadi_optimal_control)
target_link_libraries(casadi_optimal_control ${TINYXML_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/optimal_control
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
