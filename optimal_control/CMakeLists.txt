cmake_minimum_required(VERSION 2.0)

SET(OPTIMAL_CONTROL
  symbolic_ocp.hpp
  variable.hpp              variable.cpp            variable_internal.hpp             variable_internal.cpp
  variable_tools.hpp        variable_tools.cpp
  flat_ocp.hpp              flat_ocp.cpp            flat_ocp_internal.hpp             flat_ocp_internal.cpp
  xml_node.hpp              xml_node.cpp
  xml_arg.hpp               xml_arg.cpp
  ocp_tools.hpp             ocp_tools.cpp
  multiple_shooting.hpp     multiple_shooting.cpp   multiple_shooting_internal.hpp   multiple_shooting_internal.cpp
)

if(ENABLE_STATIC)
add_library(optimal_control STATIC ${OPTIMAL_CONTROL})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(optimal_control SHARED ${OPTIMAL_CONTROL})
endif(ENABLE_SHARED)
install(TARGETS optimal_control
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/optimal_control
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
