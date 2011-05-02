cmake_minimum_required(VERSION 2.0)

if(WITH_TIMEDVARIABLE)
add_definitions(-DWITH_TIMEDVARIABLE)
endif(WITH_TIMEDVARIABLE)

SET(OPTIMAL_CONTROL
  optimica_ocp.hpp          optimica_ocp.cpp
  variable.hpp              variable.cpp            variable_internal.hpp             variable_internal.cpp
  variable_tools.hpp        variable_tools.cpp
  fmi_parser.hpp            fmi_parser.cpp          fmi_parser_internal.hpp           fmi_parser_internal.cpp
  xml_node.hpp              xml_node.cpp
  xml_arg.hpp               xml_arg.cpp
  ocp_tools.hpp             ocp_tools.cpp
  multiple_shooting.hpp     multiple_shooting.cpp   multiple_shooting_internal.hpp   multiple_shooting_internal.cpp
)


add_library(optimal_control STATIC ${OPTIMAL_CONTROL})
if(ENABLE_SHARED)
add_library(optimal_control_shared SHARED ${OPTIMAL_CONTROL})
endif(ENABLE_SHARED)
