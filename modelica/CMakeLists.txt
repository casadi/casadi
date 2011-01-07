cmake_minimum_required(VERSION 2.0)

SET(MODELICA
  optimica_ocp.hpp
  optimica_ocp.cpp
  ocp_variables.hpp  
  ocp_variables.cpp
  variable.hpp
  variable.cpp
  variable_internal.hpp
  variable_internal.cpp
  variable_tools.hpp
  variable_tools.cpp
  fmi_parser.hpp
  fmi_parser.cpp
  fmi_parser_internal.hpp
  fmi_parser_internal.cpp
  xml_node.hpp
  xml_node.cpp
  xml_arg.hpp
  xml_arg.cpp
)


add_library(modelica STATIC ${MODELICA})
if(ENABLE_SHARED)
add_library(modelica_shared SHARED ${MODELICA})
endif(ENABLE_SHARED)
