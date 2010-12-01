cmake_minimum_required(VERSION 2.0)

SET(MODELICA
  optimica_ocp.hpp
  optimica_ocp.cpp
  variable.hpp
  variable.cpp
  fmi_parser.hpp
  fmi_parser.cpp
  xml_parser.hpp
  xml_parser.cpp
  xml_node.hpp
  xml_node.cpp
  xml_arg.hpp
  xml_arg.cpp
)


add_library(modelica STATIC ${MODELICA})
if(ENABLE_SHARED)
add_library(modelica_shared SHARED ${MODELICA})
endif(ENABLE_SHARED)

add_subdirectory(external_packages/tinyxml)
add_subdirectory(test)
