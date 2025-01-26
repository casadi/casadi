cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(XmlFile tinyxml
  tinyxml_interface.hpp
  tinyxml_interface.cpp
  tinyxml_interface_meta.cpp)
casadi_plugin_link_libraries(XmlFile tinyxml tinyxml2::tinyxml2)
