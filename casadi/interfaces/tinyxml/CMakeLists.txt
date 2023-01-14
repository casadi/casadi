cmake_minimum_required(VERSION 2.8.12)

casadi_plugin(XmlFile tinyxml
  tinyxml_interface.hpp
  tinyxml_interface.cpp
  tinyxml_interface_meta.cpp)
casadi_plugin_link_libraries(XmlFile tinyxml tinyxml2)

install(TARGETS tinyxml2
RUNTIME DESTINATION ${LIB_PREFIX}
LIBRARY DESTINATION ${LIB_PREFIX})