cmake_minimum_required(VERSION 3.10.2)

# Workaround to be removed when cmake version is high enough for either:
#  - PkgConfig to support absolute paths
#  - target_link_directories to be supported
if(CBC_LIBRARY_DIRS)
  link_directories(${CBC_LIBRARY_DIRS})
endif()

casadi_plugin(Conic cbc
  cbc_interface.hpp
  cbc_interface.cpp
  cbc_interface_meta.cpp)
casadi_plugin_link_libraries(Conic cbc cbc)

