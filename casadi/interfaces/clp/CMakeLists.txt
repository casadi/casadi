cmake_minimum_required(VERSION 3.10.2)

# Workaround to be removed when cmake version is high enough for either:
#  - PkgConfig to support absolute paths
#  - target_link_directories to be supported
if(CLP_LIBRARY_DIRS)
  link_directories(${CLP_LIBRARY_DIRS})
endif()

casadi_plugin(Conic clp
  clp_interface.hpp
  clp_interface.cpp
  clp_interface_meta.cpp)
casadi_plugin_link_libraries(Conic clp clp)
