cmake_minimum_required(VERSION 2.0)

SET(OPTIMAL_CONTROL
  variable.hpp              variable.cpp            variable_internal.hpp             variable_internal.cpp
  variable_tools.hpp        variable_tools.cpp
  symbolic_ocp.hpp          symbolic_ocp.cpp
  xml_node.hpp              xml_node.cpp
  ocp_tools.hpp             ocp_tools.cpp
  multiple_shooting.hpp     multiple_shooting.cpp   multiple_shooting_internal.hpp   multiple_shooting_internal.cpp
  collocation.hpp           collocation.cpp         collocation_internal.hpp         collocation_internal.cpp
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
