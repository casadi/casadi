cmake_minimum_required(VERSION 2.6)

SET(OPTIMAL_CONTROL
  variable.hpp                   variable.cpp
  symbolic_ocp.hpp               symbolic_ocp.cpp
  xml_node.hpp                   xml_node.cpp
  direct_single_shooting.hpp     direct_single_shooting.cpp     direct_single_shooting_internal.hpp     direct_single_shooting_internal.cpp
  direct_multiple_shooting.hpp   direct_multiple_shooting.cpp   direct_multiple_shooting_internal.hpp   direct_multiple_shooting_internal.cpp
  direct_collocation.hpp         direct_collocation.cpp         direct_collocation_internal.hpp         direct_collocation_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_optimal_control STATIC ${OPTIMAL_CONTROL})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_optimal_control SHARED ${OPTIMAL_CONTROL})
target_link_libraries(casadi_optimal_control
  casadi
  casadi_integration
  ${TINYXML_LIBRARIES}
)
endif(ENABLE_SHARED)
install(TARGETS casadi_optimal_control
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/optimal_control
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
