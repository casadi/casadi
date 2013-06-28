cmake_minimum_required(VERSION 2.6)
SET(INTEGRATION_SRCS
  collocation_integrator.hpp
  collocation_integrator.cpp
  collocation_integrator_internal.hpp
  collocation_integrator_internal.cpp
  rk_integrator.hpp
  rk_integrator.cpp
  rk_integrator_internal.hpp
  rk_integrator_internal.cpp
  integration_tools.cpp
  integration_tools.hpp
)

if(ENABLE_STATIC)
add_library(casadi_integration STATIC ${INTEGRATION_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_integration SHARED ${INTEGRATION_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_integration
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/integration
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
