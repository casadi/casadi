cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

# CVODES
casadi_plugin(Integrator cvodes
  sundials_interface.hpp  sundials_interface.cpp
  cvodes_interface.hpp    cvodes_interface.cpp
  cvodes_interface_meta.cpp
  )
casadi_plugin_link_libraries(Integrator cvodes ${SUNDIALS_LIBRARIES})

# IDAS
casadi_plugin(Integrator idas
  idas_interface.hpp      idas_interface.cpp
  idas_interface_meta.cpp
  )
casadi_plugin_link_libraries(Integrator idas casadi_integrator_cvodes)

# KINSOL
casadi_plugin(ImplicitFunction kinsol 
   kinsol_interface.hpp    kinsol_interface.cpp
   kinsol_interface_meta.cpp
   )
casadi_plugin_link_libraries(ImplicitFunction kinsol  ${SUNDIALS_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sundials
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

