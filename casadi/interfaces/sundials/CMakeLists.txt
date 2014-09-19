cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

# CVODES
casadi_plugin(Integrator cvodes
  sundials_internal.hpp  sundials_internal.cpp
  cvodes_internal.hpp    cvodes_internal.cpp
  cvodes_internal_meta.cpp
  )
casadi_plugin_link_libraries(Integrator cvodes ${SUNDIALS_LIBRARIES})

# IDAS
casadi_plugin(Integrator idas
  idas_internal.hpp      idas_internal.cpp
  idas_internal_meta.cpp
  )
casadi_plugin_link_libraries(Integrator idas casadi_integrator_cvodes)

# KINSOL
casadi_plugin(ImplicitFunction kinsol 
   kinsol_internal.hpp    kinsol_internal.cpp
   kinsol_internal_meta.cpp
   )
casadi_plugin_link_libraries(ImplicitFunction kinsol  ${SUNDIALS_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sundials
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

