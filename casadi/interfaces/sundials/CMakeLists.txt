cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

# CVODES
casadi_library(casadi_integrator_cvodes
  sundials_internal.hpp  sundials_internal.cpp
  cvodes_internal.hpp    cvodes_internal.cpp)
target_link_libraries(casadi_integrator_cvodes ${SUNDIALS_LIBRARIES})

# IDAS
casadi_library(casadi_integrator_idas
  idas_internal.hpp      idas_internal.cpp)
target_link_libraries(casadi_integrator_idas casadi_integrator_cvodes)

# KINSOL
casadi_library(casadi_implicitfunction_kinsol 
   kinsol_internal.hpp    kinsol_internal.cpp)
target_link_libraries(casadi_implicitfunction_kinsol ${SUNDIALS_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sundials
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

