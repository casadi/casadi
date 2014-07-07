cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

# Common for CVODES and IDAS
casadi_library(casadi_sundials_common 
  sundials_internal.hpp  sundials_internal.cpp  sundials_integrator.hpp  sundials_integrator.cpp)
target_link_libraries(casadi_sundials_common ${SUNDIALS_LIBRARIES})

# CVODES
casadi_library(casadi_integrator_cvodes
  cvodes_internal.hpp    cvodes_internal.cpp    cvodes_integrator.hpp    cvodes_integrator.cpp)
target_link_libraries(casadi_integrator_cvodes casadi_sundials_common)

# IDAS
casadi_library(casadi_integrator_idas
  idas_internal.hpp      idas_internal.cpp      idas_integrator.hpp      idas_integrator.cpp)
target_link_libraries(casadi_integrator_idas casadi_sundials_common)

# KINSOL
casadi_library(casadi_implicitfunction_kinsol 
   kinsol_internal.hpp    kinsol_internal.cpp    kinsol_solver.hpp        kinsol_solver.cpp)
target_link_libraries(casadi_implicitfunction_kinsol ${SUNDIALS_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sundials
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

