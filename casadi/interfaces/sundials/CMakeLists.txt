cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SUNDIALS_INCLUDE_DIR})

if(PRECOND_TEST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPRECOND_TEST")
endif()

set(SUNDIALS_INTERFACE_SRCS
  sundials_internal.hpp  sundials_internal.cpp  sundials_integrator.hpp  sundials_integrator.cpp
  cvodes_internal.hpp    cvodes_internal.cpp    cvodes_integrator.hpp    cvodes_integrator.cpp
  idas_internal.hpp      idas_internal.cpp      idas_integrator.hpp      idas_integrator.cpp
  kinsol_internal.hpp    kinsol_internal.cpp    kinsol_solver.hpp        kinsol_solver.cpp)

casadi_library(casadi_sundials_interface ${SUNDIALS_INTERFACE_SRCS})
target_link_libraries(casadi_sundials_interface ${SUNDIALS_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/sundials
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

