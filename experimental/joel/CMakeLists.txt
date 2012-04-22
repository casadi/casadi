find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(LIFTOPT QUIET)
find_package(OOQP QUIET)
find_package(ACADO QUIET)

IF(IPOPT_FOUND)
include_directories(${IPOPT_INCLUDE_DIR})
ENDIF(IPOPT_FOUND)

IF(SUNDIALS_FOUND)
INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})
ENDIF(SUNDIALS_FOUND)

IF(IPOPT_FOUND)

# Exaple 10.1 from Biegler
add_executable(biegler_10_1 biegler_10_1.cpp)
target_link_libraries(biegler_10_1 ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})

# Exaple 10.2 from Biegler
add_executable(biegler_10_2 biegler_10_2.cpp)
target_link_libraries(biegler_10_2 ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})

# Simple OCP example
add_executable(simple_ocp_sequential simple_ocp_sequential.cpp)
target_link_libraries(simple_ocp_sequential ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES} ${EXTRA_LIBS})

add_executable(simple_ocp_simultaneous simple_ocp_simultaneous.cpp)
target_link_libraries(simple_ocp_simultaneous ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES} ${EXTRA_LIBS})


ENDIF(IPOPT_FOUND)


IF(KNITRO_FOUND)

# Exaple 10.1 from Biegler using KNITRO
add_executable(biegler_10_1_knitro biegler_10_1_knitro.cpp)
target_link_libraries(biegler_10_1_knitro knitro_interface casadi ${KNITRO_LIBRARIES} ${CASADI_DEPENDENCIES})

ENDIF(KNITRO_FOUND)

IF(SUNDIALS_FOUND AND IPOPT_FOUND)
add_executable(convection_diffusion convection_diffusion.cpp)
target_link_libraries(convection_diffusion
  ipopt_interface
  sundials_interface
  csparse_interface
  optimal_control
  casadi
  ${IPOPT_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
  tinyxml
  csparse 
  ${CASADI_DEPENDENCIES}
  )

ENDIF(SUNDIALS_FOUND AND IPOPT_FOUND)

# Det minor from ADOL-C
add_executable(symbolic_fact symbolic_fact.cpp)
target_link_libraries(symbolic_fact csparse_interface casadi ${CASADI_DEPENDENCIES} csparse )


# Det minor from ADOL-C
add_executable(det_minor det_minor.cpp)
target_link_libraries(det_minor casadi ${CASADI_DEPENDENCIES})

# Det minor example from CasADi papper
add_executable(det_minor_paper det_minor_paper.cpp )
target_link_libraries(det_minor_paper casadi ${CASADI_DEPENDENCIES})

if(LIFTOPT_FOUND)
include_directories(${LIFTOPT_INCLUDE_DIRS})
add_executable(test_liftopt test_liftopt.cpp)
target_link_libraries(test_liftopt liftopt_interface casadi ${CASADI_DEPENDENCIES} ${LIFTOPT_LIBRARIES})

endif(LIFTOPT_FOUND)

add_executable(mx_jacobian mx_jacobian.cpp)
target_link_libraries(mx_jacobian casadi ${CASADI_DEPENDENCIES})

IF(IPOPT_FOUND)
if(OOQP_FOUND)
add_executable(test_ooqp test_ooqp.cpp )
target_link_libraries(test_ooqp ooqp_interface casadi ${OOQP_LIBRARIES} ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})
endif(OOQP_FOUND)
endif(IPOPT_FOUND)

add_executable(ticket_309 ticket_309.cpp )
target_link_libraries(ticket_309 casadi ${CASADI_DEPENDENCIES})

if(WITH_LLVM)
add_subdirectory(llvm)
endif(WITH_LLVM)

# single shooting using a RK4 integrator and SQP
IF(IPOPT_FOUND AND ACADO_FOUND AND ACADO_SHARED_LIBRARIES)
  add_executable(single_shooting_rk4 single_shooting_rk4.cpp)
  target_link_libraries(single_shooting_rk4 integration nonlinear_programming qpoases_interface ipopt_interface casadi ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES} ${ACADO_SHARED_LIBRARIES})
ENDIF(IPOPT_FOUND AND ACADO_FOUND AND ACADO_SHARED_LIBRARIES)

