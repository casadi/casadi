# Add IPOPT link directories
if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Exaple 10.1 from Biegler
if(IPOPT_FOUND)
  add_executable(biegler_10_1 biegler_10_1.cpp)
  target_link_libraries(biegler_10_1
    casadi_ipopt_interface casadi
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})
endif()

# Exaple 10.2 from Biegler
if(IPOPT_FOUND)
  add_executable(biegler_10_2 biegler_10_2.cpp)
  target_link_libraries(biegler_10_2
    casadi_ipopt_interface casadi
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})
endif()

# Simple OCP example
if(IPOPT_FOUND)
  add_executable(simple_ocp_sequential simple_ocp_sequential.cpp)
  target_link_libraries(simple_ocp_sequential
    casadi_ipopt_interface casadi
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES} ${EXTRA_LIBS})
endif()

if(IPOPT_FOUND)
  add_executable(simple_ocp_simultaneous simple_ocp_simultaneous.cpp)
  target_link_libraries(simple_ocp_simultaneous
    casadi_ipopt_interface casadi
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES} ${EXTRA_LIBS})
endif()

# Exaple 10.1 from Biegler using KNITRO
if(KNITRO_FOUND)
  add_executable(biegler_10_1_knitro biegler_10_1_knitro.cpp)
  target_link_libraries(biegler_10_1_knitro
    casadi_knitro_interface casadi
    ${KNITRO_LIBRARIES} ${CASADI_DEPENDENCIES})
endif()

if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(convection_diffusion convection_diffusion.cpp)
  target_link_libraries(convection_diffusion
    casadi_ipopt_interface casadi_sundials_interface casadi_csparse_interface casadi_optimal_control casadi
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${TINYXML_LIBRARIES} ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES})
endif()

# Det minor from ADOL-C
if(WITH_CSPARSE)
  add_executable(symbolic_fact symbolic_fact.cpp)
  target_link_libraries(symbolic_fact
    casadi_csparse_interface casadi
    ${CASADI_DEPENDENCIES} ${CSPARSE_LIBRARIES})
endif()

# Det minor from ADOL-C
add_executable(det_minor det_minor.cpp)
target_link_libraries(det_minor casadi ${CASADI_DEPENDENCIES})

# Det minor example from CasADi papper
add_executable(det_minor_paper det_minor_paper.cpp)
target_link_libraries(det_minor_paper casadi ${CASADI_DEPENDENCIES})

if(IPOPT_FOUND)
  if(OOQP_FOUND)
    add_executable(test_ooqp test_ooqp.cpp)
    target_link_libraries(test_ooqp casadi_ooqp_interface casadi ${OOQP_LIBRARIES} ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})
  endif()
endif()

add_executable(issue_367 issue_367.cpp)
target_link_libraries(issue_367 casadi ${CASADI_DEPENDENCIES})

if(WITH_LLVM)
  add_subdirectory(llvm)
endif()

# single shooting using a RK4 integrator and SQP
if(IPOPT_FOUND AND WITH_QPOASES)
  add_executable(vdp_lifted_newton vdp_lifted_newton.cpp)
  target_link_libraries(vdp_lifted_newton casadi_integration casadi_nonlinear_programming casadi_qpoases_interface casadi_ipopt_interface casadi ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES} ${QPOASES_LIBRARIES})

  add_executable(albersmeyer_toy_example albersmeyer_toy_example.cpp)
  target_link_libraries(albersmeyer_toy_example casadi_integration casadi_nonlinear_programming casadi_qpoases_interface casadi_ipopt_interface casadi ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES} ${QPOASES_LIBRARIES})
endif()

if(WITH_QPOASES AND WITH_IPOPT)
  add_executable(shallow_water shallow_water.cpp)
  target_link_libraries(shallow_water casadi_core)
endif()

if(WITH_QPOASES AND WITH_IPOPT)
  add_executable(nmpc_example nmpc_example.cpp)
  target_link_libraries(nmpc_example casadi_qpoases_interface casadi_nonlinear_programming casadi_ipopt_interface casadi ${QPOASES_LIBRARIES} ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES})
endif()

if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(issue_934 issue_934.cpp)
  target_link_libraries(issue_934 casadi_integration casadi_sundials_interface casadi_csparse_interface casadi ${CSPARSE_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES})
endif()
