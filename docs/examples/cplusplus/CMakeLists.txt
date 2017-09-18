include_directories(../../)

# Add IPOPT link directories
if(WITH_IPOPT)
  link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi)

# Callback functions
add_executable(callback callback.cpp)
target_link_libraries(callback casadi)

# Small example on how sparsity can be propagated throw a CasADi expression
add_executable(propagating_sparsity propagating_sparsity.cpp)
target_link_libraries(propagating_sparsity casadi)

# Rosenbrock problem
if(WITH_IPOPT)
  add_executable(rosenbrock rosenbrock.cpp)
  target_link_libraries(rosenbrock casadi)
endif()

# Rocket using Ipopt
if(WITH_IPOPT)
  add_executable(rocket_ipopt rocket_ipopt.cpp)
  target_link_libraries(rocket_ipopt casadi)
endif()

# Race car using Ipopt
if(WITH_IPOPT)
  add_executable(race_car race_car.cpp)
  target_link_libraries(race_car casadi)
endif()

# basic snopt example
if(WITH_SNOPT)
  add_executable(snopt_nl snopt_nl.cpp)
  target_link_libraries(snopt_nl casadi)

  add_executable(rocket_snopt rocket_snopt.cpp)
  target_link_libraries(rocket_snopt casadi)
endif()

# Rocket solved with MX + SX
if(WITH_IPOPT)
  add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
  target_link_libraries(rocket_mx_and_sx casadi)
endif()

# Parse AMPL models and solve with Ipopt
if(WITH_IPOPT)
  add_executable(ipopt_nl ipopt_nl.cpp)
  target_link_libraries(ipopt_nl casadi)

  add_executable(sqpmethod_nl sqpmethod_nl.cpp)
  target_link_libraries(sqpmethod_nl casadi)
endif()

# Parse AMPL models and solve with Worhp
if(WITH_WORHP AND WITH_IPOPT)
  add_executable(worhp_nl worhp_nl.cpp)
  target_link_libraries(worhp_nl casadi)
endif()

# Parse AMPL models and solve with KNITRO
if(WITH_KNITRO)
  add_executable(knitro_nl knitro_nl.cpp)
  target_link_libraries(knitro_nl casadi)
endif()

# rocket example using sundials and ipopt, also demonstrating linking with plugins
if(WITH_SUNDIALS AND WITH_IPOPT)
  add_executable(rocket_single_shooting rocket_single_shooting.cpp)
  target_link_libraries(rocket_single_shooting
    casadi_integrator_cvodes
    casadi_integrator_idas
    casadi_nlpsol_ipopt
    casadi_integrator_rk
    casadi_nlpsol_scpgen)
endif()

# Writing a multiple shooting code from scratch
if(WITH_SUNDIALS AND WITH_IPOPT)
  add_executable(multiple_shooting_from_scratch multiple_shooting_from_scratch.cpp)
  target_link_libraries(multiple_shooting_from_scratch casadi)
endif()

# Solve linear system of equations
add_executable(test_linsol test_linsol.cpp)
target_link_libraries(test_linsol casadi)

# Test integrators
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(sensitivity_analysis sensitivity_analysis.cpp)
  target_link_libraries(sensitivity_analysis casadi)
endif()

# Parametric NLP
if(WITH_IPOPT)
  add_executable(parametric_nlp parametric_nlp.cpp)
  target_link_libraries(parametric_nlp casadi)
endif()

# Test OpenCL and show all devices
if(WITH_OPENCL)
  add_executable(test_opencl test_opencl.cpp)
  target_link_libraries(test_opencl ${OPENCL_LIBRARIES})
endif()

if(WITH_DL AND WITH_IPOPT)
  add_executable(nlp_codegen nlp_codegen.cpp)
  target_link_libraries(nlp_codegen casadi)
endif()

if(WITH_DL AND NOT WIN32)
  add_executable(codegen_usage codegen_usage.cpp)
  target_link_libraries(codegen_usage casadi)
endif()

# Implicit Runge-Kutta integrator from scratch
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(implicit_runge-kutta implicit_runge-kutta.cpp)
  target_link_libraries(implicit_runge-kutta casadi)
endif()

if (WITH_QPOASES)
  add_executable(chain_qp chain_qp.cpp)
  target_link_libraries(chain_qp casadi)

  add_executable(lowlevel_qp lowlevel_qp.cpp)
  target_link_libraries(lowlevel_qp casadi)
endif()

if(WITH_BLOCKSQP)
  add_executable(test_blocksqp test_blocksqp.cpp)
  target_link_libraries(test_blocksqp casadi)
endif()

# DaeBuilder
add_executable(daebuilder daebuilder.cpp)
target_link_libraries(daebuilder casadi)
