include_directories(../../)

# Add IPOPT link directories
if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi_core)

# Small example on how sparsity can be propagated throw a CasADi expression
add_executable(propagating_sparsity propagating_sparsity.cpp)
target_link_libraries(propagating_sparsity casadi_core)

# Rocket using Ipopt
if(IPOPT_FOUND)
  add_executable(rocket_ipopt rocket_ipopt.cpp)
  target_link_libraries(rocket_ipopt casadi_nlpsolver_ipopt)
endif()

# basic snopt example
if(SNOPT_FOUND)
  add_executable(snopt_nl snopt_nl.cpp)
  target_link_libraries(snopt_nl
    casadi_nlpsolver_snopt)

  add_executable(rocket_snopt rocket_snopt.cpp)
  target_link_libraries(rocket_snopt casadi_nlpsolver_snopt)
endif()

# Rocket solved with MX + SX
if(IPOPT_FOUND)
  add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
  target_link_libraries(rocket_mx_and_sx casadi_nlpsolver_ipopt)
endif()

# Parse AMPL models and solve with Ipopt
if(IPOPT_FOUND)
  add_executable(ipopt_nl ipopt_nl.cpp)
  target_link_libraries(ipopt_nl
    casadi_nlpsolver_ipopt)

  add_executable(sqpmethod_nl sqpmethod_nl.cpp)
  target_link_libraries(sqpmethod_nl
    casadi_nlpsolver_ipopt casadi_nlpsolver_sqpmethod casadi_qpsolver_nlp)

endif()

if(SQIC_FOUND AND IPOPT_FOUND)
  add_executable(stabilized_sqpmethod_nl stabilized_sqpmethod_nl.cpp)
  target_link_libraries(stabilized_sqpmethod_nl
    casadi_stabilizedqpsolver_qp)
endif()

# Parse AMPL models and solve with Worhp
if(WORHP_FOUND AND IPOPT_FOUND)
  add_executable(worhp_nl worhp_nl.cpp)
  target_link_libraries(worhp_nl
    casadi_nlpsolver_worhp)
endif()

# Parse AMPL models and solve with KNITRO
if(KNITRO_FOUND)
  add_executable(knitro_nl knitro_nl.cpp)
  target_link_libraries(knitro_nl
    casadi_nlpsolver_knitro)
endif()

# rocket example using sundials and ipopt
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(rocket_single_shooting rocket_single_shooting.cpp)
  target_link_libraries(rocket_single_shooting
    casadi_integrator_cvodes casadi_integrator_idas casadi_nlpsolver_ipopt casadi_integrator_rk
    casadi_nlpsolver_scpgen casadi_implicitfunction_newton casadi_qpsolver_nlp)
endif()

# Van der pol with multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(vdp_multiple_shooting vdp_multiple_shooting.cpp)
  target_link_libraries(vdp_multiple_shooting
    casadi_integrator_cvodes casadi_nlpsolver_ipopt  casadi_linearsolver_csparse
    casadi_optimal_control casadi_integrator_oldcollocation
    casadi_implicitfunction_newton)
endif()

# Writing a multiple shooting code from scratch
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(multiple_shooting_from_scratch multiple_shooting_from_scratch.cpp)
  target_link_libraries(multiple_shooting_from_scratch
    casadi_integrator_cvodes casadi_nlpsolver_ipopt)
endif()

# CSTR multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(cstr cstr.cpp)
  target_link_libraries(cstr casadi_optimal_control)
endif()

# CSparse via CasADi
if(WITH_CSPARSE)
  add_executable(test_csparse_casadi test_csparse_casadi.cpp)
  target_link_libraries(test_csparse_casadi casadi_linearsolver_csparse)
endif()

# Test integrators
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(sensitivity_analysis sensitivity_analysis.cpp)
  target_link_libraries(sensitivity_analysis
    casadi_integrator_collocation casadi_integrator_oldcollocation casadi_integrator_rk casadi_integrator_cvodes casadi_integrator_idas casadi_implicitfunction_kinsol casadi_linearsolver_csparse)
endif()

# Parametric sensitivities with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(parametric_sensitivities parametric_sensitivities.cpp)
  target_link_libraries(parametric_sensitivities casadi_nlpsolver_ipopt)
endif()

# Reduced hessian calculations with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(reduced_hessian reduced_hessian.cpp)
  target_link_libraries(reduced_hessian casadi_nlpsolver_ipopt)
endif()

# Parametric NLP
if(IPOPT_FOUND AND USE_CXX11)
  add_executable(parametric_nlp parametric_nlp.cpp)
  target_link_libraries(parametric_nlp casadi_nlpsolver_ipopt)
endif()

# Test OpenCL and show all devices
if(WITH_OPENCL)
  add_executable(test_opencl test_opencl.cpp)
  target_link_libraries(test_opencl ${OPENCL_LIBRARIES})
endif()

if(WITH_DL AND IPOPT_FOUND)
  add_executable(nlp_codegen nlp_codegen.cpp)
  target_link_libraries(nlp_codegen casadi_nlpsolver_ipopt)
endif()

if(WITH_DL AND NOT WIN32)
  add_executable(codegen_usage codegen_usage.cpp)
  target_link_libraries(codegen_usage casadi_core)
endif()

# Implicit Runge-Kutta integrator from scratch
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(implicit_runge-kutta implicit_runge-kutta.cpp)
  target_link_libraries(implicit_runge-kutta
    casadi_integrator_cvodes casadi_linearsolver_csparse
    casadi_implicitfunction_newton casadi_qpsolver_nlp)
endif()
