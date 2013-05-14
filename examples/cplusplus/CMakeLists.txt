include_directories(../../)

# Add IPOPT link directories
if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi ${CASADI_DEPENDENCIES})

# Determinant by minor expansion using MX
add_executable(det_minor_mx det_minor_mx.cpp)
target_link_libraries(det_minor_mx casadi ${CASADI_DEPENDENCIES})

# Small example on how sparsity can be propagated throw a CasADi expression
add_executable(propagating_sparsity propagating_sparsity.cpp)
target_link_libraries(propagating_sparsity casadi ${CASADI_DEPENDENCIES})

# Rocket using Ipopt
if(IPOPT_FOUND)
  add_executable(rocket_ipopt rocket_ipopt.cpp)
  target_link_libraries(rocket_ipopt 
    casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()


# Rocket solved with MX + SX
if(IPOPT_FOUND)
  add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
  target_link_libraries(rocket_mx_and_sx 
    casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Parse AMPL models and solve with Ipopt
if(IPOPT_FOUND)
  add_executable(ipopt_nl ipopt_nl.cpp)
  target_link_libraries(ipopt_nl 
    casadi_ipopt_interface casadi_nonlinear_programming casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )

  add_executable(sqpmethod_nl sqpmethod_nl.cpp)
  target_link_libraries(sqpmethod_nl 
    casadi_ipopt_interface casadi_nonlinear_programming casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )

endif()

# Parse AMPL models and solve with Worhp
if(WORHP_FOUND AND IPOPT_FOUND)
  add_executable(worhp_nl worhp_nl.cpp)
  target_link_libraries(worhp_nl 
    casadi_worhp_interface casadi casadi_nonlinear_programming
    ${WORHP_LIBRARIES} ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES} # NOTE!!! ADDED IPOPT due to missing symbols
  )
endif()

# Parse AMPL models and solve with KNITRO
if(KNITRO_FOUND)
  add_executable(knitro_nl knitro_nl.cpp)
  target_link_libraries(knitro_nl 
    casadi_knitro_interface casadi casadi_nonlinear_programming
    ${KNITRO_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# rocket example using sundials and ipopt
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(rocket_single_shooting rocket_single_shooting.cpp)
  target_link_libraries(rocket_single_shooting 
    casadi_sundials_interface casadi_ipopt_interface casadi_integration casadi_nonlinear_programming casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Van der pol with multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(vdp_multiple_shooting vdp_multiple_shooting.cpp)
  target_link_libraries(vdp_multiple_shooting   
    casadi_sundials_interface casadi_ipopt_interface  casadi_csparse_interface
    casadi_optimal_control casadi_integration casadi_nonlinear_programming casadi
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${TINYXML_LIBRARIES} ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# Writing a multiple shooting code from scratch
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(multiple_shooting_from_scratch multiple_shooting_from_scratch.cpp)
  target_link_libraries(multiple_shooting_from_scratch 
    casadi_sundials_interface casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES} )
endif()

# CSTR multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(cstr cstr.cpp)
  target_link_libraries(cstr 
    casadi_sundials_interface casadi_ipopt_interface casadi_optimal_control casadi_csparse_interface casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${TINYXML_LIBRARIES} ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# CSparse via CasADi
if(WITH_CSPARSE)
  add_executable(test_csparse_casadi test_csparse_casadi.cpp)
  target_link_libraries(test_csparse_casadi 
    casadi_csparse_interface casadi 
    ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# Test integrators
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(sensitivity_analysis sensitivity_analysis.cpp)
  target_link_libraries(sensitivity_analysis
    casadi_integration casadi_sundials_interface casadi_csparse_interface casadi
    ${SUNDIALS_LIBRARIES} ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

if(LIFTOPT_FOUND)
  add_executable(test_liftopt test_liftopt.cpp)
  target_link_libraries(test_liftopt 
    casadi_liftopt_interface casadi 
    ${LIFTOPT_LIBRARIES} ${CASADI_DEPENDENCIES} )
endif()

# Parametric sensitivities with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(parametric_sensitivities parametric_sensitivities.cpp)
  target_link_libraries(parametric_sensitivities 
    casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Reduced hessian calculations with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(reduced_hessian reduced_hessian.cpp)
  target_link_libraries(reduced_hessian
    casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Parametric NLP
if(IPOPT_FOUND AND USE_CXX11)
  add_executable(parametric_nlp parametric_nlp.cpp)
  target_link_libraries(parametric_nlp
    casadi_ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Test OpenCL and show all devices
if(WITH_OPENCL)
  add_executable(test_opencl test_opencl.cpp)
  target_link_libraries(test_opencl ${OPENCL_LIBRARIES})
endif()

if(WITH_DL AND IPOPT_FOUND)
  add_executable(nlp_codegen nlp_codegen.cpp)
  target_link_libraries(nlp_codegen casadi_ipopt_interface casadi ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES})
endif()

if(WITH_DL)
  add_executable(codegen_usage codegen_usage.cpp)
  target_link_libraries(codegen_usage casadi ${CASADI_DEPENDENCIES})
endif()

# Implicit Runge-Kutta integrator from scratch
if(WITH_SUNDIALS AND WITH_CSPARSE)
  add_executable(implicit_runge-kutta implicit_runge-kutta.cpp)
  target_link_libraries(implicit_runge-kutta
    casadi_nonlinear_programming casadi_integration casadi_sundials_interface casadi_csparse_interface casadi
    ${SUNDIALS_LIBRARIES} ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

