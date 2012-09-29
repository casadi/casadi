include_directories(../../)

# Add IPOPT link directories
if(IPOPT_FOUND)
  link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi ${CASADI_DEPENDENCIES})

# Working with weak references
add_executable(weak_references weak_references.cpp )
target_link_libraries(weak_references casadi ${CASADI_DEPENDENCIES})

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
    ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()


# Rocket solved with MX + SX
if(IPOPT_FOUND)
  add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
  target_link_libraries(rocket_mx_and_sx 
    ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Parse AMPL models and solve with Ipopt
if(IPOPT_FOUND)
  add_executable(ipopt_nl ipopt_nl.cpp)
  target_link_libraries(ipopt_nl 
    ipopt_interface casadi_nonlinear_programming casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# rocket example using sundials and ipopt
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(rocket_single_shooting rocket_single_shooting.cpp)
  target_link_libraries(rocket_single_shooting 
    sundials_interface ipopt_interface integration casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Van der pol with multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(vdp_multiple_shooting vdp_multiple_shooting.cpp)
  target_link_libraries(vdp_multiple_shooting
    sundials_interface ipopt_interface optimal_control csparse_interface casadi
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} tinyxml ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# Writing a multiple shooting code from scratch
if(WITH_SUNDIALS AND IPOPT_FOUND)
  add_executable(multiple_shooting_from_scratch multiple_shooting_from_scratch.cpp)
  target_link_libraries(multiple_shooting_from_scratch 
    sundials_interface ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES} )
endif()

# CSTR multiple shooting
if(WITH_SUNDIALS AND IPOPT_FOUND AND WITH_CSPARSE)
  add_executable(cstr cstr.cpp)
  target_link_libraries(cstr 
    sundials_interface ipopt_interface optimal_control csparse_interface casadi 
    ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} tinyxml ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# CSparse via CasADi
if(WITH_CSPARSE)
  add_executable(test_csparse_casadi test_csparse_casadi.cpp)
  target_link_libraries(test_csparse_casadi 
    csparse_interface casadi 
    ${CSPARSE_LIBRARIES} ${CASADI_DEPENDENCIES} 
  )
endif()

# Sundials - just integrate
if(WITH_SUNDIALS AND LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
  add_executable(sundials_testing sundials_testing.cpp)
  target_link_libraries(sundials_testing
    sundials_interface lapack_interface csparse_interface casadi
    ${SUNDIALS_LIBRARIES} ${LAPACK_LIBRARIES} ${CSPARSE_LIBRARIES} ${BLAS_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

if(LIFTOPT_FOUND)
  add_executable(test_liftopt test_liftopt.cpp)
  target_link_libraries(test_liftopt 
    liftopt_interface casadi 
    ${LIFTOPT_LIBRARIES} ${CASADI_DEPENDENCIES} )
endif()

# Parametric sensitivities with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(parametric_sensitivities parametric_sensitivities.cpp)
  target_link_libraries(parametric_sensitivities 
    ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Reduced hessian calculations with sIPOPT
if(IPOPT_FOUND AND USE_CXX11 AND WITH_SIPOPT)
  add_executable(reduced_hessian reduced_hessian.cpp)
  target_link_libraries(reduced_hessian
    ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

# Parametric NLP
if(IPOPT_FOUND AND USE_CXX11)
  add_executable(parametric_nlp parametric_nlp.cpp)
  target_link_libraries(parametric_nlp
    ipopt_interface casadi 
    ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES}
  )
endif()

