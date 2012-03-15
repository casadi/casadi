find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(LIFTOPT QUIET)

include_directories(../../)


IF(IPOPT_FOUND)
  include_directories(${IPOPT_INCLUDE_DIR})
ENDIF(IPOPT_FOUND)

IF(SUNDIALS_FOUND)
  INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})
ENDIF(SUNDIALS_FOUND)

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi ${CASADI_DEPENDENCIES})

# Determinant by minor expansion using MX
add_executable(det_minor_mx det_minor_mx.cpp)
target_link_libraries(det_minor_mx casadi ${CASADI_DEPENDENCIES})

IF(IPOPT_FOUND)

# Rocket using Ipopt
add_executable(rocket_ipopt rocket_ipopt.cpp)
target_link_libraries(rocket_ipopt ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})

# Rocket solved with MX + SX
add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
target_link_libraries(rocket_mx_and_sx ipopt_interface casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})

# Parse AMPL models and solve with Ipopt
add_executable(ipopt_nl ipopt_nl.cpp)
target_link_libraries(ipopt_nl ipopt_interface nonlinear_programming casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})

ENDIF(IPOPT_FOUND)


IF(SUNDIALS_FOUND AND IPOPT_FOUND)

# rocket example using sundials and ipopt
add_executable(rocket_single_shooting rocket_single_shooting.cpp)
target_link_libraries(rocket_single_shooting sundials_interface ipopt_interface integration casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES})

# Van der pol with multiple shooting
add_executable(vdp_multiple_shooting vdp_multiple_shooting.cpp)
target_link_libraries(vdp_multiple_shooting
  sundials_interface
  ipopt_interface
  optimal_control
  csparse_interface
  casadi
  ${IPOPT_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
  tinyxml
  csparse
  ${CASADI_DEPENDENCIES} 
)

# How to write a multiple shooting code from scratch
add_executable(multiple_shooting_from_scratch multiple_shooting_from_scratch.cpp)
target_link_libraries(multiple_shooting_from_scratch sundials_interface ipopt_interface casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} ${CASADI_DEPENDENCIES} )

# Van der pol with multiple shooting
add_executable(cstr cstr.cpp)
target_link_libraries(cstr 
  sundials_interface 
  ipopt_interface 
  optimal_control 
  csparse_interface
  casadi 
  ${IPOPT_LIBRARIES} 
  ${SUNDIALS_LIBRARIES} 
  tinyxml
  csparse
  ${CASADI_DEPENDENCIES} 
)


ENDIF(SUNDIALS_FOUND AND IPOPT_FOUND)

# SuperLU via CasADi
if(WITH_SUPERLU)
IF(BLAS_FOUND)
add_executable(test_superlu_casadi test_superlu_casadi.cpp)
target_link_libraries(test_superlu_casadi superlu_interface casadi superlu ${BLAS_LIBRARIES} ${CASADI_DEPENDENCIES} )
ENDIF(BLAS_FOUND)
endif(WITH_SUPERLU)

include_directories(../../external_packages/CSparse/Include)

# CSparse via CasADi
add_executable(test_csparse_casadi test_csparse_casadi.cpp)
target_link_libraries(test_csparse_casadi csparse_interface casadi csparse ${CASADI_DEPENDENCIES} )

# Sundials - just integrate
IF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)
add_executable(sundials_testing sundials_testing.cpp)
target_link_libraries(sundials_testing
  sundials_interface
  superlu_interface 
  lapack_interface
  csparse_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  superlu
  csparse
  ${BLAS_LIBRARIES}
  ${CASADI_DEPENDENCIES}
)
ENDIF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)

if(LIFTOPT_FOUND)
include_directories(${LIFTOPT_INCLUDE_DIRS})
add_executable(test_liftopt test_liftopt.cpp)
target_link_libraries(test_liftopt liftopt_interface casadi ${LIFTOPT_LIBRARIES} ${CASADI_DEPENDENCIES} )

endif(LIFTOPT_FOUND)

