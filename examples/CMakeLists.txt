find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)

include_directories(../)

IF(IPOPT_FOUND)
include_directories(${IPOPT_INCLUDE_DIR})
ENDIF(IPOPT_FOUND)

IF(SUNDIALS_FOUND)
INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})
ENDIF(SUNDIALS_FOUND)

# Error handling
add_executable(casadi_error_handling casadi_error_handling.cpp)
target_link_libraries(casadi_error_handling casadi)

IF(IPOPT_FOUND)

# Rocket using Ipopt
add_executable(rocket_ipopt rocket_ipopt.cpp)
target_link_libraries(rocket_ipopt ipopt_interface casadi ${IPOPT_LIBRARIES})

# Rocket solved with MX + SX
add_executable(rocket_mx_and_sx rocket_mx_and_sx.cpp)
target_link_libraries(rocket_mx_and_sx ipopt_interface casadi ${IPOPT_LIBRARIES})

# Exaple 10.1 from Biegler
add_executable(biegler_10_1 biegler_10_1.cpp)
target_link_libraries(biegler_10_1 ipopt_interface casadi ${IPOPT_LIBRARIES})

# Exaple 10.2 from Biegler
add_executable(biegler_10_2 biegler_10_2.cpp)
target_link_libraries(biegler_10_2 ipopt_interface casadi ${IPOPT_LIBRARIES})
ENDIF(IPOPT_FOUND)


IF(KNITRO_FOUND)

# Exaple 10.1 from Biegler using KNITRO
add_executable(biegler_10_1_knitro biegler_10_1_knitro.cpp)
target_link_libraries(biegler_10_1_knitro knitro_interface casadi ${KNITRO_LIBRARIES})

ENDIF(KNITRO_FOUND)

IF(SUNDIALS_FOUND AND IPOPT_FOUND)
add_executable(convection_diffusion convection_diffusion.cpp)
target_link_libraries(convection_diffusion ipopt_interface sundials_interface optimal_control casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES})

# rocket example using sundials and ipopt
add_executable(rocket_ipopt_sundials rocket_ipopt_sundials.cpp)
target_link_libraries(rocket_ipopt_sundials sundials_interface ipopt_interface casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES})

# Van der pol with multiple shooting
add_executable(vdp_multiple_shooting vdp_multiple_shooting.cpp)
target_link_libraries(vdp_multiple_shooting sundials_interface ipopt_interface casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES})

# Van der pol with multiple shooting
add_executable(cstr cstr.cpp)
target_link_libraries(cstr sundials_interface ipopt_interface optimal_control casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES} tinyxml)

ENDIF(SUNDIALS_FOUND AND IPOPT_FOUND)

# SuperLU via CasADi
if(WITH_SUPERLU)
IF(BLAS_FOUND)
add_executable(test_superlu_casadi test_superlu_casadi.cpp)
target_link_libraries(test_superlu_casadi superlu_interface casadi superlu ${BLAS_LIBRARIES})
ENDIF(BLAS_FOUND)
endif(WITH_SUPERLU)

# Det minor from ADOL-C
add_executable(det_minor det_minor.cpp)
target_link_libraries(det_minor casadi)

# Sundials - just integrate
IF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)
add_executable(sundials_testing sundials_testing.cpp)
target_link_libraries(sundials_testing
  sundials_interface
  superlu_interface 
  lapack_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  superlu
  ${BLAS_LIBRARIES}
)
ENDIF(SUNDIALS_FOUND AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)





