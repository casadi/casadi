find_package(IPOPT QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)

include_directories(../)

IF(FOUND_IPOPT)
include_directories(${IPOPT_INCLUDE_DIR})
ENDIF(FOUND_IPOPT)

IF(FOUND_SUNDIALS)
INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})
ENDIF(FOUND_SUNDIALS)

IF(FOUND_IPOPT)

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
ENDIF(FOUND_IPOPT)

IF(FOUND_SUNDIALS AND FOUND_IPOPT)
add_executable(convection_diffusion convection_diffusion.cpp ocp_tools.hpp ocp_tools.cpp)
target_link_libraries(convection_diffusion ipopt_interface sundials_interface casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES})

# rocket example using sundials and ipopt
add_executable(rocket_ipopt_sundials rocket_ipopt_sundials.cpp)
target_link_libraries(rocket_ipopt_sundials sundials_interface ipopt_interface casadi ${IPOPT_LIBRARIES} ${SUNDIALS_LIBRARIES})
ENDIF(FOUND_SUNDIALS AND FOUND_IPOPT)

# SuperLU via CasADi
IF(BLAS_FOUND)
add_executable(test_superlu_casadi test_superlu_casadi.cpp)
target_link_libraries(test_superlu_casadi superlu_interface casadi superlu ${BLAS_LIBRARIES})
ENDIF(BLAS_FOUND)

# Det minor from ADOL-C
add_executable(det_minor det_minor.cpp)
target_link_libraries(det_minor casadi)

# Sundials - just integrate
IF(FOUND_SUNDIALS)
add_executable(rocket_sundials rocket_sundials.cpp)
target_link_libraries(rocket_sundials
  sundials_interface
  casadi
  ${SUNDIALS_LIBRARIES}
)
ENDIF(FOUND_SUNDIALS)





