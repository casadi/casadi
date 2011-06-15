find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(SUNDIALS QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(LIFTOPT QUIET)

IF(IPOPT_FOUND)
include_directories(${IPOPT_INCLUDE_DIR})
ENDIF(IPOPT_FOUND)

IF(SUNDIALS_FOUND)
INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})
ENDIF(SUNDIALS_FOUND)

IF(IPOPT_FOUND)

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

ENDIF(SUNDIALS_FOUND AND IPOPT_FOUND)

# Det minor from ADOL-C
add_executable(det_minor det_minor.cpp)
target_link_libraries(det_minor casadi)


if(LIFTOPT_FOUND)
include_directories(${LIFTOPT_INCLUDE_DIRS})
add_executable(test_liftopt test_liftopt.cpp)
target_link_libraries(test_liftopt liftopt_interface casadi ${LIFTOPT_LIBRARIES})

endif(LIFTOPT_FOUND)

add_executable(mx_jacobian mx_jacobian.cpp)
target_link_libraries(mx_jacobian casadi)



