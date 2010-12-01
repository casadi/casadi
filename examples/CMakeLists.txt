find_package(IPOPT QUIET)
find_package(SUNDIALS QUIET)

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

# Rocket solved with multiple shooting
add_executable(rocket_multiple_shooting rocket_multiple_shooting.cpp)
target_link_libraries(rocket_multiple_shooting ipopt_interface casadi ${IPOPT_LIBRARIES})

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
ENDIF(FOUND_SUNDIALS AND FOUND_IPOPT)


