cmake_minimum_required(VERSION 2.0)
include_directories(../../)
#SET(CMAKE_CXX_FLAGS "-O3 -std=c++0x -DWITH_CPP0X")
#SET(CMAKE_CXX_FLAGS "-O3 -fopenmp")
# SET(CMAKE_CXX_FLAGS "-O3")


find_package(IPOPT QUIET)
find_package(SUNDIALS REQUIRED QUIET)
find_package(MPI QUIET)
INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIR})

# continuously stirred tank reactor with dynamic programming
add_executable(rocket_sundials rocket_sundials.cpp)
target_link_libraries(rocket_sundials
  sundials_interface
  casadi
  ${SUNDIALS_LIBRARIES}
)

# rocket example using sundials and ipopt
IF(FOUND_IPOPT)
add_executable(rocket_ipopt_sundials rocket_ipopt_sundials.cpp)
target_link_libraries(rocket_ipopt_sundials
  sundials_interface
  ipopt_interface
  casadi
  ${MPI_LIBRARIES}
  ${IPOPT_LIBRARIES}
  ${SUNDIALS_LIBRARIES}
)
ENDIF(FOUND_IPOPT)


# add_executable(convection_diffusion convection_diffusion.cpp convection_diffusion_model.cpp convection_diffusion_model.hpp convection_diffusion_parameters.hpp)
# target_link_libraries(convection_diffusion
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES}
# )

# add_executable(convection_diffusion_simulate convection_diffusion_simulate.cpp convection_diffusion_model.cpp convection_diffusion_model.hpp convection_diffusion_parameters.hpp)
# target_link_libraries(convection_diffusion_simulate
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES}
# )

# add_executable(convection_diffusion_tabulate convection_diffusion_tabulate.cpp convection_diffusion_model.cpp convection_diffusion_model.hpp convection_diffusion_parameters.hpp)
# target_link_libraries(convection_diffusion_tabulate
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES}
# )

# add_executable(convection_diffusion_simulate_optimal convection_diffusion_simulate_optimal.cpp convection_diffusion_model.cpp convection_diffusion_model.hpp convection_diffusion_parameters.hpp)
# target_link_libraries(convection_diffusion_simulate_optimal
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES}
# )

# add_executable(van_der_pol_sundials van_der_pol_sundials.cpp)
# target_link_libraries(van_der_pol_sundials
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES})

# add_executable(chrom chrom.cpp)
# target_link_libraries(chrom
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES})

# add_executable(distillation_column distillation_column.cpp)
# target_link_libraries(distillation_column
#   integrator ocp 
#   casadi
#   ${SUNDIALS_LIBRARIES})


#XML AND SUNDIALS

# add_executable(xml_sundials xml_sundials.cpp)
# target_link_libraries(xml_sundials
#   xml
#   integrator ocp
#   casadi 
#   ${SUNDIALS_LIBRARIES}
#   tinyxml)

# continuously stirred tank reactor with dynamic programming
# add_executable(cstr_dp cstr_dp.cpp)
# target_link_libraries(cstr_dp
#   sundials_interface
#   casadi
# )
# 
# # van der pol with dynamic programming
# add_executable(van_der_pol_dp van_der_pol_dp.cpp)
# target_link_libraries(van_der_pol_dp
#   sundials_interface
#   casadi
# )
# 
# # storage optimization with dynamic programming
# add_executable(storage_dp storage_dp.cpp)
# target_link_libraries(storage_dp
#   sundials_interface
#   casadi
# )

