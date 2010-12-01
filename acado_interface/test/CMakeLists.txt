cmake_minimum_required(VERSION 2.0)
include_directories(../../)
#SET(CMAKE_CXX_FLAGS "-O3 -std=c++0x -DWITH_CPP0X")
#SET(CMAKE_CXX_FLAGS "-O3 -fopenmp")
SET(CMAKE_CXX_FLAGS "-O3")


find_package(ACADO REQUIRED)
INCLUDE_DIRECTORIES(${ACADO_INCLUDE_DIR})

# pendulum example
add_executable(pendulum pendulum.cpp)
target_link_libraries(pendulum
  acado_interface
  casadi
  ${ACADO_LIBRARIES}
)

# rocket example
add_executable(rocket_acado rocket_acado.cpp)
target_link_libraries(rocket_acado
  acado_interface
  casadi
  ${ACADO_LIBRARIES}
)




