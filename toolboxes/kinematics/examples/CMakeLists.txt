cmake_minimum_required(VERSION 2.0)
include_directories(../)
include_directories(../../../)

add_executable(example_cartpendulum cartpendulum.cpp)
target_link_libraries(
  example_cartpendulum
  kinematics
  casadi 
)

# add_executable(example_spinningtop spinningtop.cpp)
# target_link_libraries(
#   example_spinningtop
#   kinematics
#   casadi 
# )