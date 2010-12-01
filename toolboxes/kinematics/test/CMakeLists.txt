cmake_minimum_required(VERSION 2.0)
include_directories(../)
include_directories(../../../)

add_executable(test_frames test_frames.cpp)
add_executable(test_kinvec test_kinvec.cpp)
add_library(kiteode kiteode.cpp)
add_executable(kite kite.cpp)
#add_executable(kite_integrate kite_integrate.cpp)
add_executable(unittest_kinvec unittest_kinvec.cpp)

add_executable(simulator simulator.cpp)

target_link_libraries(
  test_frames
  kinematics
  casadi 
)
target_link_libraries(
  test_kinvec
  kinematics
  casadi 
)
target_link_libraries(
  kiteode
  kinematics
  casadi 
)

target_link_libraries(
  kite
  kinematics
  casadi 
)
target_link_libraries(
  unittest_kinvec
  kinematics
  casadi 
)
target_link_libraries(
  simulator
  kiteode
  kinematics
  casadi 
)

#target_link_libraries(
#  kite_integrate
 #kinematics
 # casadi 
#  kiteode
#)
