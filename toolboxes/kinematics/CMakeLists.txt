cmake_minimum_required(VERSION 2.0)

include_directories(../..)

add_library(kinematics
  frame.hpp
  frame.cpp
  frame_node.hpp
  frame_node.cpp
  kinvec.hpp
  kinvec.cpp
  kinetics.hpp
  kinetics.cpp
)

add_subdirectory(test)
add_subdirectory(examples)