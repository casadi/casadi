cmake_minimum_required(VERSION 2.8.6)

# Simple library that shows how CasADi can be extended
add_library(extending_casadi extending_casadi.hpp extending_casadi.cpp)
target_link_libraries(extending_casadi casadi)

