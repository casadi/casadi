cmake_minimum_required(VERSION 2.0)

include_directories(../..)

# Calling SuperLU directly
add_executable(test_superlu test_superlu.c)
target_link_libraries(test_superlu superlu blas)

# SuperLU via CasADi
add_executable(test_superlu_casadi test_superlu_casadi.cpp)
target_link_libraries(test_superlu_casadi superlu_interface casadi superlu blas)

