cmake_minimum_required(VERSION 2.0)

include_directories(../..)

add_executable(superlu superlu.c)
target_link_libraries(superlu superlu_casadi blas)
