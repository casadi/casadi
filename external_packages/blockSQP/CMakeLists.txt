# Joel Andersson, 2016
cmake_minimum_required(VERSION 3.0.0)
project(blockSQP CXX Fortran)

#enable_language(Fortran REQUIRED)

find_package(LAPACK)
include_directories(.)
include_directories(/Users/jaeandersson/local/include/)
link_directories(/Users/jaeandersson/local/lib/)
link_directories(/usr/local/lib/)

set(BLOCKSQP_SRCS
  blocksqp.hpp
  blocksqp_general_purpose.cpp
  blocksqp_glob.cpp
  blocksqp_hess.cpp
  blocksqp_iter.cpp
  blocksqp_main.cpp
  blocksqp_matrix.cpp
  blocksqp_options.cpp
  blocksqp_problemspec.cpp
  blocksqp_qp.cpp
  blocksqp_restoration.cpp
  blocksqp_stats.cpp
  )

add_library(blocksqp ${BLOCKSQP_SRCS})
add_executable(example1 example1.cpp)
target_link_libraries(example1 blocksqp qpOASES ma57 ${LAPACK_LIBRARIES} metis gfortran)
