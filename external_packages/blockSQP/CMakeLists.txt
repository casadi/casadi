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
  blocksqp_defs.hpp
  blocksqp_general_purpose.cpp
  blocksqp_general_purpose.hpp
  blocksqp_glob.cpp
  blocksqp_hess.cpp
  blocksqp_iter.cpp
  blocksqp_iterate.hpp
  blocksqp_main.cpp
  blocksqp_matrix.cpp
  blocksqp_matrix.hpp
  blocksqp_method.hpp
  blocksqp_options.cpp
  blocksqp_options.hpp
  blocksqp_problemspec.cpp
  blocksqp_problemspec.hpp
  blocksqp_qp.cpp
  blocksqp_restoration.cpp
  blocksqp_restoration.hpp
  blocksqp_stats.cpp
  blocksqp_stats.hpp
  )

add_library(blocksqp ${BLOCKSQP_SRCS})
add_executable(example1 example1.cpp)
target_link_libraries(example1 blocksqp qpOASES ma57 ${LAPACK_LIBRARIES} metis gfortran)
