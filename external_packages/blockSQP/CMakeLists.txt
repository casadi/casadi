# Joel Andersson, 2016
cmake_minimum_required(VERSION 3.0.0)
project(blockSQP CXX)

include_directories(${QPOASES_INCLUDE_DIR})

set(BLOCKSQP_SRCS
  blocksqp.hpp
  blocksqp.cpp
)

casadi_external_library(casadi_blocksqp ${BLOCKSQP_SRCS})
target_link_libraries(casadi_blocksqp ${QPOASES_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(blocksqp_example1 blocksqp_example1.cpp)
target_link_libraries(blocksqp_example1 casadi_blocksqp)
