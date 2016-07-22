# Joel Andersson, 2016
cmake_minimum_required(VERSION 3.0.0)
project(blockSQP CXX)

set(BLOCKSQP_SRCS
  blocksqp.hpp
  blocksqp.cpp
)

casadi_external_library(casadi_blocksqp ${BLOCKSQP_SRCS})
