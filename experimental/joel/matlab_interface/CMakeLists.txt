cmake_minimum_required(VERSION 2.0)
project(casadi_matlab_interface)

include_directories(/usr/local/matlab/extern/include)

SET(CASADI_MATLAB_SRCS
  matlab_interface.h
  matlab_interface.cxx
)
add_library(casadi_matlab SHARED ${CASADI_MATLAB_SRCS})

