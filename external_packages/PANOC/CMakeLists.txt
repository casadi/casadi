# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2018

cmake_minimum_required(VERSION 2.8.6)

project(PANOS C)

include_directories(include)
file(GLOB panoc_SRC
    "PANOC/*.c"
)
file(GLOB PANOC/constraint_SRC
    "PANOC/constraints/*.c"
)
add_library(casadi_panoc STATIC ${panoc_SRC} ${constraint_SRC})

