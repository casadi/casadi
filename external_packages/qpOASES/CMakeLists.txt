# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# KU Leuven 2012

cmake_minimum_required(VERSION 3.10.2)

project(qpOASES CXX)

if(WITH_NO_QPOASES_BANNER)
  add_definitions(-D__NO_COPYRIGHT__)
endif()

# Uncomment to bypass user-defined linear solver (requires HSL)
#add_definitions(-DSOLVER_MA27)

include_directories(include)
set(QPOASES_LIBRARIES ${QPOASES_LIBRARIES})
file(GLOB SRC src/*.cpp)
file(GLOB EXTRAS_SRC src/extras/*.cpp)

add_library(casadi_qpoases STATIC ${SRC} ${EXTRAS_SRC})
target_link_libraries(casadi_qpoases lapack ${BLAS_LIBRARIES})

# Uncomment to bypass user-defined linear solver (requires HSL)
#target_link_libraries(casadi_qpoases ${HSL_LIBRARIES})
