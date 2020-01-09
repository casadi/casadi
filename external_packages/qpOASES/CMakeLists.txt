# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012

cmake_minimum_required(VERSION 2.8.6)

project(qpOASES CXX)

# ifort on Windows
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" AND CMAKE_HOST_WIN32)
	add_definitions(-Ddgemm_=DGEMM)
	add_definitions(-Ddpotrf_=DPOTRF)
	add_definitions(-Ddtrtrs_=DTRTRS)
	add_definitions(-Ddtrcon_=DTRCON)
endif()

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
target_link_libraries(casadi_qpoases ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Uncomment to bypass user-defined linear solver (requires HSL)
#target_link_libraries(casadi_qpoases ${HSL_LIBRARIES})
