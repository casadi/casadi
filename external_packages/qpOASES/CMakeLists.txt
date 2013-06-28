# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012


cmake_minimum_required(VERSION 2.6)

PROJECT(qpOASES CXX)

INCLUDE_DIRECTORIES(include)
SET(QPOASES_LIBRARIES ${QPOASES_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
FILE(GLOB SRC src/*.cpp)
FILE(GLOB EXTRAS_SRC src/extras/*.cpp)

if(ENABLE_STATIC)
add_library(casadi_qpoases STATIC ${SRC} ${EXTRAS_SRC})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_qpoases SHARED ${SRC} ${EXTRAS_SRC})
endif(ENABLE_SHARED)
install(TARGETS casadi_qpoases
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpoases PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
