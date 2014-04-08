# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012

cmake_minimum_required(VERSION 2.6)

project(qpOASES CXX)

include_directories(include)
set(QPOASES_LIBRARIES ${QPOASES_LIBRARIES})
file(GLOB SRC src/*.cpp)
file(GLOB EXTRAS_SRC src/extras/*.cpp)

if(ENABLE_STATIC)
add_library(casadi_qpoases STATIC ${SRC} ${EXTRAS_SRC})
endif()
if(ENABLE_SHARED)
add_library(casadi_qpoases SHARED ${SRC} ${EXTRAS_SRC})
endif()
target_link_libraries(casadi_qpoases ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS casadi_qpoases
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_qpoases PROPERTIES COMPILE_FLAGS "-Wno-overloaded-virtual")
endif()
