# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012

cmake_minimum_required(VERSION 2.8.6)

project(qpOASES CXX)

# Uncomment to bypass user-defined linear solver (requires HSL)
#add_definitions(-DSOLVER_MA27)

# cmake_policy use MACOSX_RPATH=ON by default
if(COMMAND cmake_policy)
  if(CMAKE_MAJOR_VERSION GREATER 2)
    cmake_policy(SET CMP0042 NEW)
  endif()
endif()

include_directories(include)
set(QPOASES_LIBRARIES ${QPOASES_LIBRARIES})
file(GLOB SRC src/*.cpp)
file(GLOB EXTRAS_SRC src/extras/*.cpp)

casadi_external_library(casadi_qpoases ${SRC} ${EXTRAS_SRC})
target_link_libraries(casadi_qpoases ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Uncomment to bypass user-defined linear solver (requires HSL)
#target_link_libraries(casadi_qpoases ${HSL_LIBRARIES})

install(DIRECTORY include
  DESTINATION include/casadi/external_packages/qpOASES
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)
