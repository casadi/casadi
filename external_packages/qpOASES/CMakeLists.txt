# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012

cmake_minimum_required(VERSION 2.6)

project(qpOASES CXX)

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

if(ENABLE_STATIC)
add_library(casadi_qpoases STATIC ${SRC} ${EXTRAS_SRC})
endif()
if(ENABLE_SHARED)
add_library(casadi_qpoases SHARED ${SRC} ${EXTRAS_SRC})
endif()
target_link_libraries(casadi_qpoases ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS casadi_qpoases
  LIBRARY DESTINATION ${PLUGIN_INSTALL_PATH}
  ARCHIVE DESTINATION ${PLUGIN_INSTALL_PATH}
)
