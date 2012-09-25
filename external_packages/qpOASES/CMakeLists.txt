# This is a alternative version of the original CMakeLists.txt
# Joris Gillis
# K.U. Leuven 2012


cmake_minimum_required(VERSION 2.0)

PROJECT(qpOASES CXX)

INCLUDE_DIRECTORIES(include)

FILE(GLOB SRC src/*.cpp)
FILE(GLOB EXTRAS_SRC src/extras/*.cpp)

if(ENABLE_STATIC)
add_library(qpoases STATIC ${SRC} ${EXTRAS_SRC})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(qpoases SHARED ${SRC} ${EXTRAS_SRC})
endif(ENABLE_SHARED)
install(TARGETS qpoases
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
