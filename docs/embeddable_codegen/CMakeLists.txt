cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
project(casadi-standalone-example CXX)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

find_package(CASADI REQUIRED)
include_directories(${CASADI_INCLUDE_DIR})
find_package(IPOPT REQUIRED)
link_directories(${IPOPT_LIBRARY_DIRS})

# Possible CasADi core dependencies (depending on options)
set(CASADI_DEPENDENCIES ${CMAKE_DL_LIBS})
find_package(OpenCL)
if(OPENCL_FOUND)
  add_definitions(-DWITH_OPENCL)
  set(CASADI_DEPENDENCIES ${CASADI_DEPENDENCIES} ${OPENCL_LIBRARIES})
endif()

add_executable(embedding_part1 embedding_part1.cpp)
target_link_libraries(embedding_part1 ${CASADI_LIBRARIES} ${CASADI_DEPENDENCIES})

add_executable(embedding_part3 embedding_part3.cpp)
target_link_libraries(embedding_part3 ${CASADI_LIBRARIES} ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES})

