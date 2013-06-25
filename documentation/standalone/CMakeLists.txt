cmake_minimum_required(VERSION 2.6)
project(casadi-standalone-example CXX)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../../cmake_modules)

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
endif(OPENCL_FOUND)

add_executable(nlp_codegen nlp_codegen.cpp)
target_link_libraries(nlp_codegen ${CASADI_LIBRARIES} ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES})

