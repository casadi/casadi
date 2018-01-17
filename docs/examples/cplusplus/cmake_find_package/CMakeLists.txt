cmake_minimum_required(VERSION 2.8.6)

# export CMAKE_PREFIX_PATH=../casadi_install/casadi/cmake/
find_package(casadi)

add_executable(casadi_demo casadi_demo.cpp)
target_link_libraries(casadi_demo casadi)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
