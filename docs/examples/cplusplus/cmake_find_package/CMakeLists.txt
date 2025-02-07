cmake_minimum_required(VERSION 3.10.2)

project(demo)

# export CMAKE_PREFIX_PATH=../casadi_install/casadi/cmake/
find_package(casadi CONFIG REQUIRED)

add_executable(casadi_demo casadi_demo.cpp)
target_link_libraries(casadi_demo casadi::casadi)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

install(TARGETS casadi_demo RUNTIME DESTINATION . LIBRARY DESTINATION .)
