cmake_minimum_required(VERSION 3.24)
project(alpaqa_test_package CXX)

find_package(alpaqa REQUIRED)
add_executable(example src/example.cpp)
target_link_libraries(example PRIVATE alpaqa::alpaqa)
install(TARGETS example)
