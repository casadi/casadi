cmake_minimum_required(VERSION 3.21)
project(alpaqa-examples C CXX)

if (PROJECT_IS_TOP_LEVEL)
    find_package(alpaqa REQUIRED CONFIG)
endif()
    
add_subdirectory(C++)
