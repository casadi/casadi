cmake_minimum_required(VERSION 3.21)
project(alpaqa-examples)

if (PROJECT_IS_TOP_LEVEL)
    find_package(alpaqa REQUIRED CONFIG)
endif()
    
find_package(casadi)
if (TARGET casadi)
    add_subdirectory(CasADi)
endif()

if (TARGET alpaqa::cutest-loader)
    add_subdirectory(CUTEst)
endif()

add_subdirectory(Presentation)
