cmake_minimum_required(VERSION 3.21...4.1)
project(alpaqa-cmake-examples CXX)

# Locate the alpaqa package
find_package(alpaqa 1.0.0 REQUIRED)

# Compile the example code into an executable
add_executable(alpaqa-qp-solver
    "problem.cpp"
    "alpaqa-qp-solver.cpp")
# Indicate that this executable depends on the alpaqa library
target_link_libraries(alpaqa-qp-solver PRIVATE alpaqa::alpaqa)

# Enable compiler warnings
include("warnings.cmake")
add_warnings_target(warnings)
target_link_libraries(alpaqa-qp-solver PRIVATE warnings)
