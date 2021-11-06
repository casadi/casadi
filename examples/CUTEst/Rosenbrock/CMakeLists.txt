cutest_sif_problem(ROSENBR)

add_executable(cutest-rosenbrock main.cpp)
target_link_libraries(cutest-rosenbrock
    PRIVATE
        alpaqa::cutest-loader
        alpaqa::alpaqa)
add_dependencies(cutest-rosenbrock CUTEst::problem-ROSENBR)