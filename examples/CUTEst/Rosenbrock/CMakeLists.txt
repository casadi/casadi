cutest_sif_problem(ROSENBR)

add_executable(cutest-rosenbrock main.cpp)
target_link_libraries(cutest-rosenbrock
    PRIVATE
        panoc-alm::cutest-loader
        panoc-alm::panoc-alm)
add_dependencies(cutest-rosenbrock CUTEst::problem-ROSENBR)