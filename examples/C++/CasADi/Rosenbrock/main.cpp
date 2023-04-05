#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>

#include <alpaqa/casadi/CasADiProblem.hpp>

#include <filesystem>
#include <iostream>
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Find the problem to load
    fs::path so_name = ROSENBROCK_FUNC_DLL;
    if (argc > 1)
        so_name = fs::path(argv[1]);
    else if (argc > 0)
        so_name = fs::canonical(fs::path(argv[0])).parent_path() / so_name;
    std::cout << "Loading " << so_name << std::endl;

    // Load the problem
    alpaqa::CasADiProblem<config_t> problem{so_name.string()};

    // Specify the bounds
    problem.C.upperbound = vec::Constant(3, alpaqa::inf<config_t>);
    problem.C.lowerbound = vec::Constant(3, -alpaqa::inf<config_t>);
    problem.D.upperbound = vec::Constant(1, 0.);
    problem.D.lowerbound = vec::Constant(1, 0.);

    // Define the solvers to use
    using Direction   = alpaqa::LBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Direction>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.tolerance             = 1e-8; // tolerance
    almparam.dual_tolerance        = 1e-8;
    almparam.penalty_update_factor = 10;
    almparam.max_iter              = 20;
    almparam.print_interval        = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    Direction::AcceleratorParams lbfgsparam;
    lbfgsparam.memory = 10;

    // Create an ALM solver using PANOC as inner solver
    OuterSolver solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x(3);
    x << 2.5, 3.0, 0.75;
    vec y(1);
    y << 1;

    // Parameter
    problem.param(0) = 100;

    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Solve the problem
    auto stats = solver(counted_problem, x, y);

    // Print the results
    std::cout << '\n' << *counted_problem.evaluations << '\n';
    vec g(problem.m);
    problem.eval_g(x, g);
    std::cout << "status: " << stats.status << '\n'
              << "x = " << x.transpose() << '\n'
              << "y = " << y.transpose() << '\n'
              << "f = " << problem.eval_f(x) << '\n'
              << "g = " << g.transpose() << '\n'
              << "ε = " << stats.ε << '\n'
              << "δ = " << stats.δ << '\n'
              << "inner: " << stats.inner.iterations << '\n'
              << "outer: " << stats.outer_iterations << '\n'
              << std::endl;
}
