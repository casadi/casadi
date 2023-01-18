#include <alpaqa/dl/dl-problem.hpp>
#include <alpaqa/structured-panoc-alm.hpp>

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

// Double precision, same as in C
USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);

int main(int argc, char *argv[]) {
    // Find the problem to load
    fs::path so_name = DLPROBLEM_DLL;
    if (argc > 1)
        so_name = fs::path(argv[1]);
    else if (argc > 0)
        so_name = fs::canonical(fs::path(argv[0])).parent_path() / so_name;
    std::cout << "Loading " << so_name << std::endl;

    // Instantiate a problem
    alpaqa::dl::DLProblem problem{so_name.string()};
    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Specify the bounds
    length_t n           = problem.get_n();
    length_t m           = problem.get_m();
    vec b                = vec::Constant(m, -1);
    problem.C.lowerbound = vec::Constant(n, -alpaqa::inf<config_t>);
    problem.C.upperbound = vec::Constant(n, +alpaqa::inf<config_t>);
    problem.D.lowerbound = vec::Constant(m, -alpaqa::inf<config_t>);
    problem.D.upperbound = b;

    // Define the solvers to use
    using Accelerator = alpaqa::StructuredLBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Accelerator>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.ε              = 1e-8; // tolerance
    almparam.δ              = 1e-8;
    almparam.Δ              = 10; // penalty update factor
    almparam.max_iter       = 20;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    Accelerator::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 2;

    // Create an ALM solver using PANOC as inner solver
    OuterSolver solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x(2);
    x << 2, 2; // decision variables
    vec y(1);
    y << 1; // Lagrange multipliers

    // Solve the problem
    auto stats = solver(counted_problem, x, y);
    // y and x have been overwritten by the solution

    // Print the results
    std::cout << '\n' << *counted_problem.evaluations << '\n';
    std::cout << "status: " << stats.status << '\n'
              << "solver: " << solver.get_name() << '\n'
              << "f = " << problem.eval_f(x) << '\n'
              << "inner iterations: " << stats.inner.iterations << '\n'
              << "outer iterations: " << stats.outer_iterations << '\n'
              << "ε = " << stats.ε << '\n'
              << "δ = " << stats.δ << '\n'
              << "elapsed time:     "
              << std::chrono::duration<double>{stats.elapsed_time}.count()
              << " s" << '\n'
              << "x = " << x.transpose() << '\n'
              << "y = " << y.transpose() << '\n'
              << "avg τ = " << (stats.inner.sum_τ / stats.inner.count_τ) << '\n'
              << "L-BFGS rejected = " << stats.inner.lbfgs_rejected << '\n'
              << "L-BFGS failures = " << stats.inner.lbfgs_failures << '\n'
              << "Line search failures = " << stats.inner.linesearch_failures
              << '\n'
              << std::endl;
}