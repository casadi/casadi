#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>
#include <alpaqa/problem/unconstr-problem.hpp>

#include <iostream>

// Problem specification
// minimize  (a - x)² + b(y - x²)²
struct RosenbrockProblem : alpaqa::UnconstrProblem<alpaqa::DefaultConfig> {
    // Problem parameters
    real_t a = 2, b = 100;

    // Number of unknowns
    length_t get_n() const { return 2; }

    // Cost
    real_t eval_f(crvec xy) const {
        auto x = xy(0), y = xy(1);
        return sq(a - x) + b * sq(y - sq(x));
    }

    // Gradient of cost
    void eval_grad_f(crvec xy, rvec grad) const {
        auto x = xy(0), y = xy(1);
        grad(0) = (2 * x) - (2 * a) - (4 * b * x * y) + (4 * b * x * sq(x));
        grad(1) = 2 * b * (y - sq(x));
    }

    // Helper function
    static real_t sq(real_t x) { return x * x; }
};

int main() {
    USING_ALPAQA_CONFIG(RosenbrockProblem::config_t);
    // Instantiate a problem
    RosenbrockProblem problem;

    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Define the solver to use
    using Direction = alpaqa::LBFGSDirection<config_t>;
    using Solver    = alpaqa::PANOCSolver<Direction>;

    // Settings for the inner PANOC solver
    Solver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 1;
    // Settings for the L-BFGS algorithm used by PANOC
    Direction::AcceleratorParams lbfgsparam;
    lbfgsparam.memory = 2;

    // Create a PANOC solver
    Solver solver{panocparam, lbfgsparam};

    // Initial guess
    vec x = vec::Zero(2); // decision variables

    // Solve the problem
    auto stats = solver(counted_problem, {.tolerance = 1e-8}, x);
    // y and x have been overwritten by the solution

    // Print the results
    std::cout << '\n' << *counted_problem.evaluations << '\n';
    std::cout << "status: " << stats.status << '\n'
              << "f = " << problem.eval_f(x) << '\n'
              << "iterations: " << stats.iterations << '\n'
              << "ε = " << stats.ε << '\n'
              << "elapsed time:     "
              << std::chrono::duration<double>{stats.elapsed_time}.count()
              << " s" << '\n'
              << "x = " << x.transpose() << '\n'
              << "avg τ = " << (stats.sum_τ / stats.count_τ) << '\n'
              << "L-BFGS rejected = " << stats.lbfgs_rejected << '\n'
              << "L-BFGS failures = " << stats.lbfgs_failures << '\n'
              << "Line search failures = " << stats.linesearch_failures << '\n'
              << std::endl;
}