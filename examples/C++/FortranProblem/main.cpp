#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/structured-panoc-alm.hpp>

#include <iostream>

// Double precision, same as in Fortran
USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);

// External Fortran routines
extern "C" {
int64_t problem_get_num_vars(void);
int64_t problem_get_num_constr(void);
double problem_eval_f(const double *);
void problem_eval_grad_f(const double *, double *);
void problem_eval_g(const double *, double *);
void problem_eval_grad_g_prod(const double *, const double *, double *);
}

// Problem specification by inheriting from alpaqa::Problem
struct FortranProblem : alpaqa::BoxConstrProblem<config_t> {
    using alpaqa::BoxConstrProblem<config_t>::BoxConstrProblem;

    real_t eval_f(crvec x) const { return problem_eval_f(x.data()); }
    void eval_grad_f(crvec x, rvec fx) const {
        problem_eval_grad_f(x.data(), fx.data());
    }
    void eval_g(crvec x, rvec gx) const { problem_eval_g(x.data(), gx.data()); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const {
        problem_eval_grad_g_prod(x.data(), y.data(), grad_gxy.data());
    }
};

int main() {
    // Instantiate a problem
    FortranProblem problem{problem_get_num_vars(), problem_get_num_constr()};
    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Specify the bounds
    vec b                = vec::Constant(problem.m, -1);
    problem.C.lowerbound = vec::Constant(problem.n, -alpaqa::inf<config_t>);
    problem.C.upperbound = vec::Constant(problem.n, +alpaqa::inf<config_t>);
    problem.D.lowerbound = vec::Constant(problem.m, -alpaqa::inf<config_t>);
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