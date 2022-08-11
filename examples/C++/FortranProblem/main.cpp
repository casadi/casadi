#include <alpaqa/panoc-alm.hpp>

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
struct FortranProblem : alpaqa::Problem<config_t> {
    using alpaqa::Problem<config_t>::Problem;

    real_t eval_f(crvec x) const override { return problem_eval_f(x.data()); }
    void eval_grad_f(crvec x, rvec fx) const override {
        problem_eval_grad_f(x.data(), fx.data());
    }
    void eval_g(crvec x, rvec gx) const override {
        problem_eval_g(x.data(), gx.data());
    }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const override {
        problem_eval_grad_g_prod(x.data(), y.data(), grad_gxy.data());
    }
};

int main() {
    // Instantiate a problem
    FortranProblem problem{problem_get_num_vars(), problem_get_num_constr()};

    // Specify the bounds
    vec b                = vec::Constant(problem.m, -1);
    problem.C.lowerbound = vec::Constant(problem.n, -alpaqa::inf<config_t>);
    problem.C.upperbound = vec::Constant(problem.n, alpaqa::inf<config_t>);
    problem.D.lowerbound = vec::Constant(problem.m, -alpaqa::inf<config_t>);
    problem.D.upperbound = b;

    // Define the solvers to use
    using Accelerator = alpaqa::LBFGS<config_t>;
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
    Accelerator::Params lbfgsparam;
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
    auto stats = solver(problem, y, x);
    // y and x have been overwritten by the solution

    // Print the results
    std::cout << "status: " << stats.status << std::endl;
    std::cout << "inner iterations: " << stats.inner.iterations << std::endl;
    std::cout << "outer iterations: " << stats.outer_iterations << std::endl;
    std::cout << "elapsed time:     " << stats.elapsed_time.count() * 1e-6
              << 's' << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "f = " << problem.eval_f(x) << std::endl;
}