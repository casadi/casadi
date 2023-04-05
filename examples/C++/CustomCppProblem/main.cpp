#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>

#include <iostream>

int main() {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);

    // Problem specification
    // minimize  ½ xᵀQx
    //  s.t.     Ax ≤ b
    struct Problem : alpaqa::BoxConstrProblem<config_t> {
        mat Q{n, n};
        mat A{m, n};
        vec b{m};
        mutable vec Qx{n};

        Problem() : alpaqa::BoxConstrProblem<config_t>{2, 1} {
            // Initialize problem matrices
            Q << 3, -1, -1, 3;
            A << 2, 1;
            b << -1;

            // Specify the bounds
            C.lowerbound = vec::Constant(n, -alpaqa::inf<config_t>);
            C.upperbound = vec::Constant(n, +alpaqa::inf<config_t>);
            D.lowerbound = vec::Constant(m, -alpaqa::inf<config_t>);
            D.upperbound = b;
        }

        // Evaluate the cost
        real_t eval_f(crvec x) const {
            Qx.noalias() = Q * x;
            return 0.5 * x.dot(Qx);
        }
        // Evaluat the gradient of the cost
        void eval_grad_f(crvec x, rvec gr) const { gr.noalias() = Q * x; }
        // Evaluate the constraints
        void eval_g(crvec x, rvec g) const { g.noalias() = A * x; }
        // Evaluate a matrix-vector product with the gradient of the constraints
        void eval_grad_g_prod(crvec x, crvec y, rvec gr) const {
            (void)x;
            gr.noalias() = A.transpose() * y;
        }
    };

    Problem problem;

    // Wrap the problem to count the function evaluations
    auto counted_problem = alpaqa::problem_with_counters_ref(problem);

    // Define the solvers to use
    using Direction   = alpaqa::LBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Direction>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.tolerance             = 1e-8; // tolerance
    almparam.dual_tolerance        = 1e-8;
    almparam.penalty_update_factor = 10; // penalty update factor
    almparam.max_iter              = 20;
    almparam.print_interval        = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 1;
    // Settings for the L-BFGS algorithm used by PANOC
    Direction::AcceleratorParams lbfgsparam;
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