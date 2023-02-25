#include <alpaqa/casadi/CasADiProblem.hpp>
#include <alpaqa/panoc-alm.hpp>
#include <alpaqa/structured-panoc-alm.hpp>

#include <test-util/eigen-matchers.hpp>
#include <stdexcept>

TEST(ALM, casadi) {
    USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);

    // Load the problem
    alpaqa::CasADiProblem<config_t> problem{ROSENBROCK_FUNC_DLL};

    // Specify the bounds
    problem.C.lowerbound = vec::Constant(2, 0.);
    problem.C.upperbound = vec::Constant(2, 5.);
    problem.D.upperbound = vec::Constant(1, 1.);
    problem.D.lowerbound = vec::Constant(1, 0.);

    // Define the solvers to use
    using Accelerator = alpaqa::LBFGSDirection<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Accelerator>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.ε              = 1e-10; // tolerance
    almparam.δ              = 1e-10;
    almparam.Δ              = 10;
    almparam.max_iter       = 20;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    Accelerator::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 10;

    // Create an ALM solver using PANOC as inner solver
    OuterSolver solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x(2);
    x << 2, 1;
    vec y(1);
    y << 1;

    // Parameter
    problem.param(0) = 2.;

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

    vec x_expected(2);
    x_expected << +7.54681194348482909e-01, +4.63926877264482063e-01;
    vec y_expected(1);
    y_expected << +1.13829175485413714e-01;
    EXPECT_EQ(stats.status, alpaqa::SolverStatus::Converged);
    EXPECT_THAT(x, EigenAlmostEqualRel(x_expected, 1e-6));
    EXPECT_THAT(y, EigenAlmostEqualRel(y_expected, 1e-6));
}
