/** 
 * @example CasADi/Rosenbrock/main.cpp
 *
 * This example shows how to generate a problem using CasADi and how to load
 * and solve it using alpaqa.
 *
 * # Problem generation using CasADi
 * @include CasADi/Rosenbrock/codegen-rosenbrock.py
 * # Problem solution using alpaqa
 */

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/problem/wrapped-problem-with-counters.hpp>

#include <alpaqa/interop/casadi/CasADiLoader.hpp>

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

    // Load the problem (with 3 decision variables and 1 general constraint)
    auto p = alpaqa::CasADiProblem<config_t>(so_name.string(), 3, 1);

    // Specify the bounds
    p.C.upperbound = vec::Constant(3, alpaqa::inf<config_t>);
    p.C.lowerbound = vec::Constant(3, -alpaqa::inf<config_t>);
    p.D.upperbound = vec::Constant(1, 0.);
    p.D.lowerbound = vec::Constant(1, 0.);

    // Define the solvers to use
    using Accelerator = alpaqa::LBFGS<config_t>;
    using InnerSolver = alpaqa::PANOCSolver<Accelerator>;
    using OuterSolver = alpaqa::ALMSolver<InnerSolver>;

    // Settings for the outer augmented Lagrangian method
    OuterSolver::Params almparam;
    almparam.ε              = 1e-8; // tolerance
    almparam.δ              = 1e-8;
    almparam.Δ              = 10;
    almparam.max_iter       = 20;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    InnerSolver::Params panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    Accelerator::Params lbfgsparam;
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
    p.param(0) = 100;

    // Evaluation counters
    auto pc = alpaqa::with_counters(p);

    // Solve the problem
    auto stats = solver(pc, y, x);

    // Print the results
    vec g(p.m);
    p.eval_g(x, g);
    std::cout << "status: " << stats.status << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "g = " << g.transpose() << std::endl;
    std::cout << "f = " << p.eval_f(x) << std::endl;
    std::cout << "inner: " << stats.inner.iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl << std::endl;

    std::cout << pc.evaluations << std::endl;
}