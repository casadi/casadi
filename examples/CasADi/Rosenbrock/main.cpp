/** 
 * @example CasADi/Rosenbrock/main.cpp
 *
 * This example shows how to generate a problem using CasADi and how to load
 * and solve it using PANOC-ALM.
 *
 * # Problem generation using CasADi
 * @include CasADi/Rosenbrock/codegen-rosenbrock.py
 * # Problem solution using PANOC-ALM
 */

#include <panoc-alm/decl/alm.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/directions/decl/lbfgs.hpp>

#include <panoc-alm/interop/casadi/CasADiLoader.hpp>

#include <iostream>

int main(int argc, char *argv[]) {
    using pa::inf;
    using pa::vec;

    auto so_name = "examples/CasADi/Rosenbrock/librosenbrock_functions.so";
    if (argc > 1)
        so_name = argv[1];

    // Load the problem (with 3 decision variables and 1 general constraint)
    pa::Problem p = pa::load_CasADi_problem(so_name, 3, 1);

    // Specify the bounds
    p.C.upperbound = vec::Constant(3, inf);
    p.C.lowerbound = vec::Constant(3, -inf);
    p.D.upperbound = vec::Constant(1, 0.);
    p.D.lowerbound = vec::Constant(1, 0.);

    // Settings for the outer augmented Lagrangian method
    pa::ALMParams almparam;
    almparam.ε              = 1e-8; // tolerance
    almparam.δ              = 1e-8;
    almparam.Δ              = 10;
    almparam.max_iter       = 20;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    pa::PANOCParams panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    pa::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 10;

    // Create an ALM solver using PANOC as inner solver
    pa::ALMSolver<pa::PANOCSolver<>> solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x(3);
    x << 2.5, 3.0, 0.75;
    vec y(1);
    y << 1;

    // Solve the problem
    auto stats = solver(p, y, x);

    // Print the results
    vec g(p.m);
    p.g(x, g);
    std::cout << "status: " << stats.status << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "g = " << g.transpose() << std::endl;
    std::cout << "f = " << p.f(x) << std::endl;
    std::cout << "inner: " << stats.inner.iterations << std::endl;
    std::cout << "outer: " << stats.outer_iterations << std::endl;
}