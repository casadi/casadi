#include <alpaqa/alm.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>

#include <iostream>

int main() {
    using alpaqa::vec;
    using alpaqa::rvec;
    using alpaqa::crvec;
    using alpaqa::mat;
    using alpaqa::inf;

    // Problem specification
    alpaqa::Problem problem(2, 1); // # decision variables, # constraints

    // minimize  ½ xᵀHx
    //  s.t.     Ax ≤ b
    mat H(2, 2);    H <<  3, -1, 
                        -1,  3;
    mat A(1, 2);    A <<  2,  1;
    vec b(1);       b << -1;

    problem.f           = [&](crvec x)                   { return 0.5 * x.dot(H * x); };
    problem.grad_f      = [&](crvec x, rvec gr)          { gr = H * x; };
    problem.g           = [&](crvec x, rvec g)           { g  = A * x; };
    problem.grad_g_prod = [&](crvec x, crvec y, rvec gr) { gr = A.transpose() * y; };

    // Specify the bounds
    problem.C.lowerbound = vec::Constant(problem.n, -inf);
    problem.C.upperbound = vec::Constant(problem.n, inf);
    problem.D.lowerbound = vec::Constant(problem.m, -inf);
    problem.D.upperbound = b;

    // Settings for the outer augmented Lagrangian method
    alpaqa::ALMParams almparam;
    almparam.ε              = 1e-8; // tolerance
    almparam.δ              = 1e-8;
    almparam.Δ              = 10;   // penalty update factor
    almparam.max_iter       = 20;
    almparam.print_interval = 1;

    // Settings for the inner PANOC solver
    alpaqa::PANOCParams panocparam;
    panocparam.max_iter       = 500;
    panocparam.print_interval = 10;
    // Settings for the L-BFGS algorithm used by PANOC
    alpaqa::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 2;

    // Create an ALM solver using PANOC as inner solver
    alpaqa::ALMSolver<alpaqa::PANOCSolver<alpaqa::LBFGS>> solver{
        almparam,                 // params for outer solver
        {panocparam, lbfgsparam}, // inner solver
    };

    // Initial guess
    vec x(2);    x << 2, 2; // decision variables
    vec y(1);    y << 1;    // Lagrange multipliers

    // Solve the problem
    auto stats = solver(problem, y, x);
    // y and x have been overwritten by the solution

    // Print the results
    std::cout << "status: " << stats.status << std::endl;
    std::cout << "inner iterations: " << stats.inner.iterations << std::endl;
    std::cout << "outer iterations: " << stats.outer_iterations << std::endl;
    std::cout << "elapsed time:     " << stats.elapsed_time.count() * 1e-6 << 's' << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "f = " << problem.f(x) << std::endl;
}