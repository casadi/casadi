#pragma once

#include <alpaqa/decl/alm.hpp>
#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/directions/decl/lbfgs.hpp>

#include "eigen-matchers.hpp"

#include <chrono>

template <class Data>
inline alpaqa::Problem make_problem() {
    using namespace alpaqa;

    unsigned nx = Data::nx(), nu = Data::nu();
    unsigned n = nu + nx;
    unsigned m = nx;

    Box C{vec(n), vec(n)};
    C.lowerbound.topRows(nu).fill(-1);
    C.lowerbound.bottomRows(nx).fill(-inf);
    C.upperbound.topRows(nu).fill(1);
    C.upperbound.bottomRows(nx).fill(inf);
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(0);
    D.upperbound.fill(0);

    auto A = Data::A();
    auto B = Data::B();

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    x0.fill(1);

    auto f = [=](crvec ux) {
        auto u = ux.topRows(nu);
        auto x = ux.bottomRows(nx);
        return x.dot(Q * x) + u.dot(R * u);
    };
    auto grad_f = [=](crvec ux, rvec grad_f) {
        auto u                = ux.topRows(nu);
        auto x                = ux.bottomRows(nx);
        grad_f.topRows(nu)    = 2 * (R * u);
        grad_f.bottomRows(nx) = 2 * (Q * x);
    };
    auto g = [=](crvec ux, rvec g_u) {
        auto u         = ux.topRows(nu);
        auto x         = ux.bottomRows(nx);
        g_u.topRows(m) = A * x0 + B * u - x;
    };
    auto grad_g = [=](crvec ux, crvec v, rvec grad_u_v) {
        (void)ux;
        grad_u_v.topRows(nu)    = B.transpose() * v;
        grad_u_v.bottomRows(nx) = -alpaqa::mat::Identity(nx, nx) * v;
    };

    return Problem{n, m, C, D, f, grad_f, g, grad_g};
}

inline void do_test(const alpaqa::Problem &p, const alpaqa::rvec expected_sol,
                    const alpaqa::rvec expected_lagrange_multipliers) {
    using namespace alpaqa;

    ALMParams almparam;
    almparam.ε        = 1e-4;
    almparam.δ        = 1e-4;
    almparam.Δ        = 5; ///< Factor used in updating the penalty parameters
    almparam.Σ₀       = 1; ///< Initial penalty parameter
    almparam.ε₀       = 1e-1; ///< Initial tolerance on x
    almparam.θ        = 0.25;
    almparam.ρ        = 1e-1;
    almparam.M        = 1e9;
    almparam.Σ_max    = 1e9;
    almparam.max_iter = 100;

    PANOCParams panocparam;
    panocparam.Lipschitz.ε = 1e-6;
    panocparam.Lipschitz.δ = 1e-12;
    panocparam.max_iter    = 10000;

    LBFGSParams lbfgsparam;

    ALMSolver<> solver{almparam, {panocparam, lbfgsparam}};

    vec x(p.n);
    x.fill(0);
    vec y(p.m);
    y.fill(0);

    auto stats = solver(p, y, x);
    x.fill(0);
    y.fill(0);

    constexpr unsigned N = 20;
    auto begin           = std::chrono::high_resolution_clock::now();
    unsigned it          = 0;
    for (unsigned i = 0; i < N; ++i) {
        x.fill(0);
        y.fill(0);
        stats = solver(p, y, x);
        it += stats.inner.iterations;
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "y = " << y.transpose() << std::endl;
    std::cout << "f(x) = " << p.f(x) << std::endl;
    auto gx = vec(p.m);
    p.g(x, gx);
    std::cout << "g(x) = " << gx.transpose() << std::endl;
    std::cout << "Inner: " << stats.inner.iterations
              << ", Outer: " << stats.outer_iterations << std::endl;

    std::cout << "Total iterations: " << it << std::endl;
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << duration.count() / N << "µs" << std::endl;
}
