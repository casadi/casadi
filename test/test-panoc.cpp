#include "eigen-matchers.hpp"
#include <gtest/gtest.h>

#include <panoc-alm-ref/fd.hpp>
#include <panoc-alm-ref/panoc-ref.hpp>

#include <panoc-alm/inner/detail/panoc-helpers.hpp>
#include <panoc-alm/inner/decl/panoc.hpp>

using pa::inf;
using pa::Problem;
using pa::real_t;
using pa::vec;

Problem build_test_problem() {
    pa::Problem p;
    p.n            = 2;
    p.m            = 2;
    p.C.upperbound = vec::Constant(2, inf);
    p.C.lowerbound = vec::Constant(2, -inf);
    p.D.upperbound = vec::Constant(2, 350);
    p.D.lowerbound = vec::Constant(2, -1);
    p.f = [](const vec &x) { // f(x) = 1/6 x₁⁴ + 2x₂ + x₂² + 1
        return 1. / 6 * std::pow(x(0), 4) + 2 * x(1) + std::pow(x(1), 2) + 1;
    };
    p.grad_f = [](const vec &x, vec &grad) {
        grad(0) = 2. / 3 * std::pow(x(0), 3);
        grad(1) = 2 * x(1) + 2;
    };
    p.g = [](const vec &x, vec &g) {
        g(0) = 2 * x(0) + 2 * x(1);
        g(1) = 3 * std::pow(x(0), 3);
    };
    p.grad_g = [](const vec &x, const vec &y, vec &grad) {
        pa::mat jacᵀ = pa::mat::Zero(2, 2);
        jacᵀ << 2, 9 * std::pow(x(0), 2), 2, 0;
        grad = jacᵀ * y;
    };
    return p;
}

// Test the evaluation of PANOC's objective function ψ(x) and its gradient ∇ψ(x)
TEST(PANOC, calc_ψ_grad_ψ) {
    auto p = build_test_problem();

    auto f = p.f;
    auto g = [&p](const vec &x) {
        vec g(p.m);
        p.g(x, g);
        return g;
    };
    auto grad_f = [&p](const vec &x) {
        vec grad(p.n);
        p.grad_f(x, grad);
        return grad;
    };
    auto grad_g_y = [&](const vec &x, const vec &y) {
        vec grad(p.n);
        p.grad_g(x, y, grad);
        return grad;
    };

    std::scientific(std::cout);

    // Some arbitrary vectors Σ, x, y
    vec Σ(2);
    Σ << 3, 2;
    vec x(2);
    x << 5, -7;
    vec y(2);
    y << 0.3, 0.7;
    vec Σ⁻¹y = Σ.asDiagonal().inverse() * y;

    auto ψ_fun = [&p, &f, &g, &Σ, &Σ⁻¹y](const vec &x) -> real_t {
        return f(x) + 0.5 * pa::dist_squared(g(x) + Σ⁻¹y, p.D, Σ);
    };

    // Compute ψ and ∇ψ manually
    vec ζ     = g(x) + Σ⁻¹y;
    vec ẑ     = pa::project(ζ, p.D);
    vec d     = ζ - ẑ;
    vec ŷ     = Σ.asDiagonal() * d;
    real_t ψ  = f(x) + 0.5 * d.dot(ŷ);
    real_t ψ2 = ψ_fun(x);

    vec grad_ψ = grad_f(x) + grad_g_y(x, ŷ);

    vec grad_ψ_res(2);
    vec work_n(2), work_m(2);

    // Finite difference check
    vec grad_ψ_fd = pa_ref::finite_diff(ψ_fun, x);
    EXPECT_THAT(print_wrap(grad_ψ),
                EigenAlmostEqual(print_wrap(grad_ψ_fd), grad_ψ(0) * 5e-7));

    // calc_ψ_grad_ψ
    real_t ψ_res =
        pa::detail::calc_ψ_grad_ψ(p, x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(print_wrap(grad_ψ_res),
                EigenAlmostEqual(print_wrap(grad_ψ), 1e-10));
    EXPECT_FLOAT_EQ(ψ_res, ψ);
    EXPECT_FLOAT_EQ(ψ_res, ψ2);

    // calc_ψ_ŷ
    work_m.setZero();
    ψ_res = pa::detail::calc_ψ_ŷ(p, x, y, Σ, work_m);
    EXPECT_THAT(print_wrap(work_m), EigenAlmostEqual(print_wrap(ŷ), 1e-10));
    EXPECT_FLOAT_EQ(ψ_res, ψ);

    // calc_grad_ψ_from_ŷ
    grad_ψ_res.setZero();
    pa::detail::calc_grad_ψ_from_ŷ(p, x, work_m, grad_ψ_res, work_n);
    EXPECT_THAT(print_wrap(grad_ψ_res),
                EigenAlmostEqual(print_wrap(grad_ψ), 1e-10));

    // calc_grad_ψ
    grad_ψ_res.setZero();
    work_n.setZero();
    work_m.setZero();
    pa::detail::calc_grad_ψ(p, x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(print_wrap(grad_ψ_res),
                EigenAlmostEqual(print_wrap(grad_ψ), 1e-10));

    // calc_err_z
    vec err_z = g(x) - ẑ;
    vec err_z_res(2);
    pa::detail::calc_err_z(p, x, y, Σ, err_z_res);
    EXPECT_THAT(print_wrap(err_z_res),
                EigenAlmostEqual(print_wrap(err_z), 1e-10));

    // pa_ref
    ψ_res      = pa_ref::detail::eval_ψ(p, x, y, Σ);
    grad_ψ_res = pa_ref::detail::eval_grad_ψ(p, x, y, Σ);
    EXPECT_THAT(print_wrap(grad_ψ_res),
                EigenAlmostEqual(print_wrap(grad_ψ), 1e-10));
    EXPECT_FLOAT_EQ(ψ_res, ψ);
}

// Compare the reference implementation of PANOC with the optimized implementation
TEST(PANOC, ref) {
    using pa::Box;
    using pa::inf;
    using pa::NaN;
    using pa::Problem;
    using pa::real_t;
    using pa::vec;

    unsigned nu = 2;
    unsigned nx = 4;
    unsigned ns = 2;
    unsigned ny = 3;
    unsigned n  = nu + ns + ny;
    unsigned m  = 3;

    auto u = [nu](auto &&x) { return x.segment(0, nu); };
    auto s = [nu, ns](auto &&x) { return x.segment(nu, ns); };
    auto y = [nu, ns, ny](auto &&x) { return x.segment(nu + ns, ny); };

    Box C{vec(n), vec(n)};
    u(C.lowerbound).fill(-10);
    u(C.upperbound).fill(10);
    s(C.lowerbound).fill(-inf);
    s(C.upperbound).fill(inf);
    y(C.lowerbound).fill(0);
    y(C.upperbound).fill(inf);
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(-inf);
    D.upperbound.fill(0);

    real_t Ts = 0.05;

    pa::mat A = pa::mat::Identity(nx, nx);
    A(0, 2)   = Ts;
    A(1, 3)   = Ts;
    pa::mat B = pa::mat::Zero(nx, nu);
    B(2, 0)   = Ts;
    B(3, 1)   = Ts;

    auto f = [=](const vec &x, const vec &u) { return A * x + B * u; };

    using Diag = Eigen::DiagonalMatrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto Q = Diag(nx);
    Q.diagonal().fill(10);
    auto R = Diag(nu);
    R.diagonal().fill(1);

    auto x0 = vec(nx);
    x0.fill(10);

    real_t scale = 1;
    auto obj_f   = [=](const vec &ux) { return scale * s(ux)(0); };
    auto grad_f  = [=](const vec &ux, vec &grad_f) {
        (void)ux;
        grad_f.fill(0);
        s(grad_f)(0) = scale;
    };
    auto g = [=](const vec &ux, vec &g_u) {
        g_u(0) = y(ux)(0) - y(ux)(1) - s(ux)(0);
        g_u(1) = f(x0, u(ux)).dot(Q * f(x0, u(ux))) - s(ux)(1);
        g_u(2) = x0.dot(Q * x0) + u(ux).dot(R * u(ux)) -
                 (y(ux)(0) - y(ux)(1) - y(ux)(2) - s(ux)(1));
    };
    auto grad_g_mat = [=](const vec &ux) {
        pa::mat grad      = pa::mat::Zero(n, m);
        s(grad.col(0))(0) = -1;
        y(grad.col(0))(0) = 1;
        y(grad.col(0))(1) = -1;

        u(grad.col(1))    = 2 * B.transpose() * (Q * (A * x0 + B * u(ux)));
        s(grad.col(1))(1) = -1;

        u(grad.col(2))    = 2 * (R * u(ux));
        s(grad.col(2))(1) = 1;
        y(grad.col(2))(0) = -1;
        y(grad.col(2))(1) = 1;
        y(grad.col(2))(2) = 1;

        return grad;
    };
    auto grad_g = [=](const vec &ux, const vec &v, vec &grad_u_v) {
        grad_u_v = grad_g_mat(ux) * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g};

    pa::PANOCParams params;
    params.lbfgs_mem                      = 20;
    params.max_iter                       = 1000;
    params.τ_min                          = 1. / (1 << 10);
    params.update_lipschitz_in_linesearch = true;

    pa::PANOCParams params_ref = params;

    pa::PANOCSolver solver{params};
    pa_ref::PANOCSolver solver_ref{params_ref};

    vec λ     = vec::Ones(m);
    vec λ_ref = vec::Ones(m);

    vec x     = vec::Zero(n);
    vec x_ref = vec::Zero(n);

    vec err_z     = vec::Constant(m, NaN);
    vec err_z_ref = vec::Constant(m, NaN);

    vec Σ(m);
    Σ.fill(1e-2);

    real_t ε = 1e-5;

    auto stats     = solver(p, Σ, ε, x, λ, err_z);
    auto stats_ref = solver_ref(p, Σ, ε, x_ref, λ_ref, err_z_ref);

    std::cout << "Normal" << std::endl;
    std::cout << "u = " << u(x).transpose() << "\ts = " << s(x).transpose()
              << "\ty = " << y(x).transpose() << std::endl;
    std::cout << stats.iterations << std::endl;
    std::cout << stats.status << std::endl << std::endl;

    std::cout << "Ref" << std::endl;
    std::cout << "u = " << u(x_ref).transpose()
              << "\ts = " << s(x_ref).transpose()
              << "\ty = " << y(x_ref).transpose() << std::endl;
    std::cout << stats_ref.iterations << std::endl;
    std::cout << stats_ref.status << std::endl << std::endl;

    EXPECT_EQ(stats.status, pa::SolverStatus::Converged);
    EXPECT_EQ(stats.status, stats_ref.status);
    EXPECT_EQ(stats.iterations, stats_ref.iterations);
    EXPECT_EQ(stats.lbfgs_failures, stats_ref.lbfgs_failures);
    EXPECT_EQ(stats.lbfgs_rejected, stats_ref.lbfgs_rejected);
    EXPECT_THAT(print_wrap(x), EigenAlmostEqual(print_wrap(x_ref), 1e-5));
    EXPECT_THAT(print_wrap(λ), EigenAlmostEqual(print_wrap(λ_ref), 1e-5));
}
