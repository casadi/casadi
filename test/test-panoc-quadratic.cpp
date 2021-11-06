#include "eigen-matchers.hpp"
#include <gtest/gtest.h>

#include <alpaqa-ref/fd.hpp>
#include <alpaqa-ref/panoc-ref.hpp>

#include <iomanip>

TEST(PANOC, quadratic) {
    using alpaqa::Box;
    using alpaqa::crvec;
    using alpaqa::inf;
    using alpaqa::NaN;
    using alpaqa::Problem;
    using alpaqa::real_t;
    using alpaqa::rvec;
    using alpaqa::vec;

    const unsigned n = 1;
    const unsigned m = 1;

    Box C{vec(n), vec(n)};
    C.lowerbound(0) = -inf;
    C.upperbound(0) = inf;
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(6);
    D.upperbound.fill(inf);

    real_t a    = 1;
    real_t b    = -10;
    auto obj_f  = [=](crvec x) { return a * x(0) * x(0) + b * x(0); };
    auto grad_f = [=](crvec x, rvec grad_f) { grad_f(0) = 2 * a * x(0) + b; };
    auto g      = [=](crvec x, rvec g_x) { g_x(0) = x(0); };

    auto grad_g = [=]([[maybe_unused]] crvec x, crvec v, rvec grad_u_v) {
        alpaqa::mat grad = alpaqa::mat::Ones(n, m);
        grad_u_v     = grad * v;
    };
    auto g_fun = [=](crvec x) {
        vec gg(m);
        g(x, gg);
        return gg;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g, {}, {}, {}};

    alpaqa::PANOCParams params;
    params.max_iter = 10;
    params.τ_min    = 1. / 16;

    alpaqa::LBFGSParams lbfgsparam;
    lbfgsparam.memory = 20;

    pa_ref::PANOCSolver solver{params, lbfgsparam};

    vec y₀ = vec::Ones(m);
    vec y  = y₀;
    vec x  = vec(n);
    x << 1;
    vec err_z = vec::Constant(m, NaN);

    vec Σ(m);
    Σ.fill(1e10);

    real_t ε = 1e-10;

    auto stats = solver(p, Σ, ε, true, x, y, err_z);

    std::cout << std::setprecision(17);

    vec gg = g_fun(x);
    vec z  = alpaqa::project(gg + Σ.asDiagonal().inverse() * y₀, D);
    std::cout << "\n===========\n" << std::endl;
    std::cout << "f(x)     = " << obj_f(x) << std::endl;
    std::cout << "x        = " << x.transpose() << std::endl;
    std::cout << "y        = " << y.transpose() << std::endl;
    std::cout << "z        = " << z.transpose() << std::endl;
    std::cout << "g(x)     = " << gg.transpose() << std::endl;
    std::cout << "g(x) - z = " << err_z.transpose() << std::endl;
    std::cout << "Iter:   " << stats.iterations << std::endl;
    std::cout << "Status: " << stats.status << std::endl << std::endl;

    EXPECT_NEAR(x(0), 6, 2e-5);
    EXPECT_NEAR(y(0), -2, 2e-5);
    EXPECT_THAT(print_wrap(err_z), EigenAlmostEqual(gg - z, 1e-15));
}