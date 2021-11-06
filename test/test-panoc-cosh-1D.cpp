#include <gtest/gtest.h>

#include <alpaqa/decl/alm.hpp>
#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/directions/decl/lbfgs.hpp>

#include <alpaqa-ref/panoc-ref.hpp>

TEST(PANOC, cosh1D) {
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
    D.lowerbound.fill(-inf);
    D.upperbound.fill(inf);

    real_t a    = 0.5;
    auto obj_f  = [=](crvec x) { return std::cosh(a * (x(0) - 10)); };
    auto grad_f = [=](crvec x, rvec grad_f) {
        grad_f(0) = a * std::sinh(a * (x(0) - 10));
    };
    auto g = [=](crvec x, rvec g_x) { g_x(0) = x(0); };

    auto grad_g = [=]([[maybe_unused]] crvec x, crvec v, rvec grad_u_v) {
        alpaqa::mat grad = alpaqa::mat::Ones(n, m);
        grad_u_v     = grad * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g, {}, {}, {}};

    alpaqa::PANOCParams params;
    params.max_iter = 100;
    params.τ_min    = 1. / 16;
    alpaqa::LBFGSParams lbfgsparams;
    lbfgsparams.memory = 20;

    pa_ref::PANOCSolver solver{params, lbfgsparams};

    vec y     = vec::Ones(m);
    vec x     = vec::Ones(n);
    vec err_z = vec::Constant(m, NaN);

    vec Σ(m);
    Σ.fill(1);

    real_t ε = 1e-4;

    auto stats = solver(p, Σ, ε, true, x, y, err_z);

    std::cout << "\n===========\n" << std::endl;
    std::cout << "f(x)     = " << obj_f(x) << std::endl;
    std::cout << "x        = " << x.transpose() << std::endl;
    std::cout << "y        = " << y.transpose() << std::endl;
    std::cout << "g(x) - z = " << err_z.transpose() << std::endl;
    std::cout << "Iter:   " << stats.iterations << std::endl;
    std::cout << "Status: " << stats.status << std::endl << std::endl;
}