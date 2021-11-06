#include <gtest/gtest.h>

#include <alpaqa-ref/fd.hpp>
#include <alpaqa-ref/panoc-ref.hpp>

TEST(PANOC, cosh) {
    using alpaqa::Box;
    using alpaqa::crvec;
    using alpaqa::inf;
    using alpaqa::NaN;
    using alpaqa::Problem;
    using alpaqa::real_t;
    using alpaqa::rvec;
    using alpaqa::vec;

    const unsigned n = 2;
    const unsigned m = 1;

    Box C{vec(n), vec(n)};
    C.lowerbound(0) = -inf;
    C.lowerbound(1) = -inf;
    C.upperbound(0) = inf;
    C.upperbound(1) = -1;
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(-inf);
    D.upperbound.fill(0);

    real_t λ1  = 1e2;
    real_t λ2  = 1e-6;
    auto obj_f = [=](crvec x) {
        return std::cosh(λ1 * x(0)) + λ2 * x(1) * x(1);
    };
    auto grad_f = [=](crvec x, rvec grad_f) {
        grad_f(0) = λ1 * std::sinh(λ1 * x(0));
        grad_f(1) = 2 * λ2 * x(1);
    };
    auto g = [=](crvec x, rvec g_x) { g_x(0) = x(0) + x(1); };

    auto grad_g = [=]([[maybe_unused]] crvec x, crvec v, rvec grad_u_v) {
        alpaqa::mat grad = alpaqa::mat::Ones(n, m);
        grad_u_v     = grad * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g, {}, {}, {}};

    alpaqa::PANOCParams params;
    params.max_iter = 3;
    params.τ_min    = 1. / 16;
    alpaqa::LBFGSParams lbfgsparams;
    lbfgsparams.memory = 20;

    pa_ref::PANOCSolver solver{params, lbfgsparams};

    vec y = vec::Ones(m);
    vec x = vec(n);
    x << 1e-3, -1e8;
    vec z     = vec::Constant(m, NaN);
    vec err_z = vec::Constant(m, NaN);

    vec Σ(m);
    Σ.fill(20);

    real_t ε = 1e-4;

    auto stats = solver(p, Σ, ε, true, x, y, err_z);

    std::cout << "\n===========\n" << std::endl;
    std::cout << "f(x)     = " << obj_f(x) << std::endl;
    std::cout << "x        = " << x.transpose() << std::endl;
    std::cout << "y        = " << y.transpose() << std::endl;
    std::cout << "z        = " << z.transpose() << std::endl;
    std::cout << "g(x) - z = " << err_z.transpose() << std::endl;
    std::cout << "Iter:   " << stats.iterations << std::endl;
    std::cout << "Status: " << stats.status << std::endl << std::endl;
}