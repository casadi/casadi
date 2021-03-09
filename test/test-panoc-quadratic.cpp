#include <gtest/gtest.h>

#include <panoc-alm-ref/fd.hpp>
#include <panoc-alm-ref/panoc-ref.hpp>

TEST(PANOC, quadratic) {
    using pa::Box;
    using pa::inf;
    using pa::NaN;
    using pa::Problem;
    using pa::real_t;
    using pa::vec;

    const unsigned n = 1;
    const unsigned m = 1;

    Box C{vec(n), vec(n)};
    C.lowerbound(0) = -inf;
    C.upperbound(0) = inf;
    Box D{vec(m), vec(m)};
    D.lowerbound.fill(-inf);
    D.upperbound.fill(inf);

    real_t a    = 1;
    real_t b    = -10;
    auto obj_f  = [=](const vec &x) { return a * x(0) * x(0) + b * x(0); };
    auto grad_f = [=](const vec &x, vec &grad_f) {
        grad_f(0) = 2 * a * x(0) + b;
    };
    auto g = [=](const vec &x, vec &g_x) { g_x(0) = x(0); };

    auto grad_g = [=]([[maybe_unused]] const vec &x, const vec &v,
                      vec &grad_u_v) {
        pa::mat grad = pa::mat::Ones(n, m);
        grad_u_v     = grad * v;
    };

    Problem p{n, m, C, D, obj_f, grad_f, g, grad_g};

    pa::PANOCParams params;
    params.lbfgs_mem = 20;
    params.max_iter  = 1;
    params.τ_min     = 1. / 16;

    pa_ref::PANOCSolver solver{params};

    vec y = vec::Ones(m);
    vec x = vec(n);
    x << 1;
    vec err_z = vec::Constant(m, NaN);

    vec Σ(m);
    Σ.fill(1);

    real_t ε = 1e-4;

    auto stats = solver(p, Σ, ε, x, y, err_z);

    std::cout << "\n===========\n" << std::endl;
    std::cout << "f(x)     = " << obj_f(x) << std::endl;
    std::cout << "x        = " << x.transpose() << std::endl;
    std::cout << "y        = " << y.transpose() << std::endl;
    std::cout << "g(x) - z = " << err_z.transpose() << std::endl;
    std::cout << "Iter:   " << stats.iterations << std::endl;
    std::cout << "Status: " << stats.status << std::endl << std::endl;
}