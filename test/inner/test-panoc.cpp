#include <gtest/gtest.h>

#include <test-util/eigen-matchers.hpp>

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>

#include <alpaqa/config/config.hpp>
#include <alpaqa/implementation/inner/panoc-helpers.tpp>
#include <alpaqa/problem/functional-problem.hpp>

USING_ALPAQA_CONFIG(alpaqa::EigenConfigd);

using Helpers      = alpaqa::detail::PANOCHelpers<config_t>;
constexpr auto inf = alpaqa::inf<config_t>;

namespace pa_ref {
vec finite_diff(std::function<real_t(crvec)> f, crvec x) {
    const auto n = x.size();
    vec grad(n);
    vec h        = vec::Zero(n);
    const auto ε = std::sqrt(std::numeric_limits<real_t>::epsilon());
    const auto δ = std::numeric_limits<real_t>::min() / ε;
    for (index_t i = 0; i < n; ++i) {
        real_t hh        = std::abs(x(i)) > δ ? x(i) * ε : δ;
        h(i)             = hh;
        grad.coeffRef(i) = (f(x + h) - f(x)) / hh;
        h(i)             = 0;
    }
    return grad;
}
} // namespace pa_ref

auto build_test_problem() {
    alpaqa::FunctionalProblem<config_t> p{2, 2};
    p.C.upperbound = vec::Constant(2, inf);
    p.C.lowerbound = vec::Constant(2, -inf);
    p.D.upperbound = vec::Constant(2, 350);
    p.D.lowerbound = vec::Constant(2, -1);
    p.f            = [](crvec x) { // f(x) = 1/6 x₁⁴ + 2x₂ + x₂² + 1
        return 1. / 6 * std::pow(x(0), 4) + 2 * x(1) + std::pow(x(1), 2) + 1;
    };
    p.grad_f = [](crvec x, rvec grad) {
        grad(0) = 2. / 3 * std::pow(x(0), 3);
        grad(1) = 2 * x(1) + 2;
    };
    p.g = [](crvec x, rvec g) {
        g(0) = 2 * x(0) + 2 * x(1);
        g(1) = 3 * std::pow(x(0), 3);
    };
    p.grad_g_prod = [](crvec x, crvec y, rvec grad) {
        mat jacᵀ = mat::Zero(2, 2);
        jacᵀ << 2, 9 * std::pow(x(0), 2), 2, 0;
        grad = jacᵀ * y;
    };
    return p;
}

// Test the evaluation of PANOC's objective function ψ(x) and its gradient ∇ψ(x)
TEST(PANOC, calc_ψ_grad_ψ) {
    auto op = build_test_problem();
    auto p  = alpaqa::TypeErasedProblem<config_t>{op};

    auto f = [&p](crvec x) { return p.eval_f(x); };
    auto g = [&p](crvec x) {
        vec g(p.get_m());
        p.eval_g(x, g);
        return g;
    };
    auto grad_f = [&p](crvec x) {
        vec grad(p.get_n());
        p.eval_grad_f(x, grad);
        return grad;
    };
    auto grad_g_y = [&](crvec x, crvec y) {
        vec grad(p.get_n());
        p.eval_grad_g_prod(x, y, grad);
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
    vec invΣy = Σ.asDiagonal().inverse() * y;

    auto ψ_fun = [&op, &f, &g, &Σ, &invΣy](crvec x) -> real_t {
        return f(x) + 0.5 * alpaqa::dist_squared(g(x) + invΣy, op.D, Σ);
    };

    // Compute ψ and ∇ψ manually
    vec ζ     = g(x) + invΣy;
    vec ẑ     = alpaqa::project(ζ, op.D);
    vec d     = ζ - ẑ;
    vec ŷ     = Σ.asDiagonal() * d;
    real_t ψ  = f(x) + 0.5 * d.dot(ŷ);
    real_t ψ2 = ψ_fun(x);

    vec grad_ψ = grad_f(x) + grad_g_y(x, ŷ);

    vec grad_ψ_res(2);
    vec work_n(2), work_m(2);

    // Finite difference check
    vec grad_ψ_fd = pa_ref::finite_diff(ψ_fun, x);
    EXPECT_THAT(grad_ψ, EigenAlmostEqual(grad_ψ_fd, grad_ψ(0) * 5e-7));

    // calc_ψ_grad_ψ
    real_t ψ_res = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));
    EXPECT_DOUBLE_EQ(ψ_res, ψ);
    EXPECT_DOUBLE_EQ(ψ_res, ψ2);

    // calc_ψ
    work_m.setZero();
    ψ_res = p.eval_ψ(x, y, Σ, work_m);
    EXPECT_THAT(work_m, EigenAlmostEqual(ŷ, 1e-10));
    EXPECT_DOUBLE_EQ(ψ_res, ψ);

    // calc_grad_ψ_from_ŷ
    grad_ψ_res.setZero();
    p.eval_grad_L(x, work_m, grad_ψ_res, work_n);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));

    // calc_grad_ψ
    grad_ψ_res.setZero();
    work_n.setZero();
    work_m.setZero();
    p.eval_grad_ψ(x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));

    // calc_err_z
    vec err_z = g(x) - ẑ;
    vec err_z_res(2);
    Helpers::calc_err_z(alpaqa::TypeErasedProblem<config_t>{p}, x, y, Σ,
                        err_z_res);
    EXPECT_THAT(err_z_res, EigenAlmostEqual(err_z, 1e-10));
}

auto build_test_problem2() {
    alpaqa::FunctionalProblem<config_t> p{2, 2};
    p.C.upperbound = vec::Constant(2, 10);
    p.C.lowerbound = vec::Constant(2, -10);
    p.D.upperbound.resize(2);
    p.D.upperbound << -1, 0;
    p.D.lowerbound.resize(2);
    p.D.lowerbound << -inf, -inf;
    p.f = [](crvec x) { // f(x) = 1/48 x₁⁴ - 2x₁x₂ + 1/24 x₁²x₂⁴ + 10
        return 1. / 48 * std::pow(x(0), 4) - 2 * x(0) * x(1) +
               1. / 24 * std::pow(x(0), 2) * std::pow(x(1), 4) + 10;
    };
    p.grad_f = [](crvec x, rvec grad) {
        grad(0) = 1. / 12 * std::pow(x(0), 3) - 2 * x(1) +
                  1. / 12 * x(0) * std::pow(x(1), 4);
        grad(1) = -2 * x(0) + 1. / 6 * std::pow(x(0), 2) * std::pow(x(1), 3);
    };
    p.g = [](crvec x, rvec g) {
        g(0) = -4 * std::pow(x(0), 2) +
               0.25 * std::pow(x(0), 2) * std::pow(x(1), 2);
        g(1) = 0.125 * std::pow(x(0), 4) - x(0) * x(1);
    };
    p.grad_g_prod = [](crvec x, crvec y, rvec grad) {
        mat gradmat(2, 2);
        gradmat <<                                      //
            -8 * x(0) + 0.5 * x(0) * std::pow(x(1), 2), //
            0.5 * std::pow(x(0), 3) - x(1),             //
            0.5 * std::pow(x(0), 2) * x(1),             //
            -x(0);
        grad = gradmat * y;
    };
    p.grad_gi = [](crvec x, index_t i, rvec grad) {
        mat gradmat(2, 2);
        gradmat <<                                      //
            -8 * x(0) + 0.5 * x(0) * std::pow(x(1), 2), //
            0.5 * std::pow(x(0), 3) - x(1),             //
            0.5 * std::pow(x(0), 2) * x(1),             //
            -x(0);
        grad = gradmat.col(i);
    };
    p.hess_L = [](crvec x, crvec y, real_t scale, rmat H) {
        // Hessian of f
        H(0, 0) = 1. / 4 * std::pow(x(0), 2) + 1. / 12 * std::pow(x(1), 4);
        H(0, 1) = -2 + 1. / 3 * x(0) * std::pow(x(1), 3);
        H(1, 0) = H(0, 1);
        H(1, 1) = 1. / 2 * std::pow(x(0), 2) * std::pow(x(1), 2);
        H *= scale;

        mat Hg1(2, 2);
        Hg1(0, 0) = 0.5 * std::pow(x(1), 2) - 8;
        Hg1(0, 1) = x(0) * x(1);
        Hg1(1, 0) = x(0) * x(1);
        Hg1(1, 1) = 0.5 * std::pow(x(0), 2);

        mat Hg2(2, 2);
        Hg2(0, 0) = 1.5 * std::pow(x(0), 2);
        Hg2(0, 1) = -1;
        Hg2(1, 0) = -1;
        Hg2(1, 1) = 0;

        H += y(0) * Hg1;
        H += y(1) * Hg2;
    };
    return p;
}

// Test the evaluation of PANOC's augmented Lagrangian hessian
TEST(PANOC, hessian) {
    auto op = build_test_problem2();
    auto p  = alpaqa::TypeErasedProblem<config_t>{op};

    auto f = [&p](crvec x) { return p.eval_f(x); };
    auto g = [&p](crvec x) {
        vec g(p.get_m());
        p.eval_g(x, g);
        return g;
    };
    auto g1     = [&g](crvec x) { return g(x)(0); };
    auto g2     = [&g](crvec x) { return g(x)(1); };
    auto grad_f = [&p](crvec x) {
        vec grad(p.get_n());
        p.eval_grad_f(x, grad);
        return grad;
    };
    auto grad_g_y = [&](crvec x, crvec y) {
        vec grad(p.get_n());
        p.eval_grad_g_prod(x, y, grad);
        return grad;
    };
    auto grad_gi = [&p](crvec x, unsigned i) {
        vec grad(p.get_n());
        p.eval_grad_gi(x, i, grad);
        return grad;
    };
    auto grad_g1 = [&grad_gi](crvec x) { return grad_gi(x, 0); };
    auto grad_g2 = [&grad_gi](crvec x) { return grad_gi(x, 1); };

    std::scientific(std::cout);

    // Some arbitrary vectors Σ, x, y
    vec Σ(2);
    Σ << 10, 10;
    vec x(2);
    x << -0.9, 3.1;
    vec y(2);
    y << 0.3, 0.7;
    vec invΣy = Σ.asDiagonal().inverse() * y;

    auto ψ_fun = [&op, &f, &g, &Σ, &invΣy](crvec x) -> real_t {
        return f(x) + 0.5 * alpaqa::dist_squared(g(x) + invΣy, op.D, Σ);
    };

    // Compute ψ and ∇ψ manually
    vec ζ     = g(x) + invΣy;
    vec ẑ     = alpaqa::project(ζ, op.D);
    vec d     = ζ - ẑ;
    vec ŷ     = Σ.asDiagonal() * d;
    real_t ψ  = f(x) + 0.5 * d.dot(ŷ);
    real_t ψ2 = ψ_fun(x);

    vec grad_ψ = grad_f(x) + grad_g_y(x, ŷ);

    vec grad_ψ_res(2);
    vec work_n(2), work_m(2);

    // Finite difference check f
    vec grad_f_fd  = pa_ref::finite_diff(f, x);
    vec grad_f_res = grad_f(x);
    EXPECT_THAT(grad_f_res,
                EigenAlmostEqual(grad_f_fd, std::abs(grad_f_res(0)) * 1e-6));

    // Finite difference check g₁
    vec grad_g1_fd  = pa_ref::finite_diff(g1, x);
    vec grad_g1_res = grad_g1(x);
    EXPECT_THAT(grad_g1_res,
                EigenAlmostEqual(grad_g1_fd, std::abs(grad_g1_res(0)) * 1e-6));
    // Finite difference check g₂
    vec grad_g2_fd  = pa_ref::finite_diff(g2, x);
    vec grad_g2_res = grad_g2(x);
    EXPECT_THAT(grad_g2_res,
                EigenAlmostEqual(grad_g2_fd, std::abs(grad_g2_res(0)) * 1e-6));

    // Finite difference check ψ
    vec grad_ψ_fd = pa_ref::finite_diff(ψ_fun, x);
    EXPECT_THAT(grad_ψ,
                EigenAlmostEqual(grad_ψ_fd, std::abs(grad_ψ(0)) * 5e-6));

    // calc_ψ_grad_ψ
    real_t ψ_res = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));
    EXPECT_DOUBLE_EQ(ψ_res, ψ);
    EXPECT_DOUBLE_EQ(ψ_res, ψ2);

    // calc_ψ
    work_m.setZero();
    ψ_res = p.eval_ψ(x, y, Σ, work_m);
    EXPECT_THAT(work_m, EigenAlmostEqual(ŷ, 1e-10));
    EXPECT_DOUBLE_EQ(ψ_res, ψ);

    // calc_grad_ψ_from_ŷ
    grad_ψ_res.setZero();
    p.eval_grad_L(x, work_m, grad_ψ_res, work_n);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));

    // calc_grad_ψ
    grad_ψ_res.setZero();
    work_n.setZero();
    work_m.setZero();
    p.eval_grad_ψ(x, y, Σ, grad_ψ_res, work_n, work_m);
    EXPECT_THAT(grad_ψ_res, EigenAlmostEqual(grad_ψ, 1e-10));

    // calc_err_z
    vec err_z = g(x) - ẑ;
    vec err_z_res(2);
    Helpers::calc_err_z(alpaqa::TypeErasedProblem<config_t>{p}, x, y, Σ,
                        err_z_res);
    EXPECT_THAT(err_z_res, EigenAlmostEqual(err_z, 1e-10));

    // Hessian of f
    auto grad_fi = [&](index_t i) {
        return [&, i](crvec x) { return grad_f(x)(i); };
    };
    vec hess_f1_fd = pa_ref::finite_diff(grad_fi(0), x);
    vec hess_f2_fd = pa_ref::finite_diff(grad_fi(1), x);
    mat H_res(2, 2);
    indexvec inner, outer;
    p.eval_hess_L(x, vec::Zero(2), 1, inner, outer, H_res.reshaped());
    EXPECT_THAT(H_res.col(0),
                EigenAlmostEqual(hess_f1_fd, std::abs(H_res.col(0)(0)) * 5e-6));
    EXPECT_THAT(H_res.col(1),
                EigenAlmostEqual(hess_f2_fd, std::abs(H_res.col(1)(0)) * 5e-6));

    // Hessian of Lagrangian
    auto grad_Li = [&](unsigned i) {
        return [&, i](crvec x) {
            vec grad_L_res = grad_f(x) + grad_g_y(x, y);
            return grad_L_res(i);
        };
    };
    p.eval_hess_L(x, y, 1, inner, outer, H_res.reshaped());
    vec hess_L1_fd = pa_ref::finite_diff(grad_Li(0), x);
    vec hess_L2_fd = pa_ref::finite_diff(grad_Li(1), x);
    EXPECT_THAT(H_res.col(0),
                EigenAlmostEqual(hess_L1_fd, std::abs(H_res.col(0)(0)) * 5e-6));
    EXPECT_THAT(H_res.col(1),
                EigenAlmostEqual(hess_L2_fd, std::abs(H_res.col(1)(0)) * 5e-6));

    // Hessian of augmented Lagrangian
    auto grad_ψi = [&](unsigned i) {
        return [&, i](crvec x) {
            vec grad_ψ_res(2);
            p.eval_grad_ψ(x, y, Σ, grad_ψ_res, work_n, work_m);
            return grad_ψ_res(i);
        };
    };

    // TODO: enable again when implemented
    static_cast<void>(grad_ψi);
    // vec gv(2);
    // Helpers::calc_augmented_lagrangian_hessian(p, x, ŷ, y, Σ, gv, H_res,
    //                                            work_n);
    // vec hess_ψ1_fd = pa_ref::finite_diff(grad_ψi(0), x);
    // vec hess_ψ2_fd = pa_ref::finite_diff(grad_ψi(1), x);
    // EXPECT_THAT(H_res.col(0),
    //             EigenAlmostEqual(hess_ψ1_fd, std::abs(H_res.col(0)(0)) * 5e-6));
    // EXPECT_THAT(H_res.col(1),
    //             EigenAlmostEqual(hess_ψ2_fd, std::abs(H_res.col(1)(0)) * 5e-6));
}
