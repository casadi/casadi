#pragma once

#include <panoc-alm/decl/alm.hpp>

namespace pa::detail {

inline void project_y(vec &y,          // inout
                      const vec &z_lb, // in
                      const vec &z_ub, // in
                      real_t M         // in
) {
    // TODO: Handle NaN correctly
    auto max_lb = [M](real_t y, real_t z_lb) {
        real_t y_lb = z_lb == -inf ? 0 : -M;
        return std::max(y, y_lb);
    };
    auto min_ub = [M](real_t y, real_t z_ub) {
        real_t y_ub = z_ub == inf ? 0 : M;
        return std::min(y, y_ub);
    };
    y = y.binaryExpr(z_lb, max_lb).binaryExpr(z_ub, min_ub);
}

inline void update_penalty_weights(const ALMParams &params, unsigned iteration,
                                   vec &e, vec &old_e, real_t norm_e, vec &Σ) {
    if (norm_e <= params.δ)
        return;
    for (Eigen::Index i = 0; i < e.rows(); ++i) {
        if (iteration == 0 || std::abs(e(i)) > params.θ * std::abs(old_e(i))) {
            Σ(i) = std::fmin(params.Σₘₐₓ,
                             std::fmax(params.Δ * std::abs(e(i)) / norm_e, 1) *
                                 Σ(i));
        }
    }
}

inline void initialize_penalty(const Problem &p, const ALMParams &params,
                               const vec &x0, vec &Σ) {
    real_t f0 = p.f(x0);
    vec g0(p.m);
    p.g(x0, g0);
    // TODO: reuse evaluations of f ang g in PANOC?
    real_t σ = params.σ₀ * std::abs(f0) / g0.squaredNorm();
    σ        = std::max(σ, params.σ₀);
    Σ.fill(σ);
}

inline void apply_preconditioning(const Problem &problem, Problem &prec_problem,
                                  const vec &x, real_t &prec_f, vec &prec_g) {
    vec grad_f(problem.n);
    vec v = vec::Zero(problem.m);
    vec grad_g(problem.n);
    problem.grad_f(x, grad_f);
    prec_f = 1. / std::max(grad_f.lpNorm<Eigen::Infinity>(), 1.);
    prec_g.resize(problem.m);

    for (Eigen::Index i = 0; i < problem.m; ++i) {
        v(i) = 1;
        problem.grad_g(x, v, grad_g);
        v(i)      = 0;
        prec_g(i) = 1. / std::max(grad_g.lpNorm<Eigen::Infinity>(), 1.);
    }

    auto prec_f_fun      = [&](const vec &x) { return problem.f(x) * prec_f; };
    auto prec_grad_f_fun = [&](const vec &x, vec &grad_f) {
        problem.grad_f(x, grad_f);
        grad_f *= prec_f;
    };
    auto prec_g_fun = [&](const vec &x, vec &g) {
        problem.g(x, g);
        g = prec_g.asDiagonal() * g;
    };
    auto prec_grad_g_fun = [&](const vec &x, const vec &v, vec &grad_g) {
        vec prec_v = prec_g.asDiagonal() * v;
        problem.grad_g(x, prec_v, grad_g);
    };
    prec_problem = Problem{
        problem.n,
        problem.m,
        problem.C,
        Box{prec_g.asDiagonal() * problem.D.upperbound,
            prec_g.asDiagonal() * problem.D.lowerbound},
        std::move(prec_f_fun),
        std::move(prec_grad_f_fun),
        std::move(prec_g_fun),
        std::move(prec_grad_g_fun),
    };

    // if (prec_g.size() <= 10) {
    //     std::cout << "prec_g:   " << prec_g.transpose() << std::endl;
    // } else {
    //     std::cout << "‖prec_g‖: " << prec_g.norm() << std::endl;
    // }
    // std::cout << "prec_f:   " << prec_f << std::endl;
}

} // namespace pa::detail