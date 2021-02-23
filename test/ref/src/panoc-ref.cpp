#include <cassert>
#include <cmath>

#include <limits>
#include <panoc-alm-ref/panoc-ref.hpp>
#include <panoc-alm/lbfgs.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pa_ref {

using pa::dist_squared;
using pa::LBFGS;
using pa::project;

namespace detail {

vec eval_g(const Problem &p, const vec &x) {
    vec g(p.m);
    p.g(x, g);
    return g;
}

vec eval_ẑ(const Problem &p, const vec &x, const vec &y, const vec &Σ) {
    vec g = eval_g(p, x);
    return project(g + Σ.asDiagonal().inverse() * y, p.D);
}

vec eval_ŷ(const Problem &p, const vec &x, const vec &y, const vec &Σ) {
    vec g = eval_g(p, x);
    vec ẑ = eval_ẑ(p, x, y, Σ);
    return Σ.asDiagonal() * (g - ẑ) + y;
}

real_t eval_ψ(const Problem &p, const vec &x, const vec &y, const vec &Σ) {
    vec g = eval_g(p, x);
    return p.f(x) +
           0.5 * dist_squared(g + Σ.asDiagonal().inverse() * y, p.D, Σ);
}

vec eval_grad_ψ(const Problem &p, const vec &x, const vec &y, const vec &Σ) {
    vec ŷ = eval_ŷ(p, x, y, Σ);
    vec grad_f(p.n), grad_gŷ(p.n);
    p.grad_f(x, grad_f);
    p.grad_g(x, ŷ, grad_gŷ);
    return grad_f + grad_gŷ;
}

vec gradient_step(const Problem &p, const vec &x, const vec &y, const vec &Σ,
                  real_t γ) {
    return x - γ * eval_grad_ψ(p, x, y, Σ);
}

vec T_γ(const Problem &p, const vec &x, const vec &y, const vec &Σ, real_t γ) {
    return project(gradient_step(p, x, y, Σ, γ), p.C);
}

real_t eval_φ(const Problem &p, const vec &x, const vec &r, const vec &y,
              const vec &Σ, real_t γ) {
    return eval_ψ(p, x, y, Σ)                //
           + 1. / (2 * γ) * r.squaredNorm()  //
           + eval_grad_ψ(p, x, y, Σ).dot(r); //
}

real_t estimate_lipschitz(const Problem &p, const vec &x, const vec &y,
                          const vec &Σ, const PANOCParams &params) {
    // Estimate Lipschitz constant using finite difference
    vec h = (x * params.Lipschitz.ε).cwiseMax(params.Lipschitz.δ);
    // Calculate ∇ψ(x₀ + h) and ∇ψ(x₀)
    vec diff = eval_grad_ψ(p, x + h, y, Σ) - eval_grad_ψ(p, x, y, Σ);
    // Estimate Lipschitz constant
    real_t L = diff.norm() / h.norm();
    return L;
}

real_t calc_error_stop_crit(const Problem &p, const vec &xₖ, const vec &x̂ₖ,
                            const vec &y, const vec &Σ, real_t γ) {
    return ((1 / γ) * (xₖ - x̂ₖ)        //
            + eval_grad_ψ(p, x̂ₖ, y, Σ) //
            - eval_grad_ψ(p, xₖ, y, Σ))
        .lpNorm<Eigen::Infinity>();
}

bool lipschitz_check(const Problem &p, const vec &xₖ, const vec &x̂ₖ,
                     const vec &y, const vec &Σ, real_t γ, real_t L) {
    real_t ψₖ     = eval_ψ(p, xₖ, y, Σ);
    real_t ψ̂xₖ    = eval_ψ(p, x̂ₖ, y, Σ);
    real_t margin = 0; // 1e-6 * std::abs(ψₖ); // TODO: why?
    vec rₖ        = xₖ - x̂ₖ;
    return ψ̂xₖ > margin + ψₖ                            //
                     - eval_grad_ψ(p, xₖ, y, Σ).dot(rₖ) //
                     + L / 2 * rₖ.squaredNorm();
}

} // namespace detail

PANOCSolver::Stats PANOCSolver::operator()(
    const Problem &problem, ///< [in]    Problem description
    vec &x,                 ///< [inout] Decision variable @f$ x @f$
    vec &z,                 ///< [out]   Slack variable @f$ x @f$
    vec &y,                 ///< [inout] Lagrange multiplier @f$ x @f$
    vec &err_z,             ///< [out]   Slack variable error @f$ g(x) - z @f$
    const vec &Σ,           ///< [in]    Constraint weights @f$ \Sigma @f$
    real_t ε                ///< [in]    Tolerance @f$ \epsilon @f$
) {
    using namespace detail;

    const auto n = problem.n;

    LBFGS lbfgs(n, params.lbfgs_mem);

    vec xₖ = x; // Value of x at the beginning of the iteration

    real_t L = estimate_lipschitz(problem, x, y, Σ, params);
    if (L < 10 * std::numeric_limits<real_t>::epsilon())
        L = 10 * std::numeric_limits<real_t>::epsilon();
    real_t γ = params.Lipschitz.Lγ_factor / L;
    real_t σ = γ * (1 - γ * L) / 2;

    std::cout << std::scientific;

    std::cout << "Initial γ = " << γ << "\n";

    // Calculate x̂ₖ, rₖ (gradient step)
    vec x̂ₖ = T_γ(problem, xₖ, y, Σ, γ);
    // Increase L until quadratic upper bound is satisfied
    real_t old_γ = γ;
    while (lipschitz_check(problem, xₖ, x̂ₖ, y, Σ, γ, L) && γ != 0) {
        L *= 2;
        σ /= 2;
        γ /= 2;

        // Calculate x̂ₖ (with new step size)
        x̂ₖ = T_γ(problem, xₖ, y, Σ, γ);
    }
    if (old_γ != γ)
        std::cout << "      update γ = " << γ << "\n";
    vec rₖ = xₖ - x̂ₖ;

    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // Check stop condition
        real_t εₖ = calc_error_stop_crit(problem, xₖ, x̂ₖ, y, Σ, γ);
        std::cout << "      εₖ    = " << εₖ << "\n";
        if (εₖ <= ε || k == params.max_iter) {
            x     = std::move(x̂ₖ);
            z     = eval_ẑ(problem, xₖ, y, Σ);
            y     = eval_ŷ(problem, xₖ, y, Σ);
            err_z = eval_g(problem, x) - z;

            Stats s;
            s.iterations = k;
            s.converged  = εₖ <= ε;
            s.failed     = false;
            return s;
        }
        if (!std::isfinite(εₖ)) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "inf/NaN"
                         "\x1b[0m"
                      << "\n";

            std::cout << "      "
                      << "xₖ    = " << xₖ.transpose() << "\n";
            std::cout << "      "
                      << "step  = "
                      << gradient_step(problem, xₖ, y, Σ, γ).transpose()
                      << "\n";
            std::cout << "      "
                      << "x̂ₖ    = " << x̂ₖ.transpose() << "\n";
            std::cout << "      "
                      << "rₖ    = " << rₖ.transpose() << "\n";
            std::cout << "      "
                      << "∇ψ(x) = "
                      << eval_grad_ψ(problem, xₖ, y, Σ).transpose() << "\n";
            std::cout << "      "
                      << "ẑ(x)  = " << eval_ẑ(problem, xₖ, y, Σ).transpose()
                      << "\n";
            std::cout << "      "
                      << "ŷ(x)  = " << eval_ŷ(problem, xₖ, y, Σ).transpose()
                      << "\n";

            Stats s;
            s.iterations = k;
            s.converged  = false;
            s.failed     = true;
            return s;
        }

        vec dₖ(n);
        { // Calculate Newton step
            vec rₖ_tmp = rₖ;
            lbfgs.apply(1, rₖ_tmp, dₖ);
        }

        // Line search
        real_t τ = 1;
        vec xₖ₊₁, x̂ₖ₊₁, rₖ₊₁;
        real_t Lₖ₊₁, σₖ₊₁, γₖ₊₁;
        do {
            Lₖ₊₁ = L;
            σₖ₊₁ = σ;
            γₖ₊₁ = γ;

            // Calculate xₖ₊₁
            xₖ₊₁ = xₖ - (1 - τ) * rₖ - τ * dₖ;
            x̂ₖ₊₁ = T_γ(problem, xₖ₊₁, y, Σ, γ);
            while (lipschitz_check(problem, xₖ₊₁, x̂ₖ₊₁, y, Σ, γₖ₊₁, Lₖ₊₁) &&
                   γₖ₊₁ != 0) {
                Lₖ₊₁ *= 2;
                σₖ₊₁ /= 2;
                γₖ₊₁ /= 2;

                // Calculate x̂ₖ (with new step size)
                x̂ₖ₊₁ = T_γ(problem, xₖ₊₁, y, Σ, γₖ₊₁);
            }
            rₖ₊₁ = xₖ₊₁ - x̂ₖ₊₁;

            // Update τ
            τ /= 2;
        } while (eval_φ(problem, xₖ₊₁, rₖ₊₁, y, Σ, γₖ₊₁) >
                     eval_φ(problem, xₖ, rₖ, y, Σ, γ) -
                         σ / (γ * γ) * rₖ.squaredNorm() &&
                 τ >= params.τ_min);

        if (τ < params.τ_min) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "Warning: Line search failed"
                         "\x1b[0m"
                      << "\n";
            xₖ₊₁ = x̂ₖ;
            rₖ₊₁ = xₖ₊₁ - x̂ₖ₊₁;
            Lₖ₊₁ = L;
            σₖ₊₁ = σ;
            γₖ₊₁ = γ; // I know ...
        }

        // Update L-BFGS
        lbfgs.update(xₖ₊₁ - xₖ, rₖ₊₁ - rₖ);

        // Advance step
        xₖ = std::move(xₖ₊₁);
        rₖ = std::move(rₖ₊₁);
        x̂ₖ = std::move(x̂ₖ₊₁);
        L  = Lₖ₊₁;
        σ  = σₖ₊₁;
        γ  = γₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace pa_ref
