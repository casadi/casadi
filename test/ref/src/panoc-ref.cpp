#include <alpaqa-ref/panoc-ref.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace pa_ref {

using alpaqa::dist_squared;
using alpaqa::LBFGS;
using alpaqa::project;

namespace detail {

vec eval_g(const Problem &prob, crvec x) {
    vec g(prob.m);
    prob.g(x, g);
    return g;
}

vec eval_ẑ(const Problem &prob, crvec x, crvec y, crvec Σ) {
    vec g = eval_g(prob, x);
    return project(g + Σ.asDiagonal().inverse() * y, prob.D);
}

vec eval_ŷ(const Problem &prob, crvec x, crvec y, crvec Σ) {
    vec g = eval_g(prob, x);
    vec ẑ = eval_ẑ(prob, x, y, Σ);
    return Σ.asDiagonal() * (g - ẑ) + y;
}

real_t eval_ψ(const Problem &prob, crvec x, crvec y, crvec Σ) {
    vec g = eval_g(prob, x);
    return prob.f(x) +
           0.5 * dist_squared(g + Σ.asDiagonal().inverse() * y, prob.D, Σ);
}

vec eval_grad_ψ(const Problem &prob, crvec x, crvec y, crvec Σ) {
    vec ŷ = eval_ŷ(prob, x, y, Σ);
    vec grad_f(prob.n), grad_gŷ(prob.n);
    prob.grad_f(x, grad_f);
    prob.grad_g_prod(x, ŷ, grad_gŷ);
    return grad_f + grad_gŷ;
}

vec projected_gradient_step(const Problem &prob, crvec x, crvec y, crvec Σ,
                            real_t γ) {
    using binary_real_f = real_t (*)(real_t, real_t);
    return (-γ * eval_grad_ψ(prob, x, y, Σ))
        .binaryExpr(prob.C.lowerbound - x, binary_real_f(std::fmax))
        .binaryExpr(prob.C.upperbound - x, binary_real_f(std::fmin));
}

real_t eval_φ(const Problem &prob, crvec x, crvec y, crvec Σ, real_t γ) {
    vec p = projected_gradient_step(prob, x, y, Σ, γ);
    return eval_ψ(prob, x, y, Σ)                //
           + 1. / (2 * γ) * p.squaredNorm()     //
           + eval_grad_ψ(prob, x, y, Σ).dot(p); //
}

real_t estimate_lipschitz(const Problem &prob, crvec x, crvec y, crvec Σ,
                          const PANOCParams &params) {
    // Estimate Lipschitz constant using finite difference
    vec h = (x * params.Lipschitz.ε).cwiseAbs().cwiseMax(params.Lipschitz.δ);
    // Calculate ∇ψ(x₀ + h) and ∇ψ(x₀)
    vec diff = eval_grad_ψ(prob, x + h, y, Σ) - eval_grad_ψ(prob, x, y, Σ);
    // Estimate Lipschitz constant
    real_t L = diff.norm() / h.norm();
    return L;
}

real_t calc_error_stop_crit(const Problem &prob, crvec xₖ, crvec x̂ₖ, crvec y,
                            crvec Σ, real_t γ) {
    vec p = projected_gradient_step(prob, xₖ, y, Σ, γ);
    return ((1 / γ) * p                    //
            + (eval_grad_ψ(prob, xₖ, y, Σ) //
               - eval_grad_ψ(prob, x̂ₖ, y, Σ)))
        .lpNorm<Eigen::Infinity>();
}

bool lipschitz_check(const Problem &prob, crvec xₖ, crvec x̂ₖ, crvec y, crvec Σ,
                     real_t γ, real_t L, real_t rounding_threshold) {
    real_t ψₖ   = eval_ψ(prob, xₖ, y, Σ);
    real_t ψx̂ₖ  = eval_ψ(prob, x̂ₖ, y, Σ);
    vec grad_ψₖ = eval_grad_ψ(prob, xₖ, y, Σ);
    vec pₖ      = projected_gradient_step(prob, xₖ, y, Σ, γ);
    if (std::abs(grad_ψₖ.dot(pₖ) / ψₖ) <= rounding_threshold)
        return false;
    return ψx̂ₖ > ψₖ + grad_ψₖ.dot(pₖ) + L / 2 * pₖ.squaredNorm();
}

} // namespace detail

using alpaqa::SolverStatus;
using std::chrono::duration_cast;
using std::chrono::microseconds;

PANOCSolver::Stats PANOCSolver::operator()(
    const Problem &problem, ///< [in]    Problem description
    crvec Σ,                ///< [in]    Constraint weights @f$ \Sigma @f$
    real_t ε,               ///< [in]    Tolerance @f$ \varepsilon @f$
    bool
        always_overwrite_results, ///< [in] Overwrite x, y and err_z even if not converged
    rvec x,                       ///< [inout] Decision variable @f$ x @f$
    rvec y,                       ///< [inout] Lagrange multiplier @f$ x @f$
    rvec err_z ///< [out]   Slack variable error @f$ g(x) - z @f$
) {
    using namespace detail;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;

    vec xₖ = x; // Value of x at the beginning of the iteration

    real_t Lₖ = estimate_lipschitz(problem, xₖ, y, Σ, params);
    if (Lₖ < 10 * std::numeric_limits<real_t>::epsilon())
        Lₖ = 10 * std::numeric_limits<real_t>::epsilon();
    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t σₖ = γₖ * (1 - γₖ * Lₖ) / 2;

    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // Calculate x̂ₖ, pₖ (gradient step)
        vec x̂ₖ = xₖ + projected_gradient_step(problem, xₖ, y, Σ, γₖ);

        real_t old_γₖ = γₖ;
        while (lipschitz_check(problem, xₖ, x̂ₖ, y, Σ, γₖ, Lₖ,
                               params.quadratic_upperbound_tolerance_factor)) {
            assert(not params.update_lipschitz_in_linesearch || k == 0);
            Lₖ *= 2;
            σₖ /= 2;
            γₖ /= 2;

            // Calculate x̂ₖ (with new step size)
            x̂ₖ = xₖ + projected_gradient_step(problem, xₖ, y, Σ, γₖ);
        }
        // Flush L-BFGS if γ changed
        if (k > 0 && γₖ != old_γₖ)
            lbfgs.changed_γ(γₖ, old_γₖ);

        // Check stop condition
        real_t εₖ = calc_error_stop_crit(problem, xₖ, x̂ₖ, y, Σ, γₖ);

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        bool out_of_time  = time_elapsed > params.max_time;
        bool out_of_iter  = k == params.max_iter;
        bool interrupted  = stop_signal.load(std::memory_order_relaxed);
        bool not_finite   = not std::isfinite(εₖ);
        bool conv         = εₖ <= ε;
        if (conv || out_of_iter || out_of_time || not_finite || interrupted) {
            if (conv || interrupted || always_overwrite_results) {
                err_z = eval_g(problem, x̂ₖ) - eval_ẑ(problem, x̂ₖ, y, Σ);
                y     = eval_ŷ(problem, x̂ₖ, y, Σ);
                x     = std::move(x̂ₖ);
            }
            s.iterations   = k; // TODO: what do we count as an iteration?
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = conv          ? SolverStatus::Converged
                             : out_of_time ? SolverStatus::MaxTime
                             : out_of_iter ? SolverStatus::MaxIter
                             : not_finite  ? SolverStatus::NotFinite
                                           : SolverStatus::Interrupted;
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------
        vec pₖ = projected_gradient_step(problem, xₖ, y, Σ, γₖ);
        vec qₖ(n);
        real_t step_size =
            params.lbfgs_stepsize == alpaqa::LBFGSStepSize::BasedOnGradientStepSize
                ? 1
                : -1;
        if (k > 0)
            lbfgs.apply(xₖ, x̂ₖ, pₖ, step_size, qₖ);
        else
            // Initialize the L-BFGS
            lbfgs.initialize(xₖ, x̂ₖ, pₖ, eval_grad_ψ(problem, xₖ, y, Σ));

        // Line search
        real_t τ = 1;
        if (k == 0) {
            τ = 0;
        } else if (not qₖ.allFinite()) {
            τ = 0;
            ++s.lbfgs_failures;
            lbfgs.reset();
        }
        real_t Lₖ₊₁, σₖ₊₁, γₖ₊₁;
        vec xₖ₊₁(n), pₖ₊₁(n), x̂ₖ₊₁(n);
        do {
            Lₖ₊₁ = Lₖ;
            σₖ₊₁ = σₖ;
            γₖ₊₁ = γₖ;

            // Calculate xₖ₊₁
            if (τ / 2 >= params.τ_min)
                xₖ₊₁ = xₖ + (1 - τ) * pₖ + τ * qₖ;
            else
                xₖ₊₁ = x̂ₖ;

            pₖ₊₁ = projected_gradient_step(problem, xₖ₊₁, y, Σ, γₖ₊₁);
            x̂ₖ₊₁ = xₖ₊₁ + pₖ₊₁;

            if (params.update_lipschitz_in_linesearch) {
                real_t old_γₖ₊₁ = γₖ₊₁;
                while (lipschitz_check(
                    problem, xₖ₊₁, x̂ₖ₊₁, y, Σ, γₖ₊₁, Lₖ₊₁,
                    params.quadratic_upperbound_tolerance_factor)) {
                    Lₖ₊₁ *= 2;
                    σₖ₊₁ /= 2;
                    γₖ₊₁ /= 2;

                    // Calculate x̂ₖ (with new step size)
                    pₖ₊₁ = projected_gradient_step(problem, xₖ₊₁, y, Σ, γₖ₊₁);
                    x̂ₖ₊₁ = xₖ₊₁ + pₖ₊₁;
                }
                // Flush L-BFGS if γ changed
                if (γₖ₊₁ != old_γₖ₊₁)
                    lbfgs.changed_γ(γₖ₊₁, old_γₖ₊₁);
            }

            // Update τ
            τ /= 2;
        } while (eval_φ(problem, xₖ₊₁, y, Σ, γₖ₊₁) >
                     eval_φ(problem, xₖ, y, Σ, γₖ) -
                         σₖ / (γₖ * γₖ) * pₖ.squaredNorm() &&
                 τ >= params.τ_min);

        // τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
        }

        // Update L-BFGS -------------------------------------------------------
        s.lbfgs_rejected +=
            not lbfgs.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁,
                             eval_grad_ψ(problem, xₖ₊₁, y, Σ), problem.C, γₖ₊₁);

        // Advance step
        xₖ = std::move(xₖ₊₁);
        Lₖ = Lₖ₊₁;
        σₖ = σₖ₊₁;
        γₖ = γₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace pa_ref
