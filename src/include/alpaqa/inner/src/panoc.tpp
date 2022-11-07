#pragma once

#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/src/panoc-helpers.tpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/quadmath/quadmath-print.hpp>
#include <alpaqa/util/src/print.tpp>

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::nanoseconds;

template <class DirectionProviderT>
std::string PANOCSolver<DirectionProviderT>::get_name() const {
    return "PANOCSolver<" + direction_provider.get_name() + ">";
}

template <class DirectionProviderT>
auto PANOCSolver<DirectionProviderT>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Constraint weights @f$ \Sigma @f$
    crvec Σ,
    /// [in]    Tolerance @f$ \varepsilon @f$
    real_t ε,
    /// [in]    Overwrite @p x, @p y and @p err_z even if not converged
    bool always_overwrite_results,
    /// [inout] Decision variable @f$ x @f$
    rvec x,
    /// [inout] Lagrange multipliers @f$ y @f$
    rvec y,
    /// [out]   Slack variable error @f$ g(x) - z @f$
    rvec err_z) -> Stats {

    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.get_n();
    const auto m = problem.get_m();

    // Allocate vectors, init L-BFGS -------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    bool need_grad_̂ψₖ = Helpers::stop_crit_requires_grad_ψx̂(params.stop_crit);

    vec xₖ = x,    // Value of x at the beginning of the iteration
        x̂ₖ(n),     // Value of x after a projected gradient step
        xₙₑₓₜ(n),  // xₖ for next iteration
        x̂ₙₑₓₜ(n),  // x̂ₖ for next iteration
        ŷx̂ₖ(m),    // ŷ(x̂ₖ) = Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂ₙₑₓₜ(m), // ŷ(x̂ₖ) for next iteration
        pₖ(n),     // Projected gradient step pₖ = x̂ₖ - xₖ
        pₙₑₓₜ(n), // Projected gradient step pₙₑₓₜ = x̂ₙₑₓₜ - xₙₑₓₜ
        qₖ(n),                         // Newton step Hₖ pₖ
        grad_ψₖ(n),                    // ∇ψ(xₖ)
        grad_̂ψₖ(need_grad_̂ψₖ ? n : 0), // ∇ψ(x̂ₖ)
        grad_ψₙₑₓₜ(n);                 // ∇ψ(xₙₑₓₜ)

    vec work_n(n), work_m(m);

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PANOC (for readability in the main algorithm)
    auto calc_ψ_ŷ = [&problem, &y, &Σ](crvec x, rvec ŷ) {
        return problem.eval_ψ(x, y, Σ, ŷ);
    };
    auto calc_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](crvec x,
                                                              rvec grad_ψ) {
        return problem.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
    };
    auto calc_grad_ψ_from_ŷ = [&problem, &work_n](crvec x, crvec ŷ,
                                                  rvec grad_ψ) {
        problem.eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work_n);
    };
    auto calc_x̂ = [&problem](real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) {
        problem.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
    };
    auto calc_err_z = [&problem, &y, &Σ](crvec x̂, rvec err_z) {
        Helpers::calc_err_z(problem, x̂, y, Σ, err_z);
    };
    auto descent_lemma = [this, &problem, &y,
                          &Σ](crvec xₖ, real_t ψₖ, crvec grad_ψₖ, rvec x̂ₖ,
                              rvec pₖ, rvec ŷx̂ₖ, real_t &ψx̂ₖ, real_t &pₖᵀpₖ,
                              real_t &grad_ψₖᵀpₖ, real_t &Lₖ, real_t &γₖ) {
        return Helpers::descent_lemma(
            problem, params.quadratic_upperbound_tolerance_factor, params.L_max,
            xₖ, ψₖ, grad_ψₖ, y, Σ, x̂ₖ, pₖ, ŷx̂ₖ, ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ, γₖ);
    };
    std::array<char, 64> print_buf;
    auto print_real = [&](real_t x) {
        return float_to_str_vw(print_buf, x, params.print_precision);
    };
    auto print_real3 = [&](real_t x) {
        return float_to_str_vw(print_buf, x, 3);
    };
    auto print_progress = [&](unsigned k, real_t ψₖ, crvec grad_ψₖ,
                              real_t pₖᵀpₖ, crvec qₖ, real_t γₖ, real_t τₖ,
                              real_t εₖ) {
        std::cout << "[PANOC] " << std::setw(6) << k
                  << ": ψ = " << print_real(ψₖ)
                  << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())
                  << ", ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ))
                  << ", ‖q‖ = " << print_real(qₖ.norm())
                  << ", γ = " << print_real(γₖ) << ", τ = " << print_real3(τₖ)
                  << ", εₖ = " << print_real(εₖ)
                  << std::endl; // Flush for Python buffering
    };

    // Estimate Lipschitz constant ---------------------------------------------

    real_t ψₖ, Lₖ;
    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        Lₖ = Helpers::initial_lipschitz_estimate(
            problem, xₖ, y, Σ, params.Lipschitz.ε, params.Lipschitz.δ,
            params.L_min, params.L_max,
            /* in ⟹ out */ ψₖ, grad_ψₖ, x̂ₖ, grad_ψₙₑₓₜ, work_n, work_m);
    }
    // Initial Lipschitz constant provided by the user
    else {
        Lₖ = params.Lipschitz.L_0;
        // Calculate ψ(xₖ), ∇ψ(x₀)
        ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);
    }
    if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t τ  = NaN<config_t>;

    // First projected gradient step -------------------------------------------

    // Calculate x̂_0, p_0 (projected gradient step)
    calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψx̂ₖ        = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t pₖᵀpₖ      = pₖ.squaredNorm();
    // Compute forward-backward envelope
    real_t φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;

    // Make sure that we don't allocate any memory in the inner loop
    ScopedMallocBlocker mb;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Quadratic upper bound -----------------------------------------------
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            // Decrease step size until quadratic upper bound is satisfied
            real_t old_γₖ =
                descent_lemma(xₖ, ψₖ, grad_ψₖ,
                              /* in ⟹ out */ x̂ₖ, pₖ, ŷx̂ₖ,
                              /* inout */ ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ, γₖ);
            if (k > 0 && γₖ != old_γₖ) { // Flush L-BFGS if γ changed
                direction_provider.changed_γ(γₖ, old_γₖ);
            } else if (k == 0) { // Initialize L-BFGS
                ScopedMallocAllower ma;
                direction_provider.initialize(problem, y, Σ, γₖ, xₖ, x̂ₖ, pₖ,
                                              grad_ψₖ);
            }
            if (γₖ != old_γₖ)
                φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
        }
        // Calculate ∇ψ(x̂ₖ)
        if (need_grad_̂ψₖ)
            calc_grad_ψ_from_ŷ(x̂ₖ, ŷx̂ₖ, /* in ⟹ out */ grad_̂ψₖ);

        // Check stop condition ------------------------------------------------
        real_t εₖ = Helpers::calc_error_stop_crit(problem, params.stop_crit, pₖ,
                                                  γₖ, xₖ, x̂ₖ, ŷx̂ₖ, grad_ψₖ,
                                                  grad_̂ψₖ, work_n, pₙₑₓₜ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, qₖ, γₖ, τ, εₖ);
        if (progress_cb) {
            ScopedMallocAllower ma;
            progress_cb({k, xₖ, pₖ, pₖᵀpₖ, x̂ₖ, φₖ, ψₖ, grad_ψₖ, ψx̂ₖ, grad_̂ψₖ,
                         Lₖ, γₖ, τ, εₖ, Σ, y, problem, params});
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status  = Helpers::check_all_stop_conditions(
             params, time_elapsed, k, stop_signal, ε, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops of PANOC.
            // TODO: move the computation of ẑ and g(x) to ALM?
            if (stop_status == SolverStatus::Converged ||
                stop_status == SolverStatus::Interrupted ||
                always_overwrite_results) {
                calc_err_z(x̂ₖ, /* in ⟹ out */ err_z);
                x = std::move(x̂ₖ);
                y = std::move(ŷx̂ₖ);
            }
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.status       = stop_status;
            s.final_γ      = γₖ;
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------
        if (k > 0)
            direction_provider.apply(γₖ, xₖ, x̂ₖ, pₖ, grad_ψₖ,
                                     /* in ⟹ out */ qₖ);

        // Line search initialization ------------------------------------------
        τ                = 1;
        real_t σₖγₖpₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₙₑₓₜ, ψₙₑₓₜ, ψx̂ₙₑₓₜ, grad_ψₙₑₓₜᵀpₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ;
        real_t Lₙₑₓₜ, γₙₑₓₜ;
        real_t ls_cond;
        // TODO: make separate parameter
        real_t margin =
            (1 + std::abs(φₖ)) * params.quadratic_upperbound_tolerance_factor;

        // Make sure quasi-Newton step is valid
        if (k == 0) {
            τ = 0; // Always use prox step on first iteration
        } else if (not qₖ.allFinite()) {
            τ = 0;
            ++s.lbfgs_failures;
            direction_provider.reset(); // Is there anything else we can do?
        }

        // Line search loop ----------------------------------------------------
        do {
            Lₙₑₓₜ = Lₖ;
            γₙₑₓₜ = γₖ;

            // Calculate xₙₑₓₜ
            if (τ / 2 < params.τ_min) { // line search failed
                xₙₑₓₜ.swap(x̂ₖ);         // → safe prox step
                ψₙₑₓₜ = ψx̂ₖ;
                if (need_grad_̂ψₖ)
                    grad_ψₙₑₓₜ.swap(grad_̂ψₖ);
                else
                    calc_grad_ψ_from_ŷ(xₙₑₓₜ, ŷx̂ₖ, /* in ⟹ out */ grad_ψₙₑₓₜ);
            } else {        // line search didn't fail (yet)
                if (τ == 1) // → faster quasi-Newton step
                    xₙₑₓₜ = xₖ + qₖ;
                else
                    xₙₑₓₜ = xₖ + (1 - τ) * pₖ + τ * qₖ;
                // Calculate ψ(xₙₑₓₜ), ∇ψ(xₙₑₓₜ)
                ψₙₑₓₜ = calc_ψ_grad_ψ(xₙₑₓₜ, /* in ⟹ out */ grad_ψₙₑₓₜ);
            }

            // Calculate x̂ₙₑₓₜ, pₙₑₓₜ (projected gradient step in xₙₑₓₜ)
            calc_x̂(γₙₑₓₜ, xₙₑₓₜ, grad_ψₙₑₓₜ, /* in ⟹ out */ x̂ₙₑₓₜ, pₙₑₓₜ);
            // Calculate ψ(x̂ₙₑₓₜ) and ŷ(x̂ₙₑₓₜ)
            ψx̂ₙₑₓₜ = calc_ψ_ŷ(x̂ₙₑₓₜ, /* in ⟹ out */ ŷx̂ₙₑₓₜ);

            // Quadratic upper bound -------------------------------------------
            grad_ψₙₑₓₜᵀpₙₑₓₜ = grad_ψₙₑₓₜ.dot(pₙₑₓₜ);
            pₙₑₓₜᵀpₙₑₓₜ      = pₙₑₓₜ.squaredNorm();
            real_t pₙₑₓₜᵀpₙₑₓₜ_ₖ = pₙₑₓₜᵀpₙₑₓₜ; // prox step with step size γₖ

            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                (void)descent_lemma(xₙₑₓₜ, ψₙₑₓₜ, grad_ψₙₑₓₜ,
                                    /* in ⟹ out */ x̂ₙₑₓₜ, pₙₑₓₜ, ŷx̂ₙₑₓₜ,
                                    /* inout */ ψx̂ₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ,
                                    grad_ψₙₑₓₜᵀpₙₑₓₜ, Lₙₑₓₜ, γₙₑₓₜ);
            }

            // Compute forward-backward envelope
            φₙₑₓₜ = ψₙₑₓₜ + 1 / (2 * γₙₑₓₜ) * pₙₑₓₜᵀpₙₑₓₜ + grad_ψₙₑₓₜᵀpₙₑₓₜ;
            // Compute line search condition
            ls_cond = φₙₑₓₜ - (φₖ - σₖγₖpₖᵀpₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= real_t(0.5) * (1 / γₙₑₓₜ - 1 / γₖ) * pₙₑₓₜᵀpₙₑₓₜ_ₖ;

            τ /= 2;
        } while (ls_cond > margin && τ >= params.τ_min);

        // If τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
            τ = 0;
        }
        τ *= 2; // restore to the value that was actually accepted
        if (k != 0) {
            s.count_τ += 1;
            s.sum_τ += τ;
            s.τ_1_accepted += τ == 1;
        }

        // Update L-BFGS -------------------------------------------------------
        if (γₖ != γₙₑₓₜ) // Flush L-BFGS if γ changed
            direction_provider.changed_γ(γₙₑₓₜ, γₖ);

        s.lbfgs_rejected += not direction_provider.update(γₖ, γₙₑₓₜ, xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ,
                                                          grad_ψₖ, grad_ψₙₑₓₜ);

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = xₖ == xₙₑₓₜ ? no_progress + 1 : 0;

        // Advance step --------------------------------------------------------
        Lₖ = Lₙₑₓₜ;
        γₖ = γₙₑₓₜ;

        ψₖ  = ψₙₑₓₜ;
        ψx̂ₖ = ψx̂ₙₑₓₜ;
        φₖ  = φₙₑₓₜ;

        xₖ.swap(xₙₑₓₜ);
        x̂ₖ.swap(x̂ₙₑₓₜ);
        ŷx̂ₖ.swap(ŷx̂ₙₑₓₜ);
        pₖ.swap(pₙₑₓₜ);
        grad_ψₖ.swap(grad_ψₙₑₓₜ);
        grad_ψₖᵀpₖ = grad_ψₙₑₓₜᵀpₙₑₓₜ;
        pₖᵀpₖ      = pₙₑₓₜᵀpₙₑₓₜ;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa
