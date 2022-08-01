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

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

template <class DirectionProviderT>
std::string PANOCSolver<DirectionProviderT>::get_name() const {
    return "PANOCSolver<" + direction_provider.get_name() + ">";
}

template <class DirectionProviderT>
typename PANOCSolver<DirectionProviderT>::Stats
PANOCSolver<DirectionProviderT>::operator()(
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
    rvec err_z) {

    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;
    const auto m = problem.m;
    auto &&C     = problem.get_C();

    // Allocate vectors, init L-BFGS -------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    bool need_grad_̂ψₖ = Helpers::stop_crit_requires_grad_ψx̂(params.stop_crit);

    vec xₖ = x,   // Value of x at the beginning of the iteration
        x̂ₖ(n),    // Value of x after a projected gradient step
        x_kp1(n),  // xₖ for next iteration
        x̂_kp1(n),  // x̂ₖ for next iteration
        ŷx̂ₖ(m),   // ŷ(x̂ₖ) = Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂_kp1(m), // ŷ(x̂ₖ) for next iteration
        pₖ(n),    // Projected gradient step pₖ = x̂ₖ - xₖ
        p_kp1(n), // Projected gradient step p_kp1 = x̂_kp1 - x_kp1
        qₖ(n),   // Newton step Hₖ pₖ
        grad_ψₖ(n),                    // ∇ψ(xₖ)
        grad_̂ψₖ(need_grad_̂ψₖ ? n : 0), // ∇ψ(x̂ₖ)
        grad_ψ_kp1(n);                  // ∇ψ(x_kp1)

    vec work_n(n), work_m(m);

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PANOC (for readability in the main algorithm)
    auto calc_ψ_ŷ = [&problem, &y, &Σ](crvec x, rvec ŷ) {
        return Helpers::calc_ψ_ŷ(problem, x, y, Σ, ŷ);
    };
    auto calc_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](crvec x,
                                                              rvec grad_ψ) {
        return Helpers::calc_ψ_grad_ψ(problem, x, y, Σ, grad_ψ, work_n, work_m);
    };
    auto calc_grad_ψ_from_ŷ = [&problem, &work_n](crvec x, crvec ŷ,
                                                  rvec grad_ψ) {
        Helpers::calc_grad_ψ_from_ŷ(problem, x, ŷ, grad_ψ, work_n);
    };
    auto calc_x̂ = [&](real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) {
        Helpers::calc_x̂(C, γ, x, grad_ψ, x̂, p);
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
    auto print_progress = [&](unsigned k, real_t ψₖ, crvec grad_ψₖ,
                              real_t pₖᵀpₖ, real_t γₖ, real_t εₖ) {
        std::cout << "[PANOC] " << std::setw(6) << k
                  << ": ψ = " << std::setw(13) << ψₖ
                  << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                  << ", ‖p‖ = " << std::setw(13) << std::sqrt(pₖᵀpₖ)
                  << ", γ = " << std::setw(13) << γₖ
                  << ", εₖ = " << std::setw(13) << εₖ << "\r\n";
    };

    // Estimate Lipschitz constant ---------------------------------------------

    real_t ψₖ, Lₖ;
    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        Lₖ = Helpers::initial_lipschitz_estimate(
            problem, xₖ, y, Σ, params.Lipschitz.ε, params.Lipschitz.δ,
            params.L_min, params.L_max,
            /* in ⟹ out */ ψₖ, grad_ψₖ, x̂ₖ, grad_ψ_kp1, work_n, work_m);
    }
    // Initial Lipschitz constant provided by the user
    else {
        Lₖ = params.Lipschitz.L_0;
        // Calculate ψ(xₖ), ∇ψ(x_0)
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
                direction_provider.initialize(xₖ, x̂ₖ, pₖ, grad_ψₖ);
            }
            if (γₖ != old_γₖ)
                φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
        }
        // Calculate ∇ψ(x̂ₖ)
        if (need_grad_̂ψₖ)
            calc_grad_ψ_from_ŷ(x̂ₖ, ŷx̂ₖ, /* in ⟹ out */ grad_̂ψₖ);

        // Check stop condition ------------------------------------------------
        real_t εₖ = Helpers::calc_error_stop_crit(
            C, params.stop_crit, pₖ, γₖ, xₖ, x̂ₖ, ŷx̂ₖ, grad_ψₖ, grad_̂ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, γₖ, εₖ);
        if (progress_cb)
            progress_cb({k, xₖ, pₖ, pₖᵀpₖ, x̂ₖ, φₖ, ψₖ, grad_ψₖ, ψx̂ₖ, grad_̂ψₖ,
                         Lₖ, γₖ, τ, εₖ, Σ, y, problem, params});

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
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = stop_status;
            s.final_γ      = γₖ;
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------
        real_t step_size =
            params.lbfgs_stepsize == LBFGSStepSize::BasedOnGradientStepSize
                ? real_t(1)
                : real_t(-1);
        if (k > 0)
            direction_provider.apply(xₖ, x̂ₖ, pₖ, step_size,
                                     /* in ⟹ out */ qₖ);

        // Line search initialization ------------------------------------------
        τ                  = 1;
        real_t σₖγₖpₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φ_kp1, ψ_kp1, ψx̂_kp1, grad_ψ_kp1ᵀp_kp1, p_kp1ᵀp_kp1;
        real_t L_kp1, γ_kp1;
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
            L_kp1 = Lₖ;
            γ_kp1 = γₖ;

            // Calculate x_kp1
            if (τ / 2 < params.τ_min) { // line search failed
                x_kp1.swap(x̂ₖ);          // → safe prox step
                ψ_kp1 = ψx̂ₖ;
                if (need_grad_̂ψₖ)
                    grad_ψ_kp1.swap(grad_̂ψₖ);
                else
                    calc_grad_ψ_from_ŷ(x_kp1, ŷx̂ₖ, /* in ⟹ out */ grad_ψ_kp1);
            } else {        // line search didn't fail (yet)
                if (τ == 1) // → faster quasi-Newton step
                    x_kp1 = xₖ + qₖ;
                else
                    x_kp1 = xₖ + (1 - τ) * pₖ + τ * qₖ;
                // Calculate ψ(x_kp1), ∇ψ(x_kp1)
                ψ_kp1 = calc_ψ_grad_ψ(x_kp1, /* in ⟹ out */ grad_ψ_kp1);
            }

            // Calculate x̂_kp1, p_kp1 (projected gradient step in x_kp1)
            calc_x̂(γ_kp1, x_kp1, grad_ψ_kp1, /* in ⟹ out */ x̂_kp1, p_kp1);
            // Calculate ψ(x̂_kp1) and ŷ(x̂_kp1)
            ψx̂_kp1 = calc_ψ_ŷ(x̂_kp1, /* in ⟹ out */ ŷx̂_kp1);

            // Quadratic upper bound -------------------------------------------
            grad_ψ_kp1ᵀp_kp1 = grad_ψ_kp1.dot(p_kp1);
            p_kp1ᵀp_kp1      = p_kp1.squaredNorm();
            real_t p_kp1ᵀp_kp1_ₖ = p_kp1ᵀp_kp1; // prox step with step size γₖ

            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                (void)descent_lemma(x_kp1, ψ_kp1, grad_ψ_kp1,
                                    /* in ⟹ out */ x̂_kp1, p_kp1, ŷx̂_kp1,
                                    /* inout */ ψx̂_kp1, p_kp1ᵀp_kp1,
                                    grad_ψ_kp1ᵀp_kp1, L_kp1, γ_kp1);
            }

            // Compute forward-backward envelope
            φ_kp1 = ψ_kp1 + 1 / (2 * γ_kp1) * p_kp1ᵀp_kp1 + grad_ψ_kp1ᵀp_kp1;
            // Compute line search condition
            ls_cond = φ_kp1 - (φₖ - σₖγₖpₖᵀpₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= real_t(0.5) * (1 / γ_kp1 - 1 / γₖ) * p_kp1ᵀp_kp1_ₖ;

            τ /= 2;
        } while (ls_cond > margin && τ >= params.τ_min);

        // If τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
            τ = 0;
        }
        if (k != 0) {
            s.count_τ += 1;
            s.sum_τ += τ * 2;
            s.τ_1_accepted += τ * 2 == 1;
        }

        // Update L-BFGS -------------------------------------------------------
        if (γₖ != γ_kp1) // Flush L-BFGS if γ changed
            direction_provider.changed_γ(γ_kp1, γₖ);

        s.lbfgs_rejected += not direction_provider.update(xₖ, x_kp1, pₖ, p_kp1,
                                                          grad_ψ_kp1, C, γ_kp1);

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = xₖ == x_kp1 ? no_progress + 1 : 0;

        // Advance step --------------------------------------------------------
        Lₖ = L_kp1;
        γₖ = γ_kp1;

        ψₖ  = ψ_kp1;
        ψx̂ₖ = ψx̂_kp1;
        φₖ  = φ_kp1;

        xₖ.swap(x_kp1);
        x̂ₖ.swap(x̂_kp1);
        ŷx̂ₖ.swap(ŷx̂_kp1);
        pₖ.swap(p_kp1);
        grad_ψₖ.swap(grad_ψ_kp1);
        grad_ψₖᵀpₖ = grad_ψ_kp1ᵀp_kp1;
        pₖᵀpₖ      = p_kp1ᵀp_kp1;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa
