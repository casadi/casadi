#pragma once

#include <panoc-alm/inner/decl/panoc.hpp>
#include <panoc-alm/inner/detail/panoc-helpers.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

template <class DirectionProviderT>
typename PANOCSolver<DirectionProviderT>::Stats
PANOCSolver<DirectionProviderT>::operator()(
    const Problem &problem, ///< [in]    Problem description
    const vec &Σ,           ///< [in]    Constraint weights @f$ \Sigma @f$
    real_t ε,               ///< [in]    Tolerance @f$ \epsilon @f$
    vec &x,                 ///< [inout] Decision variable @f$ x @f$
    vec &y,                 ///< [inout] Lagrange multipliers @f$ y @f$
    vec &err_z              ///< [out]   Slack variable error @f$ g(x) - z @f$
) {
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;
    const auto m = problem.m;

    // Allocate vectors, init L-BFGS -------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    vec xₖ = x,       // Value of x at the beginning of the iteration
        x̂ₖ(n),        // Value of x after a projected gradient step
        xₖ₊₁(n),      // xₖ for next iteration
        x̂ₖ₊₁(n),      // x̂ₖ for next iteration
        ŷx̂ₖ(m),       // ŷ(x̂ₖ) = Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂ₖ₊₁(m),     // ŷ(x̂ₖ) for next iteration
        pₖ(n),        // pₖ = x̂ₖ - xₖ
        pₖ₊₁(n),      // pₖ₊₁ = x̂ₖ₊₁ - xₖ₊₁
        qₖ(n),        // Newton step Hₖ pₖ
        grad_ψₖ(n),   // ∇ψ(xₖ)
        grad_̂ψₖ(n),   // ∇ψ(x̂ₖ)
        grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)

    vec work_n(n), work_m(m);
    direction_provider.resize(n, params.lbfgs_mem);

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PANOC (for readability in the main algorithm)
    auto calc_ψ_ŷ = [&problem, &y, &Σ](const vec &x, vec &ŷ) {
        return detail::calc_ψ_ŷ(problem, x, y, Σ, ŷ);
    };
    auto calc_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](const vec &x,
                                                              vec &grad_ψ) {
        return detail::calc_ψ_grad_ψ(problem, x, y, Σ, grad_ψ, work_n, work_m);
    };
    auto calc_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](const vec &x,
                                                            vec &grad_ψ) {
        detail::calc_grad_ψ(problem, x, y, Σ, grad_ψ, work_n, work_m);
    };
    auto calc_grad_ψ_from_ŷ = [&problem, &work_n](const vec &x, const vec &ŷ,
                                                  vec &grad_ψ) {
        detail::calc_grad_ψ_from_ŷ(problem, x, ŷ, grad_ψ, work_n);
    };
    auto calc_x̂ = [&problem](real_t γ, const vec &x, const vec &grad_ψ, vec &x̂,
                             vec &p) {
        detail::calc_x̂(problem, γ, x, grad_ψ, x̂, p);
    };
    auto calc_err_z = [&problem, &y, &Σ](const vec &x̂, vec &err_z) {
        detail::calc_err_z(problem, x̂, y, Σ, err_z);
    };
    auto print_progress = [&](unsigned k, real_t ψₖ, const vec &grad_ψₖ,
                              real_t norm_sq_pₖ, real_t γₖ, real_t εₖ) {
        std::cout << "[PANOC] " << std::setw(6) << k
                  << ": ψ = " << std::setw(13) << ψₖ
                  << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                  << ", ‖p‖ = " << std::setw(13) << std::sqrt(norm_sq_pₖ)
                  << ", γ = " << std::setw(13) << γₖ
                  << ", εₖ = " << std::setw(13) << εₖ << "\r\n";
    };

    // Estimate Lipschitz constant ---------------------------------------------

    // Finite difference approximation of ∇²ψ in starting point
    vec h(n);
    h = (x * params.Lipschitz.ε).cwiseAbs().cwiseMax(params.Lipschitz.δ);
    x += h;

    // Calculate ∇ψ(x₀ + h)
    calc_grad_ψ(x, /* in ⟹ out */ grad_ψₖ₊₁);

    // Calculate ψ(xₖ), ∇ψ(x₀)
    real_t ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

    // Estimate Lipschitz constant
    real_t Lₖ = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    if (Lₖ < std::numeric_limits<real_t>::epsilon()) {
        Lₖ = std::numeric_limits<real_t>::epsilon();
    } else if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }

    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t σₖ = γₖ * (1 - γₖ * Lₖ) / 2;

    // First projected gradient step -------------------------------------------

    // Calculate x̂₀, p₀ (projected gradient step)
    calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψ̂xₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);

    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t norm_sq_pₖ = pₖ.squaredNorm();

    // Compute forward-backward envelope
    real_t φₖ = ψₖ + 1 / (2 * γₖ) * norm_sq_pₖ + grad_ψₖᵀpₖ;

    unsigned no_progress = 0;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Quadratic upper bound -----------------------------------------------
        // Decrease step size until quadratic upper bound is satisfied
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            real_t margin = 0; // 1e-6 * std::abs(ψₖ); // TODO: Why?
            while (ψ̂xₖ > ψₖ + margin + grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ) {
                Lₖ *= 2;
                σₖ /= 2;
                γₖ /= 2;

                // Flush L-BFGS if γ changed
                if (k > 0)
                    direction_provider.gamma_changed();

                // Calculate x̂ₖ and pₖ (with new step size)
                calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
                // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
                grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
                norm_sq_pₖ = pₖ.squaredNorm();

                // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
                ψ̂xₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);
            }
        }

        // Initialize the L-BFGS
        if (k == 0)
            direction_provider.initialize(xₖ, pₖ, grad_ψₖ, γₖ);

        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(x̂ₖ, ŷx̂ₖ, /* in ⟹ out */ grad_̂ψₖ);

        // Check stop condition ------------------------------------------------
        real_t εₖ = detail::calc_error_stop_crit(pₖ, γₖ, grad_̂ψₖ, grad_ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, norm_sq_pₖ, γₖ, εₖ);

        if (progress_cb)
            progress_cb({k, xₖ, pₖ, norm_sq_pₖ, x̂ₖ, ψₖ, grad_ψₖ, ψ̂xₖ, grad_̂ψₖ,
                         γₖ, εₖ, Σ, y, problem, params});

        auto time_elapsed    = std::chrono::steady_clock::now() - start_time;
        bool out_of_time     = time_elapsed > params.max_time;
        bool out_of_iter     = k == params.max_iter;
        bool interrupted     = stop_signal.stop_requested();
        bool not_finite      = not std::isfinite(εₖ);
        bool conv            = εₖ <= ε;
        bool max_no_progress = no_progress > params.lbfgs_mem;
        bool exit = conv || out_of_iter || out_of_time || not_finite ||
                    interrupted || max_no_progress;
        if (exit) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops of PANOC.
            // TODO: move the computation of ẑ and g(x) to ALM?
            calc_err_z(x̂ₖ, /* in ⟹ out */ err_z);
            x              = std::move(x̂ₖ);
            y              = std::move(ŷx̂ₖ);
            s.iterations   = k; // TODO: what do we count as an iteration?
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = conv              ? SolverStatus::Converged
                             : out_of_time     ? SolverStatus::MaxTime
                             : out_of_iter     ? SolverStatus::MaxIter
                             : not_finite      ? SolverStatus::NotFinite
                             : max_no_progress ? SolverStatus::NoProgress
                                               : SolverStatus::Interrupted;
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------
        if (k > 0) {
            qₖ = pₖ;
            direction_provider.apply(qₖ);
        }

        // Line search initialization ------------------------------------------
        real_t τ            = 1;
        real_t σ_norm_γ⁻¹pₖ = σₖ * norm_sq_pₖ / (γₖ * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, ψ̂xₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, norm_sq_pₖ₊₁;
        real_t Lₖ₊₁, σₖ₊₁, γₖ₊₁;
        real_t ls_cond;

        // Make sure quasi-Newton step is valid
        if (k == 0) {
            τ = 0;
        } else if (not qₖ.allFinite()) {
            τ = 0;
            ++s.lbfgs_failures;
            direction_provider.reset(); // Is there anything else we can do?
        }

        // Line search loop ----------------------------------------------------
        do {
            Lₖ₊₁ = Lₖ;
            σₖ₊₁ = σₖ;
            γₖ₊₁ = γₖ;

            // Calculate xₖ₊₁
            if (τ / 2 < params.τ_min) // line search failed
                xₖ₊₁.swap(x̂ₖ);        // safe prox step
            else                      // line search not failed (yet)
                xₖ₊₁ = xₖ + (1 - τ) * pₖ + τ * qₖ; // faster quasi-Newton step

            // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            ψₖ₊₁ = calc_ψ_grad_ψ(xₖ₊₁, /* in ⟹ out */ grad_ψₖ₊₁);
            // Calculate x̂ₖ₊₁, pₖ₊₁ (projected gradient step)
            calc_x̂(γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁);
            // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
            ψ̂xₖ₊₁ = calc_ψ_ŷ(x̂ₖ₊₁, /* in ⟹ out */ ŷx̂ₖ₊₁);

            // Quadratic upper bound -------------------------------------------
            real_t margin  = 0; // 1e-6 * std::abs(ψₖ₊₁); // TODO: Why?
            grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
            norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
            real_t norm_sq_pₖ₊₁_ₖ = norm_sq_pₖ₊₁; // prox step with step size γₖ
            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                while (ψ̂xₖ₊₁ > ψₖ₊₁ + margin + grad_ψₖ₊₁ᵀpₖ₊₁ +
                                   0.5 * Lₖ₊₁ * norm_sq_pₖ₊₁) {
                    Lₖ₊₁ *= 2;
                    σₖ₊₁ /= 2;
                    γₖ₊₁ /= 2;
                    // Flush L-BFGS if γ changed
                    direction_provider.gamma_changed();

                    // Calculate x̂ₖ₊₁ and pₖ₊₁ (with new step size)
                    calc_x̂(γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁);
                    // Calculate ∇ψ(xₖ₊₁)ᵀpₖ₊₁ and ‖pₖ₊₁‖²
                    grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
                    norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
                    // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
                    ψ̂xₖ₊₁ = calc_ψ_ŷ(x̂ₖ₊₁, /* in ⟹ out */ ŷx̂ₖ₊₁);
                }
            }

            // Compute forward-backward envelope
            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * norm_sq_pₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;

            τ /= 2;

            ls_cond = φₖ₊₁ - (φₖ - σ_norm_γ⁻¹pₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= (0.5 / γₖ₊₁ - 0.5 / γₖ) * norm_sq_pₖ₊₁_ₖ;
        } while (ls_cond > 0 && τ >= params.τ_min);

        // τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
        }

        // Update L-BFGS -------------------------------------------------------
        s.lbfgs_rejected += not direction_provider.update(
            xₖ, xₖ₊₁, pₖ, pₖ₊₁, grad_ψₖ₊₁, problem.C, γₖ₊₁);

        // Check if we made any progress
        if (no_progress > 0 || k % params.lbfgs_mem == 0)
            no_progress = xₖ == xₖ₊₁ ? no_progress + 1 : 0;

        // Advance step --------------------------------------------------------
        Lₖ = Lₖ₊₁;
        σₖ = σₖ₊₁;
        γₖ = γₖ₊₁;

        ψₖ  = ψₖ₊₁;
        ψ̂xₖ = ψ̂xₖ₊₁;
        φₖ  = φₖ₊₁;

        xₖ.swap(xₖ₊₁);
        x̂ₖ.swap(x̂ₖ₊₁);
        ŷx̂ₖ.swap(ŷx̂ₖ₊₁);
        pₖ.swap(pₖ₊₁);
        grad_ψₖ.swap(grad_ψₖ₊₁);
        grad_ψₖᵀpₖ = grad_ψₖ₊₁ᵀpₖ₊₁;
        norm_sq_pₖ = norm_sq_pₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace pa
