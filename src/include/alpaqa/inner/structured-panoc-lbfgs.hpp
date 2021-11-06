#pragma once

#include <alpaqa/inner/decl/structured-panoc-lbfgs.hpp>
#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

inline StructuredPANOCLBFGSSolver::Stats StructuredPANOCLBFGSSolver::operator()(
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

    // Allocate vectors, init L-BFGS -------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    bool need_grad_̂ψₖ = detail::stop_crit_requires_grad_̂ψₖ(params.stop_crit);

    vec xₖ = x,   // Value of x at the beginning of the iteration
        x̂ₖ(n),    // Value of x after a projected gradient step
        xₖ₊₁(n),  // xₖ for next iteration
        x̂ₖ₊₁(n),  // x̂ₖ for next iteration
        ŷx̂ₖ(m),   // ŷ(x̂ₖ) = Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂ₖ₊₁(m), // ŷ(x̂ₖ) for next iteration
        pₖ(n),    // Projected gradient step pₖ = x̂ₖ - xₖ
        pₖ₊₁(n), // Projected gradient step pₖ₊₁ = x̂ₖ₊₁ - xₖ₊₁
        qₖ(n),   // Newton step Hₖ pₖ
        grad_ψₖ(n),                    // ∇ψ(xₖ)
        grad_̂ψₖ(need_grad_̂ψₖ ? n : 0), // ∇ψ(x̂ₖ)
        grad_ψₖ₊₁(n);                  // ∇ψ(xₖ₊₁)

    vec work_n(n), work_m(m);

    vec work_n2(n);
    vec HqK(n);

    using indexvec = std::vector<vec::Index>;
    indexvec J;
    J.reserve(n);
    lbfgs.resize(n);

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress   = 0;
    unsigned last_γ_change = 0;

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PANOC (for readability in the main algorithm)
    auto calc_ψ_ŷ = [&problem, &y, &Σ](crvec x, rvec ŷ) {
        return detail::calc_ψ_ŷ(problem, x, y, Σ, ŷ);
    };
    auto calc_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](crvec x,
                                                              rvec grad_ψ) {
        return detail::calc_ψ_grad_ψ(problem, x, y, Σ, grad_ψ, work_n, work_m);
    };
    auto calc_grad_ψ_from_ŷ = [&problem, &work_n](crvec x, crvec ŷ,
                                                  rvec grad_ψ) {
        detail::calc_grad_ψ_from_ŷ(problem, x, ŷ, grad_ψ, work_n);
    };
    auto calc_x̂ = [&problem](real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) {
        detail::calc_x̂(problem, γ, x, grad_ψ, x̂, p);
    };
    auto calc_err_z = [&problem, &y, &Σ](crvec x̂, rvec err_z) {
        detail::calc_err_z(problem, x̂, y, Σ, err_z);
    };
    auto descent_lemma = [this, &problem, &y,
                          &Σ](crvec xₖ, real_t ψₖ, crvec grad_ψₖ, rvec x̂ₖ,
                              rvec pₖ, rvec ŷx̂ₖ, real_t &ψx̂ₖ, real_t &pₖᵀpₖ,
                              real_t &grad_ψₖᵀpₖ, real_t &Lₖ, real_t &γₖ) {
        return detail::descent_lemma(
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
    if (params.Lipschitz.L₀ <= 0) {
        Lₖ = detail::initial_lipschitz_estimate(
            problem, xₖ, y, Σ, params.Lipschitz.ε, params.Lipschitz.δ,
            params.L_min, params.L_max,
            /* in ⟹ out */ ψₖ, grad_ψₖ, x̂ₖ, grad_ψₖ₊₁, work_n, work_m);
    }
    // Initial Lipschitz constant provided by the user
    else {
        Lₖ = params.Lipschitz.L₀;
        // Calculate ψ(xₖ), ∇ψ(x₀)
        ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);
    }
    if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t τ  = NaN;

    // First projected gradient step -------------------------------------------

    // Calculate x̂₀, p₀ (projected gradient step)
    calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψx̂ₖ        = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t pₖᵀpₖ      = pₖ.squaredNorm();
    // Compute forward-backward envelope
    real_t φₖ   = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
    real_t nmΦₖ = φₖ;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Heuristically try to increase step size
        bool force_qub_check = false;
        bool skip_qub_check  = false;
        if (k > 0 && params.hessian_step_size_heuristic > 0 &&
            k - last_γ_change > params.hessian_step_size_heuristic) {
            detail::calc_augmented_lagrangian_hessian_prod_fd(
                problem, xₖ, y, Σ, grad_ψₖ, grad_ψₖ, HqK, work_n, work_n2,
                work_m);
            real_t η = grad_ψₖ.squaredNorm() / grad_ψₖ.dot(HqK);
            if (η > γₖ * 10)
                η = γₖ * 10;
            if (η > γₖ) {
                real_t Lₖ₊₁ = params.Lipschitz.Lγ_factor / η;
                Lₖ₊₁ = std::clamp(Lₖ₊₁, params.L_min, params.L_max);
                real_t γₖ₊₁ = params.Lipschitz.Lγ_factor / Lₖ₊₁;
                // Calculate x̂ₖ, pₖ (projected gradient step)
                calc_x̂(γₖ₊₁, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁);
                // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
                real_t grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ.dot(pₖ₊₁);
                real_t pₖ₊₁ᵀpₖ₊₁      = pₖ₊₁.squaredNorm();
                // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
                real_t ψx̂ₖ₊₁ = calc_ψ_ŷ(x̂ₖ₊₁, /* in ⟹ out */ ŷx̂ₖ₊₁);
                // Check quadratic upper bound
                real_t margin = (1 + std::abs(ψₖ)) *
                                params.quadratic_upperbound_tolerance_factor;
                constexpr bool always_accept = false;
                if (always_accept || ψx̂ₖ₊₁ - ψₖ <= grad_ψₖ₊₁ᵀpₖ₊₁ +
                                                       0.5 * Lₖ₊₁ * pₖ₊₁ᵀpₖ₊₁ +
                                                       margin) {
                    // If QUB is satisfied, accept new, larger step
                    Lₖ = Lₖ₊₁;
                    γₖ = γₖ₊₁;
                    x̂ₖ.swap(x̂ₖ₊₁);
                    pₖ.swap(pₖ₊₁);
                    grad_ψₖᵀpₖ = grad_ψₖ₊₁ᵀpₖ₊₁;
                    pₖᵀpₖ      = pₖ₊₁ᵀpₖ₊₁;
                    ψx̂ₖ        = ψx̂ₖ₊₁;
                    // Recompute forward-backward envelope
                    φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
                    // No need to re-check later
                    skip_qub_check  = not always_accept;
                    force_qub_check = always_accept;
                    last_γ_change   = k;
                } else {
                    // nothing is updated (except x̂ₖ₊₁, pₖ₊₁, ŷx̂ₖ₊₁, which will
                    // be overwritten before the next use)
                }
            }
        }

        // Quadratic upper bound -----------------------------------------------
        bool qub_check =
            k == 0 || params.update_lipschitz_in_linesearch == false;
        if ((qub_check && !skip_qub_check) || force_qub_check) {
            // Decrease step size until quadratic upper bound is satisfied
            real_t old_γₖ =
                descent_lemma(xₖ, ψₖ, grad_ψₖ,
                              /* in ⟹ out */ x̂ₖ, pₖ, ŷx̂ₖ,
                              /* inout */ ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ, γₖ);
            if (γₖ != old_γₖ) {
                last_γ_change = k;
                φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
            }
        }
        // Calculate ∇ψ(x̂ₖ)
        if (need_grad_̂ψₖ)
            calc_grad_ψ_from_ŷ(x̂ₖ, ŷx̂ₖ, /* in ⟹ out */ grad_̂ψₖ);

        // Check stop condition ------------------------------------------------
        real_t εₖ = detail::calc_error_stop_crit(
            problem.C, params.stop_crit, pₖ, γₖ, xₖ, x̂ₖ, ŷx̂ₖ, grad_ψₖ, grad_̂ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, γₖ, εₖ);
        if (progress_cb)
            progress_cb({k, xₖ, pₖ, pₖᵀpₖ, x̂ₖ, φₖ, ψₖ, grad_ψₖ, ψx̂ₖ, grad_̂ψₖ,
                         Lₖ, γₖ, τ, εₖ, Σ, y, problem, params});

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status  = detail::check_all_stop_conditions(
            params, time_elapsed, k, stop_signal, ε, εₖ, no_progress);
        if (stop_status != SolverStatus::Unknown) {
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
            return s;
        }

        // Calculate Newton step -----------------------------------------------

        if (k > 0) { // No L-BFGS estimate on first iteration → no Newton step
            J.clear();
            // Find inactive indices J
            for (vec::Index i = 0; i < n; ++i) {
                real_t gd = xₖ(i) - γₖ * grad_ψₖ(i);
                if (gd < problem.C.lowerbound(i)) {        // i ∊ J̲ ⊆ K
                    qₖ(i) = pₖ(i);                         //
                } else if (problem.C.upperbound(i) < gd) { // i ∊ J̅ ⊆ K
                    qₖ(i) = pₖ(i);                         //
                } else {                                   // i ∊ J
                    J.push_back(i);
                    qₖ(i) = params.hessian_vec ? 0 : -grad_ψₖ(i);
                }
            }

            if (not J.empty()) {     // There are inactive indices J
                if (J.size() == n) { // There are no active indices K
                    qₖ = -grad_ψₖ;
                } else if (params.hessian_vec) { // There are active indices K
                    if (params.hessian_vec_finite_differences) {
                        detail::calc_augmented_lagrangian_hessian_prod_fd(
                            problem, xₖ, y, Σ, grad_ψₖ, qₖ, HqK, work_n,
                            work_n2, work_m);
                    } else {
                        problem.hess_L_prod(xₖ, y, qₖ, HqK);
                        if (params.full_augmented_hessian) {
                            auto &g = work_m;
                            problem.g(xₖ, g);
                            for (vec::Index i = 0; i < m; ++i) {
                                real_t ζ      = g(i) + y(i) / Σ(i);
                                bool inactive = problem.D.lowerbound(i) < ζ &&
                                                ζ < problem.D.upperbound(i);
                                if (not inactive) {
                                    problem.grad_gi(xₖ, i, work_n);
                                    auto t = Σ(i) * work_n.dot(qₖ);
                                    // TODO: the dot product is more work than
                                    //       strictly necessary (only over K)
                                    for (auto j : J)
                                        HqK(j) += work_n(j) * t;
                                }
                            }
                        }
                    }

                    for (auto j : J) // Compute right-hand side of 6.1c
                        qₖ(j) = -grad_ψₖ(j) - HqK(j);
                }

                real_t stepsize = params.lbfgs_stepsize ==
                                          LBFGSStepSize::BasedOnGradientStepSize
                                      ? γₖ
                                      : -1;
                // If all indices are inactive, we can use standard L-BFGS,
                // if there are active indices, we need the specialized version
                // that only applies L-BFGS to the inactive indices
                bool success = lbfgs.apply(qₖ, stepsize, J);
                // If L-BFGS application failed, qₖ(J) still contains
                // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
                // A good alternative is to use H₀ = γI as an L-BFGS estimate.
                // This seems to be better than just falling back to a projected
                // gradient step.
                if (not success) {
                    if (J.size() == n)
                        qₖ *= γₖ;
                    else
                        for (auto j : J)
                            qₖ(j) *= γₖ;
                }
            }
        }

        // Line search initialization ------------------------------------------
        τ                  = 1;
        real_t σₖγₖ⁻¹pₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, ψx̂ₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, pₖ₊₁ᵀpₖ₊₁;
        real_t Lₖ₊₁, γₖ₊₁;
        real_t ls_cond;
        real_t w = params.nonmonotone_linesearch;
        nmΦₖ     = k == 0 ? φₖ : w * nmΦₖ + (1 - w) * φₖ;
        // TODO: make separate parameter
        real_t margin =
            (1 + std::abs(nmΦₖ)) * params.quadratic_upperbound_tolerance_factor;

        // Make sure quasi-Newton step is valid
        if (k == 0) {
            τ = 0; // Always use prox step on first iteration
        } else if (not qₖ.allFinite()) {
            τ = 0;
            ++s.lbfgs_failures;
        } else if (J.empty()) {
            τ = 0; // All constraints are active, no Newton step possible
        }

        // Line search loop ----------------------------------------------------
        unsigned last_γ_changeₖ₊₁;
        do {
            last_γ_changeₖ₊₁ = last_γ_change;
            Lₖ₊₁             = Lₖ;
            γₖ₊₁             = γₖ;

            // Calculate xₖ₊₁
            if (τ / 2 < params.τ_min) { // line search failed
                xₖ₊₁.swap(x̂ₖ);          // → safe prox step
                ψₖ₊₁ = ψx̂ₖ;
                if (need_grad_̂ψₖ)
                    grad_ψₖ₊₁.swap(grad_̂ψₖ);
                else
                    calc_grad_ψ_from_ŷ(xₖ₊₁, ŷx̂ₖ, /* in ⟹ out */ grad_ψₖ₊₁);
            } else {        // line search didn't fail (yet)
                if (τ == 1) // → faster quasi-Newton step
                    xₖ₊₁ = xₖ + qₖ;
                else
                    xₖ₊₁ = xₖ + (1 - τ) * pₖ + τ * qₖ;
                // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
                ψₖ₊₁ = calc_ψ_grad_ψ(xₖ₊₁, /* in ⟹ out */ grad_ψₖ₊₁);
            }

            // Calculate x̂ₖ₊₁, pₖ₊₁ (projected gradient step in xₖ₊₁)
            calc_x̂(γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁);
            // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
            ψx̂ₖ₊₁ = calc_ψ_ŷ(x̂ₖ₊₁, /* in ⟹ out */ ŷx̂ₖ₊₁);

            // Quadratic upper bound -------------------------------------------
            grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
            pₖ₊₁ᵀpₖ₊₁      = pₖ₊₁.squaredNorm();
            real_t pₖ₊₁ᵀpₖ₊₁_ₖ = pₖ₊₁ᵀpₖ₊₁; // prox step with step size γₖ

            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                real_t old_γₖ₊₁ = descent_lemma(
                    xₖ₊₁, ψₖ₊₁, grad_ψₖ₊₁,
                    /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁, ŷx̂ₖ₊₁,
                    /* inout */ ψx̂ₖ₊₁, pₖ₊₁ᵀpₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, Lₖ₊₁, γₖ₊₁);
                if (old_γₖ₊₁ != γₖ₊₁)
                    last_γ_changeₖ₊₁ = k;
            }

            // Compute forward-backward envelope
            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * pₖ₊₁ᵀpₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;
            // Compute line search condition
            ls_cond = φₖ₊₁ - (nmΦₖ - σₖγₖ⁻¹pₖᵀpₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= (0.5 / γₖ₊₁ - 0.5 / γₖ) * pₖ₊₁ᵀpₖ₊₁_ₖ;

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
            // std::cout << J.size() << "," << τ*2 << "\r\n";
        }

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = xₖ == xₖ₊₁ ? no_progress + 1 : 0;

        // Update L-BFGS
        const bool force = true;
        s.lbfgs_rejected += not lbfgs.update(xₖ, xₖ₊₁, grad_ψₖ, grad_ψₖ₊₁,
                                             LBFGS::Sign::Positive, force);

        // Advance step --------------------------------------------------------
        last_γ_change = last_γ_changeₖ₊₁;
        Lₖ            = Lₖ₊₁;
        γₖ            = γₖ₊₁;

        ψₖ  = ψₖ₊₁;
        ψx̂ₖ = ψx̂ₖ₊₁;
        φₖ  = φₖ₊₁;

        xₖ.swap(xₖ₊₁);
        x̂ₖ.swap(x̂ₖ₊₁);
        ŷx̂ₖ.swap(ŷx̂ₖ₊₁);
        pₖ.swap(pₖ₊₁);
        grad_ψₖ.swap(grad_ψₖ₊₁);
        grad_ψₖᵀpₖ = grad_ψₖ₊₁ᵀpₖ₊₁;
        pₖᵀpₖ      = pₖ₊₁ᵀpₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa
