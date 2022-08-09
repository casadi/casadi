#pragma once

#include <alpaqa/inner/src/panoc-helpers.tpp>
#include <alpaqa/inner/structured-panoc.hpp>
#include <alpaqa/util/max-history.hpp>

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
using std::chrono::microseconds;

template <Config Conf>
void StructuredPANOCLBFGSSolver<Conf>::compute_quasi_newton_step(
    /// [in]    Solver parameters
    const Params &params,
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Step size
    real_t γₖ,
    /// [in]    Current iterate
    crvec xₖ,
    /// [in]    Lagrange multipliers
    crvec y,
    /// [in]    Penalty weights
    crvec Σ,
    /// [in]    Gradient at current iterate
    crvec grad_ψₖ,
    /// [in]    Projected gradient
    crvec pₖ,
    /// [out]   Quasi-Newton step
    rvec qₖ,
    /// [out]   Indices of active constraints, size at least n
    indexstdvec &J,
    /// [out]   Hessian-vector product of active constraints
    rvec HqK,
    /// [inout] L-BFGS estimate to apply
    LBFGS &lbfgs,
    ///         Dimension n
    rvec work_n,
    ///         Dimension n
    rvec work_n2,
    ///         Dimension m
    rvec work_m) {

    auto n = problem.n, m = problem.m;
    J.clear();
    // Find inactive indices J
    for (index_t i = 0; i < n; ++i) {
        real_t gd = xₖ(i) - γₖ * grad_ψₖ(i);
        if (gd <= problem.get_C().lowerbound(i)) {        // i ∊ J̲ ⊆ K
            qₖ(i) = pₖ(i);                                //
        } else if (problem.get_C().upperbound(i) <= gd) { // i ∊ J̅ ⊆ K
            qₖ(i) = pₖ(i);                                //
        } else {                                          // i ∊ J
            J.push_back(i);
            qₖ(i) = params.hessian_vec ? 0 : -grad_ψₖ(i);
        }
    }

    if (not J.empty()) {     // There are inactive indices J
        if (J.size() == n) { // There are no active indices K
            qₖ = -grad_ψₖ;
        } else if (params.hessian_vec) { // There are active indices K
            if (params.hessian_vec_finite_differences) {
                Helpers::calc_augmented_lagrangian_hessian_prod_fd(
                    problem, xₖ, y, Σ, grad_ψₖ, qₖ, HqK, work_n, work_n2,
                    work_m);
            } else {
                problem.eval_hess_L_prod(xₖ, y, qₖ, HqK);
                if (params.full_augmented_hessian) {
                    auto &g = work_m;
                    problem.eval_g(xₖ, g);
                    for (index_t i = 0; i < m; ++i) {
                        real_t ζ      = g(i) + y(i) / Σ(i);
                        bool inactive = problem.get_D().lowerbound(i) < ζ &&
                                        ζ < problem.get_D().upperbound(i);
                        if (not inactive) {
                            problem.eval_grad_gi(xₖ, i, work_n);
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

        real_t stepsize =
            params.lbfgs_stepsize == LBFGSStepSize::BasedOnGradientStepSize
                ? γₖ
                : -1;
        // If all indices are inactive, we can use standard L-BFGS,
        // if there are active indices, we need the specialized version
        // that only applies L-BFGS to the inactive indices
        bool success = lbfgs.apply_masked(qₖ, stepsize, J);
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

template <Config Conf>
auto StructuredPANOCLBFGSSolver<Conf>::operator()(
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

    const auto n = problem.n;
    const auto m = problem.m;

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
        qₖ = vec::Constant(n, NaN<config_t>), // Newton step Hₖ pₖ
        grad_ψₖ(n),                           // ∇ψ(xₖ)
        grad_̂ψₖ(need_grad_̂ψₖ ? n : 0),        // ∇ψ(x̂ₖ)
        grad_ψₙₑₓₜ(n);                        // ∇ψ(xₙₑₓₜ)

    vec work_n(n), work_m(m);

    vec work_n2(n);
    vec HqK(n);

    std::vector<index_t> J;
    J.reserve(n);
    lbfgs.resize(n);

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress   = 0;
    unsigned last_γ_change = 0;

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PANOC (for readability in the main algorithm)
    auto calc_ψ_ŷ = [&problem, &y, &Σ](crvec x, rvec ŷ) {
        return problem.eval_ψ_ŷ(x, y, Σ, ŷ);
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
        Helpers::calc_x̂(problem.get_C(), γ, x, grad_ψ, x̂, p);
    };
    auto proj_grad_step = [&problem](real_t γ, crvec x, crvec grad_ψ) {
        return Helpers::projected_gradient_step(problem.get_C(), γ, x, grad_ψ);
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
    std::array<char, 64> printbuf;
    auto print_real = [&](real_t x) {
        return float_to_str_vw(printbuf, x, params.print_precision);
    };
    auto print_progress = [&](unsigned k, real_t ψₖ, crvec grad_ψₖ,
                              real_t pₖᵀpₖ, real_t γₖ, real_t εₖ) {
        std::cout << "[PANOC] " << std::setw(6) << k
                  << ": ψ = " << print_real(ψₖ)
                  << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())
                  << ", ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ))
                  << ", γ = " << print_real(γₖ) << ", εₖ = " << print_real(εₖ)
                  << "\r\n";
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

    // Calculate x̂₀, p₀ (projected gradient step)
    calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψx̂ₖ        = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t pₖᵀpₖ      = pₖ.squaredNorm();
    // Compute forward-backward envelope
    real_t φₖ   = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
    real_t nmΦₖ = φₖ;
    // History of FPR norms
    MaxHistory<real_t> fpr_buffer(params.fpr_shortcut_accept_factor //
                                      ? params.fpr_shortcut_history
                                      : 0);
    if (params.fpr_shortcut_accept_factor && params.fpr_shortcut_history)
        fpr_buffer.add(proj_grad_step(1, xₖ, grad_ψₖ).norm());

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Heuristically try to increase step size
        bool force_qub_check = false;
        bool skip_qub_check  = false;
        if (k > 0 && params.hessian_step_size_heuristic > 0 &&
            k - last_γ_change > params.hessian_step_size_heuristic) {
            Helpers::calc_augmented_lagrangian_hessian_prod_fd(
                problem, xₖ, y, Σ, grad_ψₖ, grad_ψₖ, HqK, work_n, work_n2,
                work_m);
            real_t η = grad_ψₖ.squaredNorm() / grad_ψₖ.dot(HqK);
            if (η > γₖ * 10)
                η = γₖ * 10;
            if (η > γₖ) {
                real_t Lₙₑₓₜ = params.Lipschitz.Lγ_factor / η;
                Lₙₑₓₜ = std::clamp(Lₙₑₓₜ, params.L_min, params.L_max);
                real_t γₙₑₓₜ = params.Lipschitz.Lγ_factor / Lₙₑₓₜ;
                // Calculate x̂ₖ, pₖ (projected gradient step)
                calc_x̂(γₙₑₓₜ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₙₑₓₜ, pₙₑₓₜ);
                // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
                real_t grad_ψₙₑₓₜᵀpₙₑₓₜ = grad_ψₖ.dot(pₙₑₓₜ);
                real_t pₙₑₓₜᵀpₙₑₓₜ      = pₙₑₓₜ.squaredNorm();
                // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
                real_t ψx̂ₙₑₓₜ = calc_ψ_ŷ(x̂ₙₑₓₜ, /* in ⟹ out */ ŷx̂ₙₑₓₜ);
                // Check quadratic upper bound
                real_t margin = (1 + std::abs(ψₖ)) *
                                params.quadratic_upperbound_tolerance_factor;
                constexpr bool always_accept = false;
                real_t rhs =
                    grad_ψₙₑₓₜᵀpₙₑₓₜ + real_t(0.5) * Lₙₑₓₜ * pₙₑₓₜᵀpₙₑₓₜ;
                if (always_accept || ψx̂ₙₑₓₜ - ψₖ <= rhs + margin) {
                    // If QUB is satisfied, accept new, larger step
                    Lₖ = Lₙₑₓₜ;
                    γₖ = γₙₑₓₜ;
                    x̂ₖ.swap(x̂ₙₑₓₜ);
                    pₖ.swap(pₙₑₓₜ);
                    grad_ψₖᵀpₖ = grad_ψₙₑₓₜᵀpₙₑₓₜ;
                    pₖᵀpₖ      = pₙₑₓₜᵀpₙₑₓₜ;
                    ψx̂ₖ        = ψx̂ₙₑₓₜ;
                    // Recompute forward-backward envelope
                    φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
                    // No need to re-check later
                    skip_qub_check  = not always_accept;
                    force_qub_check = always_accept;
                    last_γ_change   = k;
                } else {
                    // nothing is updated (except x̂ₙₑₓₜ, pₙₑₓₜ, ŷx̂ₙₑₓₜ, which will
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
        real_t εₖ =
            Helpers::calc_error_stop_crit(problem.get_C(), params.stop_crit, pₖ,
                                          γₖ, xₖ, x̂ₖ, ŷx̂ₖ, grad_ψₖ, grad_̂ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, γₖ, εₖ);
        if (progress_cb) {
            Eigen::Map<indexvec> mJ{J.data(), Eigen::Index(J.size())};
            progress_cb({k,  xₖ, pₖ,      pₖᵀpₖ, x̂ₖ,      qₖ,    mJ,
                         φₖ, ψₖ, grad_ψₖ, ψx̂ₖ,   grad_̂ψₖ, Lₖ,    γₖ,
                         τ,  εₖ, Σ,       y,     problem, params});
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
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = stop_status;
            return s;
        }

        // Calculate Newton step -----------------------------------------------

        if (k > 0) { // No L-BFGS estimate on first iteration → no Newton step
            compute_quasi_newton_step(params, problem, γₖ, xₖ, y, Σ, grad_ψₖ,
                                      pₖ, qₖ, J, HqK, lbfgs, work_n, work_n2,
                                      work_m);
        }

        // Line search initialization ------------------------------------------
        τ                = 1;
        real_t σₖγₖpₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₙₑₓₜ, ψₙₑₓₜ, ψx̂ₙₑₓₜ, grad_ψₙₑₓₜᵀpₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ;
        real_t Lₙₑₓₜ, γₙₑₓₜ;
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
        unsigned last_γ_changeₙₑₓₜ;
        bool fpr_shortcut = false;
        do {
            last_γ_changeₙₑₓₜ = last_γ_change;
            Lₙₑₓₜ             = Lₖ;
            γₙₑₓₜ             = γₖ;

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
                real_t old_γₙₑₓₜ =
                    descent_lemma(xₙₑₓₜ, ψₙₑₓₜ, grad_ψₙₑₓₜ,
                                  /* in ⟹ out */ x̂ₙₑₓₜ, pₙₑₓₜ, ŷx̂ₙₑₓₜ,
                                  /* inout */ ψx̂ₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ,
                                  grad_ψₙₑₓₜᵀpₙₑₓₜ, Lₙₑₓₜ, γₙₑₓₜ);
                if (old_γₙₑₓₜ != γₙₑₓₜ)
                    last_γ_changeₙₑₓₜ = k;
            }

            // Compute forward-backward envelope
            φₙₑₓₜ = ψₙₑₓₜ + 1 / (2 * γₙₑₓₜ) * pₙₑₓₜᵀpₙₑₓₜ + grad_ψₙₑₓₜᵀpₙₑₓₜ;
            // Compute line search condition
            ls_cond = φₙₑₓₜ - (nmΦₖ - σₖγₖpₖᵀpₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= real_t(0.5) * (1 / γₙₑₓₜ - 1 / γₖ) * pₙₑₓₜᵀpₙₑₓₜ_ₖ;

            real_t fpr1_new = params.fpr_shortcut_accept_factor && τ == 1
                                  ? proj_grad_step(1, xₙₑₓₜ, grad_ψₙₑₓₜ).norm()
                                  : inf<config_t>;

            fpr_shortcut =
                fpr1_new < params.fpr_shortcut_accept_factor * fpr_buffer.max();
            if (fpr_shortcut) {
                fpr_buffer.add(fpr1_new);
                ++s.fpr_shortcuts;
            }

            τ /= 2;
        } while (ls_cond > margin && !fpr_shortcut && τ >= params.τ_min);

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
            no_progress = xₖ == xₙₑₓₜ ? no_progress + 1 : 0;

        // Update L-BFGS
        const bool force = true;
        s.lbfgs_rejected += not lbfgs.update(xₖ, xₙₑₓₜ, grad_ψₖ, grad_ψₙₑₓₜ,
                                             LBFGS::Sign::Positive, force);

        // Advance step --------------------------------------------------------
        last_γ_change = last_γ_changeₙₑₓₜ;
        Lₖ            = Lₙₑₓₜ;
        γₖ            = γₙₑₓₜ;

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
