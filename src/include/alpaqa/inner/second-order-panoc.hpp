#pragma once

#include <alpaqa/inner/decl/second-order-panoc.hpp>
#include <alpaqa/inner/detail/panoc-helpers.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <Eigen/Cholesky>

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

inline SecondOrderPANOCSolver::Stats SecondOrderPANOCSolver::operator()(
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

    vec xₖ = x,   // Value of x at the beginning of the iteration
        x̂ₖ(n),    // Value of x after a projected gradient step
        xₖ₊₁(n),  // xₖ for next iteration
        x̂ₖ₊₁(n),  // x̂ₖ for next iteration
        ŷx̂ₖ(m),   // ŷ(x̂ₖ) = Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂ₖ₊₁(m), // ŷ(x̂ₖ) for next iteration
        pₖ(n),    // Projected gradient step pₖ = x̂ₖ - xₖ
        pₖ₊₁(n), // Projected gradient step pₖ₊₁ = x̂ₖ₊₁ - xₖ₊₁
        qₖ(n),   // Newton step Hₖ pₖ
        grad_ψₖ(n),   // ∇ψ(xₖ)
        grad_̂ψₖ(n),   // ∇ψ(x̂ₖ)
        grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)

    vec work_n(n), work_m(m);

    mat H(n, n);    // Hessian of the (augmented) Lagrangian function and ψ
    vec g(m),       // Value of the constraint function
        grad_gi(n); // Gradient of one of the constraints

    using indexvec  = std::vector<vec::Index>;
    using bvec      = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
    bvec inactive_C = bvec::Zero(n);

    indexvec J, K;
    J.reserve(n);
    K.reserve(n);
    vec qJ(n); // Solution of Newton Hessian system
    vec rhs(n);

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

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
            /* in ⟹ out */ ψₖ, grad_ψₖ, x̂ₖ, grad_̂ψₖ, work_n, work_m);
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

    // First projected gradient step -------------------------------------------

    // Calculate x̂₀, p₀ (projected gradient step)
    calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψx̂ₖ        = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷx̂ₖ);
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t pₖᵀpₖ      = pₖ.squaredNorm();
    // Forward-backward envelope
    real_t φₖ;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Quadratic upper bound -----------------------------------------------
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            // Decrease step size until quadratic upper bound is satisfied
            descent_lemma(xₖ, ψₖ, grad_ψₖ,
                          /* in ⟹ out */ x̂ₖ, pₖ, ŷx̂ₖ,
                          /* inout */ ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ, γₖ);
            φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
        }
        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(x̂ₖ, ŷx̂ₖ, /* in ⟹ out */ grad_̂ψₖ);

        // Check stop condition ------------------------------------------------
        real_t εₖ = detail::calc_error_stop_crit(
            problem.C, params.stop_crit, pₖ, γₖ, xₖ, x̂ₖ, ŷx̂ₖ, grad_ψₖ, grad_̂ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, γₖ, εₖ);
        if (progress_cb)
            progress_cb({k, xₖ, pₖ, pₖᵀpₖ, x̂ₖ, ψₖ, grad_ψₖ, ψx̂ₖ, grad_̂ψₖ, Lₖ,
                         γₖ, εₖ, Σ, y, problem, params});

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

        // TODO: all of this is suboptimal and ugly :(

        // TODO: write helper lambda above
        detail::calc_augmented_lagrangian_hessian(problem, xₖ, ŷx̂ₖ, y, Σ, g, H,
                                                  grad_gi);

        K.clear();
        J.clear();
        for (vec::Index i = 0; i < n; ++i) {
            real_t gd_step = xₖ(i) - γₖ * grad_ψₖ(i);
            if (gd_step < problem.C.lowerbound(i)) {
                K.push_back(i);
            } else if (problem.C.upperbound(i) < gd_step) {
                K.push_back(i);
            } else {
                J.push_back(i);
            }
        }
        qₖ                  = pₖ;
        bool newton_success = true;
        if (not J.empty()) {
            // Compute right-hand side of 6.1c
            vec::Index ri = 0;
            for (auto j : J) {
                real_t hess_q = 0;
                for (auto k : K)
                    hess_q += H(j, k) * qₖ(k);
                rhs(ri++) = -grad_ψₖ(j) - hess_q;
            }

            // Permute H to get the xJxJ block in the top left
            auto hess_Ljj = H.block(0, 0, J.size(), J.size());
            vec::Index r  = 0;
            for (auto rj : J) {
                vec::Index c = 0;
                for (auto cj : J) {
                    hess_Ljj(r, c) = H(rj, cj);
                    ++c;
                }
                ++r;
            }
            auto ldl = hess_Ljj.ldlt(); // TODO: does this allocate?
            if (ldl.isPositive()) {
                qJ.topRows(J.size()) = ldl.solve(rhs.topRows(J.size()));
            } else {
                std::cerr << "\x1b[0;31m"
                             "not PD"
                             "\x1b[0m"
                          << std::endl;
                // newton_success = false; // TODO
            }

            r = 0;
            for (auto j : J)
                qₖ(j) = qJ(r++);
        }

        // Line search initialization ------------------------------------------
        real_t τ           = 1;
        real_t σₖγₖ⁻¹pₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, ψx̂ₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, pₖ₊₁ᵀpₖ₊₁;
        real_t Lₖ₊₁, γₖ₊₁;
        real_t ls_cond;

        // Make sure quasi-Newton step is valid
        if (not newton_success || not qₖ.allFinite()) {
            τ = 0;
            ++s.newton_failures;
        }

        // Line search loop ----------------------------------------------------
        do {
            Lₖ₊₁ = Lₖ;
            γₖ₊₁ = γₖ;

            // Calculate xₖ₊₁
            if (τ / 2 < params.τ_min) { // line search failed
                xₖ₊₁.swap(x̂ₖ);          // → safe prox step
                ψₖ₊₁ = ψx̂ₖ;
                grad_ψₖ₊₁.swap(grad_̂ψₖ);
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
                descent_lemma(xₖ₊₁, ψₖ₊₁, grad_ψₖ₊₁,
                              /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁, ŷx̂ₖ₊₁,
                              /* inout */ ψx̂ₖ₊₁, pₖ₊₁ᵀpₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁,
                              Lₖ₊₁, γₖ₊₁);
            }

            // Compute forward-backward envelope
            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * pₖ₊₁ᵀpₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;
            // Compute line search condition
            ls_cond = φₖ₊₁ - (φₖ - σₖγₖ⁻¹pₖᵀpₖ);
            if (params.alternative_linesearch_cond)
                ls_cond -= (0.5 / γₖ₊₁ - 0.5 / γₖ) * pₖ₊₁ᵀpₖ₊₁_ₖ;

            τ /= 2;
        } while (ls_cond > 0 && τ >= params.τ_min);

        // τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
            τ = 0;
        }
        if (k != 0) {
            s.count_τ += 1;
            s.sum_τ += τ * 2;
            s.τ_1_accepted += τ * 2 == 1;
        }

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = xₖ == xₖ₊₁ ? no_progress + 1 : 0;

        // Advance step --------------------------------------------------------
        Lₖ = Lₖ₊₁;
        γₖ = γₖ₊₁;

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
