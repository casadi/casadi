#pragma once

#include <alpaqa/inner/decl/panoc.hpp>
#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/util/alloc.hpp>
#include <alpaqa/util/atomic_stop_signal.hpp>
#include <alpaqa/util/box.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace alpaqa {

namespace detail {

/// Estimate the Lipschitz constant of the gradient @f$ \nabla \psi @f$ using
/// finite differences.
template <class ObjFunT, class ObjGradFunT>
real_t initial_lipschitz_estimate(
    /// [in]    Objective
    const ObjFunT &ψ,
    /// [in]    Gradient of
    const ObjGradFunT &grad_ψ,
    /// [in]    Current iterate @f$ x^k @f$
    crvec xₖ,
    /// [in]    Finite difference step size relative to xₖ
    real_t ε,
    /// [in]    Minimum absolute finite difference step size
    real_t δ,
    /// [in]    Minimum allowed Lipschitz estimate.
    real_t L_min,
    /// [in]    Maximum allowed Lipschitz estimate.
    real_t L_max,
    /// [out]   @f$ \psi(x^k) @f$
    real_t &ψₖ,
    /// [out]   Gradient @f$ \nabla \psi(x^k) @f$
    rvec grad_ψₖ,
    /// [in]
    vec_allocator &alloc) {

    auto h  = (xₖ * ε).cwiseAbs().cwiseMax(δ);
    auto xh = alloc.alloc() = xₖ + h;
    real_t norm_h           = h.norm();
    auto grad_ψxh           = alloc.alloc();
    // Calculate ∇ψ(x₀ + h)
    grad_ψ(xh, grad_ψxh);
    // Calculate ψ(xₖ), ∇ψ(x₀)
    ψₖ = ψ(xₖ);
    grad_ψ(xₖ, grad_ψₖ);

    // Estimate Lipschitz constant using finite differences
    real_t L = (grad_ψxh - grad_ψₖ).norm() / norm_h;
    alloc.free(xh, grad_ψxh);
    return std::clamp(L, L_min, L_max);
}

inline void calc_x̂(real_t γ,     ///< [in]  Step size
                   crvec x,      ///< [in]  Decision variable @f$ x @f$
                   const Box &C, ///< [in]  Box constraints @f$ C @f$
                   crvec grad_ψ, ///< [in]  @f$ \nabla \psi(x^k) @f$
                   rvec x̂, ///< [out] @f$ \hat{x}^k = T_{\gamma^k}(x^k) @f$
                   rvec p  ///< [out] @f$ \hat{x}^k - x^k @f$
) {
    p = projected_gradient_step(C, γ, x, grad_ψ);
    x̂ = x + p;
}

/// Increase the estimate of the Lipschitz constant of the objective gradient
/// and decrease the step size until quadratic upper bound or descent lemma is
/// satisfied:
/// @f[ \psi(x + p) \le \psi(x) + \nabla\psi(x)^\top p + \frac{L}{2} \|p\|^2 @f]
/// The projected gradient iterate @f$ \hat x^k @f$ and step @f$ p^k @f$ are
/// updated with the new step size @f$ \gamma^k @f$, and so are other quantities
/// that depend on them, such as @f$ \nabla\psi(x^k)^\top p^k @f$ and
/// @f$ \|p\|^2 @f$. The intermediate vector @f$ \hat y(x^k) @f$ can be used to
/// compute the gradient @f$ \nabla\psi(\hat x^k) @f$ or to update the Lagrange
/// multipliers.
///
/// @return The original step size, before it was reduced by this function.
template <class ObjFunT>
real_t descent_lemma(
    /// [in]    Objective.
    const ObjFunT &ψ,
    /// [in]    Box constraints.
    const Box &C,
    /// [in]    Tolerance used to ignore rounding errors when the function
    ///         @f$ \psi(x) @f$ is relatively flat or the step size is very
    ///         small, which could cause @f$ \psi(x^k) < \psi(\hat x^k) @f$,
    ///         which is mathematically impossible but could occur in finite
    ///         precision floating point arithmetic.
    real_t rounding_tolerance,
    /// [in]    Maximum allowed Lipschitz constant estimate (prevents infinite
    ///         loop if function or derivatives are discontinuous, and keeps
    ///         step size bounded away from zero).
    real_t L_max,
    /// [in]    Current iterate @f$ x^k @f$
    crvec xₖ,
    /// [in]    Objective function @f$ \psi(x^k) @f$
    real_t ψₖ,
    /// [in]    Gradient of objective @f$ \nabla\psi(x^k) @f$
    crvec grad_ψₖ,
    /// [out]   Projected gradient iterate @f$ \hat x^k @f$
    rvec x̂ₖ,
    /// [out]   Projected gradient step @f$ p^k @f$
    rvec pₖ,
    /// [inout] Objective function @f$ \psi(\hat x^k) @f$
    real_t &ψx̂ₖ,
    /// [inout] Squared norm of the step @f$ \left\| p^k \right\|^2 @f$
    real_t &norm_sq_pₖ,
    /// [inout] Gradient of objective times step @f$ \nabla\psi(x^k)^\top p^k@f$
    real_t &grad_ψₖᵀpₖ,
    /// [inout] Lipschitz constant estimate @f$ L_{\nabla\psi}^k @f$
    real_t &Lₖ,
    /// [inout] Step size @f$ \gamma^k @f$
    real_t &γₖ) {

    real_t old_γₖ = γₖ;
    real_t margin = (1 + std::abs(ψₖ)) * rounding_tolerance;
    while (ψx̂ₖ - ψₖ > grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ + margin) {
        if (not(Lₖ * 2 <= L_max))
            break;

        Lₖ *= 2;
        γₖ /= 2;

        // Calculate x̂ₖ and pₖ (with new step size)
        calc_x̂(γₖ, xₖ, C, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
        // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
        grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
        norm_sq_pₖ = pₖ.squaredNorm();

        // Calculate ψ(x̂ₖ)
        ψx̂ₖ = ψ(x̂ₖ);
    }
    return old_γₖ;
}

inline void print_progress(unsigned k, real_t ψₖ, crvec grad_ψₖ, real_t pₖᵀpₖ,
                           real_t γₖ, real_t εₖ) {
    std::cout << "[PANOC] " << std::setw(6) << k                 //
              << ": ψ = " << std::setw(13) << ψₖ                 //
              << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()  //
              << ", ‖p‖ = " << std::setw(13) << std::sqrt(pₖᵀpₖ) //
              << ", γ = " << std::setw(13) << γₖ                 //
              << ", εₖ = " << std::setw(13) << εₖ << "\r\n";
};

} // namespace detail

using std::chrono::duration_cast;
using std::chrono::microseconds;

constexpr size_t panoc_min_alloc_size = 10;

template <class ObjFunT, class ObjGradFunT, class DirectionT>
inline PANOCStats panoc_impl(ObjFunT &ψ, ObjGradFunT &grad_ψ, const Box &C,
                             rvec x, real_t ε, const PANOCParams &params,
                             vec_allocator &alloc,
                             DirectionT &direction_provider) {
    auto start_time = std::chrono::steady_clock::now();
    PANOCStats s;

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Future parameters?
    AtomicStopSignal stop_signal;

    // Estimate Lipschitz constant ---------------------------------------------

    real_t ψₖ, Lₖ;
    auto grad_ψₖ = alloc.alloc();
    auto xₖ = alloc.alloc() = x;
    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L₀ <= 0) {
        Lₖ = detail::initial_lipschitz_estimate(
            ψ, grad_ψ, xₖ, params.Lipschitz.ε, params.Lipschitz.δ, params.L_min,
            params.L_max, ψₖ, grad_ψₖ, alloc);
    }
    // Initial Lipschitz constant provided by the user
    else {
        Lₖ = params.Lipschitz.L₀;
        // Calculate ψ(xₖ), ∇ψ(x₀)
        ψₖ = ψ(xₖ);
        grad_ψ(xₖ, grad_ψₖ);
    }
    if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t τ  = NaN;

    // First projected gradient step -------------------------------------------

    auto x̂ₖ = alloc.alloc(), pₖ = alloc.alloc(), grad_̂ψₖ = alloc.alloc();
    auto xₖ₊₁ = alloc.alloc(), qₖ = alloc.alloc(), grad_ψₖ₊₁ = alloc.alloc();
    auto x̂ₖ₊₁ = alloc.alloc(), pₖ₊₁ = alloc.alloc();

    // Calculate x̂₀, p₀ (projected gradient step)
    detail::calc_x̂(γₖ, xₖ, C, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
    // Calculate ψ(x̂ₖ)
    real_t ψx̂ₖ        = ψ(x̂ₖ);
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t pₖᵀpₖ      = pₖ.squaredNorm();
    // Compute forward-backward envelope
    real_t φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Quadratic upper bound -----------------------------------------------
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            // Decrease step size until quadratic upper bound is satisfied
            real_t old_γₖ = detail::descent_lemma(
                ψ, C, params.quadratic_upperbound_tolerance_factor,
                params.L_max, xₖ, ψₖ, grad_ψₖ,
                /* in ⟹ out */ x̂ₖ, pₖ,
                /* inout */ ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ, γₖ);
            if (k > 0 && γₖ != old_γₖ) // Flush L-BFGS if γ changed
                direction_provider.changed_γ(γₖ, old_γₖ);
            else if (k == 0) // Initialize L-BFGS
                direction_provider.initialize(xₖ, x̂ₖ, pₖ, grad_ψₖ);
            if (γₖ != old_γₖ)
                φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
        }
        // Calculate ∇ψ(x̂ₖ)
        grad_ψ(x̂ₖ, grad_̂ψₖ);
        // Check stop condition ------------------------------------------------
        real_t εₖ = detail::calc_error_stop_crit(
            C, params.stop_crit, pₖ, γₖ, xₖ, x̂ₖ, vec(), grad_ψₖ, grad_̂ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            detail::print_progress(k, ψₖ, grad_ψₖ, pₖᵀpₖ, γₖ, εₖ);

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status  = detail::check_all_stop_conditions(
            params, time_elapsed, k, stop_signal, ε, εₖ, no_progress);
        if (stop_status != SolverStatus::Unknown) {
            x              = std::move(x̂ₖ);
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = stop_status;
            alloc.free(grad_ψₖ, xₖ, x̂ₖ, pₖ, grad_̂ψₖ, xₖ₊₁, qₖ, grad_ψₖ₊₁, x̂ₖ₊₁,
                       pₖ₊₁);
            return s;
        }

        // Calculate quasi-Newton step -----------------------------------------
        real_t step_size = -1;
        if (params.lbfgs_stepsize == LBFGSStepSize::BasedOnGradientStepSize)
            step_size = 1;
        if (k > 0)
            direction_provider.apply(xₖ, x̂ₖ, pₖ, step_size,
                                     /* in ⟹ out */ qₖ);

        // Line search initialization ------------------------------------------
        τ                  = 1;
        real_t σₖγₖ⁻¹pₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, ψx̂ₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, pₖ₊₁ᵀpₖ₊₁;
        real_t Lₖ₊₁, γₖ₊₁;
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
                ψₖ₊₁ = ψ(xₖ₊₁);
                grad_ψ(xₖ₊₁, grad_ψₖ₊₁);
            }

            // Calculate x̂ₖ₊₁, pₖ₊₁ (projected gradient step in xₖ₊₁)
            detail::calc_x̂(γₖ₊₁, xₖ₊₁, C, grad_ψₖ₊₁, /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁);
            // Calculate ψ(x̂ₖ₊₁)
            ψx̂ₖ₊₁ = ψ(x̂ₖ₊₁);

            // Quadratic upper bound -------------------------------------------
            grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
            pₖ₊₁ᵀpₖ₊₁      = pₖ₊₁.squaredNorm();
            real_t pₖ₊₁ᵀpₖ₊₁_ₖ = pₖ₊₁ᵀpₖ₊₁; // prox step with step size γₖ

            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                (void)detail::descent_lemma(
                    ψ, C, params.quadratic_upperbound_tolerance_factor,
                    params.L_max, xₖ₊₁, ψₖ₊₁, grad_ψₖ₊₁,
                    /* in ⟹ out */ x̂ₖ₊₁, pₖ₊₁,
                    /* inout */ ψx̂ₖ₊₁, pₖ₊₁ᵀpₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, Lₖ₊₁, γₖ₊₁);
            }

            // Compute forward-backward envelope
            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * pₖ₊₁ᵀpₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;
            // Compute line search condition
            ls_cond = φₖ₊₁ - (φₖ - σₖγₖ⁻¹pₖᵀpₖ);
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
        }

        // Update L-BFGS -------------------------------------------------------
        if (γₖ != γₖ₊₁) // Flush L-BFGS if γ changed
            direction_provider.changed_γ(γₖ₊₁, γₖ);

        s.lbfgs_rejected += not direction_provider.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁,
                                                          grad_ψₖ₊₁, C, γₖ₊₁);

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
        pₖ.swap(pₖ₊₁);
        grad_ψₖ.swap(grad_ψₖ₊₁);
        grad_ψₖᵀpₖ = grad_ψₖ₊₁ᵀpₖ₊₁;
        pₖᵀpₖ      = pₖ₊₁ᵀpₖ₊₁;
    }
    throw std::logic_error("[PANOC] loop error");
}

template <class DirectionProviderT = LBFGS, class ObjFunT, class ObjGradFunT>
PANOCStats panoc(ObjFunT &ψ, ObjGradFunT &grad_ψ, const Box &C, rvec x,
                 real_t ε, const PANOCParams &params,
                 PANOCDirection<DirectionProviderT> direction,
                 vec_allocator &alloc) {
    if (alloc.size() < panoc_min_alloc_size)
        throw std::logic_error("Allocator too small, should be at least " +
                               std::to_string(panoc_min_alloc_size));
    const auto n = x.size();
    if (alloc.vector_size() != n)
        throw std::logic_error("Allocator size mismatch");

    return panoc_impl(ψ, grad_ψ, C, x, ε, params, alloc, direction);
}

template <class DirectionProviderT = LBFGS, class ObjFunT, class ObjGradFunT>
PANOCStats panoc(ObjFunT &ψ, ObjGradFunT &grad_ψ, const Box &C, rvec x,
                 real_t ε, const PANOCParams &params,
                 PANOCDirection<DirectionProviderT> direction) {
    vec_allocator alloc{panoc_min_alloc_size, x.size()};
    return panoc_impl(ψ, grad_ψ, C, x, ε, params, alloc, direction);
}

} // namespace alpaqa