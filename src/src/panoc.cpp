#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "private/panoc-helpers.hpp"
#include <panoc-alm/lbfgs.hpp>
#include <panoc-alm/panoc.hpp>
#include <panoc-alm/solverstatus.hpp>

namespace pa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

PANOCSolver::Stats PANOCSolver::operator()(
    const Problem &problem, ///< [in]    Problem description
    vec &x,                 ///< [inout] Decision variable @f$ x @f$
    vec &z,                 ///< [out]   Slack variable @f$ x @f$
    vec &y,                 ///< [inout] Lagrange multiplier @f$ x @f$
    vec &err_z,             ///< [out]   Slack variable error @f$ g(x) - z @f$
    const vec &Σ,           ///< [in]    Constraint weights @f$ \Sigma @f$
    real_t ε                ///< [in]    Tolerance @f$ \epsilon @f$
) {
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    using namespace detail;

    const auto n = x.size();
    const auto m = z.size();

    // TODO: allocates
    LBFGS lbfgs;
    SpecializedLBFGS slbfgs;
    params.specialized_lbfgs ? slbfgs.resize(n, params.lbfgs_mem)
                             : lbfgs.resize(n, params.lbfgs_mem);

    vec xₖ = x,       // Value of x at the beginning of the iteration
        x̂ₖ(n),        // Value of x after a projected gradient step
        xₖ₊₁(n),      // xₖ for next iteration
        x̂ₖ₊₁(n),      // x̂ₖ for next iteration
        ŷx̂ₖ(m),       // Σ (g(x̂ₖ) - ẑₖ)
        ŷx̂ₖ₊₁(m),     // ŷ(x̂ₖ) for next iteration
        pₖ(n),        // xₖ - x̂ₖ
        pₖ₊₁(n),      // xₖ₊₁ - x̂ₖ₊₁
        qₖ(n),        // Newton step Hₖ pₖ
        grad_ψₖ(n),   // ∇ψ(xₖ)
        grad_̂ψₖ(n),   // ∇ψ(x̂ₖ)
        grad_ψₖ₊₁(n); // ∇ψ(xₖ₊₁)

    vec work_n(n), work_m(m);

    // Estimate Lipschitz constant using finite difference
    vec h(n);
    h = (x * params.Lipschitz.ε).cwiseAbs().cwiseMax(params.Lipschitz.δ);
    x += h;

    // Calculate ∇ψ(x₀ + h)
    calc_grad_ψ(problem, x, y, Σ, // in
                grad_ψₖ₊₁,        // out
                work_n, work_m);  // work
    // Calculate ψ(xₖ), ∇ψ(x₀)
    real_t ψₖ = calc_ψ_grad_ψ(problem, xₖ, y, Σ, // in
                              grad_ψₖ,           // out
                              work_n, work_m);   // work

    // Estimate Lipschitz constant
    real_t Lₖ = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    if (Lₖ < std::numeric_limits<real_t>::epsilon()) {
        Lₖ = std::numeric_limits<real_t>::epsilon();
    } else if (!std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }

    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t σₖ = γₖ * (1 - γₖ * Lₖ) / 2;

    // Calculate x̂₀, p₀ (projected gradient step)
    // don't check progress here
    calc_x̂(problem, γₖ, xₖ, grad_ψₖ, // in
           x̂ₖ, pₖ);                  // out
    // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
    real_t ψ̂xₖ = calc_ψ_ŷ(problem, x̂ₖ, y, Σ, // in
                          ŷx̂ₖ);              // out

    real_t margin = 0; // 1e-6 * std::abs(ψₖ); // TODO: OpEn did this. Why?
    real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
    real_t norm_sq_pₖ = pₖ.squaredNorm();

    // Compute forward-backward envelope
    real_t φₖ = ψₖ + 1 / (2 * γₖ) * norm_sq_pₖ + grad_ψₖᵀpₖ;

    // Main loop
    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // Decrease step size until quadratic upper bound is satisfied
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            while (ψ̂xₖ > ψₖ + margin + grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ) {
                Lₖ *= 2;
                σₖ /= 2;
                γₖ /= 2;

                // Flush L-BFGS if γ changed
                if (k > 0 && params.specialized_lbfgs == false)
                    lbfgs.reset();

                // Calculate x̂ₖ and pₖ (with new step size)
                calc_x̂(problem, γₖ, xₖ, grad_ψₖ, // in
                       x̂ₖ, pₖ);                  // out
                // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
                grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
                norm_sq_pₖ = pₖ.squaredNorm();

                // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
                ψ̂xₖ = calc_ψ_ŷ(problem, x̂ₖ, y, Σ, // in
                               ŷx̂ₖ);              // out
            }
        }

        // Initialize the L-BFGS
        if (params.specialized_lbfgs == true && k == 0)
            slbfgs.initialize(xₖ, grad_ψₖ, x̂ₖ, γₖ);

        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(problem, x̂ₖ, ŷx̂ₖ, // in
                           grad_̂ψₖ,          // out
                           work_n);          // work

        // Check stop condition
        real_t εₖ = calc_error_stop_crit(pₖ, γₖ, grad_̂ψₖ, grad_ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0) {
            std::cout << "[PANOC] " << std::setw(6) << k
                      << ": ψ = " << std::setw(13) << ψₖ
                      << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                      << ", ‖p‖ = " << std::setw(13) << std::sqrt(norm_sq_pₖ)
                      << ", γ = " << std::setw(13) << γₖ
                      << ", εₖ = " << std::setw(13) << εₖ << "\r\n";
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        bool out_of_time  = time_elapsed > params.max_time;
        if (εₖ <= ε || k == params.max_iter || out_of_time) {
            if (params.print_interval > 0) {
                const Eigen::IOFormat fmt(16, 0, "\t", " ", "", "", "", "");
                std::cout << "∇ψₖ:       " << grad_ψₖ.transpose().format(fmt)
                          << std::endl;
                std::cout << "∇̂ψₖ:       " << grad_̂ψₖ.transpose().format(fmt)
                          << std::endl;
                std::cout << "∇̂ψₖ - ∇ψₖ: "
                          << (grad_̂ψₖ - grad_ψₖ).transpose().format(fmt)
                          << std::endl;
                std::cout << "p/γ:       " << (pₖ / γₖ).transpose().format(fmt)
                          << std::endl;
                std::cout << "p:         " << pₖ.transpose().format(fmt)
                          << std::endl;
                std::cout << "γ·∇ψₖ:     "
                          << (γₖ * grad_ψₖ.transpose()).format(fmt)
                          << std::endl;
                std::cout << "xl:        "
                          << problem.C.lowerbound.transpose().format(fmt)
                          << std::endl;
                std::cout << "x:         " << xₖ.transpose().format(fmt)
                          << std::endl;
                std::cout << "xu:        "
                          << problem.C.upperbound.transpose().format(fmt)
                          << std::endl;
                std::cout << "x̂:         " << x̂ₖ.transpose().format(fmt)
                          << std::endl;
                std::cout << "γ:         " << γₖ << std::endl;
            }

            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops of PANOC.
            // TODO: move the computation of ẑ and g(x) to ALM?
            calc_ẑ(problem, x̂ₖ, y, Σ, z, err_z);
            x = std::move(x̂ₖ);
            y = std::move(ŷx̂ₖ);

            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = εₖ <= ε       ? SolverStatus::Converged
                             : out_of_time ? SolverStatus::MaxTime
                                           : SolverStatus::MaxIter;
            return s;
        } else if (!std::isfinite(εₖ)) {
            std::cerr << "[PANOC] "
                         "\x1b[0;31m"
                         "inf/NaN"
                         "\x1b[0m"
                      << std::endl;
            std::cout << "[k]   " << k << std::endl;
            std::cout << "qₖ₋₁: " << qₖ.transpose() << std::endl;
            std::cout << "xₖ:   " << xₖ.transpose() << std::endl;
            std::cout << "x̂ₖ:   " << x̂ₖ.transpose() << std::endl;
            std::cout << "ŷx̂ₖ:  " << ŷx̂ₖ.transpose() << std::endl;
            std::cout << "pₖ:   " << pₖ.transpose() << std::endl;
            std::cout << "γₖ:   " << γₖ << std::endl;
            std::cout << "∇_̂ψₖ:  " << grad_̂ψₖ.transpose() << std::endl;
            std::cout << "∇ψₖ:  " << grad_ψₖ.transpose() << std::endl;

            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = SolverStatus::NotFinite;
            return s;
        } else if (stop_signal.load(std::memory_order_relaxed)) {
            calc_ẑ(problem, x̂ₖ, y, Σ, z, err_z);
            x = std::move(x̂ₖ);
            y = std::move(ŷx̂ₖ);

            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = SolverStatus::Interrupted;
            return s;
        }

        // Calculate quasi-Newton step
        if (k > 0) {
            qₖ = pₖ;
            params.specialized_lbfgs ? slbfgs.apply(qₖ) : lbfgs.apply(qₖ);
        }

        // Line search
        real_t σ_norm_γ⁻¹pₖ = σₖ * norm_sq_pₖ / (γₖ * γₖ);
        real_t φₖ₊₁, ψₖ₊₁, ψ̂xₖ₊₁, grad_ψₖ₊₁ᵀpₖ₊₁, norm_sq_pₖ₊₁;
        real_t τ = 1;
        real_t Lₖ₊₁, σₖ₊₁, γₖ₊₁;
        real_t ls_cond;

        // Make sure quasi-Newton step is valid
        if (k == 0) {
            τ = 0;
        } else if (qₖ.hasNaN()) {
            τ = 0;
            ++s.lbfgs_failures;
            params.specialized_lbfgs ? slbfgs.reset() : lbfgs.reset();
        }

        // Line search loop
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
            ψₖ₊₁ = calc_ψ_grad_ψ(problem, xₖ₊₁, y, Σ, // in
                                 grad_ψₖ₊₁,           // out
                                 work_n, work_m);     // work
            // Calculate x̂ₖ₊₁, pₖ₊₁ (projected gradient step)
            calc_x̂(problem, γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, // in
                   x̂ₖ₊₁, pₖ₊₁);                    // out
            // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
            ψ̂xₖ₊₁ = calc_ψ_ŷ(problem, x̂ₖ₊₁, y, Σ, // in
                             ŷx̂ₖ₊₁);              // out

            // TODO: OpEn did this. Why?
            real_t margin  = 0; // 1e-6 * std::abs(ψₖ);
            grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
            norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                while (ψ̂xₖ₊₁ > ψₖ₊₁ + margin + grad_ψₖ₊₁ᵀpₖ₊₁ +
                                   0.5 * Lₖ₊₁ * norm_sq_pₖ₊₁) {
                    Lₖ₊₁ *= 2;
                    σₖ₊₁ /= 2;
                    γₖ₊₁ /= 2;
                    // Flush L-BFGS if γ changed
                    if (params.specialized_lbfgs == false)
                        lbfgs.reset();

                    // Calculate x̂ₖ₊₁ and pₖ₊₁ (with new step size)
                    calc_x̂(problem, γₖ₊₁, xₖ₊₁, grad_ψₖ₊₁, // in
                           x̂ₖ₊₁, pₖ₊₁);                    // out
                    // Calculate ∇ψ(xₖ₊₁)ᵀpₖ₊₁ and ‖pₖ₊₁‖²
                    grad_ψₖ₊₁ᵀpₖ₊₁ = grad_ψₖ₊₁.dot(pₖ₊₁);
                    norm_sq_pₖ₊₁   = pₖ₊₁.squaredNorm();
                    // Calculate ψ(x̂ₖ₊₁) and ŷ(x̂ₖ₊₁)
                    ψ̂xₖ₊₁ = calc_ψ_ŷ(problem, x̂ₖ₊₁, y, Σ, // in
                                     ŷx̂ₖ₊₁);              // out
                }
            }

            // Compute forward-backward envelope
            φₖ₊₁ = ψₖ₊₁ + 1 / (2 * γₖ₊₁) * norm_sq_pₖ₊₁ + grad_ψₖ₊₁ᵀpₖ₊₁;

            τ /= 2;

            ls_cond = φₖ₊₁ - (φₖ - σ_norm_γ⁻¹pₖ);
        } while (ls_cond > 0 && τ >= params.τ_min);

        // τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min && k != 0) {
            ++s.linesearch_failures;
        }

        // Update L-BFGS
        s.lbfgs_rejected +=
            params.specialized_lbfgs
                ? !slbfgs.update(xₖ₊₁, grad_ψₖ₊₁, x̂ₖ₊₁, problem.C, γₖ₊₁)
                : !lbfgs.update(xₖ₊₁ - xₖ, pₖ - pₖ₊₁);

        // Advance step
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
