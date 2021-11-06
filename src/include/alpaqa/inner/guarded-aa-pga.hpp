#pragma once

#include <alpaqa/inner/decl/panoc-stop-crit.hpp>
#include <alpaqa/inner/detail/anderson-helpers.hpp>
#include <alpaqa/inner/detail/panoc-helpers.hpp>
#include <alpaqa/util/atomic_stop_signal.hpp>
#include <alpaqa/util/lipschitz.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace alpaqa {

struct GAAPGAParams {
    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams Lipschitz;
    /// Length of the history to keep in the limited-memory QR algorithm.
    unsigned limitedqr_mem = 10;
    /// Maximum number of inner iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);
    /// Minimum Lipschitz constant estimate.
    real_t L_min = 1e-5;
    /// Maximum Lipschitz constant estimate.
    real_t L_max = 1e20;
    /// What stopping criterion to use.
    PANOCStopCrit stop_crit = PANOCStopCrit::ApproxKKT;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;

    real_t quadratic_upperbound_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();

    /// Maximum number of iterations without any progress before giving up.
    unsigned max_no_progress = 10;

    bool full_flush_on_γ_change = true;
};

struct GAAPGAProgressInfo {
    unsigned k;
    crvec x;
    crvec p;
    real_t norm_sq_p;
    crvec x_hat;
    real_t ψ;
    crvec grad_ψ;
    real_t ψ_hat;
    crvec grad_ψ_hat;
    real_t L;
    real_t γ;
    real_t ε;
    crvec Σ;
    crvec y;
    const Problem &problem;
    const GAAPGAParams &params;
};

/// Guarded Anderson Accelerated Proximal Gradient Algorithm.
/// Vien V. Mai and Mikael Johansson, Anderson Acceleration of Proximal Gradient Methods.
/// https://arxiv.org/abs/1910.08590v2
///
/// @ingroup    grp_InnerSolvers
class GAAPGASolver {
  public:
    using Params = GAAPGAParams;

    GAAPGASolver(const Params &params) : params(params) {}

    struct Stats {
        SolverStatus status = SolverStatus::Unknown;
        real_t ε            = inf;
        std::chrono::microseconds elapsed_time;
        unsigned iterations                 = 0;
        unsigned accelerated_steps_accepted = 0;
    };

    using ProgressInfo = GAAPGAProgressInfo;

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec λ,                        // inout
                     rvec err_z);                   // out

    GAAPGASolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const { return "GAAPGA"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
};

using std::chrono::duration_cast;
using std::chrono::microseconds;

inline GAAPGASolver::Stats
GAAPGASolver::operator()(const Problem &problem,        // in
                         crvec Σ,                       // in
                         real_t ε,                      // in
                         bool always_overwrite_results, // in
                         rvec x,                        // inout
                         rvec y,                        // inout
                         rvec err_z                     // out
) {
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;
    const auto m = problem.m;

    vec xₖ = x,      // Value of x at the beginning of the iteration
        x̂ₖ(n),       // Value of x after a projected gradient step
        gₖ(n),       // <?>
        rₖ₋₁(n),     // <?>
        rₖ(n),       // <?>
        pₖ(n),       // Projected gradient step
        yₖ(n),       // Value of x after a gradient or AA step
        ŷₖ(m),       // <?>
        grad_ψₖ(n),  // ∇ψ(xₖ)
        grad_ψx̂ₖ(n); // ∇ψ(x̂ₖ)

    vec work_n(n), work_n2(n), work_m(m);

    unsigned m_AA = std::min(params.limitedqr_mem, n);
    LimitedMemoryQR qr(n, m_AA);
    mat G(n, m_AA);
    vec γ_LS(m_AA);

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within AAPGA (for readability in the main algorithm)
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
    auto print_progress = [&](unsigned k, real_t ψₖ, crvec grad_ψₖ, crvec pₖ,
                              real_t γₖ, real_t εₖ) {
        std::cout << "[AAPGA] " << std::setw(6) << k
                  << ": ψ = " << std::setw(13) << ψₖ
                  << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                  << ", ‖p‖ = " << std::setw(13) << pₖ.norm()
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
            /* in ⟹ out */ grad_ψₖ, /* work */ x̂ₖ, work_n2, work_n, work_m);
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

    rₖ₋₁     = -γₖ * grad_ψₖ;
    yₖ       = xₖ + rₖ₋₁;
    xₖ       = project(yₖ, problem.C);
    G.col(0) = yₖ;

    unsigned no_progress = 0;

    // Calculate gradient in second iterate ------------------------------------

    // Calculate ψ(x₁) and ∇ψ(x₁)
    ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

    // Main loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // From previous iteration:
        //  - xₖ
        //  - grad_ψₖ
        //  - ψₖ
        //  - rₖ₋₁
        //  - history in qr and G

        // Quadratic upper bound -----------------------------------------------

        // Projected gradient step: x̂ₖ and pₖ
        calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
        // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
        real_t ψx̂ₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷₖ);
        // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
        real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
        real_t pₖᵀpₖ      = pₖ.squaredNorm();

        real_t old_γₖ = descent_lemma(xₖ, ψₖ, grad_ψₖ, x̂ₖ, pₖ, ŷₖ, ψx̂ₖ, pₖᵀpₖ,
                                      grad_ψₖᵀpₖ, Lₖ, γₖ);

        // Flush or update Anderson buffers if step size changed
        if (γₖ != old_γₖ) {
            if (params.full_flush_on_γ_change) {
                // Save the latest function evaluation gₖ at the first index
                size_t newest_g_idx = qr.ring_tail();
                if (newest_g_idx != 0)
                    G.col(0) = G.col(newest_g_idx);
                // Flush everything else and reset indices
                qr.reset();
            } else {
                // When not near the boundaries of the feasible set,
                // r(x) = g(x) - x = Π(x - γ∇ψ(x)) - x = -γ∇ψ(x),
                // in other words, r(x) is proportional to γ, and so is Δr,
                // so when γ changes, these values have to be updated as well
                qr.scale_R(γₖ / old_γₖ);
            }
            rₖ₋₁ *= γₖ / old_γₖ;
        }

        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(x̂ₖ, ŷₖ, /* in ⟹ out */ grad_ψx̂ₖ);

        // Check stop condition ------------------------------------------------

        real_t εₖ = detail::calc_error_stop_crit(
            problem.C, params.stop_crit, pₖ, γₖ, xₖ, x̂ₖ, ŷₖ, grad_ψₖ, grad_ψx̂ₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖ, γₖ, εₖ);
        if (progress_cb)
            progress_cb({k, xₖ, pₖ, pₖᵀpₖ, x̂ₖ, ψₖ, grad_ψₖ, ψx̂ₖ, grad_ψx̂ₖ, Lₖ,
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
                y = std::move(ŷₖ);
            }
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<microseconds>(time_elapsed);
            s.status       = stop_status;
            return s;
        }

        // Standard gradient descent step
        gₖ = xₖ - γₖ * grad_ψₖ;
        rₖ = gₖ - yₖ;

        // Solve Anderson acceleration least squares problem and update history
        minimize_update_anderson(qr, G, rₖ, rₖ₋₁, gₖ, γ_LS, yₖ);

        // Project accelerated step onto feasible set
        xₖ = project(yₖ, problem.C);

        // Calculate the objective at the projected accelerated point
        real_t ψₖ₊₁  = calc_ψ_ŷ(xₖ, /* in ⟹ out */ ŷₖ);
        real_t old_ψ = ψₖ;

        // Check sufficient decrease condition for Anderson iterate
        bool sufficient_decrease;
        if (0) // From paper
            sufficient_decrease = ψₖ₊₁ <= ψₖ - 0.5 * γₖ * grad_ψₖ.squaredNorm();
        else // Since we compute ψ(x̂ₖ) we might as well pick the best one
            sufficient_decrease = ψₖ₊₁ <= ψx̂ₖ;

        if (sufficient_decrease) {
            // Accept Anderson step
            // yₖ and xₖ are already overwritten by yₑₓₜ and Π(yₑₓₜ)
            ψₖ = ψₖ₊₁;
            calc_grad_ψ_from_ŷ(xₖ, ŷₖ, /* in ⟹ out */ grad_ψₖ);
        }
        // If not satisfied, take normal projected gradient step
        else {
            yₖ.swap(gₖ);
            xₖ.swap(x̂ₖ);
            ψₖ = ψx̂ₖ;
            grad_ψₖ.swap(grad_ψx̂ₖ);
        }
        rₖ.swap(rₖ₋₁);
        s.accelerated_steps_accepted += sufficient_decrease;

        // Check if we made any progress, prevents us from exceeding the maximum
        // number of iterations doing nothing if the step size gets too small
        // TODO: is this a valid test?
        no_progress = (ψₖ == old_ψ) ? no_progress + 1 : 0;
    }
    throw std::logic_error("[AAPGA] loop error");
}

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<GAAPGASolver::Stats> {
    std::chrono::microseconds elapsed_time;
    unsigned iterations                 = 0;
    unsigned accelerated_steps_accepted = 0;
};

inline InnerStatsAccumulator<GAAPGASolver::Stats> &
operator+=(InnerStatsAccumulator<GAAPGASolver::Stats> &acc,
           const GAAPGASolver::Stats &s) {
    acc.elapsed_time += s.elapsed_time;
    acc.iterations += s.iterations;
    acc.accelerated_steps_accepted += s.accelerated_steps_accepted;
    return acc;
}

} // namespace alpaqa