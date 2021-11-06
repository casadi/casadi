#pragma once

#include <alpaqa/inner/decl/panoc-stop-crit.hpp>
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

struct PGAParams {
    /// Parameters related to the Lipschitz constant estimate and step size.
    LipschitzEstimateParams Lipschitz;
    /// Maximum number of inner iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);
    /// Minimum Lipschitz constant estimate.
    real_t L_min = 1e-5;
    /// Maximum Lipschitz constant estimate.
    real_t L_max = 1e20;
    /// What stop criterion to use.
    PANOCStopCrit stop_crit = PANOCStopCrit::ApproxKKT;

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;

    real_t quadratic_upperbound_tolerance_factor =
        10 * std::numeric_limits<real_t>::epsilon();
};

struct PGAProgressInfo {
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
    const PGAParams &params;
};

/// Standard Proximal Gradient Algorithm without any bells and whistles.
///
/// @ingroup    grp_InnerSolvers
class PGASolver {
  public:
    using Params = PGAParams;

    PGASolver(const Params &params) : params(params) {}

    struct Stats {
        SolverStatus status = SolverStatus::Unknown;
        real_t ε            = inf;
        std::chrono::microseconds elapsed_time;
        unsigned iterations = 0;
    };

    using ProgressInfo = PGAProgressInfo;

    Stats operator()(const Problem &problem,        // in
                     crvec Σ,                       // in
                     real_t ε,                      // in
                     bool always_overwrite_results, // in
                     rvec x,                        // inout
                     rvec λ,                        // inout
                     rvec err_z);                   // out

    PGASolver &
    set_progress_callback(std::function<void(const ProgressInfo &)> cb) {
        this->progress_cb = cb;
        return *this;
    }

    std::string get_name() const { return "PGA"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
    std::function<void(const ProgressInfo &)> progress_cb;
};

using std::chrono::duration_cast;
using std::chrono::microseconds;

inline PGASolver::Stats
PGASolver::operator()(const Problem &problem,        // in
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
        pₖ(n),       // Projected gradient step
        ŷₖ(m),       // <?>
        grad_ψₖ(n),  // ∇ψ(xₖ)
        grad_ψx̂ₖ(n); // ∇ψ(x̂ₖ)

    vec work_n(n), work_m(m);

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within PGA (for readability in the main algorithm)
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
        std::cout << "[PGA]   " << std::setw(6) << k
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
            /* in ⟹ out */ ψₖ, grad_ψₖ, x̂ₖ, grad_ψx̂ₖ, work_n, work_m);
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

    unsigned no_progress = 0;

    // Main loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {
        // From previous iteration:
        //  - xₖ
        //  - grad_ψₖ
        //  - ψₖ

        // Quadratic upper bound -----------------------------------------------

        // Projected gradient step: x̂ₖ and pₖ
        calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
        // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
        real_t ψx̂ₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷₖ);
        // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
        real_t grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
        real_t pₖᵀpₖ      = pₖ.squaredNorm();

        // Decrease step size until quadratic upper bound is satisfied
        descent_lemma(xₖ, ψₖ, grad_ψₖ, x̂ₖ, pₖ, ŷₖ, ψx̂ₖ, pₖᵀpₖ, grad_ψₖᵀpₖ, Lₖ,
                      γₖ);

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

        auto time_elapsed    = std::chrono::steady_clock::now() - start_time;
        bool out_of_time     = time_elapsed > params.max_time;
        bool out_of_iter     = k == params.max_iter;
        bool interrupted     = stop_signal.stop_requested();
        bool not_finite      = not std::isfinite(εₖ);
        bool conv            = εₖ <= ε;
        bool max_no_progress = no_progress > 1;
        bool exit = conv || out_of_iter || out_of_time || not_finite ||
                    interrupted || max_no_progress;
        if (exit) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops.
            // TODO: move the computation of ẑ and g(x) to ALM?
            if (conv || interrupted || always_overwrite_results) {
                calc_err_z(x̂ₖ, /* in ⟹ out */ err_z);
                x = std::move(x̂ₖ);
                y = std::move(ŷₖ);
            }
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

        if (xₖ == x̂ₖ) // TODO: this is quite expensive to do on each iteration
            ++no_progress;
        else
            no_progress = 0;

        xₖ.swap(x̂ₖ);
        grad_ψₖ.swap(grad_ψx̂ₖ);
        ψₖ = ψx̂ₖ;
    }
    throw std::logic_error("[PGA]   loop error");
}

template <class InnerSolverStats>
struct InnerStatsAccumulator;

template <>
struct InnerStatsAccumulator<PGASolver::Stats> {
    std::chrono::microseconds elapsed_time;
    unsigned iterations = 0;
};

inline InnerStatsAccumulator<PGASolver::Stats> &
operator+=(InnerStatsAccumulator<PGASolver::Stats> &acc,
           const PGASolver::Stats &s) {
    acc.elapsed_time += s.elapsed_time;
    acc.iterations += s.iterations;
    return acc;
}

} // namespace alpaqa