#pragma once

#include <cstddef>
#include <panoc-alm/inner/detail/limited-memory-qr.hpp>
#include <panoc-alm/inner/detail/panoc-helpers.hpp>
#include <panoc-alm/util/atomic_stop_signal.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <Eigen/QR> // TODO: remove, only necessary for reference implementation

namespace pa {

struct GuardedAAPGAParams {
    struct {
        /// Relative step size for initial finite difference Lipschitz estimate.
        real_t ε = 1e-6;
        /// Minimum step size for initial finite difference Lipschitz estimate.
        real_t δ = 1e-12;
        /// Factor that relates step size γ and Lipschitz constant.
        real_t Lγ_factor = 0.95;
    } Lipschitz; ///< Parameters related to the Lipschitz constant estimate
                 ///  and step size.
    /// Length of the history to keep in the limited-memory QR algorithm.
    unsigned limitedqr_mem = 10;
    /// Maximum number of inner iterations.
    unsigned max_iter = 100;
    /// Maximum duration.
    std::chrono::microseconds max_time = std::chrono::minutes(5);

    /// When to print progress. If set to zero, nothing will be printed.
    /// If set to N != 0, progress is printed every N iterations.
    unsigned print_interval = 0;
};

/// Guarded Anderson Accelerated Proximal Gradient Algorithm
/// Vien V. Mai and Mikael Johansson, Anderson Acceleration of Proximal Gradient Methods.
/// https://arxiv.org/abs/1910.08590v2
///
/// @ingroup    grp_InnerSolvers
class GuardedAAPGA {
  public:
    using Params = GuardedAAPGAParams;

    GuardedAAPGA(const Params &params) : params(params) {}

    struct Stats {
        unsigned iterations = 0;
        real_t ε            = inf;
        std::chrono::microseconds elapsed_time;
        SolverStatus status = SolverStatus::Unknown;

        unsigned linesearch_failures = 0; // TODO: unused
        unsigned lbfgs_failures      = 0; // TODO: unused
        unsigned lbfgs_rejected      = 0; // TODO: unused
    };

    Stats operator()(const Problem &problem, // in
                     const vec &Σ,           // in
                     real_t ε,               // in
                     vec &x,                 // inout
                     vec &λ,                 // inout
                     vec &err_z);            // out

    std::string get_name() const { return "GuardedAAPGA"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
};

using std::chrono::duration_cast;
using std::chrono::microseconds;

inline GuardedAAPGA::Stats
GuardedAAPGA::operator()(const Problem &problem, // in
                         const vec &Σ,           // in
                         real_t ε,               // in
                         vec &x,                 // inout
                         vec &y,                 // inout
                         vec &err_z              // out
) {
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;
    const auto m = problem.m;

    vec xₖ = x,       // Value of x at the beginning of the iteration
        x̂ₖ(n),        // Value of x after a projected gradient step
        gₖ(n),        // <?>
        rₖ₋₁(n),      // <?>
        rₖ(n),        // <?>
        pₖ(n),        // Projected gradient step
        yₖ(n),        // Value of x after a gradient or AA step
        ŷₖ(m),        // <?>
        grad_ψₖ(n),   // ∇ψ(xₖ)
        grad_ψₖ₊₁(n), // ∇ψ(xₖ₊₁)
        grad_ψx̂ₖ(n);  // ∇ψ(x̂ₖ)

    vec work_n(n), work_m(m);

    unsigned m_AA = std::min(params.limitedqr_mem, n);
    LimitedMemoryQR qr(n, m_AA);
    mat G(n, m_AA);
    vec γ_LS(m_AA);

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within AAPGA (for readability in the main algorithm)
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
                              const vec &pₖ, real_t γₖ, real_t εₖ) {
        std::cout << "[AAPGA] " << std::setw(6) << k
                  << ": ψ = " << std::setw(13) << ψₖ
                  << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                  << ", ‖p‖ = " << std::setw(13) << pₖ.norm()
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

    // Calculate ∇ψ(x₀)
    calc_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

    // Estimate Lipschitz constant
    real_t Lₖ = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    if (Lₖ < std::numeric_limits<real_t>::epsilon()) {
        Lₖ = std::numeric_limits<real_t>::epsilon();
    } else if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }

    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;

    // First projected gradient step -------------------------------------------

    rₖ₋₁     = -γₖ * grad_ψₖ;
    yₖ       = xₖ + rₖ₋₁;
    xₖ       = project(yₖ, problem.C);
    G.col(0) = yₖ;

    unsigned no_progress = 0; // TODO

    // Calculate gradient in second iterate ------------------------------------

    // Calculate ψ(x₁) and ∇ψ(x₁)
    real_t ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

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
        real_t norm_sq_pₖ = pₖ.squaredNorm();
        real_t margin     = 0; // 1e-6 * std::abs(ψₖ₊₁); // TODO: Why?

        real_t old_γₖ = γₖ;
        // Decrease step size until quadratic upper bound is satisfied
        while (ψx̂ₖ > ψₖ + margin + grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ) {
            Lₖ *= 2;
            γₖ /= 2;

            // Projected gradient step: x̂ₖ and pₖ (with new step size)
            calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
            // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
            ψx̂ₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷₖ);
            // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
            grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
            norm_sq_pₖ = pₖ.squaredNorm();
        }

        // Flush Anderson buffers if γ changed
        if (γₖ != old_γₖ) {
            qr.reset();
            // TODO: show mathematically that this is a sensible thing to do:
            rₖ₋₁ *= γₖ / old_γₖ;
        }

        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(x̂ₖ, ŷₖ, /* in ⟹ out */ grad_ψx̂ₖ);

        // Check stop condition ------------------------------------------------

        real_t εₖ = detail::calc_error_stop_crit(pₖ, γₖ, grad_ψx̂ₖ, grad_ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖ, γₖ, εₖ);

        auto time_elapsed    = std::chrono::steady_clock::now() - start_time;
        bool out_of_time     = time_elapsed > params.max_time;
        bool out_of_iter     = k == params.max_iter;
        bool interrupted     = stop_signal.stop_requested();
        bool not_finite      = not std::isfinite(εₖ);
        bool conv            = εₖ <= ε;
        bool max_no_progress = no_progress > params.limitedqr_mem;
        bool exit = conv || out_of_iter || out_of_time || not_finite ||
                    interrupted || max_no_progress;
        if (exit) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops.
            // TODO: move the computation of ẑ and g(x) to ALM?
            calc_err_z(x̂ₖ, /* in ⟹ out */ err_z);
            x              = std::move(x̂ₖ);
            y              = std::move(ŷₖ);
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

        // Standard gradient descent step
        gₖ = xₖ - γₖ * grad_ψₖ;
        rₖ = gₖ - yₖ;

        // Update QR factorization for Anderson acceleration
        if (qr.num_columns() == m_AA)
            qr.remove_column();
        qr.add_column(rₖ - rₖ₋₁);

        // Solve least squares problem Anderson acceleration
        qr.solve(rₖ, γ_LS);

#if 0
            std::cout << "Q:\n" << qr.get_Q() << std::endl;
            std::cout << "R:\n" << qr.get_R() << std::endl;

        mat G_full(n, qr.num_columns() + 1);

        vec α_vec(qr.num_columns() + 1);
        α_vec(0) = γ_LS(0);
        α_vec.segment(1, qr.num_columns() - 1) =
            γ_LS.segment(1, qr.num_columns() - 1) -
            γ_LS.segment(0, qr.num_columns() - 1);
        α_vec(qr.num_columns()) = 1 - γ_LS(qr.num_columns() - 1);

#endif

        // Compute Anderson acceleration next iterate yₑₓₜ = ∑ₙ₌₀ αₙ gₙ
        size_t g_idx = qr.ring_head();
        // α₀ = γ₀      for n = 0
        real_t α      = γ_LS(0);
        yₖ            = α * G.col(g_idx);
#if 0
        G_full.col(0) = G.col(g_idx);
#endif
        for (size_t n = 1; n < qr.num_columns(); ++n) {
            g_idx = qr.ring_next(g_idx);
            // αₙ = γₙ - γₙ₋₁       for  0 < n < mₖ
            α = γ_LS(n) - γ_LS(n - 1);
            yₖ += α * G.col(g_idx);
#if 0
            G_full.col(n) = G.col(g_idx);
#endif
        }
        // αₘ = 1 - γₘ₋₁       for n = mₖ
        α = 1 - γ_LS(qr.num_columns() - 1);
        yₖ += α * gₖ;
#if 0
        G_full.col(qr.num_columns()) = gₖ;
#endif

#if 0
        std::cout << "G:\n" << G_full << std::endl;

        std::cout << "γ_LS:  " << γ_LS.topRows(qr.num_columns()).transpose()
                  << std::endl;
        std::cout << "α_LS:  " << α_vec.transpose() << std::endl;

        std::cout << "y_ext: " << yₖ.transpose() << std::endl;
        std::cout << "y_ext: " << (G_full * α_vec).transpose() << std::endl;
#endif

        // Add the new column to G
        G.col(qr.ring_tail()) = gₖ; // TODO: avoid copy

        // Project accelerated step onto feasible set
        xₖ = project(yₖ, problem.C);

        // Calculate the objective at the projected accelerated point
        real_t ψₖ₊₁ = calc_ψ_ŷ(xₖ, /* in ⟹ out */ ŷₖ);

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
            // std::cout << "\x1b[0;32m"
            //              "Accepted"
            //              "\x1b[0m"
            //           << std::endl;
        }
        // If not satisfied, take normal projected gradient step
        else {
            yₖ = gₖ; // TODO: eliminate copy somehow
            xₖ.swap(x̂ₖ);
            ψₖ = ψx̂ₖ;
            grad_ψₖ.swap(grad_ψx̂ₖ);
            // std::cout << "\x1b[0;31m"
            //              "Rejected"
            //              "\x1b[0m"
            //           << std::endl;
        }
        rₖ.swap(rₖ₋₁);
    }
    throw std::logic_error("[AAPGA] loop error");
}

class RefGuardedAAPGA {
  public:
    using Params = GuardedAAPGAParams;

    RefGuardedAAPGA(const Params &params) : params(params) {}

    struct Stats {
        unsigned iterations = 0;
        real_t ε            = inf;
        std::chrono::microseconds elapsed_time;
        SolverStatus status = SolverStatus::Unknown;

        unsigned linesearch_failures = 0; // TODO: unused
        unsigned lbfgs_failures      = 0; // TODO: unused
        unsigned lbfgs_rejected      = 0; // TODO: unused
    };

    Stats operator()(const Problem &problem, // in
                     const vec &Σ,           // in
                     real_t ε,               // in
                     vec &x,                 // inout
                     vec &λ,                 // inout
                     vec &err_z);            // out

    std::string get_name() const { return "GuardedAAPGA"; }

    void stop() { stop_signal.stop(); }

    const Params &get_params() const { return params; }

  private:
    Params params;
    AtomicStopSignal stop_signal;
};

inline void rotate_columns_left(mat &M) {
    Eigen::PermutationMatrix<Eigen::Dynamic> P(M.cols());
    P.setIdentity();
    std::rotate(P.indices().data(), P.indices().data() + 1,
                P.indices().data() + P.size());
    M.applyOnTheRight(P);
}

inline RefGuardedAAPGA::Stats
RefGuardedAAPGA::operator()(const Problem &problem, // in
                            const vec &Σ,           // in
                            real_t ε,               // in
                            vec &x,                 // inout
                            vec &y,                 // inout
                            vec &err_z              // out
) {
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.n;
    const auto m = problem.m;

    vec xₖ = x,       // Value of x at the beginning of the iteration
        x̂ₖ(n),        // Value of x after a projected gradient step
        gₖ(n),        // <?>
        rₖ₋₁(n),      // <?>
        rₖ(n),        // <?>
        pₖ(n),        // Projected gradient step
        yₖ(n),        // Value of x after a gradient or AA step
        ŷₖ(m),        // <?>
        grad_ψₖ(n),   // ∇ψ(xₖ)
        grad_ψₖ₊₁(n), // ∇ψ(xₖ₊₁)
        grad_ψx̂ₖ(n);  // ∇ψ(x̂ₖ)

    vec work_n(n), work_m(m);

    unsigned m_AA = std::min(params.limitedqr_mem, n);
    LimitedMemoryQR qr(n, m_AA);
    mat G(n, m_AA);
    mat R(n, m_AA);
    vec γ_LS(m_AA);
    size_t stored_cols   = 0;
    size_t stored_cols_G = 1;

    // Helper functions --------------------------------------------------------

    // Wrappers for helper functions that automatically pass along any arguments
    // that are constant within AA-PGA (for readability in the main algorithm)
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
                              const vec &pₖ, real_t γₖ, real_t εₖ) {
        std::cout << "[AAPGA] " << std::setw(6) << k
                  << ": ψ = " << std::setw(13) << ψₖ
                  << ", ‖∇ψ‖ = " << std::setw(13) << grad_ψₖ.norm()
                  << ", ‖p‖ = " << std::setw(13) << pₖ.norm()
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

    // Calculate ∇ψ(x₀)
    calc_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

    // Estimate Lipschitz constant
    real_t Lₖ = (grad_ψₖ₊₁ - grad_ψₖ).norm() / h.norm();
    if (Lₖ < std::numeric_limits<real_t>::epsilon()) {
        Lₖ = std::numeric_limits<real_t>::epsilon();
    } else if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }

    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;

    // First projected gradient step -------------------------------------------

    rₖ₋₁     = -γₖ * grad_ψₖ;
    yₖ       = xₖ + rₖ₋₁;
    xₖ       = project(yₖ, problem.C);
    G.col(0) = yₖ;

    unsigned no_progress = 0; // TODO

    // Calculate gradient in second iterate ------------------------------------

    // Calculate ψ(x₁) and ∇ψ(x₁)
    real_t ψₖ = calc_ψ_grad_ψ(xₖ, /* in ⟹ out */ grad_ψₖ);

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
        real_t norm_sq_pₖ = pₖ.squaredNorm();
        real_t margin     = 0; // 1e-6 * std::abs(ψₖ₊₁); // TODO: Why?

        real_t old_γₖ = γₖ;
        // Decrease step size until quadratic upper bound is satisfied
        while (ψx̂ₖ > ψₖ + margin + grad_ψₖᵀpₖ + 0.5 * Lₖ * norm_sq_pₖ) {
            Lₖ *= 2;
            γₖ /= 2;

            // Projected gradient step: x̂ₖ and pₖ (with new step size)
            calc_x̂(γₖ, xₖ, grad_ψₖ, /* in ⟹ out */ x̂ₖ, pₖ);
            // Calculate ψ(x̂ₖ) and ŷ(x̂ₖ)
            ψx̂ₖ = calc_ψ_ŷ(x̂ₖ, /* in ⟹ out */ ŷₖ);
            // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
            grad_ψₖᵀpₖ = grad_ψₖ.dot(pₖ);
            norm_sq_pₖ = pₖ.squaredNorm();
        }

        // Flush Anderson buffers if γ changed
        if (γₖ != old_γₖ) {
            qr.reset();
            // TODO: show mathematically that this is a sensible thing to do:
            rₖ₋₁ *= γₖ / old_γₖ;
        }

        // Calculate ∇ψ(x̂ₖ)
        calc_grad_ψ_from_ŷ(x̂ₖ, ŷₖ, /* in ⟹ out */ grad_ψx̂ₖ);

        // Check stop condition ------------------------------------------------

        real_t εₖ = detail::calc_error_stop_crit(pₖ, γₖ, grad_ψx̂ₖ, grad_ψₖ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, ψₖ, grad_ψₖ, pₖ, γₖ, εₖ);

        auto time_elapsed    = std::chrono::steady_clock::now() - start_time;
        bool out_of_time     = time_elapsed > params.max_time;
        bool out_of_iter     = k == params.max_iter;
        bool interrupted     = stop_signal.stop_requested();
        bool not_finite      = not std::isfinite(εₖ);
        bool conv            = εₖ <= ε;
        bool max_no_progress = no_progress > params.limitedqr_mem;
        bool exit = conv || out_of_iter || out_of_time || not_finite ||
                    interrupted || max_no_progress;
        if (exit) {
            // TODO: We could cache g(x) and ẑ, but would that faster?
            //       It saves 1 evaluation of g per ALM iteration, but requires
            //       many extra stores in the inner loops.
            // TODO: move the computation of ẑ and g(x) to ALM?
            calc_err_z(x̂ₖ, /* in ⟹ out */ err_z);
            x              = std::move(x̂ₖ);
            y              = std::move(ŷₖ);
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

        // Standard gradient descent step
        gₖ = xₖ - γₖ * grad_ψₖ;
        rₖ = gₖ - yₖ;

        // Update QR factorization for Anderson acceleration
        if (stored_cols == m_AA) {
            rotate_columns_left(R);
            --stored_cols;
        }
        R.col(stored_cols) = rₖ - rₖ₋₁;
        ++stored_cols;

        // Solve least squares problem Anderson acceleration
        auto qr = R.block(0, 0, n, stored_cols).householderQr();
        γ_LS    = qr.solve(rₖ);

        std::cout << "Q:\n"
                  << qr.householderQ() * mat::Identity(n, stored_cols)
                  << std::endl;
        std::cout << "R:\n"
                  << qr.matrixQR()
                         .block(0, 0, stored_cols, stored_cols)
                         .triangularView<Eigen::Upper>()
                         .toDenseMatrix()
                  << std::endl;

        // Compute Anderson acceleration next iterate yₑₓₜ = ∑ₙ₌₀ αₙ gₙ
        size_t g_idx = 0;
        // α₀ = γ₀      for n = 0
        real_t α = γ_LS(0);
        yₖ       = α * G.col(g_idx);
        for (size_t n = 1; n < stored_cols; ++n) {
            ++g_idx;
            // αₙ = γₙ - γₙ₋₁       for  0 < n < mₖ
            α = γ_LS(n) - γ_LS(n - 1);
            yₖ += α * G.col(g_idx);
        }
        // αₘ = 1 - γₘ₋₁       for n = mₖ
        α = 1 - γ_LS(stored_cols - 1);
        yₖ += α * gₖ;

        mat G_full(n, stored_cols_G + 1);
        G_full.leftCols(stored_cols) = G.block(0, 0, n, stored_cols_G);
        G_full.rightCols(1)          = gₖ;
        std::cout << "G:\n" << G_full << std::endl;

        vec α_vec(stored_cols + 1);
        α_vec(0) = γ_LS(0);
        α_vec.segment(1, stored_cols - 1) =
            γ_LS.segment(1, stored_cols - 1) - γ_LS.segment(0, stored_cols - 1);
        α_vec(stored_cols) = 1 - γ_LS(stored_cols - 1);
        std::cout << "γ_LS:  " << γ_LS.topRows(stored_cols).transpose()
                  << std::endl;
        std::cout << "α_LS:  " << α_vec.transpose() << std::endl;

        std::cout << "y_ext: " << yₖ.transpose() << std::endl;
        std::cout << "y_ext: " << (G_full * α_vec).transpose() << std::endl;

        // Add the new column to G
        if (stored_cols_G == m_AA) {
            rotate_columns_left(G);
            --stored_cols_G;
        }
        G.col(stored_cols_G) = gₖ; // TODO: avoid copy
        ++stored_cols_G;

        // Project accelerated step onto feasible set
        xₖ = project(yₖ, problem.C);

        // std::cout << "yₖ:    " << yₖ.transpose() << std::endl;

        // Calculate the objective at the projected accelerated point
        real_t ψₖ₊₁ = calc_ψ_ŷ(xₖ, /* in ⟹ out */ ŷₖ);

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
            std::cout << "\x1b[0;32m"
                         "Accepted"
                         "\x1b[0m"
                      << std::endl;
        }
        // If not satisfied, take normal projected gradient step
        else {
            yₖ = gₖ; // TODO: eliminate copy somehow
            xₖ.swap(x̂ₖ);
            ψₖ = ψx̂ₖ;
            grad_ψₖ.swap(grad_ψx̂ₖ);
            std::cout << "\x1b[0;31m"
                         "Rejected"
                         "\x1b[0m"
                      << std::endl;
        }
        rₖ.swap(rₖ₋₁);
    }
    throw std::logic_error("[AAPGA] loop error");
}

} // namespace pa