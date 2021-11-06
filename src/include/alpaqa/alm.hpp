#pragma once

#include <alpaqa/detail/alm-helpers.hpp>
#include <alpaqa/util/solverstatus.hpp>

#include <iomanip>
#include <iostream>

namespace alpaqa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

template <class InnerSolverT>
typename ALMSolver<InnerSolverT>::Stats
ALMSolver<InnerSolverT>::operator()(const Problem &problem, rvec y, rvec x) {
    auto start_time = std::chrono::steady_clock::now();

    constexpr auto sigNaN = std::numeric_limits<real_t>::signaling_NaN();
    vec Σ                 = vec::Constant(problem.m, sigNaN);
    vec Σ_old             = vec::Constant(problem.m, sigNaN);
    vec error₁            = vec::Constant(problem.m, sigNaN);
    vec error₂            = vec::Constant(problem.m, sigNaN);
    real_t norm_e₁        = sigNaN;
    real_t norm_e₂        = sigNaN;

    Stats s;

    Problem prec_problem;
    real_t prec_f;
    vec prec_g;

    if (params.preconditioning)
        detail::apply_preconditioning(problem, prec_problem, x, prec_f, prec_g);
    const auto &p = params.preconditioning ? prec_problem : problem;

    // Initialize the penalty weights
    if (params.Σ₀ > 0) {
        Σ.fill(params.Σ₀);
    }
    // Initial penalty weights from problem
    else {
        detail::initialize_penalty(p, params, x, Σ);
    }

    real_t ε                   = params.ε₀;
    real_t ε_old               = sigNaN;
    real_t Δ                   = params.Δ;
    real_t ρ                   = params.ρ;
    bool first_successful_iter = true;

    for (unsigned int i = 0; i < params.max_iter; ++i) {
        // TODO: this is unnecessary when the previous iteration lowered the
        // penalty update factor.
        detail::project_y(y, p.D.lowerbound, p.D.upperbound, params.M);
        // Check if we're allowed to lower the penalty factor even further.
        bool out_of_penalty_factor_updates =
            (first_successful_iter
                 ? s.initial_penalty_reduced == params.max_num_initial_retries
                 : s.penalty_reduced == params.max_num_retries) ||
            (s.initial_penalty_reduced + s.penalty_reduced ==
             params.max_total_num_retries);
        bool out_of_iter = i + 1 == params.max_iter;
        // If this is the final iteration, or the final chance to reduce the
        // penalty update factor, the inner solver can just return its results,
        // even if it doesn't converge.
        bool overwrite_results = out_of_iter || out_of_penalty_factor_updates;

        // Inner solver
        // ------------

        // Call the inner solver to minimize the augmented lagrangian for fixed
        // Lagrange multipliers y.
        auto ps = inner_solver(p, Σ, ε, overwrite_results, x, y, error₂);
        bool inner_converged = ps.status == SolverStatus::Converged;
        // Accumulate the inner solver statistics
        s.inner_convergence_failures += not inner_converged;
        s.inner += ps;

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        bool out_of_time  = time_elapsed > params.max_time;
        bool backtrack =
            not inner_converged && not overwrite_results && not out_of_time;

        // Print statistics of current iteration
        if (params.print_interval != 0 && i % params.print_interval == 0) {
            real_t δ       = backtrack ? NaN : vec_util::norm_inf(error₂);
            auto color     = inner_converged ? "\x1b[0;32m" : "\x1b[0;31m";
            auto color_end = "\x1b[0m";
            std::cout << "[\x1b[0;34mALM\x1b[0m]   " << std::setw(5) << i
                      << ": ‖Σ‖ = " << std::setw(13) << Σ.norm()
                      << ", ‖y‖ = " << std::setw(13) << y.norm()
                      << ", δ = " << std::setw(13) << δ
                      << ", ε = " << std::setw(13) << ps.ε
                      << ", Δ = " << std::setw(13) << Δ
                      << ", status = " << color << std::setw(13) << ps.status
                      << color_end << ", iter = " << std::setw(13)
                      << ps.iterations << "\r\n";
        }

        // TODO: check penalty size?
        if (ps.status == SolverStatus::Interrupted) {
            s.ε                = ps.ε;
            s.δ                = vec_util::norm_inf(error₂);
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i + 1;
            s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
            s.status           = ps.status;
            if (params.preconditioning)
                y = prec_g.asDiagonal() * y / prec_f;
            return s;
        }

        // Backtrack and lower penalty if inner solver did not converge
        if (backtrack) {
            // This means the inner solver didn't produce a solution that
            // satisfies the required tolerance.
            // The best thing we can do now is to restore the penalty to its
            // previous value (when the inner solver did converge), then lower
            // the penalty factor, and update the penalty with this smaller
            // factor.
            // error₂ was not overwritten by the inner solver, so it still
            // contains the error from the iteration before the previous
            // successful iteration. error₁ contains the error of the last
            // successful iteration.
            if (not first_successful_iter) {
                // We have a previous Σ and error
                // Recompute penalty with smaller Δ
                Δ = std::fmax(1., Δ * params.Δ_lower);
                detail::update_penalty_weights(params, Δ, first_successful_iter,
                                               error₁, error₂, norm_e₁, norm_e₂,
                                               Σ_old, Σ);
                // Recompute the primal tolerance with larger ρ
                ρ = std::fmin(0.5, ρ * params.ρ_increase); // keep ρ <= 0.5
                ε = std::fmax(ρ * ε_old, params.ε);
                ++s.penalty_reduced;
            } else {
                // We don't have a previous Σ, simply lower the current Σ and
                // increase ε
                Σ *= params.Σ₀_lower;
                ε *= params.ε₀_increase;
                ++s.initial_penalty_reduced;
            }
        }

        // If the inner solver did converge, increase penalty
        else {
            // After this line, error₁ contains the error of the current
            // (successful) iteration, and error₂ contains the error of the
            // previous successful iteration.
            error₂.swap(error₁);
            norm_e₂ = std::exchange(norm_e₁, vec_util::norm_inf(error₁));

            // Check the termination criteria
            bool alm_converged =
                ps.ε <= params.ε && inner_converged && norm_e₁ <= params.δ;
            bool exit = alm_converged || out_of_iter || out_of_time;
            if (exit) {
                s.ε                = ps.ε;
                s.δ                = norm_e₁;
                s.norm_penalty     = Σ.norm();
                s.outer_iterations = i + 1;
                s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
                s.status           = alm_converged ? SolverStatus::Converged
                                     : out_of_time ? SolverStatus::MaxTime
                                     : out_of_iter ? SolverStatus::MaxIter
                                                   : SolverStatus::Unknown;
                if (params.preconditioning)
                    y = prec_g.asDiagonal() * y / prec_f;
                return s;
            }
            // After this line, Σ_old contains the penalty used in the current
            // (successful) iteration.
            Σ_old.swap(Σ);
            // Update Σ to contain the penalty to use on the next iteration.
            detail::update_penalty_weights(params, Δ, first_successful_iter,
                                           error₁, error₂, norm_e₁, norm_e₂,
                                           Σ_old, Σ);
            // Lower the primal tolerance for the inner solver.
            ε_old = std::exchange(ε, std::fmax(ρ * ε, params.ε));
            first_successful_iter = false;
        }
    }
    throw std::logic_error("[ALM]   loop error");
}

} // namespace alpaqa