#pragma once

#include <iomanip>
#include <iostream>
#include <panoc-alm/detail/alm-helpers.hpp>

namespace pa {

using std::chrono::duration_cast;
using std::chrono::microseconds;

template <class InnerSolverT>
typename ALMSolver<InnerSolverT>::Stats
ALMSolver<InnerSolverT>::operator()(const Problem &problem, vec &y, vec &x) {
    auto start_time = std::chrono::steady_clock::now();

    constexpr auto sigNaN = std::numeric_limits<real_t>::signaling_NaN();
    vec Σ                 = vec::Constant(problem.m, sigNaN);
    vec error             = vec::Constant(problem.m, sigNaN);
    vec error_old         = vec::Constant(problem.m, sigNaN);

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

    real_t ε = params.ε₀;

    for (unsigned int i = 0; i < params.max_iter; ++i) {
        detail::project_y(y, p.D.lowerbound, p.D.upperbound, params.M);
        auto ps = inner_solver(p, Σ, ε, x, y, error);
        s.inner_iterations += ps.iterations;
        s.inner_linesearch_failures += ps.linesearch_failures;
        s.inner_lbfgs_failures += ps.lbfgs_failures;
        s.inner_lbfgs_rejected += ps.lbfgs_rejected;
        s.inner_convergence_failures += ps.status != SolverStatus::Converged;
        real_t norm_e = vec_util::norm_inf(error);

        if (params.print_interval != 0 && i % params.print_interval == 0) {
            std::cout << "[\x1b[0;34mALM\x1b[0m]   " << std::setw(5) << i
                      << ": ‖Σ‖ = " << std::setw(13) << Σ.norm()
                      << ", ‖y‖ = " << std::setw(13) << y.norm()
                      << ", δ = " << std::setw(13) << norm_e
                      << ", ε = " << std::setw(13) << ps.ε << "\r\n";
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        // TODO: check penalty size?
        if (ps.status == SolverStatus::NotFinite ||
            ps.status == SolverStatus::Interrupted) {
            s.ε                = ps.ε;
            s.δ                = norm_e;
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i;
            s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
            s.status           = ps.status;
            if (params.preconditioning)
                y = prec_g.asDiagonal() * y / prec_f;
            return s;
        }

        bool converged = ps.ε <= params.ε &&
                         ps.status == SolverStatus::Converged &&
                         norm_e <= params.δ;
        bool out_of_time = time_elapsed > params.max_time;
        if (converged || i + 1 == params.max_iter || out_of_time) {
            s.ε                = ps.ε;
            s.δ                = norm_e;
            s.norm_penalty     = Σ.norm();
            s.outer_iterations = i + 1;
            s.elapsed_time     = duration_cast<microseconds>(time_elapsed);
            s.status           = converged     ? SolverStatus::Converged
                                 : out_of_time ? SolverStatus::MaxTime
                                               : SolverStatus::MaxIter;
            if (params.preconditioning)
                y = prec_g.asDiagonal() * y / prec_f;
            return s;
        }
        detail::update_penalty_weights(params, i, error, error_old, norm_e, Σ);
        ε = std::fmax(params.ρ * ε, params.ε);
        error_old.swap(error);
    }
    throw std::logic_error("[ALM]   loop error");
}

} // namespace pa