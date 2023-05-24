#pragma once

#include <alpaqa/inner/pantr.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <alpaqa/config/config.hpp>
#include <alpaqa/implementation/inner/panoc-helpers.tpp>
#include <alpaqa/implementation/util/print.tpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/quadmath/quadmath-print.hpp>
#include <alpaqa/util/timed.hpp>

namespace alpaqa {

template <class DirectionProviderT>
std::string PANTRSolver<DirectionProviderT>::get_name() const {
    return "PANTRSolver<" + std::string(direction.get_name()) + ">";
}

template <class DirectionProviderT>
auto PANTRSolver<DirectionProviderT>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Solve options
    const SolveOptions &opts,
    /// [inout] Decision variable @f$ x @f$
    rvec x,
    /// [inout] Lagrange multipliers @f$ y @f$
    rvec y,
    /// [in]    Constraint weights @f$ \Sigma @f$
    crvec Σ,
    /// [out]   Slack variable error @f$ g(x) - \Pi_D(g(x) + \Sigma^{-1} y) @f$
    rvec err_z) -> Stats {

    if (opts.check)
        problem.check();

    using std::chrono::nanoseconds;
    auto os         = opts.os ? opts.os : this->os;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto n = problem.get_n();
    const auto m = problem.get_m();

    // Represents an iterate in the algorithm, keeping track of some
    // intermediate values and function evaluations.
    struct Iterate {
        vec x;      //< Decision variables
        vec x̂;      //< Decision variables after proximal gradient step
        vec grad_ψ; //< Gradient of cost in x
        vec p;      //< Proximal gradient step in x
        vec ŷx̂;     //< Candidate Lagrange multipliers in x̂
        real_t ψx       = NaN<config_t>; //< Cost in x
        real_t ψx̂       = NaN<config_t>; //< Cost in x̂
        real_t γ        = NaN<config_t>; //< Step size γ
        real_t L        = NaN<config_t>; //< Lipschitz estimate L
        real_t pᵀp      = NaN<config_t>; //< Norm squared of p
        real_t grad_ψᵀp = NaN<config_t>; //< Dot product of gradient and p
        real_t hx̂       = NaN<config_t>; //< Non-smooth function value in x̂

        // @pre    @ref ψx, @ref hx̂ @ref pᵀp, @ref grad_ψᵀp
        // @return φγ
        real_t fbe() const { return ψx + hx̂ + pᵀp / (2 * γ) + grad_ψᵀp; }

        Iterate(length_t n, length_t m) : x(n), x̂(n), grad_ψ(n), p(n), ŷx̂(m) {}
    } iterates[3]{{n, m}, {n, m}, {n, m}};
    Iterate *curr = &iterates[0];
    Iterate *prox = &iterates[1];
    Iterate *cand = &iterates[2];

    bool need_grad_ψx̂ = Helpers::stop_crit_requires_grad_ψx̂(params.stop_crit);
    vec grad_ψx̂(n);
    vec work_n(n), work_m(m);
    vec q(n); // (quasi-)Newton step Hₖ pₖ
    std::chrono::nanoseconds direction_duration{};

    // Problem functions -------------------------------------------------------

    auto eval_ψ_grad_ψ = [&problem, &y, &Σ, &work_n, &work_m](Iterate &i) {
        i.ψx = problem.eval_ψ_grad_ψ(i.x, y, Σ, i.grad_ψ, work_n, work_m);
    };
    auto eval_prox_grad_step = [&problem](Iterate &i) {
        i.hx̂  = problem.eval_prox_grad_step(i.γ, i.x, i.grad_ψ, i.x̂, i.p);
        i.pᵀp = i.p.squaredNorm();
        i.grad_ψᵀp = i.p.dot(i.grad_ψ);
    };
    auto eval_ψx̂ = [&problem, &y, &Σ](Iterate &i) {
        i.ψx̂ = problem.eval_ψ(i.x̂, y, Σ, i.ŷx̂);
    };
    auto eval_grad_ψx̂ = [&problem, &work_n](Iterate &i, rvec grad_ψx̂) {
        problem.eval_grad_L(i.x̂, i.ŷx̂, grad_ψx̂, work_n);
    };

    // Helper functions --------------------------------------------------------

    auto qub_violated = [this](const Iterate &i) {
        real_t margin =
            (1 + std::abs(i.ψx)) * params.quadratic_upperbound_tolerance_factor;
        return i.ψx̂ > i.ψx + i.grad_ψᵀp + real_t(0.5) * i.L * i.pᵀp + margin;
    };
    auto backtrack_qub = [&](Iterate &i) {
        while (i.L < params.L_max && qub_violated(i)) {
            i.γ /= 2;
            i.L *= 2;
            // Compute x̂, p, ψ(x̂)
            eval_prox_grad_step(i);
            eval_ψx̂(i);
        }
    };

    // Printing ----------------------------------------------------------------

    std::array<char, 64> print_buf;
    auto print_real = [this, &print_buf](real_t x) {
        return float_to_str_vw(print_buf, x, params.print_precision);
    };
    auto print_real3 = [&print_buf](real_t x) {
        return float_to_str_vw(print_buf, x, 3);
    };

    auto print_progress_1 = [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ,
                                real_t pₖᵀpₖ, real_t γₖ, real_t εₖ, real_t Δₖ) {
        if (k == 0)
            *os << "┌─[PANTR]\n";
        else
            *os << "├─ " << std::setw(6) << k << " ──\n";
        *os << "│   φγ = " << print_real(φₖ)               //
            << ",    ψ = " << print_real(ψₖ)               //
            << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())   //
            << ",  ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ)) //
            << ",    γ = " << print_real(γₖ)               //
            << ",    Δ = " << print_real(Δₖ)               //
            << ",    ε = " << print_real(εₖ) << '\n';
    };
    auto print_progress_2 = [&](crvec qₖ, real_t ρₖ, bool accept,
                                std::chrono::nanoseconds direction_duration) {
        *os << "│  ‖q‖ = " << print_real(qₖ.norm()) //
            << ",    ρ = " << print_real3(ρₖ)       //
            << ", time = "
            << print_real3(
                   static_cast<real_t>(1e6) *
                   std::chrono::duration<real_t>(direction_duration).count())
            << " µs, "
            << (accept ? "\033[0;32maccepted\033[0m"
                       : "\033[0;35mrejected\033[0m") //
            << std::endl; // Flush for Python buffering
    };
    auto print_progress_n = [&](SolverStatus status) {
        *os << "└─ " << status << " ──"
            << std::endl; // Flush for Python buffering
    };
    auto do_progress_cb = [this, &s, &problem, &Σ, &y,
                           &opts](unsigned k, Iterate &it, crvec q,
                                  crvec grad_ψx̂, real_t Δ, real_t ρ, real_t εₖ,
                                  SolverStatus status) {
        if (!progress_cb)
            return;
        ScopedMallocAllower ma;
        alpaqa::util::Timed t{s.time_progress_callback};
        progress_cb(ProgressInfo{
            .k          = k,
            .status     = status,
            .x          = it.x,
            .p          = it.p,
            .norm_sq_p  = it.pᵀp,
            .x̂          = it.x̂,
            .φγ         = it.fbe(),
            .ψ          = it.ψx,
            .grad_ψ     = it.grad_ψ,
            .ψ_hat      = it.ψx̂,
            .grad_ψ_hat = grad_ψx̂,
            .q          = q,
            .L          = it.L,
            .γ          = it.γ,
            .Δ          = Δ,
            .ρ          = ρ,
            .ε          = εₖ,
            .Σ          = Σ,
            .y          = y,
            .outer_iter = opts.outer_iter,
            .problem    = &problem,
            .params     = &params,
        });
    };

    // Initialization ----------------------------------------------------------

    curr->x = x;

    // Estimate Lipschitz constant ---------------------------------------------

    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        curr->L = Helpers::initial_lipschitz_estimate(
            problem, curr->x, y, Σ, params.Lipschitz.ε, params.Lipschitz.δ,
            params.L_min, params.L_max,
            /* in ⟹ out */ curr->ψx, curr->grad_ψ, curr->x̂, cand->grad_ψ,
            work_n, work_m);
    }
    // Initial Lipschitz constant provided by the user
    else {
        curr->L = params.Lipschitz.L_0;
        // Calculate ψ(xₖ), ∇ψ(x₀)
        eval_ψ_grad_ψ(*curr);
    }
    if (not std::isfinite(curr->L)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    curr->γ = params.Lipschitz.Lγ_factor / curr->L;

    // First proximal gradient step --------------------------------------------

    eval_prox_grad_step(*curr);
    eval_ψx̂(*curr);
    backtrack_qub(*curr);

    // Loop data ---------------------------------------------------------------

    unsigned k            = 0; // iteration
    bool accept_candidate = false;
    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;
    // Trust radius
    real_t Δ = params.initial_radius;
    if (!std::isfinite(Δ) || Δ == 0)
        Δ = real_t(0.1) * curr->grad_ψ.norm();
    Δ = std::fmax(Δ, params.min_radius);
    // Reduction ratio
    real_t ρ = NaN<config_t>;

    // Main PANTR loop
    // =========================================================================

    ScopedMallocBlocker mb; // Don't allocate in the inner loop
    while (true) {

        // Check stopping criteria ---------------------------------------------

        // Calculate ∇ψ(x̂ₖ)
        if (need_grad_ψx̂)
            eval_grad_ψx̂(*curr, grad_ψx̂);
        bool have_grad_ψx̂ = need_grad_ψx̂;

        real_t εₖ = Helpers::calc_error_stop_crit(
            problem, params.stop_crit, curr->p, curr->γ, curr->x, curr->x̂,
            curr->ŷx̂, curr->grad_ψ, grad_ψx̂, work_n, cand->p);

        // Print progress ------------------------------------------------------

        bool do_print =
            params.print_interval != 0 && k % params.print_interval == 0;
        if (do_print)
            print_progress_1(k, curr->fbe(), curr->ψx, curr->grad_ψ, curr->pᵀp,
                             curr->γ, εₖ, Δ);

        // Return solution -----------------------------------------------------

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status  = Helpers::check_all_stop_conditions(
            params, opts, time_elapsed, k, stop_signal, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            do_progress_cb(k, *curr, null_vec<config_t>, grad_ψx̂, NaN<config_t>,
                           NaN<config_t>, εₖ, stop_status);
            bool do_final_print = params.print_interval != 0;
            if (!do_print && do_final_print)
                print_progress_1(k, curr->fbe(), curr->ψx, curr->grad_ψ,
                                 curr->pᵀp, curr->γ, εₖ, Δ);
            if (do_print || do_final_print)
                print_progress_n(stop_status);
            // Overwrite output arguments
            if (stop_status == SolverStatus::Converged ||
                stop_status == SolverStatus::Interrupted ||
                opts.always_overwrite_results) {
                auto &ŷ = curr->ŷx̂;
                if (err_z.size() > 0)
                    err_z = Σ.asDiagonal().inverse() * (ŷ - y);
                x = std::move(curr->x̂);
                y = std::move(curr->ŷx̂);
            }
            // Save statistics
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.status       = stop_status;
            s.final_γ      = curr->γ;
            s.final_ψ      = curr->ψx̂;
            s.final_h      = curr->hx̂;
            s.final_φγ     = curr->fbe();
            return s;
        }

        // Perform FBS step ----------------------------------------------------

        // x̂ₖ = xₖ + pₖ
        auto compute_FBS_step = [&] {
            assert(curr->L >= params.L_max || !qub_violated(*curr));
            // Calculate ∇ψ(x̂ₖ)
            if (not have_grad_ψx̂)
                eval_grad_ψx̂(*curr, grad_ψx̂);
            have_grad_ψx̂ = true;
            prox->x      = curr->x̂;
            prox->ψx     = curr->ψx̂;
            prox->grad_ψ.swap(grad_ψx̂);
            prox->γ = curr->γ;
            prox->L = curr->L;
            eval_ψ_grad_ψ(*prox);
            eval_prox_grad_step(*prox);
        };

        // store x̂ₖ in prox->x
        compute_FBS_step();

        // Initialize direction
        if (k == 0) {
            ScopedMallocAllower ma;
            direction.initialize(problem, y, Σ, prox->γ, prox->x, prox->x̂,
                                 prox->p, prox->grad_ψ);
        }

        // Check if x̂ₖ + q provides sufficient decrease
        auto compute_candidate_fbe = [&](crvec q) {
            // Candidate step xₖ₊₁ = x̂ₖ + q
            cand->x = prox->x + q;
            // Compute ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            eval_ψ_grad_ψ(*cand);
            cand->γ = prox->γ;
            cand->L = prox->L;
            // Compute x̂ₖ₊₁, pₖ₊₁, ψ(x̂ₖ₊₁)
            eval_prox_grad_step(*cand);

            // Quadratic upper bound in candidate point
            if (params.compute_ratio_using_new_stepsize) {
                eval_ψx̂(*cand);
                backtrack_qub(*cand);
            }
        };

        // Check ratio ρ
        auto compute_candidate_ratio = [this, prox, cand](real_t q_model) {
            real_t ϕγ      = prox->fbe();
            real_t ϕγ_next = cand->fbe();
            real_t margin  = (1 + std::abs(ϕγ)) * params.TR_tolerance_factor;
            real_t ρ       = (ϕγ - ϕγ_next + margin) / (-q_model);
            return params.ratio_approx_fbe_quadratic_model
                       ? ρ / (1 - params.Lipschitz.Lγ_factor)
                       : ρ;
        };

        // update trust radius accordingly
        auto compute_updated_radius = [this](crvec q, real_t ρ, real_t old_Δ) {
            // Very successful TR step
            if (ρ >= params.ratio_threshold_good)
                return std::max(params.radius_factor_good * q.norm(), old_Δ);
            // Successful TR step
            else if (ρ >= params.ratio_threshold_acceptable)
                return old_Δ * params.radius_factor_acceptable;
            // Unsuccessful TR step
            else
                return params.radius_factor_rejected * q.norm();
        };

        // Compute trust region direction from x̂ₖ
        auto compute_trust_region_step = [&](rvec q, real_t Δ) {
            auto t0 = std::chrono::steady_clock::now();
            real_t q_model = direction.apply(prox->γ, prox->x, prox->x̂, prox->p,
                                             prox->grad_ψ, Δ, q);
            auto t1            = std::chrono::steady_clock::now();
            direction_duration = t1 - t0;

            // Check if step is valid
            if (not q.allFinite()) {
                *os << "Direction fail: not finite" << std::endl;
                ++s.direction_failures;
                direction.reset();
                return +inf<config_t>;
            }
            if (q_model >= 0) {
                *os << "Direction fail: no decrease on model (" << q_model
                    << ')' << std::endl;
                ++s.direction_failures;
                direction.reset(); // Is there anything else we can do?
            }
            return q_model;
        };

        // Solve TR subproblem and update radius
        accept_candidate           = false;
        bool accelerated_iteration = k > 0 || direction.has_initial_direction();
        if (accelerated_iteration && !params.disable_acceleration) {
            if (auto q_model = compute_trust_region_step(q, Δ); q_model < 0) {
                compute_candidate_fbe(q);
                ρ                = compute_candidate_ratio(q_model);
                accept_candidate = ρ >= params.ratio_threshold_acceptable;
                Δ                = std::fmax(compute_updated_radius(q, ρ, Δ),
                                             params.min_radius);
            }
        }

        // Progress callback
        do_progress_cb(k, *curr, q, grad_ψx̂, Δ, ρ, εₖ, SolverStatus::Busy);

        // Accept TR step
        if (accept_candidate) {
            // Quadratic upper bound in next iterate
            if (!params.compute_ratio_using_new_stepsize) {
                eval_ψx̂(*cand);
                backtrack_qub(*cand);
            }
            // Flush L-BFGS if γ changed
            if (prox->γ != cand->γ) {
                direction.changed_γ(cand->γ, prox->γ);
                if (params.recompute_last_prox_step_after_direction_reset) {
                    std::tie(prox->γ, prox->L) = std::tie(cand->γ, cand->L);
                    eval_prox_grad_step(*prox);
                }
            }
            // update L-BFGS
            s.direction_update_rejected += not direction.update(
                prox->γ, cand->γ, prox->x, cand->x, prox->p, cand->p,
                prox->grad_ψ, cand->grad_ψ);

            if (do_print)
                print_progress_2(q, ρ, true, direction_duration);
            // Candidate becomes new iterate
            std::swap(curr, cand);
        }
        // Fall back to proximal gradient step
        else {
            if (accelerated_iteration)
                ++s.accelerated_step_rejected;
            // Quadratic upper bound in x̂ₖ
            eval_ψx̂(*prox);
            backtrack_qub(*prox);
            if (prox->γ != curr->γ) {
                direction.changed_γ(prox->γ, curr->γ);
                if (params.recompute_last_prox_step_after_direction_reset) {
                    std::tie(curr->γ, curr->L) = std::tie(prox->γ, prox->L);
                    eval_prox_grad_step(*curr);
                }
            }
            // update direction
            if (params.update_direction_on_prox_step)
                s.direction_update_rejected += not direction.update(
                    curr->γ, prox->γ, curr->x, prox->x, curr->p, prox->p,
                    curr->grad_ψ, prox->grad_ψ);
            if (do_print && accelerated_iteration)
                print_progress_2(q, ρ, false, direction_duration);
            // x̂ₖ becomes new iterate
            std::swap(curr, prox);
        }

#ifndef NDEBUG
        { // Make sure that we don't rely on any data from previous iterations,
            // reset to NaN:
            ScopedMallocAllower ma;
            *prox = {n, m};
            *cand = {n, m};
        }
#endif

        // Advance step --------------------------------------------------------
        ++k;
    }
    throw std::logic_error("[PANTR] loop error");
}

} // namespace alpaqa
