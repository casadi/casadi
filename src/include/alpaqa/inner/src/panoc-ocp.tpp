#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc-ocp/dynamics-eval.hpp>
#include <alpaqa/inner/directions/panoc-ocp/lqr.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/util/index-set.hpp>
#include <alpaqa/util/src/print.tpp>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace alpaqa {

template <Config Conf>
std::string PANOCOCPSolver<Conf>::get_name() const {
    return "PANOCOCPSolver<" + std::string(config_t::get_name()) + '>';
}

template <Config Conf>
auto PANOCOCPSolver<Conf>::operator()(
    /// [in]    Problem description
    const Problem &problem,
    /// [in]    Tolerance @f$ \varepsilon @f$
    const real_t ε,
    /// [inout] Decision variable @f$ u @f$
    rvec u) -> Stats {

    using std::chrono::nanoseconds;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto N   = problem.get_N();
    const auto nu  = problem.get_nu();
    const auto nx  = problem.get_nx();
    const auto n   = nu * N;
    const auto Nxu = n + nx * (N + 1);

    bool enable_lbfgs = params.gn_interval > 1;

    // Allocate storage --------------------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    vec xuₖ(Nxu),                 // x and u at the beginning of the iteration
        x̂uₖ(Nxu),                 // x and u after a projected gradient step
        xuₙₑₓₜ(Nxu),              // xuₖ for next iteration
        x̂uₙₑₓₜ(Nxu),              // x̂uₖ for next iteration
        uₖ(enable_lbfgs ? n : 0), // For L-BFGS update
        uₙₑₓₜ(enable_lbfgs ? n : 0), // For L-BFGS update
        quₖ(enable_lbfgs ? n : 0),   // For L-BFGS update
        pₖ(n), // Projected gradient step pₖ = x̂ₖ - xₖ
        pₙₑₓₜ(n), // Projected gradient step pₙₑₓₜ = x̂ₙₑₓₜ - xₙₑₓₜ
        qxuₖ(Nxu),       // Newton step, including states
        grad_ψₖ(n),      // ∇ψ(xₖ)
        grad_ψₙₑₓₜ(n),   // ∇ψ(xₙₑₓₜ)
        work_p(nx),      //
        work_w(nx + nu); //
    Box<config_t> U{vec::Constant(nu, NaN<config_t>),
                    vec::Constant(n, NaN<config_t>)};

    DynamicsEvaluator<config_t> eval{problem};
    detail::IndexSet<config_t> J{N, nu};
    using LQRFactor = alpaqa::StatefulLQRFactor<config_t>;
    LQRFactor lqr{{.N = N, .nx = nx, .nu = nu}};
    LBFGSParams<config_t> lbfgs_param{.memory = N}; // TODO: make configurable
    LBFGS<config_t> lbfgs{lbfgs_param, enable_lbfgs ? n : 0};

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Helper functions --------------------------------------------------------

    auto eval_proj_grad_step_box = [&U](real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                        rvec p) {
        using binary_real_f = real_t (*)(real_t, real_t);
        p                   = (-γ * grad_ψ)
                .binaryExpr(U.lowerbound - x, binary_real_f(std::fmax))
                .binaryExpr(U.upperbound - x, binary_real_f(std::fmin));
        x̂ = x + p;
    };

    auto eval_prox = [&eval_proj_grad_step_box, &eval, N,
                      nu](real_t γₖ, crvec xuₖ, crvec grad_ψₖ, rvec x̂uₖ,
                          rvec pₖ) {
        real_t pₖᵀpₖ      = 0;
        real_t grad_ψₖᵀpₖ = 0;
        for (index_t t = 0; t < N; ++t) {
            auto &&grad_ψₖ_t = grad_ψₖ.segment(t * nu, nu);
            auto &&pₖ_t      = pₖ.segment(t * nu, nu);
            eval_proj_grad_step_box(γₖ, eval.uk(xuₖ, t), grad_ψₖ_t,
                                    /* in ⟹ out */ eval.uk(x̂uₖ, t), pₖ_t);
            // Calculate ∇ψ(xₖ)ᵀpₖ and ‖pₖ‖²
            pₖᵀpₖ += pₖ_t.squaredNorm();
            grad_ψₖᵀpₖ += grad_ψₖ_t.dot(pₖ_t);
        }
        return std::make_tuple(pₖᵀpₖ, grad_ψₖᵀpₖ);
    };

    auto descent_lemma = [this, &eval, &eval_prox](
                             crvec xuₖ, real_t ψₖ, crvec grad_ψₖ, rvec x̂uₖ,
                             rvec pₖ, real_t &ψx̂ₖ, real_t &pₖᵀpₖ,
                             real_t &grad_ψₖᵀpₖ, real_t &Lₖ, real_t &γₖ) {
        const auto rounding_tolerance =
            params.quadratic_upperbound_tolerance_factor;
        real_t old_γₖ = γₖ;
        real_t margin = (1 + std::abs(ψₖ)) * rounding_tolerance;
        while (ψx̂ₖ - ψₖ > grad_ψₖᵀpₖ + real_t(0.5) * Lₖ * pₖᵀpₖ + margin) {
            if (not(Lₖ * 2 <= params.L_max))
                break;

            Lₖ *= 2;
            γₖ /= 2;

            // Calculate ûₖ and pₖ (with new step size)
            std::tie(pₖᵀpₖ, grad_ψₖᵀpₖ) = eval_prox(γₖ, xuₖ, grad_ψₖ, x̂uₖ, pₖ);

            // Calculate ψ(x̂ₖ)
            ψx̂ₖ = eval.forward(x̂uₖ);
        }
        return old_γₖ;
    };

    auto calc_error_stop_crit = [this, &eval_prox](real_t γ, crvec xuₖ,
                                                   crvec grad_ψₖ, crvec pₖ,
                                                   real_t pₖᵀpₖ, rvec work_xu,
                                                   rvec work_p) {
        switch (params.stop_crit) {
            case PANOCStopCrit::ProjGradNorm: {
                return vec_util::norm_inf(pₖ);
            }
            case PANOCStopCrit::ProjGradNorm2: {
                return std::sqrt(pₖᵀpₖ);
            }
            case PANOCStopCrit::ProjGradUnitNorm: {
                eval_prox(1, xuₖ, grad_ψₖ, work_xu, work_p);
                return vec_util::norm_inf(work_p);
            }
            case PANOCStopCrit::ProjGradUnitNorm2: {
                auto [pTp, gTp] = eval_prox(1, xuₖ, grad_ψₖ, work_xu, work_p);
                return std::sqrt(pTp);
            }
            case PANOCStopCrit::FPRNorm: {
                return vec_util::norm_inf(pₖ) / γ;
            }
            case PANOCStopCrit::FPRNorm2: {
                return std::sqrt(pₖᵀpₖ) / γ;
            }
            case PANOCStopCrit::ApproxKKT: [[fallthrough]];
            case PANOCStopCrit::ApproxKKT2: [[fallthrough]];
            case PANOCStopCrit::Ipopt: [[fallthrough]];
            case PANOCStopCrit::LBFGSBpp: [[fallthrough]];
            default:
                throw std::invalid_argument("Unsupported stopping criterion");
        }
    };

    auto check_all_stop_conditions =
        [this, ε](
            /// [in]    Time elapsed since the start of the algorithm
            auto time_elapsed,
            /// [in]    The current iteration number
            unsigned iteration,
            /// [in]    Tolerance of the current iterate
            real_t εₖ,
            /// [in]    The number of successive iterations no progress was made
            unsigned no_progress) {
            bool out_of_time     = time_elapsed > params.max_time;
            bool out_of_iter     = iteration == params.max_iter;
            bool interrupted     = stop_signal.stop_requested();
            bool not_finite      = not std::isfinite(εₖ);
            bool conv            = εₖ <= ε;
            bool max_no_progress = no_progress > params.max_no_progress;
            return conv              ? SolverStatus::Converged
                   : out_of_time     ? SolverStatus::MaxTime
                   : out_of_iter     ? SolverStatus::MaxIter
                   : not_finite      ? SolverStatus::NotFinite
                   : max_no_progress ? SolverStatus::NoProgress
                   : interrupted     ? SolverStatus::Interrupted
                                     : SolverStatus::Busy;
        };

    auto assign_interleave_xu = [N, nu, nx](crvec u, rvec xu) {
        for (index_t t = 0; t < N; ++t)
            xu.segment(t * (nx + nu) + nx, nu) = u.segment(t * nu, nu);
    };
    auto assign_extract_u = [N, nu, nx](crvec xu, rvec u) {
        for (index_t t = 0; t < N; ++t)
            u.segment(t * nu, nu) = xu.segment(t * (nx + nu) + nx, nu);
    };

    // Note: updates the stats in xuₖ, the Jacobians in eval, and
    // overwrites grad_ψₙₑₓₜ
    auto initial_lipschitz_estimate =
        [&eval, &grad_ψₙₑₓₜ, &work_p, &work_w, N, nu](
            /// [inout]    Current iterate @f$ x^k @f$
            rvec xuₖ,
            /// Whether to compute the Jacobians of the dynamics in xuₖ
            bool do_gn_step,
            /// [in]    Finite difference step size relative to xₖ
            real_t ε,
            /// [in]    Minimum absolute finite difference step size
            real_t δ,
            /// [in]    Minimum allowed Lipschitz estimate.
            real_t L_min,
            /// [in]    Maximum allowed Lipschitz estimate.
            real_t L_max,
            /// [out]   @f$ \psi(x^k) @f$
            real_t &ψ,
            /// [out]   Gradient @f$ \nabla \psi(x^k) @f$
            rvec grad_ψ,
            ///         Workspace with the same dimensions as xuₖ
            rvec work_xuₖ) {
            // Calculate ψ(x₀), ∇ψ(x₀)
            ψ = eval.forward(xuₖ);
            do_gn_step ? eval.backward_with_jac(xuₖ, grad_ψ, work_p)
                       : eval.backward(xuₖ, grad_ψ, work_p, work_w);
            // Select a small step h for finite differences
            auto h = (-grad_ψ * ε).cwiseAbs().cwiseMax(δ); // TODO: remove abs
            real_t norm_h = h.norm();
            for (index_t t = 0; t < N; ++t)
                eval.uk(work_xuₖ, t) = eval.uk(xuₖ, t) + h.segment(t * nu, nu);
            eval.forward_simulate(work_xuₖ); // needed for backwards sweep
            // Calculate ∇ψ(x₀ + h)
            eval.backward(work_xuₖ, grad_ψₙₑₓₜ, work_p, work_w);

            // Estimate Lipschitz constant using finite differences
            real_t L = (grad_ψₙₑₓₜ - grad_ψ).norm() / norm_h;
            return std::clamp(L, L_min, L_max);
        };

    std::array<char, 64> print_buf;
    auto print_real = [&](real_t x) {
        return float_to_str_vw(print_buf, x, params.print_precision);
    };
    auto print_real3 = [&](real_t x) {
        return float_to_str_vw(print_buf, x, 3);
    };
    auto print_progress = [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ,
                              real_t pₖᵀpₖ, crvec qₖ, real_t γₖ, real_t τₖ,
                              real_t εₖ, bool did_gn, length_t nJ,
                              real_t min_rcond) {
        std::cout << "[PANOC] " << std::setw(6) << k
                  << ": φγ = " << print_real(φₖ) << ", ψ = " << print_real(ψₖ)
                  << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())
                  << ", ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ))
                  << ", γ = " << print_real(γₖ) << ", εₖ = " << print_real(εₖ);
        if (k > 0)
            std::cout << ", τ = " << print_real3(τₖ)
                      << ", ‖q‖ = " << print_real(qₖ.norm())
                      << ", #J =" << std::setw(5) << nJ
                      << ", rcond = " << print_real3(min_rcond) << ", "
                      << (did_gn ? "GN" : "L-BFGS");
        std::cout << std::endl; // Flush for Python buffering
    };

    // Initialize inputs and initial state (do not simulate states yet) --------

    assign_interleave_xu(u, xuₖ);
    problem.get_x_init(xuₖ.topRows(nx));
    x̂uₖ.topRows(nx)    = xuₖ.topRows(nx);
    xuₙₑₓₜ.topRows(nx) = xuₖ.topRows(nx);
    x̂uₙₑₓₜ.topRows(nx) = xuₖ.topRows(nx);
    qxuₖ.topRows(nx).setZero();
    if (enable_lbfgs)
        uₖ = u;

    problem.get_U(U);

    bool do_gn_step     = true;
    bool have_jacobians = false; // TODO: just as a sanity check
    bool did_gn         = false;

    // Estimate Lipschitz constant ---------------------------------------------

    real_t ψₖ, Lₖ;
    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        Lₖ = initial_lipschitz_estimate(xuₖ, do_gn_step, params.Lipschitz.ε,
                                        params.Lipschitz.δ, params.L_min,
                                        params.L_max,
                                        /* in ⟹ out */ ψₖ, grad_ψₖ, xuₙₑₓₜ);
        have_jacobians = do_gn_step;
    }
    // Initial Lipschitz constant provided by the user
    else {
        Lₖ = params.Lipschitz.L_0;
        // Calculate ψ(x₀), ∇ψ(x₀)
        ψₖ = eval.forward(xuₖ);
        do_gn_step ? eval.backward_with_jac(xuₖ, grad_ψₖ, work_p)
                   : eval.backward(xuₖ, grad_ψₖ, work_p, work_w);
        have_jacobians = do_gn_step;
    }
    if (not std::isfinite(Lₖ)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    real_t γₖ = params.Lipschitz.Lγ_factor / Lₖ;
    real_t τ  = NaN<config_t>;

    // First projected gradient step -------------------------------------------

    real_t pₖᵀpₖ;
    real_t grad_ψₖᵀpₖ;

    // Calculate x̂_0, p_0 (projected gradient step)
    std::tie(pₖᵀpₖ, grad_ψₖᵀpₖ) = eval_prox(γₖ, xuₖ, grad_ψₖ, x̂uₖ, pₖ);
    // Calculate ψ(x̂ₖ)
    real_t ψx̂ₖ = eval.forward(x̂uₖ);
    // Compute forward-backward envelope
    real_t φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;

    // Make sure that we don't allocate any memory in the inner loop
    ScopedMallocBlocker mb;

    // Main PANOC loop
    // =========================================================================
    for (unsigned k = 0; k <= params.max_iter; ++k) {

        // Quadratic upper bound -----------------------------------------------
        if (k == 0 || params.update_lipschitz_in_linesearch == false) {
            // Decrease step size until quadratic upper bound is satisfied
            real_t old_γₖ = descent_lemma(xuₖ, ψₖ, grad_ψₖ, x̂uₖ, pₖ, ψx̂ₖ, pₖᵀpₖ,
                                          grad_ψₖᵀpₖ, Lₖ, γₖ);
            if (γₖ != old_γₖ)
                φₖ = ψₖ + 1 / (2 * γₖ) * pₖᵀpₖ + grad_ψₖᵀpₖ;
        }

        // Check stop condition ------------------------------------------------
        real_t εₖ =
            calc_error_stop_crit(γₖ, xuₖ, grad_ψₖ, pₖ, pₖᵀpₖ, x̂uₙₑₓₜ, pₙₑₓₜ);

        // Print progress
        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, φₖ, ψₖ, grad_ψₖ, pₖᵀpₖ, qxuₖ, γₖ, τ, εₖ, did_gn,
                           J.sizes().sum(), lqr.min_rcond);
        if (progress_cb) {
            ScopedMallocAllower ma;
            progress_cb({k, xuₖ, pₖ, pₖᵀpₖ, x̂uₖ, φₖ, ψₖ, grad_ψₖ, ψx̂ₖ, qxuₖ,
                         did_gn, J.sizes().sum(), lqr.min_rcond, Lₖ, γₖ, τ, εₖ,
                         problem, params});
        }

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status =
            check_all_stop_conditions(time_elapsed, k, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            // TODO: loop
            assign_extract_u(x̂uₖ, u);
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.status       = stop_status;
            s.final_γ      = γₖ;
            return s;
        }

        // Line search initialization ------------------------------------------
        τ                = 1;
        real_t σₖγₖpₖᵀpₖ = (1 - γₖ * Lₖ) * pₖᵀpₖ / (2 * γₖ);
        real_t φₙₑₓₜ, ψₙₑₓₜ, ψx̂ₙₑₓₜ, grad_ψₙₑₓₜᵀpₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ;
        real_t Lₙₑₓₜ, γₙₑₓₜ;
        real_t ls_cond;
        // TODO: make separate parameter
        real_t margin =
            (1 + std::abs(φₖ)) * params.quadratic_upperbound_tolerance_factor;

        // Calculate Gauss-Newton step -----------------------------------------

        did_gn = do_gn_step;
        if (do_gn_step) {
            // Forward simulation has been carried out on xuₖ and the Jacobians
            // are stored in 'eval'
            if (!have_jacobians)
                throw std::logic_error("have_jacobians");

            // Make sure that Q and R are stored as well
            eval.hessians(xuₖ);

            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui = eval.uk(xuₖ, t)(i);
                // Gradient descent step.
                real_t gs = ui - γₖ * grad_ψₖ(t * nu + i);
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub) {
                    eval.uk(qxuₖ, t)(i) = U.upperbound(i) - ui;
                    return false;
                } else if (active_lb) {
                    eval.uk(qxuₖ, t)(i) = U.lowerbound(i) - ui;
                    return false;
                } else { // Store inactive indices
                    return true;
                }
            };
            J.update(is_constr_inactive);

            auto Ak    = [&](index_t i) -> crmat { return eval.Ak(i); };
            auto Bk    = [&](index_t i) -> crmat { return eval.Bk(i); };
            auto Qk    = [&](index_t k) -> crmat { return eval.Qk(k); };
            auto Rk    = [&](index_t k) -> crmat { return eval.Rk(k); };
            auto qk    = [&](index_t k) -> crvec { return eval.qk(k); };
            auto rk    = [&](index_t k) -> crvec { return eval.rk(k); };
            auto uk_eq = [&](index_t k) -> crvec { return eval.uk(qxuₖ, k); };
            auto Jk    = [&](index_t k) -> crindexvec { return J.indices(k); };
            auto Kk    = [&](index_t k) -> crindexvec {
                return J.compl_indices(k);
            };

            lqr.factor_masked(Ak, Bk, Qk, Rk, qk, rk, uk_eq, Jk, Kk,
                              params.lqr_factor_cholesky);
            lqr.solve_masked(Ak, Bk, Jk, qxuₖ);
        } else {
            if (!enable_lbfgs)
                throw std::logic_error("enable_lbfgs");

            // Find inactive indices J
            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui     = eval.uk(xuₖ, t)(i);
                real_t grad_i = grad_ψₖ(t * nu + i);
                // Gradient descent step.
                real_t gs = ui - γₖ * grad_i;
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub) {
                    quₖ(t * nu + i) = pₖ(t * nu + i);
                    return false;
                } else if (active_lb) {
                    quₖ(t * nu + i) = pₖ(t * nu + i);
                    return false;
                } else { // Store inactive indices
                    quₖ(t * nu + i) = -grad_i;
                    return true;
                }
            };

            auto J_idx = J.indices();
            J_idx.setConstant(0x5A5A'5A5A);
            index_t nJ = 0;
            for (index_t t = 0; t < N; ++t)
                for (index_t i = 0; i < nu; ++i)
                    if (is_constr_inactive(t, i))
                        J_idx(nJ++) = t * nu + i;

            auto J_lbfgs = J_idx.topRows(nJ);

            // If all indices are inactive, we can use standard L-BFGS,
            // if there are active indices, we need the specialized version
            // that only applies L-BFGS to the inactive indices
            bool success = lbfgs.apply_masked(quₖ, γₖ, J_lbfgs);
            // If L-BFGS application failed, qₖ(J) still contains
            // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
            if (not success) {
                τ = 0;
            } else {
                // for (index_t t = N; t-- > 0;)
                //     for (index_t i = nu; i-- > 0;)
                //         qxuₖ(t * (nx + nu) + nx + i) = qxuₖ(t * (0 + nu) + nx + i);
                assign_interleave_xu(quₖ, qxuₖ);
            }
        }

        // Make sure quasi-Newton step is valid
        if (not qxuₖ.allFinite()) {
            τ = 0;
            ++s.lbfgs_failures;
            // Is there anything we can do?
            if (not did_gn)
                lbfgs.reset();
        }

        bool do_next_gn =
            params.gn_interval == 0 || ((k + 1) % params.gn_interval) == 0;
        do_gn_step = do_next_gn || (do_gn_step && params.gn_sticky);

        // Line search loop ----------------------------------------------------
        do {
            Lₙₑₓₜ = Lₖ;
            γₙₑₓₜ = γₖ;

            // Calculate xₙₑₓₜ
            if (τ / 2 < params.τ_min) { // line search failed
                xuₙₑₓₜ.swap(x̂uₖ);       // → safe prox step
                ψₙₑₓₜ = ψx̂ₖ;
                do_gn_step ? eval.backward_with_jac(xuₙₑₓₜ, grad_ψₙₑₓₜ, work_p)
                           : eval.backward(xuₙₑₓₜ, grad_ψₙₑₓₜ, work_p, work_w);
                have_jacobians = do_gn_step;
            } else {          // line search didn't fail (yet)
                if (τ == 1) { // → faster quasi-Newton step
                    xuₙₑₓₜ = xuₖ + qxuₖ;
                } else {
                    do_gn_step = do_next_gn;
                    for (index_t t = 0; t < N; ++t)
                        eval.uk(xuₙₑₓₜ, t) = eval.uk(xuₖ, t) +
                                             (1 - τ) * pₖ.segment(t * nu, nu) +
                                             τ * eval.uk(qxuₖ, t);
                }
                // Calculate ψ(xₙₑₓₜ), ∇ψ(xₙₑₓₜ)
                ψₙₑₓₜ = eval.forward(xuₙₑₓₜ);
                do_gn_step ? eval.backward_with_jac(xuₙₑₓₜ, grad_ψₙₑₓₜ, work_p)
                           : eval.backward(xuₙₑₓₜ, grad_ψₙₑₓₜ, work_p, work_w);
                have_jacobians = do_gn_step;
            }

            // Calculate x̂ₙₑₓₜ, pₙₑₓₜ (projected gradient step in xₙₑₓₜ)
            eval_prox(γₙₑₓₜ, xuₙₑₓₜ, grad_ψₙₑₓₜ, /* in ⟹ out */ x̂uₙₑₓₜ, pₙₑₓₜ);
            // Calculate ψ(x̂ₙₑₓₜ)
            ψx̂ₙₑₓₜ = eval.forward(x̂uₙₑₓₜ);

            // Quadratic upper bound -------------------------------------------
            grad_ψₙₑₓₜᵀpₙₑₓₜ = grad_ψₙₑₓₜ.dot(pₙₑₓₜ);
            pₙₑₓₜᵀpₙₑₓₜ      = pₙₑₓₜ.squaredNorm();

            if (params.update_lipschitz_in_linesearch == true) {
                // Decrease step size until quadratic upper bound is satisfied
                (void)descent_lemma(xuₙₑₓₜ, ψₙₑₓₜ, grad_ψₙₑₓₜ,
                                    /* in ⟹ out */ x̂uₙₑₓₜ, pₙₑₓₜ,
                                    /* inout */ ψx̂ₙₑₓₜ, pₙₑₓₜᵀpₙₑₓₜ,
                                    grad_ψₙₑₓₜᵀpₙₑₓₜ, Lₙₑₓₜ, γₙₑₓₜ);
            }

            // Compute forward-backward envelope
            φₙₑₓₜ = ψₙₑₓₜ + 1 / (2 * γₙₑₓₜ) * pₙₑₓₜᵀpₙₑₓₜ + grad_ψₙₑₓₜᵀpₙₑₓₜ;
            // Compute line search condition
            ls_cond = φₙₑₓₜ - (φₖ - σₖγₖpₖᵀpₖ);

            τ /= 2;
        } while (ls_cond > margin && τ >= params.τ_min);

        // If τ < τ_min the line search failed and we accepted the prox step
        if (τ < params.τ_min) {
            ++s.linesearch_failures;
            τ = 0;
        }
        τ *= 2; // restore to the value that was actually accepted
        s.count_τ += 1;
        s.sum_τ += τ;
        s.τ_1_accepted += τ == 1;

        // Update L-BFGS -------------------------------------------------------

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = xuₖ == xuₙₑₓₜ ? no_progress + 1 : 0;

        // Update L-BFGS
        if (enable_lbfgs) {
            const bool force = true;
            assign_extract_u(xuₙₑₓₜ, uₙₑₓₜ);
            if (did_gn && params.reset_lbfgs_on_gn_step) {
                lbfgs.reset();
            } else {
                s.lbfgs_rejected +=
                    not lbfgs.update(uₖ, uₙₑₓₜ, grad_ψₖ, grad_ψₙₑₓₜ,
                                     LBFGS<config_t>::Sign::Positive, force);
            }
        }

        // Advance step --------------------------------------------------------
        Lₖ = Lₙₑₓₜ;
        γₖ = γₙₑₓₜ;

        ψₖ  = ψₙₑₓₜ;
        ψx̂ₖ = ψx̂ₙₑₓₜ;
        φₖ  = φₙₑₓₜ;

        xuₖ.swap(xuₙₑₓₜ);
        x̂uₖ.swap(x̂uₙₑₓₜ);
        pₖ.swap(pₙₑₓₜ);
        grad_ψₖ.swap(grad_ψₙₑₓₜ);
        grad_ψₖᵀpₖ = grad_ψₙₑₓₜᵀpₙₑₓₜ;
        pₖᵀpₖ      = pₙₑₓₜᵀpₙₑₓₜ;
        if (enable_lbfgs)
            uₖ.swap(uₙₑₓₜ);
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa