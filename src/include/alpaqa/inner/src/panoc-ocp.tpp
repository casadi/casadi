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

namespace detail {

template <Config Conf>
void assign_interleave_xu(OCPDim<Conf> dim, crvec<Conf> u, rvec<Conf> xu) {
    for (index_t<Conf> t = 0; t < dim.N; ++t)
        xu.segment(t * (dim.nx + dim.nu) + dim.nx, dim.nu) =
            u.segment(t * dim.nu, dim.nu);
}
template <Config Conf>
void assign_interleave_xu(OCPDim<Conf> dim, crvec<Conf> x, crvec<Conf> u,
                          rvec<Conf> xu) {
    for (index_t<Conf> t = 0; t < dim.N; ++t) {
        xu.segment(t * (dim.nx + dim.nu), dim.nx) =
            x.segment(t * dim.nx, dim.nx);
        xu.segment(t * (dim.nx + dim.nu) + dim.nx, dim.nu) =
            u.segment(t * dim.nu, dim.nu);
    }
    xu.segment(dim.N * (dim.nx + dim.nu), dim.nx) =
        x.segment(dim.N * dim.nx, dim.nx);
}
template <Config Conf>
void assign_extract_u(OCPDim<Conf> dim, crvec<Conf> xu, rvec<Conf> u) {
    for (index_t<Conf> t = 0; t < dim.N; ++t)
        u.segment(t * dim.nu, dim.nu) =
            xu.segment(t * (dim.nx + dim.nu) + dim.nx, dim.nu);
}
template <Config Conf>
void assign_extract_x(OCPDim<Conf> dim, crvec<Conf> xu, rvec<Conf> x) {
    for (index_t<Conf> t = 0; t < dim.N + 1; ++t)
        x.segment(t * dim.nx, dim.nx) =
            xu.segment(t * (dim.nx + dim.nu), dim.nx);
}

template <Config Conf>
vec<Conf> extract_u(const TypeErasedControlProblem<Conf> &problem,
                    crvec<Conf> xu) {
    auto dim = problem.get_dimensions();
    vec<Conf> u(dim.N * dim.nu);
    assign_extract_u(dim, xu, u);
    return u;
}
template <Config Conf>
vec<Conf> extract_x(const TypeErasedControlProblem<Conf> &problem,
                    crvec<Conf> xu) {
    auto dim = problem.get_dimensions();
    vec<Conf> x((dim.N + 1) * (dim.nx + dim.nu));
    assign_extract_x(dim, xu, x);
    return x;
}

} // namespace detail

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::u() const -> vec {
    return detail::extract_u(problem, xu);
}

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::x() const -> vec {
    return detail::extract_x(problem, xu);
}

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::û() const -> vec {
    return detail::extract_u(problem, x̂u);
}

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::x̂() const -> vec {
    return detail::extract_x(problem, x̂u);
}

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::qu() const -> vec {
    return detail::extract_u(problem, q);
}

template <Config Conf>
auto PANOCOCPProgressInfo<Conf>::qx() const -> vec {
    return detail::extract_x(problem, q);
}

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

    bool enable_lbfgs = params.gn_interval != 1;

    // Allocate storage --------------------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    vec quₖ(enable_lbfgs ? n : 0), // For L-BFGS update
        qxuₖ(Nxu),                 // Newton step, including states
        work_p(nx),                //
        work_w(nx + nu);           //
    Box<config_t> U = Box<config_t>::NaN(nu);

    DynamicsEvaluator<config_t> eval{problem};
    detail::IndexSet<config_t> J{N, nu};
    using LQRFactor = alpaqa::StatefulLQRFactor<config_t>;
    LQRFactor lqr{{.N = N, .nx = nx, .nu = nu}};
    LBFGSParams<config_t> lbfgs_param{.memory = N}; // TODO: make configurable
    LBFGS<config_t> lbfgs{lbfgs_param, enable_lbfgs ? n : 0};

    // Functions for accessing the LQR matrices and index sets
    auto Ak    = [&](index_t i) -> crmat { return eval.Ak(i); };
    auto Bk    = [&](index_t i) -> crmat { return eval.Bk(i); };
    auto Qk    = [&](index_t k) -> crmat { return eval.Qk(k); };
    auto Rk    = [&](index_t k) -> crmat { return eval.Rk(k); };
    auto qk    = [&](index_t k) -> crvec { return eval.qk(k); };
    auto rk    = [&](index_t k) -> crvec { return eval.rk(k); };
    auto uk_eq = [&](index_t k) -> crvec { return eval.uk(qxuₖ, k); };
    auto Jk    = [&](index_t k) -> crindexvec { return J.indices(k); };
    auto Kk    = [&](index_t k) -> crindexvec { return J.compl_indices(k); };

    // Iterates ----------------------------------------------------------------

    struct Iterate {
        vec xu;     ///< Inputs u interleaved with states x
        vec xû;     ///< Inputs u interleaved with states x after prox grad
        vec grad_ψ; ///< Gradient of cost in u
        vec p;      ///< Proximal gradient step in u
        vec u;      ///< Inputs u (used for L-BFGS only)
        real_t ψu           = NaN<config_t>; ///< Cost in u
        real_t ψû           = NaN<config_t>; ///< Cost in û
        real_t γ            = NaN<config_t>; ///< Step size γ
        real_t L            = NaN<config_t>; ///< Lipschitz estimate L
        real_t pᵀp          = NaN<config_t>; ///< Norm squared of p
        real_t grad_ψᵀp     = NaN<config_t>; ///< Dot product of gradient and p
        bool have_jacobians = false;

        /// @pre    @ref ψu, @ref pᵀp, @pre grad_ψᵀp
        /// @return φγ
        real_t fbe() const { return ψu + pᵀp / (2 * γ) + grad_ψᵀp; }

        Iterate(OCPDim<config_t> dim, bool enable_lbfgs)
            : xu{dim.N * (dim.nx + dim.nu) + dim.nx},
              xû{dim.N * (dim.nx + dim.nu) + dim.nx}, grad_ψ{dim.N * dim.nu},
              p{dim.N * dim.nu}, u{enable_lbfgs ? dim.N * dim.nu : 0} {}
    } iterates[2]{
        {problem.get_dimensions(), enable_lbfgs},
        {problem.get_dimensions(), enable_lbfgs},
    };
    Iterate *curr = &iterates[0];
    Iterate *next = &iterates[1];

    // Helper functions --------------------------------------------------------

    auto eval_proj_grad_step_box = [&U](real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                        rvec p) {
        using binary_real_f = real_t (*)(real_t, real_t);
        p                   = (-γ * grad_ψ)
                .binaryExpr(U.lowerbound - x, binary_real_f(std::fmax))
                .binaryExpr(U.upperbound - x, binary_real_f(std::fmin));
        x̂ = x + p;
    };

    auto eval_prox_impl = [&eval_proj_grad_step_box, &eval, N,
                           nu](real_t γ, crvec xu, crvec grad_ψ, rvec x̂u,
                               rvec p) {
        real_t pᵀp      = 0;
        real_t grad_ψᵀp = 0;
        for (index_t t = 0; t < N; ++t) {
            auto &&grad_ψ_t = grad_ψ.segment(t * nu, nu);
            auto &&p_t      = p.segment(t * nu, nu);
            eval_proj_grad_step_box(γ, eval.uk(xu, t), grad_ψ_t,
                                    /* in ⟹ out */ eval.uk(x̂u, t), p_t);
            // Calculate ∇ψ(x)ᵀp and ‖p‖²
            pᵀp += p_t.squaredNorm();
            grad_ψᵀp += grad_ψ_t.dot(p_t);
        }
        return std::make_tuple(pᵀp, grad_ψᵀp);
    };

    auto calc_error_stop_crit = [this, &eval_prox_impl](
                                    real_t γ, crvec xuₖ, crvec grad_ψₖ,
                                    crvec pₖ, real_t pₖᵀpₖ, rvec work_xu,
                                    rvec work_p) {
        switch (params.stop_crit) {
            case PANOCStopCrit::ProjGradNorm: {
                return vec_util::norm_inf(pₖ);
            }
            case PANOCStopCrit::ProjGradNorm2: {
                return std::sqrt(pₖᵀpₖ);
            }
            case PANOCStopCrit::ProjGradUnitNorm: {
                eval_prox_impl(1, xuₖ, grad_ψₖ, work_xu, work_p);
                return vec_util::norm_inf(work_p);
            }
            case PANOCStopCrit::ProjGradUnitNorm2: {
                auto [pTp, gTp] =
                    eval_prox_impl(1, xuₖ, grad_ψₖ, work_xu, work_p);
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

    auto assign_interleave_xu = [dim{problem.get_dimensions()}](crvec u,
                                                                rvec xu) {
        detail::assign_interleave_xu(dim, u, xu);
    };
    auto assign_extract_u = [dim{problem.get_dimensions()}](crvec xu, rvec u) {
        detail::assign_extract_u(dim, xu, u);
    };

    /// @pre    @ref Iterate::γ, @ref Iterate::xu, @ref Iterate::grad_ψ
    /// @post   @ref Iterate::xû, @ref Iterate::p, @ref Iterate::pᵀp,
    ///         @ref Iterate::grad_ψᵀp
    auto eval_prox = [&eval_prox_impl](Iterate &i) {
        std::tie(i.pᵀp, i.grad_ψᵀp) =
            eval_prox_impl(i.γ, i.xu, i.grad_ψ, i.xû, i.p);
    };

    /// @pre    @ref Iterate::xu
    /// @post   @ref Iterate::ψu
    auto eval_forward = [&eval](Iterate &i) { i.ψu = eval.forward(i.xu); };
    /// @pre    @ref Iterate::xû
    /// @post   @ref Iterate::ψû
    auto eval_forward_hat = [&eval](Iterate &i) { i.ψû = eval.forward(i.xû); };

    /// @pre    @ref Iterate::xu
    /// @post   @ref Iterate::grad_ψ, @ref Iterate::have_jacobians
    auto eval_backward = [&eval, &work_p, &work_w](Iterate &i, bool with_jac) {
        with_jac ? eval.backward_with_jac(i.xu, i.grad_ψ, work_p)
                 : eval.backward(i.xu, i.grad_ψ, work_p, work_w);
        i.have_jacobians = with_jac;
    };

    auto qub_violated = [this](const Iterate &i) {
        real_t margin =
            (1 + std::abs(i.ψu)) * params.quadratic_upperbound_tolerance_factor;
        return i.ψû > i.ψu + i.grad_ψᵀp + real_t(0.5) * i.L * i.pᵀp + margin;
    };

    auto linesearch_violated = [this](const Iterate &curr,
                                      const Iterate &next) {
        real_t σ  = params.β * (1 - curr.γ * curr.L) / (2 * curr.γ);
        real_t φγ = curr.fbe();
        real_t margin = (1 + std::abs(φγ)) * params.linesearch_tolerance_factor;
        return next.fbe() > φγ - σ * curr.pᵀp + margin;
    };

    auto initial_lipschitz_estimate =
        [&eval, &eval_forward, &eval_backward, &work_p, &work_w, N, nu](
            /// Iterate, updates xu, ψ, grad_ψ, have_jacobians, L
            Iterate *it,
            /// Whether to compute the Jacobians of the dynamics in xu
            bool do_gn_step,
            /// [in]    Finite difference step size relative to x
            real_t ε,
            /// [in]    Minimum absolute finite difference step size
            real_t δ,
            /// [in]    Minimum allowed Lipschitz estimate.
            real_t L_min,
            /// [in]    Maximum allowed Lipschitz estimate.
            real_t L_max,
            ///         Workspace with the same dimensions as xu, with x_init
            rvec work_xu,
            ///         Workspace with the same dimensions as grad_ψ
            rvec work_grad_ψ) {
            // Calculate ψ(x₀), ∇ψ(x₀)
            eval_forward(*it);
            eval_backward(*it, do_gn_step);
            // Select a small step h for finite differences
            // TODO: remove abs
            auto h        = (it->grad_ψ * -ε).cwiseAbs().cwiseMax(δ);
            real_t norm_h = h.norm();
            // work_xu = xu + h
            for (index_t t = 0; t < N; ++t)
                eval.uk(work_xu, t) =
                    eval.uk(it->xu, t) + h.segment(t * nu, nu);
            // Calculate ψ(x₀ + h)
            eval.forward_simulate(work_xu); // needed for backwards sweep
            // Calculate ∇ψ(x₀ + h)
            eval.backward(work_xu, work_grad_ψ, work_p, work_w);

            // Estimate Lipschitz constant using finite differences
            it->L = (work_grad_ψ - it->grad_ψ).norm() / norm_h;
            it->L = std::clamp(it->L, L_min, L_max);
        };

    // Printing ----------------------------------------------------------------

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

    assign_interleave_xu(u, curr->xu);           // initial guess
    problem.get_x_init(curr->xu.topRows(nx));    // initial state
    curr->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xu.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
    qxuₖ.setZero();                              // set the states to zero
    if (enable_lbfgs)
        curr->u = u;

    problem.get_U(U); // input box constraints

    bool do_gn_step = params.gn_interval > 0 and !params.disable_acceleration;
    bool did_gn     = false;

    // Make sure that we don't allocate any memory in the inner loop
    ScopedMallocBlocker mb;

    // Estimate Lipschitz constant ---------------------------------------------

    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        initial_lipschitz_estimate(curr, do_gn_step, params.Lipschitz.ε,
                                   params.Lipschitz.δ, params.L_min,
                                   params.L_max, next->xu, next->grad_ψ);
    }
    // Initial Lipschitz constant provided by the user
    else {
        curr->L = params.Lipschitz.L_0;
        // Calculate ψ(x₀), ∇ψ(x₀)
        eval_forward(*curr);
        eval_backward(*curr, do_gn_step);
    }
    if (not std::isfinite(curr->L)) {
        s.status = SolverStatus::NotFinite;
        return s;
    }
    curr->γ = params.Lipschitz.Lγ_factor / curr->L;
    eval_prox(*curr);
    eval_forward_hat(*curr);

    unsigned k  = 0;
    real_t τ    = NaN<config_t>;
    length_t nJ = -1;

    // Keep track of how many successive iterations didn't update the iterate
    unsigned no_progress = 0;

    // Main PANOC loop
    // =========================================================================
    while (true) {

        // Check stop condition ------------------------------------------------

        real_t εₖ = calc_error_stop_crit(curr->γ, curr->xu, curr->grad_ψ,
                                         curr->p, curr->pᵀp, next->xû, next->p);

        // Print progress ------------------------------------------------------

        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, curr->fbe(), curr->ψu, curr->grad_ψ, curr->pᵀp,
                           qxuₖ, curr->γ, τ, εₖ, did_gn,
                           did_gn ? J.sizes().sum() : nJ, lqr.min_rcond);
        if (progress_cb) {
            ScopedMallocAllower ma;
            progress_cb({.k             = k,
                         .xu            = curr->xu,
                         .p             = curr->p,
                         .norm_sq_p     = curr->pᵀp,
                         .x̂u            = curr->xû,
                         .φγ            = curr->fbe(),
                         .ψ             = curr->ψu,
                         .grad_ψ        = curr->grad_ψ,
                         .ψ_hat         = curr->ψû,
                         .q             = qxuₖ,
                         .gn            = did_gn,
                         .nJ            = did_gn ? J.sizes().sum() : nJ,
                         .lqr_min_rcond = lqr.min_rcond,
                         .L             = curr->L,
                         .γ             = curr->γ,
                         .τ             = τ,
                         .ε             = εₖ,
                         .problem       = problem,
                         .params        = params});
        }

        // Return solution -----------------------------------------------------

        auto time_elapsed = std::chrono::steady_clock::now() - start_time;
        auto stop_status =
            check_all_stop_conditions(time_elapsed, k, εₖ, no_progress);
        if (stop_status != SolverStatus::Busy) {
            assign_extract_u(curr->xû, u);
            s.iterations    = k;
            s.ε             = εₖ;
            s.elapsed_time  = duration_cast<nanoseconds>(time_elapsed);
            s.time_forward  = eval.time.forward;
            s.time_backward = eval.time.backward;
            s.time_backward_jacobians = eval.time.backward_jacobians;
            s.time_hessians           = eval.time.hessians;
            s.status                  = stop_status;
            s.final_γ                 = curr->γ;
            return s;
        }

        // Calculate Gauss-Newton step -----------------------------------------

        real_t τ_init = 1;
        did_gn        = do_gn_step;
        if (params.disable_acceleration) {
            τ_init = 0;
        } else if (do_gn_step) {
            // Forward simulation has been carried out on xuₖ and the Jacobians
            // are stored in 'eval'
            if (!curr->have_jacobians)
                throw std::logic_error("have_jacobians");

            // Make sure that Q and R are stored as well
            eval.hessians(curr->xu);

            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui = eval.uk(curr->xu, t)(i);
                // Gradient descent step.
                real_t gs = ui - curr->γ * curr->grad_ψ(t * nu + i);
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
            { // Find active indices J
                detail::Timed t{s.time_indices};
                J.update(is_constr_inactive);
            }

            { // LQR factor
                detail::Timed t{s.time_lqr_factor};
                lqr.factor_masked(Ak, Bk, Qk, Rk, qk, rk, uk_eq, Jk, Kk,
                                  params.lqr_factor_cholesky);
            }
            { // LQR solve
                detail::Timed t{s.time_lqr_solve};
                lqr.solve_masked(Ak, Bk, Jk, qxuₖ);
            }
        } else {
            if (!enable_lbfgs)
                throw std::logic_error("enable_lbfgs");

            // Find inactive indices J
            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui     = eval.uk(curr->xu, t)(i);
                real_t grad_i = curr->grad_ψ(t * nu + i);
                // Gradient descent step.
                real_t gs = ui - curr->γ * grad_i;
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub || active_lb) {
                    quₖ(t * nu + i) = curr->p(t * nu + i);
                    return false;
                } else { // Store inactive indices
                    quₖ(t * nu + i) = -grad_i;
                    return true;
                }
            };

            auto J_idx = J.indices();
            nJ         = 0;
            {
                detail::Timed t{s.time_lbfgs_indices};
                for (index_t t = 0; t < N; ++t)
                    for (index_t i = 0; i < nu; ++i)
                        if (is_constr_inactive(t, i))
                            J_idx(nJ++) = t * nu + i;
            }
            auto J_lbfgs = J_idx.topRows(nJ);

            // If all indices are inactive, we can use standard L-BFGS,
            // if there are active indices, we need the specialized version
            // that only applies L-BFGS to the inactive indices
            bool success = [&] {
                detail::Timed t{s.time_lbfgs_apply};
                return lbfgs.apply_masked(quₖ, curr->γ, J_lbfgs);
            }();
            // If L-BFGS application failed, qₖ(J) still contains
            // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
            if (not success) {
                τ_init = 0;
            } else {
                // for (index_t t = N; t-- > 0;)
                //     for (index_t i = nu; i-- > 0;)
                //         qxuₖ(t * (nx + nu) + nx + i) = qxuₖ(t * (0 + nu) + nx + i);
                assign_interleave_xu(quₖ, qxuₖ);
            }
        }

        // Make sure quasi-Newton step is valid
        if (not qxuₖ.allFinite()) {
            τ_init = 0;
            // Is there anything we can do?
            if (not did_gn)
                lbfgs.reset();
        }
        s.lbfgs_failures += (τ_init == 0 && k > 0);

        bool do_next_gn = params.gn_interval > 0 &&
                          ((k + 1) % params.gn_interval) == 0 &&
                          !params.disable_acceleration;
        do_gn_step = do_next_gn || (do_gn_step && params.gn_sticky);

        // Line search ---------------------------------------------------------

        next->γ = curr->γ;
        next->L = curr->L;
        τ       = τ_init;

        while (!stop_signal.stop_requested()) {

            // xₖ₊₁ = xₖ + (1-τ) pₖ + τ qₖ
            if (τ == 0) {
                next->xu = curr->xû;
                next->ψu = curr->ψû;
                // Calculate ∇ψ(xₖ₊₁)
                curr->have_jacobians = false;
                eval_backward(*next, do_gn_step);
            } else {
                if (τ == 1) {
                    next->xu = curr->xu + qxuₖ;
                } else {
                    do_gn_step = do_next_gn;
                    for (index_t t = 0; t < N; ++t)
                        eval.uk(next->xu, t) =
                            eval.uk(curr->xu, t) +
                            (1 - τ) * curr->p.segment(t * nu, nu) +
                            τ * eval.uk(qxuₖ, t);
                }
                // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
                eval_forward(*next); // Not necessary for DDP
                curr->have_jacobians = false;
                eval_backward(*next, do_gn_step);
            }

            // Calculate x̂ₖ₊₁, ψ(x̂ₖ₊₁)
            eval_prox(*next);
            eval_forward_hat(*next);

            // Quadratic upper bound
            if (next->L < params.L_max && qub_violated(*next)) {
                next->γ /= 2;
                next->L *= 2;
                τ = τ_init;
                continue;
            }

            // Line search condition
            if (τ > 0 && linesearch_violated(*curr, *next)) {
                τ /= 2;
                if (τ < params.τ_min)
                    τ = 0;
                continue;
            }

            // QUB and line search satisfied
            break;
        }

        // If τ < τ_min the line search failed and we accepted the prox step
        s.linesearch_failures += (τ == 0 && τ_init > 0);
        s.τ_1_accepted += τ == 1;
        s.count_τ += 1;
        s.sum_τ += τ;

        // Check if we made any progress
        if (no_progress > 0 || k % params.max_no_progress == 0)
            no_progress = curr->xu == next->xu ? no_progress + 1 : 0;

        // Update L-BFGS -------------------------------------------------------

        if (enable_lbfgs) {
            const bool force = true;
            assign_extract_u(next->xu, next->u);
            if (did_gn && params.reset_lbfgs_on_gn_step) {
                lbfgs.reset();
            } else {
                detail::Timed t{s.time_lbfgs_update};
                s.lbfgs_rejected += not lbfgs.update(
                    curr->u, next->u, curr->grad_ψ, next->grad_ψ,
                    LBFGS<config_t>::Sign::Positive, force);
            }
        }

        // Advance step --------------------------------------------------------
        std::swap(curr, next);
        ++k;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa