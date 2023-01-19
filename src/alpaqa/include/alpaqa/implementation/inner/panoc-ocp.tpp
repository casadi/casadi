#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/implementation/util/print.tpp>
#include <alpaqa/inner/directions/panoc-ocp/lqr.hpp>
#include <alpaqa/inner/directions/panoc-ocp/ocp-vars.hpp>
#include <alpaqa/inner/panoc-ocp.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/index-set.hpp>
#include <concepts>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

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
    rvec u,
    /// [in]    Lagrange multipliers @f$ y @f$
    crvec y,
    /// [in]    Penalty factors @f$ \mu @f$
    crvec μ) -> Stats {

    using std::chrono::nanoseconds;
    auto start_time = std::chrono::steady_clock::now();
    Stats s;

    const auto N    = problem.get_N();
    const auto nu   = problem.get_nu();
    const auto nx   = problem.get_nx();
    const auto nc   = problem.get_nc();
    const auto nc_N = problem.get_nc_N();
    const auto n    = nu * N;

    bool enable_lbfgs = params.gn_interval != 1;

    // Allocate storage --------------------------------------------------------

    // TODO: the L-BFGS objects and vectors allocate on each iteration of ALM,
    //       and there are more vectors than strictly necessary.

    OCPEvaluator<config_t> eval{problem};
    auto &vars = eval.vars;
    alpaqa::detail::IndexSet<config_t> J{N, nu};
    using LQRFactor = alpaqa::StatefulLQRFactor<config_t>;
    LQRFactor lqr{{.N = N, .nx = nx, .nu = nu}};
    LBFGS<config_t> lbfgs{params.lbfgs_params, enable_lbfgs ? n : 0};
    mat jacs = vars.create_AB();
    vec qr   = vars.create_qr();

    vec q(n); // Newton step, including states
    Box<config_t> U   = Box<config_t>::NaN(nu);
    Box<config_t> D   = Box<config_t>::NaN(nc);
    Box<config_t> D_N = Box<config_t>::NaN(nc_N);

    // Workspace storage
    vec work_2x(nx * 2), work_λ(nx), work_c(std::max(nc, nc_N));
    auto work_x  = work_2x.topRows(nx);
    auto work_cN = work_c.topRows(nc_N);
    auto work_ck = work_c.topRows(nc);
    vec work_R(problem.get_R_work_size()), work_S(problem.get_S_work_size());

    // ALM
    assert((nc == 0 && nc_N == 0) || μ.size() == N + 1);
    assert(y.size() == nc * N + nc_N);

    // Functions for accessing the LQR matrices and index sets
    auto ABk = [&](index_t i) -> crmat { return vars.ABk(jacs, i); };
    auto Qk  = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](rmat out) {
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Qk input");
#endif
                alpaqa::detail::Timed t{s.time_hessians};
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                auto xk  = vars.xk(storage, k);
                if (k < N)
                    problem.eval_add_Q(k, xuk, hk, out);
                else
                    problem.eval_add_Q_N(xk, hk, out);
                if (k < N) {
                    if (nc > 0) {
                        auto ck = vars.ck(storage, k);
                        auto yk = y.segment(k * nc, k < N ? nc : nc_N);
                        auto ζ  = ck + (real_t(1) / μ(k)) * yk;
                        for (index_t i = 0; i < nc; ++i)
                            work_ck(i) = μ(k) * (ζ(i) < D.lowerbound(i) ||
                                                 ζ(i) > D.upperbound(i));
                        problem.eval_add_gn_hess_constr(k, xk, work_ck, out);
                    }
                } else {
                    if (nc_N > 0) {
                        auto ck = vars.ck(storage, k);
                        auto yk = y.segment(k * nc, k < N ? nc : nc_N);
                        auto ζ  = ck + (real_t(1) / μ(k)) * yk;
                        for (index_t i = 0; i < nc_N; ++i)
                            work_cN(i) = μ(k) * (ζ(i) < D_N.lowerbound(i) ||
                                                 ζ(i) > D_N.upperbound(i));
                        problem.eval_add_gn_hess_constr_N(xk, work_cN, out);
                    }
                }
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Qk output");
#endif
            };
        };
    };
    auto Rk = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask, rmat out) {
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Rk input");
#endif
                alpaqa::detail::Timed t{s.time_hessians};
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_R_masked(k, xuk, hk, mask, out, work_R);
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Rk output");
#endif
            };
        };
    };
    auto Sk = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask, rmat out) {
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Sk input");
#endif
                alpaqa::detail::Timed t{s.time_hessians};
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_S_masked(k, xuk, hk, mask, out, work_S);
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Sk output");
#endif
            };
        };
    };
    auto Rk_prod = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask_J, crindexvec mask_K, crvec v,
                          rvec out) {
#ifndef NDEBUG
                {
                    ScopedMallocAllower ma;
                    vec vv = v(mask_K);
                    check_finiteness(vv.reshaped(), "Rk_prod input v");
                }
                check_finiteness(out.reshaped(), "Rk_prod input");
#endif
                alpaqa::detail::Timed t{s.time_hessians};
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_R_prod_masked(k, xuk, hk, mask_J, mask_K, v,
                                               out, work_R);
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Rk_prod output");
#endif
            };
        };
    };
    auto Sk_prod = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask_K, crvec v, rvec out) {
#ifndef NDEBUG
                {
                    ScopedMallocAllower ma;
                    vec vv = v(mask_K);
                    check_finiteness(vv.reshaped(), "Sk_prod input v");
                }
                check_finiteness(out.reshaped(), "Sk_prod input");
#endif
                alpaqa::detail::Timed t{s.time_hessians};
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_S_prod_masked(k, xuk, hk, mask_K, v, out,
                                               work_S);
#ifndef NDEBUG
                check_finiteness(out.reshaped(), "Sk_prod output");
#endif
            };
        };
    };
    auto mut_qrk = [&](index_t k) -> rvec { return vars.qrk(qr, k); };
    auto mut_q_N = [&]() -> rvec { return vars.qk(qr, N); };
    auto qk      = [&](index_t k) -> crvec { return vars.qk(qr, k); };
    auto rk      = [&](index_t k) -> crvec { return vars.rk(qr, k); };
    auto uk_eq   = [&](index_t k) -> crvec { return q.segment(k * nu, nu); };
    auto Jk      = [&](index_t k) -> crindexvec { return J.indices(k); };
    auto Kk      = [&](index_t k) -> crindexvec { return J.compl_indices(k); };

    // Iterates ----------------------------------------------------------------

    // Represents an iterate in the algorithm, keeping track of some
    // intermediate values and function evaluations.
    struct Iterate {
        vec xu;     //< Inputs u interleaved with states x
        vec xû;     //< Inputs u interleaved with states x after prox grad
        vec grad_ψ; //< Gradient of cost in u
        vec p;      //< Proximal gradient step in u
        vec u;      //< Inputs u (used for L-BFGS only)
        real_t ψu       = NaN<config_t>; //< Cost in u
        real_t ψû       = NaN<config_t>; //< Cost in û
        real_t γ        = NaN<config_t>; //< Step size γ
        real_t L        = NaN<config_t>; //< Lipschitz estimate L
        real_t pᵀp      = NaN<config_t>; //< Norm squared of p
        real_t grad_ψᵀp = NaN<config_t>; //< Dot product of gradient and p

        // @pre    @ref ψu, @ref pᵀp, @pre grad_ψᵀp
        // @return φγ
        real_t fbe() const { return ψu + pᵀp / (2 * γ) + grad_ψᵀp; }

        Iterate(const OCPVariables<config_t> &vars, bool enable_lbfgs)
            : xu{vars.create()}, xû{vars.create()}, grad_ψ{vars.N * vars.nu()},
              p{vars.N * vars.nu()}, u{enable_lbfgs ? vars.N * vars.nu() : 0} {}
    } iterates[2]{
        {vars, enable_lbfgs},
        {vars, enable_lbfgs},
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

    auto eval_prox_impl = [&](real_t γ, crvec xu, crvec grad_ψ, rvec x̂u,
                              rvec p) {
        alpaqa::detail::Timed t{s.time_prox};
        real_t pᵀp      = 0;
        real_t grad_ψᵀp = 0;
        for (index_t t = 0; t < N; ++t) {
            auto &&grad_ψ_t = grad_ψ.segment(t * nu, nu);
            auto &&p_t      = p.segment(t * nu, nu);
            eval_proj_grad_step_box(γ, vars.uk(xu, t), grad_ψ_t,
                                    /* in ⟹ out */ vars.uk(x̂u, t), p_t);
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

    auto assign_interleave_xu = [&vars](crvec u, rvec xu) {
        detail::assign_interleave_xu(vars, u, xu);
    };
    auto assign_extract_u = [&vars](crvec xu, rvec u) {
        detail::assign_extract_u(vars, xu, u);
    };

    // @pre    @ref Iterate::γ, @ref Iterate::xu, @ref Iterate::grad_ψ
    // @post   @ref Iterate::xû, @ref Iterate::p, @ref Iterate::pᵀp,
    //         @ref Iterate::grad_ψᵀp
    auto eval_prox = [&](Iterate &i) {
        std::tie(i.pᵀp, i.grad_ψᵀp) =
            eval_prox_impl(i.γ, i.xu, i.grad_ψ, i.xû, i.p);
    };

    // @pre    @ref Iterate::xu
    // @post   @ref Iterate::ψu
    auto eval_forward = [&](Iterate &i) {
        alpaqa::detail::Timed t{s.time_forward};
        i.ψu = eval.forward(i.xu, D, D_N, μ, y);
    };
    // @pre    @ref Iterate::xû
    // @post   @ref Iterate::ψû
    auto eval_forward_hat = [&](Iterate &i) {
        alpaqa::detail::Timed t{s.time_forward};
        i.ψû = eval.forward(i.xû, D, D_N, μ, y);
    };

    // @pre    @ref Iterate::xu
    // @post   @ref Iterate::grad_ψ, q, q_N
    auto eval_backward = [&](Iterate &i) {
        alpaqa::detail::Timed t{s.time_backward};
        eval.backward(i.xu, i.grad_ψ, work_λ, work_x, work_c, mut_qrk, mut_q_N,
                      D, D_N, μ, y);
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
        [&](
            /// Iterate, updates xu, ψ, grad_ψ, have_jacobians, L
            Iterate *it,
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
            eval_backward(*it);
            // Select a small step h for finite differences
            auto h        = it->grad_ψ.unaryExpr([&](real_t g) {
                return g > 0 ? std::max(g * ε, δ) : std::min(g * ε, -δ);
            });
            real_t norm_h = h.norm();
            // work_xu = xu - h
            for (index_t t = 0; t < N; ++t)
                vars.uk(work_xu, t) =
                    vars.uk(it->xu, t) - h.segment(t * nu, nu);

            { // Calculate ψ(x₀ - h)
                alpaqa::detail::Timed t{s.time_forward};
                eval.forward_simulate(work_xu); // needed for backwards sweep
            }
            { // Calculate ∇ψ(x₀ + h)
                alpaqa::detail::Timed t{s.time_backward};
                eval.backward(work_xu, work_grad_ψ, work_λ, work_x, work_c,
                              mut_qrk, mut_q_N, D, D_N, μ, y);
            }
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
    [[maybe_unused]] auto print_progress =
        [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ, real_t pₖᵀpₖ,
            crvec qₖ, real_t γₖ, real_t τₖ, real_t εₖ, bool did_gn, length_t nJ,
            real_t min_rcond) {
            *os << "[PANOC] " << std::setw(6) << k
                << ": φγ = " << print_real(φₖ) << ", ψ = " << print_real(ψₖ)
                << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())
                << ", ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ))
                << ", γ = " << print_real(γₖ) << ", εₖ = " << print_real(εₖ);
            if (k > 0)
                *os << ", τ = " << print_real3(τₖ)
                    << ", ‖q‖ = " << print_real(qₖ.norm())
                    << ", #J =" << std::setw(5) << nJ
                    << ", rcond = " << print_real3(min_rcond) << ", "
                    << (did_gn ? "GN" : "L-BFGS");
            *os << std::endl; // Flush for Python buffering
        };

    auto print_progress_1 = [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ,
                                real_t pₖᵀpₖ, real_t γₖ, real_t εₖ) {
        if (k == 0)
            *os << "┌─[PANOC]\n";
        else
            *os << "├─ " << std::setw(6) << k << '\n';
        *os << "│   φγ = " << print_real(φₖ)               //
            << ",    ψ = " << print_real(ψₖ)               //
            << ", ‖∇ψ‖ = " << print_real(grad_ψₖ.norm())   //
            << ",  ‖p‖ = " << print_real(std::sqrt(pₖᵀpₖ)) //
            << ",    γ = " << print_real(γₖ)               //
            << ",    ε = " << print_real(εₖ) << '\n';
    };
    auto print_progress_2 = [&](crvec qₖ, real_t τₖ, bool did_gn, length_t nJ,
                                real_t min_rcond) {
        *os << "│  ‖q‖ = " << print_real(qₖ.norm())                       //
            << ",   #J = " << std::setw(7 + params.print_precision) << nJ //
            << ", cond = " << print_real3(real_t(1) / min_rcond)          //
            << ",    τ = " << print_real3(τₖ)                             //
            << ",    " << (did_gn ? "GN" : "L-BFGS")                      //
            << std::endl; // Flush for Python buffering
    };
    auto print_progress_n = [&](SolverStatus status) {
        *os << "└─ " << status << " ──"
            << std::endl; // Flush for Python buffering
    };

    // Initialize inputs and initial state (do not simulate states yet) --------

    assign_interleave_xu(u, curr->xu);           // initial guess
    problem.get_x_init(curr->xu.topRows(nx));    // initial state
    curr->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xu.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
    if (enable_lbfgs)
        curr->u = u;

    problem.get_U(U);     // input box constraints
    problem.get_D(D);     // general constraints
    problem.get_D_N(D_N); // general terminal constraints

    bool do_gn_step = params.gn_interval > 0 and !params.disable_acceleration;
    bool did_gn     = false;

    // Make sure that we don't allocate any memory in the inner loop
    ScopedMallocBlocker mb;

    // Estimate Lipschitz constant ---------------------------------------------

    // Finite difference approximation of ∇²ψ in starting point
    if (params.Lipschitz.L_0 <= 0) {
        initial_lipschitz_estimate(curr, params.Lipschitz.ε, params.Lipschitz.δ,
                                   params.L_min, params.L_max, next->xu,
                                   next->grad_ψ);
    }
    // Initial Lipschitz constant provided by the user
    else {
        curr->L = params.Lipschitz.L_0;
        // Calculate ψ(x₀), ∇ψ(x₀)
        eval_forward(*curr);
        eval_backward(*curr);
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
        bool do_print =
            params.print_interval != 0 && k % params.print_interval == 0;
        if (do_print)
            print_progress_1(k, curr->fbe(), curr->ψu, curr->grad_ψ, curr->pᵀp,
                             curr->γ, εₖ);
        if (progress_cb) {
            ScopedMallocAllower ma;
            alpaqa::detail::Timed t{s.time_progress_callback};
            progress_cb({.k             = k,
                         .xu            = curr->xu,
                         .p             = curr->p,
                         .norm_sq_p     = curr->pᵀp,
                         .x̂u            = curr->xû,
                         .φγ            = curr->fbe(),
                         .ψ             = curr->ψu,
                         .grad_ψ        = curr->grad_ψ,
                         .ψ_hat         = curr->ψû,
                         .q             = q,
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
            bool do_final_print = params.print_interval != 0;
            if (!do_print && do_final_print)
                print_progress_1(k, curr->fbe(), curr->ψu, curr->grad_ψ,
                                 curr->pᵀp, curr->γ, εₖ);
            if (do_print || do_final_print)
                print_progress_n(stop_status);
            assign_extract_u(curr->xû, u);
            s.iterations   = k;
            s.ε            = εₖ;
            s.elapsed_time = duration_cast<nanoseconds>(time_elapsed);
            s.time_lqr_factor -= s.time_hessians;
            s.status  = stop_status;
            s.final_γ = curr->γ;
            return s;
        }

        // Calculate Gauss-Newton step -----------------------------------------

        real_t τ_init = 1;
        did_gn        = do_gn_step;
        if (params.disable_acceleration) {
            τ_init = 0;
        } else if (do_gn_step) {
            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui = vars.uk(curr->xu, t)(i);
                // Gradient descent step.
                real_t gs = ui - curr->γ * curr->grad_ψ(t * nu + i);
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub) {
                    q(nu * t + i) = U.upperbound(i) - ui;
                    return false;
                } else if (active_lb) {
                    q(nu * t + i) = U.lowerbound(i) - ui;
                    return false;
                } else { // Store inactive indices
                    return true;
                }
            };
            { // Find active indices J
                alpaqa::detail::Timed t{s.time_indices};
                J.update(is_constr_inactive);
            }
            { // evaluate the Jacobians
                alpaqa::detail::Timed t{s.time_jacobians};
                for (index_t t = 0; t < N; ++t)
                    problem.eval_jac_f(t, vars.xk(curr->xu, t),
                                       vars.uk(curr->xu, t), vars.ABk(jacs, t));
            }
            { // LQR factor
                alpaqa::detail::Timed t{s.time_lqr_factor};
                lqr.factor_masked(ABk, Qk(curr->xu), Rk(curr->xu), Sk(curr->xu),
                                  Rk_prod(curr->xu), Sk_prod(curr->xu), qk, rk,
                                  uk_eq, Jk, Kk, params.lqr_factor_cholesky);
            }
            { // LQR solve
                alpaqa::detail::Timed t{s.time_lqr_solve};
                lqr.solve_masked(ABk, Jk, q, work_2x);
            }
        } else {
            if (!enable_lbfgs)
                throw std::logic_error("enable_lbfgs");

            // Find inactive indices J
            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui     = vars.uk(curr->xu, t)(i);
                real_t grad_i = curr->grad_ψ(t * nu + i);
                // Gradient descent step.
                real_t gs = ui - curr->γ * grad_i;
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub || active_lb) {
                    q(t * nu + i) = curr->p(t * nu + i);
                    return false;
                } else { // Store inactive indices
                    q(t * nu + i) = -grad_i;
                    return true;
                }
            };

            auto J_idx = J.indices();
            nJ         = 0;
            {
                alpaqa::detail::Timed t{s.time_lbfgs_indices};
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
                alpaqa::detail::Timed t{s.time_lbfgs_apply};
                return lbfgs.apply_masked(q, curr->γ, J_lbfgs);
            }();
            // If L-BFGS application failed, qₖ(J) still contains
            // -∇ψ(x)(J) - HqK(J) or -∇ψ(x)(J), which is not a valid step.
            if (not success)
                τ_init = 0;
        }

        // Make sure quasi-Newton step is valid
        if (not q.allFinite()) {
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

        next->γ       = curr->γ;
        next->L       = curr->L;
        τ             = τ_init;
        real_t τ_prev = -1;

        // xₖ₊₁ = xₖ + pₖ
        auto take_safe_step = [&] {
            next->xu = curr->xû;
            next->ψu = curr->ψû;
            // Calculate ∇ψ(xₖ₊₁)
            eval_backward(*next);
        };

        // xₖ₊₁ = xₖ + (1-τ) pₖ + τ qₖ
        auto take_accelerated_step = [&](real_t τ) {
            if (τ == 1) {
                for (index_t t = 0; t < N; ++t)
                    vars.uk(next->xu, t) =
                        vars.uk(curr->xu, t) + q.segment(t * nu, nu);
            } else {
                do_gn_step = do_next_gn;
                for (index_t t = 0; t < N; ++t)
                    vars.uk(next->xu, t) =
                        vars.uk(curr->xu, t) +
                        (1 - τ) * curr->p.segment(t * nu, nu) +
                        τ * q.segment(t * nu, nu);
            }
            // Calculate ψ(xₖ₊₁), ∇ψ(xₖ₊₁)
            eval_forward(*next); // Not necessary for DDP
            eval_backward(*next);
        };

        // Backtracking line search loop
        while (!stop_signal.stop_requested()) {

            // Recompute step only if τ changed
            if (τ != τ_prev) {
                τ != 0 ? take_accelerated_step(τ) : take_safe_step();
                τ_prev = τ;
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
                alpaqa::detail::Timed t{s.time_lbfgs_update};
                s.lbfgs_rejected += not lbfgs.update(
                    curr->u, next->u, curr->grad_ψ, next->grad_ψ,
                    LBFGS<config_t>::Sign::Positive, force);
            }
        }

        // Print ---------------------------------------------------------------
        if (do_print && k != 0)
            print_progress_2(q, τ, did_gn, did_gn ? J.sizes().sum() : nJ,
                             lqr.min_rcond);

        // Advance step --------------------------------------------------------
        std::swap(curr, next);
        ++k;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa