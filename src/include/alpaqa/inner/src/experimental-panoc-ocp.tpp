#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc-ocp/experimental-lqr.hpp>
#include <alpaqa/inner/experimental-panoc-ocp.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/index-set.hpp>
#include <alpaqa/util/src/print.tpp>
#include <concepts>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace alpaqa::experimental {

template <class V, class Conf>
concept VectorRefLike =
    std::convertible_to<V, rvec<Conf>> || std::convertible_to<V, crvec<Conf>>;

template <Config Conf, VectorRefLike<Conf> V>
constexpr auto const_or_mut_rvec(V &&v) {
    if constexpr (Eigen::internal::is_lvalue<std::remove_reference_t<V>>::value)
        return rvec<Conf>{v};
    else
        return crvec<Conf>{v};
}

template <Config Conf>
struct OCPVariables {
    USING_ALPAQA_CONFIG(Conf);

    OCPVariables(
        /// nx, nu, nh, nc
        const std::array<index_t, 4> &sizes,
        /// nx, nh, nc
        const std::array<index_t, 3> &sizes_N,
        /// Horizon length
        length_t N)
        : N{N} {
        std::partial_sum(sizes.begin(), sizes.end(), indices.begin());
        std::partial_sum(sizes_N.begin(), sizes_N.end(), indices_N.begin());
    }
    OCPVariables(const TypeErasedOCProblem<config_t> &prob)
        : OCPVariables{
              {prob.get_nx(), prob.get_nu(), prob.get_nh(), prob.get_nc()},
              {prob.get_nx(), prob.get_nh_N(), prob.get_nc_N()},
              prob.get_N(),
          } {}

    enum Indices {
        i_u   = 0,
        i_h   = 1,
        i_c   = 2,
        i_h_N = 0,
        i_c_N = 1,
    };
    length_t N;
    std::array<index_t, 4> indices;
    std::array<index_t, 3> indices_N;
    length_t size(size_t i) const { return indices[i + 1] - indices[i]; }
    length_t size_N(size_t i) const { return indices_N[i + 1] - indices_N[i]; }
    length_t nx() const { return indices[0]; }
    length_t nu() const { return size(i_u); }
    length_t nxu() const { return nx() + nu(); }
    length_t nh() const { return size(i_h); }
    length_t nc() const { return size(i_c); }
    length_t nx_N() const { return indices_N[0]; }
    length_t nh_N() const { return size_N(i_h_N); }
    length_t nc_N() const { return size_N(i_c_N); }

    vec create() const { return vec(N * indices.back() + indices_N.back()); }
    auto xk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(v.segment(t * indices.back(), nx()));
    }
    auto xuk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(
            v.segment(t * indices.back(), nxu()));
    }
    auto uk(VectorRefLike<config_t> auto &&v, index_t t) const {
        assert(t < N);
        return const_or_mut_rvec<config_t>(
            v.segment(t * indices.back() + indices[0], nu()));
    }
    auto hk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(v.segment(
            t * indices.back() + (t < N ? indices[i_h] : indices_N[i_h_N]),
            (t < N ? nh() : nh_N())));
    }
    auto ck(VectorRefLike<config_t> auto &&v, index_t t) const {
        return v.segment(t * indices.back() +
                             (t < N ? indices[i_c] : indices_N[i_c_N]),
                         (t < N ? nc() : nc_N()));
    }

    vec create_qr() const { return vec(N * nxu()); }
    auto qk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(v.segment(t * nxu(), nx()));
    }
    auto rk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(v.segment(t * nxu() + nx(), nu()));
    }
    auto qrk(VectorRefLike<config_t> auto &&v, index_t t) const {
        return const_or_mut_rvec<config_t>(v.segment(t * nxu(), nxu()));
    }

    mat create_AB() const { return mat(nx(), nxu() * N); }
    rmat ABk(rmat AB, index_t t) const {
        return AB.middleCols(t * nxu(), nxu());
    }
    rmat Ak(rmat AB, index_t t) const { return AB.middleCols(t * nxu(), nx()); }
    rmat Bk(rmat AB, index_t t) const {
        return AB.middleCols(t * nxu() + nx(), nu());
    }
};

template <Config Conf>
struct OCPEvaluator {
    USING_ALPAQA_CONFIG(Conf);
    using OCPVars = OCPVariables<config_t>;
    using Problem = TypeErasedOCProblem<config_t>;
    using Box     = alpaqa::Box<config_t>;
    const Problem *problem;
    OCPVars vars;
    OCPEvaluator(const Problem &problem) : problem{&problem}, vars{problem} {}

    length_t N() const { return vars.N; }

    /// @pre x0 and u initialized
    /// @post x, h and c updated
    /// @return @f$ V(u) =
    ///             \sum_{k=0}^{N-1} \ell(h_k(x_k, u_k)) + V_f(h_N(x_N)) @f$
    real_t forward(rvec storage, const Box &D, const Box &D_N, crvec μ,
                   crvec y) const {
        real_t V = 0;
        auto N   = this->N();
        auto nc  = vars.nc();
        for (index_t t = 0; t < N; ++t) {
            auto xk = vars.xk(storage, t);
            auto uk = vars.uk(storage, t);
            auto hk = vars.hk(storage, t);
            auto ck = vars.ck(storage, t);
            problem->eval_h(t, xk, uk, hk);
            V += problem->eval_l(t, hk);
            if (nc > 0) {
                problem->eval_constr(t, xk, ck);
                auto yk = y.segment(t * nc, nc);
                auto ζ  = ck + (real_t(1) / μ(t)) * yk;
                V += real_t(0.5) * μ(t) *
                     projecting_difference(ζ, D).squaredNorm();
            }
            problem->eval_f(t, xk, uk, vars.xk(storage, t + 1));
        }
        auto xN = vars.xk(storage, N);
        auto hN = vars.hk(storage, N);
        auto cN = vars.ck(storage, N);
        problem->eval_h_N(xN, hN);
        V += problem->eval_l_N(hN);
        if (nc > 0) {
            problem->eval_constr_N(xN, cN);
            auto yN = y.segment(N * nc, nc);
            auto ζ  = cN + (real_t(1) / μ(N)) * yN;
            V += real_t(0.5) * μ(N) *
                 projecting_difference(ζ, D_N).squaredNorm();
        }
        return V;
    }

    /// @pre x0 and u initialized
    /// @post x updated
    void forward_simulate(rvec storage) const {
        for (index_t t = 0; t < N(); ++t) {
            auto xk = vars.xk(storage, t);
            auto uk = vars.uk(storage, t);
            problem->eval_f(t, xk, uk, vars.xk(storage, t + 1));
        }
    }

    /// @pre x0 and u initialized
    void forward_simulate(crvec u, rvec x) const {
        assert(u.size() == vars.N * vars.nu());
        assert(x.size() == vars.nx());
        for (index_t t = 0; t < N(); ++t) {
            auto uk = u.segment(t * vars.nu(), vars.nu());
            problem->eval_f(t, x, uk, x);
        }
    }

    /// @pre x, u, h and c initialized (i.e. forward was called)
    void backward(rvec storage, rvec g, rvec λ, rvec w, rvec v, const auto &qr,
                  const Box &D, const Box &D_N, crvec μ, crvec y) const {
        auto N    = this->N();
        auto nc   = vars.nc();
        auto nc_N = vars.nc_N();
        auto nu   = vars.nu();
        auto nx   = vars.nx();
        assert((nc <= 0 && nc_N <= 0) || w.size() == nx);
        assert((nc <= 0 && nc_N <= 0) || v.size() == nc);
        auto qN = qr(N).topRows(nx);
        auto xN = vars.xk(storage, N);
        auto hN = vars.hk(storage, N);
        auto vN = v.topRows(nc_N);
        auto vk = v.topRows(nc);
        // λ ← ∇h(x)·∇l(h)
        problem->eval_q_N(xN, hN, λ);
        // λ ← ∇h(x)·∇l(h) + ∇c(x)·μ·(c(x) + μ⁻¹y - Π(c(x) + μ⁻¹y; D))
        if (nc > 0) {
            auto cN = vars.ck(storage, N);
            auto yN = y.segment(N * nc, nc_N);
            auto ζ  = cN + (real_t(1) / μ(N)) * yN;
            vN      = μ(N) * projecting_difference(ζ, D_N);
            problem->eval_grad_constr_prod_N(xN, vN, w);
            λ += w;
        }
        qN = λ;
        for (index_t t = N; t-- > 0;) {
            auto gt    = g.segment(t * nu, nu);
            auto hk    = vars.hk(storage, t);
            auto xuk   = vars.xuk(storage, t);
            auto xk    = vars.xk(storage, t);
            auto uk    = vars.uk(storage, t);
            auto &&qrk = qr(t);
            auto &&qk  = qrk.topRows(nx);
            auto &&rk  = qrk.bottomRows(nu);
            // /q\ ← /Aᵀ\ λ
            // \r/   \Bᵀ/ λ
            problem->eval_grad_f_prod(t, xk, uk, λ, qrk);
            // λ ← Aᵀλ, ∇ψ ← Bᵀλ
            λ  = qk;
            gt = rk;
            // /q\ ← ∇h(x,u)·∇l(h)
            // \r/
            problem->eval_qr(t, xuk, hk, qrk);
            // q ← ∇h(x)·∇l(h) + ∇c(x)·μ·(c(x) + μ⁻¹y - Π(c(x) + μ⁻¹y; D))
            if (nc > 0) {
                auto ck = vars.ck(storage, t);
                auto yk = y.segment(t * nc, nc);
                auto ζ  = ck + (real_t(1) / μ(t)) * yk;
                vk      = μ(t) * projecting_difference(ζ, D);
                problem->eval_grad_constr_prod(t, xk, vk, w);
                qk += w;
            }
            // λ ← q + Aᵀλ, ∇ψ ← r + Bᵀλ
            λ += qk;
            gt += rk;
        }
    }
};

namespace detail {

template <Config Conf>
void assign_interleave_xu(const OCPVariables<Conf> &dim, crvec<Conf> u,
                          rvec<Conf> storage) {
    for (index_t<Conf> t = 0; t < dim.N; ++t)
        dim.uk(storage, t) = u.segment(t * dim.nu(), dim.nu());
}
template <Config Conf>
void assign_interleave_xu(const OCPVariables<Conf> &dim, crvec<Conf> x,
                          crvec<Conf> u, rvec<Conf> storage) {
    for (index_t<Conf> t = 0; t < dim.N; ++t) {
        dim.xk(storage, t) = x.segment(t * dim.nx(), dim.nx());
        dim.uk(storage, t) = u.segment(t * dim.nu(), dim.nu());
    }
    dim.xk(storage, dim.N) = x.segment(dim.N * dim.nx(), dim.nx());
}
template <Config Conf>
void assign_extract_u(const OCPVariables<Conf> &dim, crvec<Conf> storage,
                      rvec<Conf> u) {
    for (index_t<Conf> t = 0; t < dim.N; ++t)
        u.segment(t * dim.nu(), dim.nu()) = dim.uk(storage, t);
}
template <Config Conf>
void assign_extract_x(const OCPVariables<Conf> &dim, crvec<Conf> storage,
                      rvec<Conf> x) {
    for (index_t<Conf> t = 0; t < dim.N + 1; ++t)
        x.segment(t * dim.nx(), dim.nx()) =
            storage.segment(t * (dim.nx() + dim.nu()), dim.nx());
}

template <Config Conf>
vec<Conf> extract_u(const TypeErasedOCProblem<Conf> &problem, crvec<Conf> xu) {
    OCPVariables<Conf> dim{problem};
    vec<Conf> u(dim.N * dim.nu());
    assign_extract_u(dim, xu, u);
    return u;
}
template <Config Conf>
vec<Conf> extract_x(const TypeErasedOCProblem<Conf> &problem, crvec<Conf> xu) {
    OCPVariables<Conf> dim{problem};
    vec<Conf> x((dim.N + 1) * dim.nx());
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
    using LQRFactor = alpaqa::experimental::StatefulLQRFactor<config_t>;
    LQRFactor lqr{{.N = N, .nx = nx, .nu = nu}};
    LBFGSParams<config_t> lbfgs_param{.memory = N}; // TODO: make configurable
    LBFGS<config_t> lbfgs{lbfgs_param, enable_lbfgs ? n : 0};
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
    vec μ = vec::Zero(N + 1);
    vec y = vec::Zero(nc * N + nc_N);

    // Functions for accessing the LQR matrices and index sets
    auto ABk = [&](index_t i) -> crmat { return vars.ABk(jacs, i); };
    auto Qk  = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](rmat out) {
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                auto xk  = vars.xk(storage, k);
                if (k < N)
                    problem.eval_add_Q(k, xuk, hk, out);
                else
                    problem.eval_add_Q_N(xk, hk, out);
                if (nc > 0) {
                    auto ck = vars.ck(storage, k);
                    auto yk = y.segment(k * nc, k < N ? nc : nc_N);
                    auto ζ  = ck + (real_t(1) / μ(k)) * yk;
                    if (k < N) {
                        for (index_t i = 0; i < nc; ++i)
                            work_ck(i) = μ(k) * (ζ(i) < D.lowerbound(i) ||
                                                 ζ(i) > D.upperbound(i));
                        problem.eval_add_gn_hess_constr(k, xk, work_ck, out);
                    } else {
                        for (index_t i = 0; i < nc_N; ++i)
                            work_cN(i) = μ(k) * (ζ(i) < D_N.lowerbound(i) ||
                                                 ζ(i) > D_N.upperbound(i));
                        problem.eval_add_gn_hess_constr_N(xk, work_cN, out);
                    }
                }
            };
        };
    };
    auto Rk = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask, rmat out) {
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_R_masked(k, xuk, hk, mask, out, work_R);
            };
        };
    };
    auto Sk = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask, rmat out) {
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_S_masked(k, xuk, hk, mask, out, work_S);
            };
        };
    };
    auto Rk_prod = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask_J, crindexvec mask_K, crvec v,
                          rvec out) {
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_R_prod_masked(k, xuk, hk, mask_J, mask_K, v,
                                               out, work_R);
            };
        };
    };
    auto Sk_prod = [&](rvec storage) {
        return [&, storage](index_t k) {
            return [&, k](crindexvec mask_K, crvec v, rvec out) {
                auto hk  = vars.hk(storage, k);
                auto xuk = vars.xuk(storage, k);
                problem.eval_add_S_prod_masked(k, xuk, hk, mask_K, v, out,
                                               work_S);
            };
        };
    };
    auto mut_qrk = [&](index_t k) -> rvec { return vars.qrk(qr, k); };
    auto qk      = [&](index_t k) -> crvec { return vars.qk(qr, k); };
    auto rk      = [&](index_t k) -> crvec { return vars.rk(qr, k); };
    auto uk_eq   = [&](index_t k) -> crvec { return q.segment(k * nu, nu); };
    auto Jk      = [&](index_t k) -> crindexvec { return J.indices(k); };
    auto Kk      = [&](index_t k) -> crindexvec { return J.compl_indices(k); };

    // Iterates ----------------------------------------------------------------

    struct Iterate {
        vec xu;     ///< Inputs u interleaved with states x
        vec xû;     ///< Inputs u interleaved with states x after prox grad
        vec grad_ψ; ///< Gradient of cost in u
        vec p;      ///< Proximal gradient step in u
        vec u;      ///< Inputs u (used for L-BFGS only)
        real_t ψu       = NaN<config_t>; ///< Cost in u
        real_t ψû       = NaN<config_t>; ///< Cost in û
        real_t γ        = NaN<config_t>; ///< Step size γ
        real_t L        = NaN<config_t>; ///< Lipschitz estimate L
        real_t pᵀp      = NaN<config_t>; ///< Norm squared of p
        real_t grad_ψᵀp = NaN<config_t>; ///< Dot product of gradient and p

        /// @pre    @ref ψu, @ref pᵀp, @pre grad_ψᵀp
        /// @return φγ
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
        alpaqa::detail::Timed{s.time_prox};
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

    /// @pre    @ref Iterate::γ, @ref Iterate::xu, @ref Iterate::grad_ψ
    /// @post   @ref Iterate::xû, @ref Iterate::p, @ref Iterate::pᵀp,
    ///         @ref Iterate::grad_ψᵀp
    auto eval_prox = [&](Iterate &i) {
        std::tie(i.pᵀp, i.grad_ψᵀp) =
            eval_prox_impl(i.γ, i.xu, i.grad_ψ, i.xû, i.p);
    };

    /// @pre    @ref Iterate::xu
    /// @post   @ref Iterate::ψu
    auto eval_forward = [&](Iterate &i) {
        alpaqa::detail::Timed{s.time_forward};
        i.ψu = eval.forward(i.xu, D, D_N, μ, y);
    };
    /// @pre    @ref Iterate::xû
    /// @post   @ref Iterate::ψû
    auto eval_forward_hat = [&](Iterate &i) {
        alpaqa::detail::Timed{s.time_forward};
        i.ψû = eval.forward(i.xû, D, D_N, μ, y);
    };

    /// @pre    @ref Iterate::xu
    /// @post   @ref Iterate::grad_ψ, @ref Iterate::have_jacobians
    auto eval_backward = [&](Iterate &i) {
        alpaqa::detail::Timed{s.time_backward};
        eval.backward(i.xu, i.grad_ψ, work_λ, work_x, work_c, mut_qrk, D, D_N,
                      μ, y);
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
                alpaqa::detail::Timed{s.time_forward};
                eval.forward_simulate(work_xu); // needed for backwards sweep
            }
            { // Calculate ∇ψ(x₀ + h)
                alpaqa::detail::Timed{s.time_backward};

                eval.backward(work_xu, work_grad_ψ, work_λ, work_x, work_c,
                              mut_qrk, D, D_N, μ, y);
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
    auto print_progress = [&](unsigned k, real_t φₖ, real_t ψₖ, crvec grad_ψₖ,
                              real_t pₖᵀpₖ, crvec qₖ, real_t γₖ, real_t τₖ,
                              real_t εₖ, bool did_gn, length_t nJ,
                              real_t min_rcond) {
        *os << "[PANOC] " << std::setw(6) << k << ": φγ = " << print_real(φₖ)
            << ", ψ = " << print_real(ψₖ)
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

    // Initialize inputs and initial state (do not simulate states yet) --------

    assign_interleave_xu(u, curr->xu);           // initial guess
    problem.get_x_init(curr->xu.topRows(nx));    // initial state
    curr->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xu.topRows(nx) = curr->xu.topRows(nx); // initial state
    next->xû.topRows(nx) = curr->xu.topRows(nx); // initial state
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

        if (params.print_interval != 0 && k % params.print_interval == 0)
            print_progress(k, curr->fbe(), curr->ψu, curr->grad_ψ, curr->pᵀp, q,
                           curr->γ, τ, εₖ, did_gn,
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
                alpaqa::detail::Timed{s.time_backward_jacobians};
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

        while (!stop_signal.stop_requested()) {

            // Recompute step only if τ changed
            if (τ != τ_prev) {
                // xₖ₊₁ = xₖ + (1-τ) pₖ + τ qₖ
                if (τ == 0) {
                    next->xu = curr->xû;
                    next->ψu = curr->ψû;
                    // Calculate ∇ψ(xₖ₊₁)
                    eval_backward(*next);
                } else {
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
                }
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

        // Advance step --------------------------------------------------------
        std::swap(curr, next);
        ++k;
    }
    throw std::logic_error("[PANOC] loop error");
}

} // namespace alpaqa::experimental