#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/dynamics.hpp>
#include <alpaqa/problem/problem-counters.hpp>

#include <chrono>

namespace alpaqa {

template <Config Conf>
struct DynamicsEvaluator {
    USING_ALPAQA_CONFIG(Conf);

    DynamicsEvaluator(const TypeErasedControlProblem<config_t> &problem)
        : problem{problem} {
        N         = problem.get_N();
        nu        = problem.get_nu();
        nx        = problem.get_nx();
        structure = problem.get_l_structure();
        AB.resize(nx, N * (nx + nu));
        qr.resize((nx + nu) * N + nx);
        switch (structure) {
            case CostStructure::General:
                QRS.resize(nx + nu, (nx + nu) * N);
                Q_N.resize(nx, nx);
                break;
            case CostStructure::DiagonalHessian:
                QRS.resize(nx + nu, N);
                Q_N.resize(nx, 1);
                break;
            case CostStructure::Quadratic:
                QRS.resize(nx + nu, nx + nu);
                Q_N.resize(nx, nx);
                break;
            case CostStructure::DiagonalQuadratic:
                QRS.resize(nx + nu, 1);
                Q_N.resize(nx, 1);
                break;
            default: throw std::logic_error("CostStructure");
        }
    }

    const TypeErasedControlProblem<config_t> &problem;
    CostStructure structure;
    length_t N, nx, nu;
    mat QRS;
    mat Q_N;
    mat AB;
    vec qr;

    // clang-format off
    rvec xuk(rvec xu, index_t k) const { return xu.segment(k * (nx + nu), nx + nu); }
    crvec xuk(crvec xu, index_t k) const { return xu.segment(k * (nx + nu), nx + nu); }
    rvec xuk(vec &xu, index_t k) const { return xu.segment(k * (nx + nu), nx + nu); }
    crvec xuk(const vec &xu, index_t k) const { return xu.segment(k * (nx + nu), nx + nu); }
    rvec xk(rvec xu, index_t k) const { return xu.segment(k * (nx + nu), nx); }
    crvec xk(crvec xu, index_t k) const { return xu.segment(k * (nx + nu), nx); }
    rvec xk(vec &xu, index_t k) const { return xu.segment(k * (nx + nu), nx); }
    crvec xk(const vec &xu, index_t k) const { return xu.segment(k * (nx + nu), nx); }
    rvec uk(rvec xu, index_t k) const { return xu.segment(k * (nx + nu) + nx, nu); }
    crvec uk(crvec xu, index_t k) const { return xu.segment(k * (nx + nu) + nx, nu); }
    rvec uk(vec &xu, index_t k) const { return xu.segment(k * (nx + nu) + nx, nu); }
    crvec uk(const vec &xu, index_t k) const { return xu.segment(k * (nx + nu) + nx, nu); }
    // clang-format on

    rvec qrk(index_t k) { return xuk(qr, k); }
    crvec qrk(index_t k) const { return xuk(qr, k); }
    rvec qk(index_t k) { return xk(qr, k); }
    crvec qk(index_t k) const { return xk(qr, k); }
    rvec rk(index_t k) { return uk(qr, k); }
    crvec rk(index_t k) const { return uk(qr, k); }

    rmat ABk(index_t k) { return AB.middleCols(k * (nx + nu), nx + nu); }
    crmat ABk(index_t k) const { return AB.middleCols(k * (nx + nu), nx + nu); }
    rmat Ak(index_t k) { return AB.middleCols(k * (nx + nu), nx); }
    crmat Ak(index_t k) const { return AB.middleCols(k * (nx + nu), nx); }
    rmat Bk(index_t k) { return AB.middleCols(k * (nx + nu) + nx, nu); }
    crmat Bk(index_t k) const { return AB.middleCols(k * (nx + nu) + nx, nu); }

    struct {
        std::chrono::nanoseconds forward{};
        std::chrono::nanoseconds backward{};
        std::chrono::nanoseconds backward_jacobians{};
        std::chrono::nanoseconds hessians{};
    } mutable time;

    rmat Qk(index_t k) {
        assert(k <= N);
        if (k == N)
            return Q_N;
        switch (structure) {
            case CostStructure::General:
                return QRS.middleCols((nx + nu) * k, nx).topRows(nx);
            case CostStructure::DiagonalHessian:
                return QRS.middleCols(k, 1).topRows(nx);
            case CostStructure::Quadratic: return QRS.leftCols(nx).topRows(nx);
            case CostStructure::DiagonalQuadratic:
                return QRS.leftCols(1).topRows(nx);
            default: throw std::logic_error("CostStructure");
        }
    }
    rmat Rk(index_t k) {
        assert(k < N);
        switch (structure) {
            case CostStructure::General:
                return QRS.middleCols((nx + nu) * k + nx, nu).bottomRows(nu);
            case CostStructure::DiagonalHessian:
                return QRS.middleCols(k, 1).bottomRows(nu);
            case CostStructure::Quadratic:
                return QRS.rightCols(nu).bottomRows(nu);
            case CostStructure::DiagonalQuadratic:
                return QRS.rightCols(1).bottomRows(nu);
            default: throw std::logic_error("CostStructure");
        }
    }

    /// @pre `xk(0)` and `uk(k)` for `0 <= k < N` initialized
    /// @post `xk(k)` for `1 <= k <= N` updated
    /// @return @f$ V(u) = \sum_{k=0}^{N-1} \ell(x_k, u_k) + V_f(x_N) @f$
    real_t forward(rvec xu) const {
        detail::Timed t{time.forward};
        assert(xu.size() == (nx + nu) * N + nx);
        real_t V = 0;
        for (index_t k = 0; k < N; ++k) {
            V += problem.eval_l(k, xuk(xu, k));
            problem.eval_f(k, xk(xu, k), uk(xu, k), xk(xu, k + 1));
        }
        V += problem.eval_l_N(xk(xu, N));
        return V;
    }

    /// @pre `xk(0)` and `uk(k)` for `0 <= k < N` initialized
    /// @post `xk(k)` for `1 <= k <= N` updated
    /// @return @f$ V(u) = \sum_{k=0}^{N-1} \ell(x_k, u_k) + V_f(x_N) @f$
    void forward_simulate(rvec xu) const {
        detail::Timed t{time.forward};
        assert(xu.size() == (nx + nu) * N + nx);
        for (index_t k = 0; k < N; ++k)
            problem.eval_f(k, xk(xu, k), uk(xu, k), xk(xu, k + 1));
    }

    /// @pre @ref forward() has been called
    /// @post `Ak(k)` for `0 <= k < N` updated
    /// @post `Bk(k)` for `0 <= k < N` updated
    /// @post `qk(k)` for `0 <= k <= N` updated
    /// @post `rk(k)` for `0 <= k < N` updated
    /// @param[out] g @f$ \nabla V_N(u) @f$
    /// @param p Work vector of dimension @f$ n_x @f$
    void backward_with_jac(crvec xu, rvec g, rvec p) {
        detail::Timed t{time.backward_jacobians};
        assert(xu.size() == (nx + nu) * N + nx);
        problem.eval_grad_l_N(xk(xu, N), p);
        qk(N) = p;
        for (index_t t = N; t-- > 0;) {
            auto &&gt = g.segment(t * nu, nu);
            auto &&qt = qk(t);
            problem.eval_jac_f(t, xk(xu, t), uk(xu, t), ABk(t));
            if (false) {
                // TODO: gradient could be combined with forward evaluation
                problem.eval_grad_l(t, xuk(xu, t), qrk(t));
                // TODO: avoid allocations
                gt = rk(t) + Bk(t).transpose() * p;
                if (t > 0)
                    p = qt + Ak(t).transpose() * p;
                // TODO: qt is not really used here if t == 0
            } else {
                gt.noalias() = Bk(t).transpose() * p;
                if (t > 0) {
                    qt.noalias() = Ak(t).transpose() * p;
                    p            = qt;
                }
                // TODO: gradient could be combined with forward evaluation
                problem.eval_grad_l(t, xuk(xu, t), qrk(t));
                if (t > 0)
                    p += qt;
                gt += rk(t);
            }
        }
    }

    /// @pre @ref forward() has been called
    /// @param[out] g @f$ \nabla V_N(u) @f$
    /// @param p Work vector of dimension @f$ n_x @f$
    /// @param w Work vector of dimension @f$ n_x + n_u @f$
    /// @param AB Work matrix of dimensions @f$ n_x \times (n_x + n_u) @f$
#if 0
    void backward(crvec xu, rvec g, rvec p, rvec w, rmat AB) {
        assert(xu.size() == (nx + nu) * N + nx);
        problem.eval_grad_l_N(xk(xu, N), p);
        for (index_t t = N; t-- > 0;) {
            auto &&gt = g.segment(t * nu, nu);
            auto &&At = AB.leftCols(nx);
            auto &&Bt = AB.rightCols(nu);
            problem.eval_jac_f(t, xk(xu, t), uk(xu, t), AB);
            problem.eval_grad_l(t, xuk(xu, t), w);
            gt = w.bottomRows(nu) + Bt.transpose() * p;
            if (t > 0)
                p = w.topRows(nx) + At.transpose() * p;
            // TODO: t.topRows(nx) is not really used here if k == 0
        }
    }
#else
    void backward(crvec xu, rvec g, rvec p, rvec w) {
        detail::Timed t{time.backward};
        assert(xu.size() == (nx + nu) * N + nx);
        problem.eval_grad_l_N(xk(xu, N), p);
        for (index_t t = N; t-- > 0;) {
            auto &&gt = g.segment(t * nu, nu);
            problem.eval_grad_f_prod(t, xk(xu, t), uk(xu, t), p, w);
            gt = w.bottomRows(nu);
            p  = w.topRows(nx);
            problem.eval_grad_l(t, xuk(xu, t), w);
            gt += w.bottomRows(nu);
            p += w.topRows(nx);
        }
    }
#endif

    void hessians(crvec xu) {
        detail::Timed t{time.hessians};
        assert(xu.size() == (nx + nu) * N + nx);
        switch (structure) {
            case CostStructure::General:
                for (index_t t = 0; t < N; ++t) {
                    auto &&QRS_t = QRS.middleCols((nx + nu) * t, nx + nu);
                    problem.eval_hess_l(t, xuk(xu, t), QRS_t);
                }
                break;
            case CostStructure::DiagonalHessian:
                for (index_t t = 0; t < N; ++t) {
                    auto &&QRS_t = QRS.middleCols(t, 1);
                    problem.eval_hess_l(t, xuk(xu, t), QRS_t);
                }
                break;
            case CostStructure::Quadratic: [[fallthrough]];
            case CostStructure::DiagonalQuadratic:
                problem.eval_hess_l(0, xuk(xu, 0), QRS);
                break;
            default: throw std::logic_error("CostStructure");
        }
        problem.eval_hess_l_N(xk(xu, N), Q_N);
    }
};

} // namespace alpaqa