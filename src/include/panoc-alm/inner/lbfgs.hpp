#pragma once

#include <cmath>
#include <panoc-alm/util/box.hpp>
#include <panoc-alm/util/vec.hpp>

#include <panoc-alm/inner/decl/lbfgs-fwd.hpp>
#include <panoc-alm/inner/detail/panoc-helpers.hpp>

#include <iostream>
#include <limits>

#include <Eigen/LU>

namespace pa {

struct LBFGSParams {
    struct {
        real_t α = 1;
        real_t ε = 1e-10;
    } cbfgs;
};

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm
class LBFGS {
  public:
    using Params = LBFGSParams;

    LBFGS(Params params) : params(params) {}
    LBFGS(Params params, size_t n, size_t history)
        : sto(n + 1, history * 2), params(params) {
        sto.fill(std::numeric_limits<real_t>::quiet_NaN());
    }

    size_t n() const { return sto.rows() - 1; }
    size_t history() const { return sto.cols() / 2; }
    decltype(auto) s(size_t i) { return sto.col(2 * i).topRows(n()); }
    decltype(auto) s(size_t i) const { return sto.col(2 * i).topRows(n()); }
    decltype(auto) y(size_t i) { return sto.col(2 * i + 1).topRows(n()); }
    decltype(auto) y(size_t i) const { return sto.col(2 * i + 1).topRows(n()); }
    decltype(auto) ρ(size_t i) { return sto.coeffRef(n(), 2 * i); }
    decltype(auto) ρ(size_t i) const { return sto.coeff(n(), 2 * i); }
    decltype(auto) α(size_t i) { return sto.coeffRef(n(), 2 * i + 1); }
    decltype(auto) α(size_t i) const { return sto.coeff(n(), 2 * i + 1); }

    size_t succ(size_t i) const { return i + 1 < history() ? i + 1 : 0; }

    static bool update_valid(LBFGSParams params, real_t yᵀs, real_t sᵀs,
                             real_t pᵀp) {
        // Smallest number we want to divide by without overflow
        const real_t min_divisor =
            std::sqrt(std::numeric_limits<real_t>::min());

        // Check if this L-BFGS update is accepted
        if (not std::isfinite(yᵀs))
            return false;
        if (sᵀs < min_divisor)
            return false;
        if (yᵀs < min_divisor)
            return false;

        // CBFGS condition: https://epubs.siam.org/doi/10.1137/S1052623499354242
        real_t α = params.cbfgs.α;
        real_t ε = params.cbfgs.ε;
        // Condition: yᵀs / sᵀs >= ε ‖p‖^α
        bool cbfgs_cond = yᵀs / sᵀs >= ε * std::pow(pᵀp, α / 2);
        if (not cbfgs_cond)
            return false;

        return true;
    }

    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ,
                const vec &pₖ₊₁) {
        const auto s = xₖ₊₁ - xₖ;
        const auto y = pₖ - pₖ₊₁;
        real_t yᵀs   = y.dot(s);
        real_t sᵀs   = s.squaredNorm();
        real_t pᵀp   = pₖ₊₁.squaredNorm();
        real_t ρ     = 1 / yᵀs;

        if (not update_valid(params, yᵀs, sᵀs, pᵀp))
            return false;

        // Store the new s and y vectors
        this->s(idx) = s;
        this->y(idx) = y;
        this->ρ(idx) = ρ;

        // Increment the index in the circular buffer
        idx = succ(idx);
        full |= idx == 0;

        return true;
    }

    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &grad_new, const Box &C, real_t γ_new) {
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return update(xₖ, xₖ₊₁, pₖ, pₖ₊₁);
    }

    template <class Vec>
    void apply(Vec &&q) {
        auto update1 = [&](size_t i) {
            α(i) = ρ(i) * (s(i).dot(q));
            q -= α(i) * y(i);
        };
        if (idx)
            for (size_t i = idx; i-- > 0;)
                update1(i);
        if (full)
            for (size_t i = history(); i-- > idx;)
                update1(i);

        // q = H₀ * q; // TODO: diagonal matrix H₀?

        auto update2 = [&](size_t i) {
            real_t β = ρ(i) * (y(i).dot(q));
            q += (α(i) - β) * s(i);
        };
        if (full)
            for (size_t i = idx; i < history(); ++i)
                update2(i);
        for (size_t i = 0; i < idx; ++i)
            update2(i);
    }

    void initialize(const vec &x, const vec &p, const vec &grad_ψ, real_t γ) {
        (void)x;
        (void)p;
        (void)grad_ψ;
        (void)γ;
        std::cout << "init: " << n() << ", " << history() << std::endl;
    }

    void gamma_changed() { reset(); }

    void reset() {
        idx  = 0;
        full = false;
    }

    void resize(size_t n, size_t history) {
        sto.resize(n + 1, history * 2);
        reset();
    }

  private:
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    storage_t sto;
    size_t idx = 0;
    bool full  = false;
    Params params;
};

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm that can
/// handle updates of the γ parameter.
class SpecializedLBFGS {
  public:
    using Params = LBFGSParams;

    SpecializedLBFGS(Params params) : params(params) {}
    SpecializedLBFGS(Params params, size_t n, size_t history)
        : sto(n + 1, history * 4 + 2), params(params) {
        sto.fill(std::numeric_limits<real_t>::quiet_NaN());
    }

    size_t n() const { return sto.rows() - 1; }
    size_t history() const { return (sto.cols() - 2) / 4; }
    decltype(auto) s(size_t i) { return sto.col(2 * i).topRows(n()); }
    decltype(auto) s(size_t i) const { return sto.col(2 * i).topRows(n()); }
    decltype(auto) y(size_t i) { return sto.col(2 * i + 1).topRows(n()); }
    decltype(auto) y(size_t i) const { return sto.col(2 * i + 1).topRows(n()); }
    decltype(auto) x(size_t i) {
        return sto.col(2 * history() + 2 * i).topRows(n());
    }
    decltype(auto) x(size_t i) const {
        return sto.col(2 * history() + 2 * i).topRows(n());
    }
    decltype(auto) g(size_t i) {
        return sto.col(2 * history() + 2 * i + 1).topRows(n());
    }
    decltype(auto) g(size_t i) const {
        return sto.col(2 * history() + 2 * i + 1).topRows(n());
    }
    decltype(auto) p() { return sto.col(4 * history()).topRows(n()); }
    decltype(auto) p() const { return sto.col(4 * history()).topRows(n()); }
    decltype(auto) w() { return sto.col(4 * history() + 1).topRows(n()); }
    decltype(auto) w() const { return sto.col(4 * history() + 1).topRows(n()); }
    decltype(auto) ρ(size_t i) { return sto.coeffRef(n(), 2 * i); }
    decltype(auto) ρ(size_t i) const { return sto.coeff(n(), 2 * i); }
    decltype(auto) α(size_t i) { return sto.coeffRef(n(), 2 * i + 1); }
    decltype(auto) α(size_t i) const { return sto.coeff(n(), 2 * i + 1); }

    void initialize(const vec &x₀, const vec &p, const vec &grad₀, real_t γ) {
        idx   = 0;
        full  = false;
        x(0)  = x₀;
        g(0)  = grad₀;
        old_γ = γ;
        (void)p; // TODO: remove this parameter?
    }

    size_t succ(size_t i) const { return i + 1 < history() ? i + 1 : 0; }
    size_t pred(size_t i) const { return i == 0 ? history() - 1 : i - 1; }

    /// Standard L-BFGS update without changing the step size γ.
    bool standard_update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ,
                         const vec &pₖ₊₁, const vec &gradₖ₊₁) {
        const auto s = xₖ₊₁ - xₖ;
        const auto y = pₖ - pₖ₊₁;

        real_t yᵀs = y.dot(s);
        real_t sᵀs = s.squaredNorm();
        real_t pᵀp = pₖ₊₁.squaredNorm();
        real_t ρ   = 1 / yᵀs;

        if (not LBFGS::update_valid(params, yᵀs, sᵀs, pᵀp))
            return false;

        // Store the new s and y vectors
        this->s(idx) = s;
        this->y(idx) = y;
        this->ρ(idx) = ρ;

        // Store x and the gradient
        this->x(succ(idx)) = xₖ₊₁;
        this->g(succ(idx)) = gradₖ₊₁;

        // Increment the index in the circular buffer
        idx = succ(idx);
        full |= idx == 0;

        return true;
    }

    /// L-BFGS update when changing the step size γ, recomputing everything.
    bool full_update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ_old_γ,
                     const vec &pₖ₊₁, const vec &gradₖ₊₁, const Box &C,
                     real_t γ) {
        auto &&sₖ = xₖ₊₁ - xₖ;
        auto &&yₖ = this->w(); // temporary workspace
        // Old pₖ is no longer valid, recompute with new γ
        (void)pₖ_old_γ;
        auto &&pₖ = this->p();
        pₖ        = detail::projected_gradient_step(C, γ, x(idx), g(idx));
        yₖ        = pₖ - pₖ₊₁;

        assert(x(idx) == xₖ);

        real_t yᵀs = yₖ.dot(sₖ);
        real_t sᵀs = sₖ.squaredNorm();
        real_t pᵀp = pₖ₊₁.squaredNorm();
        real_t ρₖ  = 1 / yᵀs;

        if (not LBFGS::update_valid(params, yᵀs, sᵀs, pᵀp))
            return false;

        old_γ = γ;

        // Recompute all residuals with new γ
        // yₖ = pₖ - pₖ₊₁
        // pₖ = Π(-γ∇ψ(xₖ), C - xₖ)
        size_t endidx = full ? idx : pred(0);
        for (size_t i = pred(idx); i != endidx; i = pred(i)) {
            this->y(i) = -pₖ /* i+1 */;
            pₖ = detail::projected_gradient_step(C, γ, this->x(i), this->g(i));
            this->y(i) += pₖ /* i */;
        }
        // Store the new s and y vectors
        this->s(idx) = sₖ;
        this->y(idx) = yₖ;
        this->ρ(idx) = ρₖ;

        // Store x and the gradient
        this->x(succ(idx)) = xₖ₊₁;
        this->g(succ(idx)) = gradₖ₊₁;

        // Increment the index in the circular buffer
        idx = succ(idx);
        full |= idx == 0;

        return true;
    }

    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &gradₖ₊₁, const Box &C, real_t γ) {
        return γ == old_γ ? standard_update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, gradₖ₊₁)
                          : full_update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, gradₖ₊₁, C, γ);
    }

    template <class Vec>
    void apply(Vec &&q) {
        auto update1 = [&](size_t i) {
            α(i) = ρ(i) * (s(i).dot(q));
            q -= α(i) * y(i);
        };
        if (idx)
            for (size_t i = idx; i-- > 0;)
                update1(i);
        if (full)
            for (size_t i = history(); i-- > idx;)
                update1(i);

        // q = H₀ * q; // TODO: diagonal matrix H₀?

        auto update2 = [&](size_t i) {
            real_t β = ρ(i) * (y(i).dot(q));
            q += (α(i) - β) * s(i);
        };
        if (full)
            for (size_t i = idx; i < history(); ++i)
                update2(i);
        for (size_t i = 0; i < idx; ++i)
            update2(i);
    }

    void resize(size_t n, size_t history) {
        sto.resize(n + 1, history * 4 + 2);
        sto.fill(std::numeric_limits<real_t>::quiet_NaN());
    }

    void reset() {
        x(0) = x(idx);
        g(0) = x(idx);
        idx  = 0;
        full = false;
    }

    void gamma_changed() {}

  private:
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    storage_t sto;
    size_t idx   = 0;
    bool full    = false;
    real_t old_γ = 0;
    Params params;
};

/// For tests only, to compare L-BFGS and BFGS
class BFGS {
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;
    storage_t B;

  public:
    BFGS() = default;
    BFGS(size_t n) : B(n, n) { reset(); }

    void reset() { B = storage_t::Identity(B.rows(), B.cols()); }

    template <class VecOut>
    void apply(real_t γ, vec v, VecOut &&Hv) {
        Hv = γ * B.partialPivLu().solve(v);
    }

    template <class VecS, class VecY>
    void update(const VecS &s, const VecY &y) {
        B = B + y * y.transpose() / y.dot(s) -
            (B * s) * (s.transpose() * B.transpose()) / (s.transpose() * B * s);
    }
};

} // namespace pa