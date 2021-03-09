#pragma once

#include <panoc-alm/util/box.hpp>
#include <panoc-alm/util/vec.hpp>

#include <iostream>
#include <limits>

#include <Eigen/LU>

namespace pa {

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm
class LBFGS {
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    storage_t storage;
    size_t idx = 0;
    bool full  = false;

  public:
    LBFGS() = default;
    LBFGS(size_t n, size_t history) : storage(n + 1, history * 2) {
        storage.fill(std::numeric_limits<real_t>::quiet_NaN());
    }

    size_t n() const { return storage.rows() - 1; }
    size_t history() const { return storage.cols() / 2; }
    decltype(auto) s(size_t i) { return storage.col(2 * i).topRows(n()); }
    decltype(auto) s(size_t i) const { return storage.col(2 * i).topRows(n()); }
    decltype(auto) y(size_t i) { return storage.col(2 * i + 1).topRows(n()); }
    decltype(auto) y(size_t i) const {
        return storage.col(2 * i + 1).topRows(n());
    }
    decltype(auto) ρ(size_t i) { return storage.coeffRef(n(), 2 * i); }
    decltype(auto) ρ(size_t i) const { return storage.coeff(n(), 2 * i); }
    decltype(auto) α(size_t i) { return storage.coeffRef(n(), 2 * i + 1); }
    decltype(auto) α(size_t i) const { return storage.coeff(n(), 2 * i + 1); }

    template <class VecS, class VecY>
    bool update(const VecS &s, const VecY &y) {
        auto ys        = s.dot(y);
        auto norm_sq_s = s.squaredNorm();

        if (!std::isfinite(ys)) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (yᵀs = inf/NaN)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: this means it's essentially zero (or denorm)
        else if (norm_sq_s <= std::numeric_limits<real_t>::min()) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (‖s‖² <= ε)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: this means it's essentially zero (or denorm)
        else if (std::abs(ys) <= std::numeric_limits<real_t>::min()) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (yᵀs <= ε)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: CBFGS

        this->s(idx) = s;
        this->y(idx) = y;
        this->ρ(idx) = 1. / ys;

        if (++idx >= history()) {
            idx  = 0;
            full = true;
        }
        return true;
    }

    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &grad_new, Box C, real_t γ_new) {
        // TODO: CBFGS
        // bool ret = update(xₖ₊₁ - xₖ, pₖ₊₁ - pₖ); // TODO: what's with the sign here?
        bool ret = update(xₖ₊₁ - xₖ, pₖ - pₖ₊₁);
        (void)grad_new;
        (void)C;
        (void)γ_new;
        return ret;
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
    }

    void gamma_changed() { reset(); }

    void reset() {
        idx  = 0;
        full = false;
    }

    void resize(size_t n, size_t history) {
        storage.resize(n + 1, history * 2);
        reset();
    }
};

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm that can
/// handle updates of the γ parameter.
class SpecializedLBFGS {
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    storage_t storage;
    size_t idx   = 0;
    bool full    = false;
    real_t old_γ = 0;

  public:
    SpecializedLBFGS() {
        assert(!"Fix the arguments passed by PANOC before using"); // TODO
    }
    SpecializedLBFGS(size_t n, size_t history)
        : storage(n + 1, history * 4 + 1) {
        storage.fill(std::numeric_limits<real_t>::quiet_NaN());
        assert(!"Fix the arguments passed by PANOC before using"); // TODO
    }

    size_t n() const { return storage.rows() - 1; }
    size_t history() const { return (storage.cols() - 1) / 4; }
    decltype(auto) s(size_t i) { return storage.col(2 * i).topRows(n()); }
    decltype(auto) s(size_t i) const { return storage.col(2 * i).topRows(n()); }
    decltype(auto) y(size_t i) { return storage.col(2 * i + 1).topRows(n()); }
    decltype(auto) y(size_t i) const {
        return storage.col(2 * i + 1).topRows(n());
    }
    decltype(auto) x(size_t i) {
        return storage.col(2 * history() + 2 * i).topRows(n());
    }
    decltype(auto) x(size_t i) const {
        return storage.col(2 * history() + 2 * i).topRows(n());
    }
    decltype(auto) g(size_t i) {
        return storage.col(2 * history() + 2 * i + 1).topRows(n());
    }
    decltype(auto) g(size_t i) const {
        return storage.col(2 * history() + 2 * i + 1).topRows(n());
    }
    decltype(auto) x̂() { return storage.col(4 * history()).topRows(n()); }
    decltype(auto) x̂() const { return storage.col(4 * history()).topRows(n()); }
    decltype(auto) ρ(size_t i) { return storage.coeffRef(n(), 2 * i); }
    decltype(auto) ρ(size_t i) const { return storage.coeff(n(), 2 * i); }
    decltype(auto) α(size_t i) { return storage.coeffRef(n(), 2 * i + 1); }
    decltype(auto) α(size_t i) const { return storage.coeff(n(), 2 * i + 1); }

    template <class VecX, class VecG, class VecXh>
    void initialize(const VecX &new_x, const VecG &new_g, const VecXh &new_x̂,
                    real_t γ) {
        idx   = 0;
        full  = false;
        x(0)  = new_x;
        g(0)  = new_g;
        x̂()   = new_x̂;
        old_γ = γ;
    }

    size_t succ(size_t i) const { return i + 1 < history() ? i + 1 : 0; }
    size_t succ() const { return succ(idx); }

    template <class VecX, class VecG, class VecXh>
    bool update(const VecX &new_x, const VecG &new_g, const VecXh &new_x̂,
                const Box &C, real_t γ) {
        s(idx) = new_x /* i+1 */ - x(idx);
        // If the step size γ changed
        if (γ != old_γ) {
            // Recompute all residuals with new γ
            // yₖ = sₖ + x̂ₖ - x̂ₖ₊₁
            // x̂ₖ = Π(xₖ - γ∇ψ(xₖ))
            //    = Π(xₖ - γgₖ)
            size_t i   = full ? succ(idx) : 0;
            x̂(/* i */) = project(x(i) - γ * g(i), C);
            for (; i != idx; i = succ(i)) {
                y(i)         = s(i) + x̂(/* i */);
                x̂(/* i+1 */) = project(x(succ(i)) - γ * g(succ(i)), C);
                y(i) -= x̂(/* i+1 */);
            }
            y(i)  = s(i) + x̂(/* i */) - new_x̂ /* i+1 */;
            old_γ = γ;
        } else {
            y(idx) = s(idx) + x̂(/* i */) - new_x̂ /* i+1 */;
        }
        x̂()          = new_x̂ /* i+1 */;
        x(succ(idx)) = new_x;
        g(succ(idx)) = new_g;

        auto ys        = s(idx).dot(y(idx));
        auto norm_sq_s = s(idx).squaredNorm();
        ρ(idx)         = 1. / ys;

        if (!std::isfinite(ys)) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (yᵀs = inf/NaN)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: this means it's essentially zero (or denorm)
        else if (norm_sq_s <= std::numeric_limits<real_t>::min()) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (‖s‖² <= ε)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: this means it's essentially zero (or denorm)
        else if (std::abs(ys) <= std::numeric_limits<real_t>::min()) {
            // std::cerr << "[LBFGS] "
            //              "\x1b[0;31m"
            //              "Warning: L-BFGS update failed (yᵀs <= ε)"
            //              "\x1b[0m"
            //           << std::endl;
            return false;
        }
        // TODO: CBFGS

        idx = succ(idx);
        full |= idx == 0;

        return true;
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
        storage.resize(n + 1, history * 4 + 1);
        storage.fill(std::numeric_limits<real_t>::quiet_NaN());
    }

    void reset() {
        x(0) = x(idx);
        g(0) = x(idx);
        idx  = 0;
        full = false;
    }
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