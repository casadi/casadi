#pragma once

#include "vec.hpp"
#include <iostream>

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
    void update(const VecS &s, const VecY &y) {
        auto ys        = s.dot(y);
        auto norm_sq_s = s.squaredNorm();

        if (!std::isfinite(ys)) {
            std::cerr << "[LBFGS] "
                         "\x1b[0;31m"
                         "Warning: L-BFGS update failed (yᵀs = inf/NaN)"
                         "\x1b[0m"
                      << std::endl;
            return;
        } else if (norm_sq_s <= std::numeric_limits<real_t>::min()) {
            std::cerr << "[LBFGS] "
                         "\x1b[0;31m"
                         "Warning: L-BFGS update failed (‖s‖² <= ε)"
                         "\x1b[0m"
                      << std::endl;
            return;
        } else if (std::abs(ys) <= std::numeric_limits<real_t>::min()) {
            std::cerr << "[LBFGS] "
                         "\x1b[0;31m"
                         "Warning: L-BFGS update failed (yᵀs <= ε)"
                         "\x1b[0m"
                      << std::endl;
            return;
        }
        // TODO: CBFGS

        this->s(idx)   = s;
        this->y(idx)   = y;
        this->ρ(idx)   = 1. / ys;

        if (++idx >= history()) {
            idx  = 0;
            full = true;
        }
    }

    template <class VecOut>
    void apply_work(real_t γ, vec &q, VecOut &r) {
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

        // r ← H₀q
        r = γ * q; // TODO: diagonal matrix H₀

        auto update2 = [&](size_t i) {
            real_t β = ρ(i) * (y(i).dot(r));
            r += (α(i) - β) * s(i);
        };
        if (full)
            for (size_t i = idx; i < history(); ++i)
                update2(i);
        for (size_t i = 0; i < idx; ++i)
            update2(i);
    }

    template <class VecOut>
    void apply(real_t γ, vec v, VecOut &&Hv) {
        apply_work(γ, v, Hv);
    }

    void reset() {
        idx  = 0;
        full = false;
    }

    void resize(size_t n, size_t history) {
        storage.resize(n + 1, history * 2);
        reset();
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