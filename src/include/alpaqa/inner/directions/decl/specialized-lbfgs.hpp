#pragma once

#include <alpaqa/inner/directions/decl/lbfgs.hpp>

namespace alpaqa {

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm that can
/// handle updates of the γ parameter.
/// @ingroup    grp_PANOCDirectionProviders
class SpecializedLBFGS {
  public:
    using Params = LBFGSParams;

    SpecializedLBFGS(Params params) : params(params) {}
    SpecializedLBFGS(Params params, size_t n, size_t history) : params(params) {
        resize(n, history);
    }

    /// Standard L-BFGS update without changing the step size γ.
    bool standard_update(crvec xₖ, crvec xₖ₊₁, crvec pₖ,
                         crvec pₖ₊₁, crvec gradₖ₊₁);
    /// L-BFGS update when changing the step size γ, recomputing everything.
    bool full_update(crvec xₖ, crvec xₖ₊₁, crvec pₖ_old_γ,
                     crvec pₖ₊₁, crvec gradₖ₊₁, const Box &C,
                     real_t γ);
    /// Update the inverse Hessian approximation using the new vectors xₖ₊₁
    /// and pₖ₊₁.
    bool update(crvec xₖ, crvec xₖ₊₁, crvec pₖ, crvec pₖ₊₁,
                crvec gradₖ₊₁, const Box &C, real_t γ);

    /// Apply the inverse Hessian approximation to the given vector q.
    template <class Vec>
    void apply(Vec &&q);

    /// Initialize with the starting point x₀ and the gradient in that point.
    void initialize(crvec x₀, crvec grad₀);

    /// Throw away the approximation and all previous vectors s and y.
    void reset();
    /// Re-allocate storage for a problem with a different size. 
    void resize(size_t n, size_t history);

    std::string get_name() const { return "SpecializedLBFGS"; }

    /// Get the size of the s, y, x and g vectors in the buffer.
    size_t n() const { return sto.rows() - 1; }
    /// Get the number of previous vectors s, y, x and g stored in the buffer.
    size_t history() const { return (sto.cols() - 2) / 4; }
    /// Get the next index in the circular buffer of previous s, y, x and g
    /// vectors.
    size_t succ(size_t i) const { return i + 1 < history() ? i + 1 : 0; }
    /// Get the previous index in the circular buffer of previous s, y, x and g
    /// vectors.
    size_t pred(size_t i) const { return i == 0 ? history() - 1 : i - 1; }

    auto s(size_t i) { return sto.col(2 * i).topRows(n()); }
    auto s(size_t i) const { return sto.col(2 * i).topRows(n()); }
    auto y(size_t i) { return sto.col(2 * i + 1).topRows(n()); }
    auto y(size_t i) const { return sto.col(2 * i + 1).topRows(n()); }
    auto x(size_t i) { return sto.col(2 * history() + 2 * i).topRows(n()); }
    auto x(size_t i) const {
        return sto.col(2 * history() + 2 * i).topRows(n());
    }
    auto g(size_t i) { return sto.col(2 * history() + 2 * i + 1).topRows(n()); }
    auto g(size_t i) const {
        return sto.col(2 * history() + 2 * i + 1).topRows(n());
    }
    auto p() { return sto.col(4 * history()).topRows(n()); }
    auto p() const { return sto.col(4 * history()).topRows(n()); }
    auto w() { return sto.col(4 * history() + 1).topRows(n()); }
    auto w() const { return sto.col(4 * history() + 1).topRows(n()); }
    real_t &ρ(size_t i) { return sto.coeffRef(n(), 2 * i); }
    const real_t &ρ(size_t i) const { return sto.coeff(n(), 2 * i); }
    real_t &α(size_t i) { return sto.coeffRef(n(), 2 * i + 1); }
    const real_t &α(size_t i) const { return sto.coeff(n(), 2 * i + 1); }

  private:
    using storage_t = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

    storage_t sto;
    size_t idx   = 0;
    bool full    = false;
    real_t old_γ = 0;
    Params params;
};

} // namespace alpaqa