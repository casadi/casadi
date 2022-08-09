#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace alpaqa {

/// Cautious BFGS update.
/// @see @ref LBFGSParams::cbfgs
template <Config Conf = DefaultConfig>
struct CBFGSParams {
    USING_ALPAQA_CONFIG(Conf);
    real_t α = 1;
    real_t ϵ = 0; ///< Set to zero to disable CBFGS check.
    explicit operator bool() const { return ϵ > 0; }
};

/// Parameters for the @ref LBFGS class.
template <Config Conf = DefaultConfig>
struct LBFGSParams {
    USING_ALPAQA_CONFIG(Conf);

    /// Length of the history to keep.
    length_t memory = 10;
    /// Reject update if @f$ y^\top s \le \text{min_div_fac} \cdot s^\top s @f$.
    real_t min_div_fac = std::numeric_limits<real_t>::epsilon();
    /// Reject update if @f$ s^\top s \le \text{min_abs_s} @f$.
    real_t min_abs_s =
        std::pow(std::numeric_limits<real_t>::epsilon(), real_t(2));
    /// Parameters in the cautious BFGS update condition
    /// @f[ \frac{y^\top s}{s^\top s} \ge \epsilon \| g \|^\alpha @f]
    /// @see https://epubs.siam.org/doi/10.1137/S1052623499354242
    CBFGSParams<config_t> cbfgs;
    /// If set to true, the inverse Hessian estimate should remain definite,
    /// i.e. a check is performed that rejects the update if
    /// @f$ y^\top s \le \text{min_div_fac} \cdot s^\top s @f$.
    /// If set to false, just try to prevent a singular Hessian by rejecting the
    /// update if
    /// @f$ \left| y^\top s \right| \le \text{min_div_fac} \cdot s^\top s @f$.
    bool force_pos_def = true;
};

/// Layout:
/// ~~~
///       ┌───── 2 m ─────┐
///     ┌ ┌───┬───┬───┬───┐
///     │ │   │   │   │   │
///     │ │ s │ y │ s │ y │
/// n+1 │ │   │   │   │   │
///     │ ├───┼───┼───┼───┤
///     │ │ ρ │ α │ ρ │ α │
///     └ └───┴───┴───┴───┘
/// ~~~
template <Config Conf = DefaultConfig>
struct LBFGSStorage {
    USING_ALPAQA_CONFIG(Conf);

    /// Re-allocate storage for a problem with a different size.
    void resize(length_t n, length_t history);

    /// Get the size of the s and y vectors in the buffer.
    length_t n() const { return sto.rows() - 1; }
    /// Get the number of previous vectors s and y stored in the buffer.
    length_t history() const { return sto.cols() / 2; }

    auto s(index_t i) { return sto.col(2 * i).topRows(n()); }
    auto s(index_t i) const {
        return std::as_const(sto).col(2 * i).topRows(n());
    }
    auto y(index_t i) { return sto.col(2 * i + 1).topRows(n()); }
    auto y(index_t i) const {
        return std::as_const(sto).col(2 * i + 1).topRows(n());
    }
    real_t &ρ(index_t i) { return sto.coeffRef(n(), 2 * i); }
    const real_t &ρ(index_t i) const {
        return std::as_const(sto).coeff(n(), 2 * i);
    }
    real_t &α(index_t i) { return sto.coeffRef(n(), 2 * i + 1); }
    real_t &α(index_t i) const { return sto.coeffRef(n(), 2 * i + 1); }

    using storage_t = mat;
    static_assert(!storage_t::IsRowMajor);
    mutable storage_t sto;
};

/// Limited memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm
/// @ingroup grp_Accelerators
template <Config Conf = DefaultConfig>
class LBFGS {
  public:
    USING_ALPAQA_CONFIG(Conf);

    using Params = LBFGSParams<config_t>;

    /// The sign of the vectors @f$ p @f$ passed to the @ref update method.
    enum class Sign {
        Positive, ///< @f$ p \sim \nabla \psi(x) @f$
        Negative, ///< @f$ p \sim -\nabla \psi(x) @f$
    };

    LBFGS(Params params) : params(params) {}
    LBFGS(Params params, length_t n) : params(params) { resize(n); }

    /// Check if the new vectors s and y allow for a valid BFGS update that
    /// preserves the positive definiteness of the Hessian approximation.
    static bool update_valid(const Params &params, real_t yᵀs, real_t sᵀs,
                             real_t pᵀp);

    /// Update the inverse Hessian approximation using the new vectors
    /// sₖ = xₙₑₓₜ - xₖ and yₖ = pₙₑₓₜ - pₖ.
    bool update_sy(crvec s, crvec y, real_t pₙₑₓₜᵀpₙₑₓₜ, bool forced = false);
    /// @see @ref update_sy
    bool update_sy_impl(const auto &s, const auto &y, real_t pₙₑₓₜᵀpₙₑₓₜ,
                        bool forced = false);

    /// Update the inverse Hessian approximation using the new vectors xₙₑₓₜ
    /// and pₙₑₓₜ.
    bool update(crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                Sign sign = Sign::Positive, bool forced = false);

    /// Apply the inverse Hessian approximation to the given vector q.
    /// Initial inverse Hessian approximation is set to @f$ H_0 = \gamma I @f$.
    /// If @p γ is negative, @f$ H_0 = \frac{s^\top y}{y^\top y} I @f$.
    bool apply(rvec q, real_t γ = -1) const;

    /// Apply the inverse Hessian approximation to the given vector q, applying
    /// only the columns and rows of the Hessian in the index set J.
    bool apply_masked(rvec q, real_t γ, crindexvec J);
    /// @copydoc apply_masked(rvec, real_t, crindexvec)
    bool apply_masked(rvec q, real_t γ, const std::vector<index_t> &J);
    /// @see @ref apply_masked
    bool apply_masked_impl(rvec q, real_t γ, const auto &J);

    /// Throw away the approximation and all previous vectors s and y.
    void reset();
    /// Re-allocate storage for a problem with a different size. Causes
    /// a @ref reset.
    void resize(length_t n);

    /// Scale the stored y vectors by the given factor.
    void scale_y(real_t factor);

    /// Get a string identifier for this accelerator.
    std::string get_name() const {
        return "LBFGS<" + std::string(config_t::get_name()) + '>';
    }
    /// Get the parameters.
    const Params &get_params() const { return params; }

    /// Get the size of the s and y vectors in the buffer.
    length_t n() const { return sto.n(); }
    /// Get the number of previous vectors s and y stored in the buffer.
    length_t history() const { return sto.history(); }
    /// Get the next index in the circular buffer of previous s and y vectors.
    index_t succ(index_t i) const { return i + 1 < history() ? i + 1 : 0; }
    /// Get the previous index in the circular buffer of s and y vectors.
    index_t pred(index_t i) const { return i > 0 ? i - 1 : history() - 1; }
    /// Get the number of previous s and y vectors currently stored in the
    /// buffer.
    length_t current_history() const { return full ? history() : idx; }

    auto s(index_t i) { return sto.s(i); }
    auto s(index_t i) const { return sto.s(i); }
    auto y(index_t i) { return sto.y(i); }
    auto y(index_t i) const { return sto.y(i); }
    real_t &ρ(index_t i) { return sto.ρ(i); }
    const real_t &ρ(index_t i) const { return sto.ρ(i); }
    real_t &α(index_t i) { return sto.α(i); }
    real_t &α(index_t i) const { return sto.α(i); }

    /// Iterate over the indices in the history buffer, oldest first.
    template <class F>
    void foreach_fwd(const F &fun) const {
        if (full)
            for (index_t i = idx; i < history(); ++i)
                fun(i);
        if (idx)
            for (index_t i = 0; i < idx; ++i)
                fun(i);
    }

    /// Iterate over the indices in the history buffer, newest first.
    template <class F>
    void foreach_rev(const F &fun) const {
        if (idx)
            for (index_t i = idx; i-- > 0;)
                fun(i);
        if (full)
            for (index_t i = history(); i-- > idx;)
                fun(i);
    }

  private:
    LBFGSStorage<config_t> sto;
    index_t idx = 0;
    bool full   = false;
    Params params;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, CBFGSParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, CBFGSParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, CBFGSParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, CBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, CBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, LBFGSParams, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, LBFGSParams, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, LBFGSParams, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, LBFGSParams, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, LBFGSParams, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, LBFGS, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, LBFGS, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, LBFGS, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, LBFGS, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, LBFGS, EigenConfigq);
#endif

} // namespace alpaqa