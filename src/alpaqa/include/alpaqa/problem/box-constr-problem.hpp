#pragma once

#include <alpaqa/problem/box.hpp>
#include <alpaqa/util/check-dim.hpp>

namespace alpaqa {

/// Implements common problem functions for minimization problems with box
/// constraints. Meant to be used as a base class for custom problem
/// implementations.
/// @ingroup grp_Problems
template <Config Conf>
class BoxConstrProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    /// Number of decision variables, dimension of x
    length_t n;
    /// Number of constraints, dimension of g(x) and z
    length_t m;

    /// Components of the constraint function with indices below this number are
    /// handled using a quadratic penalty method rather than using an
    /// augmented Lagrangian method. Specifically, the Lagrange multipliers for
    /// these components (which determine the shifts in ALM) are kept at zero.
    index_t penalty_alm_split = 0;

    BoxConstrProblem(length_t n, ///< Number of decision variables
                     length_t m) ///< Number of constraints
        : n{n}, m{m} {}

    BoxConstrProblem(Box C, Box D, vec l1_reg = vec(0))
        : n{C.lowerbound.size()}, m{D.lowerbound.size()}, C{std::move(C)}, D{std::move(D)},
          l1_reg{std::move(l1_reg)} {}

    void resize(length_t n, length_t m) {
        this->n = n;
        this->m = m;
        C       = Box{n};
        D       = Box{m};
    }

    BoxConstrProblem(const BoxConstrProblem &)                = default;
    BoxConstrProblem &operator=(const BoxConstrProblem &)     = default;
    BoxConstrProblem(BoxConstrProblem &&) noexcept            = default;
    BoxConstrProblem &operator=(BoxConstrProblem &&) noexcept = default;

    /// Constraints of the decision variables, @f$ x \in C @f$
    Box C{this->n};
    /// Other constraints, @f$ g(x) \in D @f$
    Box D{this->m};
    /// @f$ \ell_1 @f$ (1-norm) regularization parameter.
    /// Possible dimensions are: @f$ 0 @f$ (no regularization), @f$ 1 @f$ (a
    /// single scalar factor), or @f$ n @f$ (a different factor for each
    /// variable).
    vec l1_reg{};

    /// Number of decision variables, @ref n
    length_t get_n() const { return n; }
    /// Number of constraints, @ref m
    length_t get_m() const { return m; }

    /** Projected gradient step for rectangular box C.
      * @f[ \begin{aligned} \hat x &= \Pi_C(x - \gamma\nabla\psi(x)) \\
      * p &= \hat x - x \\
      * &= \max(\underline x - x, \;\min(-\gamma\nabla\psi(x), \overline x - x)
      * \end{aligned} @f] */
    static real_t eval_proj_grad_step_box(const Box &C, real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                          rvec p) {
        p = (-γ * grad_ψ).cwiseMax(C.lowerbound - x).cwiseMin(C.upperbound - x);
        x̂ = x + p;
        return real_t(0);
    }

    /** Proximal gradient step for rectangular box C with ℓ₁-regularization.
      * @f[ \begin{aligned} h(x) &= \|x\|_1 + \delta_C(x) \\
      * \hat x &= \prox_{\gamma h}(x - \gamma\nabla\psi(x)) \\
      * &= -\max\big(
      *         x - \overline x,
      *         \;\min\big(
      *             x - \underline x,
      *             \;\min\big(
      *                 \gamma(\nabla\psi(x) + \lambda),
      *                 \;\max\big(
      *                     \gamma(\nabla\psi(x) - \lambda),
      *                     x
      *                 \big)
      *             \big)
      *         \big)
      *     \big) \end{aligned} @f] */
    static void eval_prox_grad_step_box_l1_impl(const Box &C, const auto &λ, real_t γ, crvec x,
                                                crvec grad_ψ, rvec x̂, rvec p) {
        p = -x.cwiseMax(γ * (grad_ψ - λ))
                 .cwiseMin(γ * (grad_ψ + λ))
                 .cwiseMin(x - C.lowerbound)
                 .cwiseMax(x - C.upperbound);
        x̂ = x + p;
    }
    /// @copydoc eval_prox_grad_step_box_l1_impl
    static real_t eval_prox_grad_step_box_l1(const Box &C, const auto &λ, real_t γ, crvec x,
                                             crvec grad_ψ, rvec x̂, rvec p) {
        eval_prox_grad_step_box_l1_impl(C, λ, γ, x, grad_ψ, x̂, p);
        return vec_util::norm_1(x̂.cwiseProduct(λ));
    }

    /// @copydoc eval_prox_grad_step_box_l1_impl
    static real_t eval_prox_grad_step_box_l1_scal(const Box &C, real_t λ, real_t γ, crvec x,
                                                  crvec grad_ψ, rvec x̂, rvec p) {
        auto n     = x.size();
        auto λ_vec = vec::Constant(n, λ);
        eval_prox_grad_step_box_l1_impl(C, λ_vec, γ, x, grad_ψ, x̂, p);
        return λ * vec_util ::norm_1(x̂);
    }

    /// @see @ref TypeErasedProblem::eval_prox_grad_step
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const {
        if (l1_reg.size() == 0)
            return eval_proj_grad_step_box(C, γ, x, grad_ψ, x̂, p);
        else if (l1_reg.size() == 1)
            return eval_prox_grad_step_box_l1_scal(C, l1_reg(0), γ, x, grad_ψ, x̂, p);
        else
            return eval_prox_grad_step_box_l1(C, l1_reg, γ, x, grad_ψ, x̂, p);
    }

    /// @see @ref TypeErasedProblem::eval_proj_diff_g
    void eval_proj_diff_g(crvec z, rvec p) const { p = alpaqa::projecting_difference(z, D); }

    static void eval_proj_multipliers_box(const Box &D, rvec y, real_t M,
                                          index_t penalty_alm_split) {
        // If there's no lower bound, the multipliers can only be positive
        auto max_lb = [M](real_t y, real_t z_lb) {
            real_t y_lb = z_lb == -alpaqa::inf<config_t> ? 0 : -M;
            return std::max(y, y_lb);
        };
        // If there's no upper bound, the multipliers can only be negative
        auto min_ub = [M](real_t y, real_t z_ub) {
            real_t y_ub = z_ub == alpaqa::inf<config_t> ? 0 : M;
            return std::min(y, y_ub);
        };
        auto num_alm    = y.size() - penalty_alm_split;
        auto &&y_qpm    = y.topRows(penalty_alm_split);
        auto &&y_alm    = y.bottomRows(num_alm);
        auto &&z_alm_lb = D.lowerbound.bottomRows(num_alm);
        auto &&z_alm_ub = D.upperbound.bottomRows(num_alm);
        y_qpm.setZero();
        y_alm = y_alm.binaryExpr(z_alm_lb, max_lb).binaryExpr(z_alm_ub, min_ub);
    }

    /// @see @ref TypeErasedProblem::eval_proj_multipliers
    void eval_proj_multipliers(rvec y, real_t M) const {
        eval_proj_multipliers_box(D, y, M, penalty_alm_split);
    }

    /// @see @ref TypeErasedProblem::get_box_C
    const Box &get_box_C() const { return C; }
    /// @see @ref TypeErasedProblem::get_box_D
    const Box &get_box_D() const { return D; }

    /// @see @ref TypeErasedProblem::eval_inactive_indices_res_lna
    index_t eval_inactive_indices_res_lna(real_t γ, crvec x, crvec grad_ψ, rindexvec J) const {
        index_t nJ = 0;
        // Helper that adds i to index set J if x ∊ C
        const auto add_to_J_if_in_box_interior = [&](real_t x_fw, index_t i) {
            if (C.lowerbound(i) < x_fw && x_fw < C.upperbound(i))
                J(nJ++) = i;
        };
        // Update the index set J for the general box + l1 case
        const auto update_J_general = [&](real_t λ, real_t x_fw, index_t i) {
            if (λ == 0) {
                add_to_J_if_in_box_interior(x_fw, i);
            } else {
                if (x_fw > γ * λ)
                    add_to_J_if_in_box_interior(x_fw - γ * λ, i);
                else if (x_fw < -γ * λ)
                    add_to_J_if_in_box_interior(x_fw + γ * λ, i);
            }
        };
        const auto nλ     = l1_reg.size();
        const bool λ_is_0 = nλ == 0 || (nλ == 1 && l1_reg(0) == 0);
        // Only box constraints
        if (λ_is_0)
            for (index_t i = 0; i < n; ++i) {
                real_t x_fw = x(i) - γ * grad_ψ(i);
                add_to_J_if_in_box_interior(x_fw, i);
            }
        // Box constraints and l1
        else
            for (index_t i = 0; i < n; ++i) {
                real_t λi   = nλ == 0 ? 0 : nλ == 1 ? l1_reg(0) : l1_reg(i);
                real_t x_fw = x(i) - γ * grad_ψ(i);
                update_J_general(λi, x_fw, i);
            }
        return nJ;
    }

    /// @see @ref TypeErasedProblem::check
    void check() const {
        util::check_dim_msg<config_t>(
            C.lowerbound, n,
            "Length of problem.C.lowerbound does not match problem size problem.n");
        util::check_dim_msg<config_t>(
            C.upperbound, n,
            "Length of problem.C.upperbound does not match problem size problem.n");
        util::check_dim_msg<config_t>(
            D.lowerbound, m,
            "Length of problem.D.lowerbound does not match problem size problem.m");
        util::check_dim_msg<config_t>(
            D.upperbound, m,
            "Length of problem.D.upperbound does not match problem size problem.m");
        if (l1_reg.size() > 1)
            util::check_dim_msg<config_t>(
                l1_reg, n,
                "Length of problem.l1_reg does not match problem size problem.n, 1 or 0");
        if (penalty_alm_split < 0 || penalty_alm_split > m)
            throw std::invalid_argument("Invalid penalty_alm_split");
    }
};

} // namespace alpaqa
