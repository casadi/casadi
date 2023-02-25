#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/util/alloc-check.hpp>
#include <alpaqa/util/check-dim.hpp>
#include <alpaqa/util/not-implemented.hpp>
#include <alpaqa/util/required-method.hpp>
#include <alpaqa/util/type-erasure.hpp>
#include <chrono>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace alpaqa {

template <Config Conf>
struct ProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    template <class F>
    using optional_function_t = util::BasicVTable::optional_function_t<F, ProblemVTable>;
    template <class F>
    using optional_const_function_t =
        util::BasicVTable::optional_const_function_t<F, ProblemVTable>;

    // clang-format off

    // Required
    required_const_function_t<void(crvec z, rvec p)>
        eval_proj_diff_g;
    required_const_function_t<void(rvec y, real_t M, index_t penalty_alm_split)>
        eval_proj_multipliers;
    required_const_function_t<real_t(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p)>
        eval_prox_grad_step;
    required_const_function_t<real_t(crvec x)>
        eval_f;
    required_const_function_t<void(crvec x, rvec grad_fx)>
        eval_grad_f;
    required_const_function_t<void(crvec x, rvec gx)>
        eval_g;
    required_const_function_t<void(crvec x, crvec y, rvec grad_gxy)>
        eval_grad_g_prod;

    // Second order
    optional_const_function_t<void(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values)>
        eval_jac_g = &default_eval_jac_g;
    optional_const_function_t<length_t()>
        get_jac_g_num_nonzeros = &default_get_jac_g_num_nonzeros;
    optional_const_function_t<void(crvec x, index_t i, rvec grad_gi)>
        eval_grad_gi = &default_eval_grad_gi;
    optional_const_function_t<void(crvec x, crvec y, real_t scale, crvec v, rvec Hv)>
        eval_hess_L_prod = &default_eval_hess_L_prod;
    optional_const_function_t<void(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values)>
        eval_hess_L = &default_eval_hess_L;
    optional_const_function_t<length_t()>
        get_hess_L_num_nonzeros = &default_get_hess_L_num_nonzeros;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv)>
        eval_hess_ψ_prod = &default_eval_hess_ψ_prod;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values)>
        eval_hess_ψ = &default_eval_hess_ψ;
    optional_const_function_t<length_t()>
        get_hess_ψ_num_nonzeros = &default_get_hess_ψ_num_nonzeros;

    // Combined evaluations
    optional_const_function_t<real_t(crvec x, rvec grad_fx)>
        eval_f_grad_f = &default_eval_f_grad_f;
    optional_const_function_t<real_t(crvec x, rvec g)>
        eval_f_g = &default_eval_f_g;
    optional_const_function_t<real_t(crvec x, rvec grad_fx, rvec g)>
        eval_f_grad_f_g = &default_eval_f_grad_f_g;
    optional_const_function_t<void(crvec x, crvec y, rvec grad_f, rvec grad_gxy)>
        eval_grad_f_grad_g_prod = &default_eval_grad_f_grad_g_prod;

    // Lagrangian and augmented lagrangian evaluations
    optional_const_function_t<void(crvec x, crvec y, rvec grad_L, rvec work_n)>
        eval_grad_L = &default_eval_grad_L;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec ŷ)>
        eval_ψ = &default_eval_ψ;
    optional_const_function_t<void(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n)>
        eval_grad_ψ_from_ŷ = &default_eval_grad_ψ_from_ŷ;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_grad_ψ = &default_eval_grad_ψ;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_ψ_grad_ψ = &default_eval_ψ_grad_ψ;

    // Constraint sets
    optional_const_function_t<const Box &()>
        get_box_C = &default_get_box_C;
    optional_const_function_t<const Box &()>
        get_box_D = &default_get_box_D;

    // Check
    required_const_function_t<void()>
        check;

    // clang-format on

    length_t n, m;

    static real_t calc_ŷ_dᵀŷ(const void *self, rvec g_ŷ, crvec y, crvec Σ,
                             const ProblemVTable &vtable) {
        if (Σ.size() == 1) {
            // ζ = g(x) + Σ⁻¹y
            g_ŷ += (1 / Σ(0)) * y;
            // d = ζ - Π(ζ, D)
            vtable.eval_proj_diff_g(self, g_ŷ, g_ŷ);
            // dᵀŷ, ŷ = Σ d
            real_t dᵀŷ = Σ(0) * g_ŷ.dot(g_ŷ);
            g_ŷ *= Σ(0);
            return dᵀŷ;
        } else {
            // ζ = g(x) + Σ⁻¹y
            g_ŷ += Σ.asDiagonal().inverse() * y;
            // d = ζ - Π(ζ, D)
            vtable.eval_proj_diff_g(self, g_ŷ, g_ŷ);
            // dᵀŷ, ŷ = Σ d
            real_t dᵀŷ = 0;
            for (index_t i = 0; i < y.size(); ++i) {
                dᵀŷ += g_ŷ(i) * Σ(i) * g_ŷ(i); // TODO: vectorize
                g_ŷ(i) = Σ(i) * g_ŷ(i);
            }
            return dᵀŷ;
        }
    }

    static void default_eval_jac_g(const void *, crvec, rindexvec, rindexvec, rvec,
                                   const ProblemVTable &) {
        throw not_implemented_error("eval_jac_g");
    }
    static length_t default_get_jac_g_num_nonzeros(const void *, const ProblemVTable &) {
        return 0;
    }
    static void default_eval_grad_gi(const void *, crvec, index_t, rvec, const ProblemVTable &) {
        throw not_implemented_error("eval_grad_gi");
    }
    static void default_eval_hess_L_prod(const void *, crvec, crvec, real_t, crvec, rvec,
                                         const ProblemVTable &) {
        throw not_implemented_error("eval_hess_L_prod");
    }
    static void default_eval_hess_L(const void *, crvec, crvec, real_t, rindexvec, rindexvec, rvec,
                                    const ProblemVTable &) {
        throw not_implemented_error("eval_hess_L");
    }
    static length_t default_get_hess_L_num_nonzeros(const void *, const ProblemVTable &) {
        return 0;
    }
    static void default_eval_hess_ψ_prod(const void *self, crvec x, crvec y, crvec, real_t scale,
                                         crvec v, rvec Hv, const ProblemVTable &vtable) {
        if (y.size() == 0 && vtable.eval_hess_L_prod != default_eval_hess_L_prod)
            return vtable.eval_hess_L_prod(self, x, y, scale, v, Hv, vtable);
        throw not_implemented_error("eval_hess_ψ_prod");
    }
    static void default_eval_hess_ψ(const void *self, crvec x, crvec y, crvec, real_t scale,
                                    rindexvec inner_idx, rindexvec outer_ptr, rvec H_values,
                                    const ProblemVTable &vtable) {
        if (y.size() == 0 && vtable.eval_hess_L != default_eval_hess_L)
            return vtable.eval_hess_L(self, x, y, scale, inner_idx, outer_ptr, H_values, vtable);
        throw not_implemented_error("eval_hess_ψ");
    }
    static length_t default_get_hess_ψ_num_nonzeros(const void *, const ProblemVTable &) {
        return 0;
    }

    static real_t default_eval_f_grad_f(const void *self, crvec x, rvec grad_fx,
                                        const ProblemVTable &vtable) {
        vtable.eval_grad_f(self, x, grad_fx);
        return vtable.eval_f(self, x);
    }
    static real_t default_eval_f_g(const void *self, crvec x, rvec g, const ProblemVTable &vtable) {
        vtable.eval_g(self, x, g);
        return vtable.eval_f(self, x);
    }
    static real_t default_eval_f_grad_f_g(const void *self, crvec x, rvec grad_fx, rvec g,
                                          const ProblemVTable &vtable) {
        vtable.eval_g(self, x, g);
        return vtable.eval_f_grad_f(self, x, grad_fx, vtable);
    }
    static void default_eval_grad_f_grad_g_prod(const void *self, crvec x, crvec y, rvec grad_f,
                                                rvec grad_gxy, const ProblemVTable &vtable) {
        vtable.eval_grad_f(self, x, grad_f);
        vtable.eval_grad_g_prod(self, x, y, grad_gxy);
    }
    static void default_eval_grad_L(const void *self, crvec x, crvec y, rvec grad_L, rvec work_n,
                                    const ProblemVTable &vtable) {
        vtable.eval_grad_f_grad_g_prod(self, x, y, grad_L, work_n, vtable);
        grad_L += work_n;
    }
    static real_t default_eval_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec ŷ,
                                 const ProblemVTable &vtable) {
        if (y.size() == 0) /* [[unlikely]] */
            return vtable.eval_f(self, x);

        real_t f   = vtable.eval_f_g(self, x, ŷ, vtable);
        real_t dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable);
        // ψ(x) = f(x) + ½ dᵀŷ
        real_t ψ = f + real_t(0.5) * dᵀŷ;
        return ψ;
    }
    static void default_eval_grad_ψ_from_ŷ(const void *self, crvec x, crvec ŷ, rvec grad_ψ,
                                           rvec work_n, const ProblemVTable &vtable) {
        if (ŷ.size() == 0) /* [[unlikely]] */
            vtable.eval_grad_f(self, x, grad_ψ);
        else
            vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable);
    }
    static void default_eval_grad_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                    rvec work_n, rvec work_m, const ProblemVTable &vtable) {
        if (y.size() == 0) /* [[unlikely]] */ {
            vtable.eval_grad_f(self, x, grad_ψ);
        } else {
            vtable.eval_g(self, x, work_m);
            (void)calc_ŷ_dᵀŷ(self, work_m, y, Σ, vtable);
            vtable.eval_grad_ψ_from_ŷ(self, x, work_m, grad_ψ, work_n, vtable);
        }
    }
    static real_t default_eval_ψ_grad_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                        rvec work_n, rvec work_m, const ProblemVTable &vtable) {
        if (y.size() == 0) /* [[unlikely]] */
            return vtable.eval_f_grad_f(self, x, grad_ψ, vtable);

        auto &ŷ = work_m;
        // ψ(x) = f(x) + ½ dᵀŷ
        real_t f   = vtable.eval_f_g(self, x, ŷ, vtable);
        real_t dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable);
        real_t ψ   = f + real_t(0.5) * dᵀŷ;
        // ∇ψ(x) = ∇f(x) + ∇g(x) ŷ
        vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable);
        return ψ;
    }
    static const Box &default_get_box_C(const void *, const ProblemVTable &) {
        throw not_implemented_error("get_box_C");
    }
    static const Box &default_get_box_D(const void *, const ProblemVTable &) {
        throw not_implemented_error("get_box_D");
    }

    template <class P>
    ProblemVTable(util::VTableTypeTag<P> t) : util::BasicVTable{t} {
        auto &vtable = *this;
        assert(t.t);

        // Initialize all methods

        // Required
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_diff_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_multipliers);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_prox_grad_step);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_g_prod);
        // Second order
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_jac_g, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_jac_g_num_nonzeros, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_gi, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L_prod, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_hess_L_num_nonzeros, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ_prod, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_hess_ψ_num_nonzeros, t.t);
        // Combined evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_g, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f_g, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_f_grad_g_prod, t.t);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_L, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ_from_ŷ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ_grad_ψ, t.t);
        // Constraint set
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_C, t.t);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_D, t.t);
        // Check
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, check);

        // Dimensions
        vtable.n = t.t->get_n();
        vtable.m = t.t->get_m();
    }
    ProblemVTable() = default;
};

/// @addtogroup grp_Problems
/// @{

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedProblem : public util::TypeErased<ProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box            = alpaqa::Box<config_t>;
    using VTable         = ProblemVTable<config_t>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using TypeErased::TypeErased;

  protected:
    using TypeErased::call;
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TypeErasedProblem make(Args &&...args) {
        return TypeErased::template make<TypeErasedProblem, T>(std::forward<Args>(args)...);
    }

    /// @name Problem dimensions
    /// @{

    /// **[Required]**
    /// Number of decision variables.
    length_t get_n() const;
    /// **[Required]**
    /// Number of constraints.
    length_t get_m() const;

    /// @}

    /// @name Required cost and constraint functions
    /// @{

    /// **[Required]**
    /// Function that evaluates the cost, @f$ f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    real_t eval_f(crvec x) const;
    /// **[Required]**
    /// Function that evaluates the gradient of the cost, @f$ \nabla f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] grad_fx
    ///         Gradient of cost function @f$ \nabla f(x) \in \R^n @f$
    void eval_grad_f(crvec x, rvec grad_fx) const;
    /// **[Required]**
    /// Function that evaluates the constraints, @f$ g(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] gx
    ///         Value of the constraints @f$ g(x) \in \R^m @f$
    void eval_g(crvec x, rvec gx) const;
    /// **[Required]**
    /// Function that evaluates the gradient of the constraints times a vector,
    /// @f$ \nabla g(x)\,y = \tp{\jac_g(x)}y @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Vector @f$ y \in \R^m @f$ to multiply the gradient by
    /// @param  [out] grad_gxy
    ///         Gradient of the constraints
    ///         @f$ \nabla g(x)\,y \in \R^n @f$
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;

    /// @}

    /// @name Projections onto constraint sets and proximal mappings
    /// @{

    /// **[Required]**
    /// Function that evaluates the difference between the given point @f$ z @f$
    /// and its projection onto the constraint set @f$ D @f$.
    /// @param  [in] z
    ///         Slack variable, @f$ z \in \R^m @f$
    /// @param  [out] p
    ///         The difference relative to its projection,
    ///         @f$ p = z - \Pi_D(z) \in \R^m @f$
    /// @note   @p z and @p p can refer to the same vector.
    void eval_proj_diff_g(crvec z, rvec p) const;
    /// **[Required]**
    /// Function that projects the Lagrange multipliers for ALM.
    /// @param  [inout] y
    ///         Multipliers, @f$ y \leftarrow \Pi_Y(y) \in \R^m @f$
    /// @param  [in] M
    ///         The radius/size of the set @f$ Y @f$. See @ref ALMParams::M.
    /// @param  [in] penalty_alm_split
    ///         See @ref ALMParams::penalty_alm_split.
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const;
    /// **[Required]**
    /// Function that evaluates a proximal gradient step.
    /// @param  [in] γ
    ///         Step size, @f$ \gamma \in \R_{>0} @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] grad_ψ
    ///         Gradient of the subproblem cost, @f$ \nabla\psi(x) \in \R^n @f$
    /// @param  [out] x̂
    ///         Next proximal gradient iterate, @f$ \hat x = T_\gamma(x) =
    ///         \prox_{\gamma h}(x - \gamma\nabla\psi(x)) \in \R^n @f$
    /// @param  [out] p
    ///         The proximal gradient step,
    ///         @f$ p = \hat x - x \in \R^n @f$
    /// @return The nonsmooth function evaluated at x̂,
    ///         @f$ h(\hat x) @f$.
    /// @note   The vector @f$ p @f$ is often used in stopping criteria, so its
    ///         numerical accuracy is more important than that of @f$ \hat x @f$.
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const;

    /// @}

    /// @name Constraint sets
    /// @{

    /// **[Optional]**
    /// Get the rectangular constraint set of the decision variables,
    /// @f$ x \in C @f$.
    const Box &get_box_C() const;
    /// **[Optional]**
    /// Get the rectangular constraint set of the general constraint function,
    /// @f$ g(x) \in D @f$.
    const Box &get_box_D() const;

    /// @}

    /// @name Functions for second-order solvers
    /// @{

    /// **[Optional]**
    /// Function that evaluates the Jacobian of the constraints as a sparse
    /// matrix, @f$ \jac_g(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [inout] inner_idx
    ///         Inner indices (row indices of nonzeros).
    /// @param  [inout] outer_ptr
    ///         Outer pointers (points to the first nonzero in each column).
    /// @param  [out] J_values
    ///         Nonzero values of the Jacobian
    ///         @f$ \jac_g(x) \in \R^{m\times n} @f$
    /// If @p J_values has size zero, this function should initialize
    /// @p inner_idx and @p outer_ptr. If @p J_values is nonempty, @p inner_idx
    /// and @p outer_ptr can be assumed to be initialized, and this function
    /// should evaluate @p J_values.
    ///
    /// Required for second-order solvers only.
    void eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const;
    /// **[Optional]**
    /// Function that gets the number of nonzeros of the Jacobian of the
    /// constraints.
    ///
    /// Required for second-order solvers only.
    length_t get_jac_g_num_nonzeros() const;
    /// **[Optional]**
    /// Function that evaluates the gradient of one specific constraint,
    /// @f$ \nabla g_i(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] i
    ///         Which constraint @f$ 0 \le i \lt m @f$
    /// @param  [out] grad_gi
    ///         Gradient of the constraint
    ///         @f$ \nabla g_i(x) \in \R^n @f$
    ///
    /// Required for second-order solvers only.
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the Lagrangian multiplied by a
    /// vector,
    /// @f$ \nabla_{xx}^2L(x, y)\,v @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] scale
    ///         Scale factor for the cost function.
    /// @param  [in] v
    ///         Vector to multiply by @f$ v \in \R^n @f$
    /// @param  [out] Hv
    ///         Hessian-vector product
    ///         @f$ \nabla_{xx}^2 L(x, y)\,v \in \R^{n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the Lagrangian as a sparse matrix,
    /// @f$ \nabla_{xx}^2L(x, y) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] scale
    ///         Scale factor for the cost function.
    /// @param  [inout] inner_idx
    ///         Inner indices (row indices of nonzeros).
    /// @param  [inout] outer_ptr
    ///         Outer pointers (points to the first nonzero in each column).
    /// @param  [out] H_values
    ///         Nonzero values of the Hessian
    ///         @f$ \nabla_{xx}^2 L(x, y) \in \R^{n\times n} @f$.
    /// If @p H_values has size zero, this function should initialize
    /// @p inner_idx and @p outer_ptr. If @p H_values is nonempty, @p inner_idx
    /// and @p outer_ptr can be assumed to be initialized, and this function
    /// should evaluate @p H_values.
    ///
    /// Required for second-order solvers only.
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr,
                     rvec H_Values) const;
    /// **[Optional]**
    /// Function that gets the number of nonzeros of the Hessian of the
    /// Lagrangian.
    ///
    /// Required for second-order solvers only.
    length_t get_hess_L_num_nonzeros() const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the augmented Lagrangian
    /// multiplied by a vector,
    /// @f$ \nabla_{xx}^2L_\Sigma(x, y)\,v @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] Σ
    ///         Penalty weights @f$ \Sigma @f$
    /// @param  [in] scale
    ///         Scale factor for the cost function.
    /// @param  [in] v
    ///         Vector to multiply by @f$ v \in \R^n @f$
    /// @param  [out] Hv
    ///         Hessian-vector product
    ///         @f$ \nabla_{xx}^2 L_\Sigma(x, y)\,v \in \R^{n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const;
    /// **[Optional]**
    /// Function that evaluates the Hessian of the augmented Lagrangian,
    /// @f$ \nabla_{xx}^2L_\Sigma(x, y) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] Σ
    ///         Penalty weights @f$ \Sigma @f$
    /// @param  [in] scale
    ///         Scale factor for the cost function.
    /// @param  [inout] inner_idx
    ///         Inner indices (row indices of nonzeros).
    /// @param  [inout] outer_ptr
    ///         Outer pointers (points to the first nonzero in each column).
    /// @param  [out] H_values
    ///         Nonzero values of the Hessian
    ///         @f$ \nabla_{xx}^2 L_\Sigma(x, y) \in \R^{n\times n} @f$
    /// If @p H_values has size zero, this function should initialize
    /// @p inner_idx and @p outer_ptr. If @p H_values is nonempty, @p inner_idx
    /// and @p outer_ptr can be assumed to be initialized, and this function
    /// should evaluate @p H_values.
    ///
    /// Required for second-order solvers only.
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx,
                     rindexvec outer_ptr, rvec H_values) const;
    /// **[Optional]**
    /// Function that gets the number of nonzeros of the Hessian of the
    /// augmented Lagrangian.
    ///
    /// Required for second-order solvers only.
    length_t get_hess_ψ_num_nonzeros() const;

    /// @}

    /// @name Combined evaluations
    /// @{

    /// **[Optional]**
    /// Evaluate both @f$ f(x) @f$ and its gradient, @f$ \nabla f(x) @f$.
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    /// **[Optional]**
    /// Evaluate both @f$ f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_g(crvec x, rvec g) const;
    /// **[Optional]**
    /// Evaluate @f$ f(x) @f$, its gradient @f$ \nabla f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const;
    /// **[Optional]**
    /// Evaluate both @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$.
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const;
    /// **[Optional]**
    /// Evaluate the gradient of the Lagrangian
    /// @f$ \nabla_x L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;

    /// @}

    /// @name Augmented Lagrangian
    /// @{

    /// **[Optional]**
    /// Calculate both ψ(x) and the vector ŷ that can later be used to compute
    /// ∇ψ.
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    ///   \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \hat y = \Sigma\, \left(g(x) + \Sigma^{-1}y - \Pi_D\left(g(x)
    ///   + \Sigma^{-1}y\right)\right) @f]
    real_t eval_ψ(crvec x, ///< [in]  Decision variable @f$ x @f$
                  crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                  crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                  rvec ŷ   ///< [out] @f$ \hat y @f$
    ) const;
    /// **[Optional]**
    /// Calculate ∇ψ(x) using ŷ.
    void eval_grad_ψ_from_ŷ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                            crvec ŷ,     ///< [in]  @f$ \hat y @f$
                            rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                            rvec work_n  ///<       Dimension @f$ n @f$
    ) const;
    /// **[Optional]**
    /// Calculate the gradient ∇ψ(x).
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    void eval_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                     crvec y,     ///< [in]  Lagrange multipliers @f$ y @f$
                     crvec Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                     rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                     rvec work_n, ///<       Dimension @f$ n @f$
                     rvec work_m  ///<       Dimension @f$ m @f$
    ) const;
    /// **[Optional]**
    /// Calculate both ψ(x) and its gradient ∇ψ(x).
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    /// \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    real_t eval_ψ_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                         crvec y,     ///< [in]  Lagrange multipliers @f$ y @f$
                         crvec Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                         rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                         rvec work_n, ///<       Dimension @f$ n @f$
                         rvec work_m  ///<       Dimension @f$ m @f$
    ) const;

    /// @}

    /// @name Checks
    /// @{

    /// **[Required]**
    /// Check that the problem formulation is well-defined, the dimensions match,
    /// etc. Throws an exception if this is not the case.
    void check() const;

    /// @}

    /// @name Querying specialized implementations
    /// @{

    /// Returns true if the problem provides an implementation of
    /// @ref eval_jac_g.
    bool provides_eval_jac_g() const { return vtable.eval_jac_g != vtable.default_eval_jac_g; }
    /// Returns true if the problem provides an implementation of
    /// @ref get_jac_g_num_nonzeros.
    bool provides_get_jac_g_num_nonzeros() const {
        return vtable.get_jac_g_num_nonzeros != vtable.default_get_jac_g_num_nonzeros;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_grad_gi.
    bool provides_eval_grad_gi() const {
        return vtable.eval_grad_gi != vtable.default_eval_grad_gi;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_L_prod.
    bool provides_eval_hess_L_prod() const {
        return vtable.eval_hess_L_prod != vtable.default_eval_hess_L_prod;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_L.
    bool provides_eval_hess_L() const { return vtable.eval_hess_L != vtable.default_eval_hess_L; }
    /// Returns true if the problem provides an implementation of
    /// @ref get_hess_L_num_nonzeros.
    bool provides_get_hess_L_num_nonzeros() const {
        return vtable.get_hess_L_num_nonzeros != vtable.default_get_hess_L_num_nonzeros;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_ψ_prod.
    bool provides_eval_hess_ψ_prod() const {
        return vtable.eval_hess_ψ_prod != vtable.default_eval_hess_ψ_prod;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref eval_hess_ψ.
    bool provides_eval_hess_ψ() const { return vtable.eval_hess_ψ != vtable.default_eval_hess_ψ; }
    /// Returns true if the problem provides an implementation of
    /// @ref get_hess_ψ_num_nonzeros.
    bool provides_get_hess_ψ_num_nonzeros() const {
        return vtable.get_hess_ψ_num_nonzeros != vtable.default_get_hess_ψ_num_nonzeros;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f, false if it uses the default implementation.
    bool provides_eval_f_grad_f() const {
        return vtable.eval_f_grad_f != vtable.default_eval_f_grad_f;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_g, false if it uses the default implementation.
    bool provides_eval_f_g() const { return vtable.eval_f_g != vtable.default_eval_f_g; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f_g, false if it uses the default implementation.
    bool provides_eval_f_grad_f_g() const {
        return vtable.eval_f_grad_f_g != vtable.default_eval_f_grad_f_g;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_f_grad_g_prod, false if it uses the default implementation.
    bool provides_eval_grad_f_grad_g_prod() const {
        return vtable.eval_grad_f_grad_g_prod != vtable.default_eval_grad_f_grad_g_prod;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_L, false if it uses the default implementation.
    bool provides_eval_grad_L() const { return vtable.eval_grad_L != vtable.default_eval_grad_L; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ, false if it uses the default implementation.
    bool provides_eval_ψ() const { return vtable.eval_ψ != vtable.default_eval_ψ; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ_from_ŷ, false if it uses the default implementation.
    bool provides_eval_grad_ψ_from_ŷ() const {
        return vtable.eval_grad_ψ_from_ŷ != vtable.default_eval_grad_ψ_from_ŷ;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ, false if it uses the default implementation.
    bool provides_eval_grad_ψ() const { return vtable.eval_grad_ψ != vtable.default_eval_grad_ψ; }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ_grad_ψ, false if it uses the default implementation.
    bool provides_eval_ψ_grad_ψ() const {
        return vtable.eval_ψ_grad_ψ != vtable.default_eval_ψ_grad_ψ;
    }
    /// Returns true if the problem provides an implementation of
    /// @ref get_box_C.
    bool provides_get_box_C() const { return vtable.get_box_C != vtable.default_get_box_C; }
    /// Returns true if the problem provides an implementation of
    /// @ref get_box_D.
    bool provides_get_box_D() const { return vtable.get_box_D != vtable.default_get_box_D; }

    /// @}

    /// @name Helpers
    /// @{

    /// Given g(x), compute the intermediate results ŷ and dᵀŷ that can later be
    /// used to compute ψ(x) and ∇ψ(x).
    ///
    /// Computes the result using the following algorithm:
    /// @f[ \begin{aligned}
    ///     \zeta &= g(x) + \Sigma^{-1} y \\[]
    ///     d &= \zeta - \Pi_D(\zeta)
    ///        = \operatorname{eval\_proj\_diff\_g}(\zeta, \zeta) \\[]
    ///     \hat y &= \Sigma d \\[]
    /// \end{aligned} @f]
    /// @see @ref page_math
    ///
    /// @param[inout]   g_ŷ
    ///                 Input @f$ g(x) @f$, outputs @f$ \hat y @f$
    /// @param[in]      y
    ///                 Lagrange multipliers @f$ y @f$
    /// @param[in]      Σ
    ///                 Penalty weights @f$ \Sigma @f$
    /// @return The inner product @f$ d^\top \hat y @f$
    real_t calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ) const;

    /// @}
};

/// @}

#ifndef DOXYGEN
template <class Tref>
explicit TypeErasedProblem(Tref &&d)
    -> TypeErasedProblem<typename std::remove_cvref_t<Tref>::config_t>;

template <class Tref, class Allocator>
explicit TypeErasedProblem(Tref &&d, Allocator alloc)
    -> TypeErasedProblem<typename std::remove_cvref_t<Tref>::config_t, Allocator>;
#endif

template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_n() const -> length_t {
    return vtable.n;
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_m() const -> length_t {
    return vtable.m;
}

template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_proj_diff_g(crvec z, rvec p) const {
    return call(vtable.eval_proj_diff_g, z, p);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_proj_multipliers(rvec y, real_t M,
                                                               index_t penalty_alm_split) const {
    return call(vtable.eval_proj_multipliers, y, M, penalty_alm_split);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ,
                                                             rvec x̂, rvec p) const -> real_t {
    return call(vtable.eval_prox_grad_step, γ, x, grad_ψ, x̂, p);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f(crvec x) const -> real_t {
    return call(vtable.eval_f, x);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_f(crvec x, rvec grad_fx) const {
    return call(vtable.eval_grad_f, x, grad_fx);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_g(crvec x, rvec gx) const {
    return call(vtable.eval_g, x, gx);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const {
    return call(vtable.eval_grad_g_prod, x, y, grad_gxy);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_gi(crvec x, index_t i, rvec grad_gi) const {
    return call(vtable.eval_grad_gi, x, i, grad_gi);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_jac_g(crvec x, rindexvec inner_idx,
                                                    rindexvec outer_ptr, rvec J_values) const {
    return call(vtable.eval_jac_g, x, inner_idx, outer_ptr, J_values);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_jac_g_num_nonzeros() const -> length_t {
    return call(vtable.get_jac_g_num_nonzeros);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v,
                                                          rvec Hv) const {
    return call(vtable.eval_hess_L_prod, x, y, scale, v, Hv);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_L(crvec x, crvec y, real_t scale,
                                                     rindexvec inner_idx, rindexvec outer_ptr,
                                                     rvec H_values) const {
    return call(vtable.eval_hess_L, x, y, scale, inner_idx, outer_ptr, H_values);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_hess_L_num_nonzeros() const -> length_t {
    return call(vtable.get_hess_L_num_nonzeros);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale,
                                                          crvec v, rvec Hv) const {
    return call(vtable.eval_hess_ψ_prod, x, y, Σ, scale, v, Hv);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale,
                                                     rindexvec inner_idx, rindexvec outer_ptr,
                                                     rvec H_values) const {
    return call(vtable.eval_hess_ψ, x, y, Σ, scale, inner_idx, outer_ptr, H_values);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_hess_ψ_num_nonzeros() const -> length_t {
    return call(vtable.get_hess_ψ_num_nonzeros);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_grad_f(crvec x, rvec grad_fx) const -> real_t {
    return call(vtable.eval_f_grad_f, x, grad_fx);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_g(crvec x, rvec g) const -> real_t {
    return call(vtable.eval_f_g, x, g);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const
    -> real_t {
    return call(vtable.eval_f_grad_f_g, x, grad_fx, g);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                                                 rvec grad_gxy) const {
    return call(vtable.eval_grad_f_grad_g_prod, x, y, grad_f, grad_gxy);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_L(crvec x, crvec y, rvec grad_L,
                                                     rvec work_n) const {
    return call(vtable.eval_grad_L, x, y, grad_L, work_n);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const -> real_t {
    return call(vtable.eval_ψ, x, y, Σ, ŷ);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ,
                                                            rvec work_n) const {
    return call(vtable.eval_grad_ψ_from_ŷ, x, ŷ, grad_ψ, work_n);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                                     rvec work_n, rvec work_m) const {
    return call(vtable.eval_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                                       rvec work_n, rvec work_m) const -> real_t {
    return call(vtable.eval_ψ_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ) const -> real_t {
    return call(vtable.calc_ŷ_dᵀŷ, g_ŷ, y, Σ);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_box_C() const -> const Box & {
    return call(vtable.get_box_C);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::get_box_D() const -> const Box & {
    return call(vtable.get_box_D);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::check() const {
    return call(vtable.check);
}

/// @addtogroup grp_Problems
/// @{

/// Problem wrapper that keeps track of the number of evaluations and the run
/// time of each function.
/// You probably want to use @ref problem_with_counters or
/// @ref problem_with_counters_ref instead of instantiating this class directly.
/// @note   The evaluation counters are stored using a `std::shared_pointers`,
///         which means that different copies of a @ref ProblemWithCounters
///         instance all share the same counters. To opt out of this behavior,
///         you can use the @ref decouple_evaluations function.
template <class Problem>
struct ProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);
    using Box = typename TypeErasedProblem<config_t>::Box;

    // clang-format off
    void eval_proj_diff_g(crvec z, rvec p) const { ++evaluations->proj_diff_g; return timed(evaluations->time.proj_diff_g, std::bind(&std::remove_cvref_t<Problem>::eval_proj_diff_g, &problem, z, p)); }
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const { ++evaluations->proj_multipliers; return timed(evaluations->time.proj_multipliers, std::bind(&std::remove_cvref_t<Problem>::eval_proj_multipliers, &problem, y, M, penalty_alm_split)); }
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { ++evaluations->prox_grad_step; return timed(evaluations->time.prox_grad_step, std::bind(&std::remove_cvref_t<Problem>::eval_prox_grad_step, &problem, γ, x, grad_ψ, x̂, p)); }
    real_t eval_f(crvec x) const { ++evaluations->f; return timed(evaluations->time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, x)); }
    void eval_grad_f(crvec x, rvec grad_fx) const { ++evaluations->grad_f; return timed(evaluations->time.grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f, &problem, x, grad_fx)); }
    void eval_g(crvec x, rvec gx) const { ++evaluations->g; return timed(evaluations->time.g, std::bind(&std::remove_cvref_t<Problem>::eval_g, &problem, x, gx)); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { ++evaluations->grad_g_prod; return timed(evaluations->time.grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_g_prod, &problem, x, y, grad_gxy)); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_gi; } { ++evaluations->grad_gi; return timed(evaluations->time.grad_gi, std::bind(&std::remove_cvref_t<Problem>::eval_grad_gi, &problem, x, i, grad_gi)); }
    void eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const requires requires { &std::remove_cvref_t<Problem>::eval_jac_g; } { ++evaluations->jac_g; return timed(evaluations->time.jac_g, std::bind(&std::remove_cvref_t<Problem>::eval_jac_g, &problem, x, inner_idx, outer_ptr, J_values)); }
    length_t get_jac_g_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_jac_g_num_nonzeros; } { return problem.get_jac_g_num_nonzeros(); }
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L_prod; } { ++evaluations->hess_L_prod; return timed(evaluations->time.hess_L_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L_prod, &problem, x, y, scale, v, Hv)); }
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L; } { ++evaluations->hess_L; return timed(evaluations->time.hess_L, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L, &problem, x, y, scale, inner_idx, outer_ptr, H_values)); }
    length_t get_hess_L_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_hess_L_num_nonzeros; } { return problem.get_hess_L_num_nonzeros(); }
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ_prod; } { ++evaluations->hess_ψ_prod; return timed(evaluations->time.hess_ψ_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ_prod, &problem, x, y, Σ, scale, v, Hv)); }
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_ψ; } { ++evaluations->hess_ψ; return timed(evaluations->time.hess_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_hess_ψ, &problem, x, y, Σ, scale, inner_idx, outer_ptr, H_values)); }
    length_t get_hess_ψ_num_nonzeros() const requires requires { &std::remove_cvref_t<Problem>::get_hess_ψ_num_nonzeros; } { return problem.get_hess_ψ_num_nonzeros(); }
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f; } { ++evaluations->f_grad_f; return timed(evaluations->time.f_grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f, &problem, x, grad_fx)); }
    real_t eval_f_g(crvec x, rvec g) const requires requires { &std::remove_cvref_t<Problem>::eval_f_g; } { ++evaluations->f_g; return timed(evaluations->time.f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_g, &problem, x, g)); }
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f_g; } { ++evaluations->f_grad_f_g; return timed(evaluations->time.f_grad_f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f_g, &problem, x, grad_fx, g)); }
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod; } { ++evaluations->grad_f_grad_g_prod; return timed(evaluations->time.grad_f_grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod, &problem, x, y, grad_f, grad_gxy)); }
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_L; } { ++evaluations->grad_L; return timed(evaluations->time.grad_L, std::bind(&std::remove_cvref_t<Problem>::eval_grad_L, &problem, x, y, grad_L, work_n)); }
    real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ; } { ++evaluations->ψ; return timed(evaluations->time.ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ, &problem, x, y, Σ, ŷ)); }
    void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ; } { ++evaluations->grad_ψ_from_ŷ; return timed(evaluations->time.grad_ψ_from_ŷ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ, &problem, x, ŷ, grad_ψ, work_n)); }
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ; } { ++evaluations->grad_ψ; return timed(evaluations->time.grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ_grad_ψ; } { ++evaluations->ψ_grad_ψ; return timed(evaluations->time.ψ_grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    const Box &get_box_C() const requires requires { &std::remove_cvref_t<Problem>::get_box_C; } { return problem.get_box_C(); }
    const Box &get_box_D() const requires requires { &std::remove_cvref_t<Problem>::get_box_D; } { return problem.get_box_D(); }
    void check() const { problem.check(); }

    [[nodiscard]] bool provides_eval_grad_gi() const requires requires (Problem p) { { p.provides_eval_grad_gi() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_gi(); }
    [[nodiscard]] bool provides_eval_jac_g() const requires requires (Problem p) { { p.provides_eval_jac_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_jac_g(); }
    [[nodiscard]] bool provides_get_jac_g_num_nozeros() const requires requires (Problem p) { { p.provides_get_jac_g_num_nozeros() } -> std::convertible_to<bool>; } { return problem.provides_get_jac_g_num_nozeros(); }
    [[nodiscard]] bool provides_eval_hess_L_prod() const requires requires (Problem p) { { p.provides_eval_hess_L_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L_prod(); }
    [[nodiscard]] bool provides_eval_hess_L() const requires requires (Problem p) { { p.provides_eval_hess_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_L(); }
    [[nodiscard]] bool provides_get_hess_L_num_nozeros() const requires requires (Problem p) { { p.provides_get_hess_L_num_nozeros() } -> std::convertible_to<bool>; } { return problem.provides_get_hess_L_num_nozeros(); }
    [[nodiscard]] bool provides_eval_hess_ψ() const requires requires (Problem p) { { p.provides_eval_hess_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_hess_ψ(); }
    [[nodiscard]] bool provides_get_hess_ψ_num_nozeros() const requires requires (Problem p) { { p.provides_get_hess_ψ_num_nozeros() } -> std::convertible_to<bool>; } { return problem.provides_get_hess_ψ_num_nozeros(); }
    [[nodiscard]] bool provides_eval_f_grad_f() const requires requires (Problem p) { { p.provides_eval_f_grad_f() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_grad_f(); }
    [[nodiscard]] bool provides_eval_f_g() const requires requires (Problem p) { { p.provides_eval_f_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_g(); }
    [[nodiscard]] bool provides_eval_f_grad_f_g() const requires requires (Problem p) { { p.provides_eval_f_grad_f_g() } -> std::convertible_to<bool>; } { return problem.provides_eval_f_grad_f_g(); }
    [[nodiscard]] bool provides_eval_grad_f_grad_g_prod() const requires requires (Problem p) { { p.provides_eval_grad_f_grad_g_prod() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_f_grad_g_prod(); }
    [[nodiscard]] bool provides_eval_grad_L() const requires requires (Problem p) { { p.provides_eval_grad_L() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_L(); }
    [[nodiscard]] bool provides_eval_ψ() const requires requires (Problem p) { { p.provides_eval_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ(); }
    [[nodiscard]] bool provides_eval_grad_ψ_from_ŷ() const requires requires (Problem p) { { p.provides_eval_grad_ψ_from_ŷ() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_ψ_from_ŷ(); }
    [[nodiscard]] bool provides_eval_grad_ψ() const requires requires (Problem p) { { p.provides_eval_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_grad_ψ(); }
    [[nodiscard]] bool provides_eval_ψ_grad_ψ() const requires requires (Problem p) { { p.provides_eval_ψ_grad_ψ() } -> std::convertible_to<bool>; } { return problem.provides_eval_ψ_grad_ψ(); }
    [[nodiscard]] bool provides_get_box_C() const requires requires (Problem p) { { p.provides_get_box_C() } -> std::convertible_to<bool>; } { return problem.provides_get_box_C(); }
    [[nodiscard]] bool provides_get_box_D() const requires requires (Problem p) { { p.provides_get_box_D() } -> std::convertible_to<bool>; } { return problem.provides_get_box_D(); }
    // clang-format on

    [[nodiscard]] length_t get_n() const { return problem.get_n(); }
    [[nodiscard]] length_t get_m() const { return problem.get_m(); }

    std::shared_ptr<EvalCounter> evaluations = std::make_shared<EvalCounter>();
    Problem problem;

    template <class P>
    explicit ProblemWithCounters(P &&problem)
        requires std::is_same_v<std::remove_cvref_t<P>, std::remove_cvref_t<Problem>>
        : problem{std::forward<P>(problem)} {}
    template <class... Args>
    explicit ProblemWithCounters(std::in_place_t, Args &&...args)
        requires(!std::is_lvalue_reference_v<Problem>)
        : problem{std::forward<Args>(args)...} {}

    /// Reset all evaluation counters and timers to zero. Affects all instances
    /// that share the same evaluations. If you only want to reset the counters
    /// of this instance, use @ref decouple_evaluations first.
    void reset_evaluations() { evaluations.reset(); }
    /// Give this instance its own evaluation counters and timers, decoupling
    /// it from any other instances they might have previously been shared with.
    /// The evaluation counters and timers are preserved (a copy is made).
    void decouple_evaluations() { evaluations = std::make_shared<EvalCounter>(*evaluations); }

  private:
    template <class TimeT, class FunT>
    static decltype(auto) timed(TimeT &time, FunT &&f) {
        detail::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper has its own copy of the given problem. Making copies of the
/// wrapper also copies the underlying problem, but does not copy the evaluation
/// counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto problem_with_counters(Problem &&p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ProblemWithCounters<Prob>;
    return ProbWithCnt{std::forward<Problem>(p)};
}

/// Wraps the given problem into a @ref ProblemWithCounters and keeps track of
/// how many times each function is called, and how long these calls took.
/// The wrapper keeps only a reference to the given problem, it is the
/// responsibility of the caller to make sure that the wrapper does not outlive
/// the original problem. Making copies of the wrapper does not copy the
/// evaluation counters, all copies share the same counters.
template <class Problem>
[[nodiscard]] auto problem_with_counters_ref(Problem &p) {
    using Prob        = std::remove_cvref_t<Problem>;
    using ProbWithCnt = ProblemWithCounters<const Prob &>;
    return ProbWithCnt{p};
}

template <Config Conf>
class BoxConstrProblem {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    /// Number of decision variables, dimension of x
    length_t n;
    /// Number of constraints, dimension of g(x) and z
    length_t m;

    BoxConstrProblem(length_t n, ///< Number of decision variables
                     length_t m) ///< Number of constraints
        : n{n}, m{m} {}

    BoxConstrProblem(Box C, Box D)
        : n{C.lowerbound.size()}, m{D.lowerbound.size()}, C{std::move(C)}, D{std::move(D)} {}

    BoxConstrProblem(const BoxConstrProblem &)                = default;
    BoxConstrProblem &operator=(const BoxConstrProblem &)     = default;
    BoxConstrProblem(BoxConstrProblem &&) noexcept            = default;
    BoxConstrProblem &operator=(BoxConstrProblem &&) noexcept = default;

    /// Constraints of the decision variables, @f$ x \in C @f$
    Box C{this->n};
    /// Other constraints, @f$ g(x) \in D @f$
    Box D{this->m};

    /// Number of decision variables, @ref n
    length_t get_n() const { return n; }
    /// Number of constraints, @ref m
    length_t get_m() const { return m; }

    /// @f$ \hat x = \Pi_C(x - \gamma\nabla\psi(x)) @f$
    /// @f$ p = \hat x - x @f$
    static real_t eval_proj_grad_step_box(const Box &C, real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                                          rvec p) {
        using binary_real_f = real_t (*)(real_t, real_t);
        p                   = (-γ * grad_ψ)
                .binaryExpr(C.lowerbound - x, binary_real_f(std::fmax))
                .binaryExpr(C.upperbound - x, binary_real_f(std::fmin));
        x̂ = x + p;
        return real_t(0);
    }

    /// @see @ref TypeErasedProblem::eval_prox_grad_step
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const {
        return eval_proj_grad_step_box(C, γ, x, grad_ψ, x̂, p);
    }

    /// @see @ref TypeErasedProblem::eval_proj_diff_g
    void eval_proj_diff_g(crvec z, rvec p) const { p = alpaqa::projecting_difference(z, D); }

    static void eval_proj_multipliers_box(const Box &D, rvec y, real_t M,
                                          index_t penalty_alm_split) {
        auto max_lb = [M](real_t y, real_t z_lb) {
            real_t y_lb = z_lb == -alpaqa::inf<config_t> ? 0 : -M;
            return std::max(y, y_lb);
        };
        auto min_ub = [M](real_t y, real_t z_ub) {
            real_t y_ub = z_ub == alpaqa::inf<config_t> ? 0 : M;
            return std::min(y, y_ub);
        };
        auto num_alm    = y.size() - penalty_alm_split;
        auto &&y_alm    = y.bottomRows(num_alm);
        auto &&z_alm_lb = D.lowerbound.bottomRows(num_alm);
        auto &&z_alm_ub = D.upperbound.bottomRows(num_alm);
        y_alm           = y_alm.binaryExpr(z_alm_lb, max_lb).binaryExpr(z_alm_ub, min_ub);
    }

    /// @see @ref TypeErasedProblem::eval_proj_multipliers
    void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const {
        eval_proj_multipliers_box(D, y, M, penalty_alm_split);
    }

    /// @see @ref TypeErasedProblem::get_box_C
    const Box &get_box_C() const { return C; }
    /// @see @ref TypeErasedProblem::get_box_D
    const Box &get_box_D() const { return D; }

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
    }
};

/// Problem class that allows specifying the basic functions as C++
/// `std::function`s.
/// @ingroup grp_Problems
template <Config Conf = DefaultConfig>
class FunctionalProblem : public BoxConstrProblem<Conf> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using BoxConstrProblem<Conf>::BoxConstrProblem;

    std::function<real_t(crvec)> f;
    std::function<void(crvec, rvec)> grad_f;
    std::function<void(crvec, rvec)> g;
    std::function<void(crvec, crvec, rvec)> grad_g_prod;
    std::function<void(crvec, rmat)> jac_g;
    std::function<void(crvec, crvec, real_t, crvec, rvec)> hess_L_prod;
    std::function<void(crvec, crvec, real_t, rmat)> hess_L;
    std::function<void(crvec, crvec, crvec, real_t, crvec, rvec)> hess_ψ_prod;
    std::function<void(crvec, crvec, crvec, real_t, rmat)> hess_ψ;

    // clang-format off
    real_t eval_f(crvec x) const { ScopedMallocAllower ma; return f(x); }
    void eval_grad_f(crvec x, rvec grad_fx) const { ScopedMallocAllower ma; grad_f(x, grad_fx); }
    void eval_g(crvec x, rvec gx) const { ScopedMallocAllower ma; g(x, gx); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { ScopedMallocAllower ma; grad_g_prod(x, y, grad_gxy); }
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const { ScopedMallocAllower ma; hess_L_prod(x, y, scale, v, Hv); }
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const { ScopedMallocAllower ma; hess_ψ_prod(x, y, Σ, scale, v, Hv); }
    // clang-format on
    void eval_jac_g(crvec x, rindexvec, rindexvec, rvec J_values) const {
        ScopedMallocAllower ma;
        if (J_values.size() > 0)
            jac_g(x, J_values.reshaped(this->m, this->n));
    }
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec, rindexvec, rvec H_values) const {
        ScopedMallocAllower ma;
        if (H_values.size() > 0)
            hess_L(x, y, scale, H_values.reshaped(this->n, this->n));
    }
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec, rindexvec,
                     rvec H_values) const {
        ScopedMallocAllower ma;
        if (H_values.size() > 0)
            hess_ψ(x, y, Σ, scale, H_values.reshaped(this->n, this->n));
    }

    /// @see @ref TypeErasedProblem::provides_eval_jac_g
    bool provides_eval_jac_g() const { return bool{jac_g}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L_prod
    bool provides_eval_hess_L_prod() const { return bool{hess_L_prod}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L
    bool provides_eval_hess_L() const { return bool{hess_L}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ_prod
    bool provides_eval_hess_ψ_prod() const { return bool{hess_ψ_prod}; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ
    bool provides_eval_hess_ψ() const { return bool{hess_ψ}; }

    FunctionalProblem(const FunctionalProblem &)                = default;
    FunctionalProblem &operator=(const FunctionalProblem &)     = default;
    FunctionalProblem(FunctionalProblem &&) noexcept            = default;
    FunctionalProblem &operator=(FunctionalProblem &&) noexcept = default;
};

template <Config Conf>
void print_provided_functions(std::ostream &os, const TypeErasedProblem<Conf> &problem) {
    os << "           eval_grad_gi: " << problem.provides_eval_grad_gi() << '\n'
       << "                  jac_g: " << problem.provides_eval_jac_g() << '\n'
       << "       eval_hess_L_prod: " << problem.provides_eval_hess_L_prod() << '\n'
       << "            eval_hess_L: " << problem.provides_eval_hess_L() << '\n'
       << "       eval_hess_ψ_prod: " << problem.provides_eval_hess_ψ_prod() << '\n'
       << "            eval_hess_ψ: " << problem.provides_eval_hess_ψ() << '\n'
       << "          eval_f_grad_f: " << problem.provides_eval_f_grad_f() << '\n'
       << "               eval_f_g: " << problem.provides_eval_f_g() << '\n'
       << "        eval_f_grad_f_g: " << problem.provides_eval_f_grad_f_g() << '\n'
       << "eval_grad_f_grad_g_prod: " << problem.provides_eval_grad_f_grad_g_prod() << '\n'
       << "            eval_grad_L: " << problem.provides_eval_grad_L() << '\n'
       << "                 eval_ψ: " << problem.provides_eval_ψ() << '\n'
       << "     eval_grad_ψ_from_ŷ: " << problem.provides_eval_grad_ψ_from_ŷ() << '\n'
       << "            eval_grad_ψ: " << problem.provides_eval_grad_ψ() << '\n'
       << "          eval_ψ_grad_ψ: " << problem.provides_eval_ψ_grad_ψ() << '\n'
       << "              get_box_C: " << problem.provides_get_box_C() << '\n'
       << "              get_box_D: " << problem.provides_get_box_D() << '\n';
}

/// @}

} // namespace alpaqa