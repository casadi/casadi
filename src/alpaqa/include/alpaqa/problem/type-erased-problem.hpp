#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/export.hpp>
#include <alpaqa/problem/box.hpp>
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

/// Struct containing function pointers to all problem functions (like the
/// objective and constraint functions, with their derivatives, and more).
/// Some default implementations are available.
/// Internal struct, it is used by @ref TypeErasedProblem.
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
    required_const_function_t<void(crvec z, rvec e)>
        eval_proj_diff_g;
    required_const_function_t<void(rvec y, real_t M)>
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
    optional_const_function_t<index_t(real_t γ, crvec x, crvec grad_ψ, rindexvec J)>
        eval_inactive_indices_res_lna = default_eval_inactive_indices_res_lna;

    // Second order
    optional_const_function_t<void(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values)>
        eval_jac_g = default_eval_jac_g;
    optional_const_function_t<length_t()>
        get_jac_g_num_nonzeros = default_get_jac_g_num_nonzeros;
    optional_const_function_t<void(crvec x, index_t i, rvec grad_gi)>
        eval_grad_gi = default_eval_grad_gi;
    optional_const_function_t<void(crvec x, crvec y, real_t scale, crvec v, rvec Hv)>
        eval_hess_L_prod = default_eval_hess_L_prod;
    optional_const_function_t<void(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values)>
        eval_hess_L = default_eval_hess_L;
    optional_const_function_t<length_t()>
        get_hess_L_num_nonzeros = default_get_hess_L_num_nonzeros;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv)>
        eval_hess_ψ_prod = default_eval_hess_ψ_prod;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values)>
        eval_hess_ψ = default_eval_hess_ψ;
    optional_const_function_t<length_t()>
        get_hess_ψ_num_nonzeros = default_get_hess_ψ_num_nonzeros;

    // Combined evaluations
    optional_const_function_t<real_t(crvec x, rvec grad_fx)>
        eval_f_grad_f = default_eval_f_grad_f;
    optional_const_function_t<real_t(crvec x, rvec g)>
        eval_f_g = default_eval_f_g;
    optional_const_function_t<void(crvec x, crvec y, rvec grad_f, rvec grad_gxy)>
        eval_grad_f_grad_g_prod = default_eval_grad_f_grad_g_prod;

    // Lagrangian and augmented lagrangian evaluations
    optional_const_function_t<void(crvec x, crvec y, rvec grad_L, rvec work_n)>
        eval_grad_L = default_eval_grad_L;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec ŷ)>
        eval_ψ = default_eval_ψ;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_grad_ψ = default_eval_grad_ψ;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_ψ_grad_ψ = default_eval_ψ_grad_ψ;

    // Constraint sets
    optional_const_function_t<const Box &()>
        get_box_C = default_get_box_C;
    optional_const_function_t<const Box &()>
        get_box_D = default_get_box_D;

    // Check
    optional_const_function_t<void()>
        check = default_check;

    // clang-format on

    ALPAQA_EXPORT static real_t calc_ŷ_dᵀŷ(const void *self, rvec g_ŷ, crvec y, crvec Σ,
                                           const ProblemVTable &vtable);
    ALPAQA_EXPORT static index_t default_eval_inactive_indices_res_lna(const void *, real_t, crvec,
                                                                       crvec, rindexvec,
                                                                       const ProblemVTable &);
    ALPAQA_EXPORT static void default_eval_jac_g(const void *, crvec, rindexvec, rindexvec, rvec,
                                                 const ProblemVTable &);
    ALPAQA_EXPORT static length_t default_get_jac_g_num_nonzeros(const void *,
                                                                 const ProblemVTable &);
    ALPAQA_EXPORT static void default_eval_grad_gi(const void *, crvec, index_t, rvec,
                                                   const ProblemVTable &);
    ALPAQA_EXPORT static void default_eval_hess_L_prod(const void *, crvec, crvec, real_t, crvec,
                                                       rvec, const ProblemVTable &);
    ALPAQA_EXPORT static void default_eval_hess_L(const void *, crvec, crvec, real_t, rindexvec,
                                                  rindexvec, rvec, const ProblemVTable &);
    ALPAQA_EXPORT static length_t default_get_hess_L_num_nonzeros(const void *,
                                                                  const ProblemVTable &);
    ALPAQA_EXPORT static void default_eval_hess_ψ_prod(const void *self, crvec x, crvec y, crvec,
                                                       real_t scale, crvec v, rvec Hv,
                                                       const ProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_hess_ψ(const void *self, crvec x, crvec y, crvec,
                                                  real_t scale, rindexvec inner_idx,
                                                  rindexvec outer_ptr, rvec H_values,
                                                  const ProblemVTable &vtable);
    ALPAQA_EXPORT static length_t default_get_hess_ψ_num_nonzeros(const void *,
                                                                  const ProblemVTable &);
    ALPAQA_EXPORT static real_t default_eval_f_grad_f(const void *self, crvec x, rvec grad_fx,
                                                      const ProblemVTable &vtable);
    ALPAQA_EXPORT static real_t default_eval_f_g(const void *self, crvec x, rvec g,
                                                 const ProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_grad_f_grad_g_prod(const void *self, crvec x, crvec y,
                                                              rvec grad_f, rvec grad_gxy,
                                                              const ProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_grad_L(const void *self, crvec x, crvec y, rvec grad_L,
                                                  rvec work_n, const ProblemVTable &vtable);
    ALPAQA_EXPORT static real_t default_eval_ψ(const void *self, crvec x, crvec y, crvec Σ, rvec ŷ,
                                               const ProblemVTable &vtable);
    ALPAQA_EXPORT static void default_eval_grad_ψ(const void *self, crvec x, crvec y, crvec Σ,
                                                  rvec grad_ψ, rvec work_n, rvec work_m,
                                                  const ProblemVTable &vtable);
    ALPAQA_EXPORT static real_t default_eval_ψ_grad_ψ(const void *self, crvec x, crvec y, crvec Σ,
                                                      rvec grad_ψ, rvec work_n, rvec work_m,
                                                      const ProblemVTable &vtable);
    ALPAQA_EXPORT static const Box &default_get_box_C(const void *, const ProblemVTable &);
    ALPAQA_EXPORT static const Box &default_get_box_D(const void *, const ProblemVTable &);
    ALPAQA_EXPORT static void default_check(const void *, const ProblemVTable &);

    length_t n, m;

    template <class P>
    ProblemVTable(std::in_place_t, P &p) : util::BasicVTable{std::in_place, p} {
        auto &vtable = *this;

        // Initialize all methods

        // Required
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_diff_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_multipliers);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_prox_grad_step);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_g_prod);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_inactive_indices_res_lna, p);
        // Second order
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_jac_g, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_jac_g_num_nonzeros, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_gi, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L_prod, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_hess_L_num_nonzeros, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ_prod, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_ψ, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_hess_ψ_num_nonzeros, p);
        // Combined evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_g, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_f_grad_g_prod, p);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_L, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ_grad_ψ, p);
        // Constraint set
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_C, p);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, get_box_D, p);
        // Check
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, check, p);

        // Dimensions
        vtable.n = p.get_n();
        vtable.m = p.get_m();
    }
    ProblemVTable() = default;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ProblemVTable, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ProblemVTable, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ProblemVTable, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ProblemVTable, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(struct, ProblemVTable, EigenConfigq);
#endif

/// @addtogroup grp_Problems
/// @{

/// The main polymorphic minimization problem interface.
///
/// This class wraps the actual problem implementation class, filling in the
/// missing member functions with sensible defaults, and providing a uniform
/// interface that is used by the solvers.
///
/// The problem implementations do not inherit from an abstract base class.
/// Instead, [structural typing](https://en.wikipedia.org/wiki/Structural_type_system)
/// is used. The @ref ProblemVTable constructor uses reflection to discover
/// which member functions are provided by the problem implementation. See
/// @ref page_problem_formulations for more information, and
/// @ref C++/CustomCppProblem/main.cpp for an example.
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
    /// @param  [out] e
    ///         The difference relative to its projection,
    ///         @f$ e = z - \Pi_D(z) \in \R^m @f$
    /// @note   @p z and @p e can refer to the same vector.
    void eval_proj_diff_g(crvec z, rvec e) const;
    /// **[Required]**
    /// Function that projects the Lagrange multipliers for ALM.
    /// @param  [inout] y
    ///         Multipliers, @f$ y \leftarrow \Pi_Y(y) \in \R^m @f$
    /// @param  [in] M
    ///         The radius/size of the set @f$ Y @f$.
    ///         See @ref ALMParams::max_multiplier.
    void eval_proj_multipliers(rvec y, real_t M) const;
    /// **[Required]**
    /// Function that computes a proximal gradient step.
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
    /// **[Optional]**
    /// Function that computes the inactive indices @f$ \mathcal J(x) @f$ for
    /// the evaluation of the linear Newton approximation of the residual, as in
    /// @cite pas2022alpaqa.
    /// @param  [in] γ
    ///         Step size, @f$ \gamma \in \R_{>0} @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] grad_ψ
    ///         Gradient of the subproblem cost, @f$ \nabla\psi(x) \in \R^n @f$
    /// @param  [out] J
    ///         The indices of the components of @f$ x @f$ that are in the
    ///         index set @f$ \mathcal J(x) @f$. In ascending order, at most n.
    /// @return The number of inactive constraints, @f$ \# \mathcal J(x) @f$.
    ///
    /// For example, in the case of box constraints, we have
    /// @f[ \mathcal J(x) \defeq \defset{i \in \N_{[0, n-1]}}{\underline x_i
    /// \lt x_i - \gamma\nabla_{\!x_i}\psi(x) \lt \overline x_i}. @f]
    index_t eval_inactive_indices_res_lna(real_t γ, crvec x, crvec grad_ψ, rindexvec J) const;

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
    /// Function that gets the number of nonzeros of the sparse Jacobian of the
    /// constraints. Should return -1 for a dense Jacobian.
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
                     rvec H_values) const;
    /// **[Optional]**
    /// Function that gets the number of nonzeros of the sparse Hessian of the
    /// Lagrangian. Should return -1 for a dense Hessian.
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
    /// @default_impl   ProblemVTable::default_eval_f_grad_f
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    /// **[Optional]**
    /// Evaluate both @f$ f(x) @f$ and @f$ g(x) @f$.
    /// @default_impl   ProblemVTable::default_eval_f_g
    real_t eval_f_g(crvec x, rvec g) const;
    /// **[Optional]**
    /// Evaluate both @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$.
    /// @default_impl   ProblemVTable::default_eval_grad_f_grad_g_prod
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const;
    /// **[Optional]**
    /// Evaluate the gradient of the Lagrangian
    /// @f$ \nabla_x L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
    /// @default_impl   ProblemVTable::default_eval_grad_L
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
    /// @default_impl   ProblemVTable::default_eval_ψ
    real_t eval_ψ(crvec x, ///< [in]  Decision variable @f$ x @f$
                  crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                  crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                  rvec ŷ   ///< [out] @f$ \hat y @f$
    ) const;
    /// **[Optional]**
    /// Calculate the gradient ∇ψ(x).
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    /// @default_impl   ProblemVTable::default_eval_grad_ψ
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
    /// @default_impl   ProblemVTable::default_eval_ψ_grad_ψ
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

    /// **[Optional]**
    /// Check that the problem formulation is well-defined, the dimensions match,
    /// etc. Throws an exception if this is not the case.
    void check() const;

    /// @}

    /// @name Querying specialized implementations
    /// @{

    /// Returns true if the problem provides an implementation of
    /// @ref eval_inactive_indices_res_lna.
    bool provides_eval_inactive_indices_res_lna() const {
        return vtable.eval_inactive_indices_res_lna != vtable.default_eval_inactive_indices_res_lna;
    }
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
    /// Returns true if the problem provides an implementation of @ref check.
    bool provides_check() const { return vtable.check != vtable.default_check; }

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
void TypeErasedProblem<Conf, Allocator>::eval_proj_diff_g(crvec z, rvec e) const {
    return call(vtable.eval_proj_diff_g, z, e);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_proj_multipliers(rvec y, real_t M) const {
    return call(vtable.eval_proj_multipliers, y, M);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ,
                                                             rvec x̂, rvec p) const -> real_t {
    return call(vtable.eval_prox_grad_step, γ, x, grad_ψ, x̂, p);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_inactive_indices_res_lna(real_t γ, crvec x,
                                                                       crvec grad_ψ,
                                                                       rindexvec J) const
    -> index_t {
    return call(vtable.eval_inactive_indices_res_lna, γ, x, grad_ψ, J);
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

template <Config Conf>
void print_provided_functions(std::ostream &os, const TypeErasedProblem<Conf> &problem) {
    os << "inactive_indices_res_lna: " << problem.provides_eval_inactive_indices_res_lna() << '\n'
       << "                 grad_gi: " << problem.provides_eval_grad_gi() << '\n'
       << "                   jac_g: " << problem.provides_eval_jac_g() << '\n'
       << "             hess_L_prod: " << problem.provides_eval_hess_L_prod() << '\n'
       << "                  hess_L: " << problem.provides_eval_hess_L() << '\n'
       << "             hess_ψ_prod: " << problem.provides_eval_hess_ψ_prod() << '\n'
       << "                  hess_ψ: " << problem.provides_eval_hess_ψ() << '\n'
       << "                f_grad_f: " << problem.provides_eval_f_grad_f() << '\n'
       << "                     f_g: " << problem.provides_eval_f_g() << '\n'
       << "      grad_f_grad_g_prod: " << problem.provides_eval_grad_f_grad_g_prod() << '\n'
       << "                  grad_L: " << problem.provides_eval_grad_L() << '\n'
       << "                       ψ: " << problem.provides_eval_ψ() << '\n'
       << "                  grad_ψ: " << problem.provides_eval_grad_ψ() << '\n'
       << "                ψ_grad_ψ: " << problem.provides_eval_ψ_grad_ψ() << '\n'
       << "               get_box_C: " << problem.provides_get_box_C() << '\n'
       << "               get_box_D: " << problem.provides_get_box_D() << '\n'
       << "                   check: " << problem.provides_check() << '\n';
}

/// @}

} // namespace alpaqa