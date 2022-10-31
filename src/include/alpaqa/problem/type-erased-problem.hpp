#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/util/type-erasure.hpp>
#include "alpaqa/util/not-implemented.hpp"
#include <chrono>
#include <stdexcept>
#include <type_traits>

namespace alpaqa {

#define ALPAQA_TE_REQUIRED_METHOD(vtable, type, member)                        \
    vtable.member = util::type_erased_wrapped<&type::member>()
#define ALPAQA_TE_OPTIONAL_METHOD(vtable, type, member)                        \
    if constexpr (requires { &type::member; })                                 \
    vtable.member = util::type_erased_wrapped<                                 \
        &type::member, const std::remove_cvref_t<decltype(vtable)> &>()
#define ALPAQA_TE_DISABLED_METHOD(vtable, member, instance_p)                  \
    if constexpr (requires { (instance_p)->provides_##member(); })             \
        if (!(instance_p)->provides_##member())                                \
    vtable.member = nullptr
#define ALPAQA_TE_DEFAULT_METHOD(vtable, member, default_)                     \
    if (vtable.member == nullptr)                                              \
    vtable.member = default_

template <Config Conf>
struct ProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);

    template <class F>
    using optional_function_t =
        util::BasicVTable::optional_function_t<F, ProblemVTable>;
    template <class F>
    using optional_const_function_t =
        util::BasicVTable::optional_const_function_t<F, ProblemVTable>;

    // clang-format off

    // Required
    required_const_function_t<void(crvec z, rvec p)>
        eval_proj_diff_g;
    required_const_function_t<void(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p)>
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
    optional_const_function_t<void(crvec x, index_t i, rvec grad_gi)>
        eval_grad_gi = nullptr;
    optional_const_function_t<void(crvec x, crvec y, crvec v, rvec Hv)>
        eval_hess_L_prod = nullptr;
    optional_const_function_t<void(crvec x, crvec y, rmat H)>
        eval_hess_L = nullptr;

    // Combined evaluations
    optional_const_function_t<real_t(crvec x, rvec grad_fx)>
        eval_f_grad_f = nullptr;
    optional_const_function_t<real_t(crvec x, rvec g)>
        eval_f_g = nullptr;
    optional_const_function_t<real_t(crvec x, rvec grad_fx, rvec g)>
        eval_f_grad_f_g = nullptr;
    optional_const_function_t<void(crvec x, crvec y, rvec grad_f, rvec grad_gxy)>
        eval_grad_f_grad_g_prod = nullptr;

    // Lagrangian and augmented lagrangian evaluations
    optional_const_function_t<void(crvec x, crvec y, rvec grad_L, rvec work_n)>
        eval_grad_L = nullptr;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec ŷ)>
        eval_ψ_ŷ = nullptr;
    optional_const_function_t<void(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n)>
        eval_grad_ψ_from_ŷ = nullptr;
    optional_const_function_t<void(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_grad_ψ = nullptr;
    optional_const_function_t<real_t(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m)>
        eval_ψ_grad_ψ = nullptr;

    // clang-format on

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

    static void default_eval_grad_gi(const void *, crvec, index_t, rvec,
                                     const ProblemVTable &) {
        throw not_implemented_error("eval_grad_gi");
    }
    static void default_eval_hess_L_prod(const void *, crvec, crvec, crvec,
                                         rvec, const ProblemVTable &) {
        throw not_implemented_error("eval_hess_L_prod");
    }
    static void default_eval_hess_L(const void *, crvec x, crvec y, rmat H,
                                    const ProblemVTable &) {
        throw not_implemented_error("eval_hess_L");
    }

    static real_t default_eval_f_grad_f(const void *self, crvec x, rvec grad_fx,
                                        const ProblemVTable &vtable) {
        vtable.eval_grad_f(self, x, grad_fx);
        return vtable.eval_f(self, x);
    }
    static real_t default_eval_f_g(const void *self, crvec x, rvec g,
                                   const ProblemVTable &vtable) {
        vtable.eval_g(self, x, g);
        return vtable.eval_f(self, x);
    }
    static real_t default_eval_f_grad_f_g(const void *self, crvec x,
                                          rvec grad_fx, rvec g,
                                          const ProblemVTable &vtable) {
        vtable.eval_g(self, x, g);
        return vtable.eval_f_grad_f(self, x, grad_fx, vtable);
    }
    static void default_eval_grad_f_grad_g_prod(const void *self, crvec x,
                                                crvec y, rvec grad_f,
                                                rvec grad_gxy,
                                                const ProblemVTable &vtable) {
        vtable.eval_grad_f(self, x, grad_f);
        vtable.eval_grad_g_prod(self, x, y, grad_gxy);
    }
    static void default_eval_grad_L(const void *self, crvec x, crvec y,
                                    rvec grad_L, rvec work_n,
                                    const ProblemVTable &vtable) {
        vtable.eval_grad_f_grad_g_prod(self, x, y, grad_L, work_n, vtable);
        grad_L += work_n;
    }
    static real_t default_eval_ψ_ŷ(const void *self, crvec x, crvec y, crvec Σ,
                                   rvec ŷ, const ProblemVTable &vtable) {
        if (y.size() == 0) /* [[unlikely]] */
            return vtable.eval_f(self, x);

        real_t f   = vtable.eval_f_g(self, x, ŷ, vtable);
        real_t dᵀŷ = calc_ŷ_dᵀŷ(self, ŷ, y, Σ, vtable);
        // ψ(x) = f(x) + ½ dᵀŷ
        real_t ψ = f + real_t(0.5) * dᵀŷ;
        return ψ;
    }
    static void default_eval_grad_ψ_from_ŷ(const void *self, crvec x, crvec ŷ,
                                           rvec grad_ψ, rvec work_n,
                                           const ProblemVTable &vtable) {
        if (ŷ.size() == 0) /* [[unlikely]] */
            vtable.eval_grad_f(self, x, grad_ψ);
        else
            vtable.eval_grad_L(self, x, ŷ, grad_ψ, work_n, vtable);
    }
    static void default_eval_grad_ψ(const void *self, crvec x, crvec y, crvec Σ,
                                    rvec grad_ψ, rvec work_n, rvec work_m,
                                    const ProblemVTable &vtable) {
        if (y.size() == 0) /* [[unlikely]] */ {
            vtable.eval_grad_f(self, x, grad_ψ);
        } else {
            vtable.eval_g(self, x, work_m);
            (void)calc_ŷ_dᵀŷ(self, work_m, y, Σ, vtable);
            vtable.eval_grad_ψ_from_ŷ(self, x, work_m, grad_ψ, work_n, vtable);
        }
    }
    static real_t default_eval_ψ_grad_ψ(const void *self, crvec x, crvec y,
                                        crvec Σ, rvec grad_ψ, rvec work_n,
                                        rvec work_m,
                                        const ProblemVTable &vtable) {
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

    template <class P>
    ProblemVTable(util::VTableTypeTag<P> t) : util::BasicVTable{t} {
        auto &vtable = *this;
        assert(t.t);

        // Initialize all methods

        // Required
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_proj_diff_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_prox_grad_step);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_f);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_g);
        ALPAQA_TE_REQUIRED_METHOD(vtable, P, eval_grad_g_prod);
        // Second order
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_gi);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L_prod);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_hess_L);
        // Combined evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_g);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_f_grad_f_g);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_f_grad_g_prod);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_L);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ_ŷ);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ_from_ŷ);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_grad_ψ);
        ALPAQA_TE_OPTIONAL_METHOD(vtable, P, eval_ψ_grad_ψ);

        // Disable optional methods

        // Second order
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_grad_gi, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_hess_L_prod, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_hess_L, t.t);
        // Combined evaluations
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_f_grad_f, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_f_g, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_f_grad_f_g, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_grad_f_grad_g_prod, t.t);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_grad_L, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_ψ_ŷ, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_grad_ψ_from_ŷ, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_grad_ψ, t.t);
        ALPAQA_TE_DISABLED_METHOD(vtable, eval_ψ_grad_ψ, t.t);

        // Provide defaults

        // Second order
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_grad_gi, default_eval_grad_gi);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_hess_L_prod,
                                 default_eval_hess_L_prod);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_hess_L, default_eval_hess_L);
        // Combined evaluations
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_f_grad_f, default_eval_f_grad_f);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_f_g, default_eval_f_g);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_f_grad_f_g,
                                 default_eval_f_grad_f_g);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_grad_f_grad_g_prod,
                                 default_eval_grad_f_grad_g_prod);
        // Lagrangian and augmented lagrangian evaluations
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_grad_L, default_eval_grad_L);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_ψ_ŷ, default_eval_ψ_ŷ);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_grad_ψ_from_ŷ,
                                 default_eval_grad_ψ_from_ŷ);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_grad_ψ, default_eval_grad_ψ);
        ALPAQA_TE_DEFAULT_METHOD(vtable, eval_ψ_grad_ψ, default_eval_ψ_grad_ψ);
    }
    ProblemVTable() = default;
};

#undef ALPAQA_TE_OPTIONAL_METHOD
#undef ALPAQA_TE_REQUIRED_METHOD
#undef ALPAQA_TE_DISABLED_METHOD
#undef ALPAQA_TE_DEFAULT_METHOD

template <Config Conf     = DefaultConfig,
          class Allocator = std::allocator<std::byte>>
class TypeErasedProblem
    : public util::TypeErased<ProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = ProblemVTable<Conf>;
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
        return TypeErased::template make<TypeErasedProblem, T>(
            std::forward<Args>(args)...);
    }

    /// @name Basic functions
    /// @{

    /// Function that evaluates the difference between the given point @f$ z @f$
    /// and its projection onto the constraint set @f$ D @f$.
    /// @param  [in] z
    ///         Slack variable, @f$ z \in \R^m @f$
    /// @param  [out] p
    ///         The difference relative to its projection,
    ///         @f$ p = z - \Pi_D(z) \in \R^m @f$
    void eval_proj_diff_g(crvec z, rvec p) const;
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
    /// @note   The vector @f$ p @f$ is often used in stopping criteria, so its
    ///         numerical accuracy is more important than that of @f$ \hat x @f$.
    void eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂,
                             rvec p) const;

    /// Function that evaluates the cost, @f$ f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    real_t eval_f(crvec x) const;
    /// Function that evaluates the gradient of the cost, @f$ \nabla f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] grad_fx
    ///         Gradient of cost function @f$ \nabla f(x) \in \R^n @f$
    void eval_grad_f(crvec x, rvec grad_fx) const;
    /// Function that evaluates the constraints, @f$ g(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] gx
    ///         Value of the constraints @f$ g(x) \in \R^m @f$
    void eval_g(crvec x, rvec gx) const;
    /// Function that evaluates the gradient of the constraints times a vector,
    /// @f$ \nabla g(x)\,y @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Vector @f$ y \in \R^m @f$ to multiply the gradient by
    /// @param  [out] grad_gxy
    ///         Gradient of the constraints
    ///         @f$ \nabla g(x)\,y \in \R^n @f$
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;

    /// @}

    /// @name Functions for second-order solvers
    /// @{

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
    /// Function that evaluates the Hessian of the Lagrangian multiplied by a
    /// vector,
    /// @f$ \nabla_{xx}^2L(x, y)\,v @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [in] v
    ///         Vector to multiply by @f$ v \in \R^n @f$
    /// @param  [out] Hv
    ///         Hessian-vector product
    ///         @f$ \nabla_{xx}^2 L(x, y)\,v \in \R^{n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const;
    /// Function that evaluates the Hessian of the Lagrangian,
    /// @f$ \nabla_{xx}^2L(x, y) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Lagrange multipliers @f$ y \in \R^m @f$
    /// @param  [out] H
    ///         Hessian @f$ \nabla_{xx}^2 L(x, y) \in \R^{n\times n} @f$
    ///
    /// Required for second-order solvers only.
    void eval_hess_L(crvec x, crvec y, rmat H) const;

    /// @}

    /// @name Combined evaluations
    /// @{

    /// Evaluate both @f$ f(x) @f$ and its gradient, @f$ \nabla f(x) @f$.
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    /// Evaluate both @f$ f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_g(crvec x, rvec g) const;
    /// Evaluate @f$ f(x) @f$, its gradient @f$ \nabla f(x) @f$ and @f$ g(x) @f$.
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const;
    /// Evaluate both @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$.
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                 rvec grad_gxy) const;
    /// Evaluate the gradient of the Lagrangian
    /// @f$ \nabla_x L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;

    /// @}

    /// @name Augmented Lagrangian
    /// @{

    /// Calculate both ψ(x) and the vector ŷ that can later be used to compute
    /// ∇ψ.
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    ///   \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \hat y = \Sigma\, \left(g(x) + \Sigma^{-1}y - \Pi_D\left(g(x)
    ///   + \Sigma^{-1}y\right)\right) @f]
    real_t eval_ψ_ŷ(crvec x, ///< [in]  Decision variable @f$ x @f$
                    crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                    crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                    rvec ŷ   ///< [out] @f$ \hat y @f$
    ) const;
    /// Calculate ∇ψ(x) using ŷ.
    void eval_grad_ψ_from_ŷ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                            crvec ŷ,     ///< [in]  @f$ \hat y @f$
                            rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                            rvec work_n  ///<       Dimension @f$ n @f$
    ) const;
    /// Calculate the gradient ∇ψ(x).
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    void eval_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                     crvec y,     ///< [in]  Lagrange multipliers @f$ y @f$
                     crvec Σ,     ///< [in]  Penalty weights @f$ \Sigma @f$
                     rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                     rvec work_n, ///<       Dimension @f$ n @f$
                     rvec work_m  ///<       Dimension @f$ m @f$
    ) const;
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

    /// @name Helpers
    /// @{

    /// Given g(x), compute the intermediate results ŷ and dᵀŷ that can later be
    /// used to compute ψ(x) and ∇ψ(x).
    /// @param[inout]   g_ŷ
    ///                 Input @f$ g(x) @f$, outputs @f$ \hat y @f$
    /// @param[in]      y
    ///                 Lagrange multipliers @f$ y @f$
    /// @param[in]      Σ
    ///                 Penalty weights @f$ \Sigma @f$
    /// @return The inner product @f$ d^\top \hat y @f$
    real_t calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y, crvec Σ) const;

    /// @}

    /// @name Querying specialized implementations
    /// @{

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
    bool provides_eval_hess_L() const {
        return vtable.eval_hess_L != vtable.default_eval_hess_L;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f, false if it uses the default implementation.
    bool provides_eval_f_grad_f() const {
        return vtable.eval_f_grad_f != vtable.default_eval_f_grad_f;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_g, false if it uses the default implementation.
    bool provides_eval_f_g() const {
        return vtable.eval_f_g != vtable.default_eval_f_g;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_f_grad_f_g, false if it uses the default implementation.
    bool provides_eval_f_grad_f_g() const {
        return vtable.eval_f_grad_f_g != vtable.default_eval_f_grad_f_g;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_f_grad_g_prod, false if it uses the default implementation.
    bool provides_eval_grad_f_grad_g_prod() const {
        return vtable.eval_grad_f_grad_g_prod !=
               vtable.default_eval_grad_f_grad_g_prod;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_L, false if it uses the default implementation.
    bool provides_eval_grad_L() const {
        return vtable.eval_grad_L != vtable.default_eval_grad_L;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ_ŷ, false if it uses the default implementation.
    bool provides_eval_ψ_ŷ() const {
        return vtable.eval_ψ_ŷ != vtable.default_eval_ψ_ŷ;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ_from_ŷ, false if it uses the default implementation.
    bool provides_eval_grad_ψ_from_ŷ() const {
        return vtable.eval_grad_ψ_from_ŷ != vtable.default_eval_grad_ψ_from_ŷ;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_grad_ψ, false if it uses the default implementation.
    bool provides_eval_grad_ψ() const {
        return vtable.eval_grad_ψ != vtable.default_eval_grad_ψ;
    }
    /// Returns true if the problem provides a specialized implementation of
    /// @ref eval_ψ_grad_ψ, false if it uses the default implementation.
    bool provides_eval_ψ_grad_ψ() const {
        return vtable.eval_ψ_grad_ψ != vtable.default_eval_ψ_grad_ψ;
    }

    /// @}
};

template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_proj_diff_g(crvec z,
                                                          rvec p) const {
    return call(vtable.eval_proj_diff_g, z, p);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_prox_grad_step(real_t γ, crvec x,
                                                             crvec grad_ψ,
                                                             rvec x̂,
                                                             rvec p) const {
    return call(vtable.eval_prox_grad_step, γ, x, grad_ψ, x̂, p);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f(crvec x) const -> real_t {
    return call(vtable.eval_f, x);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_f(crvec x,
                                                     rvec grad_fx) const {
    return call(vtable.eval_grad_f, x, grad_fx);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_g(crvec x, rvec gx) const {
    return call(vtable.eval_g, x, gx);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_g_prod(crvec x, crvec y,
                                                          rvec grad_gxy) const {
    return call(vtable.eval_grad_g_prod, x, y, grad_gxy);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_gi(crvec x, index_t i,
                                                      rvec grad_gi) const {
    return call(vtable.eval_grad_gi, x, i, grad_gi);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_L_prod(crvec x, crvec y,
                                                          crvec v,
                                                          rvec Hv) const {
    return call(vtable.eval_hess_L_prod, x, y, v, Hv);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_hess_L(crvec x, crvec y,
                                                     rmat H) const {
    return call(vtable.eval_hess_L, x, y, H);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_grad_f(crvec x,
                                                       rvec grad_fx) const
    -> real_t {
    return call(vtable.eval_f_grad_f, x, grad_fx);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_g(crvec x, rvec g) const
    -> real_t {
    return call(vtable.eval_f_g, x, g);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_f_grad_f_g(crvec x, rvec grad_fx,
                                                         rvec g) const
    -> real_t {
    return call(vtable.eval_f_grad_f_g, x, grad_fx, g);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_f_grad_g_prod(
    crvec x, crvec y, rvec grad_f, rvec grad_gxy) const {
    return call(vtable.eval_grad_f_grad_g_prod, x, y, grad_f, grad_gxy);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_L(crvec x, crvec y,
                                                     rvec grad_L,
                                                     rvec work_n) const {
    return call(vtable.eval_grad_L, x, y, grad_L, work_n);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_ψ_ŷ(crvec x, crvec y, crvec Σ,
                                                  rvec ŷ) const -> real_t {
    return call(vtable.eval_ψ_ŷ, x, y, Σ, ŷ);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_ψ_from_ŷ(crvec x, crvec ŷ,
                                                            rvec grad_ψ,
                                                            rvec work_n) const {
    return call(vtable.eval_grad_ψ_from_ŷ, x, ŷ, grad_ψ, work_n);
}
template <Config Conf, class Allocator>
void TypeErasedProblem<Conf, Allocator>::eval_grad_ψ(crvec x, crvec y, crvec Σ,
                                                     rvec grad_ψ, rvec work_n,
                                                     rvec work_m) const {
    return call(vtable.eval_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::eval_ψ_grad_ψ(crvec x, crvec y,
                                                       crvec Σ, rvec grad_ψ,
                                                       rvec work_n,
                                                       rvec work_m) const
    -> real_t {
    return call(vtable.eval_ψ_grad_ψ, x, y, Σ, grad_ψ, work_n, work_m);
}
template <Config Conf, class Allocator>
auto TypeErasedProblem<Conf, Allocator>::calc_ŷ_dᵀŷ(rvec g_ŷ, crvec y,
                                                    crvec Σ) const -> real_t {
    return call(vtable.calc_ŷ_dᵀŷ, g_ŷ, y, Σ);
}

struct ProblemCounters {
    unsigned proj_diff_g{};
    unsigned prox_grad_step{};
    unsigned f{};
    unsigned grad_f{};
    unsigned g{};
    unsigned grad_g_prod{};
    unsigned grad_gi{};
    unsigned hess_L_prod{};
    unsigned hess_L{};
    unsigned f_grad_f{};
    unsigned f_g{};
    unsigned f_grad_f_g{};
    unsigned grad_f_grad_g_prod{};
    unsigned grad_L{};
    unsigned ψ_ŷ{};
    unsigned grad_ψ_from_ŷ{};
    unsigned grad_ψ{};
    unsigned ψ_grad_ψ{};

    struct Time {
        std::chrono::nanoseconds proj_diff_g{};
        std::chrono::nanoseconds prox_grad_step{};
        std::chrono::nanoseconds f{};
        std::chrono::nanoseconds grad_f{};
        std::chrono::nanoseconds g{};
        std::chrono::nanoseconds grad_g_prod{};
        std::chrono::nanoseconds grad_gi{};
        std::chrono::nanoseconds hess_L_prod{};
        std::chrono::nanoseconds hess_L{};
        std::chrono::nanoseconds f_grad_f{};
        std::chrono::nanoseconds f_g{};
        std::chrono::nanoseconds f_grad_f_g{};
        std::chrono::nanoseconds grad_f_grad_g_prod{};
        std::chrono::nanoseconds grad_L{};
        std::chrono::nanoseconds ψ_ŷ{};
        std::chrono::nanoseconds grad_ψ_from_ŷ{};
        std::chrono::nanoseconds grad_ψ{};
        std::chrono::nanoseconds ψ_grad_ψ{};
    } time{};
};

namespace detail {
template <class T>
struct Timed {
    Timed(T &time) : time(time) {
        time -= std::chrono::steady_clock::now().time_since_epoch();
    }
    ~Timed() { time += std::chrono::steady_clock::now().time_since_epoch(); }
    Timed(const Timed &)            = delete;
    Timed(Timed &&)                 = delete;
    Timed &operator=(const Timed &) = delete;
    Timed &operator=(Timed &&)      = delete;
    T &time;
};
template <class T>
Timed(T &) -> Timed<T>;
} // namespace detail

template <class Problem>
struct ProblemWithCounters {
    USING_ALPAQA_CONFIG_TEMPLATE(std::remove_cvref_t<Problem>::config_t);

    // clang-format off
    void eval_proj_diff_g(crvec z, rvec p) const { ++evaluations.proj_diff_g; return timed(evaluations.time.proj_diff_g, std::bind(&std::remove_cvref_t<Problem>::eval_proj_diff_g, &problem, z, p)); }
    void eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { ++evaluations.prox_grad_step; return timed(evaluations.time.prox_grad_step, std::bind(&std::remove_cvref_t<Problem>::eval_prox_grad_step, &problem, γ, x, grad_ψ, x̂, p)); }
    real_t eval_f(crvec x) const { ++evaluations.f; return timed(evaluations.time.f, std::bind(&std::remove_cvref_t<Problem>::eval_f, &problem, x)); }
    void eval_grad_f(crvec x, rvec grad_fx) const { ++evaluations.grad_f; return timed(evaluations.time.grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f, &problem, x, grad_fx)); }
    void eval_g(crvec x, rvec gx) const { ++evaluations.g; return timed(evaluations.time.g, std::bind(&std::remove_cvref_t<Problem>::eval_g, &problem, x, gx)); }
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { ++evaluations.grad_g_prod; return timed(evaluations.time.grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_g_prod, &problem, x, y, grad_gxy)); }
    void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_gi; } { ++evaluations.grad_gi; return timed(evaluations.time.grad_gi, std::bind(&std::remove_cvref_t<Problem>::eval_grad_gi, &problem, x, i, grad_gi)); }
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L_prod; } { ++evaluations.hess_L_prod; return timed(evaluations.time.hess_L_prod, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L_prod, &problem, x, y, v, Hv)); }
    void eval_hess_L(crvec x, crvec y, rmat H) const requires requires { &std::remove_cvref_t<Problem>::eval_hess_L; } { ++evaluations.hess_L; return timed(evaluations.time.hess_L, std::bind(&std::remove_cvref_t<Problem>::eval_hess_L, &problem, x, y, H)); }
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f; } { ++evaluations.f_grad_f; return timed(evaluations.time.f_grad_f, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f, &problem, x, grad_fx)); }
    real_t eval_f_g(crvec x, rvec g) const requires requires { &std::remove_cvref_t<Problem>::eval_f_g; } { ++evaluations.f_g; return timed(evaluations.time.f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_g, &problem, x, g)); }
    real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const requires requires { &std::remove_cvref_t<Problem>::eval_f_grad_f_g; } { ++evaluations.f_grad_f_g; return timed(evaluations.time.f_grad_f_g, std::bind(&std::remove_cvref_t<Problem>::eval_f_grad_f_g, &problem, x, grad_fx, g)); }
    void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod; } { ++evaluations.grad_f_grad_g_prod; return timed(evaluations.time.grad_f_grad_g_prod, std::bind(&std::remove_cvref_t<Problem>::eval_grad_f_grad_g_prod, &problem, x, y, grad_f, grad_gxy)); }
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_L; } { ++evaluations.grad_L; return timed(evaluations.time.grad_L, std::bind(&std::remove_cvref_t<Problem>::eval_grad_L, &problem, x, y, grad_L, work_n)); }
    real_t eval_ψ_ŷ(crvec x, crvec y, crvec Σ, rvec ŷ) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ_ŷ; } { ++evaluations.ψ_ŷ; return timed(evaluations.time.ψ_ŷ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ_ŷ, &problem, x, y, Σ, ŷ)); }
    void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ; } { ++evaluations.grad_ψ_from_ŷ; return timed(evaluations.time.grad_ψ_from_ŷ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ_from_ŷ, &problem, x, ŷ, grad_ψ, work_n)); }
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_grad_ψ; } { ++evaluations.grad_ψ; return timed(evaluations.time.grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const requires requires { &std::remove_cvref_t<Problem>::eval_ψ_grad_ψ; } { ++evaluations.ψ_grad_ψ; return timed(evaluations.time.ψ_grad_ψ, std::bind(&std::remove_cvref_t<Problem>::eval_ψ_grad_ψ, &problem, x, y, Σ, grad_ψ, work_n, work_m)); }

    template <class... Args> decltype(auto) provides_eval_grad_gi(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_grad_gi; } { return problem.provides_eval_grad_gi(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_hess_L_prod(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_hess_L_prod; } { return problem.provides_eval_hess_L_prod(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_hess_L(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_hess_L; } { return problem.provides_eval_hess_L(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_f_grad_f(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_f_grad_f; } { return problem.provides_eval_f_grad_f(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_f_g(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_f_g; } { return problem.provides_eval_f_g(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_f_grad_f_g(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_f_grad_f_g; } { return problem.provides_eval_f_grad_f_g(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_grad_f_grad_g_prod(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_grad_f_grad_g_prod; } { return problem.provides_eval_grad_f_grad_g_prod(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_grad_L(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_grad_L; } { return problem.provides_eval_grad_L(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_ψ_ŷ(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_ψ_ŷ; } { return problem.provides_eval_ψ_ŷ(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_grad_ψ_from_ŷ(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_grad_ψ_from_ŷ; } { return problem.provides_eval_grad_ψ_from_ŷ(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_grad_ψ(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_grad_ψ; } { return problem.provides_eval_grad_ψ(std::forward<Args>(args)...); }
    template <class... Args> decltype(auto) provides_eval_ψ_grad_ψ(Args... args) const requires requires { &std::remove_cvref_t<Problem>::provides_eval_ψ_grad_ψ; } { return problem.provides_eval_ψ_grad_ψ(std::forward<Args>(args)...); }
    // clang-format on

    mutable ProblemCounters evaluations;
    Problem problem;

    ProblemWithCounters(const Problem &problem) : problem(problem) {}
    ProblemWithCounters(Problem &&problem) requires(!std::is_lvalue_reference_v<Problem>)
        : problem(std::forward<Problem>(problem)) {}

  private:
    template <class TimeT, class FunT>
    static decltype(auto) timed(TimeT &time, FunT &&f) {
        detail::Timed timed{time};
        return std::forward<FunT>(f)();
    }
};

} // namespace alpaqa