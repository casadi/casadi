#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/util/not-implemented.hpp>

#include <functional>
#include <memory>
#include <type_traits>

namespace alpaqa {

/// Base class for problem definitions.
template <Config Conf = DefaultConfig>
class ProblemBase {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    /// Number of decision variables, dimension of x
    length_t n;
    /// Number of constraints, dimension of g(x) and z
    length_t m;

    ProblemBase(length_t n, ///< Number of decision variables
                length_t m) ///< Number of constraints
        : n{n}, m{m} {}

    ProblemBase(const ProblemBase &)            = default;
    ProblemBase &operator=(const ProblemBase &) = default;
    ProblemBase(ProblemBase &&)                 = default;
    ProblemBase &operator=(ProblemBase &&)      = default;

    virtual std::unique_ptr<ProblemBase> clone() const & = 0;
    virtual std::unique_ptr<ProblemBase> clone()      && = 0;
    virtual ~ProblemBase()                               = default;

    /// @name Basic functions
    /// @{

    /// Function that evaluates the cost, @f$ f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    virtual real_t eval_f(crvec x) const;
    /// Function that evaluates the gradient of the cost, @f$ \nabla f(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] grad_fx
    ///         Gradient of cost function @f$ \nabla f(x) \in \R^n @f$
    virtual void eval_grad_f(crvec x, rvec grad_fx) const;
    /// Function that evaluates the constraints, @f$ g(x) @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [out] gx
    ///         Value of the constraints @f$ g(x) \in \R^m @f$
    virtual void eval_g(crvec x, rvec gx) const;
    /// Function that evaluates the gradient of the constraints times a vector,
    /// @f$ \nabla g(x)\,y @f$
    /// @param  [in] x
    ///         Decision variable @f$ x \in \R^n @f$
    /// @param  [in] y
    ///         Vector @f$ y \in \R^m @f$ to multiply the gradient by
    /// @param  [out] grad_gxy
    ///         Gradient of the constraints
    ///         @f$ \nabla g(x)\,y \in \R^n @f$
    virtual void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;

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
    virtual void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const;
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
    virtual void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const;
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
    virtual void eval_hess_L(crvec x, crvec y, rmat H) const;

    /// @}

    /// @name Combined evaluations
    /// @{

    /// Evaluate both @f$ f(x) @f$ and its gradient, @f$ \nabla f(x) @f$.
    virtual real_t eval_f_grad_f(crvec x, rvec grad_fx) const;
    /// Evaluate both @f$ f(x) @f$ and @f$ g(x) @f$.
    virtual real_t eval_f_g(crvec x, rvec g) const;
    /// Evaluate @f$ f(x) @f$, its gradient @f$ \nabla f(x) @f$ and @f$ g(x) @f$.
    virtual real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const;
    /// Evaluate both @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$.
    virtual void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f,
                                         rvec grad_gxy) const;
    /// Evaluate the gradient of the Lagrangian
    /// @f$ \nabla_x L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
    virtual void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;

    /// @}

    /// @name Augmented Lagrangian
    /// @{

    /// Calculate both ψ(x) and the vector ŷ that can later be used to compute
    /// ∇ψ.
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    ///   \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \hat y = \Sigma\, \left(g(x) + \Sigma^{-1}y - \Pi_D\left(g(x)
    ///   + \Sigma^{-1}y\right)\right) @f]
    virtual real_t eval_ψ_ŷ(crvec x, ///< [in]  Decision variable @f$ x @f$
                            crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                            crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                            rvec ŷ   ///< [out] @f$ \hat y @f$
    ) const;
    /// Calculate ∇ψ(x) using ŷ.
    virtual void
    eval_grad_ψ_from_ŷ(crvec x,     ///< [in]  Decision variable @f$ x @f$
                       crvec ŷ,     ///< [in]  @f$ \hat y @f$
                       rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                       rvec work_n  ///<       Dimension @f$ n @f$
    ) const;
    /// Calculate the gradient ∇ψ(x).
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    virtual void eval_grad_ψ(crvec x, ///< [in]  Decision variable @f$ x @f$
                             crvec y, ///< [in]  Lagrange multipliers @f$ y @f$
                             crvec Σ, ///< [in]  Penalty weights @f$ \Sigma @f$
                             rvec grad_ψ, ///< [out] @f$ \nabla \psi(x) @f$
                             rvec work_n, ///<       Dimension @f$ n @f$
                             rvec work_m  ///<       Dimension @f$ m @f$
    ) const;
    /// Calculate both ψ(x) and its gradient ∇ψ(x).
    /// @f[ \psi(x) = f(x) + \tfrac{1}{2}
    /// \text{dist}_\Sigma^2\left(g(x) + \Sigma^{-1}y,\;D\right) @f]
    /// @f[ \nabla \psi(x) = \nabla f(x) + \nabla g(x)\,\hat y(x) @f]
    virtual real_t
    eval_ψ_grad_ψ(crvec x,     ///< [in]  Decision variable @f$ x @f$
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

    /// @name Constraint sets
    /// @{

    virtual const Box &get_C() const = 0;
    virtual const Box &get_D() const = 0;

    /// @}
};

/**
 * @brief   Problem description for minimization problems.
 * 
 * @f[ \begin{aligned}
 *  & \underset{x}{\text{minimize}}
 *  & & f(x) &&&& f : \R^n \rightarrow \R \\
 *  & \text{subject to}
 *  & & \underline{x} \le x \le \overline{x} \\
 *  &&& \underline{z} \le g(x) \le \overline{z} &&&& g :
 *  \R^n \rightarrow \R^m
 * \end{aligned} @f]
 * <br>
 * @f[ \begin{aligned}
 *  C &\defeq \defset{x\in\R^n}{\underline x \le x \le \overline x} \\
 *  D &\defeq \defset{z\in\R^m}{\underline z \le z \le \overline z} \\
 * \end{aligned} @f]
 *
 * The functions in the “Basic functions” section have to be implemented by the
 * user. Functions in the “Combined evaluations” and “Augmented Lagrangian” 
 * sections are optional, by default, they are computed by evaluating the 
 * “Basic functions”, but performance can be improved by providing functions 
 * that directly compute multiple quantities at once.
 *
 * @ingroup grp_Problems
 */
template <Config Conf = DefaultConfig>
class Problem : public ProblemBase<Conf> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using ProblemBase = alpaqa::ProblemBase<Conf>;
    using Box         = typename ProblemBase::Box;
    using ProblemBase::m;
    using ProblemBase::n;

    /// The parameter vector
    vec param;
    /// Constraints of the decision variables, @f$ x \in C @f$
    Box C{vec::Constant(this->n, +inf<Conf>),
          vec::Constant(this->n, -inf<Conf>)};
    /// Other constraints, @f$ g(x) \in D @f$
    Box D{vec::Constant(this->m, +inf<Conf>),
          vec::Constant(this->m, -inf<Conf>)};

    Problem(length_t n, ///< Number of decision variables
            length_t m, ///< Number of constraints
            vec param)  ///< Problem parameter
        : ProblemBase{n, m}, param{std::move(param)} {}

    Problem(length_t n,               ///< Number of decision variables
            length_t m,               ///< Number of constraints
            length_t p = length_t{0}) ///< Size of the problem parameter
        : Problem{n, m, vec::Constant(p, NaN<config_t>)} {}

    Problem(length_t n, ///< Number of decision variables
            length_t m, ///< Number of constraints
            vec param,  ///< Problem parameter
            Box C,      ///< Box constraints on x
            Box D)      ///< Box constraints on g(x)
        : ProblemBase{n, m}, param{std::move(param)}, C{std::move(C)},
          D{std::move(D)} {}

    Problem(length_t n, ///< Number of decision variables
            length_t m, ///< Number of constraints
            length_t p, ///< Size of the problem parameter
            Box C,      ///< Box constraints on x
            Box D)      ///< Box constraints on g(x)
        : Problem{n, m, vec::Constant(p, NaN<config_t>), std::move(C),
                  std::move(D)} {}

    Problem(length_t n, ///< Number of decision variables
            length_t m, ///< Number of constraints
            Box C,      ///< Box constraints on x
            Box D)      ///< Box constraints on g(x)
        : Problem{n, m, length_t{0}, std::move(C), std::move(D)} {}

    std::unique_ptr<ProblemBase> clone() const & override;
    std::unique_ptr<ProblemBase> clone() && override;

    void set_C(Box C) {
        assert(C.lowerbound.size() == n);
        assert(C.upperbound.size() == n);
        this->C = std::move(C);
    }
    Box &get_C() { return C; }
    const Box &get_C() const override { return C; }

    void set_D(Box D) {
        assert(D.lowerbound.size() == m);
        assert(D.upperbound.size() == m);
        this->D = std::move(D);
    }
    Box &get_D() { return D; }
    const Box &get_D() const override { return D; }

    /// @name Parameters
    /// @{

    void set_param(vec param) {
        assert(param.size() == this->param.size());
        this->param = std::move(param);
    }
    vec &get_param() { return param; }
    const vec &get_param() const { return param; }

    /// @}
};

/// @ref Problem class that allows specifying the basic functions as C++
/// `std::function`s.
/// @ingroup grp_Problems
template <Config Conf = DefaultConfig>
class FunctionalProblem : public Problem<Conf> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using Problem<Conf>::Problem;

    std::function<real_t(crvec)> f;
    std::function<void(crvec, rvec)> grad_f;
    std::function<void(crvec, rvec)> g;
    std::function<void(crvec, crvec, rvec)> grad_g_prod;
    std::function<void(crvec, index_t, rvec)> grad_gi;
    std::function<void(crvec, crvec, crvec, rvec)> hess_L_prod;
    std::function<void(crvec, crvec, rmat)> hess_L;

    real_t eval_f(crvec x) const override;
    void eval_grad_f(crvec x, rvec grad_fx) const override;
    void eval_g(crvec x, rvec gx) const override;
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const override;
    void eval_grad_gi(crvec x, index_t i, rvec gr_gi) const override;
    void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const override;
    void eval_hess_L(crvec x, crvec y, rmat H) const override;

    FunctionalProblem(const FunctionalProblem &)            = default;
    FunctionalProblem &operator=(const FunctionalProblem &) = default;
    FunctionalProblem(FunctionalProblem &&)                 = default;
    FunctionalProblem &operator=(FunctionalProblem &&)      = default;

    std::unique_ptr<ProblemBase<Conf>> clone() const & override;
    std::unique_ptr<ProblemBase<Conf>> clone() && override;
};

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ProblemBase, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ProblemBase, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ProblemBase, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ProblemBase, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, ProblemBase, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, Problem, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, Problem, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, Problem, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, Problem, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, Problem, EigenConfigq);
#endif

ALPAQA_EXPORT_EXTERN_TEMPLATE(class, FunctionalProblem, DefaultConfig);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, FunctionalProblem, EigenConfigf);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, FunctionalProblem, EigenConfigd);
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, FunctionalProblem, EigenConfigl);
#ifdef ALPAQA_WITH_QUAD_PRECISION
ALPAQA_EXPORT_EXTERN_TEMPLATE(class, FunctionalProblem, EigenConfigq);
#endif

} // namespace alpaqa
