#pragma once

#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/util/alloc-check.hpp>

namespace alpaqa {

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
    std::function<void(crvec, index_t, rvec)> grad_gi;
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
    void eval_grad_gi(crvec x, index_t i, rvec grad_gix) const { ScopedMallocAllower ma; grad_gi(x, i, grad_gix); }
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

    /// @see @ref TypeErasedProblem::provides_eval_grad_gi
    bool provides_eval_grad_gi() const { return bool{grad_gi}; }
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

} // namespace alpaqa
