#ifndef ALPAQA_PROBLEM_H
#define ALPAQA_PROBLEM_H

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>

namespace casadi {

// Forward declarations
class AlpaqaInterface;
struct AlpaqaMemory;

class AlpaqaProblem : public alpaqa::BoxConstrProblem<alpaqa::DefaultConfig> {
  public:
    AlpaqaProblem(const AlpaqaInterface& solver, AlpaqaMemory* mem);
    ~AlpaqaProblem();

    //AlpaqaProblem(const AlpaqaProblem &);
    //AlpaqaProblem &operator=(const AlpaqaProblem &);
    //AlpaqaProblem(AlpaqaProblem &&) noexcept;
    //AlpaqaProblem &operator=(AlpaqaProblem &&) noexcept;

    // clang-format off
    [[nodiscard]] real_t eval_f(crvec x) const;
    void eval_grad_f(crvec x, rvec grad_fx) const;
    real_t eval_f_grad_f(crvec x, rvec grad_fx) const; // NOLINT(*nodiscard)
    void eval_g(crvec x, rvec g) const;
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;
    void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const;
    real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const; // NOLINT(*nodiscard)
    void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const;
    [[nodiscard]] real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const;
    void eval_grad_gi(crvec x, index_t i, rvec grad_i) const;
    [[nodiscard]] length_t get_jac_g_num_nonzeros() const;
    void eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const;
    void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const;
    [[nodiscard]] length_t get_hess_L_num_nonzeros() const;
    void eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const;
    void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const;
    [[nodiscard]] length_t get_hess_ψ_num_nonzeros() const;
    void eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const;
    // clang-format on

    /// @see @ref TypeErasedProblem::provides_eval_grad_gi
    [[nodiscard]] bool provides_eval_grad_gi() const { return false; }
    /// @see @ref TypeErasedProblem::provides_eval_jac_g
    [[nodiscard]] bool provides_eval_jac_g() const { return true; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L_prod
    [[nodiscard]] bool provides_eval_hess_L_prod() const { return true; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_L
    [[nodiscard]] bool provides_eval_hess_L() const { return true; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ_prod
    [[nodiscard]] bool provides_eval_hess_ψ_prod() const { return true; }
    /// @see @ref TypeErasedProblem::provides_eval_hess_ψ
    [[nodiscard]] bool provides_eval_hess_ψ() const { return true; }

  private:
    const AlpaqaInterface& solver_;
    AlpaqaMemory* mem_;
};



} // namespace casadi

#endif /* ALPAQA_PROBLEM_H */
