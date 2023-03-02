#pragma once

#include <alpaqa/ipopt-adapter-export.h>
#include <alpaqa/problem/type-erased-problem.hpp>

#include <IpTNLP.hpp>

namespace alpaqa {

/// Based on https://coin-or.github.io/Ipopt/INTERFACES.html
class IPOPT_INTERFACE_EXPORT IpoptAdapter : public Ipopt::TNLP {
  public:
    USING_ALPAQA_CONFIG(DefaultConfig);
    using Problem = TypeErasedProblem<config_t>;
    const Problem &problem;
    vec initial_guess;
    using Index  = Ipopt::Index;
    using Number = Ipopt::Number;

    struct Results {
        Ipopt::SolverReturn status = Ipopt::SolverReturn::UNASSIGNED;
        vec solution_x, solution_y, solution_g;
        real_t solution_f, infeasibility, nlp_error;
        length_t iter_count;
    } results;

    IpoptAdapter(const Problem &problem);
    IpoptAdapter(Problem &&) = delete;

    /// @name Functions required by Ipopt
    /// @{

    bool get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
                      IndexStyleEnum &index_style) override;

    bool get_bounds_info(Index n, Number *x_l, Number *x_u, Index m,
                         Number *g_l, Number *g_u) override;

    bool get_starting_point(Index n, bool init_x, Number *x,
                            [[maybe_unused]] bool init_z,
                            [[maybe_unused]] Number *z_L,
                            [[maybe_unused]] Number *z_U,
                            [[maybe_unused]] Index m,
                            [[maybe_unused]] bool init_lambda,
                            [[maybe_unused]] Number *lambda) override;

    bool eval_f(Index n, const Number *x, [[maybe_unused]] bool new_x,
                Number &obj_value) override;

    bool eval_grad_f(Index n, const Number *x, [[maybe_unused]] bool new_x,
                     Number *grad_f) override;

    bool eval_g(Index n, const Number *x, [[maybe_unused]] bool new_x, Index m,
                Number *g) override;

    bool eval_jac_g(Index n, const Number *x, [[maybe_unused]] bool new_x,
                    Index m, Index nele_jac, Index *iRow, Index *jCol,
                    Number *values) override;

    bool eval_h(Index n, const Number *x, [[maybe_unused]] bool new_x,
                Number obj_factor, Index m, const Number *lambda,
                [[maybe_unused]] bool new_lambda, Index nele_hess, Index *iRow,
                Index *jCol, Number *values) override;

    void finalize_solution(Ipopt::SolverReturn status, Index n, const Number *x,
                           [[maybe_unused]] const Number *z_L,
                           [[maybe_unused]] const Number *z_U, Index m,
                           const Number *g, const Number *lambda,
                           Number obj_value, const Ipopt::IpoptData *ip_data,
                           Ipopt::IpoptCalculatedQuantities *ip_cq) override;

    /// @}

  private:
    struct Sparsity {
        indexvec inner_idx, outer_ptr;
    } sparsity_J, sparsity_H;

    static void set_sparsity(Index n, Index m, [[maybe_unused]] Index nele,
                             Index *iRow, Index *jCol, const Sparsity &sp);
};

} // namespace alpaqa