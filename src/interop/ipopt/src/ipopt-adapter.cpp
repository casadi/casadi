#include <alpaqa/ipopt/ipopt-adapter.hpp>

#include <IpIpoptCalculatedQuantities.hpp>

namespace alpaqa {

IpoptAdapter::IpoptAdapter(const Problem &problem) : problem(problem) {
    const length_t n = problem.get_n();
    if (length_t nnz_H = problem.get_hess_L_num_nonzeros()) {
        sparsity_H.inner_idx.resize(nnz_H);
        sparsity_H.outer_ptr.resize(n + 1);
        mvec null{nullptr, 0};
        problem.eval_hess_L(null, null, 0, sparsity_H.inner_idx,
                            sparsity_H.outer_ptr, null);
    }
    if (length_t nnz_J = problem.get_jac_g_num_nonzeros()) {
        sparsity_J.inner_idx.resize(nnz_J);
        sparsity_J.outer_ptr.resize(n + 1);
        mvec null{nullptr, 0};
        problem.eval_jac_g(null, sparsity_J.inner_idx, sparsity_J.outer_ptr,
                           null);
    }
}

bool IpoptAdapter::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
                                Index &nnz_h_lag, IndexStyleEnum &index_style) {
    n         = static_cast<Index>(problem.get_n());
    m         = static_cast<Index>(problem.get_m());
    nnz_jac_g = static_cast<Index>(sparsity_J.inner_idx.size());
    if (nnz_jac_g == 0)
        nnz_jac_g = n * m;
    nnz_h_lag = static_cast<Index>(sparsity_H.inner_idx.size());
    if (nnz_h_lag == 0)
        nnz_h_lag = n * n;
    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;
    return true;
}

bool IpoptAdapter::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m,
                                   Number *g_l, Number *g_u) {
    const auto &C = problem.get_box_C();
    mvec{x_l, n}  = C.lowerbound;
    mvec{x_u, n}  = C.upperbound;
    const auto &D = problem.get_box_D();
    mvec{g_l, m}  = D.lowerbound;
    mvec{g_u, m}  = D.upperbound;
    return true;
}

bool IpoptAdapter::get_starting_point(Index n, bool init_x, Number *x,
                                      bool init_z, Number *z_L, Number *z_U,
                                      Index m, bool init_lambda,
                                      Number *lambda) {
    if (init_x) {
        if (initial_guess.size() > 0)
            mvec{x, n} = initial_guess;
        else
            mvec{x, n}.setZero();
    }
    if (init_z) {
        if (initial_guess_bounds_multipliers_l.size() > 0)
            mvec{z_L, n} = initial_guess_bounds_multipliers_l;
        else
            mvec{z_L, n}.setZero();
        if (initial_guess_bounds_multipliers_u.size() > 0)
            mvec{z_U, n} = initial_guess_bounds_multipliers_u;
        else
            mvec{z_U, n}.setZero();
    }
    if (init_lambda) {
        if (initial_guess_multipliers.size() > 0)
            mvec{lambda, m} = initial_guess_multipliers;
        else
            mvec{lambda, m}.setZero();
    }
    return true;
}

bool IpoptAdapter::eval_f(Index n, const Number *x, [[maybe_unused]] bool new_x,
                          Number &obj_value) {
    obj_value = problem.eval_f(cmvec{x, n});
    return true;
}

bool IpoptAdapter::eval_grad_f(Index n, const Number *x,
                               [[maybe_unused]] bool new_x, Number *grad_f) {
    problem.eval_grad_f(cmvec{x, n}, mvec{grad_f, n});
    return true;
}

bool IpoptAdapter::eval_g(Index n, const Number *x, [[maybe_unused]] bool new_x,
                          Index m, Number *g) {
    problem.eval_g(cmvec{x, n}, mvec{g, m});
    return true;
}

bool IpoptAdapter::eval_jac_g(Index n, const Number *x,
                              [[maybe_unused]] bool new_x, Index m,
                              Index nele_jac, Index *iRow, Index *jCol,
                              Number *values) {
    if (!problem.provides_eval_jac_g())
        throw std::logic_error("Missing required function: eval_jac_g");
    if (values == nullptr) // Initialize sparsity
        set_sparsity(n, m, nele_jac, iRow, jCol, sparsity_J);
    else // Evaluate values
        problem.eval_jac_g(cmvec{x, n}, sparsity_J.inner_idx,
                           sparsity_J.outer_ptr, mvec{values, nele_jac});
    return true;
}

bool IpoptAdapter::eval_h(Index n, const Number *x, [[maybe_unused]] bool new_x,
                          Number obj_factor, Index m, const Number *lambda,
                          [[maybe_unused]] bool new_lambda, Index nele_hess,
                          Index *iRow, Index *jCol, Number *values) {
    if (!problem.provides_eval_hess_L())
        throw std::logic_error("Missing required function: eval_hess_L");
    // Initialize sparsity
    if (values == nullptr) {
        set_sparsity(n, n, nele_hess, iRow, jCol, sparsity_H);
    }
    // Evaluate values
    else {
        problem.eval_hess_L(cmvec{x, n}, cmvec{lambda, m}, obj_factor,
                            sparsity_H.inner_idx, sparsity_H.outer_ptr,
                            mvec{values, nele_hess});
        // For dense matrices, set lower triangle to zero
        // TODO: make this more efficient in alpaqa problem interface
        if (sparsity_H.inner_idx.size() == 0) {
            mmat H{values, n, n};
            for (Index c = 0; c < n; ++c)
                for (Index r = c + 1; r < n; ++r)
                    H(r, c) = 0;
        }
    }
    return true;
}
void IpoptAdapter::finalize_solution(Ipopt::SolverReturn status, Index n,
                                     const Number *x, const Number *z_L,
                                     const Number *z_U, Index m,
                                     const Number *g, const Number *lambda,
                                     Number obj_value,
                                     const Ipopt::IpoptData *ip_data,
                                     Ipopt::IpoptCalculatedQuantities *ip_cq) {
    results.status        = status;
    results.solution_x    = cmvec{x, n};
    results.solution_z_L  = cmvec{z_L, n};
    results.solution_z_U  = cmvec{z_U, n};
    results.solution_y    = cmvec{lambda, m};
    results.solution_g    = cmvec{g, m};
    results.solution_f    = obj_value;
    results.infeasibility = ip_cq->curr_constraint_violation();
    results.nlp_error     = ip_cq->unscaled_curr_nlp_error();
    results.iter_count    = ip_data->iter_count();
}

void IpoptAdapter::set_sparsity(Index n, Index m, [[maybe_unused]] Index nele,
                                Index *iRow, Index *jCol, const Sparsity &sp) {
    // sparse
    if (sp.inner_idx.size() > 0) {
        Index l = 0; // column major, jacobian is m×n, hessian is n×n
        for (Index c = 0; c < n; ++c) {
            auto inner_start = static_cast<Index>(sp.outer_ptr(c));
            auto inner_end   = static_cast<Index>(sp.outer_ptr(c + 1));
            for (Index i = inner_start; i < inner_end; ++i) {
                assert(l < nele);
                jCol[l] = c;
                iRow[l] = static_cast<Index>(sp.inner_idx(i));
                ++l;
            }
        }
        assert(l == nele);
    }
    // dense
    else {
        Index l = 0; // column major, jacobian is m×n, hessian is n×n
        for (Index c = 0; c < n; ++c) {
            for (Index r = 0; r < m; ++r) {
                assert(l < nele);
                iRow[l] = r;
                jCol[l] = c;
                ++l;
            }
        }
        assert(l == nele);
    }
}

} // namespace alpaqa
