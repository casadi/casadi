#include "alpaqa_problem.hpp"
#include "alpaqa_interface.hpp"
#include "casadi/core/exception.hpp"

namespace casadi {

AlpaqaProblem::AlpaqaProblem(const AlpaqaInterface& solver, AlpaqaMemory* mem)
  : alpaqa::BoxConstrProblem<alpaqa::DefaultConfig>{solver.nx_, solver.ng_},
    solver_(solver), mem_(mem) {

}

//AlpaqaProblem::AlpaqaProblem(const AlpaqaProblem &) = default;
//AlpaqaProblem& AlpaqaProblem::operator=(const AlpaqaProblem &) = default;
//AlpaqaProblem::AlpaqaProblem(AlpaqaProblem &&) noexcept = default;
//AlpaqaProblem & AlpaqaProblem::operator=(AlpaqaProblem &&) noexcept = default;
AlpaqaProblem::~AlpaqaProblem() = default;

double AlpaqaProblem::eval_f(crvec x) const {
  double obj_value;
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->res[0] = &obj_value;
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_f")==0, "Failing evaluating eval_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_f failed:" + std::string(ex.what()));
    }
  }

  return obj_value;
}

void AlpaqaProblem::eval_grad_f(crvec x, rvec grad_fx) const {
  eval_f_grad_f(x, grad_fx);
}

double AlpaqaProblem::eval_f_grad_f(crvec x, rvec grad_fx) const {
  double obj_value;
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->res[0] = &obj_value;
  mem_->res[1] = grad_fx.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_f_grad_f")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_f_grad_f failed:" + std::string(ex.what()));
    }
  }
  return obj_value;
}


template<typename T1, typename T2>
void copy(const T1* x, casadi_int n, T2* y) {
  casadi_int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

void AlpaqaProblem::eval_g(crvec x, rvec g) const {
    mem_->arg[0] = x.data();
    mem_->arg[1] = mem_->d_nlp.p;
    mem_->res[0] = g.data();
    try {
      casadi_assert(solver_.calc_function(mem_, "nlp_g")==0, "Failing evaluating eval_f_grad_f");
    } catch(KeyboardInterruptException& ex) {
      casadi_warning("KeyboardInterruptException");
      throw KeyboardInterruptException();
    } catch (std::exception& ex) {
      if (solver_.show_eval_warnings_) {
        casadi_warning("AlpaqaProblem::eval_g failed:" + std::string(ex.what()));
      }
    }
}

void AlpaqaProblem::eval_jac_g(crvec x, rindexvec inner_idx,
                                     rindexvec outer_ptr, rvec J_values) const {
  if (J_values.size()>0) {
    mem_->arg[0] = x.data();
    mem_->arg[1] = mem_->d_nlp.p;
    mem_->res[0] = J_values.data();
    try {
      casadi_assert(solver_.calc_function(mem_, "nlp_jac_g")==0, "Failing evaluating eval_f_grad_f");
    } catch(KeyboardInterruptException& ex) {
      casadi_warning("KeyboardInterruptException");
      throw KeyboardInterruptException();
    } catch (std::exception& ex) {
      if (solver_.show_eval_warnings_) {
        casadi_warning("AlpaqaProblem::eval_jac_g failed:" + std::string(ex.what()));
      }
    }
  } else {
    const Sparsity& sp = solver_.jacg_sp_;
    if (!sp.is_dense()) {
      copy(sp.row(), sp.nnz(), inner_idx.data());
      copy(sp.colind(), get_n()+1, outer_ptr.data());
    }
  }
}

void AlpaqaProblem::eval_hess_L(crvec x, crvec y, real_t scale,
                                      rindexvec inner_idx, rindexvec outer_ptr,
                                      rvec H_values) const {
  if (H_values.size()>0) {
    mem_->arg[0] = x.data();
    mem_->arg[1] = mem_->d_nlp.p;
    mem_->arg[2] = y.data();
    mem_->arg[3] = &scale;
    mem_->res[0] = H_values.data();
    try {
      casadi_assert(solver_.calc_function(mem_, "nlp_hess_L")==0, "Failing evaluating eval_f_grad_f");
    } catch(KeyboardInterruptException& ex) {
      casadi_warning("KeyboardInterruptException");
      throw KeyboardInterruptException();
    } catch (std::exception& ex) {
      if (solver_.show_eval_warnings_) {
        casadi_warning("AlpaqaProblem::eval_hess_L failed:" + std::string(ex.what()));
      }
    }
  } else {
    const Sparsity& sp = solver_.get_function("nlp_hess_L").sparsity_out(0);
    if (!sp.is_dense()) {
      copy(sp.row(), sp.nnz(), inner_idx.data());
      copy(sp.colind(), get_n()+1, outer_ptr.data());
    }
  }
}

void AlpaqaProblem::eval_hess_L_prod(crvec x, crvec y, real_t scale,
                                           crvec v, rvec Hv) const {
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = y.data();
  mem_->arg[3] = &scale;
  mem_->arg[4] = v.data();
  mem_->res[0] = Hv.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_hess_L_prod")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_hess_L_prod failed:" + std::string(ex.what()));
    }
  }
}

void AlpaqaProblem::eval_grad_g_prod(crvec, crvec, rvec) const {
  casadi_error("Not implemented");
}

void AlpaqaProblem::eval_grad_gi(crvec, index_t, rvec) const {
  casadi_error("Not implemented");
}

double AlpaqaProblem::eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const {
  double res;
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = y.data();
  mem_->arg[3] = Σ.data();
  mem_->arg[4] = this->D.lowerbound.data();
  mem_->arg[5] = this->D.upperbound.data();
  mem_->res[0] = &res;
  mem_->res[1] = ŷ.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_psi")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_psi failed:" + std::string(ex.what()));
    }
  }
  return res;
}

double AlpaqaProblem::eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec,
                                   rvec) const {
  double res;
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = y.data();
  mem_->arg[3] = Σ.data();
  mem_->arg[4] = this->D.lowerbound.data();
  mem_->arg[5] = this->D.upperbound.data();
  mem_->res[0] = &res;
  mem_->res[1] = grad_ψ.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_grad_psi")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_grad_psi failed:" + std::string(ex.what()));
    }
  }
  return res;
}

void AlpaqaProblem::eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ,
                                      rvec a, rvec b) const {
  eval_ψ_grad_ψ(x, y, Σ, grad_ψ, a, b);
}

void AlpaqaProblem::eval_grad_L(crvec x, crvec y, rvec grad_L,
                                      rvec) const {
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = y.data();
  mem_->res[0] = grad_L.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_grad_L")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_grad_L failed:" + std::string(ex.what()));
    }
  }
}

void AlpaqaProblem::eval_hess_ψ(crvec x, crvec y, crvec Σ, real_t scale,
                                      rindexvec inner_idx, rindexvec outer_ptr,
                                      rvec H_values) const {
    if (H_values.size()>0) {
      mem_->arg[0] = x.data();
      mem_->arg[1] = mem_->d_nlp.p;
      mem_->arg[2] = y.data();
      mem_->arg[3] = Σ.data();
      mem_->arg[4] = &scale;
      mem_->arg[5] = this->D.lowerbound.data();
      mem_->arg[6] = this->D.upperbound.data();
      mem_->res[0] = H_values.data();
      try {
        casadi_assert(solver_.calc_function(mem_, "nlp_hess_psi")==0, "Failing evaluating eval_f_grad_f");
      } catch(KeyboardInterruptException& ex) {
        casadi_warning("KeyboardInterruptException");
        throw KeyboardInterruptException();
      } catch (std::exception& ex) {
        if (solver_.show_eval_warnings_) {
          casadi_warning("AlpaqaProblem::eval_hess_psi failed:" + std::string(ex.what()));
        }
      }
    } else {
      const Sparsity& sp = solver_.get_function("nlp_hess_psi").sparsity_out(0);
      if (!sp.is_dense()) {
        copy(sp.row(), sp.nnz(), inner_idx.data());
        copy(sp.colind(), get_n()+1, outer_ptr.data());
      }
    }
}

void AlpaqaProblem::eval_hess_ψ_prod(crvec x, crvec y, crvec Σ,
                                           real_t scale, crvec v,
                                           rvec Hv) const {
  mem_->arg[0] = x.data();
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = y.data();
  mem_->arg[3] = Σ.data();
  mem_->arg[4] = &scale;
  mem_->arg[5] = this->D.lowerbound.data();
  mem_->arg[6] = this->D.upperbound.data();
  mem_->arg[7] = v.data();
  mem_->res[0] = Hv.data();
  try {
    casadi_assert(solver_.calc_function(mem_, "nlp_hess_psi_prod")==0, "Failing evaluating eval_f_grad_f");
  } catch(KeyboardInterruptException& ex) {
    casadi_warning("KeyboardInterruptException");
    throw KeyboardInterruptException();
  } catch (std::exception& ex) {
    if (solver_.show_eval_warnings_) {
      casadi_warning("AlpaqaProblem::eval_hess_psi_prod failed:" + std::string(ex.what()));
    }
  }
}

alpaqa::DefaultConfig::length_t AlpaqaProblem::get_hess_L_num_nonzeros() const {
    const Sparsity& sp = solver_.get_function("nlp_hess_L").sparsity_out(0);
    return sp.is_dense() ? 0 : sp.nnz();
}

alpaqa::DefaultConfig::length_t AlpaqaProblem::get_hess_ψ_num_nonzeros() const {
    const Sparsity& sp = solver_.get_function("nlp_hess_psi").sparsity_out(0);
    return sp.is_dense() ? 0 : sp.nnz();
}


alpaqa::DefaultConfig::length_t AlpaqaProblem::get_jac_g_num_nonzeros() const {
  const Sparsity& sp = solver_.jacg_sp_;
  return sp.is_dense() ? 0 : sp.nnz();
}

} // namespace casadi