//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_UNKNOWN" "1"
// C-REPLACE "SOLVER_RET_LIMITED" "2"

// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "const_cast<int*>" "(int*) "

// SYMBOL "madnlp_prob"
template<typename T1>
struct casadi_madnlp_prob {
  const casadi_nlpsol_prob<T1>* nlp;

  // 1-index cco format
  madnlp_int *nzj_i, *nzj_j;
  madnlp_int *nzh_i, *nzh_j;

  casadi_int nnz_hess_l, nnz_jac_g;

  OracleCallback nlp_hess_l;
  OracleCallback nlp_jac_g;
  OracleCallback nlp_grad_f;
  OracleCallback nlp_f;
  OracleCallback nlp_g;
};
// C-REPLACE "casadi_madnlp_prob<T1>" "struct casadi_madnlp_prob"

// SYMBOL "madnlp_data"
template<typename T1>
struct casadi_madnlp_data {
  // Problem structure
  const casadi_madnlp_prob<T1>* prob;
  // Problem structure
  casadi_nlpsol_data<T1>* nlp;

  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1 *w;

  // temporary data
  T1 *x, *lam_g, *g, *f, *jac_g, *grad_f, *hess_l, *grad_l;

  int unified_return_status;
  int success;

  //struct blasfeo_dvec v, r;
  //struct blasfeo_dmat R;

  struct MadnlpCInterface c_interface;
  struct MadnlpCStats stats;
  struct MadnlpCSolver *solver;
};
// C-REPLACE "casadi_madnlp_data<T1>" "struct casadi_madnlp_data"

// SYMBOL "madnlp_setup"
template<typename T1>
void casadi_madnlp_setup(casadi_madnlp_prob<T1>* p) {

}

// SYMBOL "madnlp_init_mem"
template<typename T1>
int madnlp_init_mem(casadi_madnlp_data<T1>* d) {
  return 0;
}

// SYMBOL "madnlp_free_mem"
template<typename T1>
void madnlp_free_mem(casadi_madnlp_data<T1>* d) {
  //Highs_destroy(d->madnlp);
}
// C-REPLACE "static_cast< casadi_madnlp_data<T1>* >" "(struct casadi_madnlp_data*)"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"
// C-REPLACE "calc_function" "casadi_oracle_call"
// C-REPLACE "casadi_error" "//casadi_error"

// SYMBOL "madnlp_eval_constr_jac"
template<typename T1>
int casadi_madnlp_eval_constr_jac(const T1* w, T1* res, void* user_data) {
  casadi_int i;
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_jac_g, d_oracle);
  //store res madnlp_data to be accessible from casadi
  d->jac_g = res;
  d->x =  (double*)w;
  return 0;
}

// SYMBOL "madnlp_eval_constr"
template<typename T1>
int casadi_madnlp_eval_constr(const T1* w, T1* res, void* user_data) {
  casadi_int i,k,column;
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_g, d_oracle);

  d->x =  (double*)w;
  d->g = res;

  return 0;
}

// SYMBOL "madnlp_eval_obj_grad"
template<typename T1>
int casadi_madnlp_eval_obj_grad(const T1* w, T1* res, void* user_data) {
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  // TODO add argument
  T1 objective_scale = 1.0;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_grad_f, d_oracle);

  d->x =  (double*)w;
  d->grad_f = res;

  //casadi_scal(p->nlp->nx, objective_scale, res);
  return 0;
}

// SYMBOL "madnlp_eval_obj"
template<typename T1>
int casadi_madnlp_eval_obj(const T1* w, T1* res, void* user_data) {
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  // TODO add argument
  T1 objective_scale = 1.0;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_f, d_oracle);

  d->x =  (double*)w;
  d->f = res;

  return 0;
}

// SYMBOL "madnlp_eval_obj"
template<typename T1>
int casadi_madnlp_eval_lag_hess(T1 objective_scale, const T1* w, const T1* lam,
                                            T1* res, void* user_data){
  casadi_int k, column, i;
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  // evaluate hessian with given parameters / inputs
  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->arg[2] = &objective_scale;
  d_oracle->arg[3] = lam;
  d_oracle->res[0] = res;

  calc_function(&d->prob->nlp_hess_l, d_oracle);
  return 0;
}

// C-REPLACE "const_cast<T1*>" "(T1*)"

// C-REPLACE "casadi_madnlp_eval_constr_jac<T1>" "casadi_madnlp_eval_constr_jac"
// C-REPLACE "casadi_madnlp_eval_obj_grad<T1>" "casadi_madnlp_eval_obj_grad"
// C-REPLACE "casadi_madnlp_eval_obj<T1>" "casadi_madnlp_eval_obj"
// C-REPLACE "casadi_madnlp_eval_constr<T1>" "casadi_madnlp_eval_constr"
// C-REPLACE "casadi_madnlp_eval_lag_hess<T1>" "casadi_madnlp_eval_lag_hess"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "madnlp_work"
template<typename T1>
void casadi_madnlp_work(const casadi_madnlp_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->nlp->nx+p->nlp->ng)); // pv

}

// SYMBOL "madnlp_init"
template<typename T1>
void casadi_madnlp_init(casadi_madnlp_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Problem structure
  const casadi_madnlp_prob<T1>* p = d->prob;
  //casadi_oracle_data<T1>* d_oracle = d->nlp->oracle;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;

  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// SYMBOL "madnlp_presolve"
template<typename T1>
void casadi_madnlp_presolve(casadi_madnlp_data<T1>* d) {
  casadi_int k, i, start, stop, nx;
  // Problem structure
  const casadi_madnlp_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  struct MadnlpCInterface* c_interface = &d->c_interface;

  c_interface->eval_constr_jac = casadi_madnlp_eval_constr_jac<T1>;
  c_interface->eval_obj_grad = casadi_madnlp_eval_obj_grad<T1>;
  c_interface->eval_obj = casadi_madnlp_eval_obj<T1>;
  c_interface->eval_constr = casadi_madnlp_eval_constr<T1>;
  c_interface->eval_lag_hess = casadi_madnlp_eval_lag_hess<T1>;
  c_interface->nw = p_nlp->nx;
  c_interface->nc = p_nlp->ng;
  c_interface->nnzo = p_nlp->nx;
  c_interface->nnzj = d->prob->nnz_jac_g;
  c_interface->nnzh = d->prob->nnz_hess_l;
  c_interface->user_data = d;

  c_interface->nzj_i = d->prob->nzj_i;
  c_interface->nzj_j = d->prob->nzj_j;
  c_interface->nzh_i = d->prob->nzh_i;
  c_interface->nzh_j = d->prob->nzh_j;

  d->solver = madnlp_c_create(&d->c_interface);

}

// SYMBOL "madnlp_c_solve"
template<typename T1>
void casadi_madnlp_solve(casadi_madnlp_data<T1>* d) {
  // Problem structure
  casadi_int k, i, column;
  const casadi_madnlp_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  d->unified_return_status = SOLVER_RET_UNKNOWN;
  d->success = 0;

  const struct MadnlpCNumericIn* in = madnlp_c_input(d->solver);
  // set initial guess
  casadi_copy(d_nlp->z, p_nlp->nx, in->x0);
  casadi_copy(d_nlp->lam + p_nlp->nx, p_nlp->ng, in->l0);
  casadi_copy(d_nlp->lbx, p_nlp->nx, in->lbx);
  casadi_copy(d_nlp->ubx, p_nlp->nx, in->ubx);
  casadi_copy(d_nlp->lbg, p_nlp->ng, in->lbg);
  casadi_copy(d_nlp->ubg, p_nlp->ng, in->ubg);

  madnlp_int ret = madnlp_c_solve(d->solver);

  const struct MadnlpCNumericOut* out = madnlp_c_output(d->solver);

  d_nlp->objective = out->obj[0];

  for (casadi_int i=0; i<p_nlp->nx; ++i) {
    // get primal solution
    d_nlp->z[i] = out->sol[i];
    // get dual solution x
    d_nlp->lam[i] = out->mul_U[i]-out->mul_L[i];
  }
  // Unpack dual solution
  for (casadi_int i=0; i<p_nlp->ng; ++i) {
    d_nlp->z[p_nlp->nx+i] = out->con[i];
    d_nlp->lam[p_nlp->nx+i] = out->mul[i];
  }

  const struct MadnlpCStats* stats = madnlp_c_get_stats(d->solver);

  printf("iter %d\n", stats->iter);

  d->stats.iter = stats->iter;
  d->stats.status = stats->status;
  d->stats.dual_feas = stats->dual_feas;
  d->stats.primal_feas = stats->primal_feas;

  if (d->stats.status==MADNLP_SOLVE_SUCCEEDED || 
      d->stats.status==MADNLP_SOLVED_TO_ACCEPTABLE_LEVEL) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
  } else if (d->stats.status==MADNLP_MAXIMUM_ITERATIONS_EXCEEDED ||
      d->stats.status==MADNLP_MAXIMUM_WALLTIME_EXCEEDED) {
    d->unified_return_status = SOLVER_RET_LIMITED;
  }

  d->success = (d->unified_return_status == SOLVER_RET_SUCCESS);

  madnlp_c_destroy(d->solver);
}