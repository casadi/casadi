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
  long *nzj_i, *nzj_j;
  long *nzh_i, *nzh_j;

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

  CNLPModel* cnlp_model;
  OptsDict* libmad_opts;
  MadNLPExecutionStats* stats;
  MadNLPSolver* solver;
};
// C-REPLACE "casadi_madnlp_data<T1>" "struct casadi_madnlp_data"

// SYMBOL "madnlp_setup"
template<typename T1>
void casadi_madnlp_setup(casadi_madnlp_prob<T1>* p) {
}

// C-REPLACE "static_cast< casadi_madnlp_data<T1>* >" "(struct casadi_madnlp_data*)"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"
// C-REPLACE "calc_function" "casadi_oracle_call"
// C-REPLACE "casadi_error" "//casadi_error"

// SYMBOL "madnlp_constr_jac_structure"
template<typename T1>
int casadi_madnlp_constr_jac_structure(libmad_int* I, libmad_int* J, void* user_data)
{
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  std::memcpy(I, p->nzj_i, p->nnz_jac_g * sizeof(libmad_int));
  std::memcpy(J, p->nzj_j, p->nnz_jac_g * sizeof(libmad_int));
  return 0;
}

// SYMBOL "madnlp_lag_hess_structure"
template<typename T1>
int casadi_madnlp_lag_hess_structure(libmad_int* I, libmad_int* J, void* user_data)
{
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  std::memcpy(I, p->nzh_i, p->nnz_hess_l * sizeof(libmad_int));
  std::memcpy(J, p->nzh_j, p->nnz_hess_l * sizeof(libmad_int));
  return 0;
}

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

  auto ret = calc_function(&d->prob->nlp_jac_g, d_oracle);
  std::cout << "jac_g: ";
  for(int ii=0; ii<3;ii++)
  {
    std::cout << res[ii] << ", ";
  }
  std::cout << std::endl;
  std::cout << "w: ";
  for(int ii=0; ii<3;ii++)
  {
    std::cout << w[ii] << ", ";
  }
  std::cout << std::endl;
  return ret;
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

  return calc_function(&d->prob->nlp_g, d_oracle);
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
  auto ret = calc_function(&d->prob->nlp_grad_f, d_oracle);

  std::cout << "jac_f: ";
  for(int ii=0; ii<3;ii++)
  {
    std::cout << res[ii] << ", ";
  }
  std::cout << std::endl;

  return ret;
  //d->x =  (double*)w;
  //d->grad_f = res;
  //casadi_scal(p->nlp->nx, objective_scale, res);
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
  return calc_function(&d->prob->nlp_f, d_oracle);
  //d->x =  (double*)w;
  //d->f = res;
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
  
  auto ret = calc_function(&d->prob->nlp_hess_l, d_oracle);
  std::cout << "w: ";
  for(int ii=0; ii<3;ii++)
  {
    std::cout << w[ii] << ", ";
  }
  std::cout << std::endl;
  std::cout << "lam: ";
  for(int ii=0; ii<1;ii++)
  {
    std::cout << lam[ii] << ", ";
  }
  std::cout << std::endl;
  std::cout << "lam_f" << objective_scale << std::endl;
  std::cout << "hess_l: ";
  for(int ii=0; ii<2;ii++)
  {
    std::cout << res[ii] << ", ";
  }
  std::cout << std::endl;
  return ret;
}

// C-REPLACE "const_cast<T1*>" "(T1*)"

// C-REPLACE "casadi_madnlp_eval_constr_jac<T1>" "casadi_madnlp_eval_constr_jac"
// C-REPLACE "casadi_madnlp_eval_obj_grad<T1>" "casadi_madnlp_eval_obj_grad"
// C-REPLACE "casadi_madnlp_eval_obj<T1>" "casadi_madnlp_eval_obj"
// C-REPLACE "casadi_madnlp_eval_constr<T1>" "casadi_madnlp_eval_constr"
// C-REPLACE "casadi_madnlp_eval_lag_hess<T1>" "casadi_madnlp_eval_lag_hess"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "madnlp_init_mem"
template<typename T1>
int casadi_madnlp_init_mem(casadi_madnlp_data<T1>* d) {
  // Problem structure
  
  return 0;
}

// SYMBOL "madnlp_free_mem"
template<typename T1>
void madnlp_free_mem(casadi_madnlp_data<T1>* d) {
  madnlp_delete_solver(d->solver);
  d->solver = nullptr;
  madnlp_delete_stats(d->stats);
  d->stats = nullptr;
  libmad_delete_options_dict(d->libmad_opts);
}

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
  libmad_nlpmodel_create(&(d->cnlp_model),
                         "Name",
                         p_nlp->nx, p_nlp->ng,
                         p->nnz_jac_g, p->nnz_hess_l,
                         casadi_madnlp_constr_jac_structure<T1>,
                         casadi_madnlp_lag_hess_structure<T1>,
                         casadi_madnlp_eval_obj<T1>,
                         casadi_madnlp_eval_constr<T1>,
                         casadi_madnlp_eval_obj_grad<T1>,
                         casadi_madnlp_eval_constr_jac<T1>,
                         casadi_madnlp_eval_lag_hess<T1>,
                         d
    );

  madnlp_create_solver(&(d->solver), d->cnlp_model, d->libmad_opts);
  std::cout << "Created Solver" << std::endl;
  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// SYMBOL "madnlp_presolve"
template<typename T1>
void casadi_madnlp_presolve(casadi_madnlp_data<T1>* d) {
}

// SYMBOL "madnlp_c_solve"
template<typename T1>
int casadi_madnlp_solve(casadi_madnlp_data<T1>* d) {
  // Problem structure
  casadi_int k, i, column;
  const casadi_madnlp_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  d->unified_return_status = SOLVER_RET_UNKNOWN;
  d->success = 0;

  std::cout << "d->cnlp_model: " << d->cnlp_model << std::endl;
  std::cout << "d_nlp: " << d_nlp << std::endl;
  // set initial guess
  libmad_nlpmodel_set_numerics(d->cnlp_model,
                               d_nlp->z, d_nlp->lam + p_nlp->nx,
                               d_nlp->lbx, d_nlp->ubx,
                               d_nlp->lbg, d_nlp->ubg);

  int ret = madnlp_solve(d->solver, d->libmad_opts, &(d->stats));

  if (ret!=0) {
    // cleanup done in free_mem
    return ret;
  }
  // get objective
  madnlp_get_obj(d->stats, &(d_nlp->objective));
  // get primal solution for x
  madnlp_get_solution(d->stats, d_nlp->z);
  // get the bound multipliers
  madnlp_get_bound_multipliers(d->stats, d_nlp->lam);
  // get the nonlinear constraint function values
  madnlp_get_constraints(d->stats, d_nlp->z + p_nlp->nx);
  // get the nonlinear constraint multipliers
  madnlp_get_multipliers(d->stats, d_nlp->lam + p_nlp->nx);

  bool success_b;
  madnlp_get_success(d->stats, &(success_b));
  if(success_b){
    d->success = 1;
  }

  if(d->success) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
  } else {
    d->unified_return_status = SOLVER_RET_LIMITED;
  }


  return 0;
}
