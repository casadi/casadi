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

//#ifdef __cplusplus
//#include "casadi/core/nlpsol_impl.hpp"
//#include "MadnlpCInterface.h"
//namespace casadi {
//#endif


// SYMBOL "madnlp_prob"
template<typename T1>
struct casadi_madnlp_prob {
  const casadi_nlpsol_prob<T1>* nlp;

  //  CCO Coordinate sparse storage
  casadi_int *nz_hess_l_i, *nz_hess_l_j;
  casadi_int *nz_jac_g_i, *nz_jaq_g_j;
  casadi_int *nz_grad_f_i, *nz_grad_f_j;

  // CSS column sparse storage
  const casadi_int *jac_g_ccs;
  const casadi_int *hess_l_ccs;
  const casadi_int *grad_f_ccs;
  // const casadi_int *sp_h, *sp_a, *sp_g;

  casadi_int nnz_hess_l, nnz_jac_g, nnz_grad_f;

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
  p->nnz_hess_l = casadi_sp_nnz(p->hess_l_ccs);
  p->nnz_jac_g = casadi_sp_nnz(p->jac_g_ccs);
  p->nnz_grad_f = casadi_sp_nnz(p->grad_f_ccs);
}

// SYMBOL "madnlp_init_mem"
template<typename T1>
int madnlp_init_mem(casadi_madnlp_data<T1>* d) {
  return 0;
}

// SYMBOL "madnlp_sparsity"
template<typename T1>
void casadi_madnlp_sparsity(const casadi_int* sp, size_t *coord_i, size_t *coord_j, T1 test) {
    // convert ccs to cco
    casadi_int ncol = sp[1];
    const casadi_int* colind = sp+2;
    const casadi_int* row = colind+ncol+1;

    for (casadi_int cc=0; cc<ncol; ++cc) {
        for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
            *coord_i++ = row[el];
            *coord_j++ = cc;
        }
    }
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

// SYMBOL "madnlp_full_eval_constr_jac"
template<typename T1>
madnlp_int casadi_madnlp_full_eval_constr_jac(const T1* w, T1* res, void* user_data) {
  casadi_int i;
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = d->g;
  d_oracle->res[1] = res;
  calc_function(&d->prob->nlp_jac_g, d_oracle);
  //store res madnlp_data to be accessible from casadi
  d->jac_g = res;
  d->x =  (double*)w;
  return 0;
}

// SYMBOL "madnlp_full_eval_constr"
template<typename T1>
madnlp_int casadi_madnlp_full_eval_constr(const T1* w, T1* res, void* user_data) {
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

  return 1; //skip
}

// SYMBOL "madnlp_full_eval_obj_grad"
template<typename T1>
madnlp_int casadi_madnlp_full_eval_obj_grad(const T1* w, T1* res, void* user_data) {
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
  return 1; // skip
}

// SYMBOL "madnlp_full_eval_obj"
template<typename T1>
madnlp_int casadi_madnlp_full_eval_obj(const T1* w, T1* res, void* user_data) {
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  // TODO add argument
  T1 objective_scale = 1.0;

  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_f, d_oracle);
  //*res *= objective_scale;

  d->x =  (double*)w;
  d->f = res;

  return 1; // skip
}

// SYMBOL "madnlp_full_eval_obj"
template<typename T1>
madnlp_int casadi_madnlp_full_eval_lag_hess(T1 objective_scale, const T1* w, const T1* lam,
                                            T1* res, void* user_data){
  casadi_int k, column, i;
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  // evaluate hessian with given parameters / inputs
  //d_oracle->arg[0] = d->x;
  //d_oracle->arg[3] = d->lam;
  d_oracle->arg[0] = w;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->arg[2] = &objective_scale;
  d_oracle->arg[3] = lam;
  //d_oracle->res[0] = d->grad_l;
  d_oracle->res[0] = res; //hess_l
  calc_function(&d->prob->nlp_hess_l, d_oracle);

  //d->x = (double*)w;
  //d->lam_g = (double*)lam;
  //d->hess_l = res;

  return 0;
}

// C-REPLACE "const_cast<T1*>" "(T1*)"


// SYMBOL "madnlp_get_bounds"
template<typename T1>
madnlp_int casadi_madnlp_get_bounds(double *lower, double *upper, const madnlp_int k, void* user_data) {
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  casadi_int nx = p->nlp->nx;

  int i=0;
  int column = 0;
  return 0;
}

// SYMBOL "madnlp_get_initial_xk"
template<typename T1>
madnlp_int casadi_madnlp_get_initial_wk(double *xk, const madnlp_int k, void* user_data) {
  casadi_madnlp_data<T1>* d = static_cast< casadi_madnlp_data<T1>* >(user_data);
  const casadi_madnlp_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  //printf("casadi_madnlp_get_initial_xk offset %lld %e\n",p->CD[k].offset_c,d_nlp->z[p->CD[k].offset_c]);
  casadi_copy(d_nlp->z+p->CD[k].offset_c, p->nx[k], xk);

  return 0;
}


// C-REPLACE "casadi_madnlp_get_nx<T1>" "casadi_madnlp_get_nx"
// C-REPLACE "casadi_madnlp_get_nu<T1>" "casadi_madnlp_get_nu"
// C-REPLACE "casadi_madnlp_get_ng<T1>" "casadi_madnlp_get_ng"
// C-REPLACE "casadi_madnlp_get_ng_ineq<T1>" "casadi_madnlp_get_ng_ineq"
// C-REPLACE "casadi_madnlp_get_horizon_length<T1>" "casadi_madnlp_get_horizon_length"
// C-REPLACE "casadi_madnlp_get_bounds<T1>" "casadi_madnlp_get_bounds"
// C-REPLACE "casadi_madnlp_get_initial_wk<T1>" "casadi_madnlp_get_initial_wk"
// C-REPLACE "casadi_madnlp_get_initial_uk<T1>" "casadi_madnlp_get_initial_uk"
// C-REPLACE "casadi_madnlp_full_eval_constr_jac<T1>" "casadi_madnlp_full_eval_constr_jac"
// C-REPLACE "casadi_madnlp_full_eval_obj_grad<T1>" "casadi_madnlp_full_eval_obj_grad"
// C-REPLACE "casadi_madnlp_full_eval_obj<T1>" "casadi_madnlp_full_eval_obj"
// C-REPLACE "casadi_madnlp_full_eval_constr<T1>" "casadi_madnlp_full_eval_constr"
// C-REPLACE "casadi_madnlp_full_eval_lag_hess<T1>" "casadi_madnlp_full_eval_lag_hess"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "madnlp_work"
template<typename T1>
void casadi_madnlp_work(const casadi_madnlp_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->nlp->nx+p->nlp->ng)); // pv

  // Persistent work vectors
  *sz_w += p->nnz_hess_l;
  *sz_w += p->nnz_jac_g;
  *sz_w += p->nnz_grad_f;
}

// SYMBOL "madnlp_init"
template<typename T1>
void casadi_madnlp_init(casadi_madnlp_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Problem structure
  const casadi_madnlp_prob<T1>* p = d->prob;
  //casadi_oracle_data<T1>* d_oracle = d->nlp->oracle;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;

  //d->lam = *w; *w += p->nlp->nx+p->nlp->ng;
  //d->a = *w; *w += casadi_sp_nnz(p->sp_a);
  //d->h = *w; *w += casadi_sp_nnz(p->sp_h);
  //d->g = *w; *w += casadi_max(p->nlp->nx,p->nlp->ng);

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

  c_interface->full_eval_constr_jac = casadi_madnlp_full_eval_constr_jac<T1>;
  c_interface->full_eval_obj_grad = casadi_madnlp_full_eval_obj_grad<T1>;
  c_interface->full_eval_obj = casadi_madnlp_full_eval_obj<T1>;
  c_interface->full_eval_constr = casadi_madnlp_full_eval_constr<T1>;
  c_interface->full_eval_lag_hess = casadi_madnlp_full_eval_lag_hess<T1>;

  c_interface->user_data = d;

  d->solver = madnlp_c_create(&d->c_interface);

  // set number of variables and  constraints
  // reading from parent nlp interface
  d->solver->nw = p->nlp->nx;
  d->solver->nc = p->nlp->ng;
  // set nnz
  d->solver->nnzo = d->prob->nnz_grad_f;
  d->solver->nnzj = d->prob->nnz_jac_g;
  d->solver->nnzh = d->prob->nnz_hess_l;

  // allocate memory
  // CoordinateSparsity vector indexes data
  // output data memory
  madnlp_c_init(d->solver);

  // TODO: Fix hack to force template resolution
  double a = 1.0;
  // transform ccs sparsity to cco
  casadi_madnlp_sparsity(d->prob->jac_g_ccs, d->solver->nzj_i, d->solver->nzj_j,a);
  casadi_madnlp_sparsity(d->prob->hess_l_ccs, d->solver->nzh_i, d->solver->nzh_j,a);
  // update sparsity indexes to 1-base
  madnlp_c_update_cco_indexes(d->solver);

  // set number bounds
  // reading from parent nlp interface
  d->solver->lbx = d->nlp->lbz;
  d->solver->ubx = d->nlp->ubz;
  d->solver->lbg = d->nlp->lbz+p->nlp->nx;
  d->solver->ubg = d->nlp->ubz+p->nlp->nx;
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

  // set initial guess
  d->solver->x0 = d_nlp->z;
  d->solver->l0 = d_nlp->lam + p_nlp->nx;

  madnlp_int ret = madnlp_c_solve(d->solver);

  d_nlp->objective = d->solver->obj[0];

  for (casadi_int i=0; i<p_nlp->nx; ++i) {
    // get primal solution
    d_nlp->z[i] = d->solver->sol[i];
    // get dual solution x
    d_nlp->lam[i] = d->solver->mul_U[i]-d->solver->mul_L[i];
  }
  for (casadi_int i=0; i<p_nlp->ng; ++i) {
    d_nlp->z[p_nlp->nx+i] = d->solver->con[i];
    d_nlp->lam[p_nlp->nx+i] = d->solver->mul[i];
  }


  if (ret==0) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
    d->success = 1;
  }

  const struct MadnlpCStats* stats = madnlp_c_get_stats(d->solver);

  d->stats.compute_sd_time = stats->compute_sd_time;
  d->stats.duinf_time = stats->duinf_time;
  d->stats.eval_hess_time = stats->eval_hess_time;
  d->stats.eval_jac_time = stats->eval_jac_time;
  d->stats.eval_cv_time = stats->eval_cv_time;
  d->stats.eval_grad_time = stats->eval_grad_time;
  d->stats.eval_obj_time = stats->eval_obj_time;
  d->stats.initialization_time = stats->initialization_time;
  d->stats.time_total = stats->time_total;
  d->stats.eval_hess_count = stats->eval_hess_count;
  d->stats.eval_jac_count = stats->eval_jac_count;
  d->stats.eval_cv_count = stats->eval_cv_count;
  d->stats.eval_grad_count = stats->eval_grad_count;
  d->stats.eval_obj_count = stats->eval_obj_count;
  d->stats.iterations_count = stats->iterations_count;
  d->stats.return_flag = stats->return_flag;

  // Unpack dual solution

  madnlp_c_destroy(d->solver);
}

//#ifdef __cplusplus
//}
//#endif

