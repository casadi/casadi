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

// SYMBOL "fatrop_mproject"
template<typename T1>
void casadi_fatrop_mproject(T1 factor, const T1* x, const casadi_int* sp_x,
                            T1* y, const casadi_int* sp_y, T1* w) {
    casadi_int ncol_y;
    const casadi_int* colind_y;
    ncol_y = sp_y[1];
    colind_y = sp_y+2;
    casadi_project(x, sp_x, y, sp_y, w);
    casadi_scal(colind_y[ncol_y], factor, y);
}

// SYMBOL "fatrop_dense_transfer"
template<typename T1>
void casadi_fatrop_dense_transfer(double factor, const T1* x,
                                    const casadi_int* sp_x, T1* y,
                                    const casadi_int* sp_y, T1* w) {
    casadi_sparsify(x, w, sp_x, 0);
    casadi_int nrow_y = sp_y[0];
    casadi_int ncol_y = sp_y[1];
    const casadi_int *colind_y = sp_y+2, *row_y = sp_y + 2 + ncol_y+1;
    /* Loop over columns of y */
    casadi_int i, el;
    for (i=0; i<ncol_y; ++i) {
        for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[nrow_y*i + row_y[el]] += factor*(*w++);
    }
}

// SYMBOL "fatrop_read_primal_data"
template<typename T1>
void casadi_fatrop_read_primal_data(const double* primal_data, T1* x, const struct FatropOcpCDims *s) {
  int k;
  for (k=0;k<s->K;++k) {
    casadi_copy(primal_data+s->ux_offs[k], s->nu[k], x+s->nx[k]+s->ux_offs[k]);
    casadi_copy(primal_data+s->ux_offs[k]+s->nu[k], s->nx[k], x+s->ux_offs[k]);
  }
}

// SYMBOL "fatrop_write_primal_data"
template<typename T1>
void casadi_fatrop_write_primal_data(const double* x, T1* primal_data, const struct FatropOcpCDims *s) {
  int k;
  for (k=0;k<s->K;++k) {
    casadi_copy(x+s->nx[k]+s->ux_offs[k], s->nu[k], primal_data+s->ux_offs[k]);
    casadi_copy(x+s->ux_offs[k], s->nx[k], primal_data+s->ux_offs[k]+s->nu[k]);
  }
}

// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"

// C-REPLACE "OracleCallback" "struct casadi_oracle_callback"
template<typename T1>
struct casadi_fatrop_prob {
  const casadi_nlpsol_prob<T1>* nlp;
  const casadi_int *nx, *nu, *ng;
  casadi_int nx_max, nu_max, nxu_max;
  // Sparsity patterns
  const casadi_int *sp_h, *sp_a;

  casadi_int nnz_h, nnz_a;

  // Sparsities
  const casadi_int *ABsp;
  const casadi_int *AB_offsets;
  const casadi_int *CDsp;
  const casadi_int *CD_offsets;
  const casadi_int *RSQsp;
  const casadi_int *RSQ_offsets;
  const casadi_int *Isp;
  const casadi_int *I_offsets;

  casadi_int N;
  casadi_ocp_block *AB, *CD, *RSQ, *I;

  OracleCallback nlp_hess_l;
  OracleCallback nlp_jac_g;
  OracleCallback nlp_grad_f;
  OracleCallback nlp_f;
  OracleCallback nlp_g;

  FatropOcpCWrite write;
  FatropOcpCFlush flush;
};
// C-REPLACE "casadi_fatrop_prob<T1>" "struct casadi_fatrop_prob"


// SYMBOL "fatrop_setup"
template<typename T1>
void casadi_fatrop_setup(casadi_fatrop_prob<T1>* p) {
  casadi_int k;
  if (p->sp_h) {
    p->nnz_h = p->sp_h[2+p->sp_h[1]];
  } else {
    p->nnz_h = 0;
  }
  p->nnz_a = p->sp_a[2+p->sp_a[1]];

  p->nx_max = 0;
  for (k=0;k<p->N+1;++k) {
    if (p->nx[k]>p->nx_max) p->nx_max = p->nx[k];
  }
  p->nu_max = 0;
  p->nxu_max = 0;
  for (k=0;k<p->N;++k) {
    if (p->nu[k]>p->nu_max) p->nu_max = p->nu[k];
    if (p->nu[k]+p->nx[k]>p->nxu_max) p->nxu_max = p->nu[k]+p->nx[k];
  }
}



// SYMBOL "fatrop_data"
template<typename T1>
struct casadi_fatrop_data {
  // Problem structure
  const casadi_fatrop_prob<T1>* prob;
  // Problem structure
  casadi_nlpsol_data<T1>* nlp;

  T1 *AB, *CD, *RSQ, *I;

  casadi_int *a_eq, *a_ineq, *a_eq_idx, *a_ineq_idx;
  casadi_int *x_eq, *x_ineq, *x_eq_idx, *x_ineq_idx;

  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;

  int unified_return_status;
  int success;

  T1 *pv, *x, *a, *g, *h, *lam;

  struct blasfeo_dvec v, r;
  struct blasfeo_dmat R;

  struct FatropOcpCInterface ocp_interface;

  struct FatropOcpCStats stats;

  struct FatropOcpCSolver *solver;
};
// C-REPLACE "casadi_fatrop_data<T1>" "struct casadi_fatrop_data"

// SYMBOL "fatrop_init_mem"
template<typename T1>
int fatrop_init_mem(casadi_fatrop_data<T1>* d) {
  return 0;
}

// SYMBOL "fatrop_free_mem"
template<typename T1>
void fatrop_free_mem(casadi_fatrop_data<T1>* d) {
  //Highs_destroy(d->fatrop);
  
}
// C-REPLACE "static_cast< casadi_fatrop_data<T1>* >" "(struct casadi_fatrop_data*)"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"
// C-REPLACE "calc_function" "casadi_oracle_call"
// C-REPLACE "casadi_error" "//casadi_error"

// SYMBOL "fatrop_full_eval_constr_jac"
template<typename T1>
fatrop_int casadi_fatrop_full_eval_constr_jac(const double* primal_data, const double* stageparams_p, const double* globalparams_p,
            struct blasfeo_dmat* BAbt_p, struct blasfeo_dmat* Ggt_p, struct blasfeo_dmat* Ggt_ineq_p, const struct FatropOcpCDims* s, void* user_data) {
  casadi_int i;
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  casadi_fatrop_read_primal_data(primal_data, d->x, s);
  d_oracle->arg[0] = d->x;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = d->g;
  d_oracle->res[1] = d->a;
  calc_function(&d->prob->nlp_jac_g, d_oracle);

  casadi_fatrop_mproject(-1.0, d->a, p->sp_a, d->AB, p->ABsp, d->pv);
  casadi_project(d->a, p->sp_a, d->CD, p->CDsp, d->pv);
  casadi_project(d->a, p->sp_a, d->I, p->Isp, d->pv);

  for (i=0;i<p->Isp[2+p->Isp[1]];++i) {
    if (d->I[0]!=1.0) {
      casadi_error("Structure mismatch: gap-closing constraints must be like this: x_{k+1}-F(xk,uk).");
    }
  }
  return 0;
}

// SYMBOL "fatrop_full_eval_contr_viol"
template<typename T1>
fatrop_int casadi_fatrop_full_eval_contr_viol(const double* primal_data, const double* stageparams_p, const double* globalparams_p,
            double* res, const struct FatropOcpCDims* s, void* user_data) {
  casadi_int i,k,column;
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  casadi_fatrop_read_primal_data(primal_data, d->x, s);
  d_oracle->arg[0] = d->x;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = d->g;
  calc_function(&d->prob->nlp_g, d_oracle);

  for (k=0;k<s->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      res[s->g_ineq_offs[k]+column] = d->g[d->a_ineq[i]];
      column++;
    }
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      res[s->g_ineq_offs[k]+column] = d->x[d->x_ineq[i]];
      column++;
    }
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      res[s->g_offs[k]+column] = d->g[d->a_eq[i]]-d->nlp->lbz[p->nlp->nx+d->a_eq[i]];
      column++;
    }
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      res[s->g_offs[k]+column] = d->x[d->x_eq[i]]-d->nlp->lbz[d->x_eq[i]];
      column++;
    }
  }
  for (k=0;k<s->K-1;++k) {
    casadi_scaled_copy(-1.0, d->g+p->AB[k].offset_r, p->nx[k+1], res+s->dyn_eq_offs[k]);
  }
  return 1; //skip
}

// SYMBOL "fatrop_full_eval_obj_grad"
template<typename T1>
fatrop_int casadi_fatrop_full_eval_obj_grad(
            double objective_scale,
            const double *primal_data,
            const double *stage_params_k,
            const double *global_params,
            double *res, const struct FatropOcpCDims* s, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  casadi_fatrop_read_primal_data(primal_data, d->x, s);
  d_oracle->arg[0] = d->x;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = d->g;
  calc_function(&d->prob->nlp_grad_f, d_oracle);

  casadi_fatrop_write_primal_data(d->g, res, s);
  casadi_scal(p->nlp->nx, objective_scale, res);
  return 1; // skip
}

// SYMBOL "fatrop_full_eval_obj"
template<typename T1>
fatrop_int casadi_fatrop_full_eval_obj(
            double objective_scale,
            const double *primal_data,
            const double *stage_params_k,
            const double *global_params,
            double *res, const struct FatropOcpCDims* s, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  casadi_fatrop_read_primal_data(primal_data, d->x, s);
  d_oracle->arg[0] = d->x;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->res[0] = res;
  calc_function(&d->prob->nlp_f, d_oracle);

  *res *= objective_scale;
  return 1; // skip
}

// SYMBOL "fatrop_full_eval_obj"
template<typename T1>
fatrop_int casadi_fatrop_full_eval_lag_hess(
            double objective_scale,
            const double *primal_data,
            const double *lam_data,
            const double *stage_params_k,
            const double *global_params,
            struct blasfeo_dmat *res, const struct FatropOcpCDims* s, void* user_data) {
  casadi_int k, column, i;
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_oracle_data<T1>* d_oracle = d_nlp->oracle;

  casadi_fatrop_read_primal_data(primal_data, d->x, s);
  for (k=0;k<s->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      d->lam[d->a_ineq[i]] = lam_data[s->g_ineq_offs[k]+column];
      column++;
    }
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      d->lam[d->a_eq[i]] = lam_data[s->g_offs[k]+column];
      column++;
    }
  }
  for (k=0;k<s->K-1;++k) {
    casadi_scaled_copy(-1.0, lam_data+s->dyn_eq_offs[k], p->nx[k+1], d->lam+p->AB[k].offset_r);
  }

  d_oracle->arg[0] = d->x;
  d_oracle->arg[1] = d_nlp->p;
  d_oracle->arg[2] = &objective_scale;
  d_oracle->arg[3] = d->lam;
  d_oracle->res[0] = d->g;
  d_oracle->res[1] = d->h;
  calc_function(&d->prob->nlp_hess_l, d_oracle);

  casadi_project(d->h, p->sp_h, d->RSQ, p->RSQsp, d->pv);

  // Note: Lagrangian is defined with lambda_dyn*(F(x_k,u_k))  with x_k+1 notably absent
  for (k=0;k<s->K-1;++k) {
    casadi_axpy(p->nx[k+1], 1.0, lam_data+s->dyn_eq_offs[k], d->g+p->CD[k+1].offset_c);
  }
  // Note: we still need to handle simple bounds
  for (k=0;k<s->K;++k) {
    column = d->a_ineq_idx[k+1]-d->a_ineq_idx[k];
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      d->g[d->x_ineq[i]] += lam_data[s->g_ineq_offs[k]+column];
      column++;
    }
    column = d->a_eq_idx[k+1]-d->a_eq_idx[k];
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      d->g[d->x_eq[i]] += lam_data[s->g_offs[k]+column];
      column++;
    }
  }

  return 0;
}

// C-REPLACE "const_cast<T1*>" "(T1*)"

// SYMBOL "fatrop_eval_BAbt"
template<typename T1>
fatrop_int casadi_fatrop_eval_BAbt(const double *states_kp1, const double *inputs_k,
    const double *states_k, const double *stage_params_k,
    const double *global_params, struct blasfeo_dmat *res, const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  blasfeo_pack_tran_dmat(p->nx[k+1], p->nx[k], d->AB+p->AB_offsets[k], p->nx[k+1], res, p->nu[k], 0);
  blasfeo_pack_tran_dmat(p->nx[k+1], p->nu[k], d->AB+p->AB_offsets[k]+p->nx[k]*p->nx[k+1], p->nx[k+1], res, 0, 0);
  blasfeo_pack_dmat(1, p->nx[k+1], const_cast<T1*>(d_nlp->lbz+p->nlp->nx+p->AB[k].offset_r), 1, res, p->nx[k]+p->nu[k], 0);

  casadi_scal(p->nx[k+1], -1.0, d->g+p->AB[k].offset_r);

  blasfeo_pack_dmat(1, p->nx[k+1], d->g+p->AB[k].offset_r, 1, res, p->nx[k]+p->nu[k], 0);

  return 0;

}

// SYMBOL "fatrop_eval_RSQrqt"
template<typename T1>
fatrop_int casadi_fatrop_eval_RSQrqt(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *lam_dyn_k,
    const double *lam_eq_k,
    const double *lam_eq_ineq_k,
    const double *stage_params_k,
    const double *global_params,
    struct blasfeo_dmat *res,
    const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;

  int n = p->nx[k]+p->nu[k];
  blasfeo_pack_dmat(p->nx[k], p->nx[k],
    d->RSQ+p->RSQ_offsets[k], n, res, p->nu[k], p->nu[k]);
  blasfeo_pack_dmat(p->nu[k], p->nu[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k]*n+p->nx[k], n, res, 0, 0);
  blasfeo_pack_dmat(p->nu[k], p->nx[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k], n, res, 0, p->nu[k]);
  blasfeo_pack_dmat(p->nx[k], p->nu[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k]*n, n, res, p->nu[k], 0);


  blasfeo_pack_dmat(1, p->nx[k], d->g+p->CD[k].offset_c, 1, res, p->nx[k]+p->nu[k], p->nu[k]);
  blasfeo_pack_dmat(1, p->nu[k], d->g+p->CD[k].offset_c+p->nx[k], 1, res, p->nx[k]+p->nu[k], 0);

  return 0;
}


// SYMBOL "fatrop_eval_Ggt"
template<typename T1>
fatrop_int casadi_fatrop_eval_Ggt(
      const double *inputs_k,
      const double *states_k,
      const double *stage_params_k,
      const double *global_params,
      struct blasfeo_dmat *res,
      const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_int i, column;

  int n_a_eq = d->a_eq_idx[k+1]-d->a_eq_idx[k];
  int n_x_eq = d->x_eq_idx[k+1]-d->x_eq_idx[k];
  int ng_eq = n_a_eq+n_x_eq;

  blasfeo_dgese(p->nx[k]+p->nu[k]+1, ng_eq, 0.0, res, 0, 0);

  column = 0;
  for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
    blasfeo_pack_tran_dmat(1, p->nx[k],
      d->CD+p->CD_offsets[k]+(d->a_eq[i]-p->CD[k].offset_r),
      p->CD[k].rows, res, p->nu[k], column);
    blasfeo_pack_tran_dmat(1, p->nu[k],
      d->CD+p->CD_offsets[k]+(d->a_eq[i]-p->CD[k].offset_r)+p->nx[k]*p->CD[k].rows,
      p->CD[k].rows, res, 0,        column);
    BLASFEO_DMATEL(res, p->nx[k]+p->nu[k], column) = d->g[d->a_eq[i]]-d->nlp->lbz[p->nlp->nx+d->a_eq[i]];
    column++;
  }
  for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
    int j = d->x_eq[i]-p->CD[k].offset_c;
    if (j>=p->nx[k]) {
      j -= p->nx[k];
    } else {
      j += p->nu[k];
    }
    BLASFEO_DMATEL(res, j, column) = 1;
    BLASFEO_DMATEL(res, p->nx[k]+p->nu[k], column) = d->x[d->x_eq[i]]-d->nlp->lbz[d->x_eq[i]];
    column++;
  }

  return 0;
}

// SYMBOL "fatrop_eval_Ggt_ineq"
template<typename T1>
fatrop_int  casadi_fatrop_eval_Ggt_ineq(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    struct blasfeo_dmat *res,
    const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_int i, column;

  int n_a_ineq = d->a_ineq_idx[k+1]-d->a_ineq_idx[k];
  int n_x_ineq = d->x_ineq_idx[k+1]-d->x_ineq_idx[k];
  int ng_ineq = n_a_ineq+n_x_ineq;

  // Ggt_ineq: [G_ineq;g_ineq]
  // G_ineq: (nu+nx by ng_ineq)
  // g_ineq: (nu+nx by 1)

  // Clear Ggt_ineq
  blasfeo_dgese(p->nx[k]+p->nu[k]+1, ng_ineq, 0.0, res, 0, 0);

  column = 0;
  for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
    blasfeo_pack_tran_dmat(1, p->nx[k],
      d->CD+p->CD_offsets[k]+(d->a_ineq[i]-p->CD[k].offset_r),
      p->CD[k].rows, res, p->nu[k], column);
    blasfeo_pack_tran_dmat(1, p->nu[k],
      d->CD+p->CD_offsets[k]+(d->a_ineq[i]-p->CD[k].offset_r)+p->nx[k]*p->CD[k].rows,
      p->CD[k].rows, res, 0,        column);
    BLASFEO_DMATEL(res, p->nx[k]+p->nu[k], column) = d->g[d->a_ineq[i]];
    column++;
  }
  for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
    int j = d->x_ineq[i]-p->CD[k].offset_c;
    if (j>=p->nx[k]) {
      j -= p->nx[k];
    } else {
      j += p->nu[k];
    }
    BLASFEO_DMATEL(res, j, column) = 1;
    BLASFEO_DMATEL(res, p->nx[k]+p->nu[k], column) = d->x[d->x_ineq[i]];
    column++;
  }

  return 0;
}

// SYMBOL "fatrop_get_nx"
template<typename T1>
fatrop_int casadi_fatrop_get_nx(const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;

  //printf("nx %lld\n", p->nx[k]);

  if (k==p->N+1) return p->nx[k-1];
  return p->nx[k];
}

// SYMBOL "fatrop_get_nu"
template<typename T1>
fatrop_int casadi_fatrop_get_nu(const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  return p->nu[k];
}

// SYMBOL "fatrop_get_ng"
template<typename T1>
fatrop_int casadi_fatrop_get_ng(const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  int ret;
  fatrop_int n_a_eq = d->a_eq_idx[k+1]-d->a_eq_idx[k];
  fatrop_int n_x_eq = d->x_eq_idx[k+1]-d->x_eq_idx[k];

  ret = n_a_eq+n_x_eq;

 /* printf("d->a_eq_idx[k] %lld\n", d->a_eq_idx[k]);
  printf("d->a_eq_idx[k+1] %lld\n", d->a_eq_idx[k+1]);
  printf("d->x_eq_idx[k] %lld\n", d->x_eq_idx[k]);
  printf("d->x_eq_idx[k+1] %lld\n", d->x_eq_idx[k+1]);


  printf("get_ng %d\n", ret);*/
  return ret;
}


// SYMBOL "fatrop_get_horizon_length"
template<typename T1>
fatrop_int casadi_fatrop_get_horizon_length(void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);

  //printf("horizon_length %lld\n", d->prob->N+1);
  return d->prob->N+1;
}

// SYMBOL "fatrop_get_ng_ineq"
template<typename T1>
fatrop_int casadi_fatrop_get_ng_ineq(const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  fatrop_int n_a_ineq = d->a_ineq_idx[k+1]-d->a_ineq_idx[k];
  fatrop_int n_x_ineq = d->x_ineq_idx[k+1]-d->x_ineq_idx[k];
  //printf("get_ng_ineq %d\n", n_a_ineq+n_x_ineq);
  return n_a_ineq+n_x_ineq;
}

// SYMBOL "fatrop_get_bounds"
template<typename T1>
fatrop_int casadi_fatrop_get_bounds(double *lower, double *upper, const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  casadi_int nx = p->nlp->nx;

  int i=0;
  int column = 0;
  for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
    lower[column] = d_nlp->lbz[nx+d->a_ineq[i]];
    upper[column] = d_nlp->ubz[nx+d->a_ineq[i]];
    column++;
  }

  for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
    lower[column] = d_nlp->lbz[d->x_ineq[i]];
    upper[column] = d_nlp->ubz[d->x_ineq[i]];
    column++;
  }

  return 0;
}

// SYMBOL "fatrop_get_initial_xk"
template<typename T1>
fatrop_int casadi_fatrop_get_initial_xk(double *xk, const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  //printf("casadi_fatrop_get_initial_xk offset %lld %e\n",p->CD[k].offset_c,d_nlp->z[p->CD[k].offset_c]);
  casadi_copy(d_nlp->z+p->CD[k].offset_c, p->nx[k], xk);

  return 0;
}

// SYMBOL "fatrop_get_initial_uk"
template<typename T1>
fatrop_int casadi_fatrop_get_initial_uk(double *uk, const fatrop_int k, void* user_data) {
  casadi_fatrop_data<T1>* d = static_cast< casadi_fatrop_data<T1>* >(user_data);
  const casadi_fatrop_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_copy(d_nlp->z+p->CD[k].offset_c+p->nx[k], p->nu[k], uk);
  return 0;
}


// SYMBOL "fatrop_work"
template<typename T1>
void casadi_fatrop_work(const casadi_fatrop_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->nlp->nx+p->nlp->ng)); // pv

  // Persistent work vectors
  *sz_w += casadi_sp_nnz(p->ABsp); // AB
  *sz_w += casadi_sp_nnz(p->CDsp); // CD
  *sz_w += casadi_sp_nnz(p->RSQsp); // RSQ
  *sz_w += casadi_sp_nnz(p->Isp); // I
  *sz_w += p->nlp->nx;  // x
  *sz_w += p->nlp->nx+p->nlp->ng; // lam
  *sz_w += casadi_sp_nnz(p->sp_a);  // a
  *sz_w += casadi_sp_nnz(p->sp_h);  // h
  *sz_w += casadi_max(p->nlp->nx,p->nlp->ng);  // g
  *sz_w += blasfeo_memsize_dvec(p->nxu_max+1)+64; // v p->nx[k]+p->nu[k]+1
  *sz_w += blasfeo_memsize_dvec(p->nx_max+p->nlp->ng)+64; // r p->nx[k+1]
  *sz_w += blasfeo_memsize_dmat(p->nxu_max, p->nxu_max)+64; // r p->nx[k]+p->nx[u]

  *sz_iw += p->N+2; // a_eq_idx
  *sz_iw += p->N+2; // a_ineq_idx
  *sz_iw += p->N+2; // x_eq_idx
  *sz_iw += p->N+2; // x_ineq_idx
  *sz_iw += p->nlp->ng; // a_eq
  *sz_iw += p->nlp->ng; // a_ineq
  *sz_iw += p->nlp->nx; // x_eq
  *sz_iw += p->nlp->nx; // x_ineq

}

// SYMBOL "fatrop_init"
template<typename T1>
void casadi_fatrop_init(casadi_fatrop_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Problem structure
  const casadi_fatrop_prob<T1>* p = d->prob;
  //casadi_oracle_data<T1>* d_oracle = d->nlp->oracle;
  //const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  
  //d->z_L = *w; *w += p_nlp->nx;
  //d->z_U = *w; *w += p_nlp->nx;

  d->AB = *w; *w += casadi_sp_nnz(p->ABsp);
  d->CD = *w; *w += casadi_sp_nnz(p->CDsp);
  d->RSQ = *w; *w += casadi_sp_nnz(p->RSQsp);
  d->I = *w; *w += casadi_sp_nnz(p->Isp);
  d->x = *w; *w += p->nlp->nx;
  d->lam = *w; *w += p->nlp->nx+p->nlp->ng;
  d->a = *w; *w += casadi_sp_nnz(p->sp_a);
  d->h = *w; *w += casadi_sp_nnz(p->sp_h);
  d->g = *w; *w += casadi_max(p->nlp->nx,p->nlp->ng);
  blasfeo_create_dvec(p->nxu_max+1, &d->v, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dvec(p->nxu_max+1)+64;
  blasfeo_create_dvec(p->nx_max+p->nlp->ng, &d->r, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dvec(p->nx_max+p->nlp->ng)+64;
  blasfeo_create_dmat(p->nxu_max, p->nxu_max, &d->R, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dmat(p->nxu_max, p->nxu_max)+64;
  
  d->a_eq_idx = *iw;   *iw += p->N+2;
  d->a_ineq_idx = *iw; *iw += p->N+2;
  d->x_eq_idx = *iw;   *iw += p->N+2;
  d->x_ineq_idx = *iw; *iw += p->N+2;
  
  d->a_eq = *iw;   *iw += p->nlp->ng;
  d->a_ineq = *iw; *iw += p->nlp->ng;
  d->x_eq = *iw;  *iw += p->nlp->nx;
  d->x_ineq = *iw; *iw += p->nlp->nx;

  d->pv = *w;

  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// C-REPLACE "casadi_fatrop_get_nx<T1>" "casadi_fatrop_get_nx"
// C-REPLACE "casadi_fatrop_get_nu<T1>" "casadi_fatrop_get_nu"
// C-REPLACE "casadi_fatrop_get_ng<T1>" "casadi_fatrop_get_ng"
// C-REPLACE "casadi_fatrop_get_ng_ineq<T1>" "casadi_fatrop_get_ng_ineq"
// C-REPLACE "casadi_fatrop_get_horizon_length<T1>" "casadi_fatrop_get_horizon_length"
// C-REPLACE "casadi_fatrop_get_bounds<T1>" "casadi_fatrop_get_bounds"
// C-REPLACE "casadi_fatrop_get_initial_xk<T1>" "casadi_fatrop_get_initial_xk"
// C-REPLACE "casadi_fatrop_get_initial_uk<T1>" "casadi_fatrop_get_initial_uk"
// C-REPLACE "casadi_fatrop_full_eval_constr_jac<T1>" "casadi_fatrop_full_eval_constr_jac"
// C-REPLACE "casadi_fatrop_full_eval_obj_grad<T1>" "casadi_fatrop_full_eval_obj_grad"
// C-REPLACE "casadi_fatrop_full_eval_obj<T1>" "casadi_fatrop_full_eval_obj"
// C-REPLACE "casadi_fatrop_full_eval_contr_viol<T1>" "casadi_fatrop_full_eval_contr_viol"
// C-REPLACE "casadi_fatrop_full_eval_lag_hess<T1>" "casadi_fatrop_full_eval_lag_hess"
// C-REPLACE "casadi_fatrop_eval_BAbt<T1>" "casadi_fatrop_eval_BAbt"
// C-REPLACE "casadi_fatrop_eval_RSQrqt<T1>" "casadi_fatrop_eval_RSQrqt"
// C-REPLACE "casadi_fatrop_eval_Ggt<T1>" "casadi_fatrop_eval_Ggt"
// C-REPLACE "casadi_fatrop_eval_Ggt_ineq<T1>" "casadi_fatrop_eval_Ggt_ineq"


// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "fatrop_presolve"
template<typename T1>
void casadi_fatrop_presolve(casadi_fatrop_data<T1>* d) {
  casadi_int k, i, start, stop, nx;
  // Problem structure
  const casadi_fatrop_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  struct FatropOcpCInterface* ocp_interface = &d->ocp_interface;

  ocp_interface->get_nx = casadi_fatrop_get_nx<T1>;
  ocp_interface->get_nu = casadi_fatrop_get_nu<T1>;
  ocp_interface->get_ng = casadi_fatrop_get_ng<T1>;
  ocp_interface->get_n_stage_params = 0;
  ocp_interface->get_n_global_params = 0;
  ocp_interface->get_default_stage_params = 0;
  ocp_interface->get_default_global_params = 0;
  ocp_interface->get_ng_ineq = casadi_fatrop_get_ng_ineq<T1>;
  ocp_interface->get_horizon_length = casadi_fatrop_get_horizon_length<T1>;

  ocp_interface->get_bounds = casadi_fatrop_get_bounds<T1>;
  ocp_interface->get_initial_xk = casadi_fatrop_get_initial_xk<T1>;
  ocp_interface->get_initial_uk = casadi_fatrop_get_initial_uk<T1>;

  ocp_interface->full_eval_constr_jac = casadi_fatrop_full_eval_constr_jac<T1>;
  ocp_interface->full_eval_obj_grad = casadi_fatrop_full_eval_obj_grad<T1>;
  ocp_interface->full_eval_obj = casadi_fatrop_full_eval_obj<T1>;
  ocp_interface->full_eval_contr_viol = casadi_fatrop_full_eval_contr_viol<T1>;
  ocp_interface->full_eval_lag_hess = casadi_fatrop_full_eval_lag_hess<T1>;

  ocp_interface->eval_BAbt = casadi_fatrop_eval_BAbt<T1>;
  ocp_interface->eval_RSQrqt = casadi_fatrop_eval_RSQrqt<T1>;
  ocp_interface->eval_Ggt = casadi_fatrop_eval_Ggt<T1>;
  ocp_interface->eval_Ggt_ineq = casadi_fatrop_eval_Ggt_ineq<T1>;
  ocp_interface->eval_rq = 0; // Computed in full_eval_obj_grad
  ocp_interface->eval_L = 0; // Computed in full_eval_obj

  nx = p_nlp->nx;

  d->a_eq_idx[0] = 0;
  d->a_ineq_idx[0] = 0;
  d->x_eq_idx[0] = 0;
  d->x_ineq_idx[0] = 0;

  // Loop over CD blocks
  for (k=0;k<p->N+1;++k) {
    d->a_eq_idx[k+1] = d->a_eq_idx[k];
    d->a_ineq_idx[k+1] = d->a_ineq_idx[k];
    start = p->CD[k].offset_r;
    stop  = p->CD[k].offset_r+p->CD[k].rows;
    for (i=start;i<stop;++i) {
      if (d_nlp->lbz[nx+i]==d_nlp->ubz[nx+i]) {
        d->a_eq[d->a_eq_idx[k+1]++] = i;
      } else {
        if (d_nlp->lbz[nx+i]==-std::numeric_limits<T1>::infinity() && d_nlp->ubz[nx+i]==std::numeric_limits<T1>::infinity()) continue;
        d->a_ineq[d->a_ineq_idx[k+1]++] = i;
      }
    }
    d->x_eq_idx[k+1] = d->x_eq_idx[k];
    d->x_ineq_idx[k+1] = d->x_ineq_idx[k];
    start = p->CD[k].offset_c;
    stop  = p->CD[k].offset_c+p->CD[k].cols;

    for (i=start;i<stop;++i) {
      if (d_nlp->lbz[i]==d_nlp->ubz[i]) {
        d->x_eq[d->x_eq_idx[k+1]++] = i;
      } else {
        if (d_nlp->lbz[i]==-std::numeric_limits<T1>::infinity() && d_nlp->ubz[i]==std::numeric_limits<T1>::infinity()) continue;
        d->x_ineq[d->x_ineq_idx[k+1]++] = i;
      }
    }
  }

  d->ocp_interface.user_data = d;

  d->solver = fatrop_ocp_c_create(&d->ocp_interface, p->write, p->flush);
}

// SYMBOL "fatrop_ocp_c_solve"
template<typename T1>
void casadi_fatrop_solve(casadi_fatrop_data<T1>* d) {
  // Problem structure
  casadi_int k, i, column;
  const casadi_fatrop_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  d->unified_return_status = SOLVER_RET_UNKNOWN;
  d->success = 0;

  fatrop_int ret = fatrop_ocp_c_solve(d->solver);

  if (ret==0) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
    d->success = 1;
  }

  const struct blasfeo_dvec* primal = fatrop_ocp_c_get_primal(d->solver);
  const struct blasfeo_dvec* dual = fatrop_ocp_c_get_dual(d->solver);
  const struct FatropOcpCDims* str = fatrop_ocp_c_get_dims(d->solver);
  const struct FatropOcpCStats* stats = fatrop_ocp_c_get_stats(d->solver);

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

  const double* primal_data = primal->pa;
  const double* dual_data = dual->pa;

  casadi_fatrop_read_primal_data(primal_data, d_nlp->z, str);
  // Unpack dual solution
  // Inequalities
  for (k=0;k<str->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      d_nlp->lam[p_nlp->nx+d->a_ineq[i]] = dual_data[str->g_ineq_offs[k]+column];
      column++;
    }
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      d_nlp->lam[d->x_ineq[i]] = dual_data[str->g_ineq_offs[k]+column];
      column++;
    }
  }
  // Equalities
  for (k=0;k<str->K;++k) {
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      d_nlp->lam[p_nlp->nx+d->a_eq[i]] = dual_data[str->g_offs[k]+column];
      column++;
    }
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      d_nlp->lam[d->x_eq[i]] = dual_data[str->g_offs[k]+column];
      column++;
    }
  }
  // Dynamics
  for (k=0;k<str->K-1;++k) {
    casadi_scaled_copy(-1.0, dual_data+str->dyn_eq_offs[k], p->nx[k+1], d_nlp->lam+p_nlp->nx+p->AB[k].offset_r);
  }

  fatrop_ocp_c_destroy(d->solver);

}