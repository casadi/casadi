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

// C-REPLACE "casadi_qp_prob<T1>" "struct casadi_qp_prob"
// C-REPLACE "casadi_qp_data<T1>" "struct casadi_qp_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "const_cast<int*>" "(int*) "


// SYMBOL "hpipm_mproject"
template<typename T1>
void casadi_hpipm_mproject(T1 factor, const T1* x, const casadi_int* sp_x,
                            T1* y, const casadi_int* sp_y, T1* w) {
    casadi_int ncol_y;
    const casadi_int* colind_y;
    ncol_y = sp_y[1];
    colind_y = sp_y+2;
    casadi_project(x, sp_x, y, sp_y, w);
    casadi_scal(colind_y[ncol_y], factor, y);
}

// SYMBOL "hpipm_dense_transfer"
template<typename T1>
void casadi_hpipm_dense_transfer(double factor, const T1* x,
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

// SYMBOL "hpipm_block"
struct casadi_hpipm_block {
    casadi_int offset_r;
    casadi_int offset_c;
    casadi_int rows;
    casadi_int cols;
};
// C-REPLACE "casadi_hpipm_block" "struct casadi_hpipm_block"

// SYMBOL "hpipm_unpack_blocks"
inline void casadi_hpipm_unpack_blocks(casadi_int N, casadi_hpipm_block* blocks, const casadi_int* packed) {
    casadi_int i;
    for (i=0;i<N;++i) {
        blocks[i].offset_r = *packed++;
        blocks[i].offset_c = *packed++;
        blocks[i].rows = *packed++;
        blocks[i].cols = *packed++;
    }
}

// SYMBOL "hpipm_ptr_block"
template<typename T1>
void casadi_hpipm_ptr_block(casadi_int N, T1** vs, T1* v, const casadi_hpipm_block* blocks, int eye) {
    casadi_int k, offset = 0;
    for(k=0;k<N;++k) {
        vs[k] = v+offset;
        if (eye) {
        offset += blocks[k].rows;
        } else {
        offset += blocks[k].rows*blocks[k].cols;
        }
    }
}

template<typename T1>
struct casadi_hpipm_prob {
  const casadi_qp_prob<T1>* qp;
  const int *nx, *nu, *ng;
  const int *nbx, *nbu, *ns;
  const int *nsbx, *nsbu, *nsg;
  // Sparsities
  const casadi_int *sp_x, *sp_ba;
  const casadi_int *Asp, *Bsp, *Csp, *Dsp;
  const casadi_int *Rsp, *Isp, *Ssp, *Qsp;
  const casadi_int *bsp;
  const casadi_int *xsp, *usp;
  const casadi_int *pisp;
  const casadi_int *theirs_xsp, *theirs_usp, *theirs_Xsp, *theirs_Usp;
  const casadi_int *lamg_gapsp, *lugsp;

  casadi_int N;
  casadi_int nx_total, nu_total, ng_total;
  casadi_hpipm_block *A, *B, *C, *D;
  casadi_hpipm_block *R, *I, *S, *Q;
  casadi_hpipm_block *b, *lug;
  casadi_hpipm_block *u, *x;
  casadi_hpipm_block *lam_ul, *lam_xl, *lam_uu, *lam_xu, *lam_cl, *lam_cu;

  T1 warm_start;
  T1 inf;

  struct d_ocp_qp_ipm_arg hpipm_options;

};
// C-REPLACE "casadi_hpipm_prob<T1>" "struct casadi_hpipm_prob"

// SYMBOL "hpipm_setup"
template<typename T1>
void casadi_hpipm_setup(casadi_hpipm_prob<T1>* p) {
  casadi_int i;

  p->nx_total = 0;
  p->nu_total = 0;
  p->ng_total = 0;
  for (i=0;i<p->N+1;++i) p->nx_total+= p->nx[i];
  for (i=0;i<p->N+1;++i) p->nu_total+= p->nu[i];
  for (i=0;i<p->N+1;++i) p->ng_total+= p->ng[i];

}



// SYMBOL "hpipm_data"
template<typename T1>
struct casadi_hpipm_data {
  // Problem structure
  const casadi_hpipm_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;
  T1 *A, *B, *C, *D;
  T1 *R, *I, *Q, *S;
  T1 *b, *b2;
  T1 *x, *q;
  T1 *u, *r;
  T1 *lg, *ug;
  T1 *pi;
  T1 *lbx, *ubx, *lbu, *ubu, *lam;

  T1 **hA, **hB, **hC, **hD;
  T1 **hR, **hI, **hQ, **hS;
  T1 **hx, **hq;
  T1 **hu, **hr;
  T1 **hlg, **hug;
  T1 **hb;

  T1 **hZl, **hZu, **hzl, **hzu, **hlls, **hlus;
  T1 **pis, **hlbx, **hubx, **hlbu, **hubu, **lams;

  int *iidxbx, *iidxbu;
  int **hidxbx, **hidxbu, **hidxs;

  int iter_count;
  int return_status;
  T1 res_stat;
  T1 res_eq;
  T1 res_ineq;
  T1 res_comp;

  T1 *pv;
};
// C-REPLACE "casadi_hpipm_data<T1>" "struct casadi_hpipm_data"


// SYMBOL "hpipm_work"
template<typename T1>
void casadi_hpipm_work(const casadi_hpipm_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->qp->nx+p->qp->na)); // pv
  // Persistent work vectors
  *sz_res += p->N; // hA
  *sz_res += p->N; // hB
  *sz_res += p->N+1; // hC
  *sz_res += p->N+1; // hD
  *sz_res += p->N+1; // hR
  *sz_res += p->N; // hI
  *sz_res += p->N+1; // hS
  *sz_res += p->N+1; // hQ
  *sz_res += p->N+1; // hx
  *sz_res += p->N+1; // hq
  *sz_res += p->N+1; // hu
  *sz_res += p->N+1; // hr
  *sz_res += p->N+1; // hlg
  *sz_res += p->N+1; // hug
  *sz_res += p->N; // hb

  *sz_res += p->N; // pis
  *sz_res += p->N+1; // hlbx
  *sz_res += p->N+1; // hubx
  *sz_res += p->N+1; // hlbu
  *sz_res += p->N+1; // hubu
  *sz_res += p->N+1; // lams
  *sz_res += p->N+1; // hidxbx
  *sz_res += p->N+1; // hidxbu
  *sz_res += p->N+1; // hidxs

  *sz_res += p->N+1; // hZl
  *sz_res += p->N+1; // hZu
  *sz_res += p->N+1; // hzl
  *sz_res += p->N+1; // hzu
  *sz_res += p->N+1; // hlls
  *sz_res += p->N+1; // hlus

  *sz_w += casadi_sp_nnz(p->Asp); // A
  *sz_w += casadi_sp_nnz(p->Bsp); // B
  *sz_w += casadi_sp_nnz(p->Csp); // C
  *sz_w += casadi_sp_nnz(p->Dsp); // D
  *sz_w += casadi_sp_nnz(p->Rsp); // R
  *sz_w += casadi_sp_nnz(p->Isp); // I
  *sz_w += casadi_sp_nnz(p->Ssp); // S
  *sz_w += casadi_sp_nnz(p->Qsp); // Q
  *sz_w += casadi_sp_nnz(p->bsp); // b
  *sz_w += casadi_sp_nnz(p->bsp); // b2
  *sz_w += casadi_sp_nnz(p->xsp);   // x
  *sz_w += casadi_sp_nnz(p->xsp);   // q
  *sz_w += casadi_sp_nnz(p->usp);   // u
  *sz_w += casadi_sp_nnz(p->usp);   // r
  *sz_w += casadi_sp_nnz(p->lugsp); // lg
  *sz_w += casadi_sp_nnz(p->lugsp); // ug
  *sz_w += casadi_sp_nnz(p->pisp); // pi
  *sz_w += p->nx_total; // lbx
  *sz_w += p->nx_total; // ubx
  *sz_w += p->nu_total; // lbu
  *sz_w += p->nu_total; // ubu
  *sz_w += p->nx_total+p->ng_total; // lam

  *sz_iw += p->nx_total; // iidxbx
  *sz_iw += p->nu_total; // iidxbu

}

// SYMBOL "hpipm_set_work"
template<typename T1>
void casadi_hpipm_set_work(casadi_hpipm_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int offset, i, k;
  
  const casadi_hpipm_prob<T1>* p = d->prob;

  d->hA = *res; *res += p->N;
  d->hB = *res; *res += p->N;
  d->hC = *res; *res += p->N + 1;
  d->hD = *res; *res += p->N + 1;
  d->hR = *res; *res += p->N + 1;
  d->hI = *res; *res += p->N;
  d->hS = *res; *res += p->N + 1;
  d->hQ = *res; *res += p->N + 1;
  d->hx = *res; *res += p->N + 1;
  d->hq = *res; *res += p->N + 1;
  d->hu = *res; *res += p->N + 1;
  d->hr = *res; *res += p->N + 1;
  d->hlg = *res; *res += p->N + 1;
  d->hug = *res; *res += p->N + 1;
  d->hb = *res; *res += p->N;
  d->pis  = *res; *res += p->N;
  d->hlbx  = *res; *res += p->N+1;
  d->hubx  = *res; *res += p->N+1;
  d->hlbu  = *res; *res += p->N+1;
  d->hubu  = *res; *res += p->N+1;
  d->lams  = *res; *res += p->N+1;
  d->hidxbx  = reinterpret_cast<int**>(*res); *res += p->N+1;
  d->hidxbu  = reinterpret_cast<int**>(*res); *res += p->N+1;
  d->hidxs = reinterpret_cast<int**>(*res); *res += p->N+1;

  d->hZl = *res; *res += p->N+1;
  d->hZu = *res; *res += p->N+1;
  d->hzl = *res; *res += p->N+1;
  d->hzu = *res; *res += p->N+1;
  d->hlls = *res; *res += p->N+1;
  d->hlus = *res; *res += p->N+1;

  d->A = *w; *w += casadi_sp_nnz(p->Asp);
  d->B = *w; *w += casadi_sp_nnz(p->Bsp);
  d->C = *w; *w += casadi_sp_nnz(p->Csp);
  d->D = *w; *w += casadi_sp_nnz(p->Dsp);
  d->R = *w; *w += casadi_sp_nnz(p->Rsp);
  d->I = *w; *w += casadi_sp_nnz(p->Isp);
  d->S = *w; *w += casadi_sp_nnz(p->Ssp);
  d->Q = *w; *w += casadi_sp_nnz(p->Qsp);
  d->b = *w; *w += casadi_sp_nnz(p->bsp);
  d->b2 = *w; *w += casadi_sp_nnz(p->bsp);
  d->x = *w; *w += casadi_sp_nnz(p->xsp);
  d->q = *w; *w += casadi_sp_nnz(p->xsp);
  d->u = *w; *w += casadi_sp_nnz(p->usp);
  d->r = *w; *w += casadi_sp_nnz(p->usp);
  d->lg = *w; *w += casadi_sp_nnz(p->lugsp);
  d->ug = *w; *w += casadi_sp_nnz(p->lugsp);
  d->pi = *w; *w += casadi_sp_nnz(p->pisp);

  d->lbx = *w; *w += p->nx_total;
  d->ubx = *w; *w += p->nx_total;
  d->lbu = *w; *w += p->nu_total;
  d->ubu = *w; *w += p->nu_total;
  d->lam = *w; *w += p->nx_total+p->ng_total;

  d->iidxbx = reinterpret_cast<int*>(*iw); *iw += p->nx_total;
  d->iidxbu = reinterpret_cast<int*>(*iw); *iw += p->nu_total;

  d->pv = *w;

  casadi_hpipm_ptr_block(p->N, d->hA, d->A, p->A, 0);
  casadi_hpipm_ptr_block(p->N, d->hB, d->B, p->B, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hC, d->C, p->C, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hD, d->D, p->D, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hR, d->R, p->R, 0);
  casadi_hpipm_ptr_block(p->N, d->hI, d->I, p->I, 1);
  casadi_hpipm_ptr_block(p->N+1, d->hS, d->S, p->S, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hQ, d->Q, p->Q, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hx, d->x, p->x, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hq, d->q, p->x, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hu, d->u, p->u, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hr, d->r, p->u, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hlg, d->lg, p->lug, 0);
  casadi_hpipm_ptr_block(p->N+1, d->hug, d->ug, p->lug, 0);
  casadi_hpipm_ptr_block(p->N, d->hb, d->b, p->b, 0);

  // Slack variables featues not currently interfaced
  for(i=0;i<p->N+1;++i) d->hZl[i] = 0;
  for(i=0;i<p->N+1;++i) d->hZu[i] = 0;
  for(i=0;i<p->N+1;++i) d->hzl[i] = 0;
  for(i=0;i<p->N+1;++i) d->hzu[i] = 0;
  for(i=0;i<p->N+1;++i) d->hlls[i] = 0;
  for(i=0;i<p->N+1;++i) d->hlus[i] = 0;

  offset = 0;
  for (k=0;k<p->N;++k) {
    d->pis[k] = d->pi+offset;
    offset += p->nx[k+1];
  }

  offset = 0;
  for (k=0;k<p->N+1;++k) {
    d->hlbx[k] = d->lbx+offset;
    d->hubx[k] = d->ubx+offset;
    offset += p->nx[k];
  }
  offset = 0;
  for (k=0;k<p->N+1;++k) {
    d->hlbu[k] = d->lbu+offset;
    d->hubu[k] = d->ubu+offset;
    offset += p->nu[k];
  }
  offset = 0;
  for (k=0;k<p->N+1;++k) {
    d->hidxbx[k] = d->iidxbx+offset;
    for (i=0;i<p->nbx[k];++i) d->hidxbx[k][i] = i;
    offset+= p->nbx[k];
  }
  offset = 0;
  for (k=0;k<p->N+1;++k) {
    d->hidxbu[k] = d->iidxbu+offset;
    for (i=0;i<p->nbu[k];++i) d->hidxbu[k][i] = i;
    offset+= p->nbu[k];
  }

}


// SYMBOL "hpipm_solve"
template<typename T1>
int casadi_hpipm_solve(casadi_hpipm_data<T1>* d, const double** arg, double** res, casadi_int* iw, double* w) {
    casadi_int k, i, j, n_row, offset;
    T1 f;
    const casadi_hpipm_prob<T1>* p = d->prob;
    const casadi_qp_prob<T1>* p_qp = p->qp;
    casadi_qp_data<T1>* d_qp = d->qp;

    int dim_size = d_ocp_qp_dim_memsize(p->N);
	  void *dim_mem = malloc(dim_size);
    
    struct d_ocp_qp_dim dim;
	  d_ocp_qp_dim_create(p->N, &dim, dim_mem);

    d_ocp_qp_dim_set_all(const_cast<int*>(p->nx), const_cast<int*>(p->nu), const_cast<int*>(p->nbx), const_cast<int*>(p->nbu),
      const_cast<int*>(p->ng), const_cast<int*>(p->nsbx), const_cast<int*>(p->nsbu), const_cast<int*>(p->nsg), &dim);

    int qp_size = d_ocp_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);

    struct d_ocp_qp qp;
    d_ocp_qp_create(&dim, &qp, qp_mem);

    // Dissect A matrix
    casadi_project(d_qp->a, p_qp->sp_a, d->A, p->Asp, d->pv);
    casadi_project(d_qp->a, p_qp->sp_a, d->B, p->Bsp, d->pv);
    casadi_project(d_qp->a, p_qp->sp_a, d->C, p->Csp, d->pv);
    casadi_project(d_qp->a, p_qp->sp_a, d->D, p->Dsp, d->pv);
    casadi_project(d_qp->a, p_qp->sp_a, d->I, p->Isp, d->pv);

    // Dissect H matrix; definition of HPIPM lacks a factor 2
    casadi_hpipm_mproject(0.5, d_qp->h, p_qp->sp_h, d->R, p->Rsp, d->pv);
    casadi_hpipm_mproject(0.5, d_qp->h, p_qp->sp_h, d->S, p->Ssp, d->pv);
    casadi_hpipm_mproject(0.5, d_qp->h, p_qp->sp_h, d->Q, p->Qsp, d->pv);

    // Dissect LBA/UBA
    casadi_hpipm_mproject(-1.0, d_qp->lba, p->sp_ba, d->b, p->bsp, d->pv);
    casadi_hpipm_mproject(-1.0, d_qp->uba, p->sp_ba, d->b2, p->bsp, d->pv);
    //casadi_assert_dev(std::equal(m->b.begin(), m->b.end(), m->b2.begin()));
    casadi_project(d_qp->lba, p->sp_ba, d->lg, p->lugsp, d->pv);
    casadi_project(d_qp->uba, p->sp_ba, d->ug, p->lugsp, d->pv);

    // Dissect LBX/UBX input
    casadi_fill(d->lbu, p->nu_total, 0.0);
    casadi_fill(d->ubu, p->nu_total, 0.0);
    casadi_fill(d->lbx, p->nx_total, 0.0);
    casadi_fill(d->ubx, p->nx_total, 0.0);

    casadi_hpipm_dense_transfer(1.0, d_qp->lbx, p->xsp, d->lbx, p->theirs_Xsp, d->pv);
    casadi_hpipm_dense_transfer(1.0, d_qp->ubx, p->xsp, d->ubx, p->theirs_Xsp, d->pv);
    casadi_hpipm_dense_transfer(1.0, d_qp->lbx, p->usp, d->lbu, p->theirs_Usp, d->pv);
    casadi_hpipm_dense_transfer(1.0, d_qp->ubx, p->usp, d->ubu, p->theirs_Usp, d->pv);

    // Dissect G
    casadi_hpipm_mproject(0.5, d_qp->g, p->sp_x, d->r, p->usp, d->pv);
    casadi_hpipm_mproject(0.5, d_qp->g, p->sp_x, d->q, p->xsp, d->pv);

    // Dissect X0
    casadi_project(d_qp->x0, p->sp_x, d->u, p->usp, d->pv);
    casadi_project(d_qp->x0, p->sp_x, d->x, p->xsp, d->pv);

    d->iter_count = -1;

    // Deal with non-unity I block
    for (k=0;k<p->N;++k) {
      n_row = p->nx[k+1];
      for (i=0;i<n_row;++i) {
        f = -1/d->hI[k][i];
        d->hb[k][i]*=f;
        for (j=0;j<p->nx[k];++j) d->hA[k][i+j*n_row]*=f;
        for (j=0;j<p->nu[k];++j) d->hB[k][i+j*n_row]*=f;
      }
    }

    // replace infinities
    casadi_clip_min(d->lbx, p->nx_total, -p->inf, 0);
    casadi_clip_min(d->lbu, p->nu_total, -p->inf, 0);
    casadi_clip_min(d->lg,  p->ng_total, -p->inf, 0); // count?
    casadi_clip_max(d->ubx, p->nx_total, p->inf, 0);
    casadi_clip_max(d->ubu, p->nu_total, p->inf, 0);
    casadi_clip_max(d->ug,  p->ng_total, p->inf, 0); // count?


    //m->fstats.at("preprocessing").toc();



    casadi_fill(d->pi, casadi_sp_nnz(p->pisp), 0.0);
    casadi_fill(d->lam, p->nx_total+p->ng_total, 0.0);

    d_ocp_qp_set_all(d->hA, d->hB, d->hb,
      d->hQ, d->hS, d->hR,
      d->hq, d->hr,
      d->hidxbx, d->hlbx, d->hubx,
      d->hidxbu, d->hlbu, d->hubu,
      d->hC, d->hD,
      d->hlg, d->hug,
      d->hZl, d->hZu, d->hzl,
      d->hzu, d->hidxs,
      d->hlls, d->hlus, &qp);

    int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);

    struct d_ocp_qp_sol qp_sol;
    d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
    void *ipm_arg_mem = malloc(ipm_arg_size);
    struct d_ocp_qp_ipm_arg myarg;
    d_ocp_qp_ipm_arg_create(&dim, &myarg, ipm_arg_mem);

    memcpy(&myarg, &p->hpipm_options, sizeof(struct d_ocp_qp_ipm_arg));

    int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &myarg);
    void *ipm_mem = malloc(ipm_size);

    struct d_ocp_qp_ipm_ws workspace;
    d_ocp_qp_ipm_ws_create(&dim, &myarg, &workspace, ipm_mem);

    // solution guess
    for (i=0; i<p->N+1; i++)	d_ocp_qp_sol_set_u(i, d->hu[i], &qp_sol);
    for (i=0; i<p->N+1; i++)	d_ocp_qp_sol_set_x(i, d->hx[i], &qp_sol);


    d_ocp_qp_print(&dim, &qp);
    //m->fstats.at("solver").tic();
		// call solver
	  d_ocp_qp_ipm_solve(&qp, &qp_sol, &myarg, &workspace);
    d_ocp_qp_ipm_get_status(&workspace, &d->return_status);
    //m->fstats.at("solver").toc();
    //m->fstats.at("postprocessing").tic();

    d_ocp_qp_ipm_get_iter(&workspace, &d->iter_count);
    d_ocp_qp_ipm_get_max_res_stat(&workspace, &d->res_stat);
	  d_ocp_qp_ipm_get_max_res_eq(&workspace, &d->res_eq);
	  d_ocp_qp_ipm_get_max_res_ineq(&workspace, &d->res_ineq);
	  d_ocp_qp_ipm_get_max_res_comp(&workspace, &d->res_comp);

	  double *stat;
    d_ocp_qp_ipm_get_stat(&workspace, &stat);
    printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
    d_print_exp_tran_mat(5, d->iter_count, stat, 5);

    printf("\nHPIPM returned with flag %i.\n", d->return_status);
    if(d->return_status == 0)
		{
        printf("\n -> QP solved!\n");
		}
	else if(d->return_status==1)
		{
        printf("\n -> Solver failed! Maximum number of iterations reached\n");
		}
	else if(d->return_status==2)
		{
        printf("\n -> Solver failed! Minimum step lenght reached\n");
		}
	else if(d->return_status==2)
		{
        printf("\n -> Solver failed! NaN in computations\n");
		}
	else
		{
        printf("\n -> Solver failed! Unknown return flag\n");
		}
	printf("\n\n");


	for(i=0; i<p->N+1; ++i) d_ocp_qp_sol_get_u(i, &qp_sol, d->hu[i]);
	for(i=0; i<p->N+1; ++i) d_ocp_qp_sol_get_x(i, &qp_sol, d->hx[i]);
	for(i=0; i<p->N; ++i) d_ocp_qp_sol_get_pi(i, &qp_sol, d->pis[i]);


  if (d_qp->lam_x) {
    offset = 0;
    for (i=0; i<p->N+1; ++i) {
      double* lam_lb = d->pv;
      double* lam_ub = d->pv+p->nbx[i]+p->nbu[i];
      d_ocp_qp_sol_get_lam_lb(i, &qp_sol, lam_lb);
      d_ocp_qp_sol_get_lam_ub(i, &qp_sol, lam_ub);
      for (k=0;k<p->nbx[i];++k) {
        d_qp->lam_x[offset+k] = (lam_ub[k+p->nbu[i]]-lam_lb[k+p->nbu[i]])*2;
      }
      offset += p->nbx[i];
      for (k=0;k<p->nbu[i];++k) {
        d_qp->lam_x[offset+k] = (lam_ub[k]-lam_lb[k])*2;
      }
      offset += p->nbu[i];
    }
  }
  
  if (d_qp->lam_a) {
    offset = 0;
    for (i=0; i<p->N+1; ++i) {
      double* lam_lg = d->pv;
      double* lam_ug = d->pv+p->ng[i];
      d_ocp_qp_sol_get_lam_lg(i, &qp_sol, lam_lg);
      d_ocp_qp_sol_get_lam_ug(i, &qp_sol, lam_ug);
      for (k=0;k<p->ng[i];++k) {
        d_qp->lam_a[offset+(i<p->N ? p->nx[i+1]: 0)+k] = (lam_ug[k]-lam_lg[k])*2;
      }
      offset += p->ng[i]+(i<p->N ? p->nx[i+1]: 0);
    }
  }

    casadi_hpipm_dense_transfer(1.0, d->x, p->theirs_Xsp, d_qp->x, p->xsp, d->pv);
    casadi_hpipm_dense_transfer(1.0, d->u, p->theirs_Usp, d_qp->x, p->usp, d->pv);


    // Deal with non-unity I block
    for (k=0;k<p->N;++k) {
      n_row = p->nx[k+1];
      for (i=0;i<n_row;++i) {
        f = -1/d->hI[k][i];
        d->pis[k][i]*=f;
      }
    }

    casadi_hpipm_dense_transfer(2.0, d->pi, p->pisp, d_qp->lam_a, p->lamg_gapsp, d->pv);

    // Construct f
    f = casadi_dot(p_qp->nx, d_qp->g, d_qp->x);
    f += 0.5*casadi_bilin(d_qp->h, p_qp->sp_h, d_qp->x, d_qp->x);

    if (d_qp->f) d_qp->f[0] = f;

    /*
    m->fstats.at("postprocessing").toc();

    */

   return d->return_status==0;
}