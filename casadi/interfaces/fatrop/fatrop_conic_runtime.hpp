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


// SYMBOL "fatrop_mproject"
template<typename T1>
void casadi_fatrop_conic_mproject(T1 factor, const T1* x, const casadi_int* sp_x,
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
void casadi_fatrop_conic_dense_transfer(double factor, const T1* x,
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

// SYMBOL "fatrop_block"
struct casadi_ocp_block {
    casadi_int offset_r;
    casadi_int offset_c;
    casadi_int rows;
    casadi_int cols;
};
// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"

// SYMBOL "fatrop_unpack_blocks"
inline void casadi_fatrop_conic_unpack_blocks(casadi_int N, casadi_ocp_block* blocks, const casadi_int* packed) {
    casadi_int i;
    for (i=0;i<N;++i) {
        blocks[i].offset_r = *packed++;
        blocks[i].offset_c = *packed++;
        blocks[i].rows = *packed++;
        blocks[i].cols = *packed++;
    }
}

// SYMBOL "fatrop_ptr_block"
template<typename T1>
void casadi_fatrop_conic_ptr_block(casadi_int N, T1** vs, T1* v, const casadi_ocp_block* blocks, int eye) {
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
struct casadi_fatrop_conic_prob {
  const casadi_qp_prob<T1>* qp;
  const int *nx, *nu, *ng;

  // Sparsities
  const casadi_int *ABsp;
  const casadi_int *AB_offsets;
  const casadi_int *CDsp;
  const casadi_int *CD_offsets;
  const casadi_int *RSQsp;
  const casadi_int *RSQ_offsets;

  casadi_int N;
  casadi_ocp_block *AB, *CD, *RSQ;

  T1 warm_start;
  T1 inf;

};
// C-REPLACE "casadi_fatrop_conic_prob<T1>" "struct casadi_fatrop_conic_prob"

// SYMBOL "fatrop_data"
template<typename T1>
struct casadi_fatrop_conic_data {
  // Problem structure
  const casadi_fatrop_conic_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;

  T1 *AB, *CD, *RSQ;

  casadi_int *a_eq, *a_ineq, *a_eq_idx, *a_ineq_idx;
  casadi_int *x_eq, *x_ineq, *x_eq_idx, *x_ineq_idx;

  int iter_count;
  const char* return_status;
  T1 res_stat;
  T1 res_eq;
  T1 res_ineq;
  T1 res_comp;

  T1 *pv;
};
// C-REPLACE "casadi_fatrop_conic_data<T1>" "struct casadi_fatrop_conic_data"


// SYMBOL "fatrop_setup"
template<typename T1>
void casadi_fatrop_conic_setup(casadi_fatrop_conic_prob<T1>* p) {


}

// SYMBOL "fatrop_solve"
template<typename T1>
int casadi_fatrop_conic_solve(casadi_fatrop_conic_data<T1>* d, const double** arg, double** res, casadi_int* iw, double* w) {
    casadi_int k, i, j, n_row, offset, start, stop;
    T1 f;
    const casadi_fatrop_conic_prob<T1>* p = d->prob;
    const casadi_qp_prob<T1>* p_qp = p->qp;
    casadi_qp_data<T1>* d_qp = d->qp;

    casadi_fatrop_conic_mproject(-1.0, d_qp->a, p_qp->sp_a, d->AB, p->ABsp, d->pv);
    casadi_project(d_qp->a, p_qp->sp_a, d->CD, p->CDsp, d->pv);

    casadi_project(d_qp->h, p_qp->sp_h, d->RSQ, p->RSQsp, d->pv);

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
        if (d_qp->lba[i]==d_qp->uba[i]) {
          d->a_eq[d->a_eq_idx[k+1]++] = i;
        } else {
          if (d_qp->lba[i]==-inf && d_qp->uba[i]==inf) continue;
          d->a_ineq[d->a_ineq_idx[k+1]++] = i;
        }
      }
      d->x_eq_idx[k+1] = d->x_eq_idx[k];
      d->x_ineq_idx[k+1] = d->x_ineq_idx[k];
      start = p->CD[k].offset_c;
      stop  = p->CD[k].offset_c+p->CD[k].cols;
      //uout() << "start" << start << "->" << stop << std::endl;
      for (i=start;i<stop;++i) {
        if (d_qp->lbx[i]==d_qp->ubx[i]) {
          d->x_eq[d->x_eq_idx[k+1]++] = i;
        } else {
          if (d_qp->lbx[i]==-inf && d_qp->ubx[i]==inf) continue;
          d->x_ineq[d->x_ineq_idx[k+1]++] = i;
        }
      }
      //uout() << "k=" << k << std::endl;
      //uout() << "a_eq" << std::vector<double>(d->a_eq+d->a_eq_idx[k], d->a_eq+d->a_eq_idx[k+1]) << std::endl;
      //uout() << "a_ineq" << std::vector<double>(d->a_ineq+d->a_ineq_idx[k], d->a_ineq+d->a_ineq_idx[k+1]) << std::endl;
      //uout() << "x_eq" << std::vector<double>(d->x_eq+d->x_eq_idx[k], d->x_eq+d->x_eq_idx[k+1]) << std::endl;
      //uout() << "x_ineq" << std::vector<double>(d->x_ineq+d->x_ineq_idx[k], d->x_ineq+d->x_ineq_idx[k+1]) << std::endl;

      //uout() << "AB=" << std::vector<double>(d->AB,d->AB+100) << std::endl;

      //uout() << "CD=" << std::vector<double>(d->CD,d->CD+100) << std::endl;

      //uout() << "RSQ=" << std::vector<double>(d->RSQ,d->RSQ+100) << std::endl;
    }


  return 0;
}

// SYMBOL "fatrop_work"
template<typename T1>
void casadi_fatrop_conic_work(const casadi_fatrop_conic_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->qp->nx+p->qp->na)); // pv
  // Persistent work vectors
  *sz_w += casadi_sp_nnz(p->ABsp); // AB
  *sz_w += casadi_sp_nnz(p->CDsp); // CD
  *sz_w += casadi_sp_nnz(p->RSQsp); // RSQ

  *sz_iw += p->N+2; // a_eq_idx
  *sz_iw += p->N+2; // a_ineq_idx
  *sz_iw += p->N+2; // x_eq_idx
  *sz_iw += p->N+2; // x_ineq_idx
  *sz_iw += p->qp->na; // a_eq
  *sz_iw += p->qp->na; // a_ineq
  *sz_iw += p->qp->nx; // x_eq
  *sz_iw += p->qp->nx; // x_ineq

}


// SYMBOL "fatrop_set_work"
template<typename T1>
void casadi_fatrop_conic_set_work(casadi_fatrop_conic_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int offset, i, k;
  
  const casadi_fatrop_conic_prob<T1>* p = d->prob;

  d->AB = *w; *w += casadi_sp_nnz(p->ABsp);
  d->CD = *w; *w += casadi_sp_nnz(p->CDsp);
  d->RSQ = *w; *w += casadi_sp_nnz(p->RSQsp);

  d->a_eq_idx = *iw;   *iw += p->N+2;
  d->a_ineq_idx = *iw; *iw += p->N+2;
  d->x_eq_idx = *iw;   *iw += p->N+2;
  d->x_ineq_idx = *iw; *iw += p->N+2;
  
  d->a_eq = *iw;   *iw += p->qp->na;
  d->a_ineq = *iw; *iw += p->qp->na;
  d->x_eq = *iw;  *iw += p->qp->nx;
  d->x_ineq = *iw; *iw += p->qp->nx;


  d->pv = *w;
}