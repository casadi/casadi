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


// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "std::numeric_limits<T1>::min()" "casadi_real_min"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "static_cast<int>" "(int) "

// SYMBOL "qp_prob"
template<typename T1>
struct casadi_qp_prob {
  // Sparsity patterns
  const casadi_int *sp_a, *sp_h;
  // Dimensions
  casadi_int nx, na, nz;
  casadi_int nnz_a, nnz_h;
};
// C-REPLACE "casadi_qp_prob<T1>" "struct casadi_qp_prob"

// SYMBOL "qp_setup"
template<typename T1>
void casadi_qp_setup(casadi_qp_prob<T1>* p) {
  p->na = p->sp_a[0];
  p->nx = p->sp_a[1];
  p->nz = p->na+p->nx;
  p->nnz_a = p->sp_a[2+p->sp_a[1]];
  p->nnz_h = p->sp_h[2+p->sp_h[1]];
}

// C-REPLACE "UnifiedReturnStatus" "int"
// C-REPLACE "bool" "int"

// SYMBOL "qp_data"
template<typename T1>
struct casadi_qp_data {
  // Problem structure
  const casadi_qp_prob<T1>* prob;

  UnifiedReturnStatus unified_return_status;
  bool success;

  // Number of iterations performed
  casadi_int iter_count;

  // QP data, pointers to arg (no allocations needed)
  const T1 *a, *h, *g, *lbx, *ubx, *lba, *uba, *x0, *lam_x0, *lam_a0;
  // QP results, pointers to res (no allocations needed)
  T1 *f, *x, *lam_x, *lam_a;

};
// C-REPLACE "casadi_qp_data<T1>" "struct casadi_qp_data"

// SYMBOL "qp_init"
template<typename T1>
void casadi_qp_init(casadi_qp_data<T1>* d, casadi_int** iw, T1** w) {
  // Local variables
  //const casadi_qp_prob<T1>* p = d->prob;
}

// SYMBOL "qp_work"
template<typename T1>
void casadi_qp_work(const casadi_qp_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_arg, sz_res
  *sz_arg = *sz_res = 0;
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
}
