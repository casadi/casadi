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
// C-REPLACE "const_cast<DAQPSettings*>" "(DAQPSettings*) "
// C-REPLACE "const_cast<T1*>" "(casadi_real*) "

template<typename T1>
struct casadi_daqp_prob {
  const casadi_qp_prob<T1>* qp;

  DAQPSettings settings;
};
// C-REPLACE "casadi_daqp_prob<T1>" "struct casadi_daqp_prob"

// SYMBOL "daqp_setup"
template<typename T1>
void casadi_daqp_setup(casadi_daqp_prob<T1>* p) {

}



// SYMBOL "daqp_data"
template<typename T1>
struct casadi_daqp_data {
  // Problem structure
  const casadi_daqp_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;

  DAQPWorkspace work;
  DAQPProblem daqp;
  DAQPResult res;

  int return_status;

};
// C-REPLACE "casadi_daqp_data<T1>" "struct casadi_daqp_data"

// SYMBOL "daqp_init_mem"
template<typename T1>
int daqp_init_mem(casadi_daqp_data<T1>* d) {
  return 0;
}

// SYMBOL "daqp_free_mem"
template<typename T1>
void daqp_free_mem(casadi_daqp_data<T1>* d) {

}

// SYMBOL "daqp_work"
template<typename T1>
void casadi_daqp_work(const casadi_daqp_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
  // Local variables
  casadi_int nx, na;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  // Get matrix number of nonzeros
  na = p_qp->na;
  nx = p_qp->nx;

  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  // Temporary work vectors
  // Persistent work vectors
  *sz_w += nx*nx; // H
  *sz_w += na*nx; // A
  *sz_w += p->qp->nz; // lbz
  *sz_w += p->qp->nz; // ubz
  *sz_w += p->qp->nz; // lam
  *sz_iw += p->qp->nz; // sense
}

// SYMBOL "daqp_init"
template<typename T1>
void casadi_daqp_init(casadi_daqp_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nx, na;
  const casadi_daqp_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  // Get matrix number of nonzeros
  na = p_qp->na;
  nx = p_qp->nx;

  d->daqp.H = *w; *w += nx*nx;
  d->daqp.A = *w; *w += na*nx;
  d->daqp.blower = *w; *w += p->qp->nz;
  d->daqp.bupper = *w; *w += p->qp->nz;
  d->res.lam = *w; *w += p->qp->nz;
  d->daqp.sense = (int*) *iw; *iw += p->qp->nz;
}


// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_LIMITED" "2"

// SYMBOL "daqp_solve"
template<typename T1>
int casadi_daqp_solve(casadi_daqp_data<T1>* d, const double** arg, double** res, casadi_int* iw, double* w) {
  // Local variables
  casadi_int i;
  int flag;
  const casadi_daqp_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;

  for (i=0;i<p_qp->nx;++i) {
    d->daqp.sense[i] = d_qp->lbx[i]==d_qp->ubx[i] ? 5 : 0;
  }
  for (i=0;i<p_qp->na;++i) {
    d->daqp.sense[p_qp->nx+i] = d_qp->lba[i]==d_qp->uba[i] ? 5 : 0;
  }

  d->daqp.n = p_qp->nx;
  d->daqp.m = p_qp->nx + p_qp->na;
  d->daqp.ms = p_qp->nx;
  d->work.settings = const_cast<DAQPSettings*>(&p->settings);
  d->res.x = d_qp->x;
  casadi_densify(d->qp->h, d->prob->qp->sp_h, d->daqp.H, 0);
  casadi_densify(d->qp->a, d->prob->qp->sp_a, d->daqp.A, 1);
  casadi_copy(d_qp->lbx, p_qp->nx, d->daqp.blower);
  casadi_copy(d_qp->lba, p_qp->na, d->daqp.blower+p_qp->nx);
  casadi_copy(d_qp->ubx, p_qp->nx, d->daqp.bupper);
  casadi_copy(d_qp->uba, p_qp->na, d->daqp.bupper+p_qp->nx);
  d->daqp.f = const_cast<T1*>(d_qp->g);

  flag = setup_daqp(&d->daqp,&d->work,&(d->res.setup_time));
  if (flag<0) return 1;
  daqp_solve(&d->res,&d->work);
  casadi_copy(d->res.lam, p_qp->nx, d_qp->lam_x);
  casadi_copy(d->res.lam+p_qp->nx, p_qp->na, d_qp->lam_a);
  *d_qp->f = d->res.fval;
  d->work.settings = 0;
  free_daqp_workspace(&d->work);
  free_daqp_ldp(&d->work);

  d->return_status = d->res.exitflag;
  d_qp->success = d->res.exitflag==EXIT_OPTIMAL;

/**
  #define EXIT_SOFT_OPTIMAL 2 
#define EXIT_OPTIMAL 1
#define EXIT_INFEASIBLE -1
#define EXIT_CYCLE -2
#define EXIT_UNBOUNDED -3
#define EXIT_ITERLIMIT -4
#define EXIT_NONCONVEX -5
#define EXIT_OVERDETERMINED_INITIAL -6
*/

  return 0;
}