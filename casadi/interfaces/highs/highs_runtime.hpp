//
//    MIT No Attribution
//
//    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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

template<typename T1>
struct casadi_highs_prob {
  const casadi_qp_prob<T1>* qp;

  const int *colinda, *rowa;
  const int *colindh, *rowh;
  const int *integrality;
};
// C-REPLACE "casadi_highs_prob<T1>" "struct casadi_highs_prob"

// SYMBOL "highs_setup"
template<typename T1>
void casadi_highs_setup(casadi_highs_prob<T1>* p) {

}



// SYMBOL "highs_data"
template<typename T1>
struct casadi_highs_data {
  // Problem structure
  const casadi_highs_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;

  int return_status;

  int simplex_iteration_count;
  int ipm_iteration_count;
  int qp_iteration_count;
  int crossover_iteration_count;
  int primal_solution_status;
  int dual_solution_status;
  int basis_validity;
  T1 mip_dual_bound;
  T1 mip_gap;
  int num_primal_infeasibilities;
  T1 max_primal_infeasibility;
  T1 sum_primal_infeasibilities;
  int num_dual_infeasibilities;
  T1 max_dual_infeasibility;
  T1 sum_dual_infeasibilities;

  void* highs;
};
// C-REPLACE "casadi_highs_data<T1>" "struct casadi_highs_data"

// SYMBOL "highs_init_mem"
template<typename T1>
int highs_init_mem(casadi_highs_data<T1>* d) {
  d->highs = Highs_create();
  return 0;
}

// SYMBOL "highs_free_mem"
template<typename T1>
void highs_free_mem(casadi_highs_data<T1>* d) {
  Highs_destroy(d->highs);
}

// SYMBOL "highs_work"
template<typename T1>
void casadi_highs_work(const casadi_highs_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
}

// SYMBOL "highs_init"
template<typename T1>
void casadi_highs_init(casadi_highs_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  

}


// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_LIMITED" "2"

// SYMBOL "highs_solve"
template<typename T1>
int casadi_highs_solve(casadi_highs_data<T1>* d, const double** arg, double** res, casadi_int* iw, double* w) {

  const casadi_highs_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;


  // Create HiGHS instance and pass problem
  int status;
  const int matrix_format = 1;
  const int sense = 1;
  const double offset = 0.0;

  status = Highs_passModel(d->highs, p_qp->nx, p_qp->na, p_qp->nnz_a, p_qp->nnz_h,
    matrix_format, matrix_format, sense, offset,
    d_qp->g, d_qp->lbx, d_qp->ubx, d_qp->lba, d_qp->uba,
    p->colinda, p->rowa, d_qp->a,
    p->colindh, p->rowh, d_qp->h,
    p->integrality);
  
  if (!(status==kHighsStatusOk || status==kHighsStatusWarning)) return 1;

  // solve incumbent model
  status = Highs_run(d->highs);

  if (!(status==kHighsStatusOk || status==kHighsStatusWarning)) return 1;

  // get primal and dual solution
  Highs_getSolution(d->highs, d_qp->x, d_qp->lam_x, 0, d_qp->lam_a);
  
  if (d_qp->lam_x) {
    casadi_scal(p_qp->nx, -1., d_qp->lam_x);
  }
  if (d_qp->lam_a) {
    casadi_scal(p_qp->na, -1., d_qp->lam_a);
  }

  if (d_qp->f) {
    *d_qp->f = Highs_getObjectiveValue(d->highs);
  }

  status = Highs_getModelStatus(d->highs);
  d->return_status = status;
  d_qp->success = status==kHighsModelStatusOptimal;

  if (status==kHighsModelStatusOptimal)
  d_qp->unified_return_status = SOLVER_RET_SUCCESS;

  if (status==kHighsModelStatusTimeLimit
        || status==kHighsModelStatusIterationLimit)
    d_qp->unified_return_status = SOLVER_RET_LIMITED;

  if (Highs_getIntInfoValue(d->highs, "simplex_iteration_count", &d->simplex_iteration_count)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "ipm_iteration_count", &d->ipm_iteration_count)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "qp_iteration_count", &d->qp_iteration_count)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "crossover_iteration_count", &d->crossover_iteration_count)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "primal_solution_status", &d->primal_solution_status)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "dual_solution_status", &d->dual_solution_status)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "basis_validity", &d->basis_validity)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "mip_dual_bound", &d->mip_dual_bound)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "mip_gap", &d->mip_gap)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "num_primal_infeasibilities", &d->num_primal_infeasibilities)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "max_primal_infeasibility", &d->max_primal_infeasibility)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "sum_primal_infeasibilities", &d->sum_primal_infeasibilities)!=kHighsStatusOk) return 1;
  if (Highs_getIntInfoValue(d->highs, "num_dual_infeasibilities", &d->num_dual_infeasibilities)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "max_dual_infeasibility", &d->max_dual_infeasibility)!=kHighsStatusOk) return 1;
  if (Highs_getDoubleInfoValue(d->highs, "sum_dual_infeasibilities", &d->sum_dual_infeasibilities)!=kHighsStatusOk) return 1;


  return 0;
}