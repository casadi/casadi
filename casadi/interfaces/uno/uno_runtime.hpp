//
//    MIT No Attribution
//
//    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//

// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"
// C-REPLACE "OracleCallback" "struct casadi_oracle_callback"
// C-REPLACE "calc_function" "casadi_oracle_call"
// C-REPLACE "static_cast< struct casadi_uno_data* >" "(struct casadi_uno_data*)"
// C-REPLACE "static_cast<uno_int>" "(uno_int) "
// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_UNKNOWN" "1"
// C-REPLACE "SOLVER_RET_LIMITED" "2"
// C-REPLACE "SOLVER_RET_NAN" "3"
// C-REPLACE "SOLVER_RET_INFEASIBLE" "4"
// C-REPLACE "SOLVER_RET_EXCEPTION" "5"

// SYMBOL "uno_prob"
template<typename T1>
struct casadi_uno_prob {
  casadi_nlpsol_prob<T1> nlp;
  const casadi_int* sp_a;
  const casadi_int* sp_h;
  const uno_int* jac_row;
  const uno_int* jac_col;
  uno_int n_jac;
  const uno_int* hess_row;
  const uno_int* hess_col;
  uno_int n_hess;
  OracleCallback nlp_f;
  OracleCallback nlp_g;
  OracleCallback nlp_grad_f;
  OracleCallback nlp_jac_g;
  OracleCallback nlp_hess_l;
  OracleCallback fwd1_nlp_grad_l;
  uno_objective_callback                    obj_cb;
  uno_objective_gradient_callback           obj_grad_cb;
  uno_constraints_callback                  constr_cb;
  uno_constraints_jacobian_callback         jac_cb;
  uno_lagrangian_hessian_callback           hess_cb;
  uno_lagrangian_hessian_operator_callback  hess_prod_cb;
};
// C-REPLACE "casadi_uno_prob<T1>" "struct casadi_uno_prob"

// SYMBOL "uno_data"
template<typename T1>
struct casadi_uno_data {
  const casadi_uno_prob<T1>* prob;
  void* solver;
  void* model;
  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;
  int unified_return_status;
  int success;
  T1 primal_infeasibility;
  T1 stationarity;
  T1 complementarity;
  uno_int iter_count;
  casadi_nlpsol_data<T1> nlp;
  casadi_oracle_data<T1> d_oracle;
};
// C-REPLACE "casadi_uno_data<T1>" "struct casadi_uno_data"

// SYMBOL "uno_obj_wrapper"
template<typename T1>
uno_int casadi_uno_obj_wrapper(uno_int n, const T1* x, T1* fval, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;
  d_oracle->arg[1] = d->nlp.p;
  d_oracle->res[0] = fval;
  return calc_function(&d->prob->nlp_f, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_obj_wrapper<T1>" "casadi_uno_obj_wrapper"

// SYMBOL "uno_obj_grad_wrapper"
template<typename T1>
uno_int casadi_uno_obj_grad_wrapper(uno_int n, const T1* x, T1* grad, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;
  d_oracle->arg[1] = d->nlp.p;
  d_oracle->res[0] = grad;
  return calc_function(&d->prob->nlp_grad_f, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_obj_grad_wrapper<T1>" "casadi_uno_obj_grad_wrapper"

// SYMBOL "uno_constr_wrapper"
template<typename T1>
uno_int casadi_uno_constr_wrapper(uno_int n, uno_int ng, const T1* x, T1* gval, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;
  d_oracle->arg[1] = d->nlp.p;
  d_oracle->res[0] = gval;
  return calc_function(&d->prob->nlp_g, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_constr_wrapper<T1>" "casadi_uno_constr_wrapper"

// SYMBOL "uno_jac_wrapper"
template<typename T1>
uno_int casadi_uno_jac_wrapper(uno_int n, uno_int nnz, const T1* x, T1* jvals, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;
  d_oracle->arg[1] = d->nlp.p;
  d_oracle->res[0] = jvals;
  return calc_function(&d->prob->nlp_jac_g, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_jac_wrapper<T1>" "casadi_uno_jac_wrapper"

// SYMBOL "uno_hess_wrapper"
template<typename T1>
uno_int casadi_uno_hess_wrapper(uno_int n, uno_int ng, uno_int nnz,
    const T1* x, T1 obj_mult, const T1* mults, T1* hvals, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;
  d_oracle->arg[1] = d->nlp.p;
  d_oracle->arg[2] = &obj_mult;
  d_oracle->arg[3] = mults;
  d_oracle->res[0] = hvals;
  return calc_function(&d->prob->nlp_hess_l, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_hess_wrapper<T1>" "casadi_uno_hess_wrapper"

// SYMBOL "uno_hess_prod_wrapper"
template<typename T1>
uno_int casadi_uno_hess_prod_wrapper(uno_int n, uno_int ng, const T1* x,
    bool evaluate_at_x, T1 obj_mult, const T1* mults, const T1* vec,
    T1* result, void* user_data) {
  casadi_uno_data<T1>* d = static_cast< casadi_uno_data<T1>* >(user_data);
  casadi_oracle_data<T1>* d_oracle = d->nlp.oracle;
  d_oracle->arg[0] = x;          // x
  d_oracle->arg[1] = d->nlp.p;   // p
  d_oracle->arg[2] = &obj_mult;  // lam:f
  d_oracle->arg[3] = mults;      // lam:g
  d_oracle->arg[4] = 0;          // out:grad:gamma:x
  d_oracle->arg[5] = vec;        // fwd:x
  d_oracle->arg[6] = 0;          // fwd:p
  d_oracle->arg[7] = 0;          // fwd:lam:f
  d_oracle->arg[8] = 0;          // fwd:lam:g
  d_oracle->res[0] = result;     // fwd:grad:gamma:x
  return calc_function(&d->prob->fwd1_nlp_grad_l, d_oracle) == 0 ? 0 : 1;
}
// C-REPLACE "casadi_uno_hess_prod_wrapper<T1>" "casadi_uno_hess_prod_wrapper"

// SYMBOL "uno_term_cb"
template<typename T1>
uno_int casadi_uno_term_cb(uno_int n, uno_int ng, const T1* primals,
    const T1* lower_mult, const T1* upper_mult, const T1* constr_mult,
    T1 obj_mult, T1 inf_pr, T1 inf_du, T1 compl_res, void* user_data) {
  return 1;
}
// C-REPLACE "casadi_uno_term_cb<T1>" "casadi_uno_term_cb"

// SYMBOL "uno_init_mem"
template<typename T1>
int casadi_uno_init_mem(casadi_uno_data<T1>* d) {
  d->solver = uno_create_solver();
  uno_set_solver_callbacks(d->solver, 0, &casadi_uno_term_cb<T1>, d);
  d->model = 0;
  d->unified_return_status = SOLVER_RET_UNKNOWN;
  d->success = 0;
  return 0;
}

// SYMBOL "uno_init_model"
template<typename T1>
void casadi_uno_init_model(casadi_uno_data<T1>* d,
    const T1* lb_g, const T1* ub_g) {
  const casadi_uno_prob<T1>* p = d->prob;
  uno_int nx = static_cast<uno_int>(p->nlp.nx);
  uno_int ng = static_cast<uno_int>(p->nlp.ng);
  // Variable bounds are set per-solve via uno_set_variables_*_bounds.
  d->model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, nx, UNO_ZERO_BASED_INDEXING);
  uno_set_user_data(d->model, d);
  uno_set_objective(d->model, UNO_MINIMIZE, p->obj_cb, p->obj_grad_cb);
  if (ng > 0) {
    uno_set_constraints(d->model, ng, p->constr_cb, lb_g, ub_g,
        p->n_jac, p->jac_row, p->jac_col, p->jac_cb);
  }
  uno_set_lagrangian_hessian(d->model, p->n_hess, UNO_UPPER_TRIANGLE,
      p->hess_row, p->hess_col, p->hess_cb);
  uno_set_lagrangian_sign_convention(d->model, UNO_MULTIPLIER_POSITIVE);
  uno_set_lagrangian_hessian_operator(d->model, p->hess_prod_cb);
}

// SYMBOL "uno_free_mem"
template<typename T1>
void casadi_uno_free_mem(casadi_uno_data<T1>* d) {
  if (d->model)  uno_destroy_model(d->model);
  if (d->solver) uno_destroy_solver(d->solver);
  d->model = 0;
  d->solver = 0;
}

// SYMBOL "uno_set_work"
template<typename T1>
void casadi_uno_set_work(casadi_uno_data<T1>* d, const T1*** arg, T1*** res,
                     casadi_int** iw, T1** w) {
  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// SYMBOL "uno_solve"
template<typename T1>
void casadi_uno_solve(casadi_uno_data<T1>* d) {
  const casadi_uno_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = &d->nlp;
  uno_int nx = static_cast<uno_int>(p->nlp.nx);
  uno_int ng = static_cast<uno_int>(p->nlp.ng);
  casadi_int i;

  if (!d->model) {
    // Codegen path: model wasn't built in init_mem (no p_nlp scope there);
    // build it now using the real constraint bounds from this first call.
    casadi_uno_init_model(d, d_nlp->lbz + nx, d_nlp->ubz + nx);
  }
  // Variable bounds aren't part of the (unconstrained) model -- set every solve.
  uno_set_variables_lower_bounds(d->model, d_nlp->lbz);
  uno_set_variables_upper_bounds(d->model, d_nlp->ubz);
  if (ng > 0) {
    uno_set_constraints_lower_bounds(d->model, d_nlp->lbz + nx);
    uno_set_constraints_upper_bounds(d->model, d_nlp->ubz + nx);
  }
  uno_set_initial_primal_iterate(d->model, d_nlp->x0);

  uno_optimize(d->solver, d->model);

  uno_get_primal_solution(d->solver, d_nlp->z);
  uno_get_constraint_dual_solution(d->solver, d_nlp->lam + nx);
  for (i = 0; i < nx; ++i) {
    // Negative sign convention: bound multipliers are added, not subtracted.
    d_nlp->lam[i] = uno_get_lower_bound_dual_solution_component(d->solver, i)
                  + uno_get_upper_bound_dual_solution_component(d->solver, i);
  }
  d_nlp->objective         = uno_get_solution_objective(d->solver);
  d->primal_infeasibility  = uno_get_solution_primal_feasibility(d->solver);
  d->stationarity          = uno_get_solution_stationarity(d->solver);
  d->complementarity       = uno_get_solution_complementarity(d->solver);
  d->iter_count            = uno_get_number_iterations(d->solver);

  uno_int opt = uno_get_optimization_status(d->solver);
  uno_int sol = uno_get_solution_status(d->solver);
  d->success = (opt == UNO_SUCCESS && sol == UNO_FEASIBLE_KKT_POINT) ? 1 : 0;
  if      (opt == UNO_EVALUATION_ERROR)            d->unified_return_status = SOLVER_RET_NAN;
  else if (opt == UNO_ITERATION_LIMIT)             d->unified_return_status = SOLVER_RET_LIMITED;
  else if (opt == UNO_TIME_LIMIT)                  d->unified_return_status = SOLVER_RET_LIMITED;
  else if (opt == UNO_ALGORITHMIC_ERROR)           d->unified_return_status = SOLVER_RET_EXCEPTION;
  else if (opt == UNO_USER_TERMINATION)            d->unified_return_status = SOLVER_RET_UNKNOWN;
  else if (sol == UNO_FEASIBLE_KKT_POINT)          d->unified_return_status = SOLVER_RET_SUCCESS;
  else if (sol == UNO_FEASIBLE_FJ_POINT)           d->unified_return_status = SOLVER_RET_SUCCESS;
  else if (sol == UNO_INFEASIBLE_STATIONARY_POINT) d->unified_return_status = SOLVER_RET_INFEASIBLE;
  else if (sol == UNO_INFEASIBLE_SMALL_STEP)       d->unified_return_status = SOLVER_RET_INFEASIBLE;
  else if (sol == UNO_FEASIBLE_SMALL_STEP)         d->unified_return_status = SOLVER_RET_LIMITED;
  else                                             d->unified_return_status = SOLVER_RET_UNKNOWN;
}
