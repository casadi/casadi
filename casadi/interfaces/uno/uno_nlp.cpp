#include "uno_nlp.hpp"
#include "uno_interface.hpp"

namespace casadi {

/* ------- Define functions to evaluate solver functions  -----------*/
UnoNlp::UnoNlp(UnoMemory* mem) : mem_(mem)
{}

uno_int UnoNlp::objective_function(const double* x, double* objective_value)
{
  mem_->arg[0] = x;
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->res[0] = objective_value;
  casadi_assert(mem_->self.calc_function(mem_, "nlp_f")==0, "Failed to evaluate objective function.");
  return 0;
}

uno_int UnoNlp::objective_function_wrapper(uno_int /*number_variables*/, const double* x, double* objective_value, void* user_data)
{
  UnoNlp* nlp = static_cast<UnoNlp*>(user_data);
  return nlp->objective_function(x, objective_value);
}

uno_int UnoNlp::constraint_functions(const double* x, double* constraint_values)
{
  mem_->arg[0] = x;
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->res[0] = constraint_values;
  casadi_assert(mem_->self.calc_function(mem_, "nlp_g")==0, "Failed to evaluate constraints.");
  return 0;
}
  
uno_int UnoNlp::constraint_functions_wrapper(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
  double* constraint_values, void* user_data)
{
  UnoNlp* nlp = static_cast<UnoNlp*>(user_data);
  return nlp->constraint_functions(x, constraint_values);
}

// Evaluate the objective gradient
uno_int UnoNlp::objective_gradient(const double* x, double* gradient)
{
  mem_->arg[0] = x;
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->res[0] = gradient;
  casadi_assert(mem_->self.calc_function(mem_, "nlp_grad_f")==0, "Failed to evaluate gradient of objective.");
  return 0;
}
uno_int UnoNlp::objective_gradient_wrapper(uno_int /*number_variables*/, const double* x, double* gradient, void* user_data)
{
  UnoNlp* nlp = static_cast<UnoNlp*>(user_data);
  return nlp->objective_gradient(x, gradient);
}

// void CasadiModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
uno_int UnoNlp::jacobian(const double* x, double* jacobian_values) {
  // Sparsity pattern of jacobian is saved in mem->self ...
  // Evaluate numerically
    mem_->arg[0] = x;
    mem_->arg[1] = mem_->d_nlp.p;
    mem_->res[0] = jacobian_values;
    // mem_->res[0] = get_ptr(casadi_tmp_constraint_jacobian);
    casadi_assert(mem_->self.calc_function(mem_, "nlp_jac_g")==0, "Failed to evaluate constraint jacobian.");

    // // Write everything into RectangularMatrix...
    // // std::cout << "Number of variables: " << this->number_variables << std::endl;
    // // std::cout << "Number of vectors in constraint jacobian: " << constraint_jacobian.size() << std::endl;
    // size_t n_rows = constraint_jacobian.size();
    // size_t n_columns = this->number_variables;
    // //clear the matrix
    // for (size_t i=0; i<n_rows; ++i){
    //   constraint_jacobian[i].clear();
    // }

    // for (size_t j=0; j<n_columns; ++j) {
    //   for (size_t k= mem_->self.jacg_sp_.colind()[j]; k< mem_->self.jacg_sp_.colind()[j + 1];++k) {
    //     const size_t i = mem_->self.jacg_sp_.row()[k];
    //     const double entry = casadi_tmp_constraint_jacobian[k];
    //     constraint_jacobian[i].insert(j, entry);
    //   }
    // }
    return 0;
}

uno_int UnoNlp::jacobian_wrapper(uno_int /*number_variables*/, uno_int /*number_jacobian_nonzeros*/, const double* x, double* jacobian_values, void* user_data)
{
  UnoNlp* nlp = static_cast<UnoNlp*>(user_data);
  return nlp->jacobian(x, jacobian_values);
}

// void CasadiModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
//        SymmetricMatrix<double>& hessian) const {
  
//   // scale by the objective sign
//   objective_multiplier *= this->objective_sign;
//   // Evaluate numerically
//   mem_->arg[0] = get_ptr(x);
//   mem_->arg[1] = mem_->d_nlp.p;
//   mem_->arg[2] = &objective_multiplier;
//   // Multipliers need to be 
//   casadi_copy(get_ptr(multipliers), this->number_constraints, get_ptr(casadi_tmp_multipliers));
//   casadi_scal(this->number_constraints, -1., get_ptr(casadi_tmp_multipliers));
//   mem_->arg[3] = get_ptr(casadi_tmp_multipliers);
//   mem_->res[0] = get_ptr(casadi_tmp_hessian);
//   casadi_assert(mem_->self.calc_function(mem_, "nlp_hess_l")==0, "Failed to evaluate Lagrangian hessian.");


//   hessian.reset();
//   // Write the hessian into Symmetric matrix ....
//   for (size_t j=0; j<this->number_variables;++j) {
//     for (size_t k=mem_->self.hesslag_sp_.colind()[j]; k< mem_->self.hesslag_sp_.colind()[j + 1];++k) {
//        const size_t i = mem_->self.hesslag_sp_.row()[k];
//        if (i <= j){
//         const double entry = casadi_tmp_hessian[k];
//         hessian.insert(entry, i, j);
//        }
//     }
//     hessian.finalize_column(j);
//  }
// }

// double CasadiModel::get_variable_lower_bound(size_t i) const {
//   return this->variables_bounds[i].lb;//mem_->d_nlp.lbz[i];
// }

// double CasadiModel::get_variable_upper_bound(size_t i) const {
//   return this->variables_bounds[i].ub;//mem_->d_nlp.ubz[i];
// }

// BoundType CasadiModel::get_variable_bound_type(size_t i) const {
//   return this->variable_status[i];
// }

// double CasadiModel::get_constraint_lower_bound(size_t j) const {
//   return this->constraint_bounds[j].lb;//mem_->d_nlp.lbz[j+this->number_variables];
// }

// double CasadiModel::get_constraint_upper_bound(size_t j) const {
//   return this->constraint_bounds[j].ub;//mem_->d_nlp.ubz[j+this->number_variables];
// }

// FunctionType CasadiModel::get_constraint_type(size_t j) const {
//   return this->constraint_type[j];
// }

// BoundType CasadiModel::get_constraint_bound_type(size_t j) const {
//   return this->constraint_status[j];
// }

// size_t CasadiModel::get_number_objective_gradient_nonzeros() const {
//   return this->number_objective_gradient_nonzeros;;
// }
// size_t CasadiModel::get_number_jacobian_nonzeros() const {
//   return this->number_jacobian_nonzeros;
// }
// size_t CasadiModel::get_number_hessian_nonzeros() const {
//   return this->number_hessian_nonzeros;
// }

// void CasadiModel::get_initial_primal_point(std::vector<double>& x) const {
//   assert(x.size() >= this->number_variables);
//   std::copy(mem_->d_nlp.x0, mem_->d_nlp.x0 + this->number_variables, begin(x));
//   // std::cout << "Vector is " << x << std::endl;
// }
// void CasadiModel::get_initial_dual_point(std::vector<double>& multipliers) const {
//   assert(multipliers.size() >= this->number_constraints);
//  std::copy(mem_->d_nlp.lam_g0, mem_->d_nlp.lam_g0 + this->number_constraints, begin(multipliers));
// }
// void CasadiModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
// // do nothing
// }

// void CasadiModel::generate_variables() {
//   // Calculate variable bounds
//   for (size_t i=0; i< this->number_variables; ++i) {
//     double lb = (mem_->d_nlp.lbx != nullptr) ? *(mem_->d_nlp.lbx+i) : -INF<double>;
//     double ub = (mem_->d_nlp.ubx != nullptr) ? *(mem_->d_nlp.ubx+i) : INF<double>;
//     if (lb == ub) {
//         WARNING << "Variable x" << i << " has identical bounds\n";
//     }
//     this->variables_bounds[i] = {lb, ub};
//     // std::cout << "lbx at" << i << ":" << lb  << std::endl;
//     // std::cout << "ubx at" << i << ":" << ub << std::endl;
//   }

//   Model::determine_bounds_types(this->variables_bounds, this->variable_status);
//   // figure out the bounded variables
//   for (size_t i=0; i< this->number_variables; ++i) {
//     const BoundType status = this->get_variable_bound_type(i);
//     if (status == BOUNDED_LOWER || status == BOUNDED_BOTH_SIDES) {
//         this->lower_bounded_variables.push_back(i);
//         if (status == BOUNDED_LOWER) {
//           this->single_lower_bounded_variables.push_back(i);
//         }
//     }
//     if (status == BOUNDED_UPPER || status == BOUNDED_BOTH_SIDES) {
//         this->upper_bounded_variables.push_back(i);
//         if (status == BOUNDED_UPPER) {
//           this->single_upper_bounded_variables.push_back(i);
//         }
//     }
//   }
// }

// const std::vector<size_t>& CasadiModel::get_linear_constraints() const {
//     return this->linear_constraints;
// }

// void CasadiModel::generate_constraints() {
// //  auto d_nlp = &this->mem_->d_nlp;
//  for (size_t i=0; i< this->number_constraints; ++i) {
//     double lb = (mem_->d_nlp.lbg != nullptr) ? *(this->mem_->d_nlp.lbg+i) : -INF<double>;
//     double ub = (mem_->d_nlp.ubg != nullptr) ? *(this->mem_->d_nlp.ubg+i) : INF<double>;
//     // std::cout << "Memory location " << d_nlp->x0 << std::endl;
//     // double lb = d_nlp->x0[0];//+this->number_variables+i;
//     // double ub = *mem_->d_nlp.ubg;//+this->number_variables+i;
//     this->constraint_bounds[i] = {lb, ub};
//     // std::cout << "lbg at" << i << ":" << lb  << std::endl;
//     // std::cout << "ubg at" << i << ":" << ub << std::endl;

//  }
//  Model::determine_bounds_types(this->constraint_bounds, this->constraint_status);
  
//  for (size_t j=0; j<this->number_constraints; ++j) {
//     if (this->get_constraint_bound_type(j) == EQUAL_BOUNDS) {
//        this->equality_constraints.push_back(j);
//     }
//     else {
//        this->inequality_constraints.push_back(j);
//     }

//   // AMPL orders the constraints based on the function type: nonlinear first, then linear
//  for (size_t j=0; j<this->number_constraints; ++j) {
//     this->constraint_type[j] = NONLINEAR;
//  }
//  }
  
// //  this->determine_constraints();
// }

// void CasadiModel::set_function_types(std::string file_name) {
//   // // allocate a temporary ASL to read Hessian sparsity pattern
//   // ASL* asl_fg = ASL_alloc(ASL_read_fg);
//   // // char* stub = getstops(file_name, option_info);
//   // //if (file_name == nullptr) {
//   // //	usage_ASL(option_info, 1);
//   // //}

//   // FILE* nl = jac0dim_ASL(asl_fg, file_name.data(), static_cast<int>(file_name.size()));
//   // // specific read function
//   // // qp_read_ASL(asl_fg, nl, ASL_findgroups);

//   // // // constraints
//   // // if (asl_fg->i.n_con_ != static_cast<int>(this->number_constraints)) {
//   // //     throw std::length_error("AMPLModel.set_function_types: inconsistent number of constraints");
//   // // }

//   this->constraint_type.reserve(this->number_constraints);
//   this->problem_type = NONLINEAR;

//   // // determine the type of each constraint and objective function
//   // // determine if the problem is nonlinear (non-quadratic objective or nonlinear constraints)
//   // this->problem_type = LINEAR;
//   // int* rowq;
//   // int* colqp;
//   // double* delsqp;
//   // for (size_t j: Range(this->number_constraints)) {
//   //     int qp = nqpcheck_ASL(asl_fg, static_cast<int>(-(j + 1)), &rowq, &colqp, &delsqp);

//   //     if (0 < qp) {
//   //       this->constraint_type[j] = QUADRATIC;
//   //       this->problem_type = NONLINEAR;
//   //     }
//   //     else if (qp == 0) {
//   //       this->constraint_type[j] = LINEAR;
//   //       this->linear_constraints.push_back(j);
//   //     }
//   //     else {
//   //       this->constraint_type[j] = NONLINEAR;
//   //       this->problem_type = NONLINEAR;
//   //     }
//   // }
//   // // objective function
//   // int qp = nqpcheck_ASL(asl_fg, 0, &rowq, &colqp, &delsqp);
//   // if (0 < qp) {
//   //     if (this->problem_type == LINEAR) {
//   //       this->problem_type = QUADRATIC;
//   //     }
//   // }
//   // else if (qp != 0) {
//   //     this->problem_type = NONLINEAR;
//   // }
//   // qp_opify_ASL(asl_fg);

//   // // deallocate memory
//   // ASL_free(&asl_fg);
// }

} //namespace casadi