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

uno_int UnoNlp::lagrangian_hessian(const double* x, double objective_multiplier, const double* multipliers, double* hessian_values)
{ 
  // scale by the objective sign
  // objective_multiplier *= this->objective_sign;
  // Evaluate numerically
  mem_->arg[0] = x;
  mem_->arg[1] = mem_->d_nlp.p;
  mem_->arg[2] = &objective_multiplier;
  // Multipliers need to be 
  // casadi_copy(get_ptr(multipliers), this->number_constraints, get_ptr(casadi_tmp_multipliers));
  // casadi_scal(this->number_constraints, -1., get_ptr(casadi_tmp_multipliers));
  mem_->arg[3] = multipliers;
  mem_->res[0] = hessian_values;
  casadi_assert(mem_->self.calc_function(mem_, "nlp_hess_l")==0, "Failed to evaluate Lagrangian hessian.");
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
  return 0;
}

uno_int UnoNlp::lagrangian_hessian_wrapper(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
            const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* user_data)
{
  UnoNlp* nlp = static_cast<UnoNlp*>(user_data);
  return nlp->lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
}

void UnoNlp::set_uno_option(void* solver, const std::string& name, const GenericType& value) {
  if (value.is_bool()) {
    uno_set_solver_bool_option(solver, name.c_str(), value.to_bool());
  } else if (value.is_int()) {
    uno_set_solver_integer_option(solver, name.c_str(), static_cast<uno_int>(value.to_int()));
  } else if (value.is_double()) {
    uno_set_solver_double_option(solver, name.c_str(), value.to_double());
  } else if (value.is_string()) {
    uno_set_solver_string_option(solver, name.c_str(), value.to_string().c_str());
  } else {
    casadi_assert(false, "Unsupported UNO option type for " + name);
  }
}

void UnoNlp::insert_casadi_options(void* solver, Dict opts)
{
    // build the (name, value) map
    Dict casadi_options = Options::sanitize(opts);
    
    // Define the preset and erase it from the options file
    std::string preset;
    for (auto&& op : casadi_options)
    {
      if (op.first=="preset")
      {
        preset = op.second.to_string();
        uno_set_solver_preset(solver, preset.c_str());
      }
      else
      {
        set_uno_option(solver, op.first, op.second);
      }
    }
}

} //namespace casadi