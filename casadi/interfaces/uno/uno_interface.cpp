/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

//Inclued UNO files
// #include "tools/Options.hpp"
// #include "tools/Timer.hpp"

// #include "optimization/Iterate.hpp"
// #include "optimization/ModelFactory.hpp"
// #include "optimization/ScaledModel.hpp"
// #include "preprocessing/Preprocessing.hpp"
// #include "linear_algebra/CSCSymmetricMatrix.hpp"

// #include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
// #include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
// #include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"

// #include "tools/Logger.hpp"

// #include "Uno.hpp"
#include "Uno_C_API.h"

// Casadi Includes
#include "uno_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include <type_traits>

extern "C" {

uno_int lagrangian_hessian(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
            const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* /*user_data*/) {
      hessian_values[0] = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
      hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0];
      hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1];
      return 0;
}

uno_int lagrangian_hessian_operator(uno_int number_variables, uno_int number_constraints, const double* x,
            bool evaluate_at_x, double objective_multiplier, const double* multipliers, const double* vector,
            double* result, void* user_data) {
      const double hessian00 = objective_multiplier*(1200*pow(x[0], 2.) - 400.*x[1] + 2.);
      const double hessian10 = -400.*objective_multiplier*x[0] - multipliers[0];
      const double hessian11 = 200.*objective_multiplier - 2.*multipliers[1];
      result[0] = hessian00*vector[0] + hessian10*vector[1];
      result[1] = hessian10*vector[0] + hessian11*vector[1];
      return 0;
}

} // extern "C"

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_UNO_EXPORT
  casadi_register_nlpsol_uno(Nlpsol::Plugin* plugin) {
    plugin->creator = UnoInterface::creator;
    plugin->name = "uno";
    plugin->doc = UnoInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &UnoInterface::options_;
    plugin->deserialize = &UnoInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_UNO_EXPORT casadi_load_nlpsol_uno() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_uno);
  }

  /*---------------------------------------------
  Constructor Destructor UnoInterface
  ----------------------------------------------*/

  UnoInterface::UnoInterface(const std::string& name, const Function& nlp) : Nlpsol(name, nlp)
  {
  }

  UnoInterface::~UnoInterface()
  {
    clear_mem();
  }

  const Options UnoInterface::options_
  = {{&Nlpsol::options_},
     {{"uno",
       {OT_DICT,
        "Options to be passed to UNO"}}
     }
  };

  UnoMemory::UnoMemory(const UnoInterface& uno_interface) : self(uno_interface), NlpsolMemory()
  {
    this->return_status = "Unset";
  }

  UnoMemory::~UnoMemory()
  {
  }

  /*----------------------------------------------------------------
  From here CasadiModel function definition
  ---------------------------------------------------------------*/

  // CasadiModel::CasadiModel(const std::string& file_name, const UnoInterface& uno_interface, UnoMemory* mem) :
  //   Model(file_name, uno_interface.nx_, uno_interface.ng_),
  //   mem_(mem),
  //   // allocate vectors
  //   casadi_tmp_gradient(this->number_variables),
  //   casadi_tmp_multipliers(this->number_constraints),
  //   casadi_tmp_constraint_jacobian(mem->self.get_function("nlp_jac_g").sparsity_out(0).nnz()),
  //   casadi_tmp_hessian(mem->self.get_function("nlp_hess_l").sparsity_out(0).nnz()),
  //   variables_bounds(this->number_variables),
  //   constraint_bounds(this->number_constraints),
  //   variable_status(this->number_variables),
  //   constraint_type(this->number_constraints),
  //   constraint_status(this->number_constraints) {
    
  //   // this->asl->i.congrd_mode = 0;

  //   // dimensions
  //   this->objective_sign = 1.;//(this->asl->i.objtype_[0] == 1) ? -1. : 1.;

  //   // variables
  //   this->generate_variables();

  //   // constraints
  //   this->equality_constraints.reserve(this->number_constraints);
  //   this->inequality_constraints.reserve(this->number_constraints);
  //   this->linear_constraints.reserve(this->number_constraints);
  //   this->generate_constraints();
  //   // this->set_function_types(file_name);

  //   // compute number of nonzeros
  //   this->number_objective_gradient_nonzeros = static_cast<size_t>(0);//static_cast<size_t>(this->asl->i.nzo_);
  //   this->number_jacobian_nonzeros = static_cast<size_t>(mem_->self.jacg_sp_.nnz());
  //   this->set_number_hessian_nonzeros();

  // }

//   void CasadiModel::set_number_hessian_nonzeros() {
//    // compute the maximum number of nonzero elements, provided that all multipliers are non-zero
//    // int (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri);
//    const int objective_number = -1;
//    const int upper_triangular = 1;
//   //  this->hessian_maximum_number_nonzeros = static_cast<size_t>((*(this->asl)->p.Sphset)(this->asl, nullptr, objective_number, 1, 1,
//   //        upper_triangular));
//    this->number_hessian_nonzeros = static_cast<size_t>(this->mem_->self.hesslag_sp_.nnz());
//    this->casadi_tmp_hessian.reserve(this->number_hessian_nonzeros);

//    // use Lagrangian scale: in AMPL, the Lagrangian is f + lambda.g, while Uno uses f - lambda.g
//    int nerror{};
//   //  lagscale_ASL(this->asl, -1., &nerror);
// }

  /*------------------------------------------
  UnoInterface function definitions
  -------------------------------------------*/

  void UnoInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="uno") {
        opts_ = op.second;
      }

    }

    // Setup NLP functions
    // create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    // Function gf_jg_fcn = create_function("nlp_gf_jg", {"x", "p"}, {"grad:f:x", "jac:g:x"});
    // jacg_sp_ = gf_jg_fcn.sparsity_out(1);

    // Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
    //                               {"hess:gamma:x:x"},
    //                               {{"gamma", {"f", "g"}}});
    // hesslag_sp_ = hess_l_fcn.sparsity_out(0);

    create_function("nlp_f", {"x", "p"}, {"f"});
    create_function("nlp_g", {"x", "p"}, {"g"});
    create_function("nlp_grad_f", {"x", "p"}, {"grad:f:x"});
    Function gf_jg_fcn = create_function("nlp_jac_g", {"x", "p"}, {"jac:g:x"});
    jacg_sp_ = gf_jg_fcn.sparsity_out(0);

    Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                  {"triu:hess:gamma:x:x"},
                                  {{"gamma", {"f", "g"}}});
    hesslag_sp_ = hess_l_fcn.sparsity_out(0);
    casadi_assert(hesslag_sp_.is_triu(), "Hessian must be upper triangular");

    // Allocate persistent memory
    alloc_w(nx_, true); // wlbx_
    alloc_w(nx_, true); // wubx_
    alloc_w(ng_, true); // wlbg_
    alloc_w(ng_, true); // wubg_
  }

  int UnoInterface::init_mem(void* mem) const {
    int return_nlpsol = Nlpsol::init_mem(mem);
    if (return_nlpsol > 0){
      printf("Return code: %d!", return_nlpsol);
    }
    auto m = static_cast<UnoMemory*>(mem);

// -------------------------------------------------------------------
    uno_int uno_major, uno_minor, uno_patch;
    uno_get_version(&uno_major, &uno_minor, &uno_patch);
    printf("Uno v%d.%d.%d\n", uno_major, uno_minor, uno_patch);

    m->uno_nlp = new UnoNlp(m);

   // Casadi model
   // memory should be freed somewhere else
  //  std::cout << "Init Memort acces" << std::endl;
  //  m->model = new CasadiModel("casadi_model", *this, m);
    // m->model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing);
    m->solver = uno_create_solver();
   return 0;

  }
//----------------------------------------------------------

  // }

  void UnoInterface::set_work(void* mem, const double**& arg, double**& res,
                                 casadi_int*& iw, double*& w) const
  {
    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);
    auto m = static_cast<UnoMemory*>(mem);
  }

  int casadi_KN_puts(const char * const str, void * const userParams) {
    std::string s(str);
    uout() << s << std::flush;
    return s.size();
  }

  void insert_casadi_options(void* solver, Dict opts) {
    // build the (name, value) map
    Dict casadi_options = Options::sanitize(opts);
    
    // Define the preset and erase it from the options file
    std::string preset;
    auto it = casadi_options.find("preset");
    if (it!=casadi_options.end())
    {
      preset = it->second.to_string();
      uno_set_solver_preset(solver, preset.c_str());
      casadi_options.erase(it);
    }
    uno_set_solver_bool_option(solver, "print_solution", true);

    // solver creation
    // uno_set_solver_bool_option(m->solver, "print_solution", true);

    // // Pass all the options to Uno
    // for (auto&& op : casadi_options) {

    //   // There might be options with a resto prefix.
    //   std::string option_name = op.first;
    //   if (startswith(option_name, "resto.")) {
    //     option_name = option_name.substr(6);
    //   }
    //   const std::string name = op.first;
    //   const std::string value = std::string(op.second);
    //   uno_options[name] = value;
    // }
  }

  // Statistics create_statistics(const Model& model, const ::Options& options) {
  //  Statistics statistics(options);
  //  statistics.add_column("iters", Statistics::int_width, options.get_int("statistics_major_column_order"));
  //  statistics.add_column("step norm", Statistics::double_width, options.get_int("statistics_step_norm_column_order"));
  //  statistics.add_column("objective", Statistics::double_width, options.get_int("statistics_objective_column_order"));
  //  if (model.is_constrained()) {
  //     statistics.add_column("primal infeas.", Statistics::double_width, options.get_int("statistics_primal_infeasibility_column_order"));
  //  }
  //  statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
  //  statistics.add_column("stationarity", Statistics::double_width, options.get_int("statistics_stationarity_column_order"));
  //  return statistics;
  // }

inline const char* return_status_string(void* solver) {
    uno_int iterate_status = uno_get_solution_status(solver);
    assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
    if (iterate_status == UNO_FEASIBLE_KKT_POINT)
    {
      return "Converged with feasible KKT point";
    }
    else if (iterate_status == UNO_FEASIBLE_FJ_POINT)
    {
      return "Converged with feasible FJ point";
    }
    else if (iterate_status == UNO_INFEASIBLE_STATIONARY_POINT)
    {
      return "Converged with infeasible stationary point";
    }
    else if (iterate_status == UNO_FEASIBLE_SMALL_STEP)
    {
      return "Terminated with feasible small step";
    }
    else if (iterate_status == UNO_INFEASIBLE_SMALL_STEP)
    {
      return "Terminated with infeasible small step";
    }
    else if (iterate_status == UNO_UNBOUNDED)
    {
      return "Terminated with unbounded problem";
    }
    else if (iterate_status == UNO_NOT_OPTIMAL)
    {
      return "Terminated with not optimal point";
    }
    else
    {
      return "Terminated with an unknown status!";
    }
  }

  inline const bool return_status_success(void* solver)
  {
    uno_int optimization_status = uno_get_optimization_status(solver);
    if (optimization_status == UNO_SUCCESS)
    {
        return true;
    }
    else
    {
        return false;
    }
  }

  int UnoInterface::solve(void* mem) const {
    auto m = static_cast<UnoMemory*>(mem);
    auto d_nlp = &m->d_nlp;

    // define preset and insert options given through casadi 
    insert_casadi_options(m->solver, opts_);

    // //   // create the statistics
    // // Statistics statistics = create_statistics(*model, uno_options);

    //     // Copy optimal constraint values to output
    //     casadi_copy(get_ptr(result.solution.evaluations.constraints), ng_, d_nlp->z + nx_);
    
    // model creation
    const uno_int base_indexing = UNO_ZERO_BASED_INDEXING;
    // variables
    // const uno_int number_variables = 2;
    double variables_lower_bounds[] = {-INFINITY, -INFINITY};
    double variables_upper_bounds[] = {0.5, INFINITY};
    // objective
    const uno_int optimization_sense = UNO_MINIMIZE;
    // constraints
    // const uno_int number_constraints = 2;
    // const uno_int number_jacobian_nonzeros = 4;
    // uno_int jacobian_row_indices[] = {0, 1, 0, 1};
    // uno_int jacobian_column_indices[] = {0, 0, 1, 1};
    double constraints_lower_bounds[] = {1., 0.};
    double constraints_upper_bounds[] = {INFINITY, INFINITY};
    const uno_int number_jacobian_nonzeros = static_cast<size_t>(this->jacg_sp_.nnz());
    std::cout << "Jacobian nnz: " << number_jacobian_nonzeros << std::endl;
    std::vector<casadi_int> row_indices = this->jacg_sp_.get_row();
    std::vector<casadi_int> column_indices = this->jacg_sp_.get_col();
    std::vector<uno_int> jacobian_row_indices(row_indices.size());
    std::vector<uno_int> jacobian_column_indices(column_indices.size());

    for (uno_int i = 0; i < number_jacobian_nonzeros; ++i) {
        jacobian_row_indices[i] = static_cast<uno_int>(row_indices[i]);
        jacobian_column_indices[i] = static_cast<uno_int>(column_indices[i]);
    }

    // Hessian
    // const uno_int number_hessian_nonzeros = 3;
    const uno_int number_hessian_nonzeros = static_cast<size_t>(this->hesslag_sp_.nnz_lower());
    std::cout << "Hessian nnz: " << number_hessian_nonzeros << std::endl;
    std::vector<casadi_int> hess_row_indices = this->hesslag_sp_.get_row();
    std::vector<casadi_int> hess_column_indices = this->hesslag_sp_.get_col();
    std::vector<uno_int> hessian_row_indices(hess_row_indices.size());
    std::vector<uno_int> hessian_column_indices(hess_column_indices.size());

    for (uno_int i = 0; i < number_hessian_nonzeros; ++i) {
        hessian_row_indices[i] = static_cast<uno_int>(row_indices[i]);
        hessian_column_indices[i] = static_cast<uno_int>(column_indices[i]);
    }
    // const char hessian_triangular_part = UNO_LOWER_TRIANGLE;
    const char hessian_triangular_part = UNO_UPPER_TRIANGLE;
    const uno_int lagrangian_sign_convention = UNO_MULTIPLIER_POSITIVE;
    // uno_int hessian_row_indices[] = {0, 1, 1};
    // uno_int hessian_column_indices[] = {0, 0, 1};
    
    
    // initial point
    double x0[] = {-2., 1.};

    // void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
    void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, nx_, variables_lower_bounds,
    variables_upper_bounds, base_indexing);
    std::cout << "Set the objective here!" << std::endl;

    UnoNlp* nlp = static_cast<UnoNlp*>(m->uno_nlp);
    
    uno_int ret;
    ret = uno_set_user_data(model, nlp);
    ret = uno_set_objective(model, optimization_sense, UnoNlp::objective_function_wrapper, UnoNlp::objective_gradient_wrapper);
    printf("uno_set_objective returned: %d\n", ret);
    // if (ret != 0) {
    //     printf("ERROR: Failed to set objective! Return code: %d\n", ret);
    //     return 1;
    // }
    
    // ret = uno_set_constraints(model, number_constraints, UnoNlp::constraint_functions_wrapper,
    //       constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
    //       jacobian_row_indices, jacobian_column_indices, UnoNlp::jacobian_wrapper);
    ret = uno_set_constraints(model, ng_, UnoNlp::constraint_functions_wrapper,
          constraints_lower_bounds, constraints_upper_bounds, number_jacobian_nonzeros,
          jacobian_row_indices.data(), jacobian_column_indices.data(), UnoNlp::jacobian_wrapper);
    printf("uno_set_constraints returned: %d\n", ret);
    // if (ret != 0) {
    //     printf("ERROR: Failed to set constraints! Return code: %d\n", ret);
    //     return 1;
    // }

    ret = uno_set_initial_primal_iterate(model, x0);
    printf("uno_set_initial_primal_iterate returned: %d\n", ret);
    // if (ret != 0) {
    //     printf("ERROR: Failed to set initial primal iterate! Return code: %d\n", ret);
    //     return 1;
    // }

    ret = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices.data(), hessian_column_indices.data(), UnoNlp::lagrangian_hessian_wrapper);
    // printf("uno_set_lagrangian_hessian returned: %d\n", ret);
    // ret = uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention);
    // printf("uno_set_lagrangian_sign_convention returned: %d\n", ret);

    // run 1: solve with no Hessian. Uno defaults to L-BFGS Hessian for NLPs
    uno_optimize(m->solver, model);
    // get the solution
    uno_int optimization_status = uno_get_optimization_status(m->solver);
    assert(optimization_status == UNO_SUCCESS);
    uno_int iterate_status = uno_get_solution_status(m->solver);
    assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
    double solution_objective = uno_get_solution_objective(m->solver);
    printf("Solution objective = %g\n", solution_objective);

    // run 2: solve with exact Hessian
    ret = uno_set_lagrangian_hessian(model, number_hessian_nonzeros, hessian_triangular_part, hessian_row_indices.data(), hessian_column_indices.data(), UnoNlp::lagrangian_hessian_wrapper);
    printf("uno_set_lagrangian_hessian returned: %d\n", ret);
    ret = uno_set_lagrangian_sign_convention(model, lagrangian_sign_convention);
    printf("uno_set_lagrangian_sign_convention returned: %d\n", ret);
    uno_optimize(m->solver, model);
    // get the solution
    optimization_status = uno_get_optimization_status(m->solver);
    assert(optimization_status == UNO_SUCCESS);
    iterate_status = uno_get_solution_status(m->solver);
    assert(iterate_status == UNO_FEASIBLE_KKT_POINT);
    solution_objective = uno_get_solution_objective(m->solver);
    printf("Solution objective = %g\n", solution_objective);
    uno_get_primal_solution(m->solver, d_nlp->z);
    // Get dual solution (constraints)
    uno_get_constraint_dual_solution(m->solver, d_nlp->lam+nx_);
    // Get dual solution (simple bounds)
    for (casadi_int i=0; i<nx_; ++i) {
      d_nlp->lam[i] = uno_get_upper_bound_dual_solution_component(m->solver, i)-uno_get_lower_bound_dual_solution_component(m->solver, i);
    }

    // Write the solution to Casadi .....
    // Negate rc to match CasADi's definition
    m->return_status = return_status_string(m->solver);
    m->success = return_status_success(m->solver);
    // Get optimal cost
    d_nlp->objective = uno_get_solution_objective(m->solver);

    const double solution_primal_feasibility = uno_get_solution_primal_feasibility(m->solver);
    printf("Primal feasibility s solution = %e\n", solution_primal_feasibility);
    const double solution_stationarity = uno_get_solution_stationarity(m->solver);
    printf("Stationarity at solution = %e\n", solution_stationarity);
    const double solution_complementarity = uno_get_solution_complementarity(m->solver);
    printf("Complementarity at solution = %e\n", solution_complementarity);

    // cleanup
    uno_destroy_solver(m->solver);
    uno_destroy_model(model);

    return 0;
  }

  Dict UnoInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    
    return stats;
  }

} // namespace casadi