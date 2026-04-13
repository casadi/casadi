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

#include "Uno_C_API.h"
#include "uno_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <type_traits>

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
    for (auto&& op : casadi_options)
    {
      if (op.first=="preset")
      {
        preset = op.second.to_string();
        uno_set_solver_preset(solver, preset.c_str());
      }
      else
      {
        UnoNlp::set_uno_option(solver, op.first, op.second);
      }
    }
  }

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

    // define preset and insert options given through  
    insert_casadi_options(m->solver, opts_);
    
    // model creation
    const uno_int base_indexing = UNO_ZERO_BASED_INDEXING;
    // objective
    const uno_int optimization_sense = UNO_MINIMIZE;
    // constraints
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
    const uno_int number_hessian_nonzeros = static_cast<size_t>(this->hesslag_sp_.nnz());
    std::cout << "Hessian nnz: " << number_hessian_nonzeros << std::endl;
    std::vector<casadi_int> hess_row_indices = this->hesslag_sp_.get_row();
    std::vector<casadi_int> hess_column_indices = this->hesslag_sp_.get_col();
    std::vector<uno_int> hessian_row_indices(hess_row_indices.size());
    std::vector<uno_int> hessian_column_indices(hess_column_indices.size());

    for (uno_int i = 0; i < number_hessian_nonzeros; ++i) {
        hessian_row_indices[i] = static_cast<uno_int>(row_indices[i]);
        hessian_column_indices[i] = static_cast<uno_int>(column_indices[i]);
    }
    const char hessian_triangular_part = UNO_UPPER_TRIANGLE;
    const uno_int lagrangian_sign_convention = UNO_MULTIPLIER_POSITIVE;    
    
    // initial point
    // double x0[] = {-2., 1.};

    void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, nx_, d_nlp->lbz, d_nlp->ubz, base_indexing);
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
    
    ret = uno_set_constraints(model, ng_, UnoNlp::constraint_functions_wrapper,
          d_nlp->lbz+nx_, d_nlp->ubz+nx_, number_jacobian_nonzeros,
          jacobian_row_indices.data(), jacobian_column_indices.data(), UnoNlp::jacobian_wrapper);
    printf("uno_set_constraints returned: %d\n", ret);
    // if (ret != 0) {
    //     printf("ERROR: Failed to set constraints! Return code: %d\n", ret);
    //     return 1;
    // }

    // ret = uno_set_initial_primal_iterate(model, x0);
    ret = uno_set_initial_primal_iterate(model, d_nlp->x0);
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