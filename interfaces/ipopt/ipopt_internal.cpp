/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "ipopt_internal.hpp"
#include "ipopt_nlp.hpp"
#include "casadi/stl_vector_tools.hpp"

using namespace std;
#include <coin/IpIpoptApplication.hpp>
namespace CasADi{

IpoptInternal::IpoptInternal(const FX& F_, const FX& G_, const FX& H_, const FX& J_, const FX& GF_) : NLPSolverInternal(F_,G_,H_,J_,GF_){
  
  // Output
  ops_["print_level"] = OT_INTEGER;
  ops_["print_user_options"] = OT_STRING;
  ops_["print_options_documentation"] = OT_STRING;
  ops_["output_file"] = OT_STRING;
  ops_["file_print_level"] = OT_INTEGER;
  ops_["option_file_name"] = OT_STRING;
    
  // Termination
  ops_["tol"] = OT_REAL;
  ops_["max_iter"] = OT_INTEGER;
  ops_["max_cpu_time"] = OT_REAL;
  ops_["dual_inf_tol"] = OT_REAL;
  ops_["constr_viol_tol"] = OT_REAL;
  ops_["compl_inf_tol"] = OT_REAL;
  ops_["acceptable_tol"] = OT_REAL;
  ops_["acceptable_iter"] = OT_INTEGER;
  ops_["acceptable_constr_viol_tol"] = OT_INTEGER;
  ops_["acceptable_dual_inf_tol"] = OT_REAL;
  ops_["acceptable_compl_inf_tol"] = OT_REAL;
  ops_["acceptable_obj_change_tol"] = OT_REAL;
  ops_["diverging_iterates_tol"] = OT_REAL;
  
  // NLP
  ops_["bound_relax_factor"] = OT_REAL;
  ops_["honor_original_bounds"] = OT_STRING;
  ops_["check_derivatives_for_naninf"] = OT_STRING;
  ops_["nlp_lower_bound_inf"] = OT_REAL;
  ops_["nlp_upper_bound_inf"] = OT_REAL;
  ops_["fixed_variable_treatment"] = OT_STRING;
  ops_["jac_c_constant"] = OT_STRING;
  ops_["jac_d_constant"] = OT_STRING;
  ops_["hessian_constant"] = OT_STRING;
  
  // Initialization
  ops_["bound_frac"] = OT_REAL;
  ops_["bound_push"] = OT_REAL;
  ops_["slack_bound_frac"] = OT_REAL;
  ops_["slack_bound_push"] = OT_REAL;
  ops_["bound_mult_init_val"] = OT_REAL;
  ops_["constr_mult_init_max"] = OT_REAL;
  ops_["bound_mult_init_method"] = OT_STRING;

  // Barrier Parameter
  ops_["mehrotra_algorithm"] = OT_STRING;
  ops_["mu_strategy"] = OT_STRING;
  ops_["mu_oracle"] = OT_STRING;
  ops_["quality_function_max_section_steps"] = OT_INTEGER;
  ops_["fixed_mu_oracle"] = OT_STRING;
  ops_["mu_init"] = OT_REAL;
  ops_["mu_max_fact"] = OT_REAL;
  ops_["mu_max"] = OT_REAL;
  ops_["mu_min"] = OT_REAL;
  ops_["mu_target"] = OT_REAL;
  ops_["barrier_tol_factor"] = OT_REAL;
  ops_["mu_linear_decrease_factor"] = OT_REAL;
  ops_["mu_superlinear_decrease_power"] = OT_REAL;
  
  // Multiplier Updates
  ops_["alpha_for_y"] = OT_STRING;
  ops_["alpha_for_y_tol"] = OT_REAL;
  ops_["recalc_y"] = OT_STRING;
  ops_["recalc_y_feas_tol"] = OT_REAL;

  // Line Search
  ops_["max_soc"] = OT_INTEGER;
  ops_["watchdog_shortened_iter_trigger"] = OT_INTEGER;
  ops_["watchdog_trial_iter_max"] = OT_INTEGER;
  ops_["accept_every_trial_step"] = OT_STRING;
  ops_["corrector_type"] = OT_STRING;

  // Warm Start
  ops_["warm_start_init_point"] = OT_STRING;
  ops_["warm_start_bound_push"] = OT_REAL;
  ops_["warm_start_bound_frac"] = OT_REAL;
  ops_["warm_start_slack_bound_frac"] = OT_REAL;
  ops_["warm_start_slack_bound_push"] = OT_REAL;
  ops_["warm_start_mult_bound_push"] = OT_REAL;
  ops_["warm_start_mult_init_max"] = OT_REAL;

  // Restoration Phase
  ops_["expect_infeasible_problem"] = OT_STRING;
  ops_["expect_infeasible_problem_ctol"] = OT_REAL;
  ops_["expect_infeasible_problem_ytol"] = OT_REAL;
  ops_["start_with_resto"] = OT_STRING;
  ops_["soft_resto_pderror_reduction_factor"] = OT_REAL;
  ops_["required_infeasibility_reduction"] = OT_REAL;
  ops_["bound_mult_reset_threshold"] = OT_REAL;
  ops_["constr_mult_reset_threshold"] = OT_REAL;
  ops_["evaluate_orig_obj_at_resto_trial"] = OT_STRING;

  // Linear Solver
  ops_["linear_solver"] = OT_STRING;
  ops_["linear_system_scaling"] = OT_STRING;
  ops_["linear_scaling_on_demand"] = OT_STRING;
  ops_["max_refinement_steps"] = OT_INTEGER;
  ops_["min_refinement_steps"] = OT_INTEGER;

// Hessian Perturbation
ops_["max_hessian_perturbation"] = OT_REAL;
ops_["min_hessian_perturbation"] = OT_REAL;
ops_["first_hessian_perturbation"] = OT_REAL;
ops_["perturb_inc_fact_first"] = OT_REAL;
ops_["perturb_inc_fact"] = OT_REAL;
ops_["perturb_dec_fact"] = OT_REAL;
ops_["jacobian_regularization_value"] = OT_REAL;

// Quasi-Newton
ops_["hessian_approximation"] = OT_STRING;
ops_["limited_memory_max_history"] = OT_INTEGER;
ops_["limited_memory_max_skipping"] = OT_INTEGER;

// Derivative Test
ops_["derivative_test"] = OT_STRING;
ops_["derivative_test_perturbation"] = OT_REAL;
ops_["derivative_test_tol"] = OT_REAL;
ops_["derivative_test_print_all"] = OT_STRING;
ops_["point_perturbation_radius"] = OT_REAL;

// MA27 Linear Solver
ops_["ma27_pivtol"] = OT_REAL;
ops_["ma27_pivtolmax"] = OT_REAL;
ops_["ma27_liw_init_factor"] = OT_REAL;
ops_["ma27_la_init_factor"] = OT_REAL;
ops_["ma27_meminc_factor"] = OT_REAL;

// MA57 Linear Solver
ops_["ma57_pivtol"] = OT_REAL;
ops_["ma57_pivtolmax"] = OT_REAL;
ops_["ma57_pre_alloc"] = OT_REAL;
ops_["ma57_pivot_order"] = OT_INTEGER;
ops_["ma57_automatic_scaling"] = OT_STRING;
ops_["ma57_block_size"] = OT_INTEGER;
ops_["ma57_node_amalgamation"] = OT_INTEGER;
ops_["ma57_small_pivot_flag" ] = OT_INTEGER;

// MUMPS Linear Solver
ops_["mumps_pivtol"] = OT_REAL;
ops_["mumps_pivtolmax"] = OT_REAL;
ops_["mumps_mem_percent"] = OT_INTEGER;
ops_["mumps_permuting_scaling"] = OT_INTEGER;
ops_["mumps_pivot_order"] = OT_INTEGER;
ops_["mumps_scaling"] = OT_INTEGER;

// Pardiso Linear Solver
ops_["pardiso_msglvl"] = OT_INTEGER;
ops_["pardiso_matching_strategy"] = OT_STRING;
ops_["pardiso_out_of_core_power"] = OT_INTEGER;

// WSMP Linear Solver
ops_["wsmp_num_threads"] = OT_INTEGER;
ops_["wsmp_ordering_option"] = OT_INTEGER;
ops_["wsmp_pivtol"] = OT_REAL;
ops_["wsmp_pivtolmax"] = OT_REAL;
ops_["wsmp_scaling"] = OT_INTEGER;
ops_["wsmp_singularity_threshold"] = OT_REAL;

// Add to options structure
for(map<string,opt_type>::const_iterator it=ops_.begin(); it!=ops_.end(); ++it)
  addOption(it->first,it->second);

  app = 0;
  userclass = 0;

  // Start the application
  app = new Ipopt::IpoptApplication();

  // Create an Ipopt user class -- need to use Ipopts spart pointer class
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = new Ipopt::SmartPtr<Ipopt::TNLP>();
  userclass = (void*)ucptr;
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;
  uc = new IpoptUserClass(this);
  
  // Limited memory hessian by default if no hessian provided
  if(H_.isNull())
    setOption("hessian_approximation","limited-memory");
}


IpoptInternal::~IpoptInternal(){
  if(app) delete app;

  // delete the smart pointer;
  if(userclass != 0){
    Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = (Ipopt::SmartPtr<Ipopt::TNLP>*)userclass;
    delete ucptr;
  }
}

void IpoptInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();

  // read options
  exact_hessian_ = !H_.isNull();
  
  if(verbose_){
    cout << "There are " << n_ << " variables and " << m_ << " constraints." << endl;
    if(exact_hessian_) std::cout << "Using exact Hessian" << std::endl;
    else             std::cout << "Using limited memory Hessian approximation" << std::endl;
  }
 
  // Pass all the options to ipopt
  for(map<string,opt_type>::const_iterator it=ops_.begin(); it!=ops_.end(); ++it)
    if(hasSetOption(it->first)){
      Option op = getOption(it->first);
      switch(it->second){
        case OT_REAL:
          app->Options()->SetNumericValue(it->first,op.toDouble());
          break;
        case OT_INTEGER:
          app->Options()->SetIntegerValue(it->first,op.toInt());
          break;
        case OT_STRING:
          app->Options()->SetStringValue(it->first,op.toString());
          break;
        default:
          throw CasadiException("Illegal type");
      }
    }
  
  // Intialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status = app->Initialize();
  if (status != Solve_Succeeded) {
    throw "Error during initialization!\n";
  }
}

void IpoptInternal::evaluate(int fsens_order, int asens_order){
  assert(fsens_order==0 && asens_order==0);

  // Get back the smart pointer
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = (Ipopt::SmartPtr<Ipopt::TNLP>*)userclass;
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;

  // Ask Ipopt to solve the problem
  Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(uc);

  if (status == Solve_Succeeded)
    std::cout << "*** The problem solved!" << std::endl;
  else
    std::cout << "*** The problem FAILED" << std::endl;

}

void IpoptInternal::finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value){
  copy(x,x+n_,result(NLP_X_OPT).begin());
  copy(z_L,z_L+n_,result(NLP_LAMBDA_LBX).begin());
  copy(z_U,z_U+n_,result(NLP_LAMBDA_UBX).begin());
  copy(lambda,lambda+m_,result(NLP_LAMBDA_OPT).begin());
  result(NLP_COST).at(0) = obj_value;
}

bool IpoptInternal::eval_h(const double* x, bool new_x, double obj_factor, const double* lambda,bool new_lambda, int nele_hess, int* iRow,int* jCol, double* values){
  log("eval_h started");
  if (values == NULL) {
    int nz=0;
    vector<int> rowind,col;
    H_.result().sparsity().getSparsityCRS(rowind,col);
    for(int r=0; r<rowind.size()-1; ++r)
      for(int el=rowind[r]; el<rowind[r+1]; ++el){
//        if(col[el]>=r){
          iRow[nz] = r;
          jCol[nz] = col[el];
          nz++;
  //      }
      }
  } else {
    // Pass input
    H_.setInput(x);
    H_.setInput(lambda,1);
    H_.setInput(obj_factor,2);

    // Evaluate
    H_.evaluate();

    // Get results
    H_.getOutput(values);
  }
  log("eval_h ok");
  return true;
}

bool IpoptInternal::eval_jac_g(int n, const double* x, bool new_x,int m, int nele_jac, int* iRow, int *jCol,double* values){
  try{
    log("eval_jac_g started");
    if (values == NULL) {
      int nz=0;
      vector<int> rowind,col;
      J_.result().sparsity().getSparsityCRS(rowind,col);
      for(int r=0; r<rowind.size()-1; ++r)
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
  //        if(col[el]>=r){
            iRow[nz] = r;
            jCol[nz] = col[el];
            nz++;
    //      }
        }
    } else {
      // Pass the argument to the function
      J_.setInput(x);
      
       // Evaluate the function
      J_.evaluate();

      // Get the output
      J_.getOutput(values);
      
      if(monitored("eval_jac_g")){
        cout << "J = " << endl;
        J_.result().printSparse();
      }
    }
    
    log("eval_jac_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_jac_g failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_f(int n, const double* x, bool new_x, double& obj_value)
{
  log("eval_f started");
  
  assert(n == n_);

  // Pass the argument to the function
  F_.setInput(x);

  // Evaluate the function
  F_.evaluate();

  // Get the result
  F_.getOutput(obj_value);

  // Printing
  if(monitored("eval_f")){
    cout << "obj_value = " << obj_value << endl;
  }

  log("eval_f ok");
  return true;
}

bool IpoptInternal::eval_g(int n, const double* x, bool new_x, int m, double* g)
{
  log("eval_g started");

  assert(n == n_);
  assert(m == m_);

  // Pass the argument to the function
  G_.setInput(x);

  // Evaluate the function and tape
  G_.evaluate();

  // Ge the result
  G_.getOutput(g);

  // Printing
  if(monitored("eval_g"))
    cout << "g = " << G_.result() << endl;
    
  log("eval_g ok");
  return true;
}

bool IpoptInternal::eval_grad_f(int n, const double* x, bool new_x, double* grad_f)
{
  log("eval_grad_f started");

  assert(n == n_);
  // Pass the argument to the function
  F_.setInput(x);
  
  // Give a seed to the function
  F_.setAdjSeed(1.0);

  // Evaluate, adjoint mode
  F_.evaluate(0,1);

  // Get the result
  F_.getAdjSens(grad_f);

  // Printing
  if(monitored("eval_grad_f")){
    cout << "grad_f = " << F_.adjSens() << endl;
  }
  
  // Check the result for regularity
  for(vector<double>::const_iterator it=F_.adjSens().begin(); it!=F_.adjSens().end(); ++it){
      if(isnan(*it) || isinf(*it)){
        log("eval_grad_f: result not regular");
        return false;
    }
  }
  
  log("eval_grad_f ok");
  return true;
}

bool IpoptInternal::get_bounds_info(int n, double* x_l, double* x_u,
                                int m, double* g_l, double* g_u)
{
  assert(n == n_);
  assert(m == m_);
  vector<double> &lbx = input(NLP_LBX);  copy(lbx.begin(),lbx.end(),x_l);
  vector<double> &ubx = input(NLP_UBX);  copy(ubx.begin(),ubx.end(),x_u);
  vector<double> &lbg = input(NLP_LBG);  copy(lbg.begin(),lbg.end(),g_l);
  vector<double> &ubg = input(NLP_UBG);  copy(ubg.begin(),ubg.end(),g_u);
  return true;
}

bool IpoptInternal::get_starting_point(int n, bool init_x, double* x,
                                   bool init_z, double* z_L, double* z_U,
                                   int m, bool init_lambda,
                                   double* lambda)
{

  // MISSING: Starting values for the dual variables
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  const vector<double> &xinit = input(NLP_X_INIT);
  copy(xinit.begin(),xinit.end(),x);
  return true;
}

void IpoptInternal::get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag)
{


  n = n_;               // number of variables
  m = m_;               // number of constraints

  // Get Jacobian sparsity pattern
  if(G_.isNull())
    nnz_jac_g = 0;
  else
    nnz_jac_g = J_.result().size();

  // Get Hessian sparsity pattern
  if(H_.isNull())
    nnz_h_lag = 0;
  else
    nnz_h_lag = H_.result().size();
}

} // namespace CasADi
