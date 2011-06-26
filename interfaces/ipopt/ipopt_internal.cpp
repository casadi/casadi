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
#include <ctime>

using namespace std;
#include <coin/IpIpoptApplication.hpp>
namespace CasADi{

IpoptInternal::IpoptInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : F_(F), G_(G), H_(H), J_(J), GF_(GF){
  addOption("pass_nonlinear_variables", OT_BOOLEAN, true);
  addOption("print_time", OT_BOOLEAN, true, "print information about execution time");
  
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
  ops_["acceptable_constr_viol_tol"] = OT_REAL;
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

  // Limited memory hessian by default if no hessian provided
  if(H_.isNull()) setOption("hessian_approximation","limited-memory");

  app = 0;
  userclass = 0;
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
  // Initialize the functions
  F_.init();
  if(!G_.isNull()) G_.init();
  n_ = F_.input(0).numel();
  m_ = G_.isNull() ? 0 : G_.output(0).numel();
  
  int n=0;
  int m=0;
  
  pn_=0;
  pm_=0;
  
  int pn=0;
  int pm=0;
  
  // Basic sanity checks
  casadi_assert_message(F_.getNumInputs()==1 or F_.getNumInputs()==2, "Wrong number of input arguments to F. Must ");
  casadi_assert_message(F_.getNumOutputs()>=1, "Wrong number of output arguments to F");
  casadi_assert_message(F_.output().scalar() && F_.output().dense(), "Output argument of F not dense scalar.");
  
  n=F_.input().numel();
  if(!G_.isNull()) {
    casadi_assert_message(G_.getNumInputs()>=1, "Wrong number of input arguments to G");
    casadi_assert_message(G_.getNumOutputs()>=1, "Wrong number of output arguments to G");
    casadi_assert_message(G_.input().numel()==n, "Inconsistent dimensions");
    m=G_.output().numel();
  }
  
  if(!H_.isNull()) {
    casadi_assert_message(H_.getNumInputs()>=3, "Wrong number of input arguments to H");
    casadi_assert_message(H_.getNumOutputs()>=1, "Wrong number of output arguments to H");
    casadi_assert_message(H_.input(0).numel()==n,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size1()==n,"Inconsistent dimensions");
    casadi_assert_message(H_.output().size2()==n,"Inconsistent dimensions");
  }
  if(!J_.isNull() && !G_.isNull()) {
    casadi_assert_message(J_.getNumInputs()>=1, "Wrong number of input arguments to J");
    casadi_assert_message(J_.getNumOutputs()>=1, "Wrong number of output arguments to J");
    casadi_assert_message(J_.input().numel()==n,"Inconsistent dimensions");
    casadi_assert_message(J_.output().size2()==n,"Inconsistent dimensions");
  }
  if(!GF_.isNull()) {
    casadi_assert_message(GF_.getNumInputs()>=1, "Wrong number of input arguments to GF");
    casadi_assert_message(GF_.getNumOutputs()>=1, "Wrong number of output arguments to GF");
    casadi_assert_message(GF_.input().numel()==n,"Inconsistent dimensions");
    casadi_assert_message(GF_.output().size1()==n,"Inconsistent dimensions");
    casadi_assert_message(GF_.output().size2()==n,"Inconsistent dimensions");
  }
  

  
  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
  }
  if(!J_.isNull()) J_.init();
  if(!H_.isNull()) H_.init();
  if(!GF_.isNull()) GF_.init();
  
  // Check if any of the functions have a second argument (i.e. to pass parameters)
  std::vector<FX> functions;
  functions.push_back(F_);
  functions.push_back(G_);
  functions.push_back(H_);
  functions.push_back(J_);
  functions.push_back(GF_);
  
  for (int k=0;k<functions.size();k++) {
    const FX &f = functions[k];
    if (f.isNull()) continue;
    if (f.getNumInputs()!=2) continue;
    pn = f.input(1).size1();
    pm = f.input(1).size2();
    
    if (pn==0 or pm==0)
     continue;
    
    if ((pn!=pn_ || pm!=pm_) && pn_!=0 && pm_!=0) {
      stringstream s;
      s << "One of your supplied functions had a second input argument, which was interpreted as a parameter of shape (" << pn_ << "x" << pm_ << ")." << std::endl;
      s << "However, another function had a second input argument of shape (" << pn << "x" << pm << ")." << std::endl;
      s << "This is inconsistent." << std::endl;
      throw CasadiException(s.str());
    }
    pn_ = pn;
    pm_ = pm;

  }
  
  // Call the init method of the base class
  NLPSolverInternal::init();
  
    // Start the application
  app = new Ipopt::IpoptApplication();

  // Create an Ipopt user class -- need to use Ipopts spart pointer class
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = new Ipopt::SmartPtr<Ipopt::TNLP>();
  userclass = (void*)ucptr;
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;
  uc = new IpoptUserClass(this);
    
  // read options
  exact_hessian_ = !(hasSetOption("hessian_approximation") && getOption("hessian_approximation")=="limited-memory");
  casadi_assert_message(!(exact_hessian_ && H_.isNull()), "No hessian has been provided");
  
  if(verbose_){
    cout << "There are " << n_ << " variables and " << m_ << " constraints." << endl;
    if(exact_hessian_) std::cout << "Using exact Hessian" << std::endl;
    else             std::cout << "Using limited memory Hessian approximation" << std::endl;
  }
 
  // Pass all the options to ipopt
  for(map<string,opt_type>::const_iterator it=ops_.begin(); it!=ops_.end(); ++it)
    if(hasSetOption(it->first)){
      GenericType op = getOption(it->first);
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

void IpoptInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);


  // Set the static parameter
  if (!F_.isNull()) {
    if (F_.getNumInputs()==2) F_.setInput(input(NLP_P),1);
  }
  if (!G_.isNull()) {
    if (G_.getNumInputs()==2) G_.setInput(input(NLP_P),1);
  }
  if (!H_.isNull()) {
    if (H_.getNumInputs()==4) H_.setInput(input(NLP_P),1);
  }
  if (!J_.isNull()) {
    if (J_.getNumInputs()==2) J_.setInput(input(NLP_P),1);
  }
  if (!GF_.isNull()) {
    if (GF_.getNumInputs()==2) GF_.setInput(input(NLP_P),1);
  }

  // Reset the counters
  t_eval_f_ = t_eval_grad_f_ = t_eval_g_ = t_eval_jac_g_ = t_eval_h_ = 0;
  
  // Get back the smart pointer
  Ipopt::SmartPtr<Ipopt::TNLP> *ucptr = (Ipopt::SmartPtr<Ipopt::TNLP>*)userclass;
  Ipopt::SmartPtr<Ipopt::TNLP> &uc = *ucptr;

  // Ask Ipopt to solve the problem
  Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(uc);
  
  if (hasOption("print_time") && bool(getOption("print_time"))) {
    // Write timings
    cout << "time spent in eval_f: " << t_eval_f_ << " s." << endl;
    cout << "time spent in eval_grad_f: " << t_eval_grad_f_ << " s." << endl;
    cout << "time spent in eval_g: " << t_eval_g_ << " s." << endl;
    cout << "time spent in eval_jac_g: " << t_eval_jac_g_ << " s." << endl;
    cout << "time spent in eval_h: " << t_eval_h_ << " s." << endl;
  }

  if (status == Solve_Succeeded)
    stats_["return_status"] = "Solve_Succeeded";
  if (status == Solved_To_Acceptable_Level)
    stats_["return_status"] = "Solved_To_Acceptable_Level";
  if (status == Infeasible_Problem_Detected)
    stats_["return_status"] = "Infeasible_Problem_Detected";
  if (status == Search_Direction_Becomes_Too_Small)
    stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
  if (status == Diverging_Iterates)
    stats_["return_status"] = "Diverging_Iterates";
  if (status == User_Requested_Stop)
    stats_["return_status"] = "User_Requested_Stop";
  if (status == Maximum_Iterations_Exceeded)
    stats_["return_status"] = "Maximum_Iterations_Exceeded";
  if (status == Restoration_Failed)
    stats_["return_status"] = "Restoration_Failed";
  if (status == Error_In_Step_Computation)
    stats_["return_status"] = "Error_In_Step_Computation";
  if (status == Not_Enough_Degrees_Of_Freedom)
    stats_["return_status"] = "Not_Enough_Degrees_Of_Freedom";
  if (status == Invalid_Problem_Definition)
    stats_["return_status"] = "Invalid_Problem_Definition";
  if (status == Invalid_Option)
    stats_["return_status"] = "Invalid_Option";
  if (status == Invalid_Number_Detected)
    stats_["return_status"] = "Invalid_Number_Detected";
  if (status == Unrecoverable_Exception)
    stats_["return_status"] = "Unrecoverable_Exception";
  if (status == NonIpopt_Exception_Thrown)
    stats_["return_status"] = "NonIpopt_Exception_Thrown";
  if (status == Insufficient_Memory)
    stats_["return_status"] = "Insufficient_Memory";

}

void IpoptInternal::finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value){
  try {
    copy(x,x+n_,output(NLP_X_OPT).begin());
    copy(z_L,z_L+n_,output(NLP_LAMBDA_LBX).begin());
    copy(z_U,z_U+n_,output(NLP_LAMBDA_UBX).begin());
    copy(lambda,lambda+m_,output(NLP_LAMBDA_OPT).begin());
    output(NLP_COST).at(0) = obj_value;
  } catch (exception& ex){
    cerr << "finalize_solution failed: " << ex.what() << endl;
  }
}

bool IpoptInternal::eval_h(const double* x, bool new_x, double obj_factor, const double* lambda,bool new_lambda, int nele_hess, int* iRow,int* jCol, double* values){
  try{
    log("eval_h started");
    double time1 = clock();
    if (values == NULL) {
      int nz=0;
      vector<int> rowind,col;
      H_.output().sparsity().getSparsityCRS(rowind,col);
      for(int r=0; r<rowind.size()-1; ++r)
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
         if(col[el]<=r){
            iRow[nz] = r;
            jCol[nz] = col[el];
            nz++;
         }
        }
    } else {
      // Pass input
      H_.setInput(x);
      H_.setInput(lambda,H_.getNumInputs()==4? 2 : 1);
      H_.setInput(obj_factor,H_.getNumInputs()==4? 3 : 2);

      // Evaluate
      H_.evaluate();

      // Get results
      H_.output().get(values,SPARSESYM);
    }
    double time2 = clock();
    t_eval_h_ += double(time2-time1)/CLOCKS_PER_SEC;
    log("eval_h ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_h failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_jac_g(int n, const double* x, bool new_x,int m, int nele_jac, int* iRow, int *jCol,double* values){
  try{
    log("eval_jac_g started");
    double time1 = clock();
    if (values == NULL) {
      int nz=0;
      vector<int> rowind,col;
      J_.output().sparsity().getSparsityCRS(rowind,col);
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
        J_.output().printSparse();
      }
    }
    
    double time2 = clock();
    t_eval_jac_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
    log("eval_jac_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_jac_g failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_f(int n, const double* x, bool new_x, double& obj_value)
{
  try {
    log("eval_f started");
    
    // Log time
    double time1 = clock();
    casadi_assert(n == n_);

    // Pass the argument to the function
    F_.setInput(x);
      
    // Evaluate the function
    F_.evaluate();

    // Get the result
    F_.getOutput(obj_value);

    // Printing
    if(monitored("eval_x"))
      cout << "x = " << F_.input() << endl;
      
    // Printing
    if(monitored("eval_f")){
      cout << "obj_value = " << obj_value << endl;
    }

    double time2 = clock();
    t_eval_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    log("eval_f ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_f failed: " << ex.what() << endl;
    return false;
  }
    
}

bool IpoptInternal::eval_g(int n, const double* x, bool new_x, int m, double* g)
{
  try {
    log("eval_g started");
    double time1 = clock();

    casadi_assert(n == n_);
    casadi_assert(m == m_);

    // Pass the argument to the function
    G_.setInput(x);

    // Evaluate the function and tape
    G_.evaluate();

    // Ge the result
    G_.getOutput(g);

    if(monitored("eval_x"))
      cout << "x = " << G_.input() << endl;
    // Printing
    if(monitored("eval_g"))
      cout << "g = " << G_.output() << endl;
      
    double time2 = clock();
    t_eval_g_ += double(time2-time1)/CLOCKS_PER_SEC;
    
    log("eval_g ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_g failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::eval_grad_f(int n, const double* x, bool new_x, double* grad_f)
{
  try {
    log("eval_grad_f started");
    double time1 = clock();
    casadi_assert(n == n_);
    
    // If no gradient function has been provided, use AD adjoint
    if(GF_.isNull()){
    
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
      
    } else {
      
      // Pass the argument to the function
      GF_.setInput(x);
      
      // Evaluate, adjoint mode
      GF_.evaluate();

      // Get the result
      GF_.getOutput(grad_f);
      
      // Printing
      if(monitored("eval_grad_f")){
        cout << "grad_f = " << GF_.output() << endl;
      }
    }
    
    double time2 = clock();
    t_eval_grad_f_ += double(time2-time1)/CLOCKS_PER_SEC;

    // Check the result for regularity
    for(int i=0; i<n; ++i){
        if(isnan(grad_f[i]) || isinf(grad_f[i])){
          log("eval_grad_f: result not regular");
          return false;
      }
    }

    log("eval_grad_f ok");
    return true;
  } catch (exception& ex){
    cerr << "eval_jac_f failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::get_bounds_info(int n, double* x_l, double* x_u,
                                int m, double* g_l, double* g_u)
{
  try {
    casadi_assert(n == n_);
    casadi_assert(m == m_);
    input(NLP_LBX).getArray(x_l,n);
    input(NLP_UBX).getArray(x_u,n);
    input(NLP_LBG).getArray(g_l,m);
    input(NLP_UBG).getArray(g_u,m);
    return true;
  } catch (exception& ex){
    cerr << "get_bounds_info failed: " << ex.what() << endl;
    return false;
  }
}

bool IpoptInternal::get_starting_point(int n, bool init_x, double* x,
                                   bool init_z, double* z_L, double* z_U,
                                   int m, bool init_lambda,
                                   double* lambda)
{
  try {
    bool warmstart = hasSetOption("warm_start_init_point") && getOption("warm_start_init_point")=="yes";
    if (warmstart) {
      input(NLP_X_INIT).getArray(x,n);
      if (init_z) {
        output(NLP_LAMBDA_LBX).getArray(z_L,n);
        output(NLP_LAMBDA_UBX).getArray(z_U,n);
      }
      if (init_lambda)
        input(NLP_LAMBDA_INIT).getArray(lambda,m);
      return true;
    } else {
      casadi_assert(init_x == true);
      casadi_assert(init_z == false);
      casadi_assert(init_lambda == false);
      input(NLP_X_INIT).getArray(x,n);
      return true;
    }
  } catch (exception& ex){
    cerr << "get_starting_point failed: " << ex.what() << endl;
    return false;
  }
}

void IpoptInternal::get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag)
{
  try {
    n = n_;               // number of variables
    m = m_;               // number of constraints

    // Get Jacobian sparsity pattern
    if(G_.isNull())
      nnz_jac_g = 0;
    else
      nnz_jac_g = J_.output().size();

    // Get Hessian sparsity pattern
    if(exact_hessian_)
      nnz_h_lag = H_.output().sparsity().sizeL();
    else
      nnz_h_lag = 0;
  } catch (exception& ex){
    cerr << "get_nlp_info failed: " << ex.what() << endl;
  }
}

int IpoptInternal::get_number_of_nonlinear_variables() const{
  try {
    if(H_.isNull() || !bool(getOption("pass_nonlinear_variables"))){
      // No Hessian has been interfaced
      return -1;
    } else {
      // Number of variables that appear nonlinearily
      int nv = 0;
      
      // Loop over the rows
      for(int i=0; i<H_.output().size1(); ++i){
        // If the row contains any non-zeros, the corresponding variable appears nonlinearily
        if(H_.output().rowind(i)!=H_.output().rowind(i+1))
          nv++;
      }
      
      // Return the number
      return nv;
    }
  } catch (exception& ex){
    cerr << "get_number_of_nonlinear_variables failed: " << ex.what() << endl;
    return -1;
  }
}

bool IpoptInternal::get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const{
  try {
    // Running index
    int el = 0;
    
    // Loop over the rows
    for(int i=0; i<H_.output().size1(); ++i){
      // If the row contains any non-zeros, the corresponding variable appears nonlinearily
      if(H_.output().rowind(i)!=H_.output().rowind(i+1)){
        pos_nonlin_vars[el++] = i;
      }
    }
    
    // Assert number and return
    casadi_assert(el==num_nonlin_vars);
    return true;
  } catch (exception& ex){
    cerr << "get_list_of_nonlinear_variables failed: " << ex.what() << endl;
    return false;
  }
}

} // namespace CasADi
