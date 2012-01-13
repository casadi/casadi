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

#include "ipopt_qp_internal.hpp"

#include "casadi/mx/mx_tools.hpp"
#include "casadi/fx/mx_function.hpp"

using namespace std;
namespace CasADi {
namespace Interfaces {

IpoptQPInternal* IpoptQPInternal::clone() const{
  // Return a deep copy
  IpoptQPInternal* node = new IpoptQPInternal(input(QP_H).sparsity(),input(QP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
IpoptQPInternal::IpoptQPInternal(const CRSSparsity& H_, const CRSSparsity &A_) : QPSolverInternal(H_,A_){
  std::cout << "Warning: IPOPT QP is highly experimental" << std::endl;

  std::map<std::string,opt_type> ops_;
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
  
}

IpoptQPInternal::~IpoptQPInternal(){ 
}

void IpoptQPInternal::evaluate(int nfdir, int nadir) {
  if (nfdir!=0 || nadir!=0) throw CasadiException("IpoptQPInternal::evaluate() not implemented for forward or backward mode");
  

  // Pass inputs of QPInternal to IpoptInternal form 
#if 0
  solver.input(NLP_P)[H_.mapping()] = input(QP_H);
  solver.input(NLP_P)[G_.mapping()] = input(QP_G);
  solver.input(NLP_P)[A_.mapping()] = input(QP_A);
#endif
  casadi_assert(0);
  
  solver.input(NLP_LBX).set(input(QP_LBX));
  solver.input(NLP_UBX).set(input(QP_UBX));
  
  solver.input(NLP_LBG).set(input(QP_LBA));
  solver.input(NLP_UBG).set(input(QP_UBA));
  
  // Delegate computation to Ipopt
  solver.evaluate();
  
  // Read the outputs from Ipopt
  output(NLP_X_OPT).set(solver.output(NLP_X_OPT));
  output(NLP_COST).set(solver.output(NLP_COST));
}

void IpoptQPInternal::init(){
  std::cout << "okay" << std::endl;
  QPSolverInternal::init();

  // Create an MX for the decision variables
  MX X("X",nx_,1);
    
  // Put H, G, A sparsities in a vector...
  std::vector< CRSSparsity > sps;
  sps.push_back(input(QP_H).sparsity());
  sps.push_back(input(QP_G).sparsity());
  sps.push_back(input(QP_A).sparsity());
  
  // So that we can pass it on to createParent
  std::pair< MX, std::vector< MX > > mypair = createParent(sps);
  
  // V groups all parameters in an MX
  MX V(mypair.first);
  std::vector< MX > variables(mypair.second);
  
  // H_, G_, A_ depend on V
  H_=variables[0];
  G_=variables[1];
  A_=variables[2];
  std::cout << "okay" << std::endl;
  // We're going to use two-argument objective and constraints to allow the use of parameters
  std::vector< MX > args;
  args.push_back(X);
  args.push_back(V);
    std::cout << "okay there" << std::endl;
  // The objective function looks exactly like a mathematical description of the NLP
  MXFunction QP_f(args, mul(trans(G_),X) + 0.5*mul(mul(trans(X),H_),X));
  QP_f.init();

  // So does the constraint function
  MXFunction QP_g(args, mul(A_,X));
  std::cout << "okay here" << std::endl;
  // Jacobian of the constraints
  MXFunction QP_j(args,A_);
  
  std:cout << (G_+mul(trans(H_),X)).dimString() << std::endl;
  // Gradient of the objective
  MXFunction QP_gf(args,G_+mul(H_,X));
  
  std::cout << "okay everyzhere" << std::endl;
  MX sigma("sigma");
    std::cout << "okay 12" << nc_ << std::endl;
  MX lambda("lambda",nc_,1);
    std::cout << "okay 4" << std::endl;
  args.push_back(lambda);
  args.push_back(sigma);
    std::cout << "okay 5" << std::endl;
  // Hessian of the Lagrangian
  MXFunction QP_h(args,H_*sigma);
    std::cout << "okay" << std::endl;
  // Generate an IpoptSolver that uses this objective and constraint
  solver = IpoptSolver(QP_f,QP_g,QP_h,QP_j,QP_gf);
  
  if (getOption("convex").toBool()) {
    setOption("mehrotra_algorithm","yes");
    setOption("mu_oracle","probing");
  }
  
  // Pass all options to the solver
  for(map<string,GenericType>::const_iterator it=dictionary().begin(); it!=dictionary().end(); ++it) {
    if (solver.hasOption(it->first)) {
     if (hasSetOption(it->first)) solver.setOption(it->first,it->second);
    }
  }
  std::cout << "okay" << std::endl;
  solver.setOption("jac_c_constant","yes");
  solver.setOption("jac_d_constant","yes");
  solver.setOption("hessian_constant","yes");
  
  solver.init();

}

} // namespace Interfaces
} // namespace CasADi

