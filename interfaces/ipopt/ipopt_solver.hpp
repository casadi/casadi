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

#ifndef IPOPT_SOLVER_HPP
#define IPOPT_SOLVER_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class IpoptInternal;
  
// List from ipopt_internal.cpp
/**
* \brief interface to IPOPT NLP solver
*
*
* available options:
\verbatim
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
\endverbatim

Input arguments:  CasADi::NLPInput\n
  
Output arguments: CasADi::NLPOutput\n
*/
class IpoptSolver : public NLPSolver {
  public:
    /// Default constructor
    IpoptSolver();

    /// \brief Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit IpoptSolver(const FX& F,         /**< F objective function: \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}]\f$*/
                         const FX& G = FX(),  /**< constraint function (default only bound constraints): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^m]\f$ */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory): \f$ [\mathbf{R}^n, \mathbf{R}^m, \mathbf{R}] \mapsto [\mathbf{R}^{n x n}]\f$ \n The third input argument for H is \f$ \sigma \f$, a scaling factor for F. */
                         const FX& J = FX(),  /**< Jacobian of G (default -> differentiate): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^{m x n}]\f$ */
                         const FX& GF = FX()  /**< Gradient of the objective function (default: adjoint mode AD on F): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^{n x n}]\f$ */
                        );

    /// Access functions of the node
    IpoptInternal* operator->();
    const IpoptInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    
};

} // namespace CasADi

#endif //IPOPT_SOLVER_HPP
