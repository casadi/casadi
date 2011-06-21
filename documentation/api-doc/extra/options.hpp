/** \class CasADi::Interfaces::LapackLUDenseInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LapackLUDense
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::KinsolInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>"Stopping criterion tolerance"</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>"Globalization strategy (\"none\" or \"linesearch\")"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::KinsolSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>"Stopping criterion tolerance"</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>"Globalization strategy (\"none\" or \"linesearch\")"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShootingInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShooting
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunctionInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td></td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symbolic_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>generate jacobian symbolically by source code transformation</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunction
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td></td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symbolic_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>generate jacobian symbolically by source code transformation</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IpoptInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
</table>
*/
/** \class CasADi::IpoptSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesIntegrator
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::JacobianInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>finite_differences</td><td>OT_BOOLEAN</td><td>false</td><td>Using finite differences instead of automatic differentiation</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Jacobian
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default", i.e. forward if n_<=m_, otherwise adjoint</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>finite_differences</td><td>OT_BOOLEAN</td><td>false</td><td>Using finite differences instead of automatic differentiation</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CplexInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td>Option()</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_cplex_problem"</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CplexSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td>Option()</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_cplex_problem"</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_time</td><td>OT_REAL</td><td>0.0</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoInterface
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_time</td><td>OT_REAL</td><td>0.0</td><td></td><td>CasADi::AcadoInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>lifted</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>optimizer</td><td>OT_STRING</td><td>"sqp"</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>lifted</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>optimizer</td><td>OT_STRING</td><td>"sqp"</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IntegratorInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Integrator
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLUInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>conditionnumber</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>diagpivotthresh</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>iterrefine</td><td>OT_STRING</td><td>"norefine"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pivotgrowth</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>printstat</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>rowperm</td><td>OT_STRING</td><td>"largediag"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symmetricmode</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLU
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>conditionnumber</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>diagpivotthresh</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>iterrefine</td><td>OT_STRING</td><td>"norefine"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pivotgrowth</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>printstat</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>rowperm</td><td>OT_STRING</td><td>"largediag"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symmetricmode</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::FXInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::FX
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>"use IDACalcIC to get consistent initial conditions"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>false</td><td>"use IDACalcIC to get consistent initial conditions"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>"IDAS scaling on cj for the user-defined linear solver module"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>"Call calc ic an extra time, with fsens=0"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>"first requested time as a fraction of the time interval"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>"maximim step size"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>"supress algebraic variables in the error testing"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasIntegrator
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>"use IDACalcIC to get consistent initial conditions"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>false</td><td>"use IDACalcIC to get consistent initial conditions"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>"IDAS scaling on cj for the user-defined linear solver module"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>"Call calc ic an extra time, with fsens=0"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>"first requested time as a fraction of the time interval"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>"maximim step size"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>"supress algebraic variables in the error testing"</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ParallelizerInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>save_corrected_input</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Parallelizer
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>save_corrected_input</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OCPSolverInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OCPSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ImplicitFunctionInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>"Stopping criterion tolerance"</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ImplicitFunction
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>"Stopping criterion tolerance"</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LinearSolverInternal
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LinearSolver
List of available options
<table>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
