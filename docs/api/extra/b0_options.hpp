/// \cond INTERNAL
/** \class casadi::RKIntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::RKIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LapackLUDenseInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>casadi::LapackLUDenseInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::LapackLUDenseInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::LapackLUDense
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>casadi::LapackLUDenseInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::LapackLUDenseInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CVodesInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>casadi::CVodesInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>casadi::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>casadi::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>casadi::FunctionInternal<br />casadi::CVodesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>casadi::CVodesInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CVodesIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>casadi::CVodesInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>casadi::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>casadi::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>casadi::FunctionInternal<br />casadi::CVodesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>casadi::CVodesInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::OOQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>casadi::OOQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>casadi::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>casadi::OOQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::OOQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>casadi::OOQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>casadi::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>casadi::OOQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QPLPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QPSOlver used to solve the LPs.</td><td>casadi::QPLPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QPSOlver</td><td>casadi::QPLPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QPLPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QPSOlver used to solve the LPs.</td><td>casadi::QPLPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QPSOlver</td><td>casadi::QPLPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::PsdIndefDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::PsdIndefDpleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::IntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Integrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LapackQRDenseInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::LapackQRDense
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ExternalFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::ExternalFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedQPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::StabilizedQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ControlSimulatorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>control_endpoint</td><td>OT_BOOLEAN</td><td>false</td><td>Include a control value at the end of the simulation domain. Used for interpolation.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>control_interpolation</td><td>OT_STRING</td><td>"none"</td><td>none|nearest|linear</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>minor_grid</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>The local grid used on each major interval, with time normalized to 1. By default, option 'nf' is used to construct a linearly spaced grid.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nf</td><td>OT_INTEGER</td><td>1</td><td>Number of minor grained integration steps per major interval. nf&gt;0 must hold. This option is not used when 'minor_grid' is provided.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>simulator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the simulator</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::ControlSimulator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>control_endpoint</td><td>OT_BOOLEAN</td><td>false</td><td>Include a control value at the end of the simulation domain. Used for interpolation.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>control_interpolation</td><td>OT_STRING</td><td>"none"</td><td>none|nearest|linear</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>minor_grid</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>The local grid used on each major interval, with time normalized to 1. By default, option 'nf' is used to construct a linearly spaced grid.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nf</td><td>OT_INTEGER</td><td>1</td><td>Number of minor grained integration steps per major interval. nf&gt;0 must hold. This option is not used when 'minor_grid' is provided.</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>simulator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the simulator</td><td>casadi::ControlSimulatorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ParallelizerInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>(serial|openmp|mpi)</td><td>casadi::ParallelizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Parallelizer
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>(serial|openmp|mpi)</td><td>casadi::ParallelizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SOCPQCQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>socp_solver</td><td>OT_SOCPSOLVER</td><td>GenericType()</td><td>The SOCPSolver used to solve the QCQPs.</td><td>casadi::SOCPQCQPInternal</td></tr>
<tr><td>socp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SOCPSOlver</td><td>casadi::SOCPQCQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SOCPQCQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>socp_solver</td><td>OT_SOCPSOLVER</td><td>GenericType()</td><td>The SOCPSolver used to solve the QCQPs.</td><td>casadi::SOCPQCQPInternal</td></tr>
<tr><td>socp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SOCPSOlver</td><td>casadi::SOCPQCQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimulatorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(initial|step)</td><td>casadi::FunctionInternal<br />casadi::SimulatorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Simulator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(initial|step)</td><td>casadi::FunctionInternal<br />casadi::SimulatorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SOCPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SOCPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::OptionsFunctionalityNode
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \endcond
/** \class casadi::OptionsFunctionality
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::XFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimpleHomotopyNLPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::HomotopyNLPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLP solver to be used by the Homotopy solver</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Homotopy solver</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>num_steps</td><td>OT_INTEGER</td><td>10</td><td>Take this many steps to go from tau=0 to tau=1.</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SimpleHomotopyNLPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::HomotopyNLPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLP solver to be used by the Homotopy solver</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Homotopy solver</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>num_steps</td><td>OT_INTEGER</td><td>10</td><td>Take this many steps to go from tau=0 to tau=1.</td><td>casadi::SimpleHomotopyNLPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SDQPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDQPSolver used to solve the SDPs.</td><td>casadi::SDQPSolverInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDQPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SDQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDQPSolver used to solve the SDPs.</td><td>casadi::SDQPSolverInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDQPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedSQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>TReta1</td><td>OT_REAL</td><td>0.8</td><td>Required predicted / actual decrease for TR increase</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>TReta2</td><td>OT_REAL</td><td>0.2</td><td>Required predicted / actual decrease for TR decrease</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>alphaMin</td><td>OT_REAL</td><td>1e-3</td><td>Used to check whether to increase rho.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.5</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>0.001</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dvMax0</td><td>OT_REAL</td><td>100</td><td>Parameter used to defined the max step length.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>eps_active</td><td>OT_REAL</td><td>1e-6</td><td>Threshold for the epsilon-active set.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gamma1</td><td>OT_REAL</td><td>2.</td><td>Trust region increase parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gamma2</td><td>OT_REAL</td><td>1.</td><td>Trust region update parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gamma3</td><td>OT_REAL</td><td>1.</td><td>Trust region decrease parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>20</td><td>Maximum number of linesearch iterations</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_time</td><td>OT_REAL</td><td>1e12</td><td>Timeout</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::StabilizedSQPInternal</td></tr>
<tr><td>muR0</td><td>OT_REAL</td><td>1e-4</td><td>Initial choice of regularization parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nu</td><td>OT_REAL</td><td>1</td><td>Parameter for primal-dual augmented Lagrangian.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>phiWeight</td><td>OT_REAL</td><td>1e-5</td><td>Weight used in pseudo-filter.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>stabilized_qp_solver</td><td>OT_STABILIZEDQPSOLVER</td><td>GenericType()</td><td>The Stabilized QP solver to be used by the SQP method</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>stabilized_qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Stabilized QP solver</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tau0</td><td>OT_REAL</td><td>1e-2</td><td>Initial parameter for the merit function optimality threshold.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for dual infeasability</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for primal infeasibility</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>yEinitial</td><td>OT_STRING</td><td>"simple"</td><td>Initial multiplier. Simple (all zero) or least (LSQ).</td><td>casadi::StabilizedSQPInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::StabilizedSQPMethod
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>TReta1</td><td>OT_REAL</td><td>0.8</td><td>Required predicted / actual decrease for TR increase</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>TReta2</td><td>OT_REAL</td><td>0.2</td><td>Required predicted / actual decrease for TR decrease</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>alphaMin</td><td>OT_REAL</td><td>1e-3</td><td>Used to check whether to increase rho.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.5</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>0.001</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dvMax0</td><td>OT_REAL</td><td>100</td><td>Parameter used to defined the max step length.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>eps_active</td><td>OT_REAL</td><td>1e-6</td><td>Threshold for the epsilon-active set.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gamma1</td><td>OT_REAL</td><td>2.</td><td>Trust region increase parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gamma2</td><td>OT_REAL</td><td>1.</td><td>Trust region update parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gamma3</td><td>OT_REAL</td><td>1.</td><td>Trust region decrease parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>20</td><td>Maximum number of linesearch iterations</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>max_time</td><td>OT_REAL</td><td>1e12</td><td>Timeout</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::StabilizedSQPInternal</td></tr>
<tr><td>muR0</td><td>OT_REAL</td><td>1e-4</td><td>Initial choice of regularization parameter</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nu</td><td>OT_REAL</td><td>1</td><td>Parameter for primal-dual augmented Lagrangian.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>phiWeight</td><td>OT_REAL</td><td>1e-5</td><td>Weight used in pseudo-filter.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>stabilized_qp_solver</td><td>OT_STABILIZEDQPSOLVER</td><td>GenericType()</td><td>The Stabilized QP solver to be used by the SQP method</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>stabilized_qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Stabilized QP solver</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tau0</td><td>OT_REAL</td><td>1e-2</td><td>Initial parameter for the merit function optimality threshold.</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for dual infeasability</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for primal infeasibility</td><td>casadi::StabilizedSQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>yEinitial</td><td>OT_STRING</td><td>"simple"</td><td>Initial multiplier. Simple (all zero) or least (LSQ).</td><td>casadi::StabilizedSQPInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SymbolicQRInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SymbolicQRInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SymbolicQRInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SymbolicQR
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SymbolicQRInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SymbolicQRInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QPStabilizerInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver used to solve the stabilized QPs.</td><td>casadi::QPStabilizerInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver instance</td><td>casadi::QPStabilizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QPStabilizer
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver used to solve the stabilized QPs.</td><td>casadi::QPStabilizerInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver instance</td><td>casadi::QPStabilizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::KnitroInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>0.0001</td><td>Feasible mode tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>1</td><td>Absolute feasible tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>0</td><td>Allow infeasible iterations</td><td>casadi::KnitroInternal</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td><td>casadi::KnitroInternal</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td><td>casadi::KnitroInternal</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td><td>casadi::KnitroInternal</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td><td>casadi::KnitroInternal</td></tr>
<tr><td>LPSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LPSolver</td><td>casadi::KnitroInternal</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td><td>casadi::KnitroInternal</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td><td>casadi::KnitroInternal</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td><td>casadi::KnitroInternal</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e-8</td><td>Maximum objective value</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td><td>casadi::KnitroInternal</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>casadi::KnitroInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::KnitroInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::KnitroSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>0.0001</td><td>Feasible mode tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>1</td><td>Absolute feasible tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>0</td><td>Allow infeasible iterations</td><td>casadi::KnitroInternal</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td><td>casadi::KnitroInternal</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td><td>casadi::KnitroInternal</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td><td>casadi::KnitroInternal</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td><td>casadi::KnitroInternal</td></tr>
<tr><td>LPSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LPSolver</td><td>casadi::KnitroInternal</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td><td>casadi::KnitroInternal</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td><td>casadi::KnitroInternal</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td><td>casadi::KnitroInternal</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e-8</td><td>Maximum objective value</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td><td>casadi::KnitroInternal</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td><td>casadi::KnitroInternal</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td><td>casadi::KnitroInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>casadi::KnitroInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::KnitroInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SnoptInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_iprint</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::SnoptInternal</td></tr>
<tr><td>_isumm</td><td>OT_INTEGER</td><td>6</td><td></td><td>casadi::SnoptInternal</td></tr>
<tr><td>_start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Warm)</td><td>casadi::SnoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td><td>casadi::SnoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_nlp|setup_nlp)</td><td>casadi::FunctionInternal<br />casadi::SnoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>casadi::SnoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SnoptSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_iprint</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::SnoptInternal</td></tr>
<tr><td>_isumm</td><td>OT_INTEGER</td><td>6</td><td></td><td>casadi::SnoptInternal</td></tr>
<tr><td>_start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Warm)</td><td>casadi::SnoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td><td>casadi::SnoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_nlp|setup_nlp)</td><td>casadi::FunctionInternal<br />casadi::SnoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>casadi::SnoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DSDPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_loglevel</td><td>OT_INTEGER</td><td>0</td><td>An integer that specifies how much logging is done on stdout.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_penalty</td><td>OT_REAL</td><td>1e5</td><td>Penality parameter lambda. Must exceed the trace of Y. This parameter heavily influences the ability of DSDP to treat linear equalities. The DSDP standard default (1e8) will make a problem with linear equality return unusable solutions.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_printlevel</td><td>OT_INTEGER</td><td>1</td><td>A printlevel of zero will disable all output. Another number indicates how often a line is printed.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_reuse</td><td>OT_INTEGER</td><td>4</td><td>Maximum on the number of times the Schur complement matrix is reused</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_rho</td><td>OT_REAL</td><td>4.0</td><td>Potential parameter. Must be &gt;=1</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_use_penalty</td><td>OT_BOOLEAN</td><td>true</td><td>Modifies the algorithm to use a penality gamma on r.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_zbar</td><td>OT_REAL</td><td>1e10</td><td>Initial upper bound on the objective of the dual problem.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td><td>casadi::DSDPInternal</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td><td>casadi::DSDPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>infinity</td><td>OT_REAL</td><td>1e30</td><td>Treat numbers higher than this as infinity</td><td>casadi::DSDPInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td><td>casadi::DSDPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td><td>casadi::DSDPInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::DSDPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_loglevel</td><td>OT_INTEGER</td><td>0</td><td>An integer that specifies how much logging is done on stdout.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_penalty</td><td>OT_REAL</td><td>1e5</td><td>Penality parameter lambda. Must exceed the trace of Y. This parameter heavily influences the ability of DSDP to treat linear equalities. The DSDP standard default (1e8) will make a problem with linear equality return unusable solutions.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_printlevel</td><td>OT_INTEGER</td><td>1</td><td>A printlevel of zero will disable all output. Another number indicates how often a line is printed.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_reuse</td><td>OT_INTEGER</td><td>4</td><td>Maximum on the number of times the Schur complement matrix is reused</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_rho</td><td>OT_REAL</td><td>4.0</td><td>Potential parameter. Must be &gt;=1</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_use_penalty</td><td>OT_BOOLEAN</td><td>true</td><td>Modifies the algorithm to use a penality gamma on r.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>_zbar</td><td>OT_REAL</td><td>1e10</td><td>Initial upper bound on the objective of the dual problem.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td><td>casadi::DSDPInternal</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td><td>casadi::DSDPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>infinity</td><td>OT_REAL</td><td>1e30</td><td>Treat numbers higher than this as infinity</td><td>casadi::DSDPInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td><td>casadi::DSDPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td><td>casadi::DSDPInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td><td>casadi::DSDPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DirectSingleShootingInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to casadi::Parallelizer</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::DirectSingleShooting
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to casadi::Parallelizer</td><td>casadi::DirectSingleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DirectMultipleShootingInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to casadi::Parallelizer</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::DirectMultipleShooting
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to casadi::Parallelizer</td><td>casadi::DirectMultipleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SDPSDQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDPSolver used to solve the SDQPs.</td><td>casadi::SDPSDQPInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDPSDQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SDPSDQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDPSolver used to solve the SDQPs.</td><td>casadi::SDPSDQPInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDPSDQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ImplicitFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::ImplicitFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::IdasInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>casadi::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>casadi::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>casadi::IdasInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>casadi::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>casadi::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td><td>casadi::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>casadi::IdasInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td><td>casadi::IdasInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup)</td><td>casadi::FunctionInternal<br />casadi::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Supress algebraic variables in the error testing</td><td>casadi::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::IdasIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>casadi::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>casadi::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>casadi::IdasInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>casadi::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>casadi::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td><td>casadi::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>casadi::IdasInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td><td>casadi::IdasInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup)</td><td>casadi::FunctionInternal<br />casadi::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Supress algebraic variables in the error testing</td><td>casadi::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::FixedStepIntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::FixedStepIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QCQPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QCQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::OCPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::OCPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::HomotopyNLPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::HomotopyNLPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::HomotopyNLPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::HomotopyNLPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SXFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>casadi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_sparsity</td><td>OT_BOOLEAN</td><td>false</td><td>Propagate sparsity patterns using just-in-time compilation to a CPU or GPU using OpenCL</td><td>casadi::SXFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SXFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>casadi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_sparsity</td><td>OT_BOOLEAN</td><td>false</td><td>Propagate sparsity patterns using just-in-time compilation to a CPU or GPU using OpenCL</td><td>casadi::SXFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NewtonImplicitInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>casadi::FunctionInternal<br />casadi::NewtonImplicitInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::NewtonImplicitSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>casadi::NewtonImplicitInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>casadi::FunctionInternal<br />casadi::NewtonImplicitInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedSQICInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::StabilizedSQICSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimpleIndefDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::SimpleIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::SimpleIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SimpleIndefDpleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::SimpleIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::SimpleIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::KinsolInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>casadi::KinsolInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>casadi::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>casadi::KinsolInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>casadi::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.</td><td>casadi::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>casadi::FunctionInternal<br />casadi::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::KinsolInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (none|linesearch)</td><td>casadi::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>casadi::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::KinsolSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>casadi::KinsolInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>casadi::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>casadi::KinsolInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>casadi::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.</td><td>casadi::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>casadi::FunctionInternal<br />casadi::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::KinsolInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (none|linesearch)</td><td>casadi::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>casadi::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CollocationIntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::CollocationIntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::CollocationIntegratorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CollocationIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::CollocationIntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::CollocationIntegratorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SQICInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SQICSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SundialsInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SundialsIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::DpleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i,i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::WorhpInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>AcceptTolFeas</td><td>OT_REAL</td><td>0.001</td><td>Tolerance for acceptable feasibility</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AcceptTolOpti</td><td>OT_REAL</td><td>0.001</td><td>Tolerance for acceptable optimality</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AlphaMinConst</td><td>OT_BOOLEAN</td><td>False</td><td>Use a constant lower bound on Armijo stepsize in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Ares</td><td>OT_INTEGERVECTOR</td><td>[42, 41, 42, 43, 44, 41, 50]</td><td>Armijo recovery strategies. Vector of size 7</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoBeta</td><td>OT_REAL</td><td>0.712</td><td>Trial stepsize decrease factor for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMaxAlpha</td><td>OT_REAL</td><td>1.0</td><td>Initial alpha for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlpha</td><td>OT_REAL</td><td>1e-06</td><td>Lower bound on alpha for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlphaRec</td><td>OT_REAL</td><td>1e-06</td><td>Lower bound on alpha for Armijo rule during recovery</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoSigma</td><td>OT_REAL</td><td>0.005</td><td>Scale factor for linearised descent check in Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AutoQPRecovery</td><td>OT_BOOLEAN</td><td>True</td><td>Enable automatic QP recovery</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSmaxblockSize</td><td>OT_INTEGER</td><td>300</td><td>Block size parameter used by certain BFGS methods</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSmethod</td><td>OT_INTEGER</td><td>0</td><td>Choose BFGS method (0: dense, 1-3: block, 100+: sparse)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSminblockSize</td><td>OT_INTEGER</td><td>300</td><td>Block size parameter used by certain BFGS methods</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSrestart</td><td>OT_INTEGER</td><td>50</td><td>Restart BFGS update after this many iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BettsFactor</td><td>OT_REAL</td><td>2.1</td><td>Update factor for Betts' Hessian regularisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BettsPoint</td><td>OT_REAL</td><td>1.0</td><td>Smallest eigenvalue of the regularised Hessian</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BoundTolFac</td><td>OT_REAL</td><td>1000.0</td><td>Factor in determining active constraints by KKT</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckFJ</td><td>OT_REAL</td><td>1e+12</td><td>Upper bound used by Fritz-John heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureDF</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of DF</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureDG</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of DG</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureHM</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of HM</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepBettsSum</td><td>OT_REAL</td><td>0.5</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepConStop</td><td>OT_REAL</td><td>1e-06</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepConvio</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMaxIter</td><td>OT_INTEGER</td><td>50</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMethod</td><td>OT_INTEGER</td><td>0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMode</td><td>OT_INTEGER</td><td>1</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepPFactor</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepPMax</td><td>OT_REAL</td><td>1000000.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepRecoveryDX</td><td>OT_BOOLEAN</td><td>False</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvBCond</td><td>OT_REAL</td><td>0.02</td><td>Block BFGS curvature condition bound</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvBFac</td><td>OT_REAL</td><td>0.3</td><td>Block BFGS curvature condition regularisation factor</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvCond</td><td>OT_REAL</td><td>0.02</td><td>BFGS Curvature condition bound</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvFac</td><td>OT_REAL</td><td>0.3</td><td>BFGS curvature condition regularisation factor</td><td>casadi::WorhpInternal</td></tr>
<tr><td>DebugMarker05</td><td>OT_INTEGER</td><td>42</td><td>Debug marker. Used to find memory alignment/padding issues</td><td>casadi::WorhpInternal</td></tr>
<tr><td>DebugMarker06</td><td>OT_INTEGER</td><td>42</td><td>Debug marker. Used to find memory alignment/padding issues</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FGtogether</td><td>OT_BOOLEAN</td><td>False</td><td>F and G cannot be evaluated separately</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FJandND</td><td>OT_BOOLEAN</td><td>False</td><td>Enable Fritz-John and non-differentiable check heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleDual</td><td>OT_BOOLEAN</td><td>False</td><td>Activate dual feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleInit</td><td>OT_BOOLEAN</td><td>False</td><td>Activate initial feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleInitTol</td><td>OT_REAL</td><td>0.001</td><td>Feasibility tolerance for no-objective feasible mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleOnly</td><td>OT_BOOLEAN</td><td>False</td><td>Activate feasible-only mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FidifEps</td><td>OT_REAL</td><td>1e-05</td><td>Finite difference perturbation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FidifHM</td><td>OT_BOOLEAN</td><td>False</td><td>Approximate Hessian by finite differences (otherwise BFGS)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterBisecAlpha</td><td>OT_BOOLEAN</td><td>True</td><td>Filter heuristic to save Armijo iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterGammaCV</td><td>OT_REAL</td><td>7.5e-06</td><td>Constraint violation decrease factor in Filter acceptance check</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterGammaF</td><td>OT_REAL</td><td>1.1e-05</td><td>Objective decrease factor in Filter acceptance check</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterIntersecAlpha</td><td>OT_BOOLEAN</td><td>True</td><td>Filter heuristic to save Armijo iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FirstDifCentral</td><td>OT_BOOLEAN</td><td>True</td><td>Use central finite difference quotient for first derivatives</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FocusOnFeas</td><td>OT_BOOLEAN</td><td>True</td><td>Enable Focus-on-Feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FocusOnFeasFactor</td><td>OT_REAL</td><td>1.36</td><td>Factor in Focus-on-Feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>GammaAlpha</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for alphamin calculation by Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>GroupMethod</td><td>OT_INTEGER</td><td>1</td><td>Select method to determine graph colouring groups</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IgnoreFilterCrit</td><td>OT_BOOLEAN</td><td>False</td><td>Activate accelerating heuristics for Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncBettsTau</td><td>OT_REAL</td><td>2.0</td><td>Increase factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncBettsTauMore</td><td>OT_REAL</td><td>100.0</td><td>Larger increase factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncreaseIWS</td><td>OT_REAL</td><td>1.0</td><td>Increase factor for estimated integer workspace requirement</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncreaseRWS</td><td>OT_REAL</td><td>1.0</td><td>Increase factor for estimated real workspace requirement</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Infty</td><td>OT_REAL</td><td>1e+20</td><td>Upper bound for numbers to be regarded as finite</td><td>casadi::WorhpInternal</td></tr>
<tr><td>InftyUnbounded</td><td>OT_REAL</td><td>1e+20</td><td>Tolerance for unboundedness detection heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>InitialLMest</td><td>OT_BOOLEAN</td><td>True</td><td>Enable initial Lagrange multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>KeepAcceptableSol</td><td>OT_BOOLEAN</td><td>True</td><td>Save acceptable solutions as fallback</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LMestQPipComTol</td><td>OT_REAL</td><td>0.003</td><td>IP complementarity tolerance in initial multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LMestQPipResTol</td><td>OT_REAL</td><td>1.0</td><td>IP residual tolerance in initial multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LinMult</td><td>OT_BOOLEAN</td><td>False</td><td>Control Lagrange multiplier update</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LogLevel</td><td>OT_INTEGER</td><td>0</td><td>Enable XML logfiles and writing interval</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LogResult</td><td>OT_INTEGER</td><td>0</td><td>Enable XML result logging and detail level</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaF</td><td>OT_REAL</td><td>0.95</td><td>Lowpass-filter update factor for objective values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaG</td><td>OT_REAL</td><td>0.95</td><td>Lowpass-filter update factor for constraint values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaMerit</td><td>OT_REAL</td><td>0.1</td><td>Lowpass-filter update factor for merit function values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassFilter</td><td>OT_BOOLEAN</td><td>True</td><td>Enable lowpass-filter termination criterion</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MAPivotThreshold</td><td>OT_REAL</td><td>1e-06</td><td>Pivoting tolerance for MA solvers</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MatrixCC</td><td>OT_BOOLEAN</td><td>False</td><td>Not to be included into a parameter file!</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxCalls</td><td>OT_INTEGER</td><td>2147483647</td><td>Upper bound to Reverse Communication calls</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxForce</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Force recovery strategy steps</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxGPart</td><td>OT_INTEGER</td><td>1</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxIter</td><td>OT_INTEGER</td><td>500</td><td>Upper bound on major iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxLScounter</td><td>OT_INTEGER</td><td>3</td><td>Control activation of Filter acceleration heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxNorm</td><td>OT_BOOLEAN</td><td>True</td><td>Select max-norm instead of 1-norm in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MeritFunction</td><td>OT_INTEGER</td><td>4</td><td>Select merit function and penalty update [0, 3..5]</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MeritGradTol</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Threshold of meritfunction gradient for increasing Hessian regularisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MinBettsTau</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Lower bound for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MoreRelax</td><td>OT_BOOLEAN</td><td>False</td><td>Introduce one relaxation variable for every constraint</td><td>casadi::WorhpInternal</td></tr>
<tr><td>NLPmethod</td><td>OT_INTEGER</td><td>1</td><td>Select (1) Meritfunction or (3) Filter globalisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>NLPprint</td><td>OT_INTEGER</td><td>2</td><td>NLP print level [-1..4]</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PairMethod</td><td>OT_INTEGER</td><td>1</td><td>Select method to determine graph colouring pairgroups</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsBar</td><td>OT_REAL</td><td>0.9</td><td>Penalty update parameter factor for MeritFunction = 3</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKFac</td><td>OT_REAL</td><td>2.0</td><td>Penalty update parameter factor for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKSequence</td><td>OT_INTEGER</td><td>2</td><td>Penalty update parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxDeltaK</td><td>OT_REAL</td><td>11.0</td><td>Max penalty for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxFac</td><td>OT_REAL</td><td>100000000.0</td><td>Max factor for increasing penalty for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdRBar</td><td>OT_REAL</td><td>2.0</td><td>Penalty update parameter for MeritFunction = 3</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PrecisionF</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>(currently unused) Relative precision of objective</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PrecisionG</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>(currently unused) Relative precision of constraints</td><td>casadi::WorhpInternal</td></tr>
<tr><td>QPscaleParam</td><td>OT_REAL</td><td>0.0</td><td>(currently unused) Scaling factor for QP</td><td>casadi::WorhpInternal</td></tr>
<tr><td>QuadraticProblem</td><td>OT_BOOLEAN</td><td>False</td><td>Not to be included into a parameter file!</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ReduceBettsTau</td><td>OT_REAL</td><td>0.3</td><td>Decrease factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RegStrategy</td><td>OT_INTEGER</td><td>1</td><td>Select Hessian regularisation strategy in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ReinitFilter</td><td>OT_BOOLEAN</td><td>False</td><td>Enables Filter-reinitialisation accelerating heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxMaxDelta</td><td>OT_REAL</td><td>0.92</td><td>Upper bound for accepting the constraint relaxation variable</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxMaxPen</td><td>OT_REAL</td><td>50000000.0</td><td>Upper bound on the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxRho</td><td>OT_REAL</td><td>6.0</td><td>Update factor for the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxStart</td><td>OT_REAL</td><td>1.0</td><td>Initial value of the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RestUntilFeas</td><td>OT_BOOLEAN</td><td>False</td><td>Do restoration until a feasible solution is found</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleConIter</td><td>OT_BOOLEAN</td><td>False</td><td>Scale constraints in every iteration</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleFacObj</td><td>OT_REAL</td><td>10.0</td><td>Value to scale large objective functions to</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleFacQP</td><td>OT_REAL</td><td>10.0</td><td>Upper bound on resulting matrix norm for QP scaling</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledFD</td><td>OT_BOOLEAN</td><td>True</td><td>Use a scaled perturbation for finite differences</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledKKT</td><td>OT_BOOLEAN</td><td>True</td><td>Scale KKT conditions</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledObj</td><td>OT_BOOLEAN</td><td>True</td><td>Scale the objective function</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledQP</td><td>OT_BOOLEAN</td><td>True</td><td>Scale some matrices handed to the QP</td><td>casadi::WorhpInternal</td></tr>
<tr><td>StartBettsTau</td><td>OT_REAL</td><td>0.1</td><td>Initial value for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingDelta</td><td>OT_REAL</td><td>0.01</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingSCV</td><td>OT_REAL</td><td>1.1</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingSF</td><td>OT_REAL</td><td>2.3</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TakeQPSol</td><td>OT_BOOLEAN</td><td>False</td><td>Evaluate QP search direction regardless of convergence</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Timeout</td><td>OT_REAL</td><td>300.0</td><td>Timeout in seconds</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolComp</td><td>OT_REAL</td><td>0.001</td><td>Complementarity tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolFeas</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolOpti</td><td>OT_REAL</td><td>1e-06</td><td>Optimality tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolWeakActive</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBig</td><td>OT_BOOLEAN</td><td>True</td><td>Enable too-big termination heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBigCV</td><td>OT_REAL</td><td>1e+25</td><td>Upper bound on constraint violation for too-big heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBigKKT</td><td>OT_REAL</td><td>1e+30</td><td>Upper bound on KKT values for too-big heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserDF</td><td>OT_BOOLEAN</td><td>True</td><td>Objective gradient values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserDG</td><td>OT_BOOLEAN</td><td>True</td><td>Jacobian values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>True</td><td>Hessian values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserHMstructure</td><td>OT_INTEGER</td><td>2</td><td>Enable automatic Hessian structure generation or checking</td><td>casadi::WorhpInternal</td></tr>
<tr><td>WeakActiveSet</td><td>OT_BOOLEAN</td><td>False</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Machine epsilon</td><td>casadi::WorhpInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>internalParChanged</td><td>OT_INTEGER</td><td>0</td><td>Counter for changed parameters. Internal use only.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />Monitor functions (eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::WorhpInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>Print information about execution time</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>7.8</td><td>IP barrier parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>2e-07</td><td>IP complementarity tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>0.88</td><td>IP fraction-to-the-boundary parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the IP method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>1e-12</td><td>IP line search minimum step size.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>2.0</td><td>The relaxation term is divided by this value if successful.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>1e-07</td><td>Maximum relaxation value.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>1e-07</td><td>Mimimum relaxation value.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>10.0</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>5e-08</td><td>IP residuals tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>True</td><td>Enable relaxation strategy when encountering an error.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>None</td><td>Select the iterative linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>None</td><td>Select preconditioner for the iterative linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>True</td><td>Enables scaling on linear solver level.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>1e-12</td><td>Tolerance for the linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>False</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>80</td><td>Imposes an upper limit on the number of minor solver iterations, i.e. for the quadratic subproblem solver.If the limit is reached before convergence,WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>None</td><td>Select the solution method used by the QP solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>0.9</td><td>NSN stepsize decrease factor.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>True</td><td>Enable gradient steps in the NSN method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>1e-06</td><td>NSN KKT tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the NSN method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>1e-11</td><td>NSN line search minimum step size.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>0.01</td><td>NSN line search slope parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>None</td><td>Controls the amount of QP solver output.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>False</td><td>Enable scaling on QP level.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>True</td><td>Use strict termination criteria in IP method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::WorhpSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>AcceptTolFeas</td><td>OT_REAL</td><td>0.001</td><td>Tolerance for acceptable feasibility</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AcceptTolOpti</td><td>OT_REAL</td><td>0.001</td><td>Tolerance for acceptable optimality</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AlphaMinConst</td><td>OT_BOOLEAN</td><td>False</td><td>Use a constant lower bound on Armijo stepsize in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Ares</td><td>OT_INTEGERVECTOR</td><td>[42, 41, 42, 43, 44, 41, 50]</td><td>Armijo recovery strategies. Vector of size 7</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoBeta</td><td>OT_REAL</td><td>0.712</td><td>Trial stepsize decrease factor for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMaxAlpha</td><td>OT_REAL</td><td>1.0</td><td>Initial alpha for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlpha</td><td>OT_REAL</td><td>1e-06</td><td>Lower bound on alpha for Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlphaRec</td><td>OT_REAL</td><td>1e-06</td><td>Lower bound on alpha for Armijo rule during recovery</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ArmijoSigma</td><td>OT_REAL</td><td>0.005</td><td>Scale factor for linearised descent check in Armijo rule</td><td>casadi::WorhpInternal</td></tr>
<tr><td>AutoQPRecovery</td><td>OT_BOOLEAN</td><td>True</td><td>Enable automatic QP recovery</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSmaxblockSize</td><td>OT_INTEGER</td><td>300</td><td>Block size parameter used by certain BFGS methods</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSmethod</td><td>OT_INTEGER</td><td>0</td><td>Choose BFGS method (0: dense, 1-3: block, 100+: sparse)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSminblockSize</td><td>OT_INTEGER</td><td>300</td><td>Block size parameter used by certain BFGS methods</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BFGSrestart</td><td>OT_INTEGER</td><td>50</td><td>Restart BFGS update after this many iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BettsFactor</td><td>OT_REAL</td><td>2.1</td><td>Update factor for Betts' Hessian regularisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BettsPoint</td><td>OT_REAL</td><td>1.0</td><td>Smallest eigenvalue of the regularised Hessian</td><td>casadi::WorhpInternal</td></tr>
<tr><td>BoundTolFac</td><td>OT_REAL</td><td>1000.0</td><td>Factor in determining active constraints by KKT</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckFJ</td><td>OT_REAL</td><td>1e+12</td><td>Upper bound used by Fritz-John heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureDF</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of DF</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureDG</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of DG</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CheckStructureHM</td><td>OT_BOOLEAN</td><td>True</td><td>Enable structural checking of HM</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepBettsSum</td><td>OT_REAL</td><td>0.5</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepConStop</td><td>OT_REAL</td><td>1e-06</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepConvio</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMaxIter</td><td>OT_INTEGER</td><td>50</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMethod</td><td>OT_INTEGER</td><td>0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepMode</td><td>OT_INTEGER</td><td>1</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepPFactor</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepPMax</td><td>OT_REAL</td><td>1000000.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CorStepRecoveryDX</td><td>OT_BOOLEAN</td><td>False</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvBCond</td><td>OT_REAL</td><td>0.02</td><td>Block BFGS curvature condition bound</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvBFac</td><td>OT_REAL</td><td>0.3</td><td>Block BFGS curvature condition regularisation factor</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvCond</td><td>OT_REAL</td><td>0.02</td><td>BFGS Curvature condition bound</td><td>casadi::WorhpInternal</td></tr>
<tr><td>CurvFac</td><td>OT_REAL</td><td>0.3</td><td>BFGS curvature condition regularisation factor</td><td>casadi::WorhpInternal</td></tr>
<tr><td>DebugMarker05</td><td>OT_INTEGER</td><td>42</td><td>Debug marker. Used to find memory alignment/padding issues</td><td>casadi::WorhpInternal</td></tr>
<tr><td>DebugMarker06</td><td>OT_INTEGER</td><td>42</td><td>Debug marker. Used to find memory alignment/padding issues</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FGtogether</td><td>OT_BOOLEAN</td><td>False</td><td>F and G cannot be evaluated separately</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FJandND</td><td>OT_BOOLEAN</td><td>False</td><td>Enable Fritz-John and non-differentiable check heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleDual</td><td>OT_BOOLEAN</td><td>False</td><td>Activate dual feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleInit</td><td>OT_BOOLEAN</td><td>False</td><td>Activate initial feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleInitTol</td><td>OT_REAL</td><td>0.001</td><td>Feasibility tolerance for no-objective feasible mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FeasibleOnly</td><td>OT_BOOLEAN</td><td>False</td><td>Activate feasible-only mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FidifEps</td><td>OT_REAL</td><td>1e-05</td><td>Finite difference perturbation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FidifHM</td><td>OT_BOOLEAN</td><td>False</td><td>Approximate Hessian by finite differences (otherwise BFGS)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterBisecAlpha</td><td>OT_BOOLEAN</td><td>True</td><td>Filter heuristic to save Armijo iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterGammaCV</td><td>OT_REAL</td><td>7.5e-06</td><td>Constraint violation decrease factor in Filter acceptance check</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterGammaF</td><td>OT_REAL</td><td>1.1e-05</td><td>Objective decrease factor in Filter acceptance check</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FilterIntersecAlpha</td><td>OT_BOOLEAN</td><td>True</td><td>Filter heuristic to save Armijo iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FirstDifCentral</td><td>OT_BOOLEAN</td><td>True</td><td>Use central finite difference quotient for first derivatives</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FocusOnFeas</td><td>OT_BOOLEAN</td><td>True</td><td>Enable Focus-on-Feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>FocusOnFeasFactor</td><td>OT_REAL</td><td>1.36</td><td>Factor in Focus-on-Feasibility mode</td><td>casadi::WorhpInternal</td></tr>
<tr><td>GammaAlpha</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for alphamin calculation by Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>GroupMethod</td><td>OT_INTEGER</td><td>1</td><td>Select method to determine graph colouring groups</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IgnoreFilterCrit</td><td>OT_BOOLEAN</td><td>False</td><td>Activate accelerating heuristics for Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncBettsTau</td><td>OT_REAL</td><td>2.0</td><td>Increase factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncBettsTauMore</td><td>OT_REAL</td><td>100.0</td><td>Larger increase factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncreaseIWS</td><td>OT_REAL</td><td>1.0</td><td>Increase factor for estimated integer workspace requirement</td><td>casadi::WorhpInternal</td></tr>
<tr><td>IncreaseRWS</td><td>OT_REAL</td><td>1.0</td><td>Increase factor for estimated real workspace requirement</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Infty</td><td>OT_REAL</td><td>1e+20</td><td>Upper bound for numbers to be regarded as finite</td><td>casadi::WorhpInternal</td></tr>
<tr><td>InftyUnbounded</td><td>OT_REAL</td><td>1e+20</td><td>Tolerance for unboundedness detection heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>InitialLMest</td><td>OT_BOOLEAN</td><td>True</td><td>Enable initial Lagrange multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>KeepAcceptableSol</td><td>OT_BOOLEAN</td><td>True</td><td>Save acceptable solutions as fallback</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LMestQPipComTol</td><td>OT_REAL</td><td>0.003</td><td>IP complementarity tolerance in initial multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LMestQPipResTol</td><td>OT_REAL</td><td>1.0</td><td>IP residual tolerance in initial multiplier estimate</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LinMult</td><td>OT_BOOLEAN</td><td>False</td><td>Control Lagrange multiplier update</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LogLevel</td><td>OT_INTEGER</td><td>0</td><td>Enable XML logfiles and writing interval</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LogResult</td><td>OT_INTEGER</td><td>0</td><td>Enable XML result logging and detail level</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaF</td><td>OT_REAL</td><td>0.95</td><td>Lowpass-filter update factor for objective values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaG</td><td>OT_REAL</td><td>0.95</td><td>Lowpass-filter update factor for constraint values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaMerit</td><td>OT_REAL</td><td>0.1</td><td>Lowpass-filter update factor for merit function values</td><td>casadi::WorhpInternal</td></tr>
<tr><td>LowPassFilter</td><td>OT_BOOLEAN</td><td>True</td><td>Enable lowpass-filter termination criterion</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MAPivotThreshold</td><td>OT_REAL</td><td>1e-06</td><td>Pivoting tolerance for MA solvers</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MatrixCC</td><td>OT_BOOLEAN</td><td>False</td><td>Not to be included into a parameter file!</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxCalls</td><td>OT_INTEGER</td><td>2147483647</td><td>Upper bound to Reverse Communication calls</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxForce</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Force recovery strategy steps</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxGPart</td><td>OT_INTEGER</td><td>1</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxIter</td><td>OT_INTEGER</td><td>500</td><td>Upper bound on major iterations</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxLScounter</td><td>OT_INTEGER</td><td>3</td><td>Control activation of Filter acceleration heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MaxNorm</td><td>OT_BOOLEAN</td><td>True</td><td>Select max-norm instead of 1-norm in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MeritFunction</td><td>OT_INTEGER</td><td>4</td><td>Select merit function and penalty update [0, 3..5]</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MeritGradTol</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Threshold of meritfunction gradient for increasing Hessian regularisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MinBettsTau</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Lower bound for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>MoreRelax</td><td>OT_BOOLEAN</td><td>False</td><td>Introduce one relaxation variable for every constraint</td><td>casadi::WorhpInternal</td></tr>
<tr><td>NLPmethod</td><td>OT_INTEGER</td><td>1</td><td>Select (1) Meritfunction or (3) Filter globalisation</td><td>casadi::WorhpInternal</td></tr>
<tr><td>NLPprint</td><td>OT_INTEGER</td><td>2</td><td>NLP print level [-1..4]</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PairMethod</td><td>OT_INTEGER</td><td>1</td><td>Select method to determine graph colouring pairgroups</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsBar</td><td>OT_REAL</td><td>0.9</td><td>Penalty update parameter factor for MeritFunction = 3</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKFac</td><td>OT_REAL</td><td>2.0</td><td>Penalty update parameter factor for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKSequence</td><td>OT_INTEGER</td><td>2</td><td>Penalty update parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxDeltaK</td><td>OT_REAL</td><td>11.0</td><td>Max penalty for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxFac</td><td>OT_REAL</td><td>100000000.0</td><td>Max factor for increasing penalty for MeritFunction = 4</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PenUpdRBar</td><td>OT_REAL</td><td>2.0</td><td>Penalty update parameter for MeritFunction = 3</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PrecisionF</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>(currently unused) Relative precision of objective</td><td>casadi::WorhpInternal</td></tr>
<tr><td>PrecisionG</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>(currently unused) Relative precision of constraints</td><td>casadi::WorhpInternal</td></tr>
<tr><td>QPscaleParam</td><td>OT_REAL</td><td>0.0</td><td>(currently unused) Scaling factor for QP</td><td>casadi::WorhpInternal</td></tr>
<tr><td>QuadraticProblem</td><td>OT_BOOLEAN</td><td>False</td><td>Not to be included into a parameter file!</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ReduceBettsTau</td><td>OT_REAL</td><td>0.3</td><td>Decrease factor for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RegStrategy</td><td>OT_INTEGER</td><td>1</td><td>Select Hessian regularisation strategy in Filter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ReinitFilter</td><td>OT_BOOLEAN</td><td>False</td><td>Enables Filter-reinitialisation accelerating heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxMaxDelta</td><td>OT_REAL</td><td>0.92</td><td>Upper bound for accepting the constraint relaxation variable</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxMaxPen</td><td>OT_REAL</td><td>50000000.0</td><td>Upper bound on the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxRho</td><td>OT_REAL</td><td>6.0</td><td>Update factor for the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RelaxStart</td><td>OT_REAL</td><td>1.0</td><td>Initial value of the constraint relaxation penalty</td><td>casadi::WorhpInternal</td></tr>
<tr><td>RestUntilFeas</td><td>OT_BOOLEAN</td><td>False</td><td>Do restoration until a feasible solution is found</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleConIter</td><td>OT_BOOLEAN</td><td>False</td><td>Scale constraints in every iteration</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleFacObj</td><td>OT_REAL</td><td>10.0</td><td>Value to scale large objective functions to</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaleFacQP</td><td>OT_REAL</td><td>10.0</td><td>Upper bound on resulting matrix norm for QP scaling</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledFD</td><td>OT_BOOLEAN</td><td>True</td><td>Use a scaled perturbation for finite differences</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledKKT</td><td>OT_BOOLEAN</td><td>True</td><td>Scale KKT conditions</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledObj</td><td>OT_BOOLEAN</td><td>True</td><td>Scale the objective function</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ScaledQP</td><td>OT_BOOLEAN</td><td>True</td><td>Scale some matrices handed to the QP</td><td>casadi::WorhpInternal</td></tr>
<tr><td>StartBettsTau</td><td>OT_REAL</td><td>0.1</td><td>Initial value for Betts' update dampening term</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingDelta</td><td>OT_REAL</td><td>0.01</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingSCV</td><td>OT_REAL</td><td>1.1</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>SwitchingSF</td><td>OT_REAL</td><td>2.3</td><td>Filter switching condition parameter</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TakeQPSol</td><td>OT_BOOLEAN</td><td>False</td><td>Evaluate QP search direction regardless of convergence</td><td>casadi::WorhpInternal</td></tr>
<tr><td>Timeout</td><td>OT_REAL</td><td>300.0</td><td>Timeout in seconds</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolComp</td><td>OT_REAL</td><td>0.001</td><td>Complementarity tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolFeas</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolOpti</td><td>OT_REAL</td><td>1e-06</td><td>Optimality tolerance</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TolWeakActive</td><td>OT_REAL</td><td>1.0</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBig</td><td>OT_BOOLEAN</td><td>True</td><td>Enable too-big termination heuristics</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBigCV</td><td>OT_REAL</td><td>1e+25</td><td>Upper bound on constraint violation for too-big heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>TooBigKKT</td><td>OT_REAL</td><td>1e+30</td><td>Upper bound on KKT values for too-big heuristic</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserDF</td><td>OT_BOOLEAN</td><td>True</td><td>Objective gradient values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserDG</td><td>OT_BOOLEAN</td><td>True</td><td>Jacobian values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>True</td><td>Hessian values supplied by caller</td><td>casadi::WorhpInternal</td></tr>
<tr><td>UserHMstructure</td><td>OT_INTEGER</td><td>2</td><td>Enable automatic Hessian structure generation or checking</td><td>casadi::WorhpInternal</td></tr>
<tr><td>WeakActiveSet</td><td>OT_BOOLEAN</td><td>False</td><td>(experimental)</td><td>casadi::WorhpInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps</td><td>OT_REAL</td><td>2.22044604925e-16</td><td>Machine epsilon</td><td>casadi::WorhpInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>internalParChanged</td><td>OT_INTEGER</td><td>0</td><td>Counter for changed parameters. Internal use only.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />Monitor functions (eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::WorhpInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>Print information about execution time</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>7.8</td><td>IP barrier parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>2e-07</td><td>IP complementarity tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>0.88</td><td>IP fraction-to-the-boundary parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the IP method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>1e-12</td><td>IP line search minimum step size.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>2.0</td><td>The relaxation term is divided by this value if successful.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>1e-07</td><td>Maximum relaxation value.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>1e-07</td><td>Mimimum relaxation value.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>10.0</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>5e-08</td><td>IP residuals tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>True</td><td>Enable relaxation strategy when encountering an error.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>None</td><td>Select the iterative linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>None</td><td>Select preconditioner for the iterative linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>True</td><td>Enables scaling on linear solver level.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>1e-12</td><td>Tolerance for the linear solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>False</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>80</td><td>Imposes an upper limit on the number of minor solver iterations, i.e. for the quadratic subproblem solver.If the limit is reached before convergence,WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>None</td><td>Select the solution method used by the QP solver.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>0.9</td><td>NSN stepsize decrease factor.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>True</td><td>Enable gradient steps in the NSN method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>1e-06</td><td>NSN KKT tolerance.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the NSN method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>1e-11</td><td>NSN line search minimum step size.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>0.01</td><td>NSN line search slope parameter.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>None</td><td>Controls the amount of QP solver output.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>False</td><td>Enable scaling on QP level.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>True</td><td>Use strict termination criteria in IP method.</td><td>casadi::WorhpInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::FunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Function
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CSparseInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CSparse
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::OldCollocationIntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the implicit solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::OldCollocationIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the implicit solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>casadi::OldCollocationIntegratorInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NLPQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSOlver used to solve the QPs.</td><td>casadi::NLPQPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSOlver</td><td>casadi::NLPQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::NLPQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSOlver used to solve the QPs.</td><td>casadi::NLPQPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSOlver</td><td>casadi::NLPQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SDPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SDPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SDPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NullspaceInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dense</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates that dense matrices can be assumed</td><td>casadi::NullspaceInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Nullspace
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dense</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates that dense matrices can be assumed</td><td>casadi::NullspaceInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SCPgenInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>gauss-newton|exact</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::SCPgenInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SCPgen
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>gauss-newton|exact</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::SCPgenInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>casadi::SCPgenInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LinearSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::LinearSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::IpoptInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_after_max_steps</td><td>OT_INTEGER</td><td>-1</td><td>Accept a trial point after maximal this number of steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td>no</td><td>Always accept the first trial step. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the complementarity conditions. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td>10000000000.0</td><td>"Acceptance" threshold for the dual infeasibility. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td>15</td><td>Number of "acceptable" iterates before triggering termination. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td>1e+20</td><td>"Acceptance" stopping criterion based on objective function change. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td>1e-06</td><td>"Acceptable" convergence tolerance (relative). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>adaptive_mu_globalization</td><td>OT_STRING</td><td>obj-constr-filter</td><td>Globalization strategy for the adaptive mu selection mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkt_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for the KKT error in the adaptive mu globalization strategies. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_fact</td><td>OT_REAL</td><td>0.9999</td><td>Sufficient decrease factor for "kkt-error" globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_iters</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of iterations requiring sufficient progress. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_monotone_init_factor</td><td>OT_REAL</td><td>0.8</td><td>Determines the initial value of the barrier parameter when switching to the monotone mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_restore_previous_iterate</td><td>OT_STRING</td><td>no</td><td>Indicates if the previous iterate should be restored if the monotone mode is entered. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_safeguard_factor</td><td>OT_REAL</td><td>0.0</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td>primal</td><td>Method to determine the step size for constraint multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td>10.0</td><td>Tolerance for switching to full equality multiplier steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_min_frac</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for the minimal step size (before switching to restoration phase). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_red_factor</td><td>OT_REAL</td><td>0.5</td><td>Fractional reduction of the trial step size in the backtracking line search. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td>10.0</td><td>Factor for mu in barrier stop test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial point to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td>constant</td><td>Initialization method for bound multipliers (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td>1.0</td><td>Initial value for the bound multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td>1000.0</td><td>Threshold for resetting bound multipliers after the restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial point to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td>1e-08</td><td>Factor for initial relaxation of the bounds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td>no</td><td>Indicates whether it is desired to check for Nan/Inf in derivative matrices (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_cup</td><td>OT_REAL</td><td>1.5</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_hat</td><td>OT_REAL</td><td>2.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_tilde</td><td>OT_REAL</td><td>5.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the complementarity conditions. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>None</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td>1000.0</td><td>Maximum allowed least-square guess of constraint multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for resetting equality and inequality multipliers after restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constraint_violation_norm_type</td><td>OT_STRING</td><td>1-norm</td><td>Norm to be used for the constraint violation in the line search. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>corrector_compl_avrg_red_fact</td><td>OT_REAL</td><td>1.0</td><td>Complementarity tolerance factor for accepting corrector step (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td>none</td><td>The type of corrector steps that should be taken (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>delta</td><td>OT_REAL</td><td>1.0</td><td>Multiplier for constraint violation in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>delta_y_max</td><td>OT_REAL</td><td>1e+12</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dependency_detection_with_rhs</td><td>OT_STRING</td><td>no</td><td>Indicates if the right hand sides of the constraints should be considered during dependency detection (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dependency_detector</td><td>OT_STRING</td><td>none</td><td>Indicates which linear solver should be used to detect linearly dependent equality constraints. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td>none</td><td>Enable derivative checker (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_first_index</td><td>OT_INTEGER</td><td>-2</td><td>Index of first quantity to be checked by derivative checker (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td>1e-08</td><td>Size of the finite difference perturbation in derivative test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td>no</td><td>Indicates whether information for all estimated derivatives should be printed. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td>0.0001</td><td>Threshold for indicating wrong derivative. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td>1e+20</td><td>Threshold for maximal value of primal iterates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td>1.0</td><td>Desired threshold for the dual infeasibility. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>epsilon_c</td><td>OT_REAL</td><td>0.01</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_min</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_penalty</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition for the penalty function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td>yes</td><td>Determines if the original objective function should be evaluated at restoration phase trial points. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td>no</td><td>Enable heuristics to quickly detect an infeasible problem. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td>0.001</td><td>Threshold for disabling "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td>100000000.0</td><td>Multiplier threshold for activating "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fast_des_fact</td><td>OT_REAL</td><td>0.1</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fast_step_computation</td><td>OT_STRING</td><td>no</td><td>Indicates if the linear system should be solved quickly. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td>5</td><td>Verbosity level for output file. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_margin_fact</td><td>OT_REAL</td><td>1e-05</td><td>Factor determining width of margin for obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_max_margin</td><td>OT_REAL</td><td>1.0</td><td>Maximum width of margin in obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_reset_trigger</td><td>OT_INTEGER</td><td>5</td><td>Number of iterations that trigger the filter reset. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>findiff_perturbation</td><td>OT_REAL</td><td>1e-07</td><td>Size of the finite difference perturbation for derivative approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td>0.0001</td><td>Size of first x-s perturbation tried. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td>average_compl</td><td>Oracle for the barrier parameter when switching to fixed mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td>make_parameter</td><td>Determines how fixed variables should be handled. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_hat</td><td>OT_REAL</td><td>0.04</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the filter margin for the barrier function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_theta</td><td>OT_REAL</td><td>1e-05</td><td>Relaxation factor in the filter margin for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_tilde</td><td>OT_REAL</td><td>4.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>exact</td><td>Indicates what Hessian information is to be used. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>hessian_approximation_space</td><td>OT_STRING</td><td>nonlinear-variables</td><td>Indicates in which subspace the Hessian information is to be approximated. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether the problem is a quadratic problem (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td>yes</td><td>Indicates whether final points should be projected into original bounds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inf_pr_output</td><td>OT_STRING</td><td>original</td><td>Determines what value is printed in the "inf_pr" output column. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all equality constraints are linear (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all inequality constraints are linear (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jacobian_approximation</td><td>OT_STRING</td><td>exact</td><td>Specifies technique to compute constraint Jacobian (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_exponent</td><td>OT_REAL</td><td>0.25</td><td>Exponent for mu in the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td>1e-08</td><td>Size of the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_d</td><td>OT_REAL</td><td>1e-05</td><td>Weight for linear damping term (to handle one-sided bounds). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_sigma</td><td>OT_REAL</td><td>10000000000.0</td><td>Factor limiting the deviation of dual variables from primal estimates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_soc</td><td>OT_REAL</td><td>0.99</td><td>Factor in the sufficient reduction rule for second order correction. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_x_dis</td><td>OT_REAL</td><td>100.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_y_dis</td><td>OT_REAL</td><td>10000.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>least_square_init_duals</td><td>OT_STRING</td><td>no</td><td>Least square initialization of all dual variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>least_square_init_primal</td><td>OT_STRING</td><td>no</td><td>Least square initialization of the primal variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_aug_solver</td><td>OT_STRING</td><td>sherman-morrison</td><td>Strategy for solving the augmented system for low-rank Hessian. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val</td><td>OT_REAL</td><td>1.0</td><td>Value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_max</td><td>OT_REAL</td><td>100000000.0</td><td>Upper bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_min</td><td>OT_REAL</td><td>1e-08</td><td>Lower bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_initialization</td><td>OT_STRING</td><td>scalar1</td><td>Initialization strategy for the limited memory quasi-Newton approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td>6</td><td>Maximum size of the history for the limited quasi-Newton Hessian approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td>2</td><td>Threshold for successive iterations where update is skipped. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_special_for_resto</td><td>OT_STRING</td><td>no</td><td>Determines if the quasi-Newton updates should be special during the restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_update_type</td><td>OT_STRING</td><td>bfgs</td><td>Quasi-Newton update formula for the limited memory approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>line_search_method</td><td>OT_STRING</td><td>filter</td><td>Globalization method used in backtracking line search (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td>yes</td><td>Flag indicating that linear scaling is only done if it seems required. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>ma27</td><td>Linear solver used for step computations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td>mc19</td><td>Method for scaling the linear system. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_ignore_singularity</td><td>OT_STRING</td><td>no</td><td>Enables MA27's ability to solve a linear system even if the matrix is singular. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Real workspace memory for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Integer workspace memory for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td>10.0</td><td>Increment factor for workspace size for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma28_pivtol</td><td>OT_REAL</td><td>0.01</td><td>Pivot tolerance for linear solver MA28. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td>yes</td><td>Controls MA57 automatic scaling (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td>16</td><td>Controls block size used by Level 3 BLAS in MA57BD (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td>16</td><td>Node amalgamation parameter (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td>5</td><td>Controls pivot order in MA57 (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td>1.05</td><td>Safety factor for work space memory allocation for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_small_pivot_flag</td><td>OT_INTEGER</td><td>0</td><td>If set to 1, then when small entries defined by CNTL(2) are detected they are removed and the corresponding pivots placed at the end of the factorization.  This can be particularly efficient if the matrix is highly rank deficient. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_nemin</td><td>OT_INTEGER</td><td>32</td><td>Node Amalgamation parameter (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_print_level</td><td>OT_INTEGER</td><td>0</td><td>Debug printing level for the linear solver MA86 (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_small</td><td>OT_REAL</td><td>1e-20</td><td>Zero Pivot Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_static</td><td>OT_REAL</td><td>0.0</td><td>Static Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_u</td><td>OT_REAL</td><td>1e-08</td><td>Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_umax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>magic_steps</td><td>OT_STRING</td><td>no</td><td>Enables magic steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum number of CPU seconds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_filter_resets</td><td>OT_INTEGER</td><td>5</td><td>Maximal allowed number of filter resets (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td>1e+20</td><td>Maximum value of regularization parameter for handling negative curvature. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>3000</td><td>Maximum number of iterations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_resto_iter</td><td>OT_INTEGER</td><td>3000000</td><td>Maximum number of successive iterations in restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of second order correction trial steps at each iteration. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_soft_resto_iters</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterations performed successively in soft restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td>no</td><td>Indicates if we want to do Mehrotra's algorithm. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_alpha_primal</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td>1e-20</td><td>Smallest perturbation of the Hessian block. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td>1</td><td>Minimum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::IpoptInternal</td></tr>
<tr><td>mu_allow_fast_monotone_decrease</td><td>OT_STRING</td><td>yes</td><td>Allow skipping of barrier problem if barrier test is already met. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td>0.1</td><td>Initial value for the barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td>0.2</td><td>Determines linear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td>1000.0</td><td>Factor for initialization of maximum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td>1e-11</td><td>Minimum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td>quality-function</td><td>Oracle for a new barrier parameter in the adaptive strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td>monotone</td><td>Update strategy for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td>1.5</td><td>Determines superlinear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td>0.0</td><td>Desired value of complementarity. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mult_diverg_feasibility_tol</td><td>OT_REAL</td><td>1e-07</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mult_diverg_y_tol</td><td>OT_REAL</td><td>100000000.0</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_dep_tol</td><td>OT_REAL</td><td>-1.0</td><td>Pivot threshold for detection of linearly dependent constraints in MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td>1000</td><td>Percentage increase in the estimated working space for MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td>7</td><td>Controls permuting and scaling in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td>7</td><td>Controls pivot order in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td>1e-06</td><td>Pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td>77</td><td>Controls scaling in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>neg_curv_test_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for heuristic to ignore wrong inertia. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>never_use_fact_cgpen_direction</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the fast Chen-Goldfarb direction (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>never_use_piecewise_penalty_ls</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the piecewise penalty method (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td>-1e+19</td><td>any bound less or equal this value will be considered -inf (i.e. not lower bounded). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_constr_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for constraint function gradient size. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_max_gradient</td><td>OT_REAL</td><td>100.0</td><td>Maximum gradient after NLP scaling. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_method</td><td>OT_STRING</td><td>gradient-based</td><td>Select the technique used for scaling the NLP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_min_value</td><td>OT_REAL</td><td>1e-08</td><td>Minimum value of gradient-based scaling values. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_obj_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for objective function gradient size. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td>1e+19</td><td>any bound greater or this value will be considered +inf (i.e. not upper bounded). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nu_inc</td><td>OT_REAL</td><td>0.0001</td><td>Increment of the penalty parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nu_init</td><td>OT_REAL</td><td>1e-06</td><td>Initial value of the penalty parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>num_linear_variables</td><td>OT_INTEGER</td><td>0</td><td>Number of linear variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>obj_max_inc</td><td>OT_REAL</td><td>5.0</td><td>Determines the upper bound on the acceptable increase of barrier objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>obj_scaling_factor</td><td>OT_REAL</td><td>1.0</td><td>Scaling factor for the objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td>File name of options file (to overwrite default). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td>File name of desired output file (leave unset for no file output). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>pardiso_iter_coarse_size</td><td>OT_INTEGER</td><td>5000</td><td>Maximum Size of Coarse Grid Matrix (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_factor</td><td>OT_REAL</td><td>0.5</td><td>dropping value for incomplete factor (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_schur</td><td>OT_REAL</td><td>0.1</td><td>dropping value for sparsify schur complement factor (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_inverse_norm_factor</td><td>OT_REAL</td><td>5000000.0</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_levels</td><td>OT_INTEGER</td><td>10</td><td>Maximum Size of Grid Levels (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_row_fill</td><td>OT_INTEGER</td><td>10000000</td><td>max fill for each row (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_relative_tol</td><td>OT_REAL</td><td>1e-06</td><td>Relative Residual Convergence (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iterative</td><td>OT_STRING</td><td>no</td><td>Switch on iterative solver in Pardiso library (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td>complete+2x2</td><td>Matching strategy to be used by Pardiso (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_max_droptol_corrections</td><td>OT_INTEGER</td><td>4</td><td>Maximal number of decreases of drop tolerance during one solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_max_iter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of Krylov-Subspace Iteration (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td>0</td><td>Pardiso message level (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td>0</td><td>Enables out-of-core variant of Pardiso (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_redo_symbolic_fact_only_if_inertia_wrong</td><td>OT_STRING</td><td>no</td><td>Toggle for handling case when elements were perturbed by Pardiso. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_repeated_perturbation_means_singular</td><td>OT_STRING</td><td>no</td><td>Interpretation of perturbed elements. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>False</td><td>n/a</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_des_fact</td><td>OT_REAL</td><td>0.2</td><td>a parameter used in penalty parameter computation (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_init_fac</td><td>OT_REAL</td><td>50.0</td><td>a parameter used to choose initial penalty parameterswhen the regularized Newton method is used. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_init_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximal value for the intial penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_init_min</td><td>OT_REAL</td><td>1.0</td><td>Minimal value for the intial penalty parameter for line search(for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_max</td><td>OT_REAL</td><td>1e+30</td><td>Maximal value for the penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_update_compl_tol</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_update_infeasibility_tol</td><td>OT_REAL</td><td>1e-09</td><td>Threshold for infeasibility in penalty parameter update test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_always_cd</td><td>OT_STRING</td><td>no</td><td>Active permanent perturbation of constraint linearization. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td>0.333333333333</td><td>Decrease factor for x-s perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td>8.0</td><td>Increase factor for x-s perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td>100.0</td><td>Increase factor for x-s perturbation for very first perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_infeasi</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_obj</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td>10.0</td><td>Maximal perturbation of an evaluation point. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_info_string</td><td>OT_STRING</td><td>no</td><td>Enables printing of additional info string at end of iteration output. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>5</td><td>Output verbosity level. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td>no</td><td>Switch to print all algorithmic options. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_options_latex_mode</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>print information about execution time</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_timing_statistics</td><td>OT_STRING</td><td>no</td><td>Switch to print timing statistics. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td>no</td><td>Print all options set by the user. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_balancing_term</td><td>OT_STRING</td><td>none</td><td>The balancing term included in the quality function for centrality. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_centrality</td><td>OT_STRING</td><td>none</td><td>The penalty term for centrality that is included in quality function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td>8</td><td>Maximum number of search steps during direct search procedure determining the optimal centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for components of the quality function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_section_qf_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for the golden section search procedure determining the optimal centering parameter (in the function value space). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_section_sigma_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for the section search procedure determining the optimal centering parameter (in sigma space). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td>no</td><td>Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility threshold for recomputation of multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>replace_bounds</td><td>OT_STRING</td><td>no</td><td>Indicates if all variable bounds should be replaced by inequality constraints (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td>0.9</td><td>Required reduction of infeasibility before leaving restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_improvement_factor</td><td>OT_REAL</td><td>0.999999999</td><td>Minimal required reduction of residual test ratio in iterative refinement. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_ratio_max</td><td>OT_REAL</td><td>1e-10</td><td>Iterative refinement tolerance (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_ratio_singular</td><td>OT_REAL</td><td>1e-05</td><td>Threshold for declaring linear system singular after failed iterative refinement. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_failure_feasibility_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for primal infeasibility to declare failure of restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_penalty_parameter</td><td>OT_REAL</td><td>1000.0</td><td>Penalty parameter in the restoration phase objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_proximity_weight</td><td>OT_REAL</td><td>1.0</td><td>Weighting factor for the proximity term in restoration phase objective. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.1</td><td>Value in penalty parameter update formula. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_max</td><td>OT_REAL</td><td>100.0</td><td>Scaling threshold for the NLP error. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_phi</td><td>OT_REAL</td><td>2.3</td><td>Exponent for linear barrier function model in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_theta</td><td>OT_REAL</td><td>1.1</td><td>Exponent for current constraint violation in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sb</td><td>OT_STRING</td><td>no</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sigma_max</td><td>OT_REAL</td><td>100.0</td><td>Maximum value of the centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sigma_min</td><td>OT_REAL</td><td>1e-06</td><td>Minimum value of the centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_corr_if_neg_curv</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step in negative curvature iteration (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_corr_in_monotone_mode</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step during monotone barrier parameter mode (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_finalize_solution_call</td><td>OT_STRING</td><td>no</td><td>Indicates if call to NLP::FinalizeSolution after optimization should be suppressed (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial slack to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial slack to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_move</td><td>OT_REAL</td><td>1.81898940355e-12</td><td>Correction size for very small slacks. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td>0.9999</td><td>Required reduction in primal-dual error in the soft restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td>no</td><td>Tells algorithm to switch to restoration phase in first iteration. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>suppress_all_output</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tau_min</td><td>OT_REAL</td><td>0.99</td><td>Lower bound on fraction-to-the-boundary parameter tau. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_min</td><td>OT_REAL</td><td>1e-06</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_min_fact</td><td>OT_REAL</td><td>0.0001</td><td>Determines constraint violation threshold in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tiny_step_tol</td><td>OT_REAL</td><td>2.22044604925e-15</td><td>Tolerance for detecting numerically insignificant steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tiny_step_y_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for quitting because of numerically insignificant steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-08</td><td>Desired convergence tolerance (relative). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>None</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>vartheta</td><td>OT_REAL</td><td>0.5</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as bound_frac for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_entire_iterate</td><td>OT_STRING</td><td>no</td><td>Tells algorithm whether to use the GetWarmStartIterate method in the NLP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td>no</td><td>Warm-start for initial point (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as mult_bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum initial value for the equality multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_same_structure</td><td>OT_STRING</td><td>no</td><td>Indicates whether a problem with a structure identical to the previous one is to be solved. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_frac for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_target_mu</td><td>OT_REAL</td><td>0.0</td><td>Unsupported! (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td>10</td><td>Number of shortened iterations that trigger the watchdog. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of watchdog iterations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_droptol</td><td>OT_REAL</td><td>0.0</td><td>Drop tolerance for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_fillin_limit</td><td>OT_REAL</td><td>0.0</td><td>Fill-in limit for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_iterative</td><td>OT_STRING</td><td>no</td><td>Switches to iterative solver in WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximal number of iterations in iterative WISMP (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_no_pivoting</td><td>OT_STRING</td><td>no</td><td>Use the static pivoting option of WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td>1</td><td>Number of threads to be used in WSMP (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(16) (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option2</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(20) (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td>0.0001</td><td>Pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td>0</td><td>Determines how the matrix is scaled by WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td>1e-18</td><td>WSMP's singularity threshold. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretent inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_write_matrix_iteration</td><td>OT_INTEGER</td><td>-1</td><td>Iteration in which the matrices are written to files. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::IpoptSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_after_max_steps</td><td>OT_INTEGER</td><td>-1</td><td>Accept a trial point after maximal this number of steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td>no</td><td>Always accept the first trial step. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the complementarity conditions. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td>10000000000.0</td><td>"Acceptance" threshold for the dual infeasibility. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td>15</td><td>Number of "acceptable" iterates before triggering termination. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td>1e+20</td><td>"Acceptance" stopping criterion based on objective function change. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td>1e-06</td><td>"Acceptable" convergence tolerance (relative). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>adaptive_mu_globalization</td><td>OT_STRING</td><td>obj-constr-filter</td><td>Globalization strategy for the adaptive mu selection mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkt_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for the KKT error in the adaptive mu globalization strategies. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_fact</td><td>OT_REAL</td><td>0.9999</td><td>Sufficient decrease factor for "kkt-error" globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_iters</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of iterations requiring sufficient progress. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_monotone_init_factor</td><td>OT_REAL</td><td>0.8</td><td>Determines the initial value of the barrier parameter when switching to the monotone mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_restore_previous_iterate</td><td>OT_STRING</td><td>no</td><td>Indicates if the previous iterate should be restored if the monotone mode is entered. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_safeguard_factor</td><td>OT_REAL</td><td>0.0</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td>primal</td><td>Method to determine the step size for constraint multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td>10.0</td><td>Tolerance for switching to full equality multiplier steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_min_frac</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for the minimal step size (before switching to restoration phase). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>alpha_red_factor</td><td>OT_REAL</td><td>0.5</td><td>Fractional reduction of the trial step size in the backtracking line search. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td>10.0</td><td>Factor for mu in barrier stop test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial point to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td>constant</td><td>Initialization method for bound multipliers (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td>1.0</td><td>Initial value for the bound multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td>1000.0</td><td>Threshold for resetting bound multipliers after the restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial point to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td>1e-08</td><td>Factor for initial relaxation of the bounds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td>no</td><td>Indicates whether it is desired to check for Nan/Inf in derivative matrices (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_cup</td><td>OT_REAL</td><td>1.5</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_hat</td><td>OT_REAL</td><td>2.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>chi_tilde</td><td>OT_REAL</td><td>5.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the complementarity conditions. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>None</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td>1000.0</td><td>Maximum allowed least-square guess of constraint multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for resetting equality and inequality multipliers after restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>constraint_violation_norm_type</td><td>OT_STRING</td><td>1-norm</td><td>Norm to be used for the constraint violation in the line search. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>corrector_compl_avrg_red_fact</td><td>OT_REAL</td><td>1.0</td><td>Complementarity tolerance factor for accepting corrector step (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td>none</td><td>The type of corrector steps that should be taken (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>delta</td><td>OT_REAL</td><td>1.0</td><td>Multiplier for constraint violation in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>delta_y_max</td><td>OT_REAL</td><td>1e+12</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dependency_detection_with_rhs</td><td>OT_STRING</td><td>no</td><td>Indicates if the right hand sides of the constraints should be considered during dependency detection (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dependency_detector</td><td>OT_STRING</td><td>none</td><td>Indicates which linear solver should be used to detect linearly dependent equality constraints. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td>none</td><td>Enable derivative checker (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_first_index</td><td>OT_INTEGER</td><td>-2</td><td>Index of first quantity to be checked by derivative checker (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td>1e-08</td><td>Size of the finite difference perturbation in derivative test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td>no</td><td>Indicates whether information for all estimated derivatives should be printed. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td>0.0001</td><td>Threshold for indicating wrong derivative. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td>1e+20</td><td>Threshold for maximal value of primal iterates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td>1.0</td><td>Desired threshold for the dual infeasibility. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>epsilon_c</td><td>OT_REAL</td><td>0.01</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_min</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_penalty</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition for the penalty function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>eta_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td>yes</td><td>Determines if the original objective function should be evaluated at restoration phase trial points. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td>no</td><td>Enable heuristics to quickly detect an infeasible problem. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td>0.001</td><td>Threshold for disabling "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td>100000000.0</td><td>Multiplier threshold for activating "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fast_des_fact</td><td>OT_REAL</td><td>0.1</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fast_step_computation</td><td>OT_STRING</td><td>no</td><td>Indicates if the linear system should be solved quickly. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td>5</td><td>Verbosity level for output file. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_margin_fact</td><td>OT_REAL</td><td>1e-05</td><td>Factor determining width of margin for obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_max_margin</td><td>OT_REAL</td><td>1.0</td><td>Maximum width of margin in obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>filter_reset_trigger</td><td>OT_INTEGER</td><td>5</td><td>Number of iterations that trigger the filter reset. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>findiff_perturbation</td><td>OT_REAL</td><td>1e-07</td><td>Size of the finite difference perturbation for derivative approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td>0.0001</td><td>Size of first x-s perturbation tried. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td>average_compl</td><td>Oracle for the barrier parameter when switching to fixed mode. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td>make_parameter</td><td>Determines how fixed variables should be handled. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_hat</td><td>OT_REAL</td><td>0.04</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the filter margin for the barrier function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_theta</td><td>OT_REAL</td><td>1e-05</td><td>Relaxation factor in the filter margin for the constraint violation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gamma_tilde</td><td>OT_REAL</td><td>4.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>exact</td><td>Indicates what Hessian information is to be used. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>hessian_approximation_space</td><td>OT_STRING</td><td>nonlinear-variables</td><td>Indicates in which subspace the Hessian information is to be approximated. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether the problem is a quadratic problem (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td>yes</td><td>Indicates whether final points should be projected into original bounds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inf_pr_output</td><td>OT_STRING</td><td>original</td><td>Determines what value is printed in the "inf_pr" output column. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all equality constraints are linear (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all inequality constraints are linear (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jacobian_approximation</td><td>OT_STRING</td><td>exact</td><td>Specifies technique to compute constraint Jacobian (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_exponent</td><td>OT_REAL</td><td>0.25</td><td>Exponent for mu in the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td>1e-08</td><td>Size of the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_d</td><td>OT_REAL</td><td>1e-05</td><td>Weight for linear damping term (to handle one-sided bounds). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_sigma</td><td>OT_REAL</td><td>10000000000.0</td><td>Factor limiting the deviation of dual variables from primal estimates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_soc</td><td>OT_REAL</td><td>0.99</td><td>Factor in the sufficient reduction rule for second order correction. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_x_dis</td><td>OT_REAL</td><td>100.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>kappa_y_dis</td><td>OT_REAL</td><td>10000.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>least_square_init_duals</td><td>OT_STRING</td><td>no</td><td>Least square initialization of all dual variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>least_square_init_primal</td><td>OT_STRING</td><td>no</td><td>Least square initialization of the primal variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_aug_solver</td><td>OT_STRING</td><td>sherman-morrison</td><td>Strategy for solving the augmented system for low-rank Hessian. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val</td><td>OT_REAL</td><td>1.0</td><td>Value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_max</td><td>OT_REAL</td><td>100000000.0</td><td>Upper bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_min</td><td>OT_REAL</td><td>1e-08</td><td>Lower bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_initialization</td><td>OT_STRING</td><td>scalar1</td><td>Initialization strategy for the limited memory quasi-Newton approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td>6</td><td>Maximum size of the history for the limited quasi-Newton Hessian approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td>2</td><td>Threshold for successive iterations where update is skipped. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_special_for_resto</td><td>OT_STRING</td><td>no</td><td>Determines if the quasi-Newton updates should be special during the restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>limited_memory_update_type</td><td>OT_STRING</td><td>bfgs</td><td>Quasi-Newton update formula for the limited memory approximation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>line_search_method</td><td>OT_STRING</td><td>filter</td><td>Globalization method used in backtracking line search (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td>yes</td><td>Flag indicating that linear scaling is only done if it seems required. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>ma27</td><td>Linear solver used for step computations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td>mc19</td><td>Method for scaling the linear system. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_ignore_singularity</td><td>OT_STRING</td><td>no</td><td>Enables MA27's ability to solve a linear system even if the matrix is singular. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Real workspace memory for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Integer workspace memory for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td>10.0</td><td>Increment factor for workspace size for MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma27_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma28_pivtol</td><td>OT_REAL</td><td>0.01</td><td>Pivot tolerance for linear solver MA28. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td>yes</td><td>Controls MA57 automatic scaling (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td>16</td><td>Controls block size used by Level 3 BLAS in MA57BD (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td>16</td><td>Node amalgamation parameter (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td>5</td><td>Controls pivot order in MA57 (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td>1.05</td><td>Safety factor for work space memory allocation for the linear solver MA57. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma57_small_pivot_flag</td><td>OT_INTEGER</td><td>0</td><td>If set to 1, then when small entries defined by CNTL(2) are detected they are removed and the corresponding pivots placed at the end of the factorization.  This can be particularly efficient if the matrix is highly rank deficient. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_nemin</td><td>OT_INTEGER</td><td>32</td><td>Node Amalgamation parameter (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_print_level</td><td>OT_INTEGER</td><td>0</td><td>Debug printing level for the linear solver MA86 (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_small</td><td>OT_REAL</td><td>1e-20</td><td>Zero Pivot Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_static</td><td>OT_REAL</td><td>0.0</td><td>Static Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_u</td><td>OT_REAL</td><td>1e-08</td><td>Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>ma86_umax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum Pivoting Threshold (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>magic_steps</td><td>OT_STRING</td><td>no</td><td>Enables magic steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum number of CPU seconds. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_filter_resets</td><td>OT_INTEGER</td><td>5</td><td>Maximal allowed number of filter resets (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td>1e+20</td><td>Maximum value of regularization parameter for handling negative curvature. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>3000</td><td>Maximum number of iterations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_resto_iter</td><td>OT_INTEGER</td><td>3000000</td><td>Maximum number of successive iterations in restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of second order correction trial steps at each iteration. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>max_soft_resto_iters</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterations performed successively in soft restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td>no</td><td>Indicates if we want to do Mehrotra's algorithm. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_alpha_primal</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td>1e-20</td><td>Smallest perturbation of the Hessian block. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td>1</td><td>Minimum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::IpoptInternal</td></tr>
<tr><td>mu_allow_fast_monotone_decrease</td><td>OT_STRING</td><td>yes</td><td>Allow skipping of barrier problem if barrier test is already met. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td>0.1</td><td>Initial value for the barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td>0.2</td><td>Determines linear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td>1000.0</td><td>Factor for initialization of maximum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td>1e-11</td><td>Minimum value for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td>quality-function</td><td>Oracle for a new barrier parameter in the adaptive strategy. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td>monotone</td><td>Update strategy for barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td>1.5</td><td>Determines superlinear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td>0.0</td><td>Desired value of complementarity. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mult_diverg_feasibility_tol</td><td>OT_REAL</td><td>1e-07</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mult_diverg_y_tol</td><td>OT_REAL</td><td>100000000.0</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_dep_tol</td><td>OT_REAL</td><td>-1.0</td><td>Pivot threshold for detection of linearly dependent constraints in MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td>1000</td><td>Percentage increase in the estimated working space for MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td>7</td><td>Controls permuting and scaling in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td>7</td><td>Controls pivot order in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td>1e-06</td><td>Pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td>77</td><td>Controls scaling in MUMPS (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>neg_curv_test_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for heuristic to ignore wrong inertia. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>never_use_fact_cgpen_direction</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the fast Chen-Goldfarb direction (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>never_use_piecewise_penalty_ls</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the piecewise penalty method (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td>-1e+19</td><td>any bound less or equal this value will be considered -inf (i.e. not lower bounded). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_constr_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for constraint function gradient size. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_max_gradient</td><td>OT_REAL</td><td>100.0</td><td>Maximum gradient after NLP scaling. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_method</td><td>OT_STRING</td><td>gradient-based</td><td>Select the technique used for scaling the NLP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_min_value</td><td>OT_REAL</td><td>1e-08</td><td>Minimum value of gradient-based scaling values. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_obj_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for objective function gradient size. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td>1e+19</td><td>any bound greater or this value will be considered +inf (i.e. not upper bounded). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nu_inc</td><td>OT_REAL</td><td>0.0001</td><td>Increment of the penalty parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>nu_init</td><td>OT_REAL</td><td>1e-06</td><td>Initial value of the penalty parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>num_linear_variables</td><td>OT_INTEGER</td><td>0</td><td>Number of linear variables (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>obj_max_inc</td><td>OT_REAL</td><td>5.0</td><td>Determines the upper bound on the acceptable increase of barrier objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>obj_scaling_factor</td><td>OT_REAL</td><td>1.0</td><td>Scaling factor for the objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td>File name of options file (to overwrite default). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td>File name of desired output file (leave unset for no file output). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>pardiso_iter_coarse_size</td><td>OT_INTEGER</td><td>5000</td><td>Maximum Size of Coarse Grid Matrix (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_factor</td><td>OT_REAL</td><td>0.5</td><td>dropping value for incomplete factor (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_schur</td><td>OT_REAL</td><td>0.1</td><td>dropping value for sparsify schur complement factor (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_inverse_norm_factor</td><td>OT_REAL</td><td>5000000.0</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_levels</td><td>OT_INTEGER</td><td>10</td><td>Maximum Size of Grid Levels (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_row_fill</td><td>OT_INTEGER</td><td>10000000</td><td>max fill for each row (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_relative_tol</td><td>OT_REAL</td><td>1e-06</td><td>Relative Residual Convergence (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_iterative</td><td>OT_STRING</td><td>no</td><td>Switch on iterative solver in Pardiso library (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td>complete+2x2</td><td>Matching strategy to be used by Pardiso (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_max_droptol_corrections</td><td>OT_INTEGER</td><td>4</td><td>Maximal number of decreases of drop tolerance during one solve. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_max_iter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of Krylov-Subspace Iteration (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td>0</td><td>Pardiso message level (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td>0</td><td>Enables out-of-core variant of Pardiso (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_redo_symbolic_fact_only_if_inertia_wrong</td><td>OT_STRING</td><td>no</td><td>Toggle for handling case when elements were perturbed by Pardiso. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_repeated_perturbation_means_singular</td><td>OT_STRING</td><td>no</td><td>Interpretation of perturbed elements. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pardiso_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>False</td><td>n/a</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_des_fact</td><td>OT_REAL</td><td>0.2</td><td>a parameter used in penalty parameter computation (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_init_fac</td><td>OT_REAL</td><td>50.0</td><td>a parameter used to choose initial penalty parameterswhen the regularized Newton method is used. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>pen_theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_init_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximal value for the intial penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_init_min</td><td>OT_REAL</td><td>1.0</td><td>Minimal value for the intial penalty parameter for line search(for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_max</td><td>OT_REAL</td><td>1e+30</td><td>Maximal value for the penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_update_compl_tol</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>penalty_update_infeasibility_tol</td><td>OT_REAL</td><td>1e-09</td><td>Threshold for infeasibility in penalty parameter update test. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_always_cd</td><td>OT_STRING</td><td>no</td><td>Active permanent perturbation of constraint linearization. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td>0.333333333333</td><td>Decrease factor for x-s perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td>8.0</td><td>Increase factor for x-s perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td>100.0</td><td>Increase factor for x-s perturbation for very first perturbation. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_infeasi</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_obj</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td>10.0</td><td>Maximal perturbation of an evaluation point. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_info_string</td><td>OT_STRING</td><td>no</td><td>Enables printing of additional info string at end of iteration output. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>5</td><td>Output verbosity level. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td>no</td><td>Switch to print all algorithmic options. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_options_latex_mode</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>print information about execution time</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_timing_statistics</td><td>OT_STRING</td><td>no</td><td>Switch to print timing statistics. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td>no</td><td>Print all options set by the user. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_balancing_term</td><td>OT_STRING</td><td>none</td><td>The balancing term included in the quality function for centrality. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_centrality</td><td>OT_STRING</td><td>none</td><td>The penalty term for centrality that is included in quality function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td>8</td><td>Maximum number of search steps during direct search procedure determining the optimal centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for components of the quality function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_section_qf_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for the golden section search procedure determining the optimal centering parameter (in the function value space). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>quality_function_section_sigma_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for the section search procedure determining the optimal centering parameter (in sigma space). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td>no</td><td>Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility threshold for recomputation of multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>replace_bounds</td><td>OT_STRING</td><td>no</td><td>Indicates if all variable bounds should be replaced by inequality constraints (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td>0.9</td><td>Required reduction of infeasibility before leaving restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_improvement_factor</td><td>OT_REAL</td><td>0.999999999</td><td>Minimal required reduction of residual test ratio in iterative refinement. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_ratio_max</td><td>OT_REAL</td><td>1e-10</td><td>Iterative refinement tolerance (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>residual_ratio_singular</td><td>OT_REAL</td><td>1e-05</td><td>Threshold for declaring linear system singular after failed iterative refinement. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_failure_feasibility_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for primal infeasibility to declare failure of restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_penalty_parameter</td><td>OT_REAL</td><td>1000.0</td><td>Penalty parameter in the restoration phase objective function. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>resto_proximity_weight</td><td>OT_REAL</td><td>1.0</td><td>Weighting factor for the proximity term in restoration phase objective. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.1</td><td>Value in penalty parameter update formula. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_max</td><td>OT_REAL</td><td>100.0</td><td>Scaling threshold for the NLP error. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_phi</td><td>OT_REAL</td><td>2.3</td><td>Exponent for linear barrier function model in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>s_theta</td><td>OT_REAL</td><td>1.1</td><td>Exponent for current constraint violation in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sb</td><td>OT_STRING</td><td>no</td><td> (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sigma_max</td><td>OT_REAL</td><td>100.0</td><td>Maximum value of the centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>sigma_min</td><td>OT_REAL</td><td>1e-06</td><td>Minimum value of the centering parameter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_corr_if_neg_curv</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step in negative curvature iteration (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_corr_in_monotone_mode</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step during monotone barrier parameter mode (unsupported!). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>skip_finalize_solution_call</td><td>OT_STRING</td><td>no</td><td>Indicates if call to NLP::FinalizeSolution after optimization should be suppressed (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial slack to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial slack to bound. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>slack_move</td><td>OT_REAL</td><td>1.81898940355e-12</td><td>Correction size for very small slacks. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td>0.9999</td><td>Required reduction in primal-dual error in the soft restoration phase. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td>no</td><td>Tells algorithm to switch to restoration phase in first iteration. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>suppress_all_output</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tau_min</td><td>OT_REAL</td><td>0.99</td><td>Lower bound on fraction-to-the-boundary parameter tau. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_min</td><td>OT_REAL</td><td>1e-06</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>theta_min_fact</td><td>OT_REAL</td><td>0.0001</td><td>Determines constraint violation threshold in the switching rule. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tiny_step_tol</td><td>OT_REAL</td><td>2.22044604925e-15</td><td>Tolerance for detecting numerically insignificant steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tiny_step_y_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for quitting because of numerically insignificant steps. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-08</td><td>Desired convergence tolerance (relative). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>None</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>casadi::IpoptInternal</td></tr>
<tr><td>vartheta</td><td>OT_REAL</td><td>0.5</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as bound_frac for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_entire_iterate</td><td>OT_STRING</td><td>no</td><td>Tells algorithm whether to use the GetWarmStartIterate method in the NLP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td>no</td><td>Warm-start for initial point (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as mult_bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum initial value for the equality multipliers. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_same_structure</td><td>OT_STRING</td><td>no</td><td>Indicates whether a problem with a structure identical to the previous one is to be solved. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_frac for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_push for the regular initializer. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warm_start_target_mu</td><td>OT_REAL</td><td>0.0</td><td>Unsupported! (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td>10</td><td>Number of shortened iterations that trigger the watchdog. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of watchdog iterations. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_droptol</td><td>OT_REAL</td><td>0.0</td><td>Drop tolerance for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_fillin_limit</td><td>OT_REAL</td><td>0.0</td><td>Fill-in limit for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_iterative</td><td>OT_STRING</td><td>no</td><td>Switches to iterative solver in WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximal number of iterations in iterative WISMP (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_no_pivoting</td><td>OT_STRING</td><td>no</td><td>Use the static pivoting option of WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td>1</td><td>Number of threads to be used in WSMP (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(16) (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option2</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(20) (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td>0.0001</td><td>Pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td>0</td><td>Determines how the matrix is scaled by WSMP. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td>1e-18</td><td>WSMP's singularity threshold. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretent inertia is correct. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
<tr><td>wsmp_write_matrix_iteration</td><td>OT_INTEGER</td><td>-1</td><td>Iteration in which the matrices are written to files. (see IPOPT documentation)</td><td>casadi::IpoptInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::LPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NLPSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::NLPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CplexInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>casadi::CplexInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>casadi::CplexInternal</td></tr>
<tr><td>dep_check</td><td>OT_STRING</td><td>"off"</td><td>Detect redundant constraints. (automatic:-1|off:0|begin:1|end:2|both:3)</td><td>casadi::CplexInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to.</td><td>casadi::CplexInternal</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format.</td><td>casadi::CplexInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>"automatic"</td><td>Determines which CPLEX algorithm to use. (automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover)</td><td>casadi::CplexInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td><td>casadi::CplexInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td><td>casadi::CplexInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td><td>casadi::CplexInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CplexSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>casadi::CplexInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>casadi::CplexInternal</td></tr>
<tr><td>dep_check</td><td>OT_STRING</td><td>"off"</td><td>Detect redundant constraints. (automatic:-1|off:0|begin:1|end:2|both:3)</td><td>casadi::CplexInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to.</td><td>casadi::CplexInternal</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format.</td><td>casadi::CplexInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>"automatic"</td><td>Determines which CPLEX algorithm to use. (automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover)</td><td>casadi::CplexInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td><td>casadi::CplexInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td><td>casadi::CplexInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td><td>casadi::CplexInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DirectCollocationInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::DirectCollocation
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::DirectCollocationInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::OCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SDPSOCPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDPSolver used to solve the SOCPs.</td><td>casadi::SDPSOCPInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDPSOCPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SDPSOCPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SOCPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>sdp_solver</td><td>OT_SDPSOLVER</td><td>GenericType()</td><td>The SDPSolver used to solve the SOCPs.</td><td>casadi::SDPSOCPInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SDPSOCPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CSparseCholeskyInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CSparseCholesky
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ImplicitFixedStepIntegratorInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::ImplicitFixedStepIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegratorInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegratorInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::MXFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::MXFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::SQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::SQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::SQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::SQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::SQPInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>casadi::SQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::SQPInternal</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::SQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::SQPInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::SQPInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::SQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::SQPInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::SQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::SQPMethod
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::SQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::SQPInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX-&gt;SX. Deprecated, use \"expand\" instead.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Use Gauss Newton Hessian approximation</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate a function for calculating the gradient of the objective.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Hessian of the Lagrangian if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Generate an exact Jacobian of the constraints if not supplied.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::SQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::SQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::SQPInternal</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>casadi::SQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::SQPInternal</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::SQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Deprecated option. Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>casadi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::SQPInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::SQPInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::SQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::SQPInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::SQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NLPSolverInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NLPImplicitInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSolver used to solve the implicit system.</td><td>casadi::NLPImplicitInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSolver</td><td>casadi::NLPImplicitInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::NLPImplicitSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSolver used to solve the implicit system.</td><td>casadi::NLPImplicitInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSolver</td><td>casadi::NLPImplicitInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QPOasesInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration. (inactive::all bounds inactive|lower::all bounds active at their lower bound|upper::all bounds active at their upper bound)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QPOasesSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration. (inactive::all bounds inactive|lower::all bounds active at their lower bound|upper::all bounds active at their upper bound)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td><td>casadi::QPOasesInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QCQPQPInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qcqp_solver</td><td>OT_QCQPSOLVER</td><td>GenericType()</td><td>The QCQPSolver used to solve the QPs.</td><td>casadi::QCQPQPInternal</td></tr>
<tr><td>qcqp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QCQPSOlver</td><td>casadi::QCQPQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::QCQPQPSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qcqp_solver</td><td>OT_QCQPSOLVER</td><td>GenericType()</td><td>The QCQPSolver used to solve the QPs.</td><td>casadi::QCQPQPInternal</td></tr>
<tr><td>qcqp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QCQPSOlver</td><td>casadi::QCQPQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CustomFunctionInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::CustomFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
