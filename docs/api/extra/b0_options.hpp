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
/** \addtogroup plugin_LinearSolver_csparsecholesky
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::CleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::CleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::CleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_CleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::CleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::CleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::CleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CollocationIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::CollocationIntegrator</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::CollocationIntegrator</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Integrator_collocation
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CondensingIndefDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
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
/** \addtogroup plugin_DpleSolver_condensing
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
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
<tr><td>integrator</td><td>OT_STRING</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::ControlSimulatorInternal</td></tr>
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
<tr><td>integrator</td><td>OT_STRING</td><td>GenericType()</td><td>An integrator creator function</td><td>casadi::ControlSimulatorInternal</td></tr>
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
/** \class casadi::CplexInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>casadi::CplexInterface</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>casadi::CplexInterface</td></tr>
<tr><td>dep_check</td><td>OT_STRING</td><td>"off"</td><td>Detect redundant constraints. (automatic:-1|off:0|begin:1|end:2|both:3)</td><td>casadi::CplexInterface</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to.</td><td>casadi::CplexInterface</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format.</td><td>casadi::CplexInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>"automatic"</td><td>Determines which CPLEX algorithm to use. (automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover)</td><td>casadi::CplexInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td><td>casadi::CplexInterface</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td><td>casadi::CplexInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td><td>casadi::CplexInterface</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_QpSolver_cplex
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td></tr>
<tr><td>dep_check</td><td>OT_STRING</td><td>"off"</td><td>Detect redundant constraints. (automatic:-1|off:0|begin:1|end:2|both:3)</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to.</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format.</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>"automatic"</td><td>Determines which CPLEX algorithm to use. (automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover)</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CsparseInterface
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
/** \addtogroup plugin_LinearSolver_csparse
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
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
/// \cond INTERNAL
/** \class casadi::CvodesInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>casadi::CvodesInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>casadi::CvodesInterface</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>casadi::CvodesInterface</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>casadi::FunctionInternal<br />casadi::CvodesInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>casadi::CvodesInterface</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Integrator_cvodes
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_DleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DleToLrDle
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_LrDleSolver_dle
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
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
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
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
/** \addtogroup general_DpleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
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
/** \class casadi::DpleToDle
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_DleSolver_dple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DpleToLrDple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_LrDpleSolver_dple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::DsdpInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_loglevel</td><td>OT_INTEGER</td><td>0</td><td>An integer that specifies how much logging is done on stdout.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_penalty</td><td>OT_REAL</td><td>1e5</td><td>Penality parameter lambda. Must exceed the trace of Y. This parameter heavily influences the ability of DSDP to treat linear equalities. The DSDP standard default (1e8) will make a problem with linear equality return unusable solutions.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_printlevel</td><td>OT_INTEGER</td><td>1</td><td>A printlevel of zero will disable all output. Another number indicates how often a line is printed.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_reuse</td><td>OT_INTEGER</td><td>4</td><td>Maximum on the number of times the Schur complement matrix is reused</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_rho</td><td>OT_REAL</td><td>4.0</td><td>Potential parameter. Must be &gt;=1</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_use_penalty</td><td>OT_BOOLEAN</td><td>true</td><td>Modifies the algorithm to use a penality gamma on r.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>_zbar</td><td>OT_REAL</td><td>1e10</td><td>Initial upper bound on the objective of the dual problem.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td><td>casadi::DsdpInterface</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td><td>casadi::DsdpInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>infinity</td><td>OT_REAL</td><td>1e30</td><td>Treat numbers higher than this as infinity</td><td>casadi::DsdpInterface</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td><td>casadi::DsdpInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td><td>casadi::DsdpInterface</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td><td>casadi::DsdpInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_SdpSolver_dsdp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>_loglevel</td><td>OT_INTEGER</td><td>0</td><td>An integer that specifies how much logging is done on stdout.</td></tr>
<tr><td>_penalty</td><td>OT_REAL</td><td>1e5</td><td>Penality parameter lambda. Must exceed the trace of Y. This parameter heavily influences the ability of DSDP to treat linear equalities. The DSDP standard default (1e8) will make a problem with linear equality return unusable solutions.</td></tr>
<tr><td>_printlevel</td><td>OT_INTEGER</td><td>1</td><td>A printlevel of zero will disable all output. Another number indicates how often a line is printed.</td></tr>
<tr><td>_reuse</td><td>OT_INTEGER</td><td>4</td><td>Maximum on the number of times the Schur complement matrix is reused</td></tr>
<tr><td>_rho</td><td>OT_REAL</td><td>4.0</td><td>Potential parameter. Must be &gt;=1</td></tr>
<tr><td>_use_penalty</td><td>OT_BOOLEAN</td><td>true</td><td>Modifies the algorithm to use a penality gamma on r.</td></tr>
<tr><td>_zbar</td><td>OT_REAL</td><td>1e10</td><td>Initial upper bound on the objective of the dual problem.</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td></tr>
<tr><td>infinity</td><td>OT_REAL</td><td>1e30</td><td>Treat numbers higher than this as infinity</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td></tr>
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
/** \class casadi::FixedSmithLrDleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iter</td><td>OT_INTEGER</td><td>100</td><td>Number of Smith iterations</td><td>casadi::FixedSmithLrDleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_LrDleSolver_fixed_smith
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>iter</td><td>OT_INTEGER</td><td>100</td><td>Number of Smith iterations</td></tr>
</table>
*/
/// \cond INTERNAL
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
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
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
/** \addtogroup general_HomotopyNlpSolver
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
/** \class casadi::IdasInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>casadi::IdasInterface</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>casadi::IdasInterface</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>casadi::IdasInterface</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>casadi::IdasInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>casadi::IdasInterface</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td><td>casadi::IdasInterface</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::IdasInterface</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>casadi::IdasInterface</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td><td>casadi::IdasInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup)</td><td>casadi::FunctionInternal<br />casadi::IdasInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Suppress algebraic variables in the error testing</td><td>casadi::IdasInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Integrator_idas
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Suppress algebraic variables in the error testing</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td></tr>
</table>
*/
/// \cond INTERNAL
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
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
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
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_ImplicitFunction
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
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
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
/** \addtogroup general_Integrator
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
/** \class casadi::IpoptInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::IpoptInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>casadi::IpoptInterface</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>casadi::IpoptInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_ipopt
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>false</td><td></td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::KinsolInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>casadi::KinsolInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>casadi::KinsolInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>casadi::KinsolInterface</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>casadi::KinsolInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.</td><td>casadi::KinsolInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>casadi::FunctionInternal<br />casadi::KinsolInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::KinsolInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (none|linesearch)</td><td>casadi::KinsolInterface</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>casadi::KinsolInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_ImplicitFunction_kinsol
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (none|linesearch)</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::KnitroInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td><td>casadi::KnitroInterface</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>1e-4</td><td>Feasible mode tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute feasible tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>1</td><td>Allow infeasible iterations</td><td>casadi::KnitroInterface</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td><td>casadi::KnitroInterface</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td><td>casadi::KnitroInterface</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td><td>casadi::KnitroInterface</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td><td>casadi::KnitroInterface</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td><td>casadi::KnitroInterface</td></tr>
<tr><td>LpSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LpSolver</td><td>casadi::KnitroInterface</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td><td>casadi::KnitroInterface</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td><td>casadi::KnitroInterface</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td><td>casadi::KnitroInterface</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e20</td><td>Maximum objective value</td><td>casadi::KnitroInterface</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td><td>casadi::KnitroInterface</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td><td>casadi::KnitroInterface</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td><td>casadi::KnitroInterface</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td><td>casadi::KnitroInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>casadi::KnitroInterface</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::KnitroInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_knitro
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>1e-4</td><td>Feasible mode tolerance</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute feasible tolerance</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>1</td><td>Allow infeasible iterations</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td></tr>
<tr><td>LpSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LpSolver</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e20</td><td>Maximum objective value</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LapackLuDense
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>casadi::LapackLuDense</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::LapackLuDense</td></tr>
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
/** \addtogroup plugin_LinearSolver_lapacklu
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LapackQrDense
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
/** \addtogroup plugin_LinearSolver_lapackqr
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LiftingIndefDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>form</td><td>OT_STRING</td><td>"A"</td><td>The form of the lifting (A:0|B:1)</td><td>casadi::LiftingIndefDpleInternal</td></tr>
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
/** \addtogroup plugin_DpleSolver_lifting
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>form</td><td>OT_STRING</td><td>"A"</td><td>The form of the lifting (A:0|B:1)</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LiftingLrDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>form</td><td>OT_STRING</td><td>"A"</td><td>The form of the lifting (A:0|B:1)</td><td>casadi::LiftingLrDpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_LrDpleSolver_lifting
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>form</td><td>OT_STRING</td><td>"A"</td><td>The form of the lifting (A:0|B:1)</td></tr>
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
/** \addtogroup general_LinearSolver
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
/** \class casadi::LpSolverInternal
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
/** \addtogroup general_LpSolver
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
/** \class casadi::LpToQp
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
/** \addtogroup plugin_LpSolver_qp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LrDleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_LrDleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LrDleToDle
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_DleSolver_lrdle
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LrDpleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_LrDpleSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::LrDpleToDple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>const_dim</td><td>OT_BOOLEAN</td><td>true</td><td>Assume constant dimension of P</td><td>casadi::DpleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DpleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
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
/** \addtogroup plugin_DpleSolver_lrdple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
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
/** \class casadi::Newton
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>casadi::Newton</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>casadi::Newton</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>casadi::Newton</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>casadi::FunctionInternal<br />casadi::Newton</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_iteration</td><td>OT_BOOLEAN</td><td>false</td><td>Print information about each iteration</td><td>casadi::Newton</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_ImplicitFunction_newton
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td></tr>
<tr><td>print_iteration</td><td>OT_BOOLEAN</td><td>false</td><td>Print information about each iteration</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::NlpSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_NlpSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
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
/** \class casadi::OldCollocationIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the implicit solver</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_STRING</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>casadi::OldCollocationIntegrator</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Integrator_oldcollocation
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the implicit solver</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td></tr>
<tr><td>startup_integrator</td><td>OT_STRING</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::OoqpInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>casadi::OoqpInterface</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>casadi::OoqpInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>casadi::OoqpInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_QpSolver_ooqp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td></tr>
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
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>psd_num_zero</td><td>OT_REAL</td><td>1e-12</td><td>Numerical zero used in Periodic Schur decomposition with slicot.This option is needed when your systems has Floquet multiplierszero or close to zero</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_DpleSolver_slicot
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td></tr>
<tr><td>psd_num_zero</td><td>OT_REAL</td><td>1e-12</td><td>Numerical zero used in Periodic Schur decomposition with slicot.This option is needed when your systems has Floquet multiplierszero or close to zero</td></tr>
</table>
*/
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
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DpleInternal</td></tr>
<tr><td>psd_num_zero</td><td>OT_REAL</td><td>1e-12</td><td>Numerical zero used in Periodic Schur decomposition with slicot.This option is needed when your systems has Floquet multiplierszero or close to zero</td><td>casadi::PsdIndefDpleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QcqpSolverInternal
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
/** \addtogroup general_QcqpSolver
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
/** \class casadi::QcqpToSocp
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
/** \addtogroup plugin_QcqpSolver_socp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpSolverInternal
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
/** \addtogroup general_QpSolver
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
/** \class casadi::QpToImplicit
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
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::ImplicitFunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_ImplicitFunction_nlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpToNlp
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
/** \addtogroup plugin_QpSolver_nlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpToQcqp
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
/** \addtogroup plugin_QpSolver_qcqp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpoasesInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration. (inactive::all bounds inactive|lower::all bounds active at their lower bound|upper::all bounds active at their upper bound)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_QpSolver_qpoases
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration. (inactive::all bounds inactive|lower::all bounds active at their lower bound|upper::all bounds active at their upper bound)</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::RkIntegrator
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
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Integrator_rk
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td></tr>
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
/** \class casadi::Scpgen
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::Scpgen</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::Scpgen</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::Scpgen</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::Scpgen</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>gauss-newton|exact</td><td>casadi::Scpgen</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::Scpgen</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::Scpgen</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>casadi::Scpgen</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::Scpgen</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>casadi::Scpgen</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::Scpgen</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>casadi::Scpgen</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::Scpgen</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::Scpgen</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>casadi::Scpgen</td></tr>
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::Scpgen</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::Scpgen</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>casadi::Scpgen</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>casadi::Scpgen</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_scpgen
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>gauss-newton|exact</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td></tr>
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SdpSolverInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_SdpSolver
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (m x m).</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SdpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SdqpSolverInternal
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
<tr><td>sdp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The SdqpSolver used to solve the SDPs.</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_SdqpSolver
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
<tr><td>sdp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The SdqpSolver used to solve the SDPs.</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SdqpToSdp
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
<tr><td>sdp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The SdqpSolver used to solve the SDPs.</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>sdp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the SDPSOlver</td><td>casadi::SdqpSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_SdqpSolver_sdp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimpleHomotopyNlp
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
<tr><td>nlp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The NLP solver to be used by the Homotopy solver</td><td>casadi::SimpleHomotopyNlp</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Homotopy solver</td><td>casadi::SimpleHomotopyNlp</td></tr>
<tr><td>num_steps</td><td>OT_INTEGER</td><td>10</td><td>Take this many steps to go from tau=0 to tau=1.</td><td>casadi::SimpleHomotopyNlp</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_HomotopyNlpSolver_simple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>nlp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The NLP solver to be used by the Homotopy solver</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Homotopy solver</td></tr>
<tr><td>num_steps</td><td>OT_INTEGER</td><td>10</td><td>Take this many steps to go from tau=0 to tau=1.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimpleIndefCleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::CleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::CleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::SimpleIndefCleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::SimpleIndefCleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::CleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_CleSolver_simple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimpleIndefDleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compressed_solve</td><td>OT_BOOLEAN</td><td>true</td><td>When a system with sparse rhs arises, compress toa smaller system with dense rhs.</td><td>casadi::SimpleIndefDleInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::DleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::SimpleIndefDleInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::SimpleIndefDleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::DleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_DleSolver_simple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>compressed_solve</td><td>OT_BOOLEAN</td><td>true</td><td>When a system with sparse rhs arises, compress toa smaller system with dense rhs.</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td></tr>
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
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::DpleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::SimpleIndefDpleInternal</td></tr>
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
/** \addtogroup plugin_DpleSolver_simple
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td></tr>
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
/** \class casadi::SmithLrDleInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eps_unstable</td><td>OT_REAL</td><td>1e-4</td><td>A margin for unstability detection</td><td>casadi::LrDleInternal</td></tr>
<tr><td>error_unstable</td><td>OT_BOOLEAN</td><td>false</td><td>Throw an exception when it is detected that Product(A_i, i=N..1) has eigenvalues greater than 1-eps_unstable</td><td>casadi::LrDleInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of iterations for the algorithm</td><td>casadi::SmithLrDleInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pos_def</td><td>OT_BOOLEAN</td><td>false</td><td>Assume P positive definite</td><td>casadi::LrDleInternal</td></tr>
<tr><td>print_iteration</td><td>OT_BOOLEAN</td><td>false</td><td>Print information about each iteration</td><td>casadi::SmithLrDleInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-12</td><td>Tolerance for satisfying the Lyapunov equation.</td><td>casadi::SmithLrDleInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_LrDleSolver_smith
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of iterations for the algorithm</td></tr>
<tr><td>print_iteration</td><td>OT_BOOLEAN</td><td>false</td><td>Print information about each iteration</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-12</td><td>Tolerance for satisfying the Lyapunov equation.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SnoptInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>_iprint</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::SnoptInterface</td></tr>
<tr><td>_isumm</td><td>OT_INTEGER</td><td>6</td><td></td><td>casadi::SnoptInterface</td></tr>
<tr><td>_start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Warm)</td><td>casadi::SnoptInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td><td>casadi::SnoptInterface</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_nlp|setup_nlp)</td><td>casadi::FunctionInternal<br />casadi::SnoptInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>casadi::SnoptInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_snopt
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>_iprint</td><td>OT_INTEGER</td><td>0</td><td></td></tr>
<tr><td>_isumm</td><td>OT_INTEGER</td><td>6</td><td></td></tr>
<tr><td>_start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Warm)</td></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SocpSolverInternal
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
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_SocpSolver
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
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SocpToSdp
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
<tr><td>ni</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Provide the size of each SOC constraint. Must sum up to N.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>print_problem</td><td>OT_BOOLEAN</td><td>false</td><td>Print out problem statement for debugging.</td><td>casadi::SocpSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_SocpSolver_sdp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SqicInterface
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
/** \addtogroup plugin_QpSolver_sqic
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Sqpmethod
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::Sqpmethod</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::Sqpmethod</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::Sqpmethod</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::Sqpmethod</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>casadi::Sqpmethod</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::Sqpmethod</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx|bfgs)</td><td>casadi::FunctionInternal<br />casadi::Sqpmethod</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::Sqpmethod</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::Sqpmethod</td></tr>
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::Sqpmethod</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::Sqpmethod</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::Sqpmethod</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::Sqpmethod</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_sqpmethod
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td></tr>
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedQpSolverInternal
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
/** \addtogroup general_StabilizedQpSolver
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
/** \class casadi::StabilizedQpToQp
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
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver used to solve the stabilized QPs.</td><td>casadi::StabilizedQpToQp</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver instance</td><td>casadi::StabilizedQpToQp</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_StabilizedQpSolver_qp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver used to solve the stabilized QPs.</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver instance</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedSqicInterface
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
/** \addtogroup plugin_StabilizedQpSolver_sqic
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::StabilizedSqp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>TReta1</td><td>OT_REAL</td><td>0.8</td><td>Required predicted / actual decrease for TR increase</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>TReta2</td><td>OT_REAL</td><td>0.2</td><td>Required predicted / actual decrease for TR decrease</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>alphaMin</td><td>OT_REAL</td><td>1e-3</td><td>Used to check whether to increase rho.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.5</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>0.001</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>dvMax0</td><td>OT_REAL</td><td>100</td><td>Parameter used to defined the max step length.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>eps_active</td><td>OT_REAL</td><td>1e-6</td><td>Threshold for the epsilon-active set.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gamma1</td><td>OT_REAL</td><td>2.</td><td>Trust region increase parameter</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>gamma2</td><td>OT_REAL</td><td>1.</td><td>Trust region update parameter</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>gamma3</td><td>OT_REAL</td><td>1.</td><td>Trust region decrease parameter</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>20</td><td>Maximum number of linesearch iterations</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>max_time</td><td>OT_REAL</td><td>1e12</td><td>Timeout</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::StabilizedSqp</td></tr>
<tr><td>muR0</td><td>OT_REAL</td><td>1e-4</td><td>Initial choice of regularization parameter</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>nu</td><td>OT_REAL</td><td>1</td><td>Parameter for primal-dual augmented Lagrangian.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>phiWeight</td><td>OT_REAL</td><td>1e-5</td><td>Weight used in pseudo-filter.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>stabilized_qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The Stabilized QP solver to be used by the SQP method</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>stabilized_qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Stabilized QP solver</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>tau0</td><td>OT_REAL</td><td>1e-2</td><td>Initial parameter for the merit function optimality threshold.</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for dual infeasability</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for primal infeasibility</td><td>casadi::StabilizedSqp</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>yEinitial</td><td>OT_STRING</td><td>"simple"</td><td>Initial multiplier. Simple (all zero) or least (LSQ).</td><td>casadi::StabilizedSqp</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_stabilizedsqp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>TReta1</td><td>OT_REAL</td><td>0.8</td><td>Required predicted / actual decrease for TR increase</td></tr>
<tr><td>TReta2</td><td>OT_REAL</td><td>0.2</td><td>Required predicted / actual decrease for TR decrease</td></tr>
<tr><td>alphaMin</td><td>OT_REAL</td><td>1e-3</td><td>Used to check whether to increase rho.</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.5</td><td>Line-search parameter, restoration factor of stepsize</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>0.001</td><td>Armijo condition, coefficient of decrease in merit</td></tr>
<tr><td>dvMax0</td><td>OT_REAL</td><td>100</td><td>Parameter used to defined the max step length.</td></tr>
<tr><td>eps_active</td><td>OT_REAL</td><td>1e-6</td><td>Threshold for the epsilon-active set.</td></tr>
<tr><td>gamma1</td><td>OT_REAL</td><td>2.</td><td>Trust region increase parameter</td></tr>
<tr><td>gamma2</td><td>OT_REAL</td><td>1.</td><td>Trust region update parameter</td></tr>
<tr><td>gamma3</td><td>OT_REAL</td><td>1.</td><td>Trust region decrease parameter</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>20</td><td>Maximum number of linesearch iterations</td></tr>
<tr><td>max_time</td><td>OT_REAL</td><td>1e12</td><td>Timeout</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td></tr>
<tr><td>muR0</td><td>OT_REAL</td><td>1e-4</td><td>Initial choice of regularization parameter</td></tr>
<tr><td>nu</td><td>OT_REAL</td><td>1</td><td>Parameter for primal-dual augmented Lagrangian.</td></tr>
<tr><td>phiWeight</td><td>OT_REAL</td><td>1e-5</td><td>Weight used in pseudo-filter.</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td></tr>
<tr><td>stabilized_qp_solver</td><td>OT_STRING</td><td>GenericType()</td><td>The Stabilized QP solver to be used by the SQP method</td></tr>
<tr><td>stabilized_qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the Stabilized QP solver</td></tr>
<tr><td>tau0</td><td>OT_REAL</td><td>1e-2</td><td>Initial parameter for the merit function optimality threshold.</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for dual infeasability</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-5</td><td>Stopping criterion for primal infeasibility</td></tr>
<tr><td>yEinitial</td><td>OT_STRING</td><td>"simple"</td><td>Initial multiplier. Simple (all zero) or least (LSQ).</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SundialsInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>expand_augmented</td><td>OT_BOOLEAN</td><td>true</td><td>If DAE callback functions are SXFunction, have augmented DAE callback function also be SXFunction.</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SymbolicQr
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SymbolicQr</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SymbolicQr</td></tr>
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
/** \addtogroup plugin_LinearSolver_symbolicqr
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::TinyXmlInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_XmlFile_tinyxml
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::WorhpInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>derivative_generator</td><td>OT_DERIVATIVEGENERATOR</td><td>GenericType()</td><td>Function that returns a derivative function given a number of forward and reverse directional derivative, overrides internal routines. Check documentation of DerivativeGenerator.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_CALLBACK</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::NlpSolverInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />Monitor functions (eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::WorhpInterface</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>worhp_p_.qp.ipBarrier</td><td>IP barrier parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>worhp_p_.qp.ipComTol</td><td>IP complementarity tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>worhp_p_.qp.ipFracBound</td><td>IP fraction-to-the-boundary parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the IP method. (LAPACK::0|MA57: only available if provided by the user:1|SuperLU::2|PARDISO: only available if provided by the user, subject to license availability:3|MUMPS: currently Linux platforms only:5|WSMP: subject to license availability:6|MA86: experimental, only available if provided by the user:7|MA97:experimental, only available if provided by the user:8)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>worhp_p_.qp.ipMinAlpha</td><td>IP line search minimum step size.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxDiv</td><td>The relaxation term is divided by this value if successful.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMax</td><td>Maximum relaxation value.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMin</td><td>Mimimum relaxation value.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMult</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>worhp_p_.qp.ipResTol</td><td>IP residuals tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.ipTryRelax</td><td>Enable relaxation strategy when encountering an error.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.lsItMaxIter</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the iterative linear solver. (none:Deactivate; use a direct linear solver.:0|CGNR::1|CGNE::2|CGS::3|BiCGSTAB::4)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select preconditioner for the iterative linear solver. (none:No preconditioner.:0|static:Static preconditioner (KKT-matrix with constant lower-right block).:1|full:Full KKT-matrix.:2)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.lsRefineMaxIter</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.lsScale</td><td>Enables scaling on linear solver level.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>worhp_p_.qp.lsTol</td><td>Tolerance for the linear solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.lsTrySimple</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.maxIter</td><td>Imposes an upper limit on the number of minor solver iterations,  i.e. for the quadratic subproblem solver.If the limit is reached before convergence, WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>GenericType()</td><td>Select the solution method used by the QP solver. (ip:Interior-Point method.:1|nsn:Nonsmooth-Newton method.:2|automatic: Prefer IP and fall back to NSN on error.:12)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>worhp_p_.qp.nsnBeta</td><td>NSN stepsize decrease factor.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.nsnGradStep</td><td>Enable gradient steps in the NSN method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>worhp_p_.qp.nsnKKT</td><td>NSN KKT tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the NSN method. (SuperLU::2|MA48: only available if provided by the user:4)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>worhp_p_.qp.nsnMinAlpha</td><td>NSN line search minimum step size.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>worhp_p_.qp.nsnSigma</td><td>NSN line search slope parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Controls the amount of QP solver output. (none:No output.:0|warn:Print warnings and errors.:1|iterations:Print iterations.:2)</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.scaleIntern</td><td>Enable scaling on QP level.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.strict</td><td>Use strict termination criteria in IP method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::NlpSolverInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_NlpSolver_worhp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>worhp_p_.qp.ipBarrier</td><td>IP barrier parameter.</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>worhp_p_.qp.ipComTol</td><td>IP complementarity tolerance.</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>worhp_p_.qp.ipFracBound</td><td>IP fraction-to-the-boundary parameter.</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the IP method. (LAPACK::0|MA57: only available if provided by the user:1|SuperLU::2|PARDISO: only available if provided by the user, subject to license availability:3|MUMPS: currently Linux platforms only:5|WSMP: subject to license availability:6|MA86: experimental, only available if provided by the user:7|MA97:experimental, only available if provided by the user:8)</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>worhp_p_.qp.ipMinAlpha</td><td>IP line search minimum step size.</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxDiv</td><td>The relaxation term is divided by this value if successful.</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMax</td><td>Maximum relaxation value.</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMin</td><td>Mimimum relaxation value.</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>worhp_p_.qp.ipRelaxMult</td><td>The relaxation term is multiplied by this value if unsuccessful.</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>worhp_p_.qp.ipResTol</td><td>IP residuals tolerance.</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.ipTryRelax</td><td>Enable relaxation strategy when encountering an error.</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.lsItMaxIter</td><td>Maximum number of iterations of the iterative linear solvers.</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the iterative linear solver. (none:Deactivate; use a direct linear solver.:0|CGNR::1|CGNE::2|CGS::3|BiCGSTAB::4)</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select preconditioner for the iterative linear solver. (none:No preconditioner.:0|static:Static preconditioner (KKT-matrix with constant lower-right block).:1|full:Full KKT-matrix.:2)</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.lsRefineMaxIter</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.lsScale</td><td>Enables scaling on linear solver level.</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>worhp_p_.qp.lsTol</td><td>Tolerance for the linear solver.</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.lsTrySimple</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>worhp_p_.qp.maxIter</td><td>Imposes an upper limit on the number of minor solver iterations,  i.e. for the quadratic subproblem solver.If the limit is reached before convergence, WORHP will activate QP recovery strategies to prevent a solver breakdown.</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>GenericType()</td><td>Select the solution method used by the QP solver. (ip:Interior-Point method.:1|nsn:Nonsmooth-Newton method.:2|automatic: Prefer IP and fall back to NSN on error.:12)</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>worhp_p_.qp.nsnBeta</td><td>NSN stepsize decrease factor.</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.nsnGradStep</td><td>Enable gradient steps in the NSN method.</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>worhp_p_.qp.nsnKKT</td><td>NSN KKT tolerance.</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the NSN method. (SuperLU::2|MA48: only available if provided by the user:4)</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>worhp_p_.qp.nsnMinAlpha</td><td>NSN line search minimum step size.</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>worhp_p_.qp.nsnSigma</td><td>NSN line search slope parameter.</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Controls the amount of QP solver output. (none:No output.:0|warn:Print warnings and errors.:1|iterations:Print iterations.:2)</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.scaleIntern</td><td>Enable scaling on QP level.</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>worhp_p_.qp.strict</td><td>Use strict termination criteria in IP method.</td></tr>
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
/** \class casadi::XmlFileInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_XmlFile
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
