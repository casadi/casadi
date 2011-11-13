/** \class CasADi::NLPSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NLPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::OOQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listends to print_level 0, 10 and 100</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::OOQPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listends to print_level 0, 10 and 100</td><td>CasADi::Interfaces::OOQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LapackQRDenseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LapackQRDense
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::FlatOCPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>eliminate_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Completely eliminate algebraic states</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>eliminate_dependent</td><td>OT_BOOLEAN</td><td>true</td><td>Eliminate variables that can be expressed as an expression of other variables</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>eliminate_dependents_with_bounds</td><td>OT_BOOLEAN</td><td>true</td><td>Verbose parsing</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>make_explicit</td><td>OT_BOOLEAN</td><td>false</td><td>Make the DAE semi-explicit</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>scale_equations</td><td>OT_BOOLEAN</td><td>false</td><td>Scale the implicit equations so that they get unity order of magnitude</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>scale_variables</td><td>OT_BOOLEAN</td><td>false</td><td>Scale the variables so that they get unity order of magnitude</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>sort_equations</td><td>OT_BOOLEAN</td><td>true</td><td>Sort the dynamic equations</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>true</td><td>Verbose parsing</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::FlatOCP
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>eliminate_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Completely eliminate algebraic states</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>eliminate_dependent</td><td>OT_BOOLEAN</td><td>true</td><td>Eliminate variables that can be expressed as an expression of other variables</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>eliminate_dependents_with_bounds</td><td>OT_BOOLEAN</td><td>true</td><td>Verbose parsing</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>make_explicit</td><td>OT_BOOLEAN</td><td>false</td><td>Make the DAE semi-explicit</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>scale_equations</td><td>OT_BOOLEAN</td><td>false</td><td>Scale the implicit equations so that they get unity order of magnitude</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>scale_variables</td><td>OT_BOOLEAN</td><td>false</td><td>Scale the variables so that they get unity order of magnitude</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>sort_equations</td><td>OT_BOOLEAN</td><td>true</td><td>Sort the dynamic equations</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>true</td><td>Verbose parsing</td><td>CasADi::OptimalControl::FlatOCPInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LapackLUDenseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LapackLUDense
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LapackLUDenseInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CollocationIntegratorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau or legendre)</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CollocationIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau or legendre)</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::KinsolInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (\"none\" or \"linesearch\")</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::KinsolSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strategy (\"none\" or \"linesearch\")</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::Sundials::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShootingInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShooting
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::OptimalControl::MultipleShootingInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SimulatorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Simulator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::QPSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::QPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::QPOasesInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Defines the amount of text output during QP solution, see Section 5.7.: \"none\", \"low\", \"medium\" or \"high\"</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::QPOasesSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Defines the amount of text output during QP solution, see Section 5.7.: \"none\", \"low\", \"medium\" or \"high\"</td><td>CasADi::Interfaces::QPOasesInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>inplace</td><td>OT_BOOLEAN</td><td>false</td><td>Evaluate with inplace operations (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>live_variables</td><td>OT_BOOLEAN</td><td>false</td><td>Reuse variables in the work vector</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SXFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>inplace</td><td>OT_BOOLEAN</td><td>false</td><td>Evaluate with inplace operations (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>live_variables</td><td>OT_BOOLEAN</td><td>false</td><td>Reuse variables in the work vector</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IpoptInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
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
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
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
/** \class CasADi::AcadoOCPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_time</td><td>OT_REAL</td><td>0.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoOCP
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td>None</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_time</td><td>OT_REAL</td><td>0.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>CasADi::Sundials::CVodesInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>CasADi::Sundials::CVodesInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>"bdf" or "adams"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>"newton" or "functional"</td><td>CasADi::Sundials::CVodesInternal</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoIntegratorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>num_algebraic</td><td>OT_INTEGER</td><td>0</td><td>Number of algebraic states</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>num_grid_points</td><td>OT_INTEGER</td><td>2</td><td>Number of uniformly distributed grid points for obtaining the solution, does not influence the integration steps</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>time_dependence</td><td>OT_BOOLEAN</td><td>true</td><td>Explicit depencency of time in the DAE</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>num_algebraic</td><td>OT_INTEGER</td><td>0</td><td>Number of algebraic states</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>num_grid_points</td><td>OT_INTEGER</td><td>2</td><td>Number of uniformly distributed grid points for obtaining the solution, does not influence the integration steps</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>time_dependence</td><td>OT_BOOLEAN</td><td>true</td><td>Explicit depencency of time in the DAE</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::JacobianInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default" (default meaning checking both)</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Jacobian
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"default"</td><td>"forward", "adjoint" or "default" (default meaning checking both)</td><td>CasADi::JacobianInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::MuscodInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>acc</td><td>OT_REAL</td><td>1e-6</td><td>accuracy of NLP solution</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>backupfile</td><td>OT_STRING</td><td>"restart.bin"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>bflag</td><td>OT_INTEGER</td><td>-1</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>cflag</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>datfile</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>eflag</td><td>OT_INTEGER</td><td>0</td><td>use two gradient evaluations per SQP iteration</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>itol</td><td>OT_REAL</td><td>1e-7</td><td>initial integration tolerance</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>levmar</td><td>OT_REAL</td><td>0.0</td><td>Levenberg Marquardt regularization of hessian</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>logfile</td><td>OT_STRING</td><td>"log.txt"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>mf</td><td>OT_INTEGER</td><td>0</td><td>maximum # of successive ``fixed'' iterations</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>mi</td><td>OT_INTEGER</td><td>100</td><td>maximum total # of SQP iterations</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"muscod_problem"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>nhtopy</td><td>OT_INTEGER</td><td>0</td><td># of homotopy steps</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>ppath</td><td>OT_STRING</td><td>"PAR"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>resfile</td><td>OT_STRING</td><td>"muscod_results.txt"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>restartfile</td><td>OT_STRING</td><td>"restart.bin"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>rfac</td><td>OT_INTEGER</td><td>0.0</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>sf</td><td>OT_INTEGER</td><td>0</td><td>index of first SQP iteration using fixed discretization</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>sflag</td><td>OT_INTEGER</td><td>0</td><td>stop after each SQP iteration</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-7</td><td>integration tolerance</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>wflag</td><td>OT_INTEGER</td><td>0</td><td>warm/cool start using final data of previous run</td><td>CasADi::MuscodInternal</td></tr>
</table>
*/
/** \class CasADi::MuscodInterface
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>acc</td><td>OT_REAL</td><td>1e-6</td><td>accuracy of NLP solution</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>backupfile</td><td>OT_STRING</td><td>"restart.bin"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>bflag</td><td>OT_INTEGER</td><td>-1</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>cflag</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>datfile</td><td>OT_STRING</td><td></td><td>None</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>eflag</td><td>OT_INTEGER</td><td>0</td><td>use two gradient evaluations per SQP iteration</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>itol</td><td>OT_REAL</td><td>1e-7</td><td>initial integration tolerance</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>levmar</td><td>OT_REAL</td><td>0.0</td><td>Levenberg Marquardt regularization of hessian</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>logfile</td><td>OT_STRING</td><td>"log.txt"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>mf</td><td>OT_INTEGER</td><td>0</td><td>maximum # of successive ``fixed'' iterations</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>mi</td><td>OT_INTEGER</td><td>100</td><td>maximum total # of SQP iterations</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"muscod_problem"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>nhtopy</td><td>OT_INTEGER</td><td>0</td><td># of homotopy steps</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>ppath</td><td>OT_STRING</td><td>"PAR"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>resfile</td><td>OT_STRING</td><td>"muscod_results.txt"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>restartfile</td><td>OT_STRING</td><td>"restart.bin"</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>rfac</td><td>OT_INTEGER</td><td>0.0</td><td></td><td>CasADi::MuscodInternal</td></tr>
<tr><td>sf</td><td>OT_INTEGER</td><td>0</td><td>index of first SQP iteration using fixed discretization</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>sflag</td><td>OT_INTEGER</td><td>0</td><td>stop after each SQP iteration</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-7</td><td>integration tolerance</td><td>CasADi::MuscodInternal</td></tr>
<tr><td>wflag</td><td>OT_INTEGER</td><td>0</td><td>warm/cool start using final data of previous run</td><td>CasADi::MuscodInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::IpoptQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::IpoptQPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>false</td><td>Specify true if you can guarantee that H will always be positive definite</td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::QPSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td></td><td></td><td>CasADi::Interfaces::IpoptQPInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lifted</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>optimizer</td><td>OT_STRING</td><td>"sqp"</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lifted</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>optimizer</td><td>OT_STRING</td><td>"sqp"</td><td></td><td>CasADi::Interfaces::LiftoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IntegratorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Integrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CplexInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td>Option()</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_cplex_problem"</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CplexSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td>Option()</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::CplexInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_cplex_problem"</td><td></td><td>CasADi::CplexInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>objsense</td><td>OT_INTEGER</td><td>CPX_MIN</td><td>optimization sense (CPX_MIN or CPX_MAX)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::XFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::XFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLUInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>conditionnumber</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>diagpivotthresh</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>iterrefine</td><td>OT_STRING</td><td>"norefine"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pivotgrowth</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>printstat</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>rowperm</td><td>OT_STRING</td><td>"largediag"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symmetricmode</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SuperLU
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>colperm</td><td>OT_STRING</td><td>"colamd"</td><td>Specifies how to permute the columns of the matrix for sparsity preservation.</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>conditionnumber</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>diagpivotthresh</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>equil</td><td>OT_BOOLEAN</td><td>true</td><td>Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>iterrefine</td><td>OT_STRING</td><td>"norefine"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pivotgrowth</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>printstat</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>rowperm</td><td>OT_STRING</td><td>"largediag"</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>symmetricmode</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_work</td><td>OT_BOOLEAN</td><td>false</td><td>keep work in memory</td><td>CasADi::SuperLUInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::FXInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::FX
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ExternalFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ExternalFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>use IDACalcIC to get consistent initial conditions</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>false</td><td>use IDACalcIC to get consistent initial conditions</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>first requested time as a fraction of the time interval</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>maximim step size</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>supress algebraic variables in the error testing</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>asens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>asens_upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>use IDACalcIC to get consistent initial conditions</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>false</td><td>use IDACalcIC to get consistent initial conditions</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>first requested time as a fraction of the time interval</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td></td><td>absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td>None</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_INTEGER</td><td>false</td><td>include the forward sensitivities in all error controls</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td></td><td>relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td></td><td>scaling factor for the components if finite differences is used</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td></td><td>specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>type of interpolation for the adjoint sensitivities ("hermite" or "polynomial")</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>is_differential</td><td>OT_INTEGERVECTOR</td><td></td><td>None</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>"gmres", "bcgstab", "tfqmr"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>maximim step size</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>supress algebraic variables in the error testing</td><td>CasADi::Sundials::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ParallelizerInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>save_corrected_input</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Parallelizer
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>serial, openmp or mpi</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>save_corrected_input</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eta</td><td>OT_REAL   </td><td>0.0001</td><td>Linesearch parameter: See Nocedal 3.4</td><td>CasADi::SQPInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of linesearch iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>mu_safety</td><td>OT_REAL   </td><td>1.1</td><td>Safety factor for linesearch mu</td><td>CasADi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SQPInternal</td></tr>
<tr><td>rho</td><td>OT_REAL   </td><td>0.5</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sigma</td><td>OT_REAL   </td><td>1.0</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tau</td><td>OT_REAL   </td><td>0.2</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>toldx</td><td>OT_REAL   </td><td>1e-12</td><td>Stopping criterion for the stepsize</td><td>CasADi::SQPInternal</td></tr>
<tr><td>tolgl</td><td>OT_REAL   </td><td>1e-12</td><td>Stopping criterion for the Lagrangian gradient</td><td>CasADi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SQPMethod
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eta</td><td>OT_REAL   </td><td>0.0001</td><td>Linesearch parameter: See Nocedal 3.4</td><td>CasADi::SQPInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Generate an exact Hessian of the Lagrangian</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of linesearch iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>mu_safety</td><td>OT_REAL   </td><td>1.1</td><td>Safety factor for linesearch mu</td><td>CasADi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SQPInternal</td></tr>
<tr><td>rho</td><td>OT_REAL   </td><td>0.5</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sigma</td><td>OT_REAL   </td><td>1.0</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tau</td><td>OT_REAL   </td><td>0.2</td><td>Linesearch parameter</td><td>CasADi::SQPInternal</td></tr>
<tr><td>toldx</td><td>OT_REAL   </td><td>1e-12</td><td>Stopping criterion for the stepsize</td><td>CasADi::SQPInternal</td></tr>
<tr><td>tolgl</td><td>OT_REAL   </td><td>1e-12</td><td>Stopping criterion for the Lagrangian gradient</td><td>CasADi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::MXFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::MXFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>topological_sorting</td><td>OT_STRING</td><td>"breadth-first"</td><td>Topological sorting algorithm: \"depth-first\" or \"breadth-first\" search</td><td>CasADi::XFunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OCPSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OCPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ImplicitFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ImplicitFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CFunctionInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CFunction
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OptionsFunctionalityNode
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
</table>
*/
/** \class CasADi::OptionsFunctionality
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
</table>
*/
/** \class CasADi::LinearSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LinearSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::GSL::GslInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::GSL::GslIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>absolute tolerence  for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"dense"</td><td>"dense", "banded" or "iterative"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_creator</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver creator function</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td>lower band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>maximum krylov subspace size</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>maximum number of steps</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nrhs</td><td>OT_INTEGER</td><td>1</td><td>number of right hand sides</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>"none", "left", "right", "both"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>should the quadratures affect the step size control</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>relative tolerence for the IVP solution</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>"simultaneous" or "staggered"</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>number of steps between two consecutive checkpoints</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>false</td><td>Stop the integrator at the end of the interval</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>start of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>end of the integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td>upper band-width of banded jacobians</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::CSparseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Interfaces::CSparse
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians: \"forward\" (only forward mode) \"reverse\" (only adjoint mode) or \"automatic\" (a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>jac_for_sens</td><td>OT_BOOLEAN</td><td>false</td><td>Create the a Jacobian function and use this to calculate forward sensitivities</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
