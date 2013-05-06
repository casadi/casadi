/** \class CasADi::NLPSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::NLPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::QPOasesInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>boundRelaxation</td><td>OT_REAL</td><td>double</td><td>Initial relaxation of bounds to start homotopy  and initial value for far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>boundTolerance</td><td>OT_REAL</td><td>double</td><td>If upper and lower bounds differ less than this tolerance, they are regarded equal, i.e. as  equality constraint.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableCholeskyRefactorisation</td><td>OT_INTEGER</td><td>int</td><td>Specifies the frequency of a full re-factorisation of projected Hessian matrix: 0: turns them off,  1: uses them at each iteration etc.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableDriftCorrection</td><td>OT_INTEGER</td><td>int</td><td>Specifies the frequency of drift corrections: 0: turns them off.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsDen</td><td>OT_REAL</td><td>double</td><td>Denominator tolerance for ratio tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsFlipping</td><td>OT_REAL</td><td>double</td><td>Tolerance of squared Cholesky diagonal factor  which triggers flipping bound.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsIterRef</td><td>OT_REAL</td><td>double</td><td>Early termination tolerance for iterative  refinement.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsLITests</td><td>OT_REAL</td><td>double</td><td>Tolerance for linear independence tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsNZCTests</td><td>OT_REAL</td><td>double</td><td>Tolerance for nonzero curvature tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsNum</td><td>OT_REAL</td><td>double</td><td>Numerator tolerance for ratio tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsRegularisation</td><td>OT_REAL</td><td>double</td><td>Scaling factor of identity matrix used for  Hessian regularisation.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>finalRamping</td><td>OT_REAL</td><td>double</td><td>Final value for ramping strategy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>growFarBounds</td><td>OT_REAL</td><td>double</td><td>Factor to grow far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialFarBounds</td><td>OT_REAL</td><td>double</td><td>Initial size for far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialRamping</td><td>OT_REAL</td><td>double</td><td>Start value for ramping strategy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration: \"inactive\": all bounds inactive, \"lower\": all bounds active at their lower bound, \"upper\": all bounds active at their upper bound.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxDualJump</td><td>OT_REAL</td><td>double</td><td>Maximum allowed jump in dual variables in  linear independence tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>maxPrimalJump</td><td>OT_REAL</td><td>double</td><td>Maximum allowed jump in primal variables in  nonzero curvature tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>numRefinementSteps</td><td>OT_INTEGER</td><td>int</td><td>Maximum number of iterative refinement steps.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>numRegularisationSteps</td><td>OT_INTEGER</td><td>int</td><td>Maximum number of successive regularisation steps.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>terminationTolerance</td><td>OT_REAL</td><td>double</td><td>Relative termination tolerance to stop homotopy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::QPOasesSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>GenericType()</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>boundRelaxation</td><td>OT_REAL</td><td>double</td><td>Initial relaxation of bounds to start homotopy  and initial value for far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>boundTolerance</td><td>OT_REAL</td><td>double</td><td>If upper and lower bounds differ less than this tolerance, they are regarded equal, i.e. as  equality constraint.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableCholeskyRefactorisation</td><td>OT_INTEGER</td><td>int</td><td>Specifies the frequency of a full re-factorisation of projected Hessian matrix: 0: turns them off,  1: uses them at each iteration etc.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableDriftCorrection</td><td>OT_INTEGER</td><td>int</td><td>Specifies the frequency of drift corrections: 0: turns them off.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables the use of  flipping bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables nonzero curvature  tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables ramping.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>BooleanType_to_bool</td><td>Enables automatic  Hessian regularisation.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsDen</td><td>OT_REAL</td><td>double</td><td>Denominator tolerance for ratio tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsFlipping</td><td>OT_REAL</td><td>double</td><td>Tolerance of squared Cholesky diagonal factor  which triggers flipping bound.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsIterRef</td><td>OT_REAL</td><td>double</td><td>Early termination tolerance for iterative  refinement.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsLITests</td><td>OT_REAL</td><td>double</td><td>Tolerance for linear independence tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsNZCTests</td><td>OT_REAL</td><td>double</td><td>Tolerance for nonzero curvature tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsNum</td><td>OT_REAL</td><td>double</td><td>Numerator tolerance for ratio tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>epsRegularisation</td><td>OT_REAL</td><td>double</td><td>Scaling factor of identity matrix used for  Hessian regularisation.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>finalRamping</td><td>OT_REAL</td><td>double</td><td>Final value for ramping strategy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>growFarBounds</td><td>OT_REAL</td><td>double</td><td>Factor to grow far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialFarBounds</td><td>OT_REAL</td><td>double</td><td>Initial size for far bounds.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialRamping</td><td>OT_REAL</td><td>double</td><td>Start value for ramping strategy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>SubjectToStatus_to_string</td><td>Initial status of bounds at first iteration: \"inactive\": all bounds inactive, \"lower\": all bounds active at their lower bound, \"upper\": all bounds active at their upper bound.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxDualJump</td><td>OT_REAL</td><td>double</td><td>Maximum allowed jump in dual variables in  linear independence tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>maxPrimalJump</td><td>OT_REAL</td><td>double</td><td>Maximum allowed jump in primal variables in  nonzero curvature tests.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>GenericType()</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>numRefinementSteps</td><td>OT_INTEGER</td><td>int</td><td>Maximum number of iterative refinement steps.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>numRegularisationSteps</td><td>OT_INTEGER</td><td>int</td><td>Maximum number of successive regularisation steps.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>PrintLevel_to_string</td><td>Defines the amount of text output during QP solution, see Section 5.7 (none|low|medium|high)</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>terminationTolerance</td><td>OT_REAL</td><td>double</td><td>Relative termination tolerance to stop homotopy.</td><td>CasADi::QPOasesInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CSparseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CSparse
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
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
<tr><td>datfile</td><td>OT_STRING</td><td></td><td></td><td>CasADi::MuscodInternal</td></tr>
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
<tr><td>datfile</td><td>OT_STRING</td><td></td><td></td><td>CasADi::MuscodInternal</td></tr>
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
/** \class CasADi::KinsolInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>CasADi::FXInternal<br />CasADi::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strateg (none|linesearch)</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::KinsolSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>CasADi::FXInternal<br />CasADi::KinsolInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>strategy</td><td>OT_STRING</td><td>"none"</td><td>Globalization strateg (none|linesearch)</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>u_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::KinsolInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>precondition an iterative solver</td><td>CasADi::KinsolInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DerivativeInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::Derivative
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>The linear solver to be used by the IP method</td><td>CasADi::IPInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::IPMethod
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>The linear solver to be used by the IP method</td><td>CasADi::IPInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::IPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::NLPImplicitInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSolver used to solve the implicit system.</td><td>CasADi::NLPImplicitInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSolver</td><td>CasADi::NLPImplicitInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NLPImplicitSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSolver used to solve the implicit system.</td><td>CasADi::NLPImplicitInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSolver</td><td>CasADi::NLPImplicitInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::WorhpInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>AcceptTolFeas</td><td>OT_REAL</td><td>worhp_p.AcceptTolFeas</td><td>Tolerance for acceptable feasibility</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AcceptTolOpti</td><td>OT_REAL</td><td>worhp_p.AcceptTolOpti</td><td>Tolerance for acceptable optimality</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AlphaMinConst</td><td>OT_BOOLEAN</td><td>worhp_p.AlphaMinConst</td><td>Use a constant lower bound on Armijo stepsize in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoBeta</td><td>OT_REAL</td><td>worhp_p.ArmijoBeta</td><td>Trial stepsize decrease factor for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMaxAlpha</td><td>OT_REAL</td><td>worhp_p.ArmijoMaxAlpha</td><td>Initial alpha for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlpha</td><td>OT_REAL</td><td>worhp_p.ArmijoMinAlpha</td><td>Lower bound on alpha for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlphaRec</td><td>OT_REAL</td><td>worhp_p.ArmijoMinAlphaRec</td><td>Lower bound on alpha for Armijo rule during recovery</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoSigma</td><td>OT_REAL</td><td>worhp_p.ArmijoSigma</td><td>Scale factor for linearised descent check in Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AutoQPRecovery</td><td>OT_BOOLEAN</td><td>worhp_p.AutoQPRecovery</td><td>Enable automatic QP recovery</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSmaxblockSize</td><td>OT_INTEGER</td><td>worhp_p.BFGSmaxblockSize</td><td>Maximum BFGS block size (depends on BFGS method)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSmethod</td><td>OT_INTEGER</td><td>worhp_p.BFGSmethod</td><td>Choose BFGS method (dense, block, sparse)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSminblockSize</td><td>OT_INTEGER</td><td>worhp_p.BFGSminblockSize</td><td>Minimum BFGS block size (depends on BFGS method)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSrestart</td><td>OT_INTEGER</td><td>worhp_p.BFGSrestart</td><td>Restart BFGS update after this many iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BettsFactor</td><td>OT_REAL</td><td>worhp_p.BettsFactor</td><td>Update factor for Betts' Hessian regularisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BettsPoint</td><td>OT_REAL</td><td>worhp_p.BettsPoint</td><td>Smallest eigenvalue of the regularised Hessian</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BoundTolFac</td><td>OT_REAL</td><td>worhp_p.BoundTolFac</td><td>Factor in determining active constraints by KKT</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckFJ</td><td>OT_REAL</td><td>worhp_p.CheckFJ</td><td>Upper bound used by Fritz-John heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureDF</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureDF</td><td>Enable structural checking of DF</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureDG</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureDG</td><td>Enable structural checking of DG</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureHM</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureHM</td><td>Enable structural checking of HM</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepBettsSum</td><td>OT_REAL</td><td>worhp_p.CorStepBettsSum</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepConStop</td><td>OT_REAL</td><td>worhp_p.CorStepConStop</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepConvio</td><td>OT_REAL</td><td>worhp_p.CorStepConvio</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMaxIter</td><td>OT_INTEGER</td><td>worhp_p.CorStepMaxIter</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMethod</td><td>OT_INTEGER</td><td>worhp_p.CorStepMethod</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMode</td><td>OT_INTEGER</td><td>worhp_p.CorStepMode</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepPFactor</td><td>OT_REAL</td><td>worhp_p.CorStepPFactor</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepPMax</td><td>OT_REAL</td><td>worhp_p.CorStepPMax</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepRecoveryDX</td><td>OT_BOOLEAN</td><td>worhp_p.CorStepRecoveryDX</td><td>Enable structural checking of HM</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvBCond</td><td>OT_REAL</td><td>worhp_p.CurvBCond</td><td>Block BFGS curvature condition bound</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvBFac</td><td>OT_REAL</td><td>worhp_p.CurvBFac</td><td>Block BFGS curvature condition regularisation factor</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvCond</td><td>OT_REAL</td><td>worhp_p.CurvCond</td><td>BFGS Curvature condition bound</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvFac</td><td>OT_REAL</td><td>worhp_p.CurvFac</td><td>BFGS curvature condition regularisation factor</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CutLength</td><td>OT_REAL</td><td>worhp_p.CutLength</td><td>Scaling factor for Cut recovery strategy</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>DebugMarker06</td><td>OT_INTEGER</td><td>worhp_p.DebugMarker06</td><td>Debug marker, only needed for ASTOS integration</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FGtogether</td><td>OT_BOOLEAN</td><td>worhp_p.FGtogether</td><td>F and G cannot be evaluated separately</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FJandND</td><td>OT_BOOLEAN</td><td>worhp_p.FJandND</td><td>Enable Fritz-John and non-differentiable check heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleDual</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleDual</td><td>Activate dual feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleInit</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleInit</td><td>Activate initial feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleInitTol</td><td>OT_REAL</td><td>worhp_p.FeasibleInitTol</td><td>Feasibility tolerance for no-objective feasible mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleOnly</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleOnly</td><td>Activate feasible-only mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FidifEps</td><td>OT_REAL</td><td>worhp_p.FidifEps</td><td>Finite difference perturbation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FidifHM</td><td>OT_BOOLEAN</td><td>worhp_p.FidifHM</td><td>Approximate Hessian by finite differences (otherwise BFGS)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterBisecAlpha</td><td>OT_BOOLEAN</td><td>worhp_p.FilterBisecAlpha</td><td>Filter heuristic to save Armijo iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterGammaCV</td><td>OT_REAL</td><td>worhp_p.FilterGammaCV</td><td>Constraint violation decrease factor in Filter acceptance check</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterGammaF</td><td>OT_REAL</td><td>worhp_p.FilterGammaF</td><td>Objective decrease factor in Filter acceptance check</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterIntersecAlpha</td><td>OT_BOOLEAN</td><td>worhp_p.FilterIntersecAlpha</td><td>Filter heuristic to save Armijo iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FirstDifCentral</td><td>OT_BOOLEAN</td><td>worhp_p.FirstDifCentral</td><td>Use central finite difference quotient for first derivatives</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FocusOnFeas</td><td>OT_BOOLEAN</td><td>worhp_p.FocusOnFeas</td><td>Enable Focus-on-Feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FocusOnFeasFactor</td><td>OT_REAL</td><td>worhp_p.FocusOnFeasFactor</td><td>Factor in Focus-on-Feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>GammaAlpha</td><td>OT_REAL</td><td>worhp_p.GammaAlpha</td><td>Safety factor for alphamin calculation by Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>GroupMethod</td><td>OT_INTEGER</td><td>worhp_p.GroupMethod</td><td>Select method to determine graph colouring groups</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IgnoreFilterCrit</td><td>OT_BOOLEAN</td><td>worhp_p.IgnoreFilterCrit</td><td>Activate accelerating heuristics for Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncBettsTau</td><td>OT_REAL</td><td>worhp_p.IncBettsTau</td><td>Increase factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncBettsTauMore</td><td>OT_REAL</td><td>worhp_p.IncBettsTauMore</td><td>Larger increase factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncreaseIWS</td><td>OT_REAL</td><td>worhp_p.IncreaseIWS</td><td>Increase factor for estimated integer workspace requirement</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncreaseRWS</td><td>OT_REAL</td><td>worhp_p.IncreaseRWS</td><td>Increase factor for estimated real workspace requirement</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Infty</td><td>OT_REAL</td><td>worhp_p.Infty</td><td>Upper bound for numbers to be regarded as finite</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>InftyUnbounded</td><td>OT_REAL</td><td>worhp_p.InftyUnbounded</td><td>Tolerance for unboundedness detection heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>InitialLMest</td><td>OT_BOOLEAN</td><td>worhp_p.InitialLMest</td><td>Enable initial Lagrange multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>KeepAcceptableSol</td><td>OT_BOOLEAN</td><td>worhp_p.KeepAcceptableSol</td><td>Save acceptable solutions as fallback</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LMestQPipComTol</td><td>OT_REAL</td><td>worhp_p.LMestQPipComTol</td><td>IP complementarity tolerance in initial multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LMestQPipResTol</td><td>OT_REAL</td><td>worhp_p.LMestQPipResTol</td><td>IP residual tolerance in initial multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LinMult</td><td>OT_BOOLEAN</td><td>worhp_p.LinMult</td><td>Control Lagrange multiplier update</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LogLevel</td><td>OT_INTEGER</td><td>worhp_p.LogLevel</td><td>Enable XML logfiles and writing interval</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LogResult</td><td>OT_INTEGER</td><td>worhp_p.LogResult</td><td>Enable XML result logging and detail level</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaF</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaF</td><td>Lowpass-filter update factor for objective values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaG</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaG</td><td>Lowpass-filter update factor for constraint values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaMerit</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaMerit</td><td>Lowpass-filter update factor for merit function values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassFilter</td><td>OT_BOOLEAN</td><td>worhp_p.LowPassFilter</td><td>Enable lowpass-filter termination criterion</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Ma57PivotThresh</td><td>OT_REAL</td><td>worhp_p.Ma57PivotThresh</td><td>Pivoting tolerance for MA57 = CNTL(1)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MatrixCC</td><td>OT_BOOLEAN</td><td>worhp_p.MatrixCC</td><td>Not to be included into a parameter file!</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxCalls</td><td>OT_INTEGER</td><td>worhp_p.MaxCalls</td><td>Upper bound to Reverse Communication calls</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxForce</td><td>OT_INTEGER</td><td>worhp_p.MaxForce</td><td>Maximum number of Force recovery strategy steps</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxGPart</td><td>OT_INTEGER</td><td>worhp_p.MaxGPart</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxIter</td><td>OT_INTEGER</td><td>worhp_p.MaxIter</td><td>Upper bound on major iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxLScounter</td><td>OT_INTEGER</td><td>worhp_p.MaxLScounter</td><td>Control activation of Filter acceleration heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxNorm</td><td>OT_BOOLEAN</td><td>worhp_p.MaxNorm</td><td>Select max-norm instead of 1-norm in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MeritFunction</td><td>OT_INTEGER</td><td>worhp_p.MeritFunction</td><td>Select merit function and penalty update [0, 3..5]</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MeritGradTol</td><td>OT_REAL</td><td>worhp_p.MeritGradTol</td><td>Threshold of meritfunction gradient for increasing Hessian regularisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MinBettsTau</td><td>OT_REAL</td><td>worhp_p.MinBettsTau</td><td>Lower bound for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MoreRelax</td><td>OT_BOOLEAN</td><td>worhp_p.MoreRelax</td><td>Introduce one relaxation variable for every constraint</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>NLPmethod</td><td>OT_INTEGER</td><td>worhp_p.NLPmethod</td><td>Select (1) Meritfunction or (3) Filter globalisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>NLPprint</td><td>OT_INTEGER</td><td>worhp_p.NLPprint</td><td>NLP print level [-1..4]</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PairMethod</td><td>OT_INTEGER</td><td>worhp_p.PairMethod</td><td>Select method to determine graph colouring pairgroups</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsBar</td><td>OT_REAL</td><td>worhp_p.PenUpdEpsBar</td><td>Penalty update parameter factor for MeritFunction = 3</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKFac</td><td>OT_REAL</td><td>worhp_p.PenUpdEpsKFac</td><td>Penalty update parameter factor for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKSequence</td><td>OT_INTEGER</td><td>worhp_p.PenUpdEpsKSequence</td><td>Penalty update parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxDeltaK</td><td>OT_REAL</td><td>worhp_p.PenUpdMaxDeltaK</td><td>Max penalty for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxFac</td><td>OT_REAL</td><td>worhp_p.PenUpdMaxFac</td><td>Max factor for increasing penalty for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdRBar</td><td>OT_REAL</td><td>worhp_p.PenUpdRBar</td><td>Penalty update parameter for MeritFunction = 3</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PrecisionF</td><td>OT_REAL</td><td>worhp_p.PrecisionF</td><td>(currently unused) Relative precision of objective</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PrecisionG</td><td>OT_REAL</td><td>worhp_p.PrecisionG</td><td>(currently unused) Relative precision of constraints</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>QPscaleParam</td><td>OT_REAL</td><td>worhp_p.QPscaleParam</td><td>(currently unused)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>QuadraticProblem</td><td>OT_BOOLEAN</td><td>worhp_p.QuadraticProblem</td><td>Not to be included into a parameter file!</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ReduceBettsTau</td><td>OT_REAL</td><td>worhp_p.ReduceBettsTau</td><td>Decrease factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RegStrategy</td><td>OT_INTEGER</td><td>worhp_p.RegStrategy</td><td>Select Hessian regularisation strategy in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ReinitFilter</td><td>OT_BOOLEAN</td><td>worhp_p.ReinitFilter</td><td>Enables Filter-reinitialisation accelerating heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxMaxDelta</td><td>OT_REAL</td><td>worhp_p.RelaxMaxDelta</td><td>Upper bound for accepting the constraint relaxation variable</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxMaxPen</td><td>OT_REAL</td><td>worhp_p.RelaxMaxPen</td><td>Upper bound on the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxRho</td><td>OT_REAL</td><td>worhp_p.RelaxRho</td><td>Update factor for the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxStart</td><td>OT_REAL</td><td>worhp_p.RelaxStart</td><td>Initial value of the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RestUntilFeas</td><td>OT_BOOLEAN</td><td>worhp_p.RestUntilFeas</td><td>Do restoration until a feasible solution is found</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleConIter</td><td>OT_BOOLEAN</td><td>worhp_p.ScaleConIter</td><td>Scale constraints in every iteration</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleFacObj</td><td>OT_REAL</td><td>worhp_p.ScaleFacObj</td><td>Value to scale large objective functions to</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleFacQP</td><td>OT_REAL</td><td>worhp_p.ScaleFacQP</td><td>Upper bound on resulting matrix norm for QP scaling</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledFD</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledFD</td><td>Use a scaled perturbation for finite differences</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledKKT</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledKKT</td><td>Scale KKT conditions</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledObj</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledObj</td><td>Scale the objective function</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledQP</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledQP</td><td>Scale some matrices handed to the QP</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>StartBettsTau</td><td>OT_REAL</td><td>worhp_p.StartBettsTau</td><td>Initial value for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingDelta</td><td>OT_REAL</td><td>worhp_p.SwitchingDelta</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingSCV</td><td>OT_REAL</td><td>worhp_p.SwitchingSCV</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingSF</td><td>OT_REAL</td><td>worhp_p.SwitchingSF</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TakeQPSol</td><td>OT_BOOLEAN</td><td>worhp_p.TakeQPSol</td><td>Evaluate QP search direction regardless of convergence</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Timeout</td><td>OT_REAL</td><td>worhp_p.Timeout</td><td>Timeout in seconds</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolComp</td><td>OT_REAL</td><td>worhp_p.TolComp</td><td>Complementarity tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolFeas</td><td>OT_REAL</td><td>worhp_p.TolFeas</td><td>Feasibility tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolOpti</td><td>OT_REAL</td><td>worhp_p.TolOpti</td><td>Optimality tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolWeakActive</td><td>OT_REAL</td><td>worhp_p.TolWeakActive</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBig</td><td>OT_BOOLEAN</td><td>worhp_p.TooBig</td><td>Enable too-big termination heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBigCV</td><td>OT_REAL</td><td>worhp_p.TooBigCV</td><td>Upper bound on constraint violation for too-big heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBigKKT</td><td>OT_REAL</td><td>worhp_p.TooBigKKT</td><td>Upper bound on KKT values for too-big heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserDF</td><td>OT_BOOLEAN</td><td>worhp_p.UserDF</td><td>Objective gradient values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserDG</td><td>OT_BOOLEAN</td><td>worhp_p.UserDG</td><td>Jacobian values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>worhp_p.UserHM</td><td>Hessian values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserHMstructure</td><td>OT_INTEGER</td><td>worhp_p.UserHMstructure</td><td>Enable automatic Hessian structure generation or checking</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>WeakActiveSet</td><td>OT_BOOLEAN</td><td>worhp_p.WeakActiveSet</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eps</td><td>OT_REAL</td><td>worhp_p.eps</td><td>Machine epsilon</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>initialised</td><td>OT_BOOLEAN</td><td>worhp_p.initialised</td><td>Automatically added initialisation flag.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>CasADi::FXInternal<br />CasADi::WorhpInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>worhp_p.qp.ipBarrier</td><td>IP barrier parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>worhp_p.qp.ipComTol</td><td>IP complementarity tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>worhp_p.qp.ipFracBound</td><td>IP fraction-to-the-boundary parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the IP method. (LAPACK::0|MA57: only available if provided by the user:1|SuperLU::2|PARDISO: only available if provided by the user, subject to license availability:3|MUMPS: currently Linux platforms only:5|WSMP: subject to license availability:6|MA86: experimental, only available if provided by the user:7|MA97:experimental, only available if provided by the user:8)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>worhp_p.qp.ipMinAlpha</td><td>IP line search minimum step size.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxDiv</td><td>The relaxation term is divided by this value if successful.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMax</td><td>Maximum relaxation value.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMin</td><td>Mimimum relaxation value.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMult</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>worhp_p.qp.ipResTol</td><td>IP residuals tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>worhp_p.qp.ipTryRelax</td><td>Enable relaxation strategy when encountering an error.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.lsItMaxIter</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the iterative linear solver. (none:Deactivate; use a direct linear solver.:0|CGNR::1|CGNE::2|CGS::3|BiCGSTAB::4)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select preconditioner for the iterative linear solver. (none:No preconditioner.:0|static:Static preconditioner (KKT-matrix with constant lower-right block).:1|full:Full KKT-matrix.:2)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.lsRefineMaxIter</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>worhp_p.qp.lsScale</td><td>Enables scaling on linear solver level.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>worhp_p.qp.lsTol</td><td>Tolerance for the linear solver.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>worhp_p.qp.lsTrySimple</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported. Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.maxIter</td><td>Imposes an upper limit on the number of minor solver iterations, i.e. for thequadratic subproblem solver. If the limit is reached before convergence,WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>GenericType()</td><td>Select the solution method used by the QP solver. (ip:Interior-Point method.:1|nsn:Nonsmooth-Newton method.:2|automatic: Prefer IP and fall back to NSN on error.:12)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>worhp_p.qp.nsnBeta</td><td>NSN stepsize decrease factor.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>worhp_p.qp.nsnGradStep</td><td>Enable gradient steps in the NSN method.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>worhp_p.qp.nsnKKT</td><td>NSN KKT tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the NSN method. (SuperLU::2|MA48: only available if provided by the user:4)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>worhp_p.qp.nsnMinAlpha</td><td>NSN line search minimum step size.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>worhp_p.qp.nsnSigma</td><td>NSN line search slope parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Controls the amount of QP solver output. (none:No output.:0|warn:Print warnings and errors.:1|iterations:Print iterations.:2)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>worhp_p.qp.scaleIntern</td><td>Enable scaling on QP level.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>worhp_p.qp.strict</td><td>Use strict termination criteria in IP method.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::WorhpSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>AcceptTolFeas</td><td>OT_REAL</td><td>worhp_p.AcceptTolFeas</td><td>Tolerance for acceptable feasibility</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AcceptTolOpti</td><td>OT_REAL</td><td>worhp_p.AcceptTolOpti</td><td>Tolerance for acceptable optimality</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AlphaMinConst</td><td>OT_BOOLEAN</td><td>worhp_p.AlphaMinConst</td><td>Use a constant lower bound on Armijo stepsize in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoBeta</td><td>OT_REAL</td><td>worhp_p.ArmijoBeta</td><td>Trial stepsize decrease factor for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMaxAlpha</td><td>OT_REAL</td><td>worhp_p.ArmijoMaxAlpha</td><td>Initial alpha for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlpha</td><td>OT_REAL</td><td>worhp_p.ArmijoMinAlpha</td><td>Lower bound on alpha for Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoMinAlphaRec</td><td>OT_REAL</td><td>worhp_p.ArmijoMinAlphaRec</td><td>Lower bound on alpha for Armijo rule during recovery</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ArmijoSigma</td><td>OT_REAL</td><td>worhp_p.ArmijoSigma</td><td>Scale factor for linearised descent check in Armijo rule</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>AutoQPRecovery</td><td>OT_BOOLEAN</td><td>worhp_p.AutoQPRecovery</td><td>Enable automatic QP recovery</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSmaxblockSize</td><td>OT_INTEGER</td><td>worhp_p.BFGSmaxblockSize</td><td>Maximum BFGS block size (depends on BFGS method)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSmethod</td><td>OT_INTEGER</td><td>worhp_p.BFGSmethod</td><td>Choose BFGS method (dense, block, sparse)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSminblockSize</td><td>OT_INTEGER</td><td>worhp_p.BFGSminblockSize</td><td>Minimum BFGS block size (depends on BFGS method)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BFGSrestart</td><td>OT_INTEGER</td><td>worhp_p.BFGSrestart</td><td>Restart BFGS update after this many iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BettsFactor</td><td>OT_REAL</td><td>worhp_p.BettsFactor</td><td>Update factor for Betts' Hessian regularisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BettsPoint</td><td>OT_REAL</td><td>worhp_p.BettsPoint</td><td>Smallest eigenvalue of the regularised Hessian</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>BoundTolFac</td><td>OT_REAL</td><td>worhp_p.BoundTolFac</td><td>Factor in determining active constraints by KKT</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckFJ</td><td>OT_REAL</td><td>worhp_p.CheckFJ</td><td>Upper bound used by Fritz-John heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureDF</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureDF</td><td>Enable structural checking of DF</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureDG</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureDG</td><td>Enable structural checking of DG</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CheckStructureHM</td><td>OT_BOOLEAN</td><td>worhp_p.CheckStructureHM</td><td>Enable structural checking of HM</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepBettsSum</td><td>OT_REAL</td><td>worhp_p.CorStepBettsSum</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepConStop</td><td>OT_REAL</td><td>worhp_p.CorStepConStop</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepConvio</td><td>OT_REAL</td><td>worhp_p.CorStepConvio</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMaxIter</td><td>OT_INTEGER</td><td>worhp_p.CorStepMaxIter</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMethod</td><td>OT_INTEGER</td><td>worhp_p.CorStepMethod</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepMode</td><td>OT_INTEGER</td><td>worhp_p.CorStepMode</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepPFactor</td><td>OT_REAL</td><td>worhp_p.CorStepPFactor</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepPMax</td><td>OT_REAL</td><td>worhp_p.CorStepPMax</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CorStepRecoveryDX</td><td>OT_BOOLEAN</td><td>worhp_p.CorStepRecoveryDX</td><td>Enable structural checking of HM</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvBCond</td><td>OT_REAL</td><td>worhp_p.CurvBCond</td><td>Block BFGS curvature condition bound</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvBFac</td><td>OT_REAL</td><td>worhp_p.CurvBFac</td><td>Block BFGS curvature condition regularisation factor</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvCond</td><td>OT_REAL</td><td>worhp_p.CurvCond</td><td>BFGS Curvature condition bound</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CurvFac</td><td>OT_REAL</td><td>worhp_p.CurvFac</td><td>BFGS curvature condition regularisation factor</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>CutLength</td><td>OT_REAL</td><td>worhp_p.CutLength</td><td>Scaling factor for Cut recovery strategy</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>DebugMarker06</td><td>OT_INTEGER</td><td>worhp_p.DebugMarker06</td><td>Debug marker, only needed for ASTOS integration</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FGtogether</td><td>OT_BOOLEAN</td><td>worhp_p.FGtogether</td><td>F and G cannot be evaluated separately</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FJandND</td><td>OT_BOOLEAN</td><td>worhp_p.FJandND</td><td>Enable Fritz-John and non-differentiable check heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleDual</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleDual</td><td>Activate dual feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleInit</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleInit</td><td>Activate initial feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleInitTol</td><td>OT_REAL</td><td>worhp_p.FeasibleInitTol</td><td>Feasibility tolerance for no-objective feasible mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FeasibleOnly</td><td>OT_BOOLEAN</td><td>worhp_p.FeasibleOnly</td><td>Activate feasible-only mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FidifEps</td><td>OT_REAL</td><td>worhp_p.FidifEps</td><td>Finite difference perturbation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FidifHM</td><td>OT_BOOLEAN</td><td>worhp_p.FidifHM</td><td>Approximate Hessian by finite differences (otherwise BFGS)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterBisecAlpha</td><td>OT_BOOLEAN</td><td>worhp_p.FilterBisecAlpha</td><td>Filter heuristic to save Armijo iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterGammaCV</td><td>OT_REAL</td><td>worhp_p.FilterGammaCV</td><td>Constraint violation decrease factor in Filter acceptance check</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterGammaF</td><td>OT_REAL</td><td>worhp_p.FilterGammaF</td><td>Objective decrease factor in Filter acceptance check</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FilterIntersecAlpha</td><td>OT_BOOLEAN</td><td>worhp_p.FilterIntersecAlpha</td><td>Filter heuristic to save Armijo iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FirstDifCentral</td><td>OT_BOOLEAN</td><td>worhp_p.FirstDifCentral</td><td>Use central finite difference quotient for first derivatives</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FocusOnFeas</td><td>OT_BOOLEAN</td><td>worhp_p.FocusOnFeas</td><td>Enable Focus-on-Feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>FocusOnFeasFactor</td><td>OT_REAL</td><td>worhp_p.FocusOnFeasFactor</td><td>Factor in Focus-on-Feasibility mode</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>GammaAlpha</td><td>OT_REAL</td><td>worhp_p.GammaAlpha</td><td>Safety factor for alphamin calculation by Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>GroupMethod</td><td>OT_INTEGER</td><td>worhp_p.GroupMethod</td><td>Select method to determine graph colouring groups</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IgnoreFilterCrit</td><td>OT_BOOLEAN</td><td>worhp_p.IgnoreFilterCrit</td><td>Activate accelerating heuristics for Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncBettsTau</td><td>OT_REAL</td><td>worhp_p.IncBettsTau</td><td>Increase factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncBettsTauMore</td><td>OT_REAL</td><td>worhp_p.IncBettsTauMore</td><td>Larger increase factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncreaseIWS</td><td>OT_REAL</td><td>worhp_p.IncreaseIWS</td><td>Increase factor for estimated integer workspace requirement</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>IncreaseRWS</td><td>OT_REAL</td><td>worhp_p.IncreaseRWS</td><td>Increase factor for estimated real workspace requirement</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Infty</td><td>OT_REAL</td><td>worhp_p.Infty</td><td>Upper bound for numbers to be regarded as finite</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>InftyUnbounded</td><td>OT_REAL</td><td>worhp_p.InftyUnbounded</td><td>Tolerance for unboundedness detection heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>InitialLMest</td><td>OT_BOOLEAN</td><td>worhp_p.InitialLMest</td><td>Enable initial Lagrange multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>KeepAcceptableSol</td><td>OT_BOOLEAN</td><td>worhp_p.KeepAcceptableSol</td><td>Save acceptable solutions as fallback</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LMestQPipComTol</td><td>OT_REAL</td><td>worhp_p.LMestQPipComTol</td><td>IP complementarity tolerance in initial multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LMestQPipResTol</td><td>OT_REAL</td><td>worhp_p.LMestQPipResTol</td><td>IP residual tolerance in initial multiplier estimate</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LinMult</td><td>OT_BOOLEAN</td><td>worhp_p.LinMult</td><td>Control Lagrange multiplier update</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LogLevel</td><td>OT_INTEGER</td><td>worhp_p.LogLevel</td><td>Enable XML logfiles and writing interval</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LogResult</td><td>OT_INTEGER</td><td>worhp_p.LogResult</td><td>Enable XML result logging and detail level</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaF</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaF</td><td>Lowpass-filter update factor for objective values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaG</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaG</td><td>Lowpass-filter update factor for constraint values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassAlphaMerit</td><td>OT_REAL</td><td>worhp_p.LowPassAlphaMerit</td><td>Lowpass-filter update factor for merit function values</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>LowPassFilter</td><td>OT_BOOLEAN</td><td>worhp_p.LowPassFilter</td><td>Enable lowpass-filter termination criterion</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Ma57PivotThresh</td><td>OT_REAL</td><td>worhp_p.Ma57PivotThresh</td><td>Pivoting tolerance for MA57 = CNTL(1)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MatrixCC</td><td>OT_BOOLEAN</td><td>worhp_p.MatrixCC</td><td>Not to be included into a parameter file!</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxCalls</td><td>OT_INTEGER</td><td>worhp_p.MaxCalls</td><td>Upper bound to Reverse Communication calls</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxForce</td><td>OT_INTEGER</td><td>worhp_p.MaxForce</td><td>Maximum number of Force recovery strategy steps</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxGPart</td><td>OT_INTEGER</td><td>worhp_p.MaxGPart</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxIter</td><td>OT_INTEGER</td><td>worhp_p.MaxIter</td><td>Upper bound on major iterations</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxLScounter</td><td>OT_INTEGER</td><td>worhp_p.MaxLScounter</td><td>Control activation of Filter acceleration heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MaxNorm</td><td>OT_BOOLEAN</td><td>worhp_p.MaxNorm</td><td>Select max-norm instead of 1-norm in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MeritFunction</td><td>OT_INTEGER</td><td>worhp_p.MeritFunction</td><td>Select merit function and penalty update [0, 3..5]</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MeritGradTol</td><td>OT_REAL</td><td>worhp_p.MeritGradTol</td><td>Threshold of meritfunction gradient for increasing Hessian regularisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MinBettsTau</td><td>OT_REAL</td><td>worhp_p.MinBettsTau</td><td>Lower bound for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>MoreRelax</td><td>OT_BOOLEAN</td><td>worhp_p.MoreRelax</td><td>Introduce one relaxation variable for every constraint</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>NLPmethod</td><td>OT_INTEGER</td><td>worhp_p.NLPmethod</td><td>Select (1) Meritfunction or (3) Filter globalisation</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>NLPprint</td><td>OT_INTEGER</td><td>worhp_p.NLPprint</td><td>NLP print level [-1..4]</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PairMethod</td><td>OT_INTEGER</td><td>worhp_p.PairMethod</td><td>Select method to determine graph colouring pairgroups</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsBar</td><td>OT_REAL</td><td>worhp_p.PenUpdEpsBar</td><td>Penalty update parameter factor for MeritFunction = 3</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKFac</td><td>OT_REAL</td><td>worhp_p.PenUpdEpsKFac</td><td>Penalty update parameter factor for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdEpsKSequence</td><td>OT_INTEGER</td><td>worhp_p.PenUpdEpsKSequence</td><td>Penalty update parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxDeltaK</td><td>OT_REAL</td><td>worhp_p.PenUpdMaxDeltaK</td><td>Max penalty for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdMaxFac</td><td>OT_REAL</td><td>worhp_p.PenUpdMaxFac</td><td>Max factor for increasing penalty for MeritFunction = 4</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PenUpdRBar</td><td>OT_REAL</td><td>worhp_p.PenUpdRBar</td><td>Penalty update parameter for MeritFunction = 3</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PrecisionF</td><td>OT_REAL</td><td>worhp_p.PrecisionF</td><td>(currently unused) Relative precision of objective</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>PrecisionG</td><td>OT_REAL</td><td>worhp_p.PrecisionG</td><td>(currently unused) Relative precision of constraints</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>QPscaleParam</td><td>OT_REAL</td><td>worhp_p.QPscaleParam</td><td>(currently unused)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>QuadraticProblem</td><td>OT_BOOLEAN</td><td>worhp_p.QuadraticProblem</td><td>Not to be included into a parameter file!</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ReduceBettsTau</td><td>OT_REAL</td><td>worhp_p.ReduceBettsTau</td><td>Decrease factor for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RegStrategy</td><td>OT_INTEGER</td><td>worhp_p.RegStrategy</td><td>Select Hessian regularisation strategy in Filter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ReinitFilter</td><td>OT_BOOLEAN</td><td>worhp_p.ReinitFilter</td><td>Enables Filter-reinitialisation accelerating heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxMaxDelta</td><td>OT_REAL</td><td>worhp_p.RelaxMaxDelta</td><td>Upper bound for accepting the constraint relaxation variable</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxMaxPen</td><td>OT_REAL</td><td>worhp_p.RelaxMaxPen</td><td>Upper bound on the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxRho</td><td>OT_REAL</td><td>worhp_p.RelaxRho</td><td>Update factor for the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RelaxStart</td><td>OT_REAL</td><td>worhp_p.RelaxStart</td><td>Initial value of the constraint relaxation penalty</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>RestUntilFeas</td><td>OT_BOOLEAN</td><td>worhp_p.RestUntilFeas</td><td>Do restoration until a feasible solution is found</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleConIter</td><td>OT_BOOLEAN</td><td>worhp_p.ScaleConIter</td><td>Scale constraints in every iteration</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleFacObj</td><td>OT_REAL</td><td>worhp_p.ScaleFacObj</td><td>Value to scale large objective functions to</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaleFacQP</td><td>OT_REAL</td><td>worhp_p.ScaleFacQP</td><td>Upper bound on resulting matrix norm for QP scaling</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledFD</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledFD</td><td>Use a scaled perturbation for finite differences</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledKKT</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledKKT</td><td>Scale KKT conditions</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledObj</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledObj</td><td>Scale the objective function</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ScaledQP</td><td>OT_BOOLEAN</td><td>worhp_p.ScaledQP</td><td>Scale some matrices handed to the QP</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>StartBettsTau</td><td>OT_REAL</td><td>worhp_p.StartBettsTau</td><td>Initial value for Betts' update dampening term</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingDelta</td><td>OT_REAL</td><td>worhp_p.SwitchingDelta</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingSCV</td><td>OT_REAL</td><td>worhp_p.SwitchingSCV</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>SwitchingSF</td><td>OT_REAL</td><td>worhp_p.SwitchingSF</td><td>Filter switching condition parameter</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TakeQPSol</td><td>OT_BOOLEAN</td><td>worhp_p.TakeQPSol</td><td>Evaluate QP search direction regardless of convergence</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>Timeout</td><td>OT_REAL</td><td>worhp_p.Timeout</td><td>Timeout in seconds</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolComp</td><td>OT_REAL</td><td>worhp_p.TolComp</td><td>Complementarity tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolFeas</td><td>OT_REAL</td><td>worhp_p.TolFeas</td><td>Feasibility tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolOpti</td><td>OT_REAL</td><td>worhp_p.TolOpti</td><td>Optimality tolerance</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TolWeakActive</td><td>OT_REAL</td><td>worhp_p.TolWeakActive</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBig</td><td>OT_BOOLEAN</td><td>worhp_p.TooBig</td><td>Enable too-big termination heuristics</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBigCV</td><td>OT_REAL</td><td>worhp_p.TooBigCV</td><td>Upper bound on constraint violation for too-big heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>TooBigKKT</td><td>OT_REAL</td><td>worhp_p.TooBigKKT</td><td>Upper bound on KKT values for too-big heuristic</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserDF</td><td>OT_BOOLEAN</td><td>worhp_p.UserDF</td><td>Objective gradient values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserDG</td><td>OT_BOOLEAN</td><td>worhp_p.UserDG</td><td>Jacobian values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>worhp_p.UserHM</td><td>Hessian values supplied by caller</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>UserHMstructure</td><td>OT_INTEGER</td><td>worhp_p.UserHMstructure</td><td>Enable automatic Hessian structure generation or checking</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>WeakActiveSet</td><td>OT_BOOLEAN</td><td>worhp_p.WeakActiveSet</td><td>(experimental)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eps</td><td>OT_REAL</td><td>worhp_p.eps</td><td>Machine epsilon</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>initialised</td><td>OT_BOOLEAN</td><td>worhp_p.initialised</td><td>Automatically added initialisation flag.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>CasADi::FXInternal<br />CasADi::WorhpInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>worhp_p.qp.ipBarrier</td><td>IP barrier parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>worhp_p.qp.ipComTol</td><td>IP complementarity tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>worhp_p.qp.ipFracBound</td><td>IP fraction-to-the-boundary parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the IP method. (LAPACK::0|MA57: only available if provided by the user:1|SuperLU::2|PARDISO: only available if provided by the user, subject to license availability:3|MUMPS: currently Linux platforms only:5|WSMP: subject to license availability:6|MA86: experimental, only available if provided by the user:7|MA97:experimental, only available if provided by the user:8)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>worhp_p.qp.ipMinAlpha</td><td>IP line search minimum step size.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxDiv</td><td>The relaxation term is divided by this value if successful.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMax</td><td>Maximum relaxation value.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMin</td><td>Mimimum relaxation value.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>worhp_p.qp.ipRelaxMult</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>worhp_p.qp.ipResTol</td><td>IP residuals tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>worhp_p.qp.ipTryRelax</td><td>Enable relaxation strategy when encountering an error.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.lsItMaxIter</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the iterative linear solver. (none:Deactivate; use a direct linear solver.:0|CGNR::1|CGNE::2|CGS::3|BiCGSTAB::4)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select preconditioner for the iterative linear solver. (none:No preconditioner.:0|static:Static preconditioner (KKT-matrix with constant lower-right block).:1|full:Full KKT-matrix.:2)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.lsRefineMaxIter</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>worhp_p.qp.lsScale</td><td>Enables scaling on linear solver level.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>worhp_p.qp.lsTol</td><td>Tolerance for the linear solver.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>worhp_p.qp.lsTrySimple</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported. Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>worhp_p.qp.maxIter</td><td>Imposes an upper limit on the number of minor solver iterations, i.e. for thequadratic subproblem solver. If the limit is reached before convergence,WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>GenericType()</td><td>Select the solution method used by the QP solver. (ip:Interior-Point method.:1|nsn:Nonsmooth-Newton method.:2|automatic: Prefer IP and fall back to NSN on error.:12)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>worhp_p.qp.nsnBeta</td><td>NSN stepsize decrease factor.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>worhp_p.qp.nsnGradStep</td><td>Enable gradient steps in the NSN method.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>worhp_p.qp.nsnKKT</td><td>NSN KKT tolerance.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>GenericType()</td><td>Select the direct linear solver used by the NSN method. (SuperLU::2|MA48: only available if provided by the user:4)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>worhp_p.qp.nsnMinAlpha</td><td>NSN line search minimum step size.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>worhp_p.qp.nsnSigma</td><td>NSN line search slope parameter.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>GenericType()</td><td>Controls the amount of QP solver output. (none:No output.:0|warn:Print warnings and errors.:1|iterations:Print iterations.:2)</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>worhp_p.qp.scaleIntern</td><td>Enable scaling on QP level.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>worhp_p.qp.strict</td><td>Use strict termination criteria in IP method.</td><td>CasADi::WorhpInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::SimulatorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(initial|step)</td><td>CasADi::FXInternal<br />CasADi::SimulatorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(initial|step)</td><td>CasADi::FXInternal<br />CasADi::SimulatorInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LapackLUDenseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::LapackLUDenseInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LapackLUDense
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LapackLUDenseInternal</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::LapackLUDenseInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CVodesInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>CasADi::FXInternal<br />CasADi::CVodesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::CVodesIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>CasADi::FXInternal<br />CasADi::CVodesInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>CasADi::CVodesInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DSDPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DSDPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>dualTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>gapTol</td><td>OT_REAL</td><td>1e-8</td><td>Convergence criterion based on distance between primal and dual objective</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxIter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of iterations</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>primalTol</td><td>OT_REAL</td><td>1e-4</td><td>Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>stepTol</td><td>OT_REAL</td><td>5e-2</td><td>Terminate the solver if the step length in the primal is below this tolerance.</td><td>CasADi::DSDPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DirectSingleShootingInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to CasADi::Parallelizer</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DirectSingleShooting
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to CasADi::Parallelizer</td><td>CasADi::DirectSingleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LiftedSQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eta</td><td>OT_REAL</td><td>0.0001</td><td>Linesearch parameter: See Nocedal 3.4</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"BFGS"</td><td>BFGS|exact</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of linesearch iterations</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp)</td><td>CasADi::FXInternal<br />CasADi::LiftedSQPInternal</td></tr>
<tr><td>mu_safety</td><td>OT_REAL</td><td>1.1</td><td>Safety factor for linesearch mu</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>num_lifted</td><td>OT_INTEGER</td><td>0</td><td>Number of variables to lift</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.5</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>sigma</td><td>OT_REAL</td><td>1.0</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tau</td><td>OT_REAL</td><td>0.2</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>toldx</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion for the stepsize</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>tolgl</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion for the Lagrangian gradient</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::LiftedSQP
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>eta</td><td>OT_REAL</td><td>0.0001</td><td>Linesearch parameter: See Nocedal 3.4</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"BFGS"</td><td>BFGS|exact</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of SQP iterations</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>100</td><td>Maximum number of linesearch iterations</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp)</td><td>CasADi::FXInternal<br />CasADi::LiftedSQPInternal</td></tr>
<tr><td>mu_safety</td><td>OT_REAL</td><td>1.1</td><td>Safety factor for linesearch mu</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>num_lifted</td><td>OT_INTEGER</td><td>0</td><td>Number of variables to lift</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.5</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>sigma</td><td>OT_REAL</td><td>1.0</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tau</td><td>OT_REAL</td><td>0.2</td><td>Linesearch parameter</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>toldx</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion for the stepsize</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>tolgl</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion for the Lagrangian gradient</td><td>CasADi::LiftedSQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::SCPgenInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"limited-memory"</td><td>limited-memory|exact</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>CasADi::FXInternal<br />CasADi::SCPgenInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::SCPgen
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"limited-memory"</td><td>limited-memory|exact</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>CasADi::FXInternal<br />CasADi::SCPgenInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>CasADi::SCPgenInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::OOQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::OOQPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>CasADi::OOQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SymbolicQRInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>CasADi::SymbolicQRInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>CasADi::SymbolicQRInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SymbolicQR
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>CasADi::SymbolicQRInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>CasADi::SymbolicQRInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
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
<tr><td>accept_after_max_steps</td><td>OT_INTEGER</td><td>-1</td><td>Accept a trial point after maximal this number of steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td>no</td><td>Always accept the first trial step. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the complementarity conditions. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td>10000000000.0</td><td>"Acceptance" threshold for the dual infeasibility. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td>15</td><td>Number of "acceptable" iterates before triggering termination. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td>1e+20</td><td>"Acceptance" stopping criterion based on objective function change. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td>1e-06</td><td>"Acceptable" convergence tolerance (relative). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>automatic</td><td>How to calculate the Jacobians.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_globalization</td><td>OT_STRING</td><td>obj-constr-filter</td><td>Globalization strategy for the adaptive mu selection mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkt_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for the KKT error in the adaptive mu globalization strategies. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_fact</td><td>OT_REAL</td><td>0.9999</td><td>Sufficient decrease factor for "kkt-error" globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_iters</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of iterations requiring sufficient progress. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_monotone_init_factor</td><td>OT_REAL</td><td>0.8</td><td>Determines the initial value of the barrier parameter when switching to the monotone mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_restore_previous_iterate</td><td>OT_STRING</td><td>no</td><td>Indicates if the previous iterate should be restored if the monotone mode is entered. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_safeguard_factor</td><td>OT_REAL</td><td>0.0</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td>primal</td><td>Method to determine the step size for constraint multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td>10.0</td><td>Tolerance for switching to full equality multiplier steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_min_frac</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for the minimal step size (before switching to restoration phase). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_red_factor</td><td>OT_REAL</td><td>0.5</td><td>Fractional reduction of the trial step size in the backtracking line search. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td>10.0</td><td>Factor for mu in barrier stop test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial point to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td>constant</td><td>Initialization method for bound multipliers (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td>1.0</td><td>Initial value for the bound multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td>1000.0</td><td>Threshold for resetting bound multipliers after the restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial point to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td>1e-08</td><td>Factor for initial relaxation of the bounds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td>no</td><td>Indicates whether it is desired to check for Nan/Inf in derivative matrices (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_cup</td><td>OT_REAL</td><td>1.5</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_hat</td><td>OT_REAL</td><td>2.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_tilde</td><td>OT_REAL</td><td>5.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the complementarity conditions. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td>1000.0</td><td>Maximum allowed least-square guess of constraint multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for resetting equality and inequality multipliers after restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constraint_violation_norm_type</td><td>OT_STRING</td><td>1-norm</td><td>Norm to be used for the constraint violation in the line search. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_compl_avrg_red_fact</td><td>OT_REAL</td><td>1.0</td><td>Complementarity tolerance factor for accepting corrector step (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td>none</td><td>The type of corrector steps that should be taken (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>delta</td><td>OT_REAL</td><td>1.0</td><td>Multiplier for constraint violation in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>delta_y_max</td><td>OT_REAL</td><td>1e+12</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dependency_detection_with_rhs</td><td>OT_STRING</td><td>no</td><td>Indicates if the right hand sides of the constraints should be considered during dependency detection (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dependency_detector</td><td>OT_STRING</td><td>none</td><td>Indicates which linear solver should be used to detect linearly dependent equality constraints. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td>none</td><td>Enable derivative checker (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_first_index</td><td>OT_INTEGER</td><td>-2</td><td>Index of first quantity to be checked by derivative checker (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td>1e-08</td><td>Size of the finite difference perturbation in derivative test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td>no</td><td>Indicates whether information for all estimated derivatives should be printed. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td>0.0001</td><td>Threshold for indicating wrong derivative. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td>1e+20</td><td>Threshold for maximal value of primal iterates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td>1.0</td><td>Desired threshold for the dual infeasibility. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>epsilon_c</td><td>OT_REAL</td><td>0.01</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_min</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_penalty</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition for the penalty function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td>yes</td><td>Determines if the original objective function should be evaluated at restoration phase trial points. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>False</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>None</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>None</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td>no</td><td>Enable heuristics to quickly detect an infeasible problem. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td>0.001</td><td>Threshold for disabling "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td>100000000.0</td><td>Multiplier threshold for activating "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fast_des_fact</td><td>OT_REAL</td><td>0.1</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fast_step_computation</td><td>OT_STRING</td><td>no</td><td>Indicates if the linear system should be solved quickly. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td>5</td><td>Verbosity level for output file. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_margin_fact</td><td>OT_REAL</td><td>1e-05</td><td>Factor determining width of margin for obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_max_margin</td><td>OT_REAL</td><td>1.0</td><td>Maximum width of margin in obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_reset_trigger</td><td>OT_INTEGER</td><td>5</td><td>Number of iterations that trigger the filter reset. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>findiff_perturbation</td><td>OT_REAL</td><td>1e-07</td><td>Size of the finite difference perturbation for derivative approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td>0.0001</td><td>Size of first x-s perturbation tried. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td>average_compl</td><td>Oracle for the barrier parameter when switching to fixed mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td>make_parameter</td><td>Determines how fixed variables should be handled. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_hat</td><td>OT_REAL</td><td>0.04</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the filter margin for the barrier function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_theta</td><td>OT_REAL</td><td>1e-05</td><td>Relaxation factor in the filter margin for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_tilde</td><td>OT_REAL</td><td>4.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>False</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>None</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>None</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>None</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>None</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>exact</td><td>Indicates what Hessian information is to be used. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation_space</td><td>OT_STRING</td><td>nonlinear-variables</td><td>Indicates in which subspace the Hessian information is to be approximated. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether the problem is a quadratic problem (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td>yes</td><td>Indicates whether final points should be projected into original bounds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>False</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>inf_pr_output</td><td>OT_STRING</td><td>original</td><td>Determines what value is printed in the "inf_pr" output column. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td></td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>False</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all equality constraints are linear (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all inequality constraints are linear (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_approximation</td><td>OT_STRING</td><td>exact</td><td>Specifies technique to compute constraint Jacobian (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>None</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_exponent</td><td>OT_REAL</td><td>0.25</td><td>Exponent for mu in the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td>1e-08</td><td>Size of the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_d</td><td>OT_REAL</td><td>1e-05</td><td>Weight for linear damping term (to handle one-sided bounds). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_sigma</td><td>OT_REAL</td><td>10000000000.0</td><td>Factor limiting the deviation of dual variables from primal estimates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_soc</td><td>OT_REAL</td><td>0.99</td><td>Factor in the sufficient reduction rule for second order correction. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_x_dis</td><td>OT_REAL</td><td>100.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_y_dis</td><td>OT_REAL</td><td>10000.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>least_square_init_duals</td><td>OT_STRING</td><td>no</td><td>Least square initialization of all dual variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>least_square_init_primal</td><td>OT_STRING</td><td>no</td><td>Least square initialization of the primal variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_aug_solver</td><td>OT_STRING</td><td>sherman-morrison</td><td>Strategy for solving the augmented system for low-rank Hessian. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val</td><td>OT_REAL</td><td>1.0</td><td>Value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_max</td><td>OT_REAL</td><td>100000000.0</td><td>Upper bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_min</td><td>OT_REAL</td><td>1e-08</td><td>Lower bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_initialization</td><td>OT_STRING</td><td>scalar1</td><td>Initialization strategy for the limited memory quasi-Newton approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td>6</td><td>Maximum size of the history for the limited quasi-Newton Hessian approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td>2</td><td>Threshold for successive iterations where update is skipped. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_special_for_resto</td><td>OT_STRING</td><td>no</td><td>Determines if the quasi-Newton updates should be special during the restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_update_type</td><td>OT_STRING</td><td>bfgs</td><td>Quasi-Newton update formula for the limited memory approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>line_search_method</td><td>OT_STRING</td><td>filter</td><td>Globalization method used in backtracking line search (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td>yes</td><td>Flag indicating that linear scaling is only done if it seems required. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>ma27</td><td>Linear solver used for step computations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td>mc19</td><td>Method for scaling the linear system. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_ignore_singularity</td><td>OT_STRING</td><td>no</td><td>Enables MA27's ability to solve a linear system even if the matrix is singular. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Real workspace memory for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Integer workspace memory for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td>10.0</td><td>Increment factor for workspace size for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma28_pivtol</td><td>OT_REAL</td><td>0.01</td><td>Pivot tolerance for linear solver MA28. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td>yes</td><td>Controls MA57 automatic scaling (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td>16</td><td>Controls block size used by Level 3 BLAS in MA57BD (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td>16</td><td>Node amalgamation parameter (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td>5</td><td>Controls pivot order in MA57 (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td>1.05</td><td>Safety factor for work space memory allocation for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_small_pivot_flag</td><td>OT_INTEGER</td><td>0</td><td>If set to 1, then when small entries defined by CNTL(2) are detected they are removed and the corresponding pivots placed at the end of the factorization.  This can be particularly efficient if the matrix is highly rank deficient. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_nemin</td><td>OT_INTEGER</td><td>32</td><td>Node Amalgamation parameter (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_print_level</td><td>OT_INTEGER</td><td>0</td><td>Debug printing level for the linear solver MA86 (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_small</td><td>OT_REAL</td><td>1e-20</td><td>Zero Pivot Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_static</td><td>OT_REAL</td><td>0.0</td><td>Static Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_u</td><td>OT_REAL</td><td>1e-08</td><td>Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_umax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>magic_steps</td><td>OT_STRING</td><td>no</td><td>Enables magic steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum number of CPU seconds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_filter_resets</td><td>OT_INTEGER</td><td>5</td><td>Maximal allowed number of filter resets (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td>1e+20</td><td>Maximum value of regularization parameter for handling negative curvature. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>3000</td><td>Maximum number of iterations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>64</td><td>Allow "number_of_adj_dir" to grow until it reaches this number</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>64</td><td>Allow "number_of_fwd_dir" to grow until it reaches this number</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_resto_iter</td><td>OT_INTEGER</td><td>3000000</td><td>Maximum number of successive iterations in restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of second order correction trial steps at each iteration. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soft_resto_iters</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterations performed successively in soft restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td>no</td><td>Indicates if we want to do Mehrotra's algorithm. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_alpha_primal</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td>1e-20</td><td>Smallest perturbation of the Hessian block. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td>1</td><td>Minimum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f)</td><td>CasADi::FXInternal<br />CasADi::IpoptInternal</td></tr>
<tr><td>mu_allow_fast_monotone_decrease</td><td>OT_STRING</td><td>yes</td><td>Allow skipping of barrier problem if barrier test is already met. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td>0.1</td><td>Initial value for the barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td>0.2</td><td>Determines linear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td>1000.0</td><td>Factor for initialization of maximum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td>1e-11</td><td>Minimum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td>quality-function</td><td>Oracle for a new barrier parameter in the adaptive strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td>monotone</td><td>Update strategy for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td>1.5</td><td>Determines superlinear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td>0.0</td><td>Desired value of complementarity. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mult_diverg_feasibility_tol</td><td>OT_REAL</td><td>1e-07</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mult_diverg_y_tol</td><td>OT_REAL</td><td>100000000.0</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_dep_tol</td><td>OT_REAL</td><td>-1.0</td><td>Pivot threshold for detection of linearly dependent constraints in MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td>1000</td><td>Percentage increase in the estimated working space for MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td>7</td><td>Controls permuting and scaling in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td>7</td><td>Controls pivot order in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td>1e-06</td><td>Pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td>77</td><td>Controls scaling in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>unnamed_shared_object</td><td>n/a</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>neg_curv_test_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for heuristic to ignore wrong inertia. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>never_use_fact_cgpen_direction</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the fast Chen-Goldfarb direction (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>never_use_piecewise_penalty_ls</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the piecewise penalty method (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td>-1e+19</td><td>any bound less or equal this value will be considered -inf (i.e. not lower bounded). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_constr_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for constraint function gradient size. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_max_gradient</td><td>OT_REAL</td><td>100.0</td><td>Maximum gradient after NLP scaling. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_method</td><td>OT_STRING</td><td>gradient-based</td><td>Select the technique used for scaling the NLP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_min_value</td><td>OT_REAL</td><td>1e-08</td><td>Minimum value of gradient-based scaling values. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_obj_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for objective function gradient size. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td>1e+19</td><td>any bound greater or this value will be considered +inf (i.e. not upper bounded). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nu_inc</td><td>OT_REAL</td><td>0.0001</td><td>Increment of the penalty parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nu_init</td><td>OT_REAL</td><td>1e-06</td><td>Initial value of the penalty parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>num_linear_variables</td><td>OT_INTEGER</td><td>0</td><td>Number of linear variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>False</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>False</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>obj_max_inc</td><td>OT_REAL</td><td>5.0</td><td>Determines the upper bound on the acceptable increase of barrier objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>obj_scaling_factor</td><td>OT_REAL</td><td>1.0</td><td>Scaling factor for the objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td>File name of options file (to overwrite default). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td>File name of desired output file (leave unset for no file output). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>None</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_coarse_size</td><td>OT_INTEGER</td><td>5000</td><td>Maximum Size of Coarse Grid Matrix (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_factor</td><td>OT_REAL</td><td>0.5</td><td>dropping value for incomplete factor (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_schur</td><td>OT_REAL</td><td>0.1</td><td>dropping value for sparsify schur complement factor (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_inverse_norm_factor</td><td>OT_REAL</td><td>5000000.0</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_levels</td><td>OT_INTEGER</td><td>10</td><td>Maximum Size of Grid Levels (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_row_fill</td><td>OT_INTEGER</td><td>10000000</td><td>max fill for each row (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_relative_tol</td><td>OT_REAL</td><td>1e-06</td><td>Relative Residual Convergence (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iterative</td><td>OT_STRING</td><td>no</td><td>Switch on iterative solver in Pardiso library (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td>complete+2x2</td><td>Matching strategy to be used by Pardiso (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_max_droptol_corrections</td><td>OT_INTEGER</td><td>4</td><td>Maximal number of decreases of drop tolerance during one solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_max_iter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of Krylov-Subspace Iteration (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td>0</td><td>Pardiso message level (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td>0</td><td>Enables out-of-core variant of Pardiso (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_redo_symbolic_fact_only_if_inertia_wrong</td><td>OT_STRING</td><td>no</td><td>Toggle for handling case when elements were perturbed by Pardiso. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_repeated_perturbation_means_singular</td><td>OT_STRING</td><td>no</td><td>Interpretation of perturbed elements. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_des_fact</td><td>OT_REAL</td><td>0.2</td><td>a parameter used in penalty parameter computation (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_init_fac</td><td>OT_REAL</td><td>50.0</td><td>a parameter used to choose initial penalty parameterswhen the regularized Newton method is used. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_init_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximal value for the intial penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_init_min</td><td>OT_REAL</td><td>1.0</td><td>Minimal value for the intial penalty parameter for line search(for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_max</td><td>OT_REAL</td><td>1e+30</td><td>Maximal value for the penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_update_compl_tol</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_update_infeasibility_tol</td><td>OT_REAL</td><td>1e-09</td><td>Threshold for infeasibility in penalty parameter update test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_always_cd</td><td>OT_STRING</td><td>no</td><td>Active permanent perturbation of constraint linearization. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td>0.333333333333</td><td>Decrease factor for x-s perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td>8.0</td><td>Increase factor for x-s perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td>100.0</td><td>Increase factor for x-s perturbation for very first perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_infeasi</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_obj</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td>10.0</td><td>Maximal perturbation of an evaluation point. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_info_string</td><td>OT_STRING</td><td>no</td><td>Enables printing of additional info string at end of iteration output. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>5</td><td>Output verbosity level. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td>no</td><td>Switch to print all algorithmic options. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_latex_mode</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_timing_statistics</td><td>OT_STRING</td><td>no</td><td>Switch to print timing statistics. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td>no</td><td>Print all options set by the user. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_balancing_term</td><td>OT_STRING</td><td>none</td><td>The balancing term included in the quality function for centrality. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_centrality</td><td>OT_STRING</td><td>none</td><td>The penalty term for centrality that is included in quality function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td>8</td><td>Maximum number of search steps during direct search procedure determining the optimal centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for components of the quality function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_section_qf_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for the golden section search procedure determining the optimal centering parameter (in the function value space). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_section_sigma_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for the section search procedure determining the optimal centering parameter (in sigma space). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td>no</td><td>Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility threshold for recomputation of multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>True</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>replace_bounds</td><td>OT_STRING</td><td>no</td><td>Indicates if all variable bounds should be replaced by inequality constraints (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td>0.9</td><td>Required reduction of infeasibility before leaving restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_improvement_factor</td><td>OT_REAL</td><td>0.999999999</td><td>Minimal required reduction of residual test ratio in iterative refinement. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_ratio_max</td><td>OT_REAL</td><td>1e-10</td><td>Iterative refinement tolerance (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_ratio_singular</td><td>OT_REAL</td><td>1e-05</td><td>Threshold for declaring linear system singular after failed iterative refinement. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_failure_feasibility_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for primal infeasibility to declare failure of restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_penalty_parameter</td><td>OT_REAL</td><td>1000.0</td><td>Penalty parameter in the restoration phase objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_proximity_weight</td><td>OT_REAL</td><td>1.0</td><td>Weighting factor for the proximity term in restoration phase objective. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.1</td><td>Value in penalty parameter update formula. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_max</td><td>OT_REAL</td><td>100.0</td><td>Scaling threshold for the NLP error. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_phi</td><td>OT_REAL</td><td>2.3</td><td>Exponent for linear barrier function model in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_theta</td><td>OT_REAL</td><td>1.1</td><td>Exponent for current constraint violation in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sb</td><td>OT_STRING</td><td>no</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sigma_max</td><td>OT_REAL</td><td>100.0</td><td>Maximum value of the centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sigma_min</td><td>OT_REAL</td><td>1e-06</td><td>Minimum value of the centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_corr_if_neg_curv</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step in negative curvature iteration (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_corr_in_monotone_mode</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step during monotone barrier parameter mode (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_finalize_solution_call</td><td>OT_STRING</td><td>no</td><td>Indicates if call to NLP::FinalizeSolution after optimization should be suppressed (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial slack to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial slack to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_move</td><td>OT_REAL</td><td>1.81898940355e-12</td><td>Correction size for very small slacks. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td>0.9999</td><td>Required reduction in primal-dual error in the soft restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>True</td><td>function is sparse</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>None</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td>no</td><td>Tells algorithm to switch to restoration phase in first iteration. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>False</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>suppress_all_output</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tau_min</td><td>OT_REAL</td><td>0.99</td><td>Lower bound on fraction-to-the-boundary parameter tau. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_min</td><td>OT_REAL</td><td>1e-06</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_min_fact</td><td>OT_REAL</td><td>0.0001</td><td>Determines constraint violation threshold in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tiny_step_tol</td><td>OT_REAL</td><td>2.22044604925e-15</td><td>Tolerance for detecting numerically insignificant steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tiny_step_y_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for quitting because of numerically insignificant steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-08</td><td>Desired convergence tolerance (relative). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>None</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>vartheta</td><td>OT_REAL</td><td>0.5</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>False</td><td>verbose evaluation -- for debugging</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as bound_frac for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_entire_iterate</td><td>OT_STRING</td><td>no</td><td>Tells algorithm whether to use the GetWarmStartIterate method in the NLP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td>no</td><td>Warm-start for initial point (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as mult_bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum initial value for the equality multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_same_structure</td><td>OT_STRING</td><td>no</td><td>Indicates whether a problem with a structure identical to the previous one is to be solved. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_frac for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_target_mu</td><td>OT_REAL</td><td>0.0</td><td>Unsupported! (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>False</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td>10</td><td>Number of shortened iterations that trigger the watchdog. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of watchdog iterations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_droptol</td><td>OT_REAL</td><td>0.0</td><td>Drop tolerance for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_fillin_limit</td><td>OT_REAL</td><td>0.0</td><td>Fill-in limit for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_iterative</td><td>OT_STRING</td><td>no</td><td>Switches to iterative solver in WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximal number of iterations in iterative WISMP (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_no_pivoting</td><td>OT_STRING</td><td>no</td><td>Use the static pivoting option of WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td>1</td><td>Number of threads to be used in WSMP (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(16) (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option2</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(20) (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td>0.0001</td><td>Pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td>0</td><td>Determines how the matrix is scaled by WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td>1e-18</td><td>WSMP's singularity threshold. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretent inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_write_matrix_iteration</td><td>OT_INTEGER</td><td>-1</td><td>Iteration in which the matrices are written to files. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
</table>
*/
/** \class CasADi::IpoptSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>accept_after_max_steps</td><td>OT_INTEGER</td><td>-1</td><td>Accept a trial point after maximal this number of steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>accept_every_trial_step</td><td>OT_STRING</td><td>no</td><td>Always accept the first trial step. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_compl_inf_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the complementarity conditions. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_constr_viol_tol</td><td>OT_REAL</td><td>0.01</td><td>"Acceptance" threshold for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_dual_inf_tol</td><td>OT_REAL</td><td>10000000000.0</td><td>"Acceptance" threshold for the dual infeasibility. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_iter</td><td>OT_INTEGER</td><td>15</td><td>Number of "acceptable" iterates before triggering termination. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_obj_change_tol</td><td>OT_REAL</td><td>1e+20</td><td>"Acceptance" stopping criterion based on objective function change. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>acceptable_tol</td><td>OT_REAL</td><td>1e-06</td><td>"Acceptable" convergence tolerance (relative). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>automatic</td><td>How to calculate the Jacobians.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_globalization</td><td>OT_STRING</td><td>obj-constr-filter</td><td>Globalization strategy for the adaptive mu selection mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkt_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for the KKT error in the adaptive mu globalization strategies. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_fact</td><td>OT_REAL</td><td>0.9999</td><td>Sufficient decrease factor for "kkt-error" globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_kkterror_red_iters</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of iterations requiring sufficient progress. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_monotone_init_factor</td><td>OT_REAL</td><td>0.8</td><td>Determines the initial value of the barrier parameter when switching to the monotone mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_restore_previous_iterate</td><td>OT_STRING</td><td>no</td><td>Indicates if the previous iterate should be restored if the monotone mode is entered. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>adaptive_mu_safeguard_factor</td><td>OT_REAL</td><td>0.0</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y</td><td>OT_STRING</td><td>primal</td><td>Method to determine the step size for constraint multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_for_y_tol</td><td>OT_REAL</td><td>10.0</td><td>Tolerance for switching to full equality multiplier steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_min_frac</td><td>OT_REAL</td><td>0.05</td><td>Safety factor for the minimal step size (before switching to restoration phase). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>alpha_red_factor</td><td>OT_REAL</td><td>0.5</td><td>Fractional reduction of the trial step size in the backtracking line search. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>barrier_tol_factor</td><td>OT_REAL</td><td>10.0</td><td>Factor for mu in barrier stop test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial point to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_method</td><td>OT_STRING</td><td>constant</td><td>Initialization method for bound multipliers (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_init_val</td><td>OT_REAL</td><td>1.0</td><td>Initial value for the bound multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_mult_reset_threshold</td><td>OT_REAL</td><td>1000.0</td><td>Threshold for resetting bound multipliers after the restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial point to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>bound_relax_factor</td><td>OT_REAL</td><td>1e-08</td><td>Factor for initial relaxation of the bounds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>check_derivatives_for_naninf</td><td>OT_STRING</td><td>no</td><td>Indicates whether it is desired to check for Nan/Inf in derivative matrices (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_cup</td><td>OT_REAL</td><td>1.5</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_hat</td><td>OT_REAL</td><td>2.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>chi_tilde</td><td>OT_REAL</td><td>5.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>compl_inf_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the complementarity conditions. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>con_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_init_max</td><td>OT_REAL</td><td>1000.0</td><td>Maximum allowed least-square guess of constraint multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_mult_reset_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for resetting equality and inequality multipliers after restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constr_viol_tol</td><td>OT_REAL</td><td>0.0001</td><td>Desired threshold for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>constraint_violation_norm_type</td><td>OT_STRING</td><td>1-norm</td><td>Norm to be used for the constraint violation in the line search. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_compl_avrg_red_fact</td><td>OT_REAL</td><td>1.0</td><td>Complementarity tolerance factor for accepting corrector step (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>corrector_type</td><td>OT_STRING</td><td>none</td><td>The type of corrector steps that should be taken (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>delta</td><td>OT_REAL</td><td>1.0</td><td>Multiplier for constraint violation in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>delta_y_max</td><td>OT_REAL</td><td>1e+12</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dependency_detection_with_rhs</td><td>OT_STRING</td><td>no</td><td>Indicates if the right hand sides of the constraints should be considered during dependency detection (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dependency_detector</td><td>OT_STRING</td><td>none</td><td>Indicates which linear solver should be used to detect linearly dependent equality constraints. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test</td><td>OT_STRING</td><td>none</td><td>Enable derivative checker (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_first_index</td><td>OT_INTEGER</td><td>-2</td><td>Index of first quantity to be checked by derivative checker (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_perturbation</td><td>OT_REAL</td><td>1e-08</td><td>Size of the finite difference perturbation in derivative test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_print_all</td><td>OT_STRING</td><td>no</td><td>Indicates whether information for all estimated derivatives should be printed. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>derivative_test_tol</td><td>OT_REAL</td><td>0.0001</td><td>Threshold for indicating wrong derivative. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>diverging_iterates_tol</td><td>OT_REAL</td><td>1e+20</td><td>Threshold for maximal value of primal iterates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>dual_inf_tol</td><td>OT_REAL</td><td>1.0</td><td>Desired threshold for the dual infeasibility. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>epsilon_c</td><td>OT_REAL</td><td>0.01</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_min</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_penalty</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition for the penalty function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>eta_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the Armijo condition. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>evaluate_orig_obj_at_resto_trial</td><td>OT_STRING</td><td>yes</td><td>Determines if the original objective function should be evaluated at restoration phase trial points. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>False</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>None</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>None</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem</td><td>OT_STRING</td><td>no</td><td>Enable heuristics to quickly detect an infeasible problem. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ctol</td><td>OT_REAL</td><td>0.001</td><td>Threshold for disabling "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>expect_infeasible_problem_ytol</td><td>OT_REAL</td><td>100000000.0</td><td>Multiplier threshold for activating "expect_infeasible_problem" option. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fast_des_fact</td><td>OT_REAL</td><td>0.1</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fast_step_computation</td><td>OT_STRING</td><td>no</td><td>Indicates if the linear system should be solved quickly. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>file_print_level</td><td>OT_INTEGER</td><td>5</td><td>Verbosity level for output file. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_margin_fact</td><td>OT_REAL</td><td>1e-05</td><td>Factor determining width of margin for obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_max_margin</td><td>OT_REAL</td><td>1.0</td><td>Maximum width of margin in obj-constr-filter adaptive globalization strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>filter_reset_trigger</td><td>OT_INTEGER</td><td>5</td><td>Number of iterations that trigger the filter reset. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>findiff_perturbation</td><td>OT_REAL</td><td>1e-07</td><td>Size of the finite difference perturbation for derivative approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>first_hessian_perturbation</td><td>OT_REAL</td><td>0.0001</td><td>Size of first x-s perturbation tried. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_mu_oracle</td><td>OT_STRING</td><td>average_compl</td><td>Oracle for the barrier parameter when switching to fixed mode. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>fixed_variable_treatment</td><td>OT_STRING</td><td>make_parameter</td><td>Determines how fixed variables should be handled. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_hat</td><td>OT_REAL</td><td>0.04</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_phi</td><td>OT_REAL</td><td>1e-08</td><td>Relaxation factor in the filter margin for the barrier function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_theta</td><td>OT_REAL</td><td>1e-05</td><td>Relaxation factor in the filter margin for the constraint violation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gamma_tilde</td><td>OT_REAL</td><td>4.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>False</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>None</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>None</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>None</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>None</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>exact</td><td>Indicates what Hessian information is to be used. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_approximation_space</td><td>OT_STRING</td><td>nonlinear-variables</td><td>Indicates in which subspace the Hessian information is to be approximated. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>hessian_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether the problem is a quadratic problem (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>honor_original_bounds</td><td>OT_STRING</td><td>yes</td><td>Indicates whether final points should be projected into original bounds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>False</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>inf_pr_output</td><td>OT_STRING</td><td>original</td><td>Determines what value is printed in the "inf_pr" output column. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td></td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>False</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_c_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all equality constraints are linear (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jac_d_constant</td><td>OT_STRING</td><td>no</td><td>Indicates whether all inequality constraints are linear (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_approximation</td><td>OT_STRING</td><td>exact</td><td>Specifies technique to compute constraint Jacobian (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>None</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_exponent</td><td>OT_REAL</td><td>0.25</td><td>Exponent for mu in the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>jacobian_regularization_value</td><td>OT_REAL</td><td>1e-08</td><td>Size of the regularization for rank-deficient constraint Jacobians. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_d</td><td>OT_REAL</td><td>1e-05</td><td>Weight for linear damping term (to handle one-sided bounds). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_sigma</td><td>OT_REAL</td><td>10000000000.0</td><td>Factor limiting the deviation of dual variables from primal estimates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_soc</td><td>OT_REAL</td><td>0.99</td><td>Factor in the sufficient reduction rule for second order correction. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_x_dis</td><td>OT_REAL</td><td>100.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>kappa_y_dis</td><td>OT_REAL</td><td>10000.0</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>least_square_init_duals</td><td>OT_STRING</td><td>no</td><td>Least square initialization of all dual variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>least_square_init_primal</td><td>OT_STRING</td><td>no</td><td>Least square initialization of the primal variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_aug_solver</td><td>OT_STRING</td><td>sherman-morrison</td><td>Strategy for solving the augmented system for low-rank Hessian. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val</td><td>OT_REAL</td><td>1.0</td><td>Value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_max</td><td>OT_REAL</td><td>100000000.0</td><td>Upper bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_init_val_min</td><td>OT_REAL</td><td>1e-08</td><td>Lower bound on value for B0 in low-rank update. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_initialization</td><td>OT_STRING</td><td>scalar1</td><td>Initialization strategy for the limited memory quasi-Newton approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_history</td><td>OT_INTEGER</td><td>6</td><td>Maximum size of the history for the limited quasi-Newton Hessian approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_max_skipping</td><td>OT_INTEGER</td><td>2</td><td>Threshold for successive iterations where update is skipped. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_special_for_resto</td><td>OT_STRING</td><td>no</td><td>Determines if the quasi-Newton updates should be special during the restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>limited_memory_update_type</td><td>OT_STRING</td><td>bfgs</td><td>Quasi-Newton update formula for the limited memory approximation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>line_search_method</td><td>OT_STRING</td><td>filter</td><td>Globalization method used in backtracking line search (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_scaling_on_demand</td><td>OT_STRING</td><td>yes</td><td>Flag indicating that linear scaling is only done if it seems required. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>ma27</td><td>Linear solver used for step computations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>linear_system_scaling</td><td>OT_STRING</td><td>mc19</td><td>Method for scaling the linear system. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_ignore_singularity</td><td>OT_STRING</td><td>no</td><td>Enables MA27's ability to solve a linear system even if the matrix is singular. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_la_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Real workspace memory for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_liw_init_factor</td><td>OT_REAL</td><td>5.0</td><td>Integer workspace memory for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_meminc_factor</td><td>OT_REAL</td><td>10.0</td><td>Increment factor for workspace size for MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA27. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma27_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma28_pivtol</td><td>OT_REAL</td><td>0.01</td><td>Pivot tolerance for linear solver MA28. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_automatic_scaling</td><td>OT_STRING</td><td>yes</td><td>Controls MA57 automatic scaling (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_block_size</td><td>OT_INTEGER</td><td>16</td><td>Controls block size used by Level 3 BLAS in MA57BD (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_node_amalgamation</td><td>OT_INTEGER</td><td>16</td><td>Node amalgamation parameter (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivot_order</td><td>OT_INTEGER</td><td>5</td><td>Controls pivot order in MA57 (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtol</td><td>OT_REAL</td><td>1e-08</td><td>Pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pivtolmax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum pivot tolerance for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_pre_alloc</td><td>OT_REAL</td><td>1.05</td><td>Safety factor for work space memory allocation for the linear solver MA57. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma57_small_pivot_flag</td><td>OT_INTEGER</td><td>0</td><td>If set to 1, then when small entries defined by CNTL(2) are detected they are removed and the corresponding pivots placed at the end of the factorization.  This can be particularly efficient if the matrix is highly rank deficient. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_nemin</td><td>OT_INTEGER</td><td>32</td><td>Node Amalgamation parameter (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_print_level</td><td>OT_INTEGER</td><td>0</td><td>Debug printing level for the linear solver MA86 (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_small</td><td>OT_REAL</td><td>1e-20</td><td>Zero Pivot Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_static</td><td>OT_REAL</td><td>0.0</td><td>Static Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_u</td><td>OT_REAL</td><td>1e-08</td><td>Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>ma86_umax</td><td>OT_REAL</td><td>0.0001</td><td>Maximum Pivoting Threshold (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>magic_steps</td><td>OT_STRING</td><td>no</td><td>Enables magic steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_cpu_time</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum number of CPU seconds. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_filter_resets</td><td>OT_INTEGER</td><td>5</td><td>Maximal allowed number of filter resets (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_hessian_perturbation</td><td>OT_REAL</td><td>1e+20</td><td>Maximum value of regularization parameter for handling negative curvature. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>3000</td><td>Maximum number of iterations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>64</td><td>Allow "number_of_adj_dir" to grow until it reaches this number</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>64</td><td>Allow "number_of_fwd_dir" to grow until it reaches this number</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_refinement_steps</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_resto_iter</td><td>OT_INTEGER</td><td>3000000</td><td>Maximum number of successive iterations in restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soc</td><td>OT_INTEGER</td><td>4</td><td>Maximum number of second order correction trial steps at each iteration. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>max_soft_resto_iters</td><td>OT_INTEGER</td><td>10</td><td>Maximum number of iterations performed successively in soft restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mehrotra_algorithm</td><td>OT_STRING</td><td>no</td><td>Indicates if we want to do Mehrotra's algorithm. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_alpha_primal</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_hessian_perturbation</td><td>OT_REAL</td><td>1e-20</td><td>Smallest perturbation of the Hessian block. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>min_refinement_steps</td><td>OT_INTEGER</td><td>1</td><td>Minimum number of iterative refinement steps per linear system solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f)</td><td>CasADi::FXInternal<br />CasADi::IpoptInternal</td></tr>
<tr><td>mu_allow_fast_monotone_decrease</td><td>OT_STRING</td><td>yes</td><td>Allow skipping of barrier problem if barrier test is already met. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_init</td><td>OT_REAL</td><td>0.1</td><td>Initial value for the barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_linear_decrease_factor</td><td>OT_REAL</td><td>0.2</td><td>Determines linear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_max_fact</td><td>OT_REAL</td><td>1000.0</td><td>Factor for initialization of maximum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_min</td><td>OT_REAL</td><td>1e-11</td><td>Minimum value for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_oracle</td><td>OT_STRING</td><td>quality-function</td><td>Oracle for a new barrier parameter in the adaptive strategy. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_strategy</td><td>OT_STRING</td><td>monotone</td><td>Update strategy for barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_superlinear_decrease_power</td><td>OT_REAL</td><td>1.5</td><td>Determines superlinear decrease rate of barrier parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mu_target</td><td>OT_REAL</td><td>0.0</td><td>Desired value of complementarity. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mult_diverg_feasibility_tol</td><td>OT_REAL</td><td>1e-07</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mult_diverg_y_tol</td><td>OT_REAL</td><td>100000000.0</td><td>tolerance for deciding if the multipliers are diverging (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_dep_tol</td><td>OT_REAL</td><td>-1.0</td><td>Pivot threshold for detection of linearly dependent constraints in MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_mem_percent</td><td>OT_INTEGER</td><td>1000</td><td>Percentage increase in the estimated working space for MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_permuting_scaling</td><td>OT_INTEGER</td><td>7</td><td>Controls permuting and scaling in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivot_order</td><td>OT_INTEGER</td><td>7</td><td>Controls pivot order in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtol</td><td>OT_REAL</td><td>1e-06</td><td>Pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver MUMPS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>mumps_scaling</td><td>OT_INTEGER</td><td>77</td><td>Controls scaling in MUMPS (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>unnamed_shared_object</td><td>n/a</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>neg_curv_test_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for heuristic to ignore wrong inertia. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>never_use_fact_cgpen_direction</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the fast Chen-Goldfarb direction (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>never_use_piecewise_penalty_ls</td><td>OT_STRING</td><td>no</td><td>Toggle to switch off the piecewise penalty method (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_lower_bound_inf</td><td>OT_REAL</td><td>-1e+19</td><td>any bound less or equal this value will be considered -inf (i.e. not lower bounded). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_constr_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for constraint function gradient size. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_max_gradient</td><td>OT_REAL</td><td>100.0</td><td>Maximum gradient after NLP scaling. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_method</td><td>OT_STRING</td><td>gradient-based</td><td>Select the technique used for scaling the NLP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_min_value</td><td>OT_REAL</td><td>1e-08</td><td>Minimum value of gradient-based scaling values. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_scaling_obj_target_gradient</td><td>OT_REAL</td><td>0.0</td><td>Target value for objective function gradient size. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nlp_upper_bound_inf</td><td>OT_REAL</td><td>1e+19</td><td>any bound greater or this value will be considered +inf (i.e. not upper bounded). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nu_inc</td><td>OT_REAL</td><td>0.0001</td><td>Increment of the penalty parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>nu_init</td><td>OT_REAL</td><td>1e-06</td><td>Initial value of the penalty parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>num_linear_variables</td><td>OT_INTEGER</td><td>0</td><td>Number of linear variables (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>False</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>False</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>obj_max_inc</td><td>OT_REAL</td><td>5.0</td><td>Determines the upper bound on the acceptable increase of barrier objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>obj_scaling_factor</td><td>OT_REAL</td><td>1.0</td><td>Scaling factor for the objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>option_file_name</td><td>OT_STRING</td><td></td><td>File name of options file (to overwrite default). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>output_file</td><td>OT_STRING</td><td></td><td>File name of desired output file (leave unset for no file output). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>None</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_coarse_size</td><td>OT_INTEGER</td><td>5000</td><td>Maximum Size of Coarse Grid Matrix (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_factor</td><td>OT_REAL</td><td>0.5</td><td>dropping value for incomplete factor (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_dropping_schur</td><td>OT_REAL</td><td>0.1</td><td>dropping value for sparsify schur complement factor (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_inverse_norm_factor</td><td>OT_REAL</td><td>5000000.0</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_levels</td><td>OT_INTEGER</td><td>10</td><td>Maximum Size of Grid Levels (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_max_row_fill</td><td>OT_INTEGER</td><td>10000000</td><td>max fill for each row (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iter_relative_tol</td><td>OT_REAL</td><td>1e-06</td><td>Relative Residual Convergence (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_iterative</td><td>OT_STRING</td><td>no</td><td>Switch on iterative solver in Pardiso library (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_matching_strategy</td><td>OT_STRING</td><td>complete+2x2</td><td>Matching strategy to be used by Pardiso (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_max_droptol_corrections</td><td>OT_INTEGER</td><td>4</td><td>Maximal number of decreases of drop tolerance during one solve. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_max_iter</td><td>OT_INTEGER</td><td>500</td><td>Maximum number of Krylov-Subspace Iteration (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_msglvl</td><td>OT_INTEGER</td><td>0</td><td>Pardiso message level (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_out_of_core_power</td><td>OT_INTEGER</td><td>0</td><td>Enables out-of-core variant of Pardiso (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_redo_symbolic_fact_only_if_inertia_wrong</td><td>OT_STRING</td><td>no</td><td>Toggle for handling case when elements were perturbed by Pardiso. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_repeated_perturbation_means_singular</td><td>OT_STRING</td><td>no</td><td>Interpretation of perturbed elements. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pardiso_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretend inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_des_fact</td><td>OT_REAL</td><td>0.2</td><td>a parameter used in penalty parameter computation (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_init_fac</td><td>OT_REAL</td><td>50.0</td><td>a parameter used to choose initial penalty parameterswhen the regularized Newton method is used. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>pen_theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_init_max</td><td>OT_REAL</td><td>100000.0</td><td>Maximal value for the intial penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_init_min</td><td>OT_REAL</td><td>1.0</td><td>Minimal value for the intial penalty parameter for line search(for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_max</td><td>OT_REAL</td><td>1e+30</td><td>Maximal value for the penalty parameter (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_update_compl_tol</td><td>OT_REAL</td><td>10.0</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>penalty_update_infeasibility_tol</td><td>OT_REAL</td><td>1e-09</td><td>Threshold for infeasibility in penalty parameter update test. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_always_cd</td><td>OT_STRING</td><td>no</td><td>Active permanent perturbation of constraint linearization. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_dec_fact</td><td>OT_REAL</td><td>0.333333333333</td><td>Decrease factor for x-s perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact</td><td>OT_REAL</td><td>8.0</td><td>Increase factor for x-s perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>perturb_inc_fact_first</td><td>OT_REAL</td><td>100.0</td><td>Increase factor for x-s perturbation for very first perturbation. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_infeasi</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>piecewisepenalty_gamma_obj</td><td>OT_REAL</td><td>1e-13</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>point_perturbation_radius</td><td>OT_REAL</td><td>10.0</td><td>Maximal perturbation of an evaluation point. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_info_string</td><td>OT_STRING</td><td>no</td><td>Enables printing of additional info string at end of iteration output. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>5</td><td>Output verbosity level. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_documentation</td><td>OT_STRING</td><td>no</td><td>Switch to print all algorithmic options. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_options_latex_mode</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_timing_statistics</td><td>OT_STRING</td><td>no</td><td>Switch to print timing statistics. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>print_user_options</td><td>OT_STRING</td><td>no</td><td>Print all options set by the user. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_balancing_term</td><td>OT_STRING</td><td>none</td><td>The balancing term included in the quality function for centrality. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_centrality</td><td>OT_STRING</td><td>none</td><td>The penalty term for centrality that is included in quality function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_max_section_steps</td><td>OT_INTEGER</td><td>8</td><td>Maximum number of search steps during direct search procedure determining the optimal centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_norm_type</td><td>OT_STRING</td><td>2-norm-squared</td><td>Norm used for components of the quality function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_section_qf_tol</td><td>OT_REAL</td><td>0.0</td><td>Tolerance for the golden section search procedure determining the optimal centering parameter (in the function value space). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>quality_function_section_sigma_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for the section search procedure determining the optimal centering parameter (in sigma space). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y</td><td>OT_STRING</td><td>no</td><td>Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>recalc_y_feas_tol</td><td>OT_REAL</td><td>1e-06</td><td>Feasibility threshold for recomputation of multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>True</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>replace_bounds</td><td>OT_STRING</td><td>no</td><td>Indicates if all variable bounds should be replaced by inequality constraints (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>required_infeasibility_reduction</td><td>OT_REAL</td><td>0.9</td><td>Required reduction of infeasibility before leaving restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_improvement_factor</td><td>OT_REAL</td><td>0.999999999</td><td>Minimal required reduction of residual test ratio in iterative refinement. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_ratio_max</td><td>OT_REAL</td><td>1e-10</td><td>Iterative refinement tolerance (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>residual_ratio_singular</td><td>OT_REAL</td><td>1e-05</td><td>Threshold for declaring linear system singular after failed iterative refinement. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_failure_feasibility_threshold</td><td>OT_REAL</td><td>0.0</td><td>Threshold for primal infeasibility to declare failure of restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_penalty_parameter</td><td>OT_REAL</td><td>1000.0</td><td>Penalty parameter in the restoration phase objective function. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>resto_proximity_weight</td><td>OT_REAL</td><td>1.0</td><td>Weighting factor for the proximity term in restoration phase objective. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>rho</td><td>OT_REAL</td><td>0.1</td><td>Value in penalty parameter update formula. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_max</td><td>OT_REAL</td><td>100.0</td><td>Scaling threshold for the NLP error. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_phi</td><td>OT_REAL</td><td>2.3</td><td>Exponent for linear barrier function model in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>s_theta</td><td>OT_REAL</td><td>1.1</td><td>Exponent for current constraint violation in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sb</td><td>OT_STRING</td><td>no</td><td> (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sigma_max</td><td>OT_REAL</td><td>100.0</td><td>Maximum value of the centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sigma_min</td><td>OT_REAL</td><td>1e-06</td><td>Minimum value of the centering parameter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_corr_if_neg_curv</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step in negative curvature iteration (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_corr_in_monotone_mode</td><td>OT_STRING</td><td>yes</td><td>Skip the corrector step during monotone barrier parameter mode (unsupported!). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>skip_finalize_solution_call</td><td>OT_STRING</td><td>no</td><td>Indicates if call to NLP::FinalizeSolution after optimization should be suppressed (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_frac</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum relative distance from the initial slack to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_bound_push</td><td>OT_REAL</td><td>0.01</td><td>Desired minimum absolute distance from the initial slack to bound. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>slack_move</td><td>OT_REAL</td><td>1.81898940355e-12</td><td>Correction size for very small slacks. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>soft_resto_pderror_reduction_factor</td><td>OT_REAL</td><td>0.9999</td><td>Required reduction in primal-dual error in the soft restoration phase. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>True</td><td>function is sparse</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>None</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>start_with_resto</td><td>OT_STRING</td><td>no</td><td>Tells algorithm to switch to restoration phase in first iteration. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>False</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>suppress_all_output</td><td>OT_STRING</td><td>no</td><td>Undocumented (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tau_min</td><td>OT_REAL</td><td>0.99</td><td>Lower bound on fraction-to-the-boundary parameter tau. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_max_fact</td><td>OT_REAL</td><td>10000.0</td><td>Determines upper bound for constraint violation in the filter. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_min</td><td>OT_REAL</td><td>1e-06</td><td>LIFENG WRITES THIS. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>theta_min_fact</td><td>OT_REAL</td><td>0.0001</td><td>Determines constraint violation threshold in the switching rule. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tiny_step_tol</td><td>OT_REAL</td><td>2.22044604925e-15</td><td>Tolerance for detecting numerically insignificant steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tiny_step_y_tol</td><td>OT_REAL</td><td>0.01</td><td>Tolerance for quitting because of numerically insignificant steps. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1e-08</td><td>Desired convergence tolerance (relative). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>None</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>var_string_md</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>vartheta</td><td>OT_REAL</td><td>0.5</td><td>a parameter used to check if the fast direction can be used asthe line search direction (for Chen-Goldfarb line search). (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>False</td><td>verbose evaluation -- for debugging</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as bound_frac for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_entire_iterate</td><td>OT_STRING</td><td>no</td><td>Tells algorithm whether to use the GetWarmStartIterate method in the NLP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_init_point</td><td>OT_STRING</td><td>no</td><td>Warm-start for initial point (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as mult_bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_mult_init_max</td><td>OT_REAL</td><td>1000000.0</td><td>Maximum initial value for the equality multipliers. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_same_structure</td><td>OT_STRING</td><td>no</td><td>Indicates whether a problem with a structure identical to the previous one is to be solved. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_frac</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_frac for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_slack_bound_push</td><td>OT_REAL</td><td>0.001</td><td>same as slack_bound_push for the regular initializer. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warm_start_target_mu</td><td>OT_REAL</td><td>0.0</td><td>Unsupported! (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>False</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_shortened_iter_trigger</td><td>OT_INTEGER</td><td>10</td><td>Number of shortened iterations that trigger the watchdog. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>watchdog_trial_iter_max</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of watchdog iterations. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_droptol</td><td>OT_REAL</td><td>0.0</td><td>Drop tolerance for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_inexact_fillin_limit</td><td>OT_REAL</td><td>0.0</td><td>Fill-in limit for inexact factorization preconditioner in WISMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_iterative</td><td>OT_STRING</td><td>no</td><td>Switches to iterative solver in WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximal number of iterations in iterative WISMP (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_no_pivoting</td><td>OT_STRING</td><td>no</td><td>Use the static pivoting option of WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_num_threads</td><td>OT_INTEGER</td><td>1</td><td>Number of threads to be used in WSMP (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(16) (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_ordering_option2</td><td>OT_INTEGER</td><td>1</td><td>Determines how ordering is done in WSMP (IPARM(20) (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtol</td><td>OT_REAL</td><td>0.0001</td><td>Pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_pivtolmax</td><td>OT_REAL</td><td>0.1</td><td>Maximum pivot tolerance for the linear solver WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_scaling</td><td>OT_INTEGER</td><td>0</td><td>Determines how the matrix is scaled by WSMP. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_singularity_threshold</td><td>OT_REAL</td><td>1e-18</td><td>WSMP's singularity threshold. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_skip_inertia_check</td><td>OT_STRING</td><td>no</td><td>Always pretent inertia is correct. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
<tr><td>wsmp_write_matrix_iteration</td><td>OT_INTEGER</td><td>-1</td><td>Iteration in which the matrices are written to files. (see IPOPT documentation)</td><td>CasADi::IpoptInternal</td></tr>
</table>
*/
/** \class CasADi::SDPSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SDPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>calc_dual</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if dual should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>calc_p</td><td>OT_BOOLEAN</td><td>true</td><td>Indicate if the P-part of primal solution should be allocated and calculated. You may want to avoid calculating this variable for problems with n large, as is always dense (n x n).</td><td>CasADi::SDPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::AcadoOCPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
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
<tr><td>absolute_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>auto_init</td><td>OT_BOOLEAN</td><td>false</td><td>initialize differential and angebraic states by a forward integration</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>dynamic_sensitivity</td><td>OT_STRING</td><td></td><td>forward_sensitivities or backward_sensitivities</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator</td><td>OT_STRING</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>integrator_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>kkt_tolerance</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_integrator_steps</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_num_iterations</td><td>OT_INTEGER</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_shooting_nodes</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>periodic_bounds</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>print_level</td><td>OT_STRING</td><td>"low"</td><td>"none", "low", "medium", "high", "debug"</td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>relaxation_parameter</td><td>OT_REAL</td><td></td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>start_time</td><td>OT_REAL</td><td>0.0</td><td></td><td>CasADi::AcadoOCPInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IdasInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>CasADi::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>CasADi::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>CasADi::IdasInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>CasADi::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>CasADi::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td><td>CasADi::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>CasADi::IdasInternal</td></tr>
<tr><td>init_z</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the algebraic states</td><td>CasADi::IdasInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td><td>CasADi::IdasInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup)</td><td>CasADi::FXInternal<br />CasADi::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Supress algebraic variables in the error testing</td><td>CasADi::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::IdasIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::IdasInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>CasADi::IdasInternal</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>CasADi::IdasInternal</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>CasADi::IdasInternal</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>CasADi::IdasInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>extra_fsens_calc_ic</td><td>OT_BOOLEAN</td><td>false</td><td>Call calc ic an extra time, with fsens=0</td><td>CasADi::IdasInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>first_time</td><td>OT_REAL</td><td>GenericType()</td><td>First requested time as a fraction of the time interval</td><td>CasADi::IdasInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstolv</td><td>OT_REALVECTOR</td><td></td><td></td><td>CasADi::IdasInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>CasADi::IdasInternal</td></tr>
<tr><td>init_z</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the algebraic states</td><td>CasADi::IdasInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_step_size</td><td>OT_REAL</td><td>0</td><td>Maximim step size</td><td>CasADi::IdasInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(correctInitialConditions|res|resS|resB|rhsQB|bjacB|jtimesB|psetupB|psolveB|psetup)</td><td>CasADi::FXInternal<br />CasADi::IdasInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Supress algebraic variables in the error testing</td><td>CasADi::IdasInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NewtonImplicitInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>CasADi::FXInternal<br />CasADi::NewtonImplicitInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NewtonImplicitSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>CasADi::NewtonImplicitInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>CasADi::FXInternal<br />CasADi::NewtonImplicitInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>num_algebraic</td><td>OT_INTEGER</td><td>0</td><td>Number of algebraic states</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>num_grid_points</td><td>OT_INTEGER</td><td>2</td><td>Number of uniformly distributed grid points for obtaining the solution, does not influence the integration steps</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>time_dependence</td><td>OT_BOOLEAN</td><td>true</td><td>Explicit depencency of time in the DAE</td><td>CasADi::AcadoIntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>num_algebraic</td><td>OT_INTEGER</td><td>0</td><td>Number of algebraic states</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>num_grid_points</td><td>OT_INTEGER</td><td>2</td><td>Number of uniformly distributed grid points for obtaining the solution, does not influence the integration steps</td><td>CasADi::AcadoIntegratorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>time_dependence</td><td>OT_BOOLEAN</td><td>true</td><td>Explicit depencency of time in the DAE</td><td>CasADi::AcadoIntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>just_in_time</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_sparsity</td><td>OT_BOOLEAN</td><td>false</td><td>Propagate sparsity patterns using just-in-time compilation to a CPU or GPU using OpenCL</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>just_in_time</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>just_in_time_sparsity</td><td>OT_BOOLEAN</td><td>false</td><td>Propagate sparsity patterns using just-in-time compilation to a CPU or GPU using OpenCL</td><td>CasADi::SXFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>CasADi::CplexInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>CasADi::CplexInternal</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to. Default: qp.dat</td><td>CasADi::CplexInternal</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format. Default: false</td><td>CasADi::CplexInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td><td>CasADi::CplexInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td><td>CasADi::CplexInternal</td></tr>
</table>
*/
/** \class CasADi::CplexSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>CasADi::CplexInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>CasADi::CplexInternal</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to. Default: qp.dat</td><td>CasADi::CplexInternal</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format. Default: false</td><td>CasADi::CplexInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>CasADi::CplexInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>simplex_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of simplex iterations.</td><td>CasADi::CplexInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol</td><td>OT_REAL</td><td>1E-6</td><td>Tolerance of solver</td><td>CasADi::CplexInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warm_start</td><td>OT_BOOLEAN</td><td>false</td><td>Use warm start with simplex methods (affects only the simplex methods).</td><td>CasADi::CplexInternal</td></tr>
</table>
*/
/** \class CasADi::LapackQRDenseInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LapackQRDense
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>0.0001</td><td>Feasible mode tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>1</td><td>Absolute feasible tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>0</td><td>Allow infeasible iterations</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>LPSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LPSolver</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e-8</td><td>Maximum objective value</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::KnitroInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>CasADi::FXInternal<br />CasADi::KnitroInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::KnitroSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>BarRule</td><td>OT_INTEGER</td><td>0</td><td>Barrier Rule</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Debug</td><td>OT_INTEGER</td><td>0</td><td>Debug level</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Delta</td><td>OT_REAL</td><td>1.0</td><td>Initial region scaling factor</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasModeTol</td><td>OT_REAL</td><td>0.0001</td><td>Feasible mode tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasTol</td><td>OT_REAL</td><td>1e-6</td><td>Feasible tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>FeasTolAbs</td><td>OT_REAL</td><td>1</td><td>Absolute feasible tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Feasible</td><td>OT_BOOLEAN</td><td>0</td><td>Allow infeasible iterations</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>GradOpt</td><td>OT_INTEGER</td><td>1</td><td>Gradient calculation method</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>HessOpt</td><td>OT_INTEGER</td><td>1</td><td>Hessian calculation method</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>HonorBnds</td><td>OT_BOOLEAN</td><td>0</td><td>Enforce bounds</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>InitPt</td><td>OT_BOOLEAN</td><td>0</td><td>Use initial point strategy</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>LPSolver</td><td>OT_BOOLEAN</td><td>0</td><td>Use LPSolver</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>LmSize</td><td>OT_INTEGER</td><td>10</td><td>Memory pairsize limit</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>MaxCgIt</td><td>OT_INTEGER</td><td>0</td><td>Maximum conjugate gradient iterations</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>MaxIt</td><td>OT_INTEGER</td><td>10000</td><td>Iteration limit</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Mu</td><td>OT_REAL</td><td>0.1</td><td>Initial barrier parameter</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Multistart</td><td>OT_BOOLEAN</td><td>0</td><td>Use multistart</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>NewPoint</td><td>OT_BOOLEAN</td><td>0</td><td>Select new-point feature</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ObjRange</td><td>OT_REAL</td><td>1e-8</td><td>Maximum objective value</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OptTol</td><td>OT_REAL</td><td>1e-6</td><td>Relative optimality tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OptTolAbs</td><td>OT_REAL</td><td>0</td><td>Absolute optimality tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>OutLev</td><td>OT_INTEGER</td><td>2</td><td>Log output level</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Pivot</td><td>OT_REAL</td><td>1e-8</td><td>Initial pivot threshold</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Scale</td><td>OT_BOOLEAN</td><td>1</td><td>Perform scaling</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ShiftInit</td><td>OT_BOOLEAN</td><td>1</td><td>Interior-point shifting initial point</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>Soc</td><td>OT_INTEGER</td><td>1</td><td>Second order correction</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>XTol</td><td>OT_REAL</td><td>1e-15</td><td>Relative solution change tolerance</td><td>CasADi::KnitroInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>CasADi::KnitroInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>CasADi::FXInternal<br />CasADi::KnitroInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::CollocationIntegratorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>hotstart</td><td>OT_BOOLEAN</td><td>true</td><td>Initialize the trajectory at the previous solution</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver</td><td>OT_IMPLICITFUNCTION</td><td>GenericType()</td><td>An implicit function solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quadrature_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>An linear solver to solver the quadrature equations</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>quadrature_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the quadrature solver</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>startup_integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An ODE/DAE integrator that can be used to generate a startup trajectory</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>startup_integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the startup integrator</td><td>CasADi::CollocationIntegratorInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>CasADi::SQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>CasADi::SQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"limited-memory"</td><td>limited-memory|exact</td><td>CasADi::SQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>CasADi::SQPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>CasADi::SQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>CasADi::FXInternal<br />CasADi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>CasADi::SQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>CasADi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::SQPMethod
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>CasADi::SQPInternal</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>CasADi::SQPInternal</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the objective function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>expand_g</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expand the constraint function in terms of scalar operations, i.e. MX->SX</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>gauss_newton</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use Gauss Newton Hessian approximation</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_gradient</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate a function for calculating the gradient of the objective</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_hessian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Hessian of the Lagrangian if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>generate_jacobian</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Generate an exact Jacobian of the constraints if not supplied</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"limited-memory"</td><td>limited-memory|exact</td><td>CasADi::SQPInternal</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FX</td><td>FX()</td><td>A function that will be called at each iteration. Input scheme is the same as NLPSolver's output scheme. Output is scalar.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>CasADi::SQPInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>maxiter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>maxiter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>CasADi::SQPInternal</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>CasADi::SQPInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>CasADi::FXInternal<br />CasADi::SQPInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parametric</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Expect F, G, H, J to have an additional input argument appended at the end, denoting fixed parameters.</td><td>CasADi::NLPSolverInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver</td><td>OT_QPSOLVER</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>CasADi::SQPInternal</td></tr>
<tr><td>qp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>CasADi::SQPInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>CasADi::SQPInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>CasADi::SQPInternal</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>CasADi::SQPInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>CasADi::NLPSolverInternal</td></tr>
</table>
*/
/** \class CasADi::FXInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>(serial|openmp|mpi)</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>(serial|openmp|mpi)</td><td>CasADi::ParallelizerInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>save_corrected_input</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::ParallelizerInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DirectMultipleShootingInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to CasADi::Parallelizer</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DirectMultipleShooting
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>GenericType()</td><td>Passed on to CasADi::Parallelizer</td><td>CasADi::DirectMultipleShootingInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::RKIntegratorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>4</td><td>Order of the interpolating polynomials</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::RKIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>expand_f</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the ODE/DAE residual function in an SX graph</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>expand_q</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the quadrature function in an SX graph</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>4</td><td>Order of the interpolating polynomials</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>CasADi::RKIntegratorInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NLPQPInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSOlver used to solve the QPs.</td><td>CasADi::NLPQPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSOlver</td><td>CasADi::NLPQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::NLPQPSolver
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>The NLPSOlver used to solve the QPs.</td><td>CasADi::NLPQPInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLPSOlver</td><td>CasADi::NLPQPInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>CasADi::ImplicitFunctionInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SundialsInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::SundialsIntegrator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-8</td><td>Absolute tolerence  for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>abstolB</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the adjoint sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>adj_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating adjoint directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>fwd_via_sct</td><td>OT_BOOLEAN</td><td>true</td><td>Generate new functions for calculating forward directional derivatives</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>linear_solver</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solverB</td><td>OT_LINEARSOLVER</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>CasADi::IntegratorInternal</td></tr>
<tr><td>upper_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded Jacobian (estimations)</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>upper_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Upper band-width of banded jacobians for backward integration [default: equal to upper_bandwidth]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditioner</td><td>OT_BOOLEAN</td><td>false</td><td>Precondition an iterative solver</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>use_preconditionerB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Precondition an iterative solver for the backwards problem [default: equal to use_preconditioner]</td><td>CasADi::SundialsInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ControlSimulatorInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>control_endpoint</td><td>OT_BOOLEAN</td><td>false</td><td>Include a control value at the end of the simulation domain. Used for interpolation.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>control_interpolation</td><td>OT_STRING</td><td>"none"</td><td>none|nearest|linear</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>minor_grid</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>The local grid used on each major interval, with time normalized to 1. By default, option 'nf' is used to construct a linearly spaced grid.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nf</td><td>OT_INTEGER</td><td>1</td><td>Number of minor grained integration steps per major interval. nf>0 must hold. This option is not used when 'minor_grid' is provided.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>simulator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the simulator</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::ControlSimulator
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>control_endpoint</td><td>OT_BOOLEAN</td><td>false</td><td>Include a control value at the end of the simulation domain. Used for interpolation.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>control_interpolation</td><td>OT_STRING</td><td>"none"</td><td>none|nearest|linear</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>integrator</td><td>OT_INTEGRATOR</td><td>GenericType()</td><td>An integrator creator function</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>integrator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the integrator</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>minor_grid</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>The local grid used on each major interval, with time normalized to 1. By default, option 'nf' is used to construct a linearly spaced grid.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nf</td><td>OT_INTEGER</td><td>1</td><td>Number of minor grained integration steps per major interval. nf>0 must hold. This option is not used when 'minor_grid' is provided.</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>simulator_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the simulator</td><td>CasADi::ControlSimulatorInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::LinearSolverInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
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
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>trans</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>CasADi::LinearSolverInternal</td></tr>
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
/** \class CasADi::DirectCollocationInternal
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
/** \class CasADi::DirectCollocation
\n
\par
<table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_mode</td><td>OT_STRING</td><td>"automatic"</td><td>How to calculate the Jacobians. (forward: only forward mode|reverse: only adjoint mode|automatic: a heuristic decides which is more appropriate)</td><td>CasADi::FXInternal</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>final_time</td><td>OT_REAL</td><td>1.0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate wether statistics must be gathered</td><td>CasADi::FXInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>jacobian_generator</td><td>OT_JACOBIANGENERATOR</td><td>GenericType()</td><td>Function pointer that returns a Jacobian function given a set of desired Jacobian blocks, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_adj_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_adj_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>max_number_of_fwd_dir</td><td>OT_INTEGER</td><td>optimized_num_dir</td><td>Allow \"number_of_fwd_dir\" to grow until it reaches this number</td><td>CasADi::FXInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>CasADi::FXInternal</td></tr>
<tr><td>name</td><td>OT_STRING</td><td>"unnamed_shared_object"</td><td>name of the object</td><td>CasADi::OptionsFunctionalityNode</td></tr>
<tr><td>nlp_solver</td><td>OT_NLPSOLVER</td><td>GenericType()</td><td>An NLPSolver creator function</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>nlp_solver_options</td><td>OT_DICTIONARY</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>CasADi::DirectCollocationInternal</td></tr>
<tr><td>number_of_adj_dir</td><td>OT_INTEGER</td><td>1</td><td>number of adjoint derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_fwd_dir</td><td>OT_INTEGER</td><td>1</td><td>number of forward derivatives to be calculated simultanously</td><td>CasADi::FXInternal</td></tr>
<tr><td>number_of_grid_points</td><td>OT_INTEGER</td><td>20</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>number_of_parameters</td><td>OT_INTEGER</td><td>0</td><td></td><td>CasADi::OCPSolverInternal</td></tr>
<tr><td>numeric_hessian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Hessians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>numeric_jacobian</td><td>OT_BOOLEAN</td><td>false</td><td>Calculate Jacobians numerically (using directional derivatives) rather than with the built-in method</td><td>CasADi::FXInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparse</td><td>OT_BOOLEAN</td><td>true</td><td>function is sparse</td><td>CasADi::FXInternal</td></tr>
<tr><td>sparsity_generator</td><td>OT_SPARSITYGENERATOR</td><td>GenericType()</td><td>Function that provides sparsity for a given input output block, overrides internal routines</td><td>CasADi::FXInternal</td></tr>
<tr><td>store_jacobians</td><td>OT_BOOLEAN</td><td>false</td><td>keep references to generated Jacobians in order to avoid generating identical Jacobians multiple times</td><td>CasADi::FXInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>CasADi::FXInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>verbose evaluation -- for debugging</td><td>CasADi::FXInternal</td></tr>
</table>
*/
