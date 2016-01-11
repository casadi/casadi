/// \cond INTERNAL
/** \class casadi::CSparseCholeskyInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Linsol_csparsecholesky
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CallbackInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \class casadi::Callback
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ClangCompiler
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>flags</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Compile flags for the JIT compiler. Default: None</td><td>casadi::ClangCompiler</td></tr>
<tr><td>include_path</td><td>OT_STRING</td><td>""</td><td>Include paths for the JIT compiler. The include directory shipped with CasADi will be automatically appended.</td><td>casadi::ClangCompiler</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Compiler_clang
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>flags</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Compile flags for the JIT compiler. Default: None</td></tr>
<tr><td>include_path</td><td>OT_STRING</td><td>""</td><td>Include paths for the JIT compiler. The include directory shipped with CasADi will be automatically appended.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CollocationIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td><td>casadi::CollocationIntegrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td><td>casadi::CollocationIntegrator</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
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
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td></tr>
<tr><td>collocation_scheme</td><td>OT_STRING</td><td>"radau"</td><td>Collocation scheme (radau|legendre)</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td></tr>
<tr><td>interpolation_order</td><td>OT_INTEGER</td><td>3</td><td>Order of the interpolating polynomials</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CommonExternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::CompilerInternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \endcond
/** \addtogroup general_Compiler
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CplexInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>barrier_maxiter</td><td>OT_INTEGER</td><td>2100000000</td><td>Maximum number of barrier iterations.</td><td>casadi::CplexInterface</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>convex</td><td>OT_BOOLEAN</td><td>true</td><td>Indicates if the QP is convex or not (affects only the barrier method).</td><td>casadi::CplexInterface</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>dep_check</td><td>OT_STRING</td><td>"off"</td><td>Detect redundant constraints. (automatic:-1|off:0|begin:1|end:2|both:3)</td><td>casadi::CplexInterface</td></tr>
<tr><td>dump_filename</td><td>OT_STRING</td><td>"qp.dat"</td><td>The filename to dump to.</td><td>casadi::CplexInterface</td></tr>
<tr><td>dump_to_file</td><td>OT_BOOLEAN</td><td>false</td><td>Dumps QP to file in CPLEX format.</td><td>casadi::CplexInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
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
/** \addtogroup plugin_Qpsol_cplex
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Linsol_csparse
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::Integrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable CVodes internal warning messages</td><td>casadi::CvodesInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_all_at_once</td><td>OT_BOOLEAN</td><td>true</td><td>Calculate all right hand sides of the sensitivity equations at once</td><td>casadi::CvodesInterface</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::Integrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_multistep_method</td><td>OT_STRING</td><td>"bdf"</td><td>Integrator scheme (bdf|adams)</td><td>casadi::CvodesInterface</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(res|resB|resQB|reset|psetupB|djacB)</td><td>casadi::FunctionInternal<br />casadi::CvodesInterface</td></tr>
<tr><td>nonlinear_solver_iteration</td><td>OT_STRING</td><td>"newton"</td><td>(newton|functional)</td><td>casadi::CvodesInterface</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::Integrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::Integrator</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::Integrator</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::Integrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::Integrator</td></tr>
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
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td></tr>
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
/** \class casadi::External
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::FixedStepIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::FixedStepIntegrator</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::GenericExternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::GurobiInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>vtype</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Type of variables, Each entry can be one of: 'continuous', 'binary', 'integer', 'semicont', 'semiint'</td><td>casadi::GurobiInterface</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Qpsol_gurobi
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>vtype</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Type of variables, Each entry can be one of: 'continuous', 'binary', 'integer', 'semicont', 'semiint'</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::Integrator</td></tr>
<tr><td>calc_ic</td><td>OT_BOOLEAN</td><td>true</td><td>Use IDACalcIC to get consistent initial conditions.</td><td>casadi::IdasInterface</td></tr>
<tr><td>calc_icB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use IDACalcIC to get consistent initial conditions for backwards system [default: equal to calc_ic].</td><td>casadi::IdasInterface</td></tr>
<tr><td>cj_scaling</td><td>OT_BOOLEAN</td><td>false</td><td>IDAS scaling on cj for the user-defined linear solver module</td><td>casadi::IdasInterface</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable IDAS internal warning messages</td><td>casadi::IdasInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
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
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::Integrator</td></tr>
<tr><td>init_xdot</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Initial values for the state derivatives</td><td>casadi::IdasInterface</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
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
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::Integrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::Integrator</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::Integrator</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>suppress_algebraic</td><td>OT_BOOLEAN</td><td>false</td><td>Suppress algebraic variables in the error testing</td><td>casadi::IdasInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::Integrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::Integrator</td></tr>
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
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::ImplicitFixedStepIntegrator</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::ImplicitToNlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::Rootfinder</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jacobian_function</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function object for calculating the Jacobian (autogenerated by default)</td><td>casadi::Rootfinder</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"csparse"</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::Rootfinder</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::Rootfinder</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nlpsol</td><td>OT_STRING</td><td>GenericType()</td><td>Name of solver.</td><td>casadi::ImplicitToNlp</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Rootfinder_nlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>nlpsol</td><td>OT_STRING</td><td>GenericType()</td><td>Name of solver.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Integrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::Integrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::Integrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::Integrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::Integrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::Integrator</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::Integrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::Integrator</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::IpoptInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>con_integer_md</td><td>OT_DICT</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICT</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>con_string_md</td><td>OT_DICT</td><td>None</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ipopt</td><td>OT_DICT</td><td>None</td><td>Options to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::IpoptInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>False</td><td>n/a</td><td>casadi::IpoptInterface</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>print information about execution time</td><td>casadi::IpoptInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>var_integer_md</td><td>OT_DICT</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICT</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>var_string_md</td><td>OT_DICT</td><td>None</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td><td>casadi::IpoptInterface</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_ipopt
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>con_integer_md</td><td>OT_DICT</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about constraints to be passed to IPOPT</td></tr>
<tr><td>con_numeric_md</td><td>OT_DICT</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about constraints to be passed to IPOPT</td></tr>
<tr><td>con_string_md</td><td>OT_DICT</td><td>None</td><td>String metadata (a dictionary with lists of strings) about constraints to be passed to IPOPT</td></tr>
<tr><td>ipopt</td><td>OT_DICT</td><td>None</td><td>Options to be passed to IPOPT</td></tr>
<tr><td>pass_nonlinear_variables</td><td>OT_BOOLEAN</td><td>False</td><td>n/a</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>print information about execution time</td></tr>
<tr><td>var_integer_md</td><td>OT_DICT</td><td>None</td><td>Integer metadata (a dictionary with lists of integers) about variables to be passed to IPOPT</td></tr>
<tr><td>var_numeric_md</td><td>OT_DICT</td><td>None</td><td>Numeric metadata (a dictionary with lists of reals) about variables to be passed to IPOPT</td></tr>
<tr><td>var_string_md</td><td>OT_DICT</td><td>None</td><td>String metadata (a dictionary with lists of strings) about variables to be passed to IPOPT</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Jit
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>hess</td><td>OT_STRING</td><td>GenericType()</td><td>Function body for Hessian</td><td>casadi::Jit</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac</td><td>OT_STRING</td><td>GenericType()</td><td>Function body for Jacobian</td><td>casadi::Jit</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::KernelSum
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::KinsolInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion tolerance</td><td>casadi::KinsolInterface</td></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::Rootfinder</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>disable_internal_warnings</td><td>OT_BOOLEAN</td><td>false</td><td>Disable KINSOL internal warning messages</td><td>casadi::KinsolInterface</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>f_scale</td><td>OT_REALVECTOR</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>gmres|bcgstab|tfqmr</td><td>casadi::KinsolInterface</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jacobian_function</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function object for calculating the Jacobian (autogenerated by default)</td><td>casadi::Rootfinder</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"csparse"</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::Rootfinder</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::Rootfinder</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>dense|banded|iterative|user_defined</td><td>casadi::KinsolInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td></td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.</td><td>casadi::KinsolInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>0</td><td></td><td>casadi::KinsolInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_djac)</td><td>casadi::FunctionInternal<br />casadi::KinsolInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
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
/** \addtogroup plugin_Rootfinder_kinsol
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>contype</td><td>OT_INTEGERVECTOR</td><td></td><td></td><td>casadi::KnitroInterface</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::KnitroInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_knitro
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
/** \class casadi::LapackLu
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>allow_equilibration_failure</td><td>OT_BOOLEAN</td><td>false</td><td></td><td>casadi::LapackLu</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>equilibration</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::LapackLu</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Linsol_lapacklu
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
/** \class casadi::LapackQr
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Linsol_lapackqr
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Linsol
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::MXFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::MapBase
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>Computational strategy for parallelization (serial|openmp)</td><td>casadi::MapBase</td></tr>
<tr><td>reduced_inputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain inputs</td><td>casadi::MapBase</td></tr>
<tr><td>reduced_outputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain outputs</td><td>casadi::MapBase</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::MapReduce
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>Computational strategy for parallelization (serial|openmp)</td><td>casadi::MapReduce</td></tr>
<tr><td>reduced_inputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain inputs</td><td>casadi::MapReduce</td></tr>
<tr><td>reduced_outputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain outputs</td><td>casadi::MapReduce</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::MapSerial
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>parallelization</td><td>OT_STRING</td><td>"serial"</td><td>Computational strategy for parallelization (serial|openmp)</td><td>casadi::MapSerial</td></tr>
<tr><td>reduced_inputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain inputs</td><td>casadi::MapSerial</td></tr>
<tr><td>reduced_outputs</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Reduction for certain outputs</td><td>casadi::MapSerial</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::Mapaccum
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::Newton
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>abstol</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on max(|F|)</td><td>casadi::Newton</td></tr>
<tr><td>abstolStep</td><td>OT_REAL</td><td>1e-12</td><td>Stopping criterion tolerance on step size</td><td>casadi::Newton</td></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::Rootfinder</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jacobian_function</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function object for calculating the Jacobian (autogenerated by default)</td><td>casadi::Rootfinder</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"csparse"</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::Rootfinder</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::Rootfinder</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>1000</td><td>Maximum number of Newton iterations to perform before returning.</td><td>casadi::Newton</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(step|stepsize|J|F|normF)</td><td>casadi::FunctionInternal<br />casadi::Newton</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print_iteration</td><td>OT_BOOLEAN</td><td>false</td><td>Print information about each iteration</td><td>casadi::Newton</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Rootfinder_newton
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
/** \class casadi::Nlpsol
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::OoqpInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>artol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setArTol to OOQP</td><td>casadi::OoqpInterface</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>mutol</td><td>OT_REAL</td><td>1e-8</td><td>tolerance as provided with setMuTol to OOQP</td><td>casadi::OoqpInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print_level</td><td>OT_INTEGER</td><td>0</td><td>Print level. OOQP listens to print_level 0, 10 and 100</td><td>casadi::OoqpInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Qpsol_ooqp
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
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \endcond
/** \class casadi::OptionsFunctionality
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpToNlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nlpsol</td><td>OT_STRING</td><td>GenericType()</td><td>Name of solver.</td><td>casadi::QpToNlp</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Qpsol_nlp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>nlpsol</td><td>OT_STRING</td><td>GenericType()</td><td>Name of solver.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::QpoasesInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>None</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>boundRelaxation</td><td>OT_REAL</td><td>10000.0</td><td>Initial relaxation of bounds to start homotopy  and initial value for far bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>boundTolerance</td><td>OT_REAL</td><td>2.221e-10</td><td>If upper and lower bounds differ less than this tolerance, they are regarded equal, i.e. as  equality constraint.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>enableCholeskyRefactorisation</td><td>OT_INTEGER</td><td>0</td><td>Specifies the frequency of a full re-factorisation of projected Hessian matrix: 0: turns them off,  1: uses them at each iteration etc.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableDriftCorrection</td><td>OT_INTEGER</td><td>1</td><td>Specifies the frequency of drift corrections: 0: turns them off.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>False</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>True</td><td>Enables the use of  far bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>True</td><td>Enables the use of  flipping bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>False</td><td>Enables condition-hardened  (but more expensive) LI test.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>True</td><td>Enables nonzero curvature  tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>True</td><td>Enables ramping.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>False</td><td>Enables automatic  Hessian regularisation.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsDen</td><td>OT_REAL</td><td>2.221e-13</td><td>Denominator tolerance for ratio tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsFlipping</td><td>OT_REAL</td><td>2.221e-13</td><td>Tolerance of squared Cholesky diagonal factor  which triggers flipping bound.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsIterRef</td><td>OT_REAL</td><td>2.221e-14</td><td>Early termination tolerance for iterative  refinement.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsLITests</td><td>OT_REAL</td><td>2.221e-11</td><td>Tolerance for linear independence tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsNZCTests</td><td>OT_REAL</td><td>6.663e-13</td><td>Tolerance for nonzero curvature tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsNum</td><td>OT_REAL</td><td>-2.221e-13</td><td>Numerator tolerance for ratio tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>epsRegularisation</td><td>OT_REAL</td><td>1.1105e-12</td><td>Scaling factor of identity matrix used for  Hessian regularisation.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>finalRamping</td><td>OT_REAL</td><td>1.0</td><td>Final value for ramping strategy.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>growFarBounds</td><td>OT_REAL</td><td>1000.0</td><td>Factor to grow far bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>initialFarBounds</td><td>OT_REAL</td><td>1000000.0</td><td>Initial size for far bounds.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>initialRamping</td><td>OT_REAL</td><td>0.5</td><td>Start value for ramping strategy.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>lower</td><td>Initial status of bounds at first iteration.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>maxDualJump</td><td>OT_REAL</td><td>100000000.0</td><td>Maximum allowed jump in dual variables in  linear independence tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>maxPrimalJump</td><td>OT_REAL</td><td>100000000.0</td><td>Maximum allowed jump in primal variables in  nonzero curvature tests.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>None</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>numRefinementSteps</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of iterative refinement steps.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>numRegularisationSteps</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of successive regularisation steps.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>medium</td><td>Defines the amount of text output during QP solution, see Section 5.7</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>terminationTolerance</td><td>OT_REAL</td><td>2.221e-09</td><td>Relative termination tolerance to stop homotopy.</td><td>casadi::QpoasesInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Qpsol_qpoases
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>CPUtime</td><td>OT_REAL</td><td>None</td><td>The maximum allowed CPU time in seconds for the whole initialisation (and the actually required one on output). Disabled if unset.</td></tr>
<tr><td>boundRelaxation</td><td>OT_REAL</td><td>10000.0</td><td>Initial relaxation of bounds to start homotopy  and initial value for far bounds.</td></tr>
<tr><td>boundTolerance</td><td>OT_REAL</td><td>2.221e-10</td><td>If upper and lower bounds differ less than this tolerance, they are regarded equal, i.e. as  equality constraint.</td></tr>
<tr><td>enableCholeskyRefactorisation</td><td>OT_INTEGER</td><td>0</td><td>Specifies the frequency of a full re-factorisation of projected Hessian matrix: 0: turns them off,  1: uses them at each iteration etc.</td></tr>
<tr><td>enableDriftCorrection</td><td>OT_INTEGER</td><td>1</td><td>Specifies the frequency of drift corrections: 0: turns them off.</td></tr>
<tr><td>enableEqualities</td><td>OT_BOOLEAN</td><td>False</td><td>Specifies whether equalities should be treated  as always active (True) or not (False)</td></tr>
<tr><td>enableFarBounds</td><td>OT_BOOLEAN</td><td>True</td><td>Enables the use of  far bounds.</td></tr>
<tr><td>enableFlippingBounds</td><td>OT_BOOLEAN</td><td>True</td><td>Enables the use of  flipping bounds.</td></tr>
<tr><td>enableFullLITests</td><td>OT_BOOLEAN</td><td>False</td><td>Enables condition-hardened  (but more expensive) LI test.</td></tr>
<tr><td>enableNZCTests</td><td>OT_BOOLEAN</td><td>True</td><td>Enables nonzero curvature  tests.</td></tr>
<tr><td>enableRamping</td><td>OT_BOOLEAN</td><td>True</td><td>Enables ramping.</td></tr>
<tr><td>enableRegularisation</td><td>OT_BOOLEAN</td><td>False</td><td>Enables automatic  Hessian regularisation.</td></tr>
<tr><td>epsDen</td><td>OT_REAL</td><td>2.221e-13</td><td>Denominator tolerance for ratio tests.</td></tr>
<tr><td>epsFlipping</td><td>OT_REAL</td><td>2.221e-13</td><td>Tolerance of squared Cholesky diagonal factor  which triggers flipping bound.</td></tr>
<tr><td>epsIterRef</td><td>OT_REAL</td><td>2.221e-14</td><td>Early termination tolerance for iterative  refinement.</td></tr>
<tr><td>epsLITests</td><td>OT_REAL</td><td>2.221e-11</td><td>Tolerance for linear independence tests.</td></tr>
<tr><td>epsNZCTests</td><td>OT_REAL</td><td>6.663e-13</td><td>Tolerance for nonzero curvature tests.</td></tr>
<tr><td>epsNum</td><td>OT_REAL</td><td>-2.221e-13</td><td>Numerator tolerance for ratio tests.</td></tr>
<tr><td>epsRegularisation</td><td>OT_REAL</td><td>1.1105e-12</td><td>Scaling factor of identity matrix used for  Hessian regularisation.</td></tr>
<tr><td>finalRamping</td><td>OT_REAL</td><td>1.0</td><td>Final value for ramping strategy.</td></tr>
<tr><td>growFarBounds</td><td>OT_REAL</td><td>1000.0</td><td>Factor to grow far bounds.</td></tr>
<tr><td>initialFarBounds</td><td>OT_REAL</td><td>1000000.0</td><td>Initial size for far bounds.</td></tr>
<tr><td>initialRamping</td><td>OT_REAL</td><td>0.5</td><td>Start value for ramping strategy.</td></tr>
<tr><td>initialStatusBounds</td><td>OT_STRING</td><td>lower</td><td>Initial status of bounds at first iteration.</td></tr>
<tr><td>maxDualJump</td><td>OT_REAL</td><td>100000000.0</td><td>Maximum allowed jump in dual variables in  linear independence tests.</td></tr>
<tr><td>maxPrimalJump</td><td>OT_REAL</td><td>100000000.0</td><td>Maximum allowed jump in primal variables in  nonzero curvature tests.</td></tr>
<tr><td>nWSR</td><td>OT_INTEGER</td><td>None</td><td>The maximum number of working set recalculations to be performed during the initial homotopy. Default is 5(nx + nc)</td></tr>
<tr><td>numRefinementSteps</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of iterative refinement steps.</td></tr>
<tr><td>numRegularisationSteps</td><td>OT_INTEGER</td><td>0</td><td>Maximum number of successive regularisation steps.</td></tr>
<tr><td>printLevel</td><td>OT_STRING</td><td>medium</td><td>Defines the amount of text output during QP solution, see Section 5.7</td></tr>
<tr><td>terminationTolerance</td><td>OT_REAL</td><td>2.221e-09</td><td>Relative termination tolerance to stop homotopy.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Qpsol
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::RkIntegrator
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::FixedStepIntegrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::FixedStepIntegrator</td></tr>
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
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::Rootfinder
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>constraints</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Constrain the unknowns. 0 (default): no constraint on ui, 1: ui &gt;= 0.0, -1: ui &lt;= 0.0, 2: ui &gt; 0.0, -2: ui &lt; 0.0.</td><td>casadi::Rootfinder</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>implicit_input</td><td>OT_INTEGER</td><td>0</td><td>Index of the input that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>implicit_output</td><td>OT_INTEGER</td><td>0</td><td>Index of the output that corresponds to the actual root-finding</td><td>casadi::Rootfinder</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jacobian_function</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function object for calculating the Jacobian (autogenerated by default)</td><td>casadi::Rootfinder</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>"csparse"</td><td>User-defined linear solver class. Needed for sensitivities.</td><td>casadi::Rootfinder</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver.</td><td>casadi::Rootfinder</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SXFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>just_in_time_opencl</td><td>OT_BOOLEAN</td><td>false</td><td>Just-in-time compilation for numeric evaluation using OpenCL (experimental)</td><td>casadi::SXFunction</td></tr>
<tr><td>just_in_time_sparsity</td><td>OT_BOOLEAN</td><td>false</td><td>Propagate sparsity patterns using just-in-time compilation to a CPU or GPU using OpenCL</td><td>casadi::SXFunction</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::Scpgen
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::Scpgen</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::Scpgen</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::Scpgen</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>gauss-newton|exact</td><td>casadi::Scpgen</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::Scpgen</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::Scpgen</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>1</td><td>Maximum number of linesearch iterations</td><td>casadi::Scpgen</td></tr>
<tr><td>merit_memsize</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::Scpgen</td></tr>
<tr><td>merit_start</td><td>OT_REAL</td><td>1e-8</td><td>Lower bound for the merit function parameter</td><td>casadi::Scpgen</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx)</td><td>casadi::FunctionInternal<br />casadi::Scpgen</td></tr>
<tr><td>name_x</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Names of the variables.</td><td>casadi::Scpgen</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::Scpgen</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::Scpgen</td></tr>
<tr><td>print_x</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Which variables to print.</td><td>casadi::Scpgen</td></tr>
<tr><td>qpsol</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::Scpgen</td></tr>
<tr><td>qpsol_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::Scpgen</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td><td>casadi::Scpgen</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td><td>casadi::Scpgen</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td><td>casadi::Scpgen</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_scpgen
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1e-4</td><td>Armijo condition, coefficient of decrease in merit</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td></tr>
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
<tr><td>qpsol</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td></tr>
<tr><td>qpsol_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the QP solver</td></tr>
<tr><td>reg_threshold</td><td>OT_REAL</td><td>1e-8</td><td>Threshold for the regularization.</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td></tr>
<tr><td>tol_pr_step</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for the step size</td></tr>
<tr><td>tol_reg</td><td>OT_REAL</td><td>1e-11</td><td>Stopping criterion for regularization</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::ShellCompiler
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc"</td><td>Compiler command</td><td>casadi::ShellCompiler</td></tr>
<tr><td>compiler_setup</td><td>OT_STRING</td><td>"-fPIC -shared"</td><td>Compiler setup command. Intended to be fixed. The 'flag' option is the prefered way to set custom flags.</td><td>casadi::ShellCompiler</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>flags</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Compile flags for the JIT compiler. Default: None</td><td>casadi::ShellCompiler</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Compiler_shell
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc"</td><td>Compiler command</td></tr>
<tr><td>compiler_setup</td><td>OT_STRING</td><td>"-fPIC -shared"</td><td>Compiler setup command. Intended to be fixed. The 'flag' option is the prefered way to set custom flags.</td></tr>
<tr><td>flags</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Compile flags for the JIT compiler. Default: None</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SimplifiedExternal
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SnoptInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td><td>casadi::SnoptInterface</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print file</td><td>OT_STRING</td><td></td><td></td><td>casadi::SnoptInterface</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td><td>casadi::SnoptInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>specs file</td><td>OT_STRING</td><td></td><td></td><td>casadi::SnoptInterface</td></tr>
<tr><td>start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Basis|Warm)</td><td>casadi::SnoptInterface</td></tr>
<tr><td>summary</td><td>OT_BOOLEAN</td><td>true</td><td></td><td>casadi::SnoptInterface</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_snopt
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>detect_linear</td><td>OT_BOOLEAN</td><td>true</td><td>Make an effort to treat linear constraints and linear variables specially.</td></tr>
<tr><td>print file</td><td>OT_STRING</td><td></td><td></td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>print information about execution time</td></tr>
<tr><td>specs file</td><td>OT_STRING</td><td></td><td></td></tr>
<tr><td>start</td><td>OT_STRING</td><td>"Cold"</td><td>(Cold|Basis|Warm)</td></tr>
<tr><td>summary</td><td>OT_BOOLEAN</td><td>true</td><td></td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::SqicInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(lp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Qpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Qpsol_sqic
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>beta</td><td>OT_REAL</td><td>0.8</td><td>Line-search parameter, restoration factor of stepsize</td><td>casadi::Sqpmethod</td></tr>
<tr><td>c1</td><td>OT_REAL</td><td>1E-4</td><td>Armijo condition, coefficient of decrease in merit</td><td>casadi::Sqpmethod</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hessian_approximation</td><td>OT_STRING</td><td>"exact"</td><td>limited-memory|exact</td><td>casadi::Sqpmethod</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>lbfgs_memory</td><td>OT_INTEGER</td><td>10</td><td>Size of L-BFGS memory.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>max_iter</td><td>OT_INTEGER</td><td>50</td><td>Maximum number of SQP iterations</td><td>casadi::Sqpmethod</td></tr>
<tr><td>max_iter_ls</td><td>OT_INTEGER</td><td>3</td><td>Maximum number of linesearch iterations</td><td>casadi::Sqpmethod</td></tr>
<tr><td>merit_memory</td><td>OT_INTEGER</td><td>4</td><td>Size of memory to store history of merit function values</td><td>casadi::Sqpmethod</td></tr>
<tr><td>min_step_size</td><td>OT_REAL</td><td>1e-10</td><td>The size (inf-norm) of the step size should not become smaller than this.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />(eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx|bfgs)</td><td>casadi::FunctionInternal<br />casadi::Sqpmethod</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print_header</td><td>OT_BOOLEAN</td><td>true</td><td>Print the header with problem statistics</td><td>casadi::Sqpmethod</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>true</td><td>Print information about execution time</td><td>casadi::Sqpmethod</td></tr>
<tr><td>qpsol</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td><td>casadi::Sqpmethod</td></tr>
<tr><td>qpsol_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the QP solver</td><td>casadi::Sqpmethod</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td><td>casadi::Sqpmethod</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td><td>casadi::Sqpmethod</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td><td>casadi::Sqpmethod</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_sqpmethod
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
<tr><td>qpsol</td><td>OT_STRING</td><td>GenericType()</td><td>The QP solver to be used by the SQP method</td></tr>
<tr><td>qpsol_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the QP solver</td></tr>
<tr><td>regularize</td><td>OT_BOOLEAN</td><td>false</td><td>Automatic regularization of Lagrange Hessian.</td></tr>
<tr><td>tol_du</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for dual infeasability</td></tr>
<tr><td>tol_pr</td><td>OT_REAL</td><td>1e-6</td><td>Stopping criterion for primal infeasibility</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>augmented_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed down to the augmented integrator, if one is constructed.</td><td>casadi::Integrator</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>exact_jacobian</td><td>OT_BOOLEAN</td><td>true</td><td>Use exact Jacobian information for the forward integration</td><td>casadi::SundialsInterface</td></tr>
<tr><td>exact_jacobianB</td><td>OT_BOOLEAN</td><td>GenericType()</td><td>Use exact Jacobian information for the backward integration [default: equal to exact_jacobian]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>finite_difference_fsens</td><td>OT_BOOLEAN</td><td>false</td><td>Use finite differences to approximate the forward sensitivity equations (if AD is not available)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_abstol</td><td>OT_REAL</td><td>GenericType()</td><td>Absolute tolerence for the forward sensitivity solution [default: equal to abstol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_err_con</td><td>OT_BOOLEAN</td><td>true</td><td>include the forward sensitivities in all error controls</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_reltol</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the forward sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_scaling_factors</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Scaling factor for the components if finite differences is used</td><td>casadi::SundialsInterface</td></tr>
<tr><td>fsens_sensitiviy_parameters</td><td>OT_INTEGERVECTOR</td><td>GenericType()</td><td>Specifies which components will be used when estimating the sensitivity equations</td><td>casadi::SundialsInterface</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grid</td><td>OT_REALVECTOR</td><td>GenericType()</td><td>Time grid</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver</td><td>OT_STRING</td><td>GenericType()</td><td>An implicit function solver</td><td>casadi::Integrator</td></tr>
<tr><td>implicit_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the NLP Solver</td><td>casadi::Integrator</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>interpolation_type</td><td>OT_STRING</td><td>"hermite"</td><td>Type of interpolation for the adjoint sensitivities (hermite|polynomial)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solver</td><td>OT_STRING</td><td>"gmres"</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>iterative_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>(gmres|bcgstab|tfqmr)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>linear_solver</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solverB</td><td>OT_STRING</td><td>GenericType()</td><td>A custom linear solver creator function for backwards integration [default: equal to linear_solver]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_optionsB</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the linear solver for backwards integration [default: equal to linear_solver_options]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_type</td><td>OT_STRING</td><td>"dense"</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>linear_solver_typeB</td><td>OT_STRING</td><td>GenericType()</td><td>(user_defined|dense|banded|iterative)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidth</td><td>OT_INTEGER</td><td>GenericType()</td><td>Lower band-width of banded Jacobian (estimations)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>lower_bandwidthB</td><td>OT_INTEGER</td><td>GenericType()</td><td>lower band-width of banded jacobians for backward integration [default: equal to lower_bandwidth]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylov</td><td>OT_INTEGER</td><td>10</td><td>Maximum Krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_krylovB</td><td>OT_INTEGER</td><td>GenericType()</td><td>Maximum krylov subspace size</td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_multistep_order</td><td>OT_INTEGER</td><td>5</td><td></td><td>casadi::SundialsInterface</td></tr>
<tr><td>max_num_steps</td><td>OT_INTEGER</td><td>10000</td><td>Maximum number of integrator steps</td><td>casadi::SundialsInterface</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>number_of_finite_elements</td><td>OT_INTEGER</td><td>20</td><td>Number of finite elements</td><td>casadi::Integrator</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_t0</td><td>OT_BOOLEAN</td><td>false</td><td>Output the state at the initial time</td><td>casadi::Integrator</td></tr>
<tr><td>pretype</td><td>OT_STRING</td><td>"none"</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>pretypeB</td><td>OT_STRING</td><td>GenericType()</td><td>(none|left|right|both)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>print_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Print out statistics after integration</td><td>casadi::Integrator</td></tr>
<tr><td>quad_err_con</td><td>OT_BOOLEAN</td><td>false</td><td>Should the quadratures affect the step size control</td><td>casadi::SundialsInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>reltol</td><td>OT_REAL</td><td>1e-6</td><td>Relative tolerence for the IVP solution</td><td>casadi::SundialsInterface</td></tr>
<tr><td>reltolB</td><td>OT_REAL</td><td>GenericType()</td><td>Relative tolerence for the adjoint sensitivity solution [default: equal to reltol]</td><td>casadi::SundialsInterface</td></tr>
<tr><td>sensitivity_method</td><td>OT_STRING</td><td>"simultaneous"</td><td>(simultaneous|staggered)</td><td>casadi::SundialsInterface</td></tr>
<tr><td>steps_per_checkpoint</td><td>OT_INTEGER</td><td>20</td><td>Number of steps between two consecutive checkpoints</td><td>casadi::SundialsInterface</td></tr>
<tr><td>stop_at_end</td><td>OT_BOOLEAN</td><td>true</td><td>Stop the integrator at the end of the interval</td><td>casadi::SundialsInterface</td></tr>
<tr><td>t0</td><td>OT_REAL</td><td>0.0</td><td>Beginning of the time horizon</td><td>casadi::Integrator</td></tr>
<tr><td>tf</td><td>OT_REAL</td><td>1.0</td><td>End of the time horizon</td><td>casadi::Integrator</td></tr>
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
/** \class casadi::Switch
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
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
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>codegen</td><td>OT_BOOLEAN</td><td>false</td><td>C-code generation</td><td>casadi::SymbolicQr</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"gcc -fPIC -O2"</td><td>Compiler command to be used for compiling generated code</td><td>casadi::SymbolicQr</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Linsol_symbolicqr
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
/** \class casadi::WorhpInterface
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>Ares</td><td>OT_INTEGERVECTOR</td><td>None</td><td>Armijo recovery strategies</td><td>casadi::WorhpInterface</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>True</td><td>Use exact Hessian</td><td>casadi::WorhpInterface</td></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)<br />(qp)</td><td>casadi::OptionsFunctionalityNode<br />casadi::Nlpsol</td></tr>
<tr><td>eval_errors_fatal</td><td>OT_BOOLEAN</td><td>false</td><td>When errors occur during evaluation of f,g,...,stop the iterations</td><td>casadi::Nlpsol</td></tr>
<tr><td>expand</td><td>OT_BOOLEAN</td><td>false</td><td>Expand the NLP function in terms of scalar operations, i.e. MX-&gt;SX</td><td>casadi::Nlpsol</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>grad_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the objective (column, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the gradient of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>grad_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated gradient of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Hessian of the Lagrangian (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>hess_lag_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Hessian of the Lagrangian.</td><td>casadi::Nlpsol</td></tr>
<tr><td>ignore_check_vec</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, the input shape of F will not be checked.</td><td>casadi::Nlpsol</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>iteration_callback</td><td>OT_FUNCTION</td><td>GenericType()</td><td>A function that will be called at each iteration with the solver as input. Check documentation of Callback.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_ignore_errors</td><td>OT_BOOLEAN</td><td>false</td><td>If set to true, errors thrown by iteration_callback will be ignored.</td><td>casadi::Nlpsol</td></tr>
<tr><td>iteration_callback_step</td><td>OT_INTEGER</td><td>1</td><td>Only call the callback function every few iterations.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the objective (sparse row, autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_f_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the objective.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g</td><td>OT_FUNCTION</td><td>GenericType()</td><td>Function for calculating the Jacobian of the constraints (autogenerated by default)</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_g_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options for the autogenerated Jacobian of the constraints.</td><td>casadi::Nlpsol</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)<br />Monitor functions (eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h)</td><td>casadi::FunctionInternal<br />casadi::WorhpInterface</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>Print information about execution time</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>None</td><td>IP barrier parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>None</td><td>IP complementarity tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>None</td><td>IP fraction-to-the-boundary parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the IP method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>None</td><td>IP line search minimum step size.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>None</td><td>The relaxation term is divided by this value if successful.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>None</td><td>Maximum relaxation value.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>None</td><td>Mimimum relaxation value.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>None</td><td>The relaxation term is multiplied by this value if unsuccessful.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>None</td><td>IP residuals tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>None</td><td>Enable relaxation strategy when encountering an error.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>None</td><td>Maximum number of iterations of the iterative linear solvers.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>None</td><td>Select the iterative linear solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>None</td><td>Select preconditioner for the iterative linear solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>None</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>None</td><td>Enables scaling on linear solver level.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>None</td><td>Tolerance for the linear solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>None</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>None</td><td>Imposes an upper limit on the number of minor solver iterations,  i.e. for the quadratic subproblem solver.If the limit is reached before convergence, WORHP will activate QP recovery strategies to prevent a solver breakdown.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>None</td><td>Select the solution method used by the QP solver.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>None</td><td>NSN stepsize decrease factor.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>None</td><td>Enable gradient steps in the NSN method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>None</td><td>NSN KKT tolerance.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the NSN method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>None</td><td>NSN line search minimum step size.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>None</td><td>NSN line search slope parameter.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>None</td><td>Controls the amount of QP solver output.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>None</td><td>Enable scaling on QP level.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>None</td><td>Use strict termination criteria in IP method.</td><td>casadi::WorhpInterface</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose_init</td><td>OT_BOOLEAN</td><td>false</td><td>Print out timing information about the different stages of initialization</td><td>casadi::Nlpsol</td></tr>
<tr><td>warn_initial_bounds</td><td>OT_BOOLEAN</td><td>false</td><td>Warn if the initial guess does not satisfy LBX and UBX</td><td>casadi::Nlpsol</td></tr>
</table>
*/
/// \endcond
/** \addtogroup plugin_Nlpsol_worhp
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td>Ares</td><td>OT_INTEGERVECTOR</td><td>None</td><td>Armijo recovery strategies</td></tr>
<tr><td>UserHM</td><td>OT_BOOLEAN</td><td>True</td><td>Use exact Hessian</td></tr>
<tr><td>print_time</td><td>OT_BOOLEAN</td><td>True</td><td>Print information about execution time</td></tr>
<tr><td>qp_ipBarrier</td><td>OT_REAL</td><td>None</td><td>IP barrier parameter.</td></tr>
<tr><td>qp_ipComTol</td><td>OT_REAL</td><td>None</td><td>IP complementarity tolerance.</td></tr>
<tr><td>qp_ipFracBound</td><td>OT_REAL</td><td>None</td><td>IP fraction-to-the-boundary parameter.</td></tr>
<tr><td>qp_ipLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the IP method.</td></tr>
<tr><td>qp_ipMinAlpha</td><td>OT_REAL</td><td>None</td><td>IP line search minimum step size.</td></tr>
<tr><td>qp_ipRelaxDiv</td><td>OT_REAL</td><td>None</td><td>The relaxation term is divided by this value if successful.</td></tr>
<tr><td>qp_ipRelaxMax</td><td>OT_REAL</td><td>None</td><td>Maximum relaxation value.</td></tr>
<tr><td>qp_ipRelaxMin</td><td>OT_REAL</td><td>None</td><td>Mimimum relaxation value.</td></tr>
<tr><td>qp_ipRelaxMult</td><td>OT_REAL</td><td>None</td><td>The relaxation term is multiplied by this value if unsuccessful.</td></tr>
<tr><td>qp_ipResTol</td><td>OT_REAL</td><td>None</td><td>IP residuals tolerance.</td></tr>
<tr><td>qp_ipTryRelax</td><td>OT_BOOLEAN</td><td>None</td><td>Enable relaxation strategy when encountering an error.</td></tr>
<tr><td>qp_lsItMaxIter</td><td>OT_INTEGER</td><td>None</td><td>Maximum number of iterations of the iterative linear solvers.</td></tr>
<tr><td>qp_lsItMethod</td><td>OT_STRING</td><td>None</td><td>Select the iterative linear solver.</td></tr>
<tr><td>qp_lsItPrecondMethod</td><td>OT_STRING</td><td>None</td><td>Select preconditioner for the iterative linear solver.</td></tr>
<tr><td>qp_lsRefineMaxIter</td><td>OT_INTEGER</td><td>None</td><td>Maximum number of iterative refinement steps of the direct linear solvers.</td></tr>
<tr><td>qp_lsScale</td><td>OT_BOOLEAN</td><td>None</td><td>Enables scaling on linear solver level.</td></tr>
<tr><td>qp_lsTol</td><td>OT_REAL</td><td>None</td><td>Tolerance for the linear solver.</td></tr>
<tr><td>qp_lsTrySimple</td><td>OT_BOOLEAN</td><td>None</td><td>Some matrices can be solved without calling a linear equation solver.Currently only diagonal matrices are supported.Non-diagonal matrices will besolved with the chosen linear equation solver.</td></tr>
<tr><td>qp_maxIter</td><td>OT_INTEGER</td><td>None</td><td>Imposes an upper limit on the number of minor solver iterations,  i.e. for the quadratic subproblem solver.If the limit is reached before convergence, WORHP will activate QP recovery strategies to prevent a solver breakdown.</td></tr>
<tr><td>qp_method</td><td>OT_STRING</td><td>None</td><td>Select the solution method used by the QP solver.</td></tr>
<tr><td>qp_nsnBeta</td><td>OT_REAL</td><td>None</td><td>NSN stepsize decrease factor.</td></tr>
<tr><td>qp_nsnGradStep</td><td>OT_BOOLEAN</td><td>None</td><td>Enable gradient steps in the NSN method.</td></tr>
<tr><td>qp_nsnKKT</td><td>OT_REAL</td><td>None</td><td>NSN KKT tolerance.</td></tr>
<tr><td>qp_nsnLsMethod</td><td>OT_STRING</td><td>None</td><td>Select the direct linear solver used by the NSN method.</td></tr>
<tr><td>qp_nsnMinAlpha</td><td>OT_REAL</td><td>None</td><td>NSN line search minimum step size.</td></tr>
<tr><td>qp_nsnSigma</td><td>OT_REAL</td><td>None</td><td>NSN line search slope parameter.</td></tr>
<tr><td>qp_printLevel</td><td>OT_STRING</td><td>None</td><td>Controls the amount of QP solver output.</td></tr>
<tr><td>qp_scaleIntern</td><td>OT_BOOLEAN</td><td>None</td><td>Enable scaling on QP level.</td></tr>
<tr><td>qp_strict</td><td>OT_BOOLEAN</td><td>None</td><td>Use strict termination criteria in IP method.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::XFunction
\n
\par
<a name='options'></a><table>
<caption>List of available options</caption>
<tr><th>Id</th><th>Type</th><th>Default</th><th>Description</th><th>Used in</th></tr>
<tr><td>ad_weight</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for derivative calculation.When there is an option of either using forward or reverse mode directional derivatives, the condition ad_weight*nf&lt;=(1-ad_weight)*na is used where nf and na are estimates of the number of forward/reverse mode directional derivatives needed. By default, ad_weight is calculated automatically, but this can be overridden by setting this option. In particular, 0 means forcing forward mode and 1 forcing reverse mode. Leave unset for (class specific) heuristics.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>ad_weight_sp</td><td>OT_REAL</td><td>GenericType()</td><td>Weighting factor for sparsity pattern calculation calculation.Overrides default behavior. Set to 0 and 1 to force forward and reverse mode respectively. Cf. option \"ad_weight\".</td><td>casadi::FunctionInternal</td></tr>
<tr><td>compiler</td><td>OT_STRING</td><td>"clang"</td><td>Just-in-time compiler plugin to be used.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>defaults_recipes</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Changes default options according to a given recipe (low-level)</td><td>casadi::OptionsFunctionalityNode</td></tr>
<tr><td>gather_stats</td><td>OT_BOOLEAN</td><td>false</td><td>Flag to indicate whether statistics must be gathered</td><td>casadi::FunctionInternal</td></tr>
<tr><td>input_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom input scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>inputs_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when the numerical values of the inputs don't make sense</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jac_penalty</td><td>OT_REAL</td><td>2</td><td>When requested for a number of forward/reverse directions,   it may be cheaper to compute first the full jacobian and then multiply with seeds, rather than obtain the requested directions in a straightforward manner. Casadi uses a heuristic to decide which is cheaper. A high value of 'jac_penalty' makes it less likely for the heurstic to chose the full Jacobian strategy. The special value -1 indicates never to use the full Jacobian strategy</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit</td><td>OT_BOOLEAN</td><td>false</td><td>Use just-in-time compiler to speed up the evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>jit_options</td><td>OT_DICT</td><td>GenericType()</td><td>Options to be passed to the jit compiler.</td><td>casadi::FunctionInternal</td></tr>
<tr><td>monitor</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Monitors to be activated (inputs|outputs)</td><td>casadi::FunctionInternal</td></tr>
<tr><td>output_scheme</td><td>OT_STRINGVECTOR</td><td>GenericType()</td><td>Custom output scheme</td><td>casadi::FunctionInternal</td></tr>
<tr><td>regularity_check</td><td>OT_BOOLEAN</td><td>true</td><td>Throw exceptions when NaN or Inf appears during evaluation</td><td>casadi::FunctionInternal</td></tr>
<tr><td>user_data</td><td>OT_VOIDPTR</td><td>GenericType()</td><td>A user-defined field that can be used to identify the function or pass additional information</td><td>casadi::FunctionInternal</td></tr>
<tr><td>verbose</td><td>OT_BOOLEAN</td><td>false</td><td>Verbose evaluation -- for debugging</td><td>casadi::FunctionInternal</td></tr>
</table>
*/
/// \endcond
