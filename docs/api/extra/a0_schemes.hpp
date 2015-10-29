/** \defgroup scheme_IntegratorOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 6) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td></td><td>Differential state at the final time.</td></tr>
<tr><td>INTEGRATOR_QF</td><td></td><td>Quadrature state at the final time.</td></tr>
<tr><td>INTEGRATOR_ZF</td><td></td><td>Algebraic variable at the final time.</td></tr>
<tr><td>INTEGRATOR_RXF</td><td></td><td>Backward differential state at the initial time.</td></tr>
<tr><td>INTEGRATOR_RQF</td><td></td><td>Backward quadrature state at the initial time.</td></tr>
<tr><td>INTEGRATOR_RZF</td><td></td><td>Backward algebraic variable at the initial time.</td></tr>
</table>
*/
/** \defgroup scheme_HessLagOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::HessLagOutput  (HESSLAG_NUM_OUT = 5) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HESSLAG_HESS</td><td></td><td>Hessian of the Lagrangian.</td></tr>
<tr><td>HESSLAG_F</td><td></td><td>Objective function.</td></tr>
<tr><td>HESSLAG_G</td><td></td><td>Constraint function.</td></tr>
<tr><td>HESSLAG_GRAD_X</td><td></td><td>Gradient of the Lagrangian with respect to x.</td></tr>
<tr><td>HESSLAG_GRAD_P</td><td></td><td>Gradient of the Lagrangian with respect to p.</td></tr>
</table>
*/
/** \defgroup scheme_LinsolInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::LinsolInput  (LINSOL_NUM_IN = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_A</td><td></td><td>The square matrix A: sparse, (n x n)</td></tr>
<tr><td>LINSOL_B</td><td></td><td>The right-hand-side matrix b: dense, (n x m)</td></tr>
</table>
*/
/** \defgroup scheme_QpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::QpSolverInput  (QP_SOLVER_NUM_IN = 9) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_H</td><td></td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_SOLVER_G</td><td></td><td>The vector g: dense, (n x 1)</td></tr>
<tr><td>QP_SOLVER_A</td><td></td><td>The matrix A: sparse, (nc x n) - product with x must be dense.</td></tr>
<tr><td>QP_SOLVER_LBA</td><td></td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_SOLVER_UBA</td><td></td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_SOLVER_LBX</td><td></td><td>dense, (n x 1)</td></tr>
<tr><td>QP_SOLVER_UBX</td><td></td><td>dense, (n x 1)</td></tr>
<tr><td>QP_SOLVER_X0</td><td></td><td>dense, (n x 1)</td></tr>
<tr><td>QP_SOLVER_LAM_X0</td><td></td><td>dense</td></tr>
</table>
*/
/** \defgroup scheme_RDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::RDAEInput  (RDAE_NUM_IN = 7) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>RDAE_RX</td><td></td><td>Backward differential state.</td></tr>
<tr><td>RDAE_RZ</td><td></td><td>Backward algebraic state.</td></tr>
<tr><td>RDAE_RP</td><td></td><td>Backward parameter vector.</td></tr>
<tr><td>RDAE_X</td><td></td><td>Forward differential state.</td></tr>
<tr><td>RDAE_Z</td><td></td><td>Forward algebraic state.</td></tr>
<tr><td>RDAE_P</td><td></td><td>Parameter vector.</td></tr>
<tr><td>RDAE_T</td><td></td><td>Explicit time dependence.</td></tr>
</table>
*/
/** \defgroup scheme_NLPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NLPOutput  (NL_NUM_OUT = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_F</td><td></td><td>Objective function.</td></tr>
<tr><td>NL_G</td><td></td><td>Constraint function.</td></tr>
</table>
*/
/** \defgroup scheme_QpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::QpSolverOutput  (QP_SOLVER_NUM_OUT = 4) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_X</td><td></td><td>The primal solution.</td></tr>
<tr><td>QP_SOLVER_COST</td><td></td><td>The optimal cost.</td></tr>
<tr><td>QP_SOLVER_LAM_A</td><td></td><td>The dual solution corresponding to linear bounds.</td></tr>
<tr><td>QP_SOLVER_LAM_X</td><td></td><td>The dual solution corresponding to simple bounds.</td></tr>
</table>
*/
/** \defgroup scheme_DAEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::DAEInput  (DAE_NUM_IN = 4) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_X</td><td></td><td>Differential state.</td></tr>
<tr><td>DAE_Z</td><td></td><td>Algebraic state.</td></tr>
<tr><td>DAE_P</td><td></td><td>Parameter.</td></tr>
<tr><td>DAE_T</td><td></td><td>Explicit time dependence.</td></tr>
</table>
*/
/** \defgroup scheme_DAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::DAEOutput  (DAE_NUM_OUT = 3) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_ODE</td><td></td><td>Right hand side of the implicit ODE.</td></tr>
<tr><td>DAE_ALG</td><td></td><td>Right hand side of algebraic equations.</td></tr>
<tr><td>DAE_QUAD</td><td></td><td>Right hand side of quadratures equations.</td></tr>
</table>
*/
/** \defgroup scheme_GradFInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::GradFInput  (GRADF_NUM_IN = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_X</td><td></td><td>Decision variable.</td></tr>
<tr><td>GRADF_P</td><td></td><td>Fixed parameter.</td></tr>
</table>
*/
/** \defgroup scheme_GradFOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::GradFOutput  (GRADF_NUM_OUT = 3) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_GRAD</td><td></td><td>Jacobian of the constraints.</td></tr>
<tr><td>GRADF_F</td><td></td><td>Objective function.</td></tr>
<tr><td>GRADF_G</td><td></td><td>Constraint function.</td></tr>
</table>
*/
/** \defgroup scheme_QPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::QPStruct  ( = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_STRUCT_H</td><td></td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_STRUCT_A</td><td></td><td>The matrix A: sparse, (nc x n) - product with x must be dense.</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::IntegratorInput  (INTEGRATOR_NUM_IN = 6) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td></td><td>Differential state at the initial time.</td></tr>
<tr><td>INTEGRATOR_P</td><td></td><td>Parameters.</td></tr>
<tr><td>INTEGRATOR_Z0</td><td></td><td>Initial guess for the algebraic variable.</td></tr>
<tr><td>INTEGRATOR_RX0</td><td></td><td>Backward differential state at the final time.</td></tr>
<tr><td>INTEGRATOR_RP</td><td></td><td>Backward parameter vector.</td></tr>
<tr><td>INTEGRATOR_RZ0</td><td></td><td>Initial guess for the backwards algebraic variable.</td></tr>
</table>
*/
/** \defgroup scheme_NLPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NLPInput  (NL_NUM_IN = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_X</td><td></td><td>Decision variable.</td></tr>
<tr><td>NL_P</td><td></td><td>Fixed parameter.</td></tr>
</table>
*/
/** \defgroup scheme_RDAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::RDAEOutput  (RDAE_NUM_OUT = 3) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>RDAE_ODE</td><td></td><td>Right hand side of ODE.</td></tr>
<tr><td>RDAE_ALG</td><td></td><td>Right hand side of algebraic equations.</td></tr>
<tr><td>RDAE_QUAD</td><td></td><td>Right hand side of quadratures.</td></tr>
</table>
*/
/** \defgroup scheme_NlpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NlpSolverOutput  (NLP_SOLVER_NUM_OUT = 6) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X</td><td></td><td>Decision variables at the optimal solution (nx x 1)</td></tr>
<tr><td>NLP_SOLVER_F</td><td></td><td>Cost function value at the optimal solution (1 x 1)</td></tr>
<tr><td>NLP_SOLVER_G</td><td></td><td>Constraints function at the optimal solution (ng x 1)</td></tr>
<tr><td>NLP_SOLVER_LAM_X</td><td></td><td>Lagrange multipliers for bounds on X at the solution (nx x 1)</td></tr>
<tr><td>NLP_SOLVER_LAM_G</td><td></td><td>Lagrange multipliers for bounds on G at the solution (ng x 1)</td></tr>
<tr><td>NLP_SOLVER_LAM_P</td><td></td><td>Lagrange multipliers for bounds on P at the solution (np x 1)</td></tr>
</table>
*/
/** \defgroup scheme_LinsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LinsolOutput  (LINSOL_NUM_OUT = 1) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_X</td><td></td><td>Solution to the linear system of equations.</td></tr>
</table>
*/
/** \defgroup scheme_NlpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NlpSolverInput  (NLP_SOLVER_NUM_IN = 8) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X0</td><td></td><td>Decision variables, initial guess (nx x 1)</td></tr>
<tr><td>NLP_SOLVER_P</td><td></td><td>Value of fixed parameters (np x 1)</td></tr>
<tr><td>NLP_SOLVER_LBX</td><td></td><td>Decision variables lower bound (nx x 1), default -inf.</td></tr>
<tr><td>NLP_SOLVER_UBX</td><td></td><td>Decision variables upper bound (nx x 1), default +inf.</td></tr>
<tr><td>NLP_SOLVER_LBG</td><td></td><td>Constraints lower bound (ng x 1), default -inf.</td></tr>
<tr><td>NLP_SOLVER_UBG</td><td></td><td>Constraints upper bound (ng x 1), default +inf.</td></tr>
<tr><td>NLP_SOLVER_LAM_X0</td><td></td><td>Lagrange multipliers for bounds on X, initial guess (nx x 1)</td></tr>
<tr><td>NLP_SOLVER_LAM_G0</td><td></td><td>Lagrange multipliers for bounds on G, initial guess (ng x 1)</td></tr>
</table>
*/
/** \defgroup scheme_JacGOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::JacGOutput  (JACG_NUM_OUT = 3) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_JAC</td><td></td><td>Jacobian of the constraints.</td></tr>
<tr><td>JACG_F</td><td></td><td>Objective function.</td></tr>
<tr><td>JACG_G</td><td></td><td>Constraint function.</td></tr>
</table>
*/
/** \defgroup scheme_JacGInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::JacGInput  (JACG_NUM_IN = 2) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_X</td><td></td><td>Decision variable.</td></tr>
<tr><td>JACG_P</td><td></td><td>Fixed parameter.</td></tr>
</table>
*/
/** \defgroup scheme_HessLagInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::HessLagInput  (HESSLAG_NUM_IN = 4) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HESSLAG_X</td><td></td><td>Decision variable.</td></tr>
<tr><td>HESSLAG_P</td><td></td><td>Fixed parameter.</td></tr>
<tr><td>HESSLAG_LAM_F</td><td></td><td>Multiplier for f. Just a scalar factor for the objective that the NLP solver might use to scale the objective</td></tr>
<tr><td>HESSLAG_LAM_G</td><td></td><td>Multiplier for g.</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CvodesInterface
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::IdasInterface
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SundialsInterface
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
