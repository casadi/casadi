/** \defgroup scheme_IntegratorOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 6) [integratorOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td>xf</td><td>Differential state at the final time .</td></tr>
<tr><td>INTEGRATOR_QF</td><td>qf</td><td>Quadrature state at the final time .</td></tr>
<tr><td>INTEGRATOR_ZF</td><td>zf</td><td>Algebraic variable at the final time .</td></tr>
<tr><td>INTEGRATOR_RXF</td><td>rxf</td><td>Backward differential state at the initial time .</td></tr>
<tr><td>INTEGRATOR_RQF</td><td>rqf</td><td>Backward quadrature state at the initial time .</td></tr>
<tr><td>INTEGRATOR_RZF</td><td>rzf</td><td>Backward algebraic variable at the initial time .</td></tr>
</table>
*/
/** \defgroup scheme_HessLagOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::HessLagOutput  (HESSLAG_NUM_OUT = 5) [hessLagOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HESSLAG_HESS</td><td>hess</td><td>Hessian of the Lagrangian .</td></tr>
<tr><td>HESSLAG_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>HESSLAG_G</td><td>g</td><td>Constraint function .</td></tr>
<tr><td>HESSLAG_GRAD_X</td><td>grad_x</td><td>Gradient of the Lagrangian with respect to x .</td></tr>
<tr><td>HESSLAG_GRAD_P</td><td>grad_p</td><td>Gradient of the Lagrangian with respect to p .</td></tr>
</table>
*/
/** \defgroup scheme_LinsolInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::LinsolInput  (LINSOL_NUM_IN = 2) [linsolIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_A</td><td>A</td><td>The square matrix A: sparse, (n x n). .</td></tr>
<tr><td>LINSOL_B</td><td>B</td><td>The right-hand-side matrix b: dense, (n x m) .</td></tr>
</table>
*/
/** \defgroup scheme_QpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::QpSolverInput  (QP_SOLVER_NUM_IN = 9) [qpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>QP_SOLVER_G</td><td>g</td><td>The vector g: dense, (n x 1) .</td></tr>
<tr><td>QP_SOLVER_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
<tr><td>QP_SOLVER_LBA</td><td>lba</td><td>dense, (nc x 1) </td></tr>
<tr><td>QP_SOLVER_UBA</td><td>uba</td><td>dense, (nc x 1) </td></tr>
<tr><td>QP_SOLVER_LBX</td><td>lbx</td><td>dense, (n x 1) </td></tr>
<tr><td>QP_SOLVER_UBX</td><td>ubx</td><td>dense, (n x 1) </td></tr>
<tr><td>QP_SOLVER_X0</td><td>x0</td><td>dense, (n x 1) </td></tr>
<tr><td>QP_SOLVER_LAM_X0</td><td>lam_x0</td><td>dense </td></tr>
</table>
*/
/** \defgroup scheme_NlpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NlpSolverOutput  (NLP_SOLVER_NUM_OUT = 6) [nlpSolverOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X</td><td>x</td><td>Decision variables at the optimal solution (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_F</td><td>f</td><td>Cost function value at the optimal solution (1 x 1) .</td></tr>
<tr><td>NLP_SOLVER_G</td><td>g</td><td>Constraints function at the optimal solution (ng x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_X</td><td>lam_x</td><td>Lagrange multipliers for bounds on X at the solution (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_G</td><td>lam_g</td><td>Lagrange multipliers for bounds on G at the solution (ng x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_P</td><td>lam_p</td><td>Lagrange multipliers for bounds on P at the solution (np x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_RDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::RDAEInput  (RDAE_NUM_IN = 7) [rdaeIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>RDAE_RX</td><td>rx</td><td>Backward differential state .</td></tr>
<tr><td>RDAE_RZ</td><td>rz</td><td>Backward algebraic state .</td></tr>
<tr><td>RDAE_RP</td><td>rp</td><td>Backward parameter vector .</td></tr>
<tr><td>RDAE_X</td><td>x</td><td>Forward differential state .</td></tr>
<tr><td>RDAE_Z</td><td>z</td><td>Forward algebraic state .</td></tr>
<tr><td>RDAE_P</td><td>p</td><td>Parameter vector .</td></tr>
<tr><td>RDAE_T</td><td>t</td><td>Explicit time dependence .</td></tr>
</table>
*/
/** \defgroup scheme_NLPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NLPOutput  (NL_NUM_OUT = 2) [nlpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>NL_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_DAEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::DAEInput  (DAE_NUM_IN = 4) [daeIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_X</td><td>x</td><td>Differential state .</td></tr>
<tr><td>DAE_Z</td><td>z</td><td>Algebraic state .</td></tr>
<tr><td>DAE_P</td><td>p</td><td>Parameter .</td></tr>
<tr><td>DAE_T</td><td>t</td><td>Explicit time dependence .</td></tr>
</table>
*/
/** \defgroup scheme_DAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::DAEOutput  (DAE_NUM_OUT = 3) [daeOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_ODE</td><td>ode</td><td>Right hand side of the implicit ODE .</td></tr>
<tr><td>DAE_ALG</td><td>alg</td><td>Right hand side of algebraic equations .</td></tr>
<tr><td>DAE_QUAD</td><td>quad</td><td>Right hand side of quadratures equations .</td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Input scheme: casadi::InputOutputScheme  ( = 20) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Output scheme: casadi::InputOutputScheme  ( = 20) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
</table>
*/
/** \defgroup scheme_GradFInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::GradFInput  (GRADF_NUM_IN = 2) [gradFIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>GRADF_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_GradFOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::GradFOutput  (GRADF_NUM_OUT = 3) [gradFOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_GRAD</td><td>grad</td><td>Jacobian of the constraints .</td></tr>
<tr><td>GRADF_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>GRADF_G</td><td>g</td><td>Constraint function .</td></tr>
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
/** \defgroup scheme_NLPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NLPInput  (NL_NUM_IN = 2) [nlpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>NL_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::IntegratorInput  (INTEGRATOR_NUM_IN = 6) [integratorIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>x0</td><td>Differential state at the initial time .</td></tr>
<tr><td>INTEGRATOR_P</td><td>p</td><td>Parameters .</td></tr>
<tr><td>INTEGRATOR_Z0</td><td>z0</td><td>Initial guess for the algebraic variable .</td></tr>
<tr><td>INTEGRATOR_RX0</td><td>rx0</td><td>Backward differential state at the final time .</td></tr>
<tr><td>INTEGRATOR_RP</td><td>rp</td><td>Backward parameter vector .</td></tr>
<tr><td>INTEGRATOR_RZ0</td><td>rz0</td><td>Initial guess for the backwards algebraic variable .</td></tr>
</table>
*/
/** \defgroup scheme_StabilizedQpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::StabilizedQpSolverInput  (STABILIZED_QP_SOLVER_NUM_IN = 12) [stabilizedQpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>STABILIZED_QP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>STABILIZED_QP_SOLVER_G</td><td>g</td><td>The vector g: dense, (n x 1) .</td></tr>
<tr><td>STABILIZED_QP_SOLVER_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
<tr><td>STABILIZED_QP_SOLVER_LBA</td><td>lba</td><td>dense, (nc x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_UBA</td><td>uba</td><td>dense, (nc x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_LBX</td><td>lbx</td><td>dense, (n x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_UBX</td><td>ubx</td><td>dense, (n x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_X0</td><td>x0</td><td>dense, (n x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_LAM_X0</td><td>lam_x0</td><td>dense </td></tr>
<tr><td>STABILIZED_QP_SOLVER_MUR</td><td>muR</td><td>dense (1 x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_MUE</td><td>muE</td><td>dense (nc x 1) </td></tr>
<tr><td>STABILIZED_QP_SOLVER_MU</td><td>mu</td><td>dense (nc x 1) </td></tr>
</table>
*/
/** \defgroup scheme_RDAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::RDAEOutput  (RDAE_NUM_OUT = 3) [rdaeOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>RDAE_ODE</td><td>ode</td><td>Right hand side of ODE. .</td></tr>
<tr><td>RDAE_ALG</td><td>alg</td><td>Right hand side of algebraic equations. .</td></tr>
<tr><td>RDAE_QUAD</td><td>quad</td><td>Right hand side of quadratures. .</td></tr>
</table>
*/
/** \defgroup scheme_QpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::QpSolverOutput  (QP_SOLVER_NUM_OUT = 4) [qpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>QP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>QP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>QP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
</table>
*/
/** \defgroup scheme_LinsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LinsolOutput  (LINSOL_NUM_OUT = 1) [linsolOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_X</td><td>X</td><td>Solution to the linear system of equations .</td></tr>
</table>
*/
/** \defgroup scheme_NlpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NlpSolverInput  (NLP_SOLVER_NUM_IN = 8) [nlpSolverIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X0</td><td>x0</td><td>Decision variables, initial guess (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_P</td><td>p</td><td>Value of fixed parameters (np x 1) .</td></tr>
<tr><td>NLP_SOLVER_LBX</td><td>lbx</td><td>Decision variables lower bound (nx x 1), default -inf .</td></tr>
<tr><td>NLP_SOLVER_UBX</td><td>ubx</td><td>Decision variables upper bound (nx x 1), default +inf .</td></tr>
<tr><td>NLP_SOLVER_LBG</td><td>lbg</td><td>Constraints lower bound (ng x 1), default -inf .</td></tr>
<tr><td>NLP_SOLVER_UBG</td><td>ubg</td><td>Constraints upper bound (ng x 1), default +inf .</td></tr>
<tr><td>NLP_SOLVER_LAM_X0</td><td>lam_x0</td><td>Lagrange multipliers for bounds on X, initial guess (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_G0</td><td>lam_g0</td><td>Lagrange multipliers for bounds on G, initial guess (ng x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_JacGOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::JacGOutput  (JACG_NUM_OUT = 3) [jacGOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_JAC</td><td>jac</td><td>Jacobian of the constraints .</td></tr>
<tr><td>JACG_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>JACG_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_JacGInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::JacGInput  (JACG_NUM_IN = 2) [jacGIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>JACG_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_HessLagInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::HessLagInput  (HESSLAG_NUM_IN = 4) [hessLagIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HESSLAG_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>HESSLAG_P</td><td>p</td><td>Fixed parameter .</td></tr>
<tr><td>HESSLAG_LAM_F</td><td>lam_f</td><td>Multiplier for f. Just a scalar factor for the objective that the NLP solver might use to scale the objective. </td></tr>
<tr><td>HESSLAG_LAM_G</td><td>lam_g</td><td>Multiplier for g .</td></tr>
</table>
*/
/// \cond INTERNAL
/** \class casadi::CSparseCholeskyInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::CollocationIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::CplexInterface
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::CsparseInterface
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
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
/** \class casadi::FixedStepIntegrator
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
/** \class casadi::ImplicitFixedStepIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::IntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/** \addtogroup general_Integrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \cond INTERNAL
/** \class casadi::IpoptInterface
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::KnitroInterface
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LapackLuDense
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LapackQrDense
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LinearSolverInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
/** \addtogroup general_LinearSolver
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \cond INTERNAL
/** \class casadi::NlpSolverInternal
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/** \addtogroup general_NlpSolver
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \cond INTERNAL
/** \class casadi::OldCollocationIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::OoqpInterface
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::QpSolverInternal
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/** \addtogroup general_QpSolver
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \cond INTERNAL
/** \class casadi::QpToNlp
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::QpoasesInterface
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::RkIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::Scpgen
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimulatorInternal
\n
\par
@copydoc scheme_IntegratorInput
*/
/// \endcond
/** \class casadi::Simulator
\n
\par
@copydoc scheme_IntegratorInput
*/
/// \cond INTERNAL
/** \class casadi::SnoptInterface
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SqicInterface
\n
\par
@copydoc scheme_QpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::Sqpmethod
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::StabilizedQpSolverInternal
\n
\par
@copydoc scheme_StabilizedQpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/** \addtogroup general_StabilizedQpSolver
\n
\par
@copydoc scheme_StabilizedQpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \cond INTERNAL
/** \class casadi::StabilizedQpToQp
\n
\par
@copydoc scheme_StabilizedQpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::StabilizedSqicInterface
\n
\par
@copydoc scheme_StabilizedQpSolverInput
<br/>
@copydoc scheme_QpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::StabilizedSqp
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
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
/// \cond INTERNAL
/** \class casadi::SymbolicQr
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::WorhpInterface
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
