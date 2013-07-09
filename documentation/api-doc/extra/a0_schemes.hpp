/** \defgroup scheme_IntegratorOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td>Differential state at the final time [xf].</td></tr>
<tr><td>INTEGRATOR_QF</td><td>Quadrature state at the final time [qf].</td></tr>
<tr><td>INTEGRATOR_RXF</td><td>Backward differential state at the initial time [rxf].</td></tr>
<tr><td>INTEGRATOR_RQF</td><td>Backward quadrature state at the initial time [rqf].</td></tr>
</table>
*/
/** \defgroup scheme_QCQPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::QCQPSolverInput  (QCQP_SOLVER_NUM_IN = 12) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_H</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. [h].</td></tr>
<tr><td>QCQP_SOLVER_G</td><td>The vector g: dense, (n x 1) [g].</td></tr>
<tr><td>QCQP_SOLVER_P</td><td>The vertical stack of all Pi. Each Pi is sparse (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. [p].</td></tr>
<tr><td>QCQP_SOLVER_Q</td><td>The vertical stack of all qi: dense, (nq n x 1) [q].</td></tr>
<tr><td>QCQP_SOLVER_R</td><td>The vertical stack of all scalars ri (nq x 1) [r].</td></tr>
<tr><td>QCQP_SOLVER_A</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. [a].</td></tr>
<tr><td>QCQP_SOLVER_LBA</td><td>dense, (nc x 1) [lba]</td></tr>
<tr><td>QCQP_SOLVER_UBA</td><td>dense, (nc x 1) [uba]</td></tr>
<tr><td>QCQP_SOLVER_LBX</td><td>dense, (n x 1) [lbx]</td></tr>
<tr><td>QCQP_SOLVER_UBX</td><td>dense, (n x 1) [ubx]</td></tr>
<tr><td>QCQP_SOLVER_X0</td><td>dense, (n x 1) [x0]</td></tr>
<tr><td>QCQP_SOLVER_LAM_X0</td><td>dense [lam_x0]</td></tr>
</table>
*/
/** \defgroup scheme_HessLagOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::HessLagOutput  (HESSLAG_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>HESSLAG_HESS</td><td>Hessian of the Lagrangian [hess].</td></tr>
<tr><td>HESSLAG_F</td><td>Objective function [f].</td></tr>
<tr><td>HESSLAG_G</td><td>Constraint function [g].</td></tr>
<tr><td>HESSLAG_GRAD_X</td><td>Gradient of the Lagrangian with respect to x [grad_x].</td></tr>
<tr><td>HESSLAG_GRAD_P</td><td>Gradient of the Lagrangian with respect to p [grad_p].</td></tr>
</table>
*/
/** \defgroup scheme_LinsolInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::LinsolInput  (LINSOL_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LINSOL_A</td><td>The square matrix A: sparse, (n x n). [A].</td></tr>
<tr><td>LINSOL_B</td><td>The right-hand-side matrix b: dense, (n x m) [B].</td></tr>
<tr><td>LINSOL_T</td><td></td></tr>
</table>
*/
/** \defgroup scheme_SOCPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SOCPOutput  (SOCP_SOLVER_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_X</td><td>The primal solution (n x 1) [x].</td></tr>
<tr><td>SOCP_SOLVER_COST</td><td>The primal optimal cost (1 x 1) [cost].</td></tr>
<tr><td>SOCP_SOLVER_LAM_A</td><td>The dual solution corresponding to the linear constraints (nc x 1) [lam_a].</td></tr>
<tr><td>SOCP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds (n x 1) [lam_x].</td></tr>
</table>
*/
/** \defgroup scheme_NLPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::NLPSolverOutput  (NLP_SOLVER_NUM_OUT = 6) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X</td><td>Decision variables at the optimal solution (nx x 1) [x].</td></tr>
<tr><td>NLP_SOLVER_F</td><td>Cost function value at the optimal solution (1 x 1) [f].</td></tr>
<tr><td>NLP_SOLVER_G</td><td>Constraints function at the optimal solution (ng x 1) [g].</td></tr>
<tr><td>NLP_SOLVER_LAM_X</td><td>Lagrange multipliers for bounds on X at the solution (nx x 1) [lam_x].</td></tr>
<tr><td>NLP_SOLVER_LAM_G</td><td>Lagrange multipliers for bounds on G at the solution (ng x 1) [lam_g].</td></tr>
<tr><td>NLP_SOLVER_LAM_P</td><td>Lagrange multipliers for bounds on P at the solution (np x 1) [lam_p].</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_Input
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ACADO_Input  (ACADO_NUM_IN = 17) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>ACADO_X_GUESS</td><td>Initial guess for x (default: 0) [x_guess].</td></tr>
<tr><td>ACADO_U_GUESS</td><td>Initial guess for u (default: 0) [u_guess].</td></tr>
<tr><td>ACADO_P_GUESS</td><td>Initial guess for p (default: 0) [p_guess].</td></tr>
<tr><td>ACADO_LBX</td><td>Lower bound on x (default: -infinity) [lbx].</td></tr>
<tr><td>ACADO_UBX</td><td>Upper bound on x (default: infinity) [ubx].</td></tr>
<tr><td>ACADO_LBX0</td><td>Lower bound on x0 (default: -infinity) [lbx0].</td></tr>
<tr><td>ACADO_UBX0</td><td>Upper bound on x0 (default: infinity) [ubx0].</td></tr>
<tr><td>ACADO_LBXF</td><td>Lower bound on xf (default: -infinity) [lbxf].</td></tr>
<tr><td>ACADO_UBXF</td><td>Upper bound on xf (default: infinity) [ubxf].</td></tr>
<tr><td>ACADO_LBU</td><td>Lower bound on u (default: -infinity) [lbu].</td></tr>
<tr><td>ACADO_UBU</td><td>Upper bound on u (default: infinity) [ubu].</td></tr>
<tr><td>ACADO_LBP</td><td>Lower bound on p (default: -infinity) [lbp].</td></tr>
<tr><td>ACADO_UBP</td><td>Upper bound on p (default: infinity) [ubp].</td></tr>
<tr><td>ACADO_LBC</td><td>Lower bound on the path constraint function (default: -infinity) [lbc].</td></tr>
<tr><td>ACADO_UBC</td><td>Upper bound on the path constraint function (default: infinity) [ubc].</td></tr>
<tr><td>ACADO_LBR</td><td>Lower bound on the initial constraint function (default: 0) [lbr].</td></tr>
<tr><td>ACADO_UBR</td><td>Upper bound on the initial constraint function (default: 0) [ubr].</td></tr>
</table>
*/
/** \defgroup scheme_SDPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SDPInput  (SDP_SOLVER_NUM_IN = 8) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDP_SOLVER_F</td><td>The vertical stack of all matrices F_i: ( nm x m) [f].</td></tr>
<tr><td>SDP_SOLVER_C</td><td>The vector c: ( n x 1) [c].</td></tr>
<tr><td>SDP_SOLVER_G</td><td>The matrix G: ( m x m) [g].</td></tr>
<tr><td>SDP_SOLVER_A</td><td>The matrix A: ( nc x n) [a].</td></tr>
<tr><td>SDP_SOLVER_LBA</td><td>Lower bounds on Ax ( nc x 1) [lba].</td></tr>
<tr><td>SDP_SOLVER_UBA</td><td>Upper bounds on Ax ( nc x 1) [uba].</td></tr>
<tr><td>SDP_SOLVER_LBX</td><td>Lower bounds on x ( n x 1 ) [lbx].</td></tr>
<tr><td>SDP_SOLVER_UBX</td><td>Upper bounds on x ( n x 1 ) [ubx].</td></tr>
</table>
*/
/** \defgroup scheme_LPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::LPSolverInput  (LP_SOLVER_NUM_IN = 6) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LP_SOLVER_C</td><td>The vector c: dense (n x 1) [c].</td></tr>
<tr><td>LP_SOLVER_A</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. [a].</td></tr>
<tr><td>LP_SOLVER_LBA</td><td>dense, (nc x 1) [lba]</td></tr>
<tr><td>LP_SOLVER_UBA</td><td>dense, (nc x 1) [uba]</td></tr>
<tr><td>LP_SOLVER_LBX</td><td>dense, (n x 1) [lbx]</td></tr>
<tr><td>LP_SOLVER_UBX</td><td>dense, (n x 1) [ubx]</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_FCN_Input
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ACADO_FCN_Input  (ACADO_FCN_NUM_IN = 6) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>ACADO_FCN_T</td><td>Time [t].</td></tr>
<tr><td>ACADO_FCN_XD</td><td>Differential state [xd].</td></tr>
<tr><td>ACADO_FCN_XA</td><td>Algebraic state [xa].</td></tr>
<tr><td>ACADO_FCN_U</td><td>Control input [u].</td></tr>
<tr><td>ACADO_FCN_P</td><td>Parameter [p].</td></tr>
<tr><td>ACADO_FCN_XDOT</td><td>Differential state derivative [xdot].</td></tr>
</table>
*/
/** \defgroup scheme_RDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::RDAEInput  (RDAE_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>RDAE_RX</td><td>Backward differential state [rx].</td></tr>
<tr><td>RDAE_RZ</td><td>Backward algebraic state [rz].</td></tr>
<tr><td>RDAE_RP</td><td>Backward parameter vector [rp].</td></tr>
<tr><td>RDAE_X</td><td>Forward differential state [x].</td></tr>
<tr><td>RDAE_Z</td><td>Forward algebraic state [z].</td></tr>
<tr><td>RDAE_P</td><td>Parameter vector [p].</td></tr>
<tr><td>RDAE_T</td><td>Explicit time dependence [t].</td></tr>
</table>
*/
/** \defgroup scheme_MUSCOD_FCN_Output
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::MUSCOD_FCN_Output  (MUSCOD_FCN_NUM_OUT = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MUSCOD_FCN_RHS</td><td></td></tr>
<tr><td>MUSCOD_FCN_RES</td><td></td></tr>
</table>
*/
/** \defgroup scheme_NLPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::NLPOutput  (NL_NUM_OUT = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NL_F</td><td>Objective function [f].</td></tr>
<tr><td>NL_G</td><td>Constraint function [g].</td></tr>
</table>
*/
/** \defgroup scheme_SOCPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SOCPInput  (SOCP_SOLVER_NUM_IN = 10) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_G</td><td>The vertical stack of all matrices Gi: ( N x n) [g].</td></tr>
<tr><td>SOCP_SOLVER_H</td><td>The vertical stack of all vectors hi: ( N x 1) [h].</td></tr>
<tr><td>SOCP_SOLVER_E</td><td>The vertical stack of all vectors ei: ( nm x 1) [e].</td></tr>
<tr><td>SOCP_SOLVER_F</td><td>The vertical stack of all scalars fi: ( m x 1) [f].</td></tr>
<tr><td>SOCP_SOLVER_C</td><td>The vector c: ( n x 1) [c].</td></tr>
<tr><td>SOCP_SOLVER_A</td><td>The matrix A: ( nc x n) [a].</td></tr>
<tr><td>SOCP_SOLVER_LBA</td><td>Lower bounds on Ax ( nc x 1) [lba].</td></tr>
<tr><td>SOCP_SOLVER_UBA</td><td>Upper bounds on Ax ( nc x 1) [uba].</td></tr>
<tr><td>SOCP_SOLVER_LBX</td><td>Lower bounds on x ( n x 1 ) [lbx].</td></tr>
<tr><td>SOCP_SOLVER_UBX</td><td>Upper bounds on x ( n x 1 ) [ubx].</td></tr>
</table>
*/
/** \defgroup scheme_SDPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SDPOutput  (SDP_SOLVER_NUM_OUT = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDP_SOLVER_X</td><td>The primal solution (n x 1) - may be used as initial guess [x].</td></tr>
<tr><td>SDP_SOLVER_P</td><td>The solution P (m x m) - may be used as initial guess [p].</td></tr>
<tr><td>SDP_SOLVER_DUAL</td><td>The dual solution (m x m) - may be used as initial guess [dual].</td></tr>
<tr><td>SDP_SOLVER_COST</td><td>The primal optimal cost (1 x 1) [cost].</td></tr>
<tr><td>SDP_SOLVER_DUAL_COST</td><td>The dual optimal cost (1 x 1) [dual_cost].</td></tr>
<tr><td>SDP_SOLVER_LAM_A</td><td>The dual solution corresponding to the linear constraints (nc x 1) [lam_a].</td></tr>
<tr><td>SDP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds (n x 1) [lam_x].</td></tr>
</table>
*/
/** \defgroup scheme_NLPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::NLPSolverInput  (NLP_SOLVER_NUM_IN = 8) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X0</td><td>Decision variables, initial guess (nx x 1) [x0].</td></tr>
<tr><td>NLP_SOLVER_P</td><td>Value of fixed parameters (np x 1) [p].</td></tr>
<tr><td>NLP_SOLVER_LBX</td><td>Decision variables lower bound (nx x 1), default -inf [lbx].</td></tr>
<tr><td>NLP_SOLVER_UBX</td><td>Decision variables upper bound (nx x 1), default +inf [ubx].</td></tr>
<tr><td>NLP_SOLVER_LBG</td><td>Constraints lower bound (ng x 1), default -inf [lbg].</td></tr>
<tr><td>NLP_SOLVER_UBG</td><td>Constraints upper bound (ng x 1), default +inf [ubg].</td></tr>
<tr><td>NLP_SOLVER_LAM_X0</td><td>Lagrange multipliers for bounds on X, initial guess (nx x 1) [lam_x0].</td></tr>
<tr><td>NLP_SOLVER_LAM_G0</td><td>Lagrange multipliers for bounds on G, initial guess (ng x 1) [lam_g0].</td></tr>
</table>
*/
/** \defgroup scheme_DAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_X</td><td>Differential state [x].</td></tr>
<tr><td>DAE_Z</td><td>Algebraic state [z].</td></tr>
<tr><td>DAE_P</td><td>Parameter [p].</td></tr>
<tr><td>DAE_T</td><td>Explicit time dependence [t].</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_Output
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::ACADO_Output  (ACADO_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>ACADO_X_OPT</td><td>Optimal states [x_opt].</td></tr>
<tr><td>ACADO_U_OPT</td><td>Optimal control inputs [u_opt].</td></tr>
<tr><td>ACADO_P_OPT</td><td>Optimal parameters [p_opt].</td></tr>
<tr><td>ACADO_COST</td><td>Optimal cost [cost].</td></tr>
</table>
*/
/** \defgroup scheme_DAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_ODE</td><td>Right hand side of the implicit ODE [ode].</td></tr>
<tr><td>DAE_ALG</td><td>Right hand side of algebraic equations [alg].</td></tr>
<tr><td>DAE_QUAD</td><td>Right hand side of quadratures equations [quad].</td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::InputOutputScheme  ( = 43) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SCHEME_ACADO_Input</td><td></td></tr>
<tr><td>SCHEME_ACADO_Output</td><td></td></tr>
<tr><td>SCHEME_ACADO_FCN_Input</td><td></td></tr>
<tr><td>SCHEME_ControlledDAEInput</td><td></td></tr>
<tr><td>SCHEME_ControlSimulatorInput</td><td></td></tr>
<tr><td>SCHEME_DAEInput</td><td></td></tr>
<tr><td>SCHEME_DAEOutput</td><td></td></tr>
<tr><td>SCHEME_RDAEInput</td><td></td></tr>
<tr><td>SCHEME_RDAEOutput</td><td></td></tr>
<tr><td>SCHEME_IntegratorInput</td><td></td></tr>
<tr><td>SCHEME_IntegratorOutput</td><td></td></tr>
<tr><td>SCHEME_LinsolInput</td><td></td></tr>
<tr><td>SCHEME_LinsolOutput</td><td></td></tr>
<tr><td>SCHEME_LPSolverInput</td><td></td></tr>
<tr><td>SCHEME_LPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_LPStruct</td><td></td></tr>
<tr><td>SCHEME_NLPInput</td><td></td></tr>
<tr><td>SCHEME_NLPOutput</td><td></td></tr>
<tr><td>SCHEME_GradFInput</td><td></td></tr>
<tr><td>SCHEME_GradFOutput</td><td></td></tr>
<tr><td>SCHEME_JacGInput</td><td></td></tr>
<tr><td>SCHEME_JacGOutput</td><td></td></tr>
<tr><td>SCHEME_HessLagInput</td><td></td></tr>
<tr><td>SCHEME_HessLagOutput</td><td></td></tr>
<tr><td>SCHEME_NLPSolverInput</td><td></td></tr>
<tr><td>SCHEME_NLPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_MayerInput</td><td></td></tr>
<tr><td>SCHEME_OCPInput</td><td></td></tr>
<tr><td>SCHEME_OCPOutput</td><td></td></tr>
<tr><td>SCHEME_QCQPSolverInput</td><td></td></tr>
<tr><td>SCHEME_QCQPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_QCQPStruct</td><td></td></tr>
<tr><td>SCHEME_QPSolverInput</td><td></td></tr>
<tr><td>SCHEME_QPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_QPStruct</td><td></td></tr>
<tr><td>SCHEME_SDPInput</td><td></td></tr>
<tr><td>SCHEME_SDPOutput</td><td></td></tr>
<tr><td>SCHEME_SDPStruct</td><td></td></tr>
<tr><td>SCHEME_SDQPInput</td><td></td></tr>
<tr><td>SCHEME_SDQPOutput</td><td></td></tr>
<tr><td>SCHEME_SDQPStruct</td><td></td></tr>
<tr><td>SCHEME_SOCPInput</td><td></td></tr>
<tr><td>SCHEME_SOCPOutput</td><td></td></tr>
<tr><td>SCHEME_SOCPStruct</td><td></td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::InputOutputScheme  ( = 43) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SCHEME_ACADO_Input</td><td></td></tr>
<tr><td>SCHEME_ACADO_Output</td><td></td></tr>
<tr><td>SCHEME_ACADO_FCN_Input</td><td></td></tr>
<tr><td>SCHEME_ControlledDAEInput</td><td></td></tr>
<tr><td>SCHEME_ControlSimulatorInput</td><td></td></tr>
<tr><td>SCHEME_DAEInput</td><td></td></tr>
<tr><td>SCHEME_DAEOutput</td><td></td></tr>
<tr><td>SCHEME_RDAEInput</td><td></td></tr>
<tr><td>SCHEME_RDAEOutput</td><td></td></tr>
<tr><td>SCHEME_IntegratorInput</td><td></td></tr>
<tr><td>SCHEME_IntegratorOutput</td><td></td></tr>
<tr><td>SCHEME_LinsolInput</td><td></td></tr>
<tr><td>SCHEME_LinsolOutput</td><td></td></tr>
<tr><td>SCHEME_LPSolverInput</td><td></td></tr>
<tr><td>SCHEME_LPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_LPStruct</td><td></td></tr>
<tr><td>SCHEME_NLPInput</td><td></td></tr>
<tr><td>SCHEME_NLPOutput</td><td></td></tr>
<tr><td>SCHEME_GradFInput</td><td></td></tr>
<tr><td>SCHEME_GradFOutput</td><td></td></tr>
<tr><td>SCHEME_JacGInput</td><td></td></tr>
<tr><td>SCHEME_JacGOutput</td><td></td></tr>
<tr><td>SCHEME_HessLagInput</td><td></td></tr>
<tr><td>SCHEME_HessLagOutput</td><td></td></tr>
<tr><td>SCHEME_NLPSolverInput</td><td></td></tr>
<tr><td>SCHEME_NLPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_MayerInput</td><td></td></tr>
<tr><td>SCHEME_OCPInput</td><td></td></tr>
<tr><td>SCHEME_OCPOutput</td><td></td></tr>
<tr><td>SCHEME_QCQPSolverInput</td><td></td></tr>
<tr><td>SCHEME_QCQPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_QCQPStruct</td><td></td></tr>
<tr><td>SCHEME_QPSolverInput</td><td></td></tr>
<tr><td>SCHEME_QPSolverOutput</td><td></td></tr>
<tr><td>SCHEME_QPStruct</td><td></td></tr>
<tr><td>SCHEME_SDPInput</td><td></td></tr>
<tr><td>SCHEME_SDPOutput</td><td></td></tr>
<tr><td>SCHEME_SDPStruct</td><td></td></tr>
<tr><td>SCHEME_SDQPInput</td><td></td></tr>
<tr><td>SCHEME_SDQPOutput</td><td></td></tr>
<tr><td>SCHEME_SDQPStruct</td><td></td></tr>
<tr><td>SCHEME_SOCPInput</td><td></td></tr>
<tr><td>SCHEME_SOCPOutput</td><td></td></tr>
<tr><td>SCHEME_SOCPStruct</td><td></td></tr>
</table>
*/
/** \defgroup scheme_GradFInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::GradFInput  (GRADF_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>GRADF_X</td><td>Decision variable [x].</td></tr>
<tr><td>GRADF_P</td><td>Fixed parameter [p].</td></tr>
</table>
*/
/** \defgroup scheme_LPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::LPSolverOutput  (LP_SOLVER_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LP_SOLVER_X</td><td>The primal solution [x].</td></tr>
<tr><td>LP_SOLVER_COST</td><td>The optimal cost [cost].</td></tr>
<tr><td>LP_SOLVER_LAM_A</td><td>The dual solution corresponding to linear bounds [lam_a].</td></tr>
<tr><td>LP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds [lam_x].</td></tr>
</table>
*/
/** \defgroup scheme_GradFOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::GradFOutput  (GRADF_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>GRADF_GRAD</td><td>Jacobian of the constraints [grad].</td></tr>
<tr><td>GRADF_F</td><td>Objective function [f].</td></tr>
<tr><td>GRADF_G</td><td>Constraint function [g].</td></tr>
</table>
*/
/** \defgroup scheme_SDQPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SDQPInput  (SDQP_SOLVER_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDQP_SOLVER_H</td><td>The matrix H: sparse ( n x n) [h].</td></tr>
<tr><td>SDQP_SOLVER_C</td><td>The vector c: ( n x 1) [c].</td></tr>
<tr><td>SDQP_SOLVER_F</td><td>The vertical stack of all matrices F_i: ( nm x m) [f].</td></tr>
<tr><td>SDQP_SOLVER_G</td><td>The matrix G: ( m x m) [g].</td></tr>
<tr><td>SDQP_SOLVER_A</td><td>The matrix A: ( nc x n) [a].</td></tr>
<tr><td>SDQP_SOLVER_LBA</td><td>Lower bounds on Ax ( nc x 1) [lba].</td></tr>
<tr><td>SDQP_SOLVER_UBA</td><td>Upper bounds on Ax ( nc x 1) [uba].</td></tr>
<tr><td>SDQP_SOLVER_LBX</td><td>Lower bounds on x ( n x 1 ) [lbx].</td></tr>
<tr><td>SDQP_SOLVER_UBX</td><td>Upper bounds on x ( n x 1 ) [ubx].</td></tr>
</table>
*/
/** \defgroup scheme_OCPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::OCPInput  (OCP_NUM_IN = 13) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_LBX</td><td>States lower bounds (nx x (ns+1)) [lbx].</td></tr>
<tr><td>OCP_UBX</td><td>States upper bounds (nx x (ns+1)) [ubx].</td></tr>
<tr><td>OCP_X_INIT</td><td>States initial guess (nx x (ns+1)) [x_init].</td></tr>
<tr><td>OCP_LBU</td><td>Controls lower bounds (nu x ns) [lbu].</td></tr>
<tr><td>OCP_UBU</td><td>Controls upper bounds (nu x ns) [ubu].</td></tr>
<tr><td>OCP_U_INIT</td><td>Controls initial guess (nu x ns) [u_init].</td></tr>
<tr><td>OCP_LBP</td><td>Parameters lower bounds (np x 1) [lbp].</td></tr>
<tr><td>OCP_UBP</td><td>Parameters upper bounds (np x 1) [ubp].</td></tr>
<tr><td>OCP_P_INIT</td><td>Parameters initial guess (np x 1) [p_init].</td></tr>
<tr><td>OCP_LBH</td><td>Point constraint lower bound (nh x (ns+1)) [lbh].</td></tr>
<tr><td>OCP_UBH</td><td>Point constraint upper bound (nh x (ns+1)) [ubh].</td></tr>
<tr><td>OCP_LBG</td><td>Lower bound for the coupling constraints [lbg].</td></tr>
<tr><td>OCP_UBG</td><td>Upper bound for the coupling constraints [ubg].</td></tr>
</table>
*/
/** \defgroup scheme_QCQPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::QCQPSolverOutput  (QCQP_SOLVER_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_X</td><td>The primal solution [x].</td></tr>
<tr><td>QCQP_SOLVER_COST</td><td>The optimal cost [cost].</td></tr>
<tr><td>QCQP_SOLVER_LAM_A</td><td>The dual solution corresponding to linear bounds [lam_a].</td></tr>
<tr><td>QCQP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds [lam_x].</td></tr>
</table>
*/
/** \defgroup scheme_MayerInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::MayerInput  (MAYER_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MAYER_X</td><td>States at the end of integration (nx x 1) [x].</td></tr>
<tr><td>MAYER_P</td><td>Problem parameters (np x 1) [p].</td></tr>
</table>
*/
/** \defgroup scheme_ControlledDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ControlledDAEInput  (CONTROL_DAE_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>CONTROL_DAE_T</td><td>Global physical time. (1-by-1) [t].</td></tr>
<tr><td>CONTROL_DAE_X</td><td>State vector (dimension nx-by-1). Should have same amount of non-zeros as DAEOutput:DAE_RES [x].</td></tr>
<tr><td>CONTROL_DAE_Z</td><td>Algebraic state vector (dimension np-by-1). [z].</td></tr>
<tr><td>CONTROL_DAE_P</td><td>Parameter vector (dimension np-by-1). [p].</td></tr>
<tr><td>CONTROL_DAE_U</td><td>Control vector (dimension nu-by-1). [u].</td></tr>
<tr><td>CONTROL_DAE_U_INTERP</td><td>Control vector, linearly interpolated (dimension nu-by-1). [u_interp].</td></tr>
<tr><td>CONTROL_DAE_X_MAJOR</td><td>State vector (dimension nx-by-1) at the last major time-step [x_major].</td></tr>
<tr><td>CONTROL_DAE_T0</td><td>Time at start of control interval (1-by-1) [t0].</td></tr>
<tr><td>CONTROL_DAE_TF</td><td>Time at end of control interval (1-by-1) [tf].</td></tr>
</table>
*/
/** \defgroup scheme_NLPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::NLPInput  (NL_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NL_X</td><td>Decision variable [x].</td></tr>
<tr><td>NL_P</td><td>Fixed parameter [p].</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::IntegratorInput  (INTEGRATOR_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>Differential state at the initial time [x0].</td></tr>
<tr><td>INTEGRATOR_P</td><td>Parameters [p].</td></tr>
<tr><td>INTEGRATOR_RX0</td><td>Backward differential state at the final time [rx0].</td></tr>
<tr><td>INTEGRATOR_RP</td><td>Backward parameter vector [rp].</td></tr>
</table>
*/
/** \defgroup scheme_QPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::QPSolverInput  (QP_SOLVER_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_SOLVER_H</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. [h].</td></tr>
<tr><td>QP_SOLVER_G</td><td>The vector g: dense, (n x 1) [g].</td></tr>
<tr><td>QP_SOLVER_A</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. [a].</td></tr>
<tr><td>QP_SOLVER_LBA</td><td>dense, (nc x 1) [lba]</td></tr>
<tr><td>QP_SOLVER_UBA</td><td>dense, (nc x 1) [uba]</td></tr>
<tr><td>QP_SOLVER_LBX</td><td>dense, (n x 1) [lbx]</td></tr>
<tr><td>QP_SOLVER_UBX</td><td>dense, (n x 1) [ubx]</td></tr>
<tr><td>QP_SOLVER_X0</td><td>dense, (n x 1) [x0]</td></tr>
<tr><td>QP_SOLVER_LAM_X0</td><td>dense [lam_x0]</td></tr>
</table>
*/
/** \defgroup scheme_OCPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory [x_opt].</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory [u_opt].</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters [p_opt].</td></tr>
<tr><td>OCP_COST</td><td>Objective/cost function for optimal solution (1 x 1) [cost].</td></tr>
</table>
*/
/** \defgroup scheme_RDAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::RDAEOutput  (RDAE_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>RDAE_ODE</td><td>Right hand side of ODE. [ode].</td></tr>
<tr><td>RDAE_ALG</td><td>Right hand side of algebraic equations. [alg].</td></tr>
<tr><td>RDAE_QUAD</td><td>Right hand side of quadratures. [quad].</td></tr>
</table>
*/
/** \defgroup scheme_JacGOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::JacGOutput  (JACG_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>JACG_JAC</td><td>Jacobian of the constraints [jac].</td></tr>
<tr><td>JACG_F</td><td>Objective function [f].</td></tr>
<tr><td>JACG_G</td><td>Constraint function [g].</td></tr>
</table>
*/
/** \defgroup scheme_LinsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::LinsolOutput  (LINSOL_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LINSOL_X</td><td>Solution to the linear system of equations [X].</td></tr>
</table>
*/
/** \defgroup scheme_SDQPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SDQPOutput  (SDQP_SOLVER_NUM_OUT = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDQP_SOLVER_X</td><td>The primal solution (n x 1) - may be used as initial guess [x].</td></tr>
<tr><td>SDQP_SOLVER_P</td><td>The solution P (m x m) - may be used as initial guess [p].</td></tr>
<tr><td>SDQP_SOLVER_DUAL</td><td>The dual solution (m x m) - may be used as initial guess [dual].</td></tr>
<tr><td>SDQP_SOLVER_COST</td><td>The primal optimal cost (1 x 1) [cost].</td></tr>
<tr><td>SDQP_SOLVER_DUAL_COST</td><td>The dual optimal cost (1 x 1) [dual_cost].</td></tr>
<tr><td>SDQP_SOLVER_LAM_A</td><td>The dual solution corresponding to the linear constraints (nc x 1) [lam_a].</td></tr>
<tr><td>SDQP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds (n x 1) [lam_x].</td></tr>
</table>
*/
/** \defgroup scheme_ControlSimulatorInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ControlSimulatorInput  (CONTROLSIMULATOR_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>CONTROLSIMULATOR_X0</td><td>Differential or algebraic state at t0 (dimension nx-by-1) [x0].</td></tr>
<tr><td>CONTROLSIMULATOR_P</td><td>Parameters that are fixed over the entire horizon (dimension np-by-1) [p].</td></tr>
<tr><td>CONTROLSIMULATOR_U</td><td>Parameters that change over the integration intervals (dimension (ns-1)-by-nu) [u].</td></tr>
</table>
*/
/** \defgroup scheme_JacGInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::JacGInput  (JACG_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>JACG_X</td><td>Decision variable [x].</td></tr>
<tr><td>JACG_P</td><td>Fixed parameter [p].</td></tr>
</table>
*/
/** \defgroup scheme_HessLagInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::HessLagInput  (HESSLAG_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>HESSLAG_X</td><td>Decision variable [x].</td></tr>
<tr><td>HESSLAG_P</td><td>Fixed parameter [p].</td></tr>
<tr><td>HESSLAG_LAM_F</td><td>Multiplier for f [lam_f].</td></tr>
<tr><td>HESSLAG_LAM_G</td><td>Multiplier for g [lam_g].</td></tr>
</table>
*/
/** \defgroup scheme_MUSCOD_FCN_Input
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::MUSCOD_FCN_Input  (MUSCOD_FCN_NUM_IN = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MUSCOD_FCN_T</td><td></td></tr>
<tr><td>MUSCOD_FCN_XD</td><td></td></tr>
<tr><td>MUSCOD_FCN_XA</td><td></td></tr>
<tr><td>MUSCOD_FCN_U</td><td></td></tr>
<tr><td>MUSCOD_FCN_P</td><td></td></tr>
</table>
*/
/** \defgroup scheme_QPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::QPSolverOutput  (QP_SOLVER_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_SOLVER_X</td><td>The primal solution [x].</td></tr>
<tr><td>QP_SOLVER_COST</td><td>The optimal cost [cost].</td></tr>
<tr><td>QP_SOLVER_LAM_A</td><td>The dual solution corresponding to linear bounds [lam_a].</td></tr>
<tr><td>QP_SOLVER_LAM_X</td><td>The dual solution corresponding to simple bounds [lam_x].</td></tr>
</table>
*/
/** \class CasADi::NLPSolverInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::NLPSolver
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::QPOasesInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::QPOasesSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::CSparseInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::CSparse
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::SimulatorInternal
\n
\par
@copydoc scheme_IntegratorInput
*/
/** \class CasADi::Simulator
\n
\par
@copydoc scheme_IntegratorInput
*/
/** \class CasADi::QCQPQPInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::QCQPQPSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::SDPSDQPInternal
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/** \class CasADi::SDPSDQPSolver
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/** \class CasADi::LapackLUDenseInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::LapackLUDense
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::CVodesInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::CVodesIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::QPSolverInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::QPSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::DSDPInternal
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/** \class CasADi::DSDPSolver
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/** \class CasADi::DirectSingleShootingInternal
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::DirectSingleShooting
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::SCPgenInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::SCPgen
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::IpoptInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::IpoptSolver
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::SDPSolverInternal
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/** \class CasADi::SDPSolver
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/** \class CasADi::IdasInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::IdasIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::AcadoIntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::AcadoIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::QCQPSolverInternal
\n
\par
@copydoc scheme_QCQPSolverInput
<br/>
@copydoc scheme_QCQPSolverOutput
*/
/** \class CasADi::QCQPSolver
\n
\par
@copydoc scheme_QCQPSolverInput
<br/>
@copydoc scheme_QCQPSolverOutput
*/
/** \class CasADi::CSparseCholeskyInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::CSparseCholesky
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::CplexInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::CplexSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::LapackQRDenseInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::LapackQRDense
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::KnitroInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::KnitroSolver
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::QPLPInternal
\n
\par
@copydoc scheme_LPSolverInput
<br/>
@copydoc scheme_LPSolverOutput
*/
/** \class CasADi::QPLPSolver
\n
\par
@copydoc scheme_LPSolverInput
<br/>
@copydoc scheme_LPSolverOutput
*/
/** \class CasADi::IntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::Integrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::SQPInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::SQPMethod
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::SOCPQCQPInternal
\n
\par
@copydoc scheme_QCQPSolverInput
<br/>
@copydoc scheme_QCQPSolverOutput
*/
/** \class CasADi::SOCPQCQPSolver
\n
\par
@copydoc scheme_QCQPSolverInput
<br/>
@copydoc scheme_QCQPSolverOutput
*/
/** \class CasADi::DirectMultipleShootingInternal
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::DirectMultipleShooting
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::SymbolicQRInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::SymbolicQR
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::SOCPSolverInternal
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/** \class CasADi::SOCPSolver
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/** \class CasADi::SDQPSolverInternal
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/** \class CasADi::SDQPSolver
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/** \class CasADi::RKIntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::RKIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::NLPQPInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::NLPQPSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::OCPSolverInternal
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::OCPSolver
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::LPSolverInternal
\n
\par
@copydoc scheme_LPSolverInput
<br/>
@copydoc scheme_LPSolverOutput
*/
/** \class CasADi::LPSolver
\n
\par
@copydoc scheme_LPSolverInput
<br/>
@copydoc scheme_LPSolverOutput
*/
/** \class CasADi::OOQPInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::OOQPSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::SundialsInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::SundialsIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::ControlSimulatorInternal
\n
\par
@copydoc scheme_ControlSimulatorInput
*/
/** \class CasADi::ControlSimulator
\n
\par
@copydoc scheme_ControlSimulatorInput
*/
/** \class CasADi::LinearSolverInternal
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::LinearSolver
\n
\par
@copydoc scheme_LinsolInput
<br/>
@copydoc scheme_LinsolOutput
*/
/** \class CasADi::CollocationIntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::CollocationIntegrator
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::SDPSOCPInternal
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/** \class CasADi::SDPSOCPSolver
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/** \class CasADi::DirectCollocationInternal
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::DirectCollocation
\n
\par
@copydoc scheme_OCPInput
<br/>
@copydoc scheme_OCPOutput
*/
/** \class CasADi::WorhpInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::WorhpSolver
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
