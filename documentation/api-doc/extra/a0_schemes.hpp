/** \defgroup scheme_IntegratorOutput
<table>
<caption>Output scheme: CasADi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td>Differential state at the final time [xf].</td></tr>
<tr><td>INTEGRATOR_QF</td><td>Quadrature state at the final time [qf].</td></tr>
<tr><td>INTEGRATOR_RXF</td><td>Backward differential state at the initial time [rxf].</td></tr>
<tr><td>INTEGRATOR_RQF</td><td>Backward quadrature state at the initial time [rqf].</td></tr>
</table>
*/
/** \defgroup scheme_SDPOutput
<table>
<caption>Output scheme: CasADi::SDPOutput  (SDP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDP_PRIMAL</td><td>The primal solution (m x 1) - may be used as initial guess [primal].</td></tr>
<tr><td>SDP_PRIMAL_P</td><td>The solution P (n x n) - may be used as initial guess [p].</td></tr>
<tr><td>SDP_DUAL</td><td>The dual solution (n x n) - may be used as initial guess [dual].</td></tr>
<tr><td>SDP_PRIMAL_COST</td><td>The primal optimal cost (1 x 1) [primal_cost].</td></tr>
<tr><td>SDP_DUAL_COST</td><td>The dual optimal cost (1 x 1) [dual_cost].</td></tr>
</table>
*/
/** \defgroup scheme_QPInput
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. [h].</td></tr>
<tr><td>QP_G</td><td>The vector G: dense, (nx x 1) [g].</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense. [a].</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1) [lba]</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1) [uba]</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1) [lbx]</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1) [ubx]</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1) [x_init]</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td>dense [lambda_init]</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_Input
<table>
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
<table>
<caption>Input scheme: CasADi::SDPInput  (SDP_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>SDP_A</td><td>The vertical stack of all matrices A_i: ( nm x n) [a].</td></tr>
<tr><td>SDP_B</td><td>The vector b: ( m x 1) [b].</td></tr>
<tr><td>SDP_C</td><td>The matrix C: ( n x n) [c].</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_FCN_Input
<table>
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
<table>
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
/** \defgroup scheme_LOFunInputs
<table>
<caption>Input scheme: CasADi::LOFunInputs  (LO_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LO_U</td><td></td></tr>
<tr><td>LO_LAMBDA</td><td></td></tr>
</table>
*/
/** \defgroup scheme_QPOutput
<table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_PRIMAL</td><td>The primal solution [primal].</td></tr>
<tr><td>QP_COST</td><td>The optimal cost [cost].</td></tr>
<tr><td>QP_LAMBDA_A</td><td>The dual solution corresponding to linear bounds [lambda_a].</td></tr>
<tr><td>QP_LAMBDA_X</td><td>The dual solution corresponding to simple bounds [lambda_x].</td></tr>
</table>
*/
/** \defgroup scheme_NLPOutput
<table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_SOLVER_NUM_OUT = 6) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X</td><td>Decision variables for optimal solution (nx x 1) [x].</td></tr>
<tr><td>NLP_SOLVER_F</td><td>Objective/cost function for optimal solution (1 x 1) [f].</td></tr>
<tr><td>NLP_SOLVER_LAM_G</td><td>Lagrange multipliers associated with G at the solution (ng x 1) [lam_g].</td></tr>
<tr><td>NLP_SOLVER_LAM_X</td><td>Lagrange multipliers associated with bounds on X at the solution (nx x 1) [lam_x].</td></tr>
<tr><td>NLP_SOLVER_LAM_P</td><td>Lagrange multipliers associated with the parameters (np x 1) [lam_p].</td></tr>
<tr><td>NLP_SOLVER_G</td><td>The constraints evaluated at the optimal solution (ng x 1) [g].</td></tr>
</table>
*/
/** \defgroup scheme_DAEInput
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_X</td><td>Differential state [x].</td></tr>
<tr><td>DAE_Z</td><td>Algebraic state [z].</td></tr>
<tr><td>DAE_P</td><td>Parameter [p].</td></tr>
<tr><td>DAE_T</td><td>Explicit time dependence [t].</td></tr>
</table>
*/
/** \defgroup scheme_ACADO_Output
<table>
<caption>Output scheme: CasADi::ACADO_Output  (ACADO_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>ACADO_X_OPT</td><td>Optimal states [x_opt].</td></tr>
<tr><td>ACADO_U_OPT</td><td>Optimal control inputs [u_opt].</td></tr>
<tr><td>ACADO_P_OPT</td><td>Optimal parameters [p_opt].</td></tr>
<tr><td>ACADO_COST</td><td>Optimal cost [cost].</td></tr>
</table>
*/
/** \defgroup scheme_DAEOutput
<table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_ODE</td><td>Right hand side of the implicit ODE [ode].</td></tr>
<tr><td>DAE_ALG</td><td>Right hand side of algebraic equations [alg].</td></tr>
<tr><td>DAE_QUAD</td><td>Right hand side of quadratures equations [quad].</td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<table>
<caption>Input scheme: CasADi::InputOutputScheme  ( = 20) </caption>
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
<tr><td>SCHEME_NLPInput</td><td></td></tr>
<tr><td>SCHEME_NLPOutput</td><td></td></tr>
<tr><td>SCHEME_MayerInput</td><td></td></tr>
<tr><td>SCHEME_OCPInput</td><td></td></tr>
<tr><td>SCHEME_OCPOutput</td><td></td></tr>
<tr><td>SCHEME_QPInput</td><td></td></tr>
<tr><td>SCHEME_QPOutput</td><td></td></tr>
<tr><td>SCHEME_SDPInput</td><td></td></tr>
<tr><td>SCHEME_SDPOutput</td><td></td></tr>
<tr><td>SCHEME_unknown</td><td></td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<table>
<caption>Output scheme: CasADi::InputOutputScheme  ( = 20) </caption>
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
<tr><td>SCHEME_NLPInput</td><td></td></tr>
<tr><td>SCHEME_NLPOutput</td><td></td></tr>
<tr><td>SCHEME_MayerInput</td><td></td></tr>
<tr><td>SCHEME_OCPInput</td><td></td></tr>
<tr><td>SCHEME_OCPOutput</td><td></td></tr>
<tr><td>SCHEME_QPInput</td><td></td></tr>
<tr><td>SCHEME_QPOutput</td><td></td></tr>
<tr><td>SCHEME_SDPInput</td><td></td></tr>
<tr><td>SCHEME_SDPOutput</td><td></td></tr>
<tr><td>SCHEME_unknown</td><td></td></tr>
</table>
*/
/** \defgroup scheme_MayerInput
<table>
<caption>Input scheme: CasADi::MayerInput  (MAYER_NUM_IN = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MAYER_X</td><td>States at the end of integration (nx x 1) [x].</td></tr>
<tr><td>MAYER_P</td><td>Problem parameters (np x 1) [p].</td></tr>
</table>
*/
/** \defgroup scheme_ControlledDAEInput
<table>
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
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_SOLVER_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X0</td><td>Decision variables initial guess (nx x 1) [x0].</td></tr>
<tr><td>NLP_SOLVER_LBX</td><td>Decision variables lower bound (nx x 1), default -inf [lbx].</td></tr>
<tr><td>NLP_SOLVER_UBX</td><td>Decision variables upper bound (nx x 1), default +inf [ubx].</td></tr>
<tr><td>NLP_SOLVER_LBG</td><td>Constraints lower bound (ng x 1), default -inf [lbg].</td></tr>
<tr><td>NLP_SOLVER_UBG</td><td>Constraints upper bound (ng x 1), default +inf [ubg].</td></tr>
<tr><td>NLP_SOLVER_LAM_G0</td><td>Lagrange multipliers associated with G, initial guess (ng x 1) [lam_g0].</td></tr>
<tr><td>NLP_SOLVER_P</td><td>Parameters on which the objective and constraints might depend (np x 1) [p].</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<table>
<caption>Input scheme: CasADi::IntegratorInput  (INTEGRATOR_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>Differential state at the initial time [x0].</td></tr>
<tr><td>INTEGRATOR_P</td><td>Parameters [p].</td></tr>
<tr><td>INTEGRATOR_RX0</td><td>Backward differential state at the final time [rx0].</td></tr>
<tr><td>INTEGRATOR_RP</td><td>Backward parameter vector [rp].</td></tr>
</table>
*/
/** \defgroup scheme_OCPOutput
<table>
<caption>Output scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory [x_opt].</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory [u_opt].</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters [p_opt].</td></tr>
<tr><td>OCP_COST</td><td>Objective/cost function for optimal solution (1 x 1) [cost].</td></tr>
</table>
*/
/** \defgroup scheme_RDAEOutput
<table>
<caption>Output scheme: CasADi::RDAEOutput  (RDAE_NUM_OUT = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>RDAE_ODE</td><td>Right hand side of ODE. [ode].</td></tr>
<tr><td>RDAE_ALG</td><td>Right hand side of algebraic equations. [alg].</td></tr>
<tr><td>RDAE_QUAD</td><td>Right hand side of quadratures. [quad].</td></tr>
</table>
*/
/** \defgroup scheme_MUSCOD_FCN_Input
<table>
<caption>Input scheme: CasADi::MUSCOD_FCN_Input  (MUSCOD_FCN_NUM_IN = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MUSCOD_FCN_T</td><td></td></tr>
<tr><td>MUSCOD_FCN_XD</td><td></td></tr>
<tr><td>MUSCOD_FCN_XA</td><td></td></tr>
<tr><td>MUSCOD_FCN_U</td><td></td></tr>
<tr><td>MUSCOD_FCN_P</td><td></td></tr>
</table>
*/
/** \defgroup scheme_OCPInput
<table>
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
/** \defgroup scheme_ControlSimulatorInput
<table>
<caption>Input scheme: CasADi::ControlSimulatorInput  (CONTROLSIMULATOR_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>CONTROLSIMULATOR_X0</td><td>Differential or algebraic state at t0 (dimension nx-by-1) [x0].</td></tr>
<tr><td>CONTROLSIMULATOR_P</td><td>Parameters that are fixed over the entire horizon (dimension np-by-1) [p].</td></tr>
<tr><td>CONTROLSIMULATOR_U</td><td>Parameters that change over the integration intervals (dimension (ns-1)-by-nu) [u].</td></tr>
</table>
*/
/** \defgroup scheme_MUSCOD_FCN_Output
<table>
<caption>Output scheme: CasADi::MUSCOD_FCN_Output  (MUSCOD_FCN_NUM_OUT = 2) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>MUSCOD_FCN_RHS</td><td></td></tr>
<tr><td>MUSCOD_FCN_RES</td><td></td></tr>
</table>
*/
/** \defgroup scheme_LOFunOutputs
<table>
<caption>Output scheme: CasADi::LOFunOutputs  (LO_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>LO_OBJRES</td><td></td></tr>
<tr><td>LO_EQ</td><td></td></tr>
<tr><td>LO_INEQ</td><td></td></tr>
<tr><td>LO_OBJ</td><td></td></tr>
<tr><td>LO_LAGFCN</td><td></td></tr>
</table>
*/
/** \class CasADi::LiftoptInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::LiftoptSolver
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::NLPSolverInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::NLPSolver
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::QPOasesInternal
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::QPOasesSolver
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::IPInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::IPMethod
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::WorhpInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::WorhpSolver
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
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
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::QPSolver
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
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
/** \class CasADi::LiftedSQPInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::LiftedSQP
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::SCPgenInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::SCPgen
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::OOQPInternal
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::OOQPSolver
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::IpoptInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::IpoptSolver
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
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
/** \class CasADi::CplexInternal
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::CplexSolver
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::KnitroInternal
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::KnitroSolver
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
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
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
*/
/** \class CasADi::SQPMethod
\n
\par
@copydoc scheme_NLPInput
<br/>
@copydoc scheme_NLPOutput
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
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
*/
/** \class CasADi::NLPQPSolver
\n
\par
@copydoc scheme_QPInput
<br/>
@copydoc scheme_QPOutput
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
