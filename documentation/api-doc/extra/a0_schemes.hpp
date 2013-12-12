/** \defgroup scheme_IntegratorOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 5) [integratorOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td>xf</td><td>Differential state at the final time .</td></tr>
<tr><td>INTEGRATOR_QF</td><td>qf</td><td>Quadrature state at the final time .</td></tr>
<tr><td>INTEGRATOR_RXF</td><td>rxf</td><td>Backward differential state at the initial time .</td></tr>
<tr><td>INTEGRATOR_RQF</td><td>rqf</td><td>Backward quadrature state at the initial time .</td></tr>
</table>
*/
/** \defgroup scheme_QCQPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::QCQPSolverInput  (QCQP_SOLVER_NUM_IN = 13) [qcqpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. .</td></tr>
<tr><td>QCQP_SOLVER_G</td><td>g</td><td>The vector g: dense, (n x 1) .</td></tr>
<tr><td>QCQP_SOLVER_P</td><td>p</td><td>The vertical stack of all Pi. Each Pi is sparse (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. .</td></tr>
<tr><td>QCQP_SOLVER_Q</td><td>q</td><td>The vertical stack of all qi: dense, (nq n x 1) .</td></tr>
<tr><td>QCQP_SOLVER_R</td><td>r</td><td>The vertical stack of all scalars ri (nq x 1) .</td></tr>
<tr><td>QCQP_SOLVER_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
<tr><td>QCQP_SOLVER_LBA</td><td>lba</td><td>dense, (nc x 1) </td></tr>
<tr><td>QCQP_SOLVER_UBA</td><td>uba</td><td>dense, (nc x 1) </td></tr>
<tr><td>QCQP_SOLVER_LBX</td><td>lbx</td><td>dense, (n x 1) </td></tr>
<tr><td>QCQP_SOLVER_UBX</td><td>ubx</td><td>dense, (n x 1) </td></tr>
<tr><td>QCQP_SOLVER_X0</td><td>x0</td><td>dense, (n x 1) </td></tr>
<tr><td>QCQP_SOLVER_LAM_X0</td><td>lam_x0</td><td>dense </td></tr>
</table>
*/
/** \defgroup scheme_HessLagOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::HessLagOutput  (HESSLAG_NUM_OUT = 6) [hessLagOut]</caption>
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
<caption>Input scheme: CasADi::LinsolInput  (LINSOL_NUM_IN = 4) [linsolIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_A</td><td>A</td><td>The square matrix A: sparse, (n x n). .</td></tr>
<tr><td>LINSOL_B</td><td>B</td><td>The right-hand-side matrix b: dense, (m x n) .</td></tr>
</table>
*/
/** \defgroup scheme_SOCPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SOCPOutput  (SOCP_SOLVER_NUM_OUT = 5) [socpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_X</td><td>x</td><td>The primal solution (n x 1) .</td></tr>
<tr><td>SOCP_SOLVER_COST</td><td>cost</td><td>The primal optimal cost (1 x 1) .</td></tr>
<tr><td>SOCP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to the linear constraints (nc x 1) .</td></tr>
<tr><td>SOCP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds (n x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_NLPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::NLPSolverOutput  (NLP_SOLVER_NUM_OUT = 7) [nlpSolverOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLP_SOLVER_X</td><td>x</td><td>Decision variables at the optimal solution (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_F</td><td>f</td><td>Cost function value at the optimal solution (1 x 1) .</td></tr>
<tr><td>NLP_SOLVER_G</td><td>g</td><td>Constraints function at the optimal solution (ng x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_X</td><td>lam_x</td><td>Lagrange multipliers for bounds on X at the solution (nx x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_G</td><td>lam_g</td><td>Lagrange multipliers for bounds on G at the solution (ng x 1) .</td></tr>
<tr><td>NLP_SOLVER_LAM_P</td><td>lam_p</td><td>Lagrange multipliers for bounds on P at the solution (np x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_SDPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SDPInput  (SDP_SOLVER_NUM_IN = 9) [sdpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDP_SOLVER_F</td><td>f</td><td>The vertical stack of all matrices F_i: ( nm x m) .</td></tr>
<tr><td>SDP_SOLVER_C</td><td>c</td><td>The vector c: ( n x 1) .</td></tr>
<tr><td>SDP_SOLVER_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDP_SOLVER_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
<tr><td>SDP_SOLVER_LBA</td><td>lba</td><td>Lower bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDP_SOLVER_UBA</td><td>uba</td><td>Upper bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDP_SOLVER_LBX</td><td>lbx</td><td>Lower bounds on x ( n x 1 ) .</td></tr>
<tr><td>SDP_SOLVER_UBX</td><td>ubx</td><td>Upper bounds on x ( n x 1 ) .</td></tr>
</table>
*/
/** \defgroup scheme_LPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::LPSolverInput  (LP_SOLVER_NUM_IN = 7) [lpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LP_SOLVER_C</td><td>c</td><td>The vector c: dense (n x 1) .</td></tr>
<tr><td>LP_SOLVER_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
<tr><td>LP_SOLVER_LBA</td><td>lba</td><td>dense, (nc x 1) </td></tr>
<tr><td>LP_SOLVER_UBA</td><td>uba</td><td>dense, (nc x 1) </td></tr>
<tr><td>LP_SOLVER_LBX</td><td>lbx</td><td>dense, (n x 1) </td></tr>
<tr><td>LP_SOLVER_UBX</td><td>ubx</td><td>dense, (n x 1) </td></tr>
</table>
*/
/** \defgroup scheme_RDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::RDAEInput  (RDAE_NUM_IN = 8) [rdaeIn]</caption>
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
<caption>Output scheme: CasADi::NLPOutput  (NL_NUM_OUT = 3) [nlpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>NL_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_SOCPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SOCPInput  (SOCP_SOLVER_NUM_IN = 11) [socpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_G</td><td>g</td><td>The vertical stack of all matrices Gi: ( N x n) .</td></tr>
<tr><td>SOCP_SOLVER_H</td><td>h</td><td>The vertical stack of all vectors hi: ( N x 1) .</td></tr>
<tr><td>SOCP_SOLVER_E</td><td>e</td><td>The vertical stack of all vectors ei: ( nm x 1) .</td></tr>
<tr><td>SOCP_SOLVER_F</td><td>f</td><td>The vertical stack of all scalars fi: ( m x 1) .</td></tr>
<tr><td>SOCP_SOLVER_C</td><td>c</td><td>The vector c: ( n x 1) .</td></tr>
<tr><td>SOCP_SOLVER_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
<tr><td>SOCP_SOLVER_LBA</td><td>lba</td><td>Lower bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SOCP_SOLVER_UBA</td><td>uba</td><td>Upper bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SOCP_SOLVER_LBX</td><td>lbx</td><td>Lower bounds on x ( n x 1 ) .</td></tr>
<tr><td>SOCP_SOLVER_UBX</td><td>ubx</td><td>Upper bounds on x ( n x 1 ) .</td></tr>
</table>
*/
/** \defgroup scheme_SDPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SDPOutput  (SDP_SOLVER_NUM_OUT = 8) [sdpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDP_SOLVER_X</td><td>x</td><td>The primal solution (n x 1) - may be used as initial guess .</td></tr>
<tr><td>SDP_SOLVER_P</td><td>p</td><td>The solution P (m x m) - may be used as initial guess .</td></tr>
<tr><td>SDP_SOLVER_DUAL</td><td>dual</td><td>The dual solution (m x m) - may be used as initial guess .</td></tr>
<tr><td>SDP_SOLVER_COST</td><td>cost</td><td>The primal optimal cost (1 x 1) .</td></tr>
<tr><td>SDP_SOLVER_DUAL_COST</td><td>dual_cost</td><td>The dual optimal cost (1 x 1) .</td></tr>
<tr><td>SDP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to the linear constraints (nc x 1) .</td></tr>
<tr><td>SDP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds (n x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_NLPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::NLPSolverInput  (NLP_SOLVER_NUM_IN = 9) [nlpSolverIn]</caption>
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
/** \defgroup scheme_DAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 5) [daeIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_X</td><td>x</td><td>Differential state .</td></tr>
<tr><td>DAE_Z</td><td>z</td><td>Algebraic state .</td></tr>
<tr><td>DAE_P</td><td>p</td><td>Parameter .</td></tr>
<tr><td>DAE_T</td><td>t</td><td>Explicit time dependence .</td></tr>
</table>
*/
/** \defgroup scheme_DAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 4) [daeOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DAE_ODE</td><td>ode</td><td>Right hand side of the implicit ODE .</td></tr>
<tr><td>DAE_ALG</td><td>alg</td><td>Right hand side of algebraic equations .</td></tr>
<tr><td>DAE_QUAD</td><td>quad</td><td>Right hand side of quadratures equations .</td></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::InputOutputScheme  ( = 42) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::InputOutputScheme  ( = 42) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
</table>
*/
/** \defgroup scheme_StabilizedQPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::StabilizedQPSolverInput  (STABILIZED_QP_SOLVER_NUM_IN = 13) [stabilizedQpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>STABILIZED_QP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. .</td></tr>
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
/** \defgroup scheme_GradFInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::GradFInput  (GRADF_NUM_IN = 3) [gradFIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>GRADF_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_LPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::LPSolverOutput  (LP_SOLVER_NUM_OUT = 5) [lpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>LP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>LP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>LP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
</table>
*/
/** \defgroup scheme_GradFOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::GradFOutput  (GRADF_NUM_OUT = 4) [gradFOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>GRADF_GRAD</td><td>grad</td><td>Jacobian of the constraints .</td></tr>
<tr><td>GRADF_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>GRADF_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_SDQPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::SDQPInput  (SDQP_SOLVER_NUM_IN = 10) [sdqpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDQP_SOLVER_H</td><td>h</td><td>The matrix H: sparse ( n x n) .</td></tr>
<tr><td>SDQP_SOLVER_C</td><td>c</td><td>The vector c: ( n x 1) .</td></tr>
<tr><td>SDQP_SOLVER_F</td><td>f</td><td>The vertical stack of all matrices F_i: ( nm x m) .</td></tr>
<tr><td>SDQP_SOLVER_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDQP_SOLVER_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
<tr><td>SDQP_SOLVER_LBA</td><td>lba</td><td>Lower bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDQP_SOLVER_UBA</td><td>uba</td><td>Upper bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDQP_SOLVER_LBX</td><td>lbx</td><td>Lower bounds on x ( n x 1 ) .</td></tr>
<tr><td>SDQP_SOLVER_UBX</td><td>ubx</td><td>Upper bounds on x ( n x 1 ) .</td></tr>
</table>
*/
/** \defgroup scheme_OCPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::OCPInput  (OCP_NUM_IN = 14) [ocpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>OCP_LBX</td><td>lbx</td><td>States lower bounds (nx x (ns+1)) .</td></tr>
<tr><td>OCP_UBX</td><td>ubx</td><td>States upper bounds (nx x (ns+1)) .</td></tr>
<tr><td>OCP_X_INIT</td><td>x_init</td><td>States initial guess (nx x (ns+1)) .</td></tr>
<tr><td>OCP_LBU</td><td>lbu</td><td>Controls lower bounds (nu x ns) .</td></tr>
<tr><td>OCP_UBU</td><td>ubu</td><td>Controls upper bounds (nu x ns) .</td></tr>
<tr><td>OCP_U_INIT</td><td>u_init</td><td>Controls initial guess (nu x ns) .</td></tr>
<tr><td>OCP_LBP</td><td>lbp</td><td>Parameters lower bounds (np x 1) .</td></tr>
<tr><td>OCP_UBP</td><td>ubp</td><td>Parameters upper bounds (np x 1) .</td></tr>
<tr><td>OCP_P_INIT</td><td>p_init</td><td>Parameters initial guess (np x 1) .</td></tr>
<tr><td>OCP_LBH</td><td>lbh</td><td>Point constraint lower bound (nh x (ns+1)) .</td></tr>
<tr><td>OCP_UBH</td><td>ubh</td><td>Point constraint upper bound (nh x (ns+1)) .</td></tr>
<tr><td>OCP_LBG</td><td>lbg</td><td>Lower bound for the coupling constraints .</td></tr>
<tr><td>OCP_UBG</td><td>ubg</td><td>Upper bound for the coupling constraints .</td></tr>
</table>
*/
/** \defgroup scheme_QCQPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::QCQPSolverOutput  (QCQP_SOLVER_NUM_OUT = 5) [qcqpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>QCQP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>QCQP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>QCQP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
</table>
*/
/** \defgroup scheme_MayerInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::MayerInput  (MAYER_NUM_IN = 3) [mayerIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>MAYER_X</td><td>x</td><td>States at the end of integration (nx x 1) .</td></tr>
<tr><td>MAYER_P</td><td>p</td><td>Problem parameters (np x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_ControlledDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ControlledDAEInput  (CONTROL_DAE_NUM_IN = 10) [controldaeIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONTROL_DAE_T</td><td>t</td><td>Global physical time. (1-by-1) .</td></tr>
<tr><td>CONTROL_DAE_X</td><td>x</td><td>State vector (dimension nx-by-1). Should have same amount of non-zeros as DAEOutput:DAE_RES .</td></tr>
<tr><td>CONTROL_DAE_Z</td><td>z</td><td>Algebraic state vector (dimension np-by-1). .</td></tr>
<tr><td>CONTROL_DAE_P</td><td>p</td><td>Parameter vector (dimension np-by-1). .</td></tr>
<tr><td>CONTROL_DAE_U</td><td>u</td><td>Control vector (dimension nu-by-1). .</td></tr>
<tr><td>CONTROL_DAE_U_INTERP</td><td>u_interp</td><td>Control vector, linearly interpolated (dimension nu-by-1). .</td></tr>
<tr><td>CONTROL_DAE_X_MAJOR</td><td>x_major</td><td>State vector (dimension nx-by-1) at the last major time-step .</td></tr>
<tr><td>CONTROL_DAE_T0</td><td>t0</td><td>Time at start of control interval (1-by-1) .</td></tr>
<tr><td>CONTROL_DAE_TF</td><td>tf</td><td>Time at end of control interval (1-by-1) .</td></tr>
</table>
*/
/** \defgroup scheme_NLPInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::NLPInput  (NL_NUM_IN = 3) [nlpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>NL_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::IntegratorInput  (INTEGRATOR_NUM_IN = 5) [integratorIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>x0</td><td>Differential state at the initial time .</td></tr>
<tr><td>INTEGRATOR_P</td><td>p</td><td>Parameters .</td></tr>
<tr><td>INTEGRATOR_RX0</td><td>rx0</td><td>Backward differential state at the final time .</td></tr>
<tr><td>INTEGRATOR_RP</td><td>rp</td><td>Backward parameter vector .</td></tr>
</table>
*/
/** \defgroup scheme_QPSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::QPSolverInput  (QP_SOLVER_NUM_IN = 10) [qpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. .</td></tr>
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
/** \defgroup scheme_OCPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 5) [ocpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>x_opt</td><td>Optimal state trajectory .</td></tr>
<tr><td>OCP_U_OPT</td><td>u_opt</td><td>Optimal control trajectory .</td></tr>
<tr><td>OCP_P_OPT</td><td>p_opt</td><td>Optimal parameters .</td></tr>
<tr><td>OCP_COST</td><td>cost</td><td>Objective/cost function for optimal solution (1 x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_RDAEOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::RDAEOutput  (RDAE_NUM_OUT = 4) [rdaeOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>RDAE_ODE</td><td>ode</td><td>Right hand side of ODE. .</td></tr>
<tr><td>RDAE_ALG</td><td>alg</td><td>Right hand side of algebraic equations. .</td></tr>
<tr><td>RDAE_QUAD</td><td>quad</td><td>Right hand side of quadratures. .</td></tr>
</table>
*/
/** \defgroup scheme_JacGOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::JacGOutput  (JACG_NUM_OUT = 4) [jacGOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_JAC</td><td>jac</td><td>Jacobian of the constraints .</td></tr>
<tr><td>JACG_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>JACG_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_LinsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::LinsolOutput  (LINSOL_NUM_OUT = 2) [linsolOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_X</td><td>X</td><td>Solution to the linear system of equations .</td></tr>
</table>
*/
/** \defgroup scheme_SDQPOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::SDQPOutput  (SDQP_SOLVER_NUM_OUT = 8) [sdqpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDQP_SOLVER_X</td><td>x</td><td>The primal solution (n x 1) - may be used as initial guess .</td></tr>
<tr><td>SDQP_SOLVER_P</td><td>p</td><td>The solution P (m x m) - may be used as initial guess .</td></tr>
<tr><td>SDQP_SOLVER_DUAL</td><td>dual</td><td>The dual solution (m x m) - may be used as initial guess .</td></tr>
<tr><td>SDQP_SOLVER_COST</td><td>cost</td><td>The primal optimal cost (1 x 1) .</td></tr>
<tr><td>SDQP_SOLVER_DUAL_COST</td><td>dual_cost</td><td>The dual optimal cost (1 x 1) .</td></tr>
<tr><td>SDQP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to the linear constraints (nc x 1) .</td></tr>
<tr><td>SDQP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds (n x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_ControlSimulatorInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::ControlSimulatorInput  (CONTROLSIMULATOR_NUM_IN = 4) [controlsimulatorIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONTROLSIMULATOR_X0</td><td>x0</td><td>Differential or algebraic state at t0 (dimension nx-by-1) .</td></tr>
<tr><td>CONTROLSIMULATOR_P</td><td>p</td><td>Parameters that are fixed over the entire horizon (dimension np-by-1) .</td></tr>
<tr><td>CONTROLSIMULATOR_U</td><td>u</td><td>Parameters that change over the integration intervals (dimension (ns-1)-by-nu) .</td></tr>
</table>
*/
/** \defgroup scheme_JacGInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::JacGInput  (JACG_NUM_IN = 3) [jacGIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>JACG_P</td><td>p</td><td>Fixed parameter .</td></tr>
</table>
*/
/** \defgroup scheme_HessLagInput
<a name='schemes'></a><table>
<caption>Input scheme: CasADi::HessLagInput  (HESSLAG_NUM_IN = 5) [hessLagIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HESSLAG_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>HESSLAG_P</td><td>p</td><td>Fixed parameter .</td></tr>
<tr><td>HESSLAG_LAM_F</td><td>lam_f</td><td>Multiplier for f. Just a scalar factor for the objective that the NLP solver might use to scale the objective. .</td></tr>
<tr><td>HESSLAG_LAM_G</td><td>lam_g</td><td>Multiplier for g .</td></tr>
</table>
*/
/** \defgroup scheme_QPSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: CasADi::QPSolverOutput  (QP_SOLVER_NUM_OUT = 5) [qpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>QP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>QP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>QP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
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
/** \class CasADi::SQICInternal
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::SQICSolver
\n
\par
@copydoc scheme_QPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::SnoptInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::SnoptSolver
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
/** \class CasADi::QPStabilizerInternal
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::QPStabilizer
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
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
/** \class CasADi::StabilizedSQPInternal
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
*/
/** \class CasADi::StabilizedSQPMethod
\n
\par
@copydoc scheme_NLPSolverInput
<br/>
@copydoc scheme_NLPSolverOutput
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
/** \class CasADi::GenericIntegratorInternal
\n
\par
@copydoc scheme_IntegratorInput
<br/>
@copydoc scheme_IntegratorOutput
*/
/** \class CasADi::GenericIntegrator
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
/** \class CasADi::StabilizedSQICInternal
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::StabilizedSQICSolver
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
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
/** \class CasADi::StabilizedQPSolverInternal
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
*/
/** \class CasADi::StabilizedQPSolver
\n
\par
@copydoc scheme_StabilizedQPSolverInput
<br/>
@copydoc scheme_QPSolverOutput
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
