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
/** \defgroup scheme_SDPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::SDPOutput  (SDP_SOLVER_NUM_OUT = 7) [sdpOut]</caption>
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
/** \defgroup scheme_LPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::LPStruct  ( = 1) [lpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LP_STRUCT_A</td><td>a</td><td>The matrix A: sparse .</td></tr>
</table>
*/
/** \defgroup scheme_SOCPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::SOCPStruct  ( = 2) [socpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SOCP_STRUCT_G</td><td>g</td><td>The horizontal stack of all matrices Gi: ( n x N) .</td></tr>
<tr><td>SOCP_STRUCT_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
</table>
*/
/** \defgroup scheme_DLEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::DLEOutput  (DLE_NUM_OUT = 1) [dleOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DLE_P</td><td>p</td><td>P matrix .</td></tr>
</table>
*/
/** \defgroup scheme_LR_DPLEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::LR_DPLEInput  (LR_DPLE_NUM_IN = 4) [lrdpleIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_DPLE_A</td><td>a</td><td>A matrices (horzcat when const_dim, blkdiag otherwise) .</td></tr>
<tr><td>LR_DPLE_V</td><td>v</td><td>V matrices (horzcat when const_dim, blkdiag otherwise) .</td></tr>
<tr><td>LR_DPLE_C</td><td>c</td><td>C matrix .</td></tr>
<tr><td>LR_DPLE_H</td><td>h</td><td>H matrix: horizontal stack of all Hi .</td></tr>
</table>
*/
/** \defgroup scheme_SDPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::SDPInput  (SDP_SOLVER_NUM_IN = 8) [sdpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDP_SOLVER_F</td><td>f</td><td>The horizontal stack of all matrices F_i: ( m x nm) .</td></tr>
<tr><td>SDP_SOLVER_C</td><td>c</td><td>The vector c: ( n x 1) .</td></tr>
<tr><td>SDP_SOLVER_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDP_SOLVER_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
<tr><td>SDP_SOLVER_LBA</td><td>lba</td><td>Lower bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDP_SOLVER_UBA</td><td>uba</td><td>Upper bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDP_SOLVER_LBX</td><td>lbx</td><td>Lower bounds on x ( n x 1 ) .</td></tr>
<tr><td>SDP_SOLVER_UBX</td><td>ubx</td><td>Upper bounds on x ( n x 1 ) .</td></tr>
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
/** \defgroup scheme_DLEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::DLEInput  (DLE_NUM_IN = 2) [dleIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DLE_A</td><td>a</td><td>A matrix .</td></tr>
<tr><td>DLE_V</td><td>v</td><td>V matrix .</td></tr>
</table>
*/
/** \defgroup scheme_QcqpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::QcqpSolverOutput  (QCQP_SOLVER_NUM_OUT = 4) [qcqpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>QCQP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>QCQP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>QCQP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
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
/** \defgroup scheme_DPLEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::DPLEInput  (DPLE_NUM_IN = 2) [dpleIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DPLE_A</td><td>a</td><td>A matrices (horzcat when const_dim, blkdiag otherwise) .</td></tr>
<tr><td>DPLE_V</td><td>v</td><td>V matrices (horzcat when const_dim, blkdiag otherwise) .</td></tr>
</table>
*/
/** \defgroup scheme_CLEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::CLEOutput  (CLE_NUM_OUT = 1) [cleOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CLE_P</td><td>p</td><td>Lyapunov matrix .</td></tr>
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
/** \defgroup scheme_NLPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NLPOutput  (NL_NUM_OUT = 2) [nlpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>NL_G</td><td>g</td><td>Constraint function .</td></tr>
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
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Input scheme: casadi::InputOutputScheme  ( = 54) []</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
</table>
*/
/** \defgroup scheme_InputOutputScheme
<a name='schemes'></a><table>
<caption>Output scheme: casadi::InputOutputScheme  ( = 54) []</caption>
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
/** \defgroup scheme_DpleVecStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::DpleVecStruct  ( = 2) [dpleStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>Dple_STRUCT_A</td><td>a</td><td>Sparsities for A_i, blkdiag form .</td></tr>
<tr><td>Dple_STRUCT_V</td><td>v</td><td>Sparsities for V_i, blkdiag form .</td></tr>
</table>
*/
/** \defgroup scheme_QPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::QPStruct  ( = 2) [qpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QP_STRUCT_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>QP_STRUCT_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
</table>
*/
/** \defgroup scheme_SDQPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::SDQPOutput  (SDQP_SOLVER_NUM_OUT = 7) [sdqpOut]</caption>
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
/** \defgroup scheme_SOCPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::SOCPOutput  (SOCP_SOLVER_NUM_OUT = 4) [socpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_X</td><td>x</td><td>The primal solution (n x 1) .</td></tr>
<tr><td>SOCP_SOLVER_COST</td><td>cost</td><td>The primal optimal cost (1 x 1) .</td></tr>
<tr><td>SOCP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to the linear constraints (nc x 1) .</td></tr>
<tr><td>SOCP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds (n x 1) .</td></tr>
</table>
*/
/** \defgroup scheme_SDQPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::SDQPInput  (SDQP_SOLVER_NUM_IN = 9) [sdqpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDQP_SOLVER_H</td><td>h</td><td>The matrix H: sparse ( n x n) .</td></tr>
<tr><td>SDQP_SOLVER_C</td><td>c</td><td>The vector c: ( n x 1) .</td></tr>
<tr><td>SDQP_SOLVER_F</td><td>f</td><td>The horizontal stack of all matrices F_i: ( m x nm) .</td></tr>
<tr><td>SDQP_SOLVER_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDQP_SOLVER_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
<tr><td>SDQP_SOLVER_LBA</td><td>lba</td><td>Lower bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDQP_SOLVER_UBA</td><td>uba</td><td>Upper bounds on Ax ( nc x 1) .</td></tr>
<tr><td>SDQP_SOLVER_LBX</td><td>lbx</td><td>Lower bounds on x ( n x 1 ) .</td></tr>
<tr><td>SDQP_SOLVER_UBX</td><td>ubx</td><td>Upper bounds on x ( n x 1 ) .</td></tr>
</table>
*/
/** \defgroup scheme_ControlledDAEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::ControlledDAEInput  (CONTROL_DAE_NUM_IN = 9) [controldaeIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONTROL_DAE_T</td><td>t</td><td>Global physical time. (1-by-1) .</td></tr>
<tr><td>CONTROL_DAE_X</td><td>x</td><td>State vector (dimension nx-by-1). Should have the same amount of non-zeros as DAEOutput:DAE_RES </td></tr>
<tr><td>CONTROL_DAE_Z</td><td>z</td><td>Algebraic state vector (dimension np-by-1). .</td></tr>
<tr><td>CONTROL_DAE_P</td><td>p</td><td>Parameter vector (dimension np-by-1). .</td></tr>
<tr><td>CONTROL_DAE_U</td><td>u</td><td>Control vector (dimension nu-by-1). .</td></tr>
<tr><td>CONTROL_DAE_U_INTERP</td><td>u_interp</td><td>Control vector, linearly interpolated (dimension nu-by-1). .</td></tr>
<tr><td>CONTROL_DAE_X_MAJOR</td><td>x_major</td><td>State vector (dimension nx-by-1) at the last major time-step .</td></tr>
<tr><td>CONTROL_DAE_T0</td><td>t0</td><td>Time at start of control interval (1-by-1) .</td></tr>
<tr><td>CONTROL_DAE_TF</td><td>tf</td><td>Time at end of control interval (1-by-1) .</td></tr>
</table>
*/
/** \defgroup scheme_DleStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::DleStruct  ( = 2) [dleStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>Dle_STRUCT_A</td><td>a</td><td>The matrix A .</td></tr>
<tr><td>Dle_STRUCT_V</td><td>v</td><td>The matrix V .</td></tr>
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
/** \defgroup scheme_SDQPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::SDQPStruct  ( = 4) [sdqpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDQP_STRUCT_H</td><td>h</td><td>The matrix H: sparse ( n x n) .</td></tr>
<tr><td>SDQP_STRUCT_F</td><td>f</td><td>The horizontal stack of all matrices F_i: ( m x nm) .</td></tr>
<tr><td>SDQP_STRUCT_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDQP_STRUCT_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
</table>
*/
/** \defgroup scheme_CLEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::CLEInput  (CLE_NUM_IN = 2) [cleIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CLE_A</td><td>a</td><td>A matrix .</td></tr>
<tr><td>CLE_V</td><td>v</td><td>V matrix .</td></tr>
</table>
*/
/** \defgroup scheme_LR_DPLEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LR_DPLEOutput  (LR_DPLE_NUM_OUT = 1) [lrdpleOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_DPLE_Y</td><td>y</td><td>Lyapunov matrix (horzcat when const_dim, blkdiag otherwise) (Cholesky of P if pos_def) .</td></tr>
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
/** \defgroup scheme_JacGOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::JacGOutput  (JACG_NUM_OUT = 3) [jacGOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>JACG_JAC</td><td>jac</td><td>Jacobian of the constraints .</td></tr>
<tr><td>JACG_F</td><td>f</td><td>Objective function .</td></tr>
<tr><td>JACG_G</td><td>g</td><td>Constraint function .</td></tr>
</table>
*/
/** \defgroup scheme_QCQPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::QCQPStruct  ( = 3) [qcqpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QCQP_STRUCT_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>QCQP_STRUCT_P</td><td>p</td><td>The horizontal stack of all Pi. Each Pi is sparse (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>QCQP_STRUCT_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
</table>
*/
/** \defgroup scheme_LinsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LinsolOutput  (LINSOL_NUM_OUT = 1) [linsolOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LINSOL_X</td><td>X</td><td>Solution to the linear system of equations .</td></tr>
</table>
*/
/** \defgroup scheme_QcqpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::QcqpSolverInput  (QCQP_SOLVER_NUM_IN = 12) [qcqpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>QCQP_SOLVER_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
<tr><td>QCQP_SOLVER_G</td><td>g</td><td>The vector g: dense, (n x 1) .</td></tr>
<tr><td>QCQP_SOLVER_P</td><td>p</td><td>The horizontal stack of all Pi. Each Pi is sparse (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. </td></tr>
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
/** \defgroup scheme_CleStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::CleStruct  ( = 3) [cleStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>Cle_STRUCT_A</td><td>a</td><td>The matrix A .</td></tr>
<tr><td>Cle_STRUCT_V</td><td>v</td><td>The matrix V .</td></tr>
<tr><td>Cle_STRUCT_C</td><td>c</td><td>The matrix C (defaults to unity) .</td></tr>
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
/** \defgroup scheme_ControlSimulatorInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::ControlSimulatorInput  (CONTROLSIMULATOR_NUM_IN = 3) [controlsimulatorIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONTROLSIMULATOR_X0</td><td>x0</td><td>Differential or algebraic state at t0 (dimension nx-by-1) .</td></tr>
<tr><td>CONTROLSIMULATOR_P</td><td>p</td><td>Parameters that are fixed over the entire horizon (dimension np-by-1) .</td></tr>
<tr><td>CONTROLSIMULATOR_U</td><td>u</td><td>Parameters that change over the integration intervals (dimension nu-by-(ns-1)) .</td></tr>
</table>
*/
/** \defgroup scheme_SOCPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::SOCPInput  (SOCP_SOLVER_NUM_IN = 10) [socpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SOCP_SOLVER_G</td><td>g</td><td>The horizontal stack of all matrices Gi: ( n x N) .</td></tr>
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
/** \defgroup scheme_LpSolverOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LpSolverOutput  (LP_SOLVER_NUM_OUT = 4) [lpOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LP_SOLVER_X</td><td>x</td><td>The primal solution .</td></tr>
<tr><td>LP_SOLVER_COST</td><td>cost</td><td>The optimal cost .</td></tr>
<tr><td>LP_SOLVER_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds .</td></tr>
<tr><td>LP_SOLVER_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds .</td></tr>
</table>
*/
/** \defgroup scheme_SDPStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::SDPStruct  ( = 3) [sdpStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SDP_STRUCT_F</td><td>f</td><td>The horizontal stack of all matrices F_i: ( m x nm) .</td></tr>
<tr><td>SDP_STRUCT_G</td><td>g</td><td>The matrix G: ( m x m) .</td></tr>
<tr><td>SDP_STRUCT_A</td><td>a</td><td>The matrix A: ( nc x n) .</td></tr>
</table>
*/
/** \defgroup scheme_LrDpleVecStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::LrDpleVecStruct  ( = 4) [lrdpleStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_Dple_STRUCT_A</td><td>a</td><td>Sparsities for A_i, blkdiag form .</td></tr>
<tr><td>LR_Dple_STRUCT_V</td><td>v</td><td>Sparsities for V_i, blkdiag form .</td></tr>
<tr><td>LR_Dple_STRUCT_C</td><td>c</td><td>Sparsities for C_i (defaults to unity), blkdiag form .</td></tr>
<tr><td>LR_Dple_STRUCT_H</td><td>h</td><td>Sparsities for H_i (defaults to unity), blkdiag form .</td></tr>
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
/** \defgroup scheme_LR_DLEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::LR_DLEOutput  (LR_DLE_NUM_OUT = 1) [lrdleOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_DLE_Y</td><td>y</td><td>Y matrix, blkdiag form .</td></tr>
</table>
*/
/** \defgroup scheme_LpSolverInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::LpSolverInput  (LP_SOLVER_NUM_IN = 6) [lpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LP_SOLVER_C</td><td>c</td><td>The vector c: dense (n x 1) .</td></tr>
<tr><td>LP_SOLVER_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense. .</td></tr>
<tr><td>LP_SOLVER_LBA</td><td>lba</td><td>dense, (nc x 1) </td></tr>
<tr><td>LP_SOLVER_UBA</td><td>uba</td><td>dense, (nc x 1) </td></tr>
<tr><td>LP_SOLVER_LBX</td><td>lbx</td><td>dense, (n x 1) </td></tr>
<tr><td>LP_SOLVER_UBX</td><td>ubx</td><td>dense, (n x 1) </td></tr>
</table>
*/
/** \defgroup scheme_LrDleStruct
<a name='schemes'></a><table>
<caption>Struct scheme: casadi::LrDleStruct  ( = 4) [lrdleStruct]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_DLE_STRUCT_A</td><td>a</td><td>The matrix A .</td></tr>
<tr><td>LR_DLE_STRUCT_V</td><td>v</td><td>The matrix V .</td></tr>
<tr><td>LR_DLE_STRUCT_C</td><td>c</td><td>The matrix C (defaults to unity) .</td></tr>
<tr><td>LR_DLE_STRUCT_H</td><td>h</td><td>H matrix: horizontal stack of all Hi .</td></tr>
</table>
*/
/** \defgroup scheme_DPLEOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::DPLEOutput  (DPLE_NUM_OUT = 1) [dpleOut]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DPLE_P</td><td>p</td><td>Lyapunov matrix (horzcat when const_dim, blkdiag otherwise) (Cholesky of P if pos_def) .</td></tr>
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
/** \defgroup scheme_LR_DLEInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::LR_DLEInput  (LR_DLE_NUM_IN = 4) [lrdleIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>LR_DLE_A</td><td>a</td><td>A matrix .</td></tr>
<tr><td>LR_DLE_V</td><td>v</td><td>V matrix .</td></tr>
<tr><td>LR_DLE_C</td><td>c</td><td>C matrix .</td></tr>
<tr><td>LR_DLE_H</td><td>h</td><td>H matrix: horizontal stack of all Hi .</td></tr>
</table>
*/
/** \defgroup scheme_HNLPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::HNLPInput  (HNL_NUM_IN = 3) [hnlpIn]</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>HNL_X</td><td>x</td><td>Decision variable .</td></tr>
<tr><td>HNL_P</td><td>p</td><td>Fixed parameter .</td></tr>
<tr><td>HNL_TAU</td><td>tau</td><td>Homotopy parameter .</td></tr>
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
/** \class casadi::CleInternal
\n
\par
@copydoc scheme_CLEInput
<br/>
@copydoc scheme_CLEOutput
*/
/// \endcond
/** \addtogroup general_CleSolver
\n
\par
@copydoc scheme_CLEInput
<br/>
@copydoc scheme_CLEOutput
*/
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
/** \class casadi::CondensingIndefDpleInternal
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::ControlSimulatorInternal
\n
\par
@copydoc scheme_ControlSimulatorInput
*/
/// \endcond
/** \class casadi::ControlSimulator
\n
\par
@copydoc scheme_ControlSimulatorInput
*/
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
/** \class casadi::DleInternal
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \endcond
/** \addtogroup general_DleSolver
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \cond INTERNAL
/** \class casadi::DleToLrDle
\n
\par
@copydoc scheme_LR_DLEInput
<br/>
@copydoc scheme_LR_DLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::DpleInternal
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \endcond
/** \addtogroup general_DpleSolver
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \cond INTERNAL
/** \class casadi::DpleToDle
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::DpleToLrDple
\n
\par
@copydoc scheme_LR_DPLEInput
<br/>
@copydoc scheme_LR_DPLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::DsdpInterface
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::FixedSmithDleInternal
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::FixedSmithLrDleInternal
\n
\par
@copydoc scheme_LR_DLEInput
<br/>
@copydoc scheme_LR_DLEOutput
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
/** \class casadi::HomotopyNLPInternal
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/** \addtogroup general_HomotopyNlpSolver
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
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
/** \class casadi::LiftingIndefDpleInternal
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LiftingLrDpleInternal
\n
\par
@copydoc scheme_LR_DPLEInput
<br/>
@copydoc scheme_LR_DPLEOutput
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
/** \class casadi::LpSolverInternal
\n
\par
@copydoc scheme_LpSolverInput
<br/>
@copydoc scheme_LpSolverOutput
*/
/// \endcond
/** \addtogroup general_LpSolver
\n
\par
@copydoc scheme_LpSolverInput
<br/>
@copydoc scheme_LpSolverOutput
*/
/// \cond INTERNAL
/** \class casadi::LpToQp
\n
\par
@copydoc scheme_LpSolverInput
<br/>
@copydoc scheme_LpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LrDleInternal
\n
\par
@copydoc scheme_LR_DLEInput
<br/>
@copydoc scheme_LR_DLEOutput
*/
/// \endcond
/** \addtogroup general_LrDleSolver
\n
\par
@copydoc scheme_LR_DLEInput
<br/>
@copydoc scheme_LR_DLEOutput
*/
/// \cond INTERNAL
/** \class casadi::LrDleToDle
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::LrDpleInternal
\n
\par
@copydoc scheme_LR_DPLEInput
<br/>
@copydoc scheme_LR_DPLEOutput
*/
/// \endcond
/** \addtogroup general_LrDpleSolver
\n
\par
@copydoc scheme_LR_DPLEInput
<br/>
@copydoc scheme_LR_DPLEOutput
*/
/// \cond INTERNAL
/** \class casadi::LrDpleToDple
\n
\par
@copydoc scheme_LR_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \endcond
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
/** \class casadi::PsdIndefDpleInternal
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \endcond
/** \class casadi::PsdIndefDpleSolver
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
*/
/// \cond INTERNAL
/** \class casadi::QcqpSolverInternal
\n
\par
@copydoc scheme_QcqpSolverInput
<br/>
@copydoc scheme_QcqpSolverOutput
*/
/// \endcond
/** \addtogroup general_QcqpSolver
\n
\par
@copydoc scheme_QcqpSolverInput
<br/>
@copydoc scheme_QcqpSolverOutput
*/
/// \cond INTERNAL
/** \class casadi::QcqpToSocp
\n
\par
@copydoc scheme_QcqpSolverInput
<br/>
@copydoc scheme_QcqpSolverOutput
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
/** \class casadi::QpToQcqp
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
/** \class casadi::SdpSolverInternal
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/// \endcond
/** \addtogroup general_SdpSolver
\n
\par
@copydoc scheme_SDPInput
<br/>
@copydoc scheme_SDPOutput
*/
/// \cond INTERNAL
/** \class casadi::SdqpSolverInternal
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/// \endcond
/** \addtogroup general_SdqpSolver
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/// \cond INTERNAL
/** \class casadi::SdqpToSdp
\n
\par
@copydoc scheme_SDQPInput
<br/>
@copydoc scheme_SDQPOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimpleHomotopyNlp
\n
\par
@copydoc scheme_NlpSolverInput
<br/>
@copydoc scheme_NlpSolverOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimpleIndefCleInternal
\n
\par
@copydoc scheme_CLEInput
<br/>
@copydoc scheme_CLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimpleIndefDleInternal
\n
\par
@copydoc scheme_DLEInput
<br/>
@copydoc scheme_DLEOutput
*/
/// \endcond
/// \cond INTERNAL
/** \class casadi::SimpleIndefDpleInternal
\n
\par
@copydoc scheme_DPLEInput
<br/>
@copydoc scheme_DPLEOutput
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
/** \class casadi::SmithLrDleInternal
\n
\par
@copydoc scheme_LR_DLEInput
<br/>
@copydoc scheme_LR_DLEOutput
*/
/// \endcond
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
/** \class casadi::SocpSolverInternal
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/// \endcond
/** \addtogroup general_SocpSolver
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
*/
/// \cond INTERNAL
/** \class casadi::SocpToSdp
\n
\par
@copydoc scheme_SOCPInput
<br/>
@copydoc scheme_SOCPOutput
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
