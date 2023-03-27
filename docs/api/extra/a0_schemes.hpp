/** \defgroup scheme_ConicInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::ConicInput  (CONIC_NUM_IN = 12)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONIC_H</td><td>h</td><td>The square matrix H: sparse, (n x n). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>CONIC_G</td><td>g</td><td>The vector g: dense, (n x 1)</td></tr>
<tr><td>CONIC_A</td><td>a</td><td>The matrix A: sparse, (nc x n) - product with x must be dense.</td></tr>
<tr><td>CONIC_LBA</td><td>lba</td><td>dense, (nc x 1)</td></tr>
<tr><td>CONIC_UBA</td><td>uba</td><td>dense, (nc x 1)</td></tr>
<tr><td>CONIC_LBX</td><td>lbx</td><td>dense, (n x 1)</td></tr>
<tr><td>CONIC_UBX</td><td>ubx</td><td>dense, (n x 1)</td></tr>
<tr><td>CONIC_X0</td><td>x0</td><td>dense, (n x 1)</td></tr>
<tr><td>CONIC_LAM_X0</td><td>lam_x0</td><td>dense</td></tr>
<tr><td>CONIC_LAM_A0</td><td>lam_a0</td><td>dense</td></tr>
<tr><td>CONIC_Q</td><td>q</td><td>The matrix Q: sparse symmetric, (np^2 x n)</td></tr>
<tr><td>CONIC_P</td><td>p</td><td>The matrix P: sparse symmetric, (np x np)</td></tr>
</table>
*/
/** \defgroup scheme_ConicOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::ConicOutput  (CONIC_NUM_OUT = 4)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>CONIC_X</td><td>x</td><td>The primal solution.</td></tr>
<tr><td>CONIC_COST</td><td>cost</td><td>The optimal cost.</td></tr>
<tr><td>CONIC_LAM_A</td><td>lam_a</td><td>The dual solution corresponding to linear bounds.</td></tr>
<tr><td>CONIC_LAM_X</td><td>lam_x</td><td>The dual solution corresponding to simple bounds.</td></tr>
</table>
*/
/** \defgroup scheme_DpleInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::DpleInput  (DPLE_NUM_IN = 2)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DPLE_A</td><td>a</td><td>A matrices (horzcat when const_dim, diagcat otherwise) [a].</td></tr>
<tr><td>DPLE_V</td><td>v</td><td>V matrices (horzcat when const_dim, diagcat otherwise) [v].</td></tr>
</table>
*/
/** \defgroup scheme_DpleOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::DpleOutput  (DPLE_NUM_OUT = 1)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>DPLE_P</td><td>p</td><td>Lyapunov matrix (horzcat when const_dim, diagcat otherwise) (Cholesky of P if pos_def) [p].</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::IntegratorInput  (INTEGRATOR_NUM_IN = 7)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>x0</td><td>Differential state at the initial time.</td></tr>
<tr><td>INTEGRATOR_Z0</td><td>z0</td><td>Initial guess for the algebraic variable at the initial time.</td></tr>
<tr><td>INTEGRATOR_P</td><td>p</td><td>Parameters.</td></tr>
<tr><td>INTEGRATOR_U</td><td>u</td><td>Piecewise constant control, a new control interval starts at each output time.</td></tr>
<tr><td>INTEGRATOR_ADJ_XF</td><td>adj_xf</td><td>Adjoint seeds corresponding to the states at the output times.</td></tr>
<tr><td>INTEGRATOR_ADJ_ZF</td><td>adj_zf</td><td>Adjoint seeds corresponding to the algebraic variables at the output times.</td></tr>
<tr><td>INTEGRATOR_ADJ_QF</td><td>adj_qf</td><td>Adjoint seeds corresponding to the quadratures at the output times.</td></tr>
</table>
*/
/** \defgroup scheme_IntegratorOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::IntegratorOutput  (INTEGRATOR_NUM_OUT = 7)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>INTEGRATOR_XF</td><td>xf</td><td>Differential state at all output times.</td></tr>
<tr><td>INTEGRATOR_ZF</td><td>zf</td><td>Algebraic variable at all output times.</td></tr>
<tr><td>INTEGRATOR_QF</td><td>qf</td><td>Quadrature state at all output times.</td></tr>
<tr><td>INTEGRATOR_ADJ_X0</td><td>adj_x0</td><td>Adjoint sensitivities corresponding to the initial state.</td></tr>
<tr><td>INTEGRATOR_ADJ_Z0</td><td>adj_z0</td><td>Adjoint sensitivities corresponding to the algebraic variable guess.</td></tr>
<tr><td>INTEGRATOR_ADJ_P</td><td>adj_p</td><td>Adjoint sensitivities corresponding to the parameter vector.</td></tr>
<tr><td>INTEGRATOR_ADJ_U</td><td>adj_u</td><td>Adjoint sensitivities corresponding to the control vector.</td></tr>
</table>
*/
/** \defgroup scheme_NLPInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NLPInput  (NL_NUM_IN = 2)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_X</td><td>x</td><td>Decision variable.</td></tr>
<tr><td>NL_P</td><td>p</td><td>Fixed parameter.</td></tr>
</table>
*/
/** \defgroup scheme_NLPOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NLPOutput  (NL_NUM_OUT = 2)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NL_F</td><td>f</td><td>Objective function.</td></tr>
<tr><td>NL_G</td><td>g</td><td>Constraint function.</td></tr>
</table>
*/
/** \defgroup scheme_NlpsolInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::NlpsolInput  (NLPSOL_NUM_IN = 8)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLPSOL_X0</td><td>x0</td><td>Decision variables, initial guess (nx x 1)</td></tr>
<tr><td>NLPSOL_P</td><td>p</td><td>Value of fixed parameters (np x 1)</td></tr>
<tr><td>NLPSOL_LBX</td><td>lbx</td><td>Decision variables lower bound (nx x 1), default -inf.</td></tr>
<tr><td>NLPSOL_UBX</td><td>ubx</td><td>Decision variables upper bound (nx x 1), default +inf.</td></tr>
<tr><td>NLPSOL_LBG</td><td>lbg</td><td>Constraints lower bound (ng x 1), default -inf.</td></tr>
<tr><td>NLPSOL_UBG</td><td>ubg</td><td>Constraints upper bound (ng x 1), default +inf.</td></tr>
<tr><td>NLPSOL_LAM_X0</td><td>lam_x0</td><td>Lagrange multipliers for bounds on X, initial guess (nx x 1)</td></tr>
<tr><td>NLPSOL_LAM_G0</td><td>lam_g0</td><td>Lagrange multipliers for bounds on G, initial guess (ng x 1)</td></tr>
</table>
*/
/** \defgroup scheme_NlpsolOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::NlpsolOutput  (NLPSOL_NUM_OUT = 6)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>NLPSOL_X</td><td>x</td><td>Decision variables at the optimal solution (nx x 1)</td></tr>
<tr><td>NLPSOL_F</td><td>f</td><td>Cost function value at the optimal solution (1 x 1)</td></tr>
<tr><td>NLPSOL_G</td><td>g</td><td>Constraints function at the optimal solution (ng x 1)</td></tr>
<tr><td>NLPSOL_LAM_X</td><td>lam_x</td><td>Lagrange multipliers for bounds on X at the solution (nx x 1)</td></tr>
<tr><td>NLPSOL_LAM_G</td><td>lam_g</td><td>Lagrange multipliers for bounds on G at the solution (ng x 1)</td></tr>
<tr><td>NLPSOL_LAM_P</td><td>lam_p</td><td>Lagrange multipliers for bounds on P at the solution (np x 1)</td></tr>
</table>
*/
/** \defgroup scheme_RootfinderInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::RootfinderInput  (ROOTFINDER_NUM_IN = 2)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>ROOTFINDER_X0</td><td>x0</td><td>Initial guess for the solution.</td></tr>
<tr><td>ROOTFINDER_P</td><td>p</td><td>Parameters.</td></tr>
</table>
*/
/** \defgroup scheme_RootfinderOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::RootfinderOutput  (ROOTFINDER_NUM_OUT = 1)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>ROOTFINDER_X</td><td>x</td><td>Solution to the system of equations.</td></tr>
</table>
*/
/** \defgroup scheme_SimulatorInput
<a name='schemes'></a><table>
<caption>Input scheme: casadi::SimulatorInput  (SIMULATOR_NUM_IN = 4)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SIMULATOR_X0</td><td>x0</td><td>Differential state at the initial time.</td></tr>
<tr><td>SIMULATOR_U</td><td>u</td><td>Controls.</td></tr>
<tr><td>SIMULATOR_Z0</td><td>z0</td><td>Initial guess for the algebraic variable at the initial time.</td></tr>
<tr><td>SIMULATOR_P</td><td>p</td><td>Parameters.</td></tr>
</table>
*/
/** \defgroup scheme_SimulatorOutput
<a name='schemes'></a><table>
<caption>Output scheme: casadi::SimulatorOutput  (SIMULATOR_NUM_OUT = 2)</caption>
<tr><th>Full name</th><th>Short</th><th>Description</th></tr>
<tr><td>SIMULATOR_X</td><td>x</td><td>Differential state.</td></tr>
<tr><td>SIMULATOR_Z</td><td>z</td><td>Algebraic variable.</td></tr>
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
/** \class casadi::CvodesSimulator
\n
\par
@copydoc scheme_SimulatorInput
<br/>
@copydoc scheme_SimulatorOutput
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
/// \cond INTERNAL
/** \class casadi::SundialsSimulator
\n
\par
@copydoc scheme_SimulatorInput
<br/>
@copydoc scheme_SimulatorOutput
*/
/// \endcond
