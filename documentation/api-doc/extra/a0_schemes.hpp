/** \class CasADi::NLPSolverInternal
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::NLPSolver
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::Interfaces::OOQPInternal
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::Interfaces::OOQPSolver
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShootingInternal
<table>
<caption>Input scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory.</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory.</td></tr>
<tr><td>OCP_XP_OPT</td><td>Optimal state derivative trajectory.</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters.</td></tr>
</table>
*/
/** \class CasADi::OptimalControl::MultipleShooting
<table>
<caption>Input scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory.</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory.</td></tr>
<tr><td>OCP_XP_OPT</td><td>Optimal state derivative trajectory.</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters.</td></tr>
</table>
*/
/** \class CasADi::SimulatorInternal
<table>
<caption>Input scheme: CasADi::IntegratorInput  (INTEGRATOR_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>Differential or algebraic state at t0 (dimension nx-by-1)</td></tr>
<tr><td>INTEGRATOR_P</td><td>Parameters p (dimension np-by-1)</td></tr>
<tr><td>INTEGRATOR_XP0</td><td>State derivative at t0 (dimension nx-by-1) This input may be changed during an IDASIntegrator::evaluate()</td></tr>
</table>
*/
/** \class CasADi::Simulator
<table>
<caption>Input scheme: CasADi::IntegratorInput  (INTEGRATOR_NUM_IN = 3) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>INTEGRATOR_X0</td><td>Differential or algebraic state at t0 (dimension nx-by-1)</td></tr>
<tr><td>INTEGRATOR_P</td><td>Parameters p (dimension np-by-1)</td></tr>
<tr><td>INTEGRATOR_XP0</td><td>State derivative at t0 (dimension nx-by-1) This input may be changed during an IDASIntegrator::evaluate()</td></tr>
</table>
*/
/** \class CasADi::QPSolverInternal
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::QPSolver
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::Interfaces::QPOasesInternal
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::Interfaces::QPOasesSolver
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::IpoptInternal
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::IpoptSolver
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesInternal
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::Sundials::CVodesIntegrator
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::Interfaces::IpoptQPInternal
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::Interfaces::IpoptQPSolver
<table>
<caption>Input scheme: CasADi::QPInput  (QP_NUM_IN = 9) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_H</td><td>The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.</td></tr>
<tr><td>QP_G</td><td>The column vector G: dense, (nx x 1)</td></tr>
<tr><td>QP_A</td><td>The matrix A: sparse, (nc x nx) - product with x must be dense.</td></tr>
<tr><td>QP_LBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_UBA</td><td>dense, (nc x 1)</td></tr>
<tr><td>QP_LBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_UBX</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_X_INIT</td><td>dense, (nx x 1)</td></tr>
<tr><td>QP_LAMBDA_INIT</td><td></td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::QPOutput  (QP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>QP_X_OPT</td><td>The optimal value of x as calculated with evaluate()</td></tr>
<tr><td>QP_COST</td><td>The value of the cost function as calculated with evaluate()</td></tr>
<tr><td>QP_LAMBDA_OPT</td><td></td></tr>
<tr><td>QP_LAMBDA_LBX</td><td></td></tr>
<tr><td>QP_LAMBDA_UBX</td><td></td></tr>
</table>
*/
/** \class CasADi::KnitroInternal
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::KnitroSolver
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptInternal
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::Interfaces::LiftoptSolver
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::IntegratorInternal
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::Integrator
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::CplexInternal
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::CplexSolver
<table>
<caption>Input scheme: CasADi::NLPInput  (NLP_NUM_IN = 7) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_INIT</td><td>Decision variables initial guess.</td></tr>
<tr><td>NLP_LBX</td><td>Decision variables lower bound.</td></tr>
<tr><td>NLP_UBX</td><td>Decision variables upper bound.</td></tr>
<tr><td>NLP_LBG</td><td>Constraints lower bound.</td></tr>
<tr><td>NLP_UBG</td><td>Constraints upper bound.</td></tr>
<tr><td>NLP_LAMBDA_INIT</td><td>Lambda multipliers initial guess.</td></tr>
<tr><td>NLP_P</td><td>Static parameters on which the objective and constraints might depend.</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::NLPOutput  (NLP_NUM_OUT = 5) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>NLP_X_OPT</td><td>Decision variables for optimal solution.</td></tr>
<tr><td>NLP_COST</td><td>Objective/cost function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_OPT</td><td>Lambda multipliers function for optimal solution.</td></tr>
<tr><td>NLP_LAMBDA_LBX</td><td>Lower bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
<tr><td>NLP_LAMBDA_UBX</td><td>Upper bound multipliers for optimal solution When in warm start mode, this output will be used as input</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasInternal
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::Sundials::IdasIntegrator
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::OCPSolverInternal
<table>
<caption>Input scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory.</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory.</td></tr>
<tr><td>OCP_XP_OPT</td><td>Optimal state derivative trajectory.</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters.</td></tr>
</table>
*/
/** \class CasADi::OCPSolver
<table>
<caption>Input scheme: CasADi::OCPOutput  (OCP_NUM_OUT = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>OCP_X_OPT</td><td>Optimal state trajectory.</td></tr>
<tr><td>OCP_U_OPT</td><td>Optimal control trajectory.</td></tr>
<tr><td>OCP_XP_OPT</td><td>Optimal state derivative trajectory.</td></tr>
<tr><td>OCP_P_OPT</td><td>Optimal parameters.</td></tr>
</table>
*/
/** \class CasADi::GSL::GslInternal
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
/** \class CasADi::GSL::GslIntegrator
<table>
<caption>Input scheme: CasADi::DAEInput  (DAE_NUM_IN = 4) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_T</td><td>Time. (1-by-1)</td></tr>
<tr><td>DAE_Y</td><td>State vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
<tr><td>DAE_P</td><td>Parameter vector (matrix).</td></tr>
<tr><td>DAE_YDOT</td><td>State derivative vector (matrix). Should have same amount of non-zeros as ODEOutput:DAE_RES</td></tr>
</table>
<br/><table>
<caption>Output scheme: CasADi::DAEOutput  (DAE_NUM_OUT = 1) </caption>
<tr><th>Name</th><th>Description</th></tr>
<tr><td>DAE_RES</td><td>Right hand side of ODE. Should have same amount of non-zeros as ODEInput:ODE_Y</td></tr>
</table>
*/
