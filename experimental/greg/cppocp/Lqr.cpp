// Lqr.cpp
// Greg Horn
// Casadi 2011

#include "Lqr.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace CasADi;
using namespace std;

Lqr::Lqr(Ode & _ode, double t0_, double tf_, int N_, SX (*cost_)(map<string,SX> state, map<string,SX> action, int timestep, int _N)) : ode(_ode) , V_0(N_) , V_x(N_) , V_xx(N_) , u_feedforward(N_) , u_feedback_gain(N_) , x_trajectory(N_) , u_trajectory(N_), cost_0(N_), cost_x(N_), cost_u(N_), cost_xx(N_), cost_xu(N_), cost_uu(N_)
{
	ode.init();

	t0 = t0_;
	tf = tf_;
	N = N_;

	costFcnExt = cost_;

	// initialize trajectory
	for (int k=0; k<N; k++){
		V_0.at(k)               = DMatrix(        1,        1, 0.0);
		V_x.at(k)               = DMatrix( ode.nx(),        1, 0.0);
		V_xx.at(k)              = DMatrix( ode.nx(), ode.nx(), 0.0);
		x_trajectory.at(k)      = DMatrix( ode.nx(),        1, 0.0);

		cost_0.at(k)  = DMatrix(        1,        1, 0.0);
		cost_x.at(k)  = DMatrix( ode.nx(),        1, 0.0);
		cost_u.at(k)  = DMatrix( ode.nu(),        1, 0.0);
		cost_xx.at(k) = DMatrix( ode.nx(), ode.nx(), 0.0);
		cost_xu.at(k) = DMatrix( ode.nx(), ode.nu(), 0.0);
		cost_uu.at(k) = DMatrix( ode.nu(), ode.nu(), 0.0);

		u_trajectory.at(k)    = DMatrix( ode.nu(),        1, 0.0);
		u_feedforward.at(k)   = DMatrix( ode.nu(),        1, 0.0);
		u_feedback_gain.at(k) = DMatrix( ode.nu(), ode.nx(), 0.0);
	}

	setupBackwardSweepFunction();
	setupCostFunctions();
}

Lqr::~Lqr() {}


void Lqr::setupCostFunctions()
{
	// inputs
	SXMatrix xk( create_symbolic( "xk", ode.nx(), 1 ) );
	SXMatrix uk( create_symbolic( "uk", ode.nu(), 1 ) );

	vector<SXMatrix> costInputs(NUM_COST_INPUTS);
	costInputs.at(IDX_COST_INPUTS_X_K) = xk;
	costInputs.at(IDX_COST_INPUTS_U_K) = uk;

	map<string,SX> xkMap = ode.getStateMap(xk);
	map<string,SX> ukMap = ode.getActionMap(uk);

	for (int k=0; k<N; k++){

		// function
		SXMatrix cost_0_k(costFcnExt( xkMap, ukMap, k, N ));

		// jacobian
		SXMatrix cost_x_k = gradient( cost_0_k, xk );
		SXMatrix cost_u_k = gradient( cost_0_k, uk );

		// hessian
		SXMatrix cost_xx_k = jacobian( cost_x_k, xk );
		SXMatrix cost_xu_k = jacobian( cost_x_k, uk ); // == jacobian( cost_u, x ).trans()
		SXMatrix cost_uu_k = jacobian( cost_u_k, uk );

		// workaround bug where size() != size1()*size2()
		makeDense(cost_0_k);
		makeDense(cost_x_k);
		makeDense(cost_u_k);
		makeDense(cost_xx_k);
		makeDense(cost_xu_k);
		makeDense(cost_uu_k);

		// outputs
		vector<SXMatrix> costOutputs(NUM_COST_OUTPUTS);
		costOutputs.at(IDX_COST_OUTPUTS_COST_0_K)  = cost_0_k;
		costOutputs.at(IDX_COST_OUTPUTS_COST_X_K)  = cost_x_k;
		costOutputs.at(IDX_COST_OUTPUTS_COST_U_K)  = cost_u_k;
		costOutputs.at(IDX_COST_OUTPUTS_COST_XX_K) = cost_xx_k;
		costOutputs.at(IDX_COST_OUTPUTS_COST_XU_K) = cost_xu_k;
		costOutputs.at(IDX_COST_OUTPUTS_COST_UU_K) = cost_uu_k;
		
		// sx function
		costFunctions.push_back( SXFunction( costInputs, costOutputs ) );
		costFunctions.at(k).init();
	}
}


void Lqr::setupBackwardSweepFunction()
{
	/*************** inputs **************/
	SXMatrix xk = create_symbolic("xk", ode.nx());
	SXMatrix uk = create_symbolic("uk", ode.nu());

	SXMatrix cost_0_k(  create_symbolic(  "cost_0",        1,        1 ) );
	SXMatrix cost_x_k(  create_symbolic(  "cost_x", ode.nx(),        1 ) );
	SXMatrix cost_u_k(  create_symbolic(  "cost_u", ode.nu(),        1 ) );
	SXMatrix cost_xx_k( create_symbolic( "cost_xx", ode.nx(), ode.nx() ) );
	SXMatrix cost_xu_k( create_symbolic( "cost_xu", ode.nx(), ode.nu() ) );
	SXMatrix cost_uu_k( create_symbolic( "cost_uu", ode.nu(), ode.nu() ) );

	SXMatrix V_0_kp1(  create_symbolic(  "V_0_kp1",        1,        1) );
	SXMatrix V_x_kp1(  create_symbolic(  "V_x_kp1", ode.nx(),        1) );
	SXMatrix V_xx_kp1( create_symbolic( "V_xx_kp1", ode.nx(), ode.nx()) );

	vector<SXMatrix> backwardSweepInputs(NUM_BACKWARD_SWEEP_INPUTS);
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_X_K)       = xk;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_U_K)       = uk;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_0_K)  = cost_0_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_X_K)  = cost_x_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_U_K)  = cost_u_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K) = cost_xx_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K) = cost_xu_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K) = cost_uu_k;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1)   = V_0_kp1;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1)   = V_x_kp1;
	backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1)  = V_xx_kp1;


	/**************** dynamics *********************/
	// dummy params for now
	map<string,SX> dummyParams;
	double dt = (tf - t0)/(N - 1);
	cout << "SWITCH BACK TO RK4 WHEN YOU GET IT WORKING\n";
	SXMatrix f = ode.rk4Step( xk, uk, uk, dummyParams, t0, dt); // timestep dependent: f(x,u,__t__) - same w cost
	// SXMatrix f = ode.eulerStep( xk, uk, dummyParams, SX(t0), SX(dt)); // timestep dependent:  f(x,u,__t__) - same w cost
	SXMatrix f_x = jacobian( f, xk );
	SXMatrix f_u = jacobian( f, uk );


	/**************** Q function *********************/
	// Q_0 = 1/2*f'*V_xx*f + V_x_kp1'*f + cost_0 + V_0_kp1
	SXMatrix Q_0 = prod( f.trans(), prod(V_xx_kp1, f) )/2 // 1/2*f'*V_xx*f (quadratic term)
		+ prod(V_x_kp1.trans(), f) // V_x_kp1'*f (linear term)
		+ cost_0_k + V_0_kp1; // cost_0 + V_0_kp1 (constant term)

	// Q_x = cost_x + f_x'*( V_x_kp1 + V_xx_kp1*f )
	SXMatrix Q_x = cost_x_k + prod(f_x.trans(), V_x_kp1 + prod( V_xx_kp1, f ) );
	
	// Q_u = cost_u + f_u'*( V_x_kp1 + V_xx_kp1*f )
	SXMatrix Q_u = cost_u_k + prod(f_u.trans(), V_x_kp1 + prod( V_xx_kp1, f ) );
	
	// Q_xx = cost_xx + f_x'*V_xx_kp1*f_x
	SXMatrix Q_xx = cost_xx_k + prod( f_x.trans(), prod( V_xx_kp1, f_x ) );

	// Q_xu = cost_xu + f_x'*V_xx_kp1*f_u
	SXMatrix Q_xu = cost_xu_k + prod( f_x.trans(), prod( V_xx_kp1, f_u ) );

	// Q_uu = cost_uu + f_u'*V_xx_kp1*f_u
	SXMatrix Q_uu = cost_uu_k + prod( f_u.trans(), prod( V_xx_kp1, f_u ) );


	/************** optimal control *******************/
	SXMatrix Q_uu_inv = inv( Q_uu );
	SXMatrix u_feedforward_k   = -prod( Q_uu_inv, Q_u );
	SXMatrix u_feedback_gain_k = -prod( Q_uu_inv, Q_xu.trans() );


	/************** value function propogation ********/
	SXMatrix V_0_k  = Q_0  - prod( Q_u.trans(), prod( Q_uu_inv, Q_u ) );
	SXMatrix V_x_k  = Q_x  - prod( Q_xu, prod( Q_uu_inv.trans(), Q_u ) );
	SXMatrix V_xx_k = Q_xx - prod( Q_xu, prod( Q_uu_inv, Q_xu.trans() ) );


	/*************** functions ****************/
	// workaround bug where size() != size1()*size2()
	makeDense( u_feedforward_k );
	makeDense( u_feedback_gain_k );
	makeDense( V_0_k );
	makeDense( V_x_k );
	makeDense( V_xx_k );

	vector<SXMatrix> backwardSweepOutputs(NUM_BACKWARD_SWEEP_OUTPUTS);
	backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDFORWARD_K)   = u_feedforward_k;
	backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDBACK_GAIN_K) = u_feedback_gain_k;
	backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K)             = V_0_k;
	backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K)             = V_x_k;
	backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K)            = V_xx_k;

	backwardSweepFcn = SXFunction( backwardSweepInputs, backwardSweepOutputs );
	backwardSweepFcn.init();

	// workaround bug where size() != size1()*size2()
	makeDense( f );
	makeDense( f_x );
	makeDense( f_u );
	vector<SXMatrix> dynamicsOutputs(3);
	dynamicsOutputs.at(0) = f;
	dynamicsOutputs.at(1) = f_x;
	dynamicsOutputs.at(2) = f_u;

	dynamicsFcn = SXFunction( backwardSweepInputs, dynamicsOutputs );
	dynamicsFcn.init();

	// workaround bug where size() != size1()*size2()
	makeDense( Q_0 );
	makeDense( Q_x );
	makeDense( Q_u );
	makeDense( Q_xx );
	makeDense( Q_xu );
	makeDense( Q_uu );

	vector<SXMatrix> qOutputs(6);
	qOutputs.at(0) = Q_0;
	qOutputs.at(1) = Q_x;
	qOutputs.at(2) = Q_u;
	qOutputs.at(3) = Q_xx;
	qOutputs.at(4) = Q_xu;
	qOutputs.at(5) = Q_uu;

	qFcn = SXFunction( backwardSweepInputs, qOutputs );
	qFcn.init();
}

void Lqr::runBackwardSweep()
{
	cout << "here yo\n";
	/******************* V(N-1) == Cost(N-1) ******************/
	// set inputs
	costFunctions.at( N-1 ).setInput( x_trajectory.at( N-1 ), IDX_COST_INPUTS_X_K );
	costFunctions.at( N-1 ).setInput( u_trajectory.at( N-1 ), IDX_COST_INPUTS_U_K );

	// evaluate
	costFunctions.at( N-1 ).evaluate();

	cout << "getting output\n";

	// get outputs
	costFunctions.at( N-1 ).getOutput(  V_0.at( N-1 ), IDX_COST_OUTPUTS_COST_0_K  );
	costFunctions.at( N-1 ).getOutput(  V_x.at( N-1 ), IDX_COST_OUTPUTS_COST_X_K  );
	costFunctions.at( N-1 ).getOutput( V_xx.at( N-1 ), IDX_COST_OUTPUTS_COST_XX_K );

	cout << "done yo\n";
	/*********** run backward sweep ************/
	for (int k = N-2; k >= 0; k--)
		takeBackwardStep(k);

	// // print value function along trajectory
	// for (int k = N-1; k >= 0; k--){
	// 	cout << endl;
	// 	cout << "V_0.at("  << k << "):\n" << V_0.at(k)  << endl << endl;
	// 	cout << "V_x.at("  << k << "):\n" << V_x.at(k)  << endl << endl;
	// 	cout << "V_xx.at(" << k << "):\n" << V_xx.at(k) << endl << endl;
	// }
}

void Lqr::runForwardSweep()
{
	// for (int k = 0; k < N-1; k++)
	// 	takeForwardStep(k);
}

void Lqr::takeBackwardStep(int timestep)
{
	/******************* evaluate cost function ******************/
	// set inputs
	costFunctions.at(timestep).setInput( x_trajectory.at( timestep ), IDX_COST_INPUTS_X_K );
	costFunctions.at(timestep).setInput( u_trajectory.at( timestep ), IDX_COST_INPUTS_U_K );

	// evaluate
	costFunctions.at(timestep).evaluate();

	// get outputs
	costFunctions.at(timestep).getOutput(  cost_0.at(timestep), IDX_COST_OUTPUTS_COST_0_K  );
	costFunctions.at(timestep).getOutput(  cost_x.at(timestep), IDX_COST_OUTPUTS_COST_X_K  );
	costFunctions.at(timestep).getOutput(  cost_u.at(timestep), IDX_COST_OUTPUTS_COST_U_K  );
	costFunctions.at(timestep).getOutput( cost_xx.at(timestep), IDX_COST_OUTPUTS_COST_XX_K );
	costFunctions.at(timestep).getOutput( cost_xu.at(timestep), IDX_COST_OUTPUTS_COST_XU_K );
	costFunctions.at(timestep).getOutput( cost_uu.at(timestep), IDX_COST_OUTPUTS_COST_UU_K );


	/*************** evaluate backward sweep function ************/
	// set inputs
	backwardSweepFcn.setInput( x_trajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_X_K       );
	backwardSweepFcn.setInput( u_trajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_U_K       );
	backwardSweepFcn.setInput(       cost_0.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_0_K  );
	backwardSweepFcn.setInput(       cost_x.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_X_K  );
	backwardSweepFcn.setInput(       cost_u.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_U_K  );
	backwardSweepFcn.setInput(      cost_xx.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K );
	backwardSweepFcn.setInput(      cost_xu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K );
	backwardSweepFcn.setInput(      cost_uu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K );
	backwardSweepFcn.setInput(          V_0.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1  );
	backwardSweepFcn.setInput(          V_x.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1  );
	backwardSweepFcn.setInput(         V_xx.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1 );

	// set outputs
	// backwardSweepFcn.setOutput(   u_feedforward.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDFORWARD_K   );
	// backwardSweepFcn.setOutput( u_feedback_gain.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDBACK_GAIN_K );
	// backwardSweepFcn.setOutput(             V_0.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K             );
	// backwardSweepFcn.setOutput(             V_x.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K             );
	// backwardSweepFcn.setOutput(            V_xx.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K            );

	// evaluate
	backwardSweepFcn.evaluate();

	// get outputs
	backwardSweepFcn.getOutput(   u_feedforward.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDFORWARD_K   );
	backwardSweepFcn.getOutput( u_feedback_gain.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_U_FEEDBACK_GAIN_K );
	backwardSweepFcn.getOutput(             V_0.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K             );
	backwardSweepFcn.getOutput(             V_x.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K             );
	backwardSweepFcn.getOutput(            V_xx.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K            );
}
