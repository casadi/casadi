// Lqr.cpp
// Greg Horn
// Casadi 2011

#include "Lqr.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace CasADi;
using namespace std;

Lqr::Lqr(Ode & _ode, double t0_, double tf_, int N_, SX (*cost_)(map<string,SX> state, map<string,SX> action)) : ode(_ode) , V_0(N_) , V_x(N_) , V_xx(N_) , u_feedforward(N_-1) , u_feedback_gain(N_-1) , x_trajectory(N_) , u_trajectory(N_-1)
{
	ode.locked = 1;
	t0 = t0_;
	tf = tf_;
	N = N_;

	costFcnExt = cost_;

	// initialize trajectory
	for (int k=0; k<N; k++){
		V_0.at(k) =              DMatrix(        1,        1, 1.0);
		V_x.at(k) =              DMatrix( ode.nx(),        1, 2.0);
		V_xx.at(k) =             DMatrix( ode.nx(), ode.nx(), 3.0);
		x_trajectory.at(k) =     DMatrix( ode.nx(),        1, 4.0);

		if (k < N - 1){
			u_trajectory.at(k) =     DMatrix( ode.nu(),        1, 5.0);
			u_feedforward.at(k) =    DMatrix( ode.nu(),        1, 6.0);
			u_feedback_gain.at(k) =  DMatrix( ode.nu(), ode.nx(), 7.0);
		}

	}

	setupFunctions();
}

Lqr::~Lqr() {}

void Lqr::setupFunctions()
{
	// inputs
	SXMatrix xk = create_symbolic("xk", ode.nx());
	SXMatrix uk = create_symbolic("uk", ode.nu());

	vector<SXMatrix> inputs_xu(2);
	inputs_xu.at(0) = xk;
	inputs_xu.at(1) = uk;

	map<string,SX> xkMap = ode.getStateMap(xk);
	map<string,SX> ukMap = ode.getActionMap(uk);

	/**************** cost *********************/
	// function
	SXMatrix cost(costFcnExt( xkMap, ukMap ));

	// jacobian
	SXMatrix cost_x = gradient( cost, xk );
	SXMatrix cost_u = gradient( cost, uk );

	// hessian
	SXMatrix cost_xx = jacobian( cost_x, xk );
	SXMatrix cost_xu = jacobian( cost_x, uk ); // == jacobian( cost_u, x ).trans()
	SXMatrix cost_uu = jacobian( cost_u, uk );

	/**************** dynamics *********************/
	// dummy params for now
	map<string,SX> dummyParams;
	double dt = (tf - t0)/(N - 1);
	cout << "SWITCH BACK TO RK4 WHEN YOU GET IT WORKING\n";
	//SXMatrix f = ode.rk4Step( xk, uk, uk, dummyParams, t0, dt); // timestep dependent: f(x,u,__t__) - same w cost
	SXMatrix f = ode.eulerStep( xk, uk, dummyParams, SX(t0), SX(dt)); // timestep dependent:  f(x,u,__t__) - same w cost
	SXMatrix f_x = jacobian( f, xk );
	SXMatrix f_u = jacobian( f, uk );


	/**************** Q function *********************/
	SXMatrix V_0_kp1( create_symbolic("V_0_kp1", 1, 1) );
	SXMatrix V_x_kp1( create_symbolic("V_x_kp1", ode.nx(), 1) );
	SXMatrix V_xx_kp1( create_symbolic("V_xx_kp1", ode.nx(), ode.nx()) );

	// Q_0 = 1/2*f'*V_xx*f + V_x_kp1'*f + cost + V_0_kp1
	SXMatrix Q_0 = prod( f.trans(), prod(V_xx_kp1, f) )/2 // 1/2*f'*V_xx*f (quadratic term)
		+ prod(V_x_kp1.trans(), f) // V_x_kp1'*f (linear term)
		+ cost + V_0_kp1; // cost + V_0_kp1 (constant term)

	// Q_x = cost_x + f_x'*( V_x_kp1 + V_xx_kp1*f )
	SXMatrix Q_x = cost_x + prod(f_x.trans(), V_x_kp1 + prod( V_xx_kp1, f ) );
	
	// Q_u = cost_u + f_u'*( V_x_kp1 + V_xx_kp1*f )
	SXMatrix Q_u = cost_u + prod(f_u.trans(), V_x_kp1 + prod( V_xx_kp1, f ) );
	
	// Q_xx = cost_xx + f_x'*V_xx_kp1*f_x
	SXMatrix Q_xx = cost_xx + prod( f_x.trans(), prod( V_xx_kp1, f_x ) );

	// Q_xu = cost_xu + f_x'*V_xx_kp1*f_u
	SXMatrix Q_xu = cost_xu + prod( f_x.trans(), prod( V_xx_kp1, f_u ) );

	// Q_uu = cost_uu + f_u'*V_xx_kp1*f_u
	SXMatrix Q_uu = cost_uu + prod( f_u.trans(), prod( V_xx_kp1, f_u ) );


	/************** optimal control *******************/
	SXMatrix Q_uu_inv = inv( Q_uu );
	SXMatrix u_feedforward_k   = -prod( Q_uu_inv, Q_u );
	SXMatrix u_feedback_gain_k = -prod( Q_uu_inv, Q_xu.trans() );


	/************** value function propogation ********/
	SXMatrix V_0_k  = Q_0  - prod( Q_u.trans(), prod( Q_uu_inv, Q_u ) );
	SXMatrix V_x_k  = Q_x  - prod( Q_xu, prod( Q_uu_inv.trans(), Q_u ) );
	SXMatrix V_xx_k = Q_xx - prod( Q_xu, prod( Q_uu_inv, Q_xu.trans() ) );


	/*************** functions ****************/
	vector<SXMatrix> inputs(LQR_NUM_INPUTS);
	inputs.at(IDX_INPUTS_X_K)      = xk;
	inputs.at(IDX_INPUTS_U_K)      = uk;
	inputs.at(IDX_INPUTS_V_0_KP1)  = V_0_kp1;
	inputs.at(IDX_INPUTS_V_X_KP1)  = V_x_kp1;
	inputs.at(IDX_INPUTS_V_XX_KP1) = V_xx_kp1;

	vector<SXMatrix> ilqr_outputs(LQR_NUM_OUTPUTS);
	ilqr_outputs.at(IDX_LQR_OUTPUTS_U_FEEDFORWARD_K)   = u_feedforward_k;
	ilqr_outputs.at(IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K) = u_feedback_gain_k;
	ilqr_outputs.at(IDX_LQR_OUTPUTS_V_0_K)             = V_0_k;
	ilqr_outputs.at(IDX_LQR_OUTPUTS_V_X_K)             = V_x_k;
	ilqr_outputs.at(IDX_LQR_OUTPUTS_V_XX_K)            = V_xx_k;

	ilqr_fcn = SXFunction( inputs, ilqr_outputs );
	ilqr_fcn.init();

	// vector<SXMatrix> cost_outputs(6);
	// cost_outputs.at(0) = cost;
	// cost_outputs.at(1) = cost_x;
	// cost_outputs.at(2) = cost_u;
	// cost_outputs.at(3) = cost_xx;
	// cost_outputs.at(4) = cost_xu;
	// cost_outputs.at(5) = cost_uu;

	// cost_fcn = SXFunction( inputs, cost_outputs );
	// cost_fcn.init();

	// vector<SXMatrix> dynamics_outputs(3);
	// dynamics_outputs.at(0) = f;
	// dynamics_outputs.at(1) = f_x;
	// dynamics_outputs.at(2) = f_u;

	// dynamics_fcn = SXFunction( inputs, dynamics_outputs );
	// dynamics_fcn.init();

	// vector<SXMatrix> Q_outputs(6);
	// Q_outputs.at(0) = Q_0;
	// Q_outputs.at(1) = Q_x;
	// Q_outputs.at(2) = Q_u;
	// Q_outputs.at(3) = Q_xx;
	// Q_outputs.at(4) = Q_xu;
	// Q_outputs.at(5) = Q_uu;

	// Q_fcn = SXFunction( inputs, Q_outputs );
	// Q_fcn.init();

	// cout << "cost:\n" << cost << endl;
	// cout << "cost_x:\n" << cost_x << endl;
	// cout << "cost_u:\n" << cost_u << endl;
	// cout << "cost_xx:\n" << cost_xx << endl;
	// cout << "cost_xu:\n" << cost_xu << endl;
	// cout << "cost_uu:\n" << cost_uu << endl;


	// cout << "ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDFORWARD_K).size1():   " << ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDFORWARD_K).size1() << endl;
	// cout << "ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDFORWARD_K).size2():   " << ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDFORWARD_K).size2() << endl;
	// cout << "ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K).size1(): " << ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K).size1() << endl;
	// cout << "ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K).size2(): " << ilqr_fcn.outputSX(IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K).size2() << endl;


	// cout << "u_feedback_gain_k:\n";
	// cout << u_feedback_gain_k << endl << endl;
	// cout << "u_feedback_gain_k.size1(): " << u_feedback_gain_k.size1() << endl;
	// cout << "u_feedback_gain_k.size2(): " << u_feedback_gain_k.size2() << endl;
}

void Lqr::runBackwardSweep()
{
	/*********** fill in random values for V(N-1) ************/
	V_0.at(N-1)[0,0] = 1.5;

	for (int k=0; k<ode.nx(); k++)
		V_x.at(N-1)[0,k] = 1.7*k+0.3;

	int counter = 0;
	for (int k=0; k<ode.nx(); k++)
		for (int j=0; j<ode.nx(); j++){
			V_xx.at(N-1).at(counter) = 0.5 + (j+1)*(k+1)*0.7;
			counter++;
		}

	// cout << endl;
	// cout << "V_0.at(N-1):\n" << V_0.at(N-1) << endl << endl;
	// cout << "V_x.at(N-1):\n" << V_x.at(N-1) << endl << endl;
	// cout << "V_xx.at(N-1):\n" << V_xx.at(N-1) << endl << endl;

	/*********** run backward sweep ************/
	for (int k = N-2; k >= 0; k--)
		takeBackwardStep(k);
}

void Lqr::runForwardSweep()
{
	// for (int k = 0; k < N-1; k++)
	// 	takeForwardStep(k);
}

void Lqr::takeBackwardStep(int timestep)
{
	// set inputs
	ilqr_fcn.setInput( x_trajectory.at( timestep     ), IDX_INPUTS_X_K      );
	ilqr_fcn.setInput( u_trajectory.at( timestep     ), IDX_INPUTS_U_K      );
	ilqr_fcn.setInput(          V_0.at( timestep + 1 ), IDX_INPUTS_V_0_KP1  );
	ilqr_fcn.setInput(          V_x.at( timestep + 1 ), IDX_INPUTS_V_X_KP1  );
	ilqr_fcn.setInput(         V_xx.at( timestep + 1 ), IDX_INPUTS_V_XX_KP1 );

	// set outputs
	ilqr_fcn.setOutput(   u_feedforward.at(timestep), IDX_LQR_OUTPUTS_U_FEEDFORWARD_K   );
	ilqr_fcn.setOutput( u_feedback_gain.at(timestep), IDX_LQR_OUTPUTS_U_FEEDBACK_GAIN_K );
	ilqr_fcn.setOutput(             V_0.at(timestep), IDX_LQR_OUTPUTS_V_0_K             );
	ilqr_fcn.setOutput(             V_x.at(timestep), IDX_LQR_OUTPUTS_V_X_K             );
	ilqr_fcn.setOutput(            V_xx.at(timestep), IDX_LQR_OUTPUTS_V_XX_K            );

	// evaluate
	ilqr_fcn.evaluate();
}
