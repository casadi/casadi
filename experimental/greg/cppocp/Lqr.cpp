/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// Lqr.cpp
// Greg Horn
// Casadi 2011

#include "Lqr.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace casadi;
using namespace std;

Lqr::Lqr(Ode & _ode, double t0_, double tf_, int N_, SX (*cost_)(map<string,SX> state, map<string,SX> action, int timestep, int _N)) : ode(_ode)
{
     ode.init();
     ode.setupIntegrators();

     t0 = t0_;
     tf = tf_;
     N = N_;

     costFcnExt = cost_;

     // initialize trajectory
     for (int k=0; k<N; k++){
	  V_0.push_back(   DMatrix(        1,        1, 0.0) );
	  V_x.push_back(   DMatrix( ode.nx(),        1, 0.0) );
	  V_xx.push_back(  DMatrix( ode.nx(), ode.nx(), 0.0) );

	  cost_0.push_back(  DMatrix(        1,        1, 0.0) );
	  cost_x.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  cost_u.push_back(  DMatrix( ode.nu(),        1, 0.0) );
	  cost_xx.push_back( DMatrix( ode.nx(), ode.nx(), 0.0) );
	  cost_xu.push_back( DMatrix( ode.nx(), ode.nu(), 0.0) );
	  cost_uu.push_back( DMatrix( ode.nu(), ode.nu(), 0.0) );

	  xTrajectory.push_back(         DMatrix( ode.nx(),        1, 0.0) );
	  xNominalTrajectory.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  uTrajectory.push_back(         DMatrix( ode.nu(),        1, 0.0) );
	  uOpenLoop.push_back(           DMatrix( ode.nu(),        1, 0.0) );
	  feedbackGain.push_back(        DMatrix( ode.nu(), ode.nx(), 0.0) );

	  Q0Trajectory.push_back(  DMatrix(        1,        1, 0.0) );
	  QxTrajectory.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  QuTrajectory.push_back(  DMatrix( ode.nu(),        1, 0.0) );
	  QxxTrajectory.push_back( DMatrix( ode.nx(), ode.nx(), 0.0) );
	  QxuTrajectory.push_back( DMatrix( ode.nx(), ode.nu(), 0.0) );
	  QuuTrajectory.push_back( DMatrix( ode.nu(), ode.nu(), 0.0) );
		
	  f0Trajectory.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  functionTrajectory.push_back(  DMatrix( ode.nx(), ode.nx(), 0.0) );
	  fuTrajectory.push_back(  DMatrix( ode.nx(), ode.nu(), 0.0) );

	  // initialize action bounds
	  vector<double> ubAction_;
	  vector<double> lbAction_;
	  for (int j=0; j<ode.nu(); j++){
	       ubAction_.push_back(1e30);
	       lbAction_.push_back(-1e30);
	  }
	  ubAction.push_back( ubAction_ );
	  lbAction.push_back( lbAction_ );
     }


     stateRegularization  = DMatrix( ode.nx(), ode.nx(), 0.0 );
     actionRegularization = DMatrix( ode.nu(), ode.nu(), 0.0 );

     setupBackwardSweepFunction();
     setupCostFunctions();
}

Lqr::~Lqr() {}


void Lqr::setupCostFunctions()
{
     // inputs
     SX xk( ssym( "xk", ode.nx(), 1 ) );
     SX uk( ssym( "uk", ode.nu(), 1 ) );

     vector<SX> costInputs(NUM_COST_INPUTS);
     costInputs.at(IDX_COST_INPUTS_X_K) = xk;
     costInputs.at(IDX_COST_INPUTS_U_K) = uk;

     map<string,SX> xkMap = ode.getStateMap(xk);
     map<string,SX> ukMap = ode.getActionMap(uk);

     for (int k=0; k<N; k++){

	  // function
	  SX cost_0_k(costFcnExt( xkMap, ukMap, k, N ));

	  // jacobian
	  SX cost_x_k = gradient( cost_0_k, xk );
	  SX cost_u_k = gradient( cost_0_k, uk );

	  // hessian
	  SX cost_xx_k = jacobian( cost_x_k, xk );
	  SX cost_xu_k = jacobian( cost_x_k, uk ); // == jacobian( cost_u, x ).trans()
	  SX cost_uu_k = jacobian( cost_u_k, uk );

	  simplify(cost_0_k);
	  simplify(cost_x_k);
	  simplify(cost_u_k);
	  simplify(cost_xx_k);
	  simplify(cost_xu_k);
	  simplify(cost_uu_k);

	  // workaround bug where size() != size1()*size2()
	  makeDense(cost_0_k);
	  makeDense(cost_x_k);
	  makeDense(cost_u_k);
	  makeDense(cost_xx_k);
	  makeDense(cost_xu_k);
	  makeDense(cost_uu_k);

	  // outputs
	  vector<SX> costOutputs(NUM_COST_OUTPUTS);
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
     SX xk = ssym("xk", ode.nx());
     SX uk = ssym("uk", ode.nu());

     SX cost_0_k(  ssym(  "cost_0",        1,        1 ) );
     SX cost_x_k(  ssym(  "cost_x", ode.nx(),        1 ) );
     SX cost_u_k(  ssym(  "cost_u", ode.nu(),        1 ) );
     SX cost_xx_k( ssym( "cost_xx", ode.nx(), ode.nx() ) );
     SX cost_xu_k( ssym( "cost_xu", ode.nx(), ode.nu() ) );
     SX cost_uu_k( ssym( "cost_uu", ode.nu(), ode.nu() ) );

     SX V_0_kp1(  ssym(  "V_0_kp1",        1,        1) );
     SX V_x_kp1(  ssym(  "V_x_kp1", ode.nx(),        1) );
     SX V_xx_kp1( ssym( "V_xx_kp1", ode.nx(), ode.nx()) );

     vector<SX> backwardSweepInputs(NUM_BACKWARD_SWEEP_INPUTS);
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
     SX f = ode.rk4Step( xk, uk, uk, dummyParams, t0, dt); // timestep dependent: f(x,u,__t__) - same as cost
     //SX f = ode.eulerStep( xk, uk, dummyParams, SX(t0), SX(dt)); // timestep dependent:  f(x,u,__t__) - same as cost
     SX f_x = jacobian( f, xk );
     SX f_u = jacobian( f, uk );


     /**************** Q function *********************/
     // Q_0 = cost_0 + V_0_kp1
     SX Q_0 = cost_0_k + V_0_kp1;

     // Q_x = cost_x + f_x'*V_x_kp1
     SX Q_x = cost_x_k + mul( f_x.trans(), V_x_kp1 );
	
     // Q_u = cost_u + f_u'* V_x_kp1
     SX Q_u = cost_u_k + mul( f_u.trans(), V_x_kp1 );
	
     // Q_xx = cost_xx + f_x'*V_xx_kp1*f_x
     SX Q_xx = cost_xx_k + mul( f_x.trans(), mul( V_xx_kp1, f_x ) );

     // Q_xu = cost_xu + f_x'*V_xx_kp1*f_u
     SX Q_xu = cost_xu_k + mul( f_x.trans(), mul( V_xx_kp1, f_u ) );

     // Q_uu = cost_uu + f_u'*V_xx_kp1*f_u
     SX Q_uu = cost_uu_k + mul( f_u.trans(), mul( V_xx_kp1, f_u ) );


     /************** optimal control *******************/
     SX Q_uu_inv = inv( Q_uu );
     SX u_feedforward_k = -mul( Q_uu_inv, Q_u );
     SX feedbackGain_k  = -mul( Q_uu_inv, Q_xu.trans() );


     /************** value function propogation ********/
     SX V_0_k  = Q_0  - mul( Q_u.trans(), mul( Q_uu_inv, Q_u ) );
     SX V_x_k  = Q_x  - mul( Q_xu, mul( Q_uu_inv.trans(), Q_u ) );
     SX V_xx_k = Q_xx - mul( Q_xu, mul( Q_uu_inv, Q_xu.trans() ) );


     /*************** backwardSweepFcn ****************/
     // workaround bug where size() != size1()*size2()
     makeDense( u_feedforward_k );
     makeDense( feedbackGain_k );
     makeDense( V_0_k );
     makeDense( V_x_k );
     makeDense( V_xx_k );

     vector<SX> backwardSweepOutputs(NUM_BACKWARD_SWEEP_OUTPUTS);
     backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_U_OPEN_LOOP_K)   = u_feedforward_k + uk;
     backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_FEEDBACK_GAIN_K) = feedbackGain_k;
     backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K)           = V_0_k;
     backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K)           = V_x_k;
     backwardSweepOutputs.at(IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K)          = V_xx_k;

     backwardSweepFcn = SXFunction( backwardSweepInputs, backwardSweepOutputs );
     backwardSweepFcn.init();

     /*************** dynamicsFcn ****************/
     // workaround bug where size() != size1()*size2()
     makeDense( f );
     makeDense( f_x );
     makeDense( f_u );
     vector<SX> dynamicsOutputs(3);
     dynamicsOutputs.at(0) = f;
     dynamicsOutputs.at(1) = f_x;
     dynamicsOutputs.at(2) = f_u;

     dynamicsFcn = SXFunction( backwardSweepInputs, dynamicsOutputs );
     dynamicsFcn.init();

     /*************** qFcn ****************/
     // workaround bug where size() != size1()*size2()
     makeDense( Q_0 );
     makeDense( Q_x );
     makeDense( Q_u );
     makeDense( Q_xx );
     makeDense( Q_xu );
     makeDense( Q_uu );

     // outputs
     vector<SX> qOutputs(NUM_COST_OUTPUTS);
     qOutputs.at(IDX_COST_OUTPUTS_COST_0_K)  = Q_0;
     qOutputs.at(IDX_COST_OUTPUTS_COST_X_K)  = Q_x;
     qOutputs.at(IDX_COST_OUTPUTS_COST_U_K)  = Q_u;
     qOutputs.at(IDX_COST_OUTPUTS_COST_XX_K) = Q_xx;
     qOutputs.at(IDX_COST_OUTPUTS_COST_XU_K) = Q_xu;
     qOutputs.at(IDX_COST_OUTPUTS_COST_UU_K) = Q_uu;
	
     // sx function
     qFcn = SXFunction( backwardSweepInputs, qOutputs );
     qFcn.init();
}

void Lqr::runBackwardSweep()
{
     xNominalTrajectory.at(N-1) = xTrajectory.at(N-1);

     /******************* V(N-1) == Cost(N-1) ******************/
     // set inputs
     costFunctions.at( N-1 ).setInput( xTrajectory.at( N-1 ), IDX_COST_INPUTS_X_K );
     costFunctions.at( N-1 ).setInput( uTrajectory.at( N-1 ), IDX_COST_INPUTS_U_K );

     // evaluate
     costFunctions.at( N-1 ).evaluate();

     // get outputs
     costFunctions.at( N-1 ).getOutput(  V_0.at( N-1 ), IDX_COST_OUTPUTS_COST_0_K  );
     costFunctions.at( N-1 ).getOutput(  V_x.at( N-1 ), IDX_COST_OUTPUTS_COST_X_K  );
     costFunctions.at( N-1 ).getOutput( V_xx.at( N-1 ), IDX_COST_OUTPUTS_COST_XX_K );

     /*********** run backward sweep ************/
     for (int k = N-2; k >= 0; k--)
	  takeBackwardStep(k);

     // print stuff along trajectory
     // for (int k = N-1; k >= 0; k--){
     // 	cout << "========================= timestep: " << k << " ======================\n";
     // 	cout << "---------- cost: ----------\n";
     // 	cout << " cost_0.at(" << k << "): "  <<  cost_0.at(k) << endl;
     // 	cout << " cost_x.at(" << k << "): "  <<  cost_x.at(k) << endl;
     // 	cout << " cost_u.at(" << k << "): "  <<  cost_u.at(k) << endl;
     // 	cout << "cost_xx.at(" << k << "):\n" << cost_xx.at(k);
     // 	cout << "cost_xu.at(" << k << "): "  << cost_xu.at(k) << endl;
     // 	cout << "cost_uu.at(" << k << "): "  << cost_uu.at(k) << endl;
     // 	cout << endl;

     // 	cout << "---------- dynamics: --------\n";
     // 	cout << "f0Trajectory.at(" << k << "): "  <<  f0Trajectory.at(k) << endl;
     // 	cout << "functionTrajectory.at(" << k << "):\n" <<  functionTrajectory.at(k);
     // 	cout << "fuTrajectory.at(" << k << "): "  <<  fuTrajectory.at(k) << endl;
     // 	cout << endl;

     // 	cout << "---------- Q: ----------\n";
     // 	cout << " Q0Trajectory(" << k << "): "  <<  Q0Trajectory.at(k) << endl;
     // 	cout << " QxTrajectory(" << k << "): "  <<  QxTrajectory.at(k) << endl;
     // 	cout << " QuTrajectory(" << k << "): "  <<  QuTrajectory.at(k) << endl;
     // 	cout << "QxxTrajectory(" << k << "):\n" << QxxTrajectory.at(k);
     // 	cout << "QxuTrajectory(" << k << "): "  << QxuTrajectory.at(k) << endl;
     // 	cout << "QuuTrajectory(" << k << "): "  << QuuTrajectory.at(k) << endl;
     // 	cout << endl;

     // 	cout << "---------- V(" << k << "):-------\n";
     // 	cout << "V_0.at("  << k << "):\n" << V_0.at(k)  << endl << endl;
     // 	cout << "V_x.at("  << k << "):\n" << V_x.at(k)  << endl << endl;
     // 	cout << "V_xx.at(" << k << "):\n" << V_xx.at(k) << endl << endl;
     // 	cout << endl;
     // }
}

void Lqr::takeBackwardStep(int timestep)
{
     /******************* evaluate cost function ******************/
     // set inputs
     costFunctions.at(timestep).setInput( xTrajectory.at( timestep ), IDX_COST_INPUTS_X_K );
     costFunctions.at(timestep).setInput( uTrajectory.at( timestep ), IDX_COST_INPUTS_U_K );

     // evaluate
     costFunctions.at(timestep).evaluate();

     // get outputs
     costFunctions.at(timestep).getOutput(  cost_0.at(timestep), IDX_COST_OUTPUTS_COST_0_K  );
     costFunctions.at(timestep).getOutput(  cost_x.at(timestep), IDX_COST_OUTPUTS_COST_X_K  );
     costFunctions.at(timestep).getOutput(  cost_u.at(timestep), IDX_COST_OUTPUTS_COST_U_K  );
     costFunctions.at(timestep).getOutput( cost_xx.at(timestep), IDX_COST_OUTPUTS_COST_XX_K );
     costFunctions.at(timestep).getOutput( cost_xu.at(timestep), IDX_COST_OUTPUTS_COST_XU_K );
     costFunctions.at(timestep).getOutput( cost_uu.at(timestep), IDX_COST_OUTPUTS_COST_UU_K );

     cost_xx.at(timestep) += stateRegularization;
     cost_uu.at(timestep) += actionRegularization;

#ifdef EVAL_DYNAMICS_FUNCTION
     {
	  /******************* evaluate dynamics functions (for debugging/plotting only) ******************/
	  // set inputs
	  dynamicsFcn.setInput( xTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_X_K       );
	  dynamicsFcn.setInput( uTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_U_K       );
	  dynamicsFcn.setInput(      cost_0.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_0_K  );
	  dynamicsFcn.setInput(      cost_x.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_X_K  );
	  dynamicsFcn.setInput(      cost_u.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_U_K  );
	  dynamicsFcn.setInput(     cost_xx.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K );
	  dynamicsFcn.setInput(     cost_xu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K );
	  dynamicsFcn.setInput(     cost_uu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K );
	  dynamicsFcn.setInput(         V_0.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1   );
	  dynamicsFcn.setInput(         V_x.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1   );
	  dynamicsFcn.setInput(        V_xx.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1  );

	  // evaluate
	  dynamicsFcn.evaluate();

	  // get outputs
	  dynamicsFcn.getOutput( f0Trajectory.at(timestep), 0 );
	  dynamicsFcn.getOutput( functionTrajectory.at(timestep), 1 );
	  dynamicsFcn.getOutput( fuTrajectory.at(timestep), 2 );
     }
#endif
 
#ifdef EVAL_Q_FUNCTION
     {
	  /******************* evaluate Q function (for debugging/plotting only) ******************/
	  // set inputs
	  qFcn.setInput( xTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_X_K       );
	  qFcn.setInput( uTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_U_K       );
	  qFcn.setInput(      cost_0.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_0_K  );
	  qFcn.setInput(      cost_x.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_X_K  );
	  qFcn.setInput(      cost_u.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_U_K  );
	  qFcn.setInput(     cost_xx.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K );
	  qFcn.setInput(     cost_xu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K );
	  qFcn.setInput(     cost_uu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K );
	  qFcn.setInput(         V_0.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1   );
	  qFcn.setInput(         V_x.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1   );
	  qFcn.setInput(        V_xx.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1  );

	  // evaluate
	  qFcn.evaluate();

	  // get outputs
	  qFcn.getOutput(  Q0Trajectory.at(timestep), IDX_COST_OUTPUTS_COST_0_K  );
	  qFcn.getOutput(  QxTrajectory.at(timestep), IDX_COST_OUTPUTS_COST_X_K  );
	  qFcn.getOutput(  QuTrajectory.at(timestep), IDX_COST_OUTPUTS_COST_U_K  );
	  qFcn.getOutput( QxxTrajectory.at(timestep), IDX_COST_OUTPUTS_COST_XX_K );
	  qFcn.getOutput( QxuTrajectory.at(timestep), IDX_COST_OUTPUTS_COST_XU_K );
	  qFcn.getOutput( QuuTrajectory.at(timestep), IDX_COST_OUTPUTS_COST_UU_K );
     }
#endif


     /*************** evaluate backward sweep function ************/
     // set inputs
     backwardSweepFcn.setInput( xTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_X_K       );
     backwardSweepFcn.setInput( uTrajectory.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_U_K       );
     backwardSweepFcn.setInput(      cost_0.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_0_K  );
     backwardSweepFcn.setInput(      cost_x.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_X_K  );
     backwardSweepFcn.setInput(      cost_u.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_U_K  );
     backwardSweepFcn.setInput(     cost_xx.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XX_K );
     backwardSweepFcn.setInput(     cost_xu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_XU_K );
     backwardSweepFcn.setInput(     cost_uu.at( timestep     ), IDX_BACKWARD_SWEEP_INPUTS_COST_UU_K );
     backwardSweepFcn.setInput(         V_0.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_0_KP1   );
     backwardSweepFcn.setInput(         V_x.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_X_KP1   );
     backwardSweepFcn.setInput(        V_xx.at( timestep + 1 ), IDX_BACKWARD_SWEEP_INPUTS_V_XX_KP1  );

     // evaluate
     backwardSweepFcn.evaluate();

     // get outputs
     backwardSweepFcn.getOutput(    uOpenLoop.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_U_OPEN_LOOP_K   );
     backwardSweepFcn.getOutput( feedbackGain.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_FEEDBACK_GAIN_K );
     backwardSweepFcn.getOutput(          V_0.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_0_K           );
     backwardSweepFcn.getOutput(          V_x.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_X_K           );
     backwardSweepFcn.getOutput(         V_xx.at(timestep), IDX_BACKWARD_SWEEP_OUTPUTS_V_XX_K          );

     // nominal trajectory (for feedback purposes)
     xNominalTrajectory.at( timestep ) = xTrajectory.at( timestep );
}


void Lqr::runForwardSweep()
{
     for (int k = 0; k < N-1; k++)
	  takeForwardStep(k);
}

void Lqr::takeForwardStep(int timestep)
{
     xNominalTrajectory.at( timestep + 1 ) = xTrajectory.at( timestep + 1 );
	
     // open loop
     uTrajectory.at(timestep) = uOpenLoop.at(timestep);
     // add feedback
     uTrajectory.at(timestep) += mul( feedbackGain.at(timestep), xTrajectory.at(timestep) - xNominalTrajectory.at(timestep) );

     // bound actions
#define BOUND( var, lb, ub ){			\
	  if ( var > ub )			\
	       var = ub;			\
	  if ( var < lb )			\
	       var = lb;			\
     }
     for (int k=0; k<ode.nu(); k++)
	  BOUND( uTrajectory.at(timestep).at(k), lbAction.at(timestep).at(k), ubAction.at(timestep).at(k) );

     double dt = (tf - t0)/(N - 1.0);
     xTrajectory.at(timestep + 1) = ode.rk4Step( xTrajectory.at(timestep),
						 uTrajectory.at(timestep),
						 t0 + timestep*dt,
						 dt);

     // xTrajectory.at(timestep + 1) = ode.eulerStep( xTrajectory.at(timestep),
     // 											  uTrajectory.at(timestep),
     // 											  dummyParams,
     // 											  t0 + timestep*dt,
     // 											  dt);
}


void Lqr::boundAction( vector<double> lb_, vector<double> ub_ )
{
     for (int k=0; k<N; k++)
	  boundAction( lb_, ub_, k );

}

void Lqr::boundAction( vector<double> lb_, vector<double> ub_, int timestep )
{
     lbAction.at(timestep) = lb_;
     ubAction.at(timestep) = ub_;
}
