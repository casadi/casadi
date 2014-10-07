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

// Ddp.cpp
// Greg Horn
// Casadi 2011

#include "Ddp.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace casadi;
using namespace std;

Ddp::Ddp(Ode & _ode, double t0_, double tf_, int N_, SX (*cost_)(map<string,SX> state, map<string,SX> action, int timestep, int _N)) : ode(_ode)
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

	  Q0Trajectory.push_back(  DMatrix(        1,        1, 0.0) );
	  QxTrajectory.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  QuTrajectory.push_back(  DMatrix( ode.nu(),        1, 0.0) );
	  QxxTrajectory.push_back( DMatrix( ode.nx(), ode.nx(), 0.0) );
	  QxuTrajectory.push_back( DMatrix( ode.nx(), ode.nu(), 0.0) );
	  QuuTrajectory.push_back( DMatrix( ode.nu(), ode.nu(), 0.0) );

	  xTrajectory.push_back(         DMatrix( ode.nx(),        1, 0.0) );
	  xNominalTrajectory.push_back(  DMatrix( ode.nx(),        1, 0.0) );
	  uTrajectory.push_back(         DMatrix( ode.nu(),        1, 0.0) );
	  uOpenLoop.push_back(           DMatrix( ode.nu(),        1, 0.0) );
	  feedbackGain.push_back(        DMatrix( ode.nu(), ode.nx(), 0.0) );

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
     setupQFunctions();
}

Ddp::~Ddp() {}


void Ddp::setupQFunctions()
{
     /*************** inputs **************/
     SX xk = ssym("xk", ode.nx());
     SX uk = ssym("uk", ode.nu());

     SX V_0_kp1(  ssym(  "V_0_kp1",        1,        1) );
     SX V_x_kp1(  ssym(  "V_x_kp1", ode.nx(),        1) );
     SX V_xx_kp1( ssym( "V_xx_kp1", ode.nx(), ode.nx()) );

     vector<SX> qInputs(NUM_Q_INPUTS);
     qInputs.at(IDX_Q_INPUTS_X_K)      = xk;
     qInputs.at(IDX_Q_INPUTS_U_K)      = uk;
     qInputs.at(IDX_Q_INPUTS_V_0_KP1)  = V_0_kp1;
     qInputs.at(IDX_Q_INPUTS_V_X_KP1)  = V_x_kp1;
     qInputs.at(IDX_Q_INPUTS_V_XX_KP1) = V_xx_kp1;

     SX dx = ssym("dx", ode.nx());
     SX du = ssym("du", ode.nu());

     /**************** dynamics *********************/
     // dummy params for now
     map<string,SX> dummyParams;
     double dt = (tf - t0)/(N - 1);
     SX f = ode.rk4Step( xk + dx, uk + du, uk + du, dummyParams, t0, dt); // timestep dependent: f(x,u,__t__) - same as cost
     SX f0 = ode.rk4Step( xk, uk, uk, dummyParams, t0, dt); // timestep dependent: f(x,u,__t__) - same as cost
     // SX f = ode.eulerStep( xk + dx, uk + du, dummyParams, SX(t0), SX(dt)); // timestep dependent:  f(x,u,__t__) - same as cost
     // SX f0 = ode.eulerStep( xk, uk, dummyParams, SX(t0), SX(dt)); // timestep dependent:  f(x,u,__t__) - same as cost

     /**************** Q function *********************/
     SX xk_p_dx( xk + dx );
     SX uk_p_du( uk + du );

     map<string,SX> xkMap = ode.getStateMap(  xk_p_dx );
     map<string,SX> ukMap = ode.getActionMap( uk_p_du );

     // map<string,SX> xkMap = ode.getStateMap(  xk );
     // map<string,SX> ukMap = ode.getActionMap( uk );

     SX df = f - f0;

     SX dxZeros( zerosSX( ode.nx(), 1 ) );
     SX duZeros( zerosSX( ode.nu(), 1 ) );

     makeDense( dxZeros );
     makeDense( duZeros );

     for (int k=0; k<N; k++){

	  // function
	  SX Q_0_k(costFcnExt( xkMap, ukMap, k, N ) + V_0_kp1 + mul( V_x_kp1.trans(), df )
			 + mul( df.trans(), mul( V_xx_kp1, df ) )/2 );

	  // jacobian
	  SX Q_x_k =  gradient( Q_0_k, dx );
	  SX Q_u_k =  gradient( Q_0_k, du );

	  // hessian
	  SX Q_xx_k = jacobian( Q_x_k, dx );
	  SX Q_xu_k = jacobian( Q_x_k, du ); // == jacobian( Q_u, x ).trans()
	  SX Q_uu_k = jacobian( Q_u_k, du );

	  Q_0_k  = substitute(  Q_0_k, dx, dxZeros );
	  Q_x_k  = substitute(  Q_x_k, dx, dxZeros );
	  Q_u_k  = substitute(  Q_u_k, dx, dxZeros );
	  Q_xx_k = substitute( Q_xx_k, dx, dxZeros );
	  Q_xu_k = substitute( Q_xu_k, dx, dxZeros );
	  Q_uu_k = substitute( Q_uu_k, dx, dxZeros );

	  Q_0_k  = substitute(  Q_0_k, du, duZeros );
	  Q_x_k  = substitute(  Q_x_k, du, duZeros );
	  Q_u_k  = substitute(  Q_u_k, du, duZeros );
	  Q_xx_k = substitute( Q_xx_k, du, duZeros );
	  Q_xu_k = substitute( Q_xu_k, du, duZeros );
	  Q_uu_k = substitute( Q_uu_k, du, duZeros );

	  simplify(Q_0_k);
	  simplify(Q_x_k);
	  simplify(Q_u_k);
	  simplify(Q_xx_k);
	  simplify(Q_xu_k);
	  simplify(Q_uu_k);

	  // workaround bug where size() != size1()*size2()
	  makeDense(Q_0_k);
	  makeDense(Q_x_k);
	  makeDense(Q_u_k);
	  makeDense(Q_xx_k);
	  makeDense(Q_xu_k);
	  makeDense(Q_uu_k);

	  // outputs
	  vector<SX> qOutputs(NUM_Q_OUTPUTS);
	  qOutputs.at(IDX_Q_OUTPUTS_Q_0_K)  = Q_0_k;
	  qOutputs.at(IDX_Q_OUTPUTS_Q_X_K)  = Q_x_k;
	  qOutputs.at(IDX_Q_OUTPUTS_Q_U_K)  = Q_u_k;
	  qOutputs.at(IDX_Q_OUTPUTS_Q_XX_K) = Q_xx_k;
	  qOutputs.at(IDX_Q_OUTPUTS_Q_XU_K) = Q_xu_k;
	  qOutputs.at(IDX_Q_OUTPUTS_Q_UU_K) = Q_uu_k;
		
	  // sx function
	  qFunctions.push_back( SXFunction( qInputs, qOutputs ) );
	  qFunctions.at(k).init();
     }
}

void Ddp::setupBackwardSweepFunction()
{
     /************ inputs ************/
     SX uk(   ssym(   "uk", ode.nu(),        1 ) );

     SX Q_0(  ssym(  "Q_0",        1,        1 ) );
     SX Q_x(  ssym(  "Q_x", ode.nx(),        1 ) );
     SX Q_u(  ssym(  "Q_u", ode.nu(),        1 ) );
     SX Q_xx( ssym( "Q_xx", ode.nx(), ode.nx() ) );
     SX Q_xu( ssym( "Q_xu", ode.nx(), ode.nu() ) );
     SX Q_uu( ssym( "Q_uu", ode.nu(), ode.nu() ) );


     vector<SX> backwardSweepInputs(NUM_BACKWARD_SWEEP_INPUTS);
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_U_K)    = uk;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_0_K)  = Q_0;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_X_K)  = Q_x;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_U_K)  = Q_u;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_XX_K) = Q_xx;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_XU_K) = Q_xu;
     backwardSweepInputs.at(IDX_BACKWARD_SWEEP_INPUTS_Q_UU_K) = Q_uu;

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
}

void Ddp::runBackwardSweep()
{
     xNominalTrajectory.at(N-1) = xTrajectory.at(N-1);

     /******************* V(N-1) == Cost(N-1) ******************/
     // set inputs
     DMatrix V_0_zero(         1,        1, 0.0 );
     DMatrix V_x_zero(  ode.nx(),        1, 0.0 );
     DMatrix V_xx_zero( ode.nx(), ode.nx(), 0.0 );

     qFunctions.at( N-1 ).setInput( xTrajectory.at(N-1), IDX_Q_INPUTS_X_K      );
     qFunctions.at( N-1 ).setInput( uTrajectory.at(N-1), IDX_Q_INPUTS_U_K      );
     qFunctions.at( N-1 ).setInput(            V_0_zero, IDX_Q_INPUTS_V_0_KP1  );
     qFunctions.at( N-1 ).setInput(            V_x_zero, IDX_Q_INPUTS_V_X_KP1  );
     qFunctions.at( N-1 ).setInput(           V_xx_zero, IDX_Q_INPUTS_V_XX_KP1 );

     // evaluate
     qFunctions.at( N-1 ).evaluate();

     // get outputs
     qFunctions.at( N-1 ).getOutput(  V_0.at( N-1 ), IDX_Q_OUTPUTS_Q_0_K  );
     qFunctions.at( N-1 ).getOutput(  V_x.at( N-1 ), IDX_Q_OUTPUTS_Q_X_K  );
     qFunctions.at( N-1 ).getOutput( V_xx.at( N-1 ), IDX_Q_OUTPUTS_Q_XX_K );

     /*********** run backward sweep ************/
     for (int k = N-2; k >= 0; k--)
	  takeBackwardStep(k);
}

void Ddp::takeBackwardStep(int timestep)
{
     /******************* evaluate Q function ******************/
     // set inputs
     qFunctions.at(timestep).setInput( xTrajectory.at( timestep     ), IDX_Q_INPUTS_X_K      );
     qFunctions.at(timestep).setInput( uTrajectory.at( timestep     ), IDX_Q_INPUTS_U_K      );
     qFunctions.at(timestep).setInput(         V_0.at( timestep + 1 ), IDX_Q_INPUTS_V_0_KP1  );
     qFunctions.at(timestep).setInput(         V_x.at( timestep + 1 ), IDX_Q_INPUTS_V_X_KP1  );
     qFunctions.at(timestep).setInput(        V_xx.at( timestep + 1 ), IDX_Q_INPUTS_V_XX_KP1 );

     // evaluate
     qFunctions.at(timestep).evaluate();

     // get outputs
     qFunctions.at(timestep).getOutput(  Q0Trajectory.at(timestep), IDX_Q_OUTPUTS_Q_0_K  );
     qFunctions.at(timestep).getOutput(  QxTrajectory.at(timestep), IDX_Q_OUTPUTS_Q_X_K  );
     qFunctions.at(timestep).getOutput(  QuTrajectory.at(timestep), IDX_Q_OUTPUTS_Q_U_K  );
     qFunctions.at(timestep).getOutput( QxxTrajectory.at(timestep), IDX_Q_OUTPUTS_Q_XX_K );
     qFunctions.at(timestep).getOutput( QxuTrajectory.at(timestep), IDX_Q_OUTPUTS_Q_XU_K );
     qFunctions.at(timestep).getOutput( QuuTrajectory.at(timestep), IDX_Q_OUTPUTS_Q_UU_K );

     QxxTrajectory.at(timestep) += stateRegularization;
     QuuTrajectory.at(timestep) += actionRegularization;

     /*************** evaluate backward sweep function ************/
     // set inputs
     backwardSweepFcn.setInput(   uTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_U_K    );
     backwardSweepFcn.setInput(  Q0Trajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_0_K  );
     backwardSweepFcn.setInput(  QxTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_X_K  );
     backwardSweepFcn.setInput(  QuTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_U_K  );
     backwardSweepFcn.setInput( QxxTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_XX_K );
     backwardSweepFcn.setInput( QxuTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_XU_K );
     backwardSweepFcn.setInput( QuuTrajectory.at( timestep ), IDX_BACKWARD_SWEEP_INPUTS_Q_UU_K );

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


void Ddp::runForwardSweep()
{
     for (int k = 0; k < N-1; k++)
	  takeForwardStep(k);
}

void Ddp::takeForwardStep(int timestep)
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


void Ddp::boundAction( vector<double> lb_, vector<double> ub_ )
{
     for (int k=0; k<N; k++)
	  boundAction( lb_, ub_, k );

}

void Ddp::boundAction( vector<double> lb_, vector<double> ub_, int timestep )
{
     lbAction.at(timestep) = lb_;
     ubAction.at(timestep) = ub_;
}
