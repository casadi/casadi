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

// ddpDoubleCartpole.cpp

#include <iostream>
#include <cstdlib>

#include <core/std_vector_tools.hpp>
#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>

#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"
#include "Ddp.hpp"

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace casadi;
using namespace std;

#define SQR(sqr_me) ((sqr_me)*(sqr_me))

void
dxdt(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
	double g  = 9.8;
	double L1 = 0.1;
	double L2 = 0.1;
	double m0 = 0.1;
	double mp = 0.03;
	double m1 = mp;
	double m2 = mp;
     
	double d1 = m0 + m1 + m2;
	double d2 = (0.5*m1 + m2)*L1;
	double d3 = 0.5 * m2 * L2;
	double d4 = ( m1/3 + m2 )*SQR(L1);
	double d5 = 0.5 * m2 * L1 * L2;
	double d6 = m2 * SQR(L2)/3;
	double f1 = (0.5*m1 + m2) * L1 * g;
	double f2 = 0.5 * m2 * L2 * g;

	// th0  = y(1);
	// th1  = y(2);
	// th2  = y(3);
	// th0d = y(4);
	// th1d = y(5);
	// th2d = y(6);

	SX th0  = state["th0"];
	SX th1  = state["th1"];
	SX th2  = state["th2"];
	SX th0d = state["th0d"];
	SX th1d = state["th1d"];
	SX th2d = state["th2d"];

	// D = [          d1,     d2*cos(th1),     d3*cos(th2);
	// 	     d2*cos(th1),              d4, d5*cos(th1-th2);
	// 	     d3*cos(th2), d5*cos(th1-th2),              d6;];

	SX D( zerosSX(3,3) );
	makeDense(D);
	D[0,0] =          d1;   D[0,1] =     d2*cos(th1);   D[0,2] =     d3*cos(th2);
	D[1,0] = d2*cos(th1);   D[1,1] =              d4;   D[1,2] = d5*cos(th1-th2);
	D[2,0] = d3*cos(th2);   D[2,1] = d5*cos(th1-th2);   D[2,2] =              d6;

	// C = [0,     -d2*sin(th1)*th1d,    -d3*sin(th2)*th2d;
	// 	    0,                     0, d5*sin(th1-th2)*th2d;
	// 	    0, -d5*sin(th1-th2)*th1d,                    0;];
	SX C( zerosSX(3,3) );
	makeDense(C);
	C[0,0] = 0;   C[0,1] =      -d2*sin(th1)*th1d;   C[0,2] =     -d3*sin(th2)*th2d;
	C[1,0] = 0;   C[1,1] =                      0;   C[1,2] =  d5*sin(th1-th2)*th2d;
	C[2,0] = 0;   C[2,1] =  -d5*sin(th1-th2)*th1d;   C[2,2] =                     0;

	// G = [0; -f1*sin(th1); -f2*sin(th2);];
	SX G( zerosSX(3,1) );
	makeDense(G);
	G.at(0) = 0;
	G.at(1) = -f1*sin(th1);
	G.at(2) = -f2*sin(th2);

	// H = [1;0;0;];
	SX H( zerosSX(3,1) );
	makeDense(H);
	H.at(0) = 1;
	H.at(1) = 0;
	H.at(2) = 0;

	// dy(1:3) = y(4:6);
	xDot["th0"] = th0d;
	xDot["th1"] = th1d;
	xDot["th2"] = th2d;

	// dy(4:6) = D\( - C*y(4:6) - G + H*u );
	SX vel( zerosSX(3,1) );
	makeDense(vel);
	vel.at(0) = th0d;
	vel.at(1) = th1d;
	vel.at(2) = th2d;
	SX accel = mul( inv(D), - mul( C, vel ) - G + mul( H, SX(action["u"]) ) );

	simplify(accel.at(0));
	simplify(accel.at(1));
	simplify(accel.at(2));

	xDot["th0d"] = accel.at(0);
	xDot["th1d"] = accel.at(1);
	xDot["th2d"] = accel.at(2);

	// cout << th0 << endl;
	// cout << th1 << endl;
	// cout << th2 << endl;

	// cout << th0d << endl;
	// cout << th1d << endl;
	// cout << th2d << endl;

	// cout << accel.at(0) << endl;
	// cout << accel.at(1) << endl;
	// cout << accel.at(2) << endl;

	outputs["cart_x"] = th0;
	outputs["cart_y"] = 0;
	outputs["bob0_x"] = th0 + L1*sin(th1);
	outputs["bob0_y"] = - L1*cos(th1);
	outputs["bob1_x"] = th0 + L1*sin(th1) + L2*sin(th2);
	outputs["bob1_y"] = - L1*cos(th1) - L2*cos(th2);
}


Ode
getOde()
{
	Ode ode("double_cartpole");
	ode.addState("th0");
	ode.addState("th1");
	ode.addState("th2");
	ode.addState("th0d");
	ode.addState("th1d");
	ode.addState("th2d");

	ode.addOutput( "cart_x" );
	ode.addOutput( "cart_y" );
	ode.addOutput( "bob0_x" );
	ode.addOutput( "bob0_y" );
	ode.addOutput( "bob1_x" );
	ode.addOutput( "bob1_y" );

	ode.addAction("u");

	ode.dxdt = &dxdt;

	return ode;
}

SX
cost(map<string,SX> state,
	 map<string,SX> action,
	 int timestep,
	 int N)
{
	SX costRet;
	if (timestep == N-1){
		costRet  = 1e4*SQR( state["th0"] );
		costRet += 100*SQR( state["th1"] );
		costRet += 100*SQR( state["th2"] );
		costRet += 100*SQR( state["dth0"] );
		costRet +=     SQR( state["th1d"] );
		costRet +=     SQR( state["th2d"] );
		return costRet;
	}

	costRet = 25*(SQR( state["th0"] ) + SQR( state["th1"] ) + SQR( state["th2"] ) + SQR( state["dth0"] ));
	costRet += 0.01*( SQR( state["th1d"] )  +SQR( state["th2d"] ));
	costRet += 0.1*SQR(action["u"]);

	// // add barrier function
	// double uUb =  10.1;
	// double uLb = -10.1;
	// double mu = 1.0;
	// SX uBarrierUb = -mu*log(  uUb - action["u"] );
	// SX uBarrierLb = -mu*log( -uLb + action["u"] );
	// costRet += uBarrierUb + uBarrierLb;

	return costRet;
}


int
main()
{
	/************* get solution **************/
	double trackLength = 4;
  
	Ode ode = getOde();
	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");
  
	MultipleShooting & ms = ocp.addMultipleShooting("double_cartpole", ode, 0.0, tEnd, 60);

	// cost function
	// SX xf      = ms.getState(      "x", ms.N-1);
	// SX thetaf  = ms.getState(  "theta", ms.N-1);
	// SX vthetaf = ms.getState( "vtheta", ms.N-1);
	ocp.objFun = tEnd;
	for (int k=0; k<ms.N-1; k++){
		SX u_k   = ms.getAction( "u", k );
		SX u_kp1 = ms.getAction( "u", k+1 );
		// ocp.objFun += 3*SQR( u_k - u_kp1 );
		ocp.objFun += 0.01*SQR( u_k );
	}

	// bounds
	ocp.boundParam("tEnd", 3, 20);
	ocp.setParamGuess("tEnd", 6);
	
	ms.boundStateAction(     "th0", -trackLength/2, trackLength/2);
	ms.boundStateAction(     "th1", -4*M_PI, 4*M_PI);
	ms.boundStateAction(     "th2", -4*M_PI, 4*M_PI);
	ms.boundStateAction(    "th0d", -20, 20);
	ms.boundStateAction(    "th1d", -200, 200);
	ms.boundStateAction(    "th2d", -200, 200);

	ms.boundStateAction("u", -100, 100 );

	// initial guess
	for (int k=0; k<ms.N; k++){
		double zero_to_one = k/(ms.N - 1.0);
		ms.setStateActionGuess( "th0", 0.0, k);
		ms.setStateActionGuess( "th1", (1 - zero_to_one)*M_PI, k);
		ms.setStateActionGuess( "th2", (1 - zero_to_one)*M_PI, k);

		ms.setStateActionGuess( "th0d", 0.0, k);
		ms.setStateActionGuess( "th1d", -M_PI/6, k);
		ms.setStateActionGuess( "th2d", -M_PI/6, k);

		//ms.setStateActionGuess( "u", sin(2*M_PI*zero_to_one), k);
	}

	/// initial conditions
	ms.boundStateAction(  "th0", 0.0, 0.0, 0);
	ms.boundStateAction(  "th1", M_PI-0.1, M_PI-0.1, 0);
	ms.boundStateAction(  "th2", M_PI-0.1, M_PI-0.1, 0);
	ms.boundStateAction( "th0d", 0.0, 0.0, 0);
	ms.boundStateAction( "th1d", 0.0, 0.0, 0);
	ms.boundStateAction( "th2d", 0.0, 0.0, 0);
	ms.boundStateAction(    "u", 0.0, 0.0, 0);

	// final conditions
	ms.boundStateAction(  "th0", 0.0, 0.0, ms.N - 1);
	ms.boundStateAction(  "th1", 0.0, 0.0, ms.N - 1);
	ms.boundStateAction(  "th2", 0.0, 0.0, ms.N - 1);
	ms.boundStateAction( "th0d", 0.0, 0.0, ms.N - 1);
	ms.boundStateAction( "th1d", 0.0, 0.0, ms.N - 1);
	ms.boundStateAction( "th2d", 0.0, 0.0, ms.N - 1);

	// solve
	SnoptInterface si(ocp);
	
	si.run();
	ocp.writeOctaveOutput("double_cartpole_multiple_shooting_out");
	
	// Print the optimal cost
	cout << "optimal objFcn: " << si.F[0] << endl;
	
	
	/************** run ddp ****************/
	cout << "\nrunning DDP\n";
	double t0 = 0;
	double tf = ocp.getParamSolution("tEnd");
 	Ddp ddp(ode, t0, tf, ms.N, &cost);

	// regularization
	ddp.stateRegularization[0,0]  = 1;
	ddp.stateRegularization[1,1]  = 1;
	ddp.stateRegularization[2,2]  = 1;
	ddp.stateRegularization[3,3]  = 1;
	ddp.actionRegularization[0,0] = 1;

	// action bounding
	vector<double> ub(1);
	vector<double> lb(1);
	ub.at(0) = 10;
	lb.at(0) = -10;
	ddp.boundAction( lb, ub );

	// load ms solution into ddp
	for (int k=0; k<ms.N; k++)
		ddp.xTrajectory.at(k) = ocp.getStateSolution(k);
	for (int k=0; k<ms.N-1; k++)
		ddp.uTrajectory.at(k) = ocp.getActionSolution(k);
	// for (int k=0; k<ms.N-1; k++)
	// 	ddp.uTrajectory.at(k).at(0) += 2.2;

	// run ddp
	for (int k=0; k<0; k++){
		ddp.runBackwardSweep();
		ddp.runForwardSweep();
		cout << "ddp iter: " << k << ", value: " << ddp.V_0.at(0) << endl;
	}

	ocp.setStates(  ddp.xTrajectory );
	ocp.setActions( ddp.uTrajectory );

	ocp.writeOctaveOutput("double_cartpole_ddp_out");
	
	cout << "successful finish\n";
	return 0;
}
