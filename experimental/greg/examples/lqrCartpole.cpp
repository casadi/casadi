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

// lqrCartpole.cpp

#include <iostream>
#include <cstdlib>

#include <core/std_vector_tools.hpp>
#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>

#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"
#include "Lqr.hpp"

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace casadi;
using namespace std;

#define SQR(sqr_me) ((sqr_me)*(sqr_me))

void
dxdt(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
	// constants
	double g = 9.8;    // acceleration due to gravity
	double l = 2.2;
	double mc = 2;
	double mp = 1;

	SX x = state["x"];
	SX theta = state["theta"];
	SX vx = state["vx"];
	SX vtheta = state["vtheta"];
	SX u = action["u"];

	SX ax = 1/(mc+mp*sin(theta)*sin(theta))*(u+mp*sin(theta)*(l*vtheta*vtheta+g*cos(theta)));
	SX atheta = 1/(l*(mc+mp*sin(theta)*sin(theta)))*(-u*cos(theta) - mp*l*vtheta*vtheta*cos(theta)*sin(theta) - (mc+mp)*g*sin(theta));

	xDot["x"] = vx;
	xDot["theta"] = vtheta;
	xDot["vx"] = ax;
	xDot["vtheta"] = atheta;

	outputs["cart_x"] = x;
	outputs["cart_y"] = 0;
	outputs["bob_x"]  = x + l*sin(theta);
	outputs["bob_y"]  = - l*cos(theta);
}

Ode
getOde()
{
	Ode ode("cartpole");
	ode.addState("x");
	ode.addState("theta");
	ode.addState("vx");
	ode.addState("vtheta");
	// ode.addState("u");

	ode.addOutput( "cart_x" );
	ode.addOutput( "cart_y" );
	ode.addOutput( "bob_x" );
	ode.addOutput( "bob_y" );

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
		costRet = N*(10*SQR(state["x"]) + 10*SQR(state["theta"] - M_PI) + SQR(state["vx"]) + SQR(state["vtheta"]));
		return costRet;
	}

	//	costRet = 2*SQR(state["x"]) + 3*SQR(state["theta"]) + 4*SQR(state["vx"]) + 5*SQR(state["vtheta"]) + 6*SQR(action["u"]);
	costRet = 0.01*SQR(action["u"]);

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
  
	MultipleShooting & ms = ocp.addMultipleShooting("cartpole", ode, 0.0, tEnd, 60);

	// cost function
	SX xf      = ms.getState(      "x", ms.N-1);
	SX thetaf  = ms.getState(  "theta", ms.N-1);
	SX vthetaf = ms.getState( "vtheta", ms.N-1);

	//	ocp.objFun = tEnd  + 50*cos(thetaf) + 5*vthetaf*vthetaf;
	ocp.objFun = ms.N*cos(thetaf) + 5*vthetaf*vthetaf;
	for (int k=0; k<ms.N-1; k++){
		SX u_k   = ms.getAction( "u", k );
		SX u_kp1 = ms.getAction( "u", k+1 );
		// ocp.objFun += 3*SQR( u_k - u_kp1 );
		ocp.objFun += 0.01*SQR( u_k );
	}

	// bounds
	ocp.boundParam("tEnd", 6, 6);
	
	ms.boundStateAction(      "x", -trackLength/2, trackLength/2);
	ms.boundStateAction(     "vx", -22, 22);
	ms.boundStateAction(  "theta", -50, 50);
	ms.boundStateAction( "vtheta", -50, 50);

	// ms.boundStateAction("u",-20,20);

	/// initial conditions
	ms.boundStateAction(      "x", 0.0, 0.0, 0);
	ms.boundStateAction(  "theta", 0.1, 0.1, 0);
	ms.boundStateAction(     "vx", 0.0, 0.0, 0);
	ms.boundStateAction( "vtheta", 0.0, 0.0, 0);
	ms.boundStateAction(      "u", 0.0, 0.0, 0);

	// final conditions
	ms.boundStateAction(      "x", 0.0, 0.0, ms.N-1);
	ms.boundStateAction(     "vx", 0.0, 0.0, ms.N-1);
	ms.boundStateAction( "vtheta", 0.0, 0.0, ms.N-1);
	ms.boundStateAction(  "theta", 3.14159, 3.14159, ms.N-1);
	ms.boundStateAction(      "u", 0.0, 0.0, ms.N-1);

	// initial guess
	for (int k=0; k<ms.N; k++){
		double zero_to_one = k/(ms.N - 1.0);
		ms.setStateActionGuess( "u", sin(2*M_PI*zero_to_one), k);
	}

	// solve
	SnoptInterface si(ocp);
	
	si.run();
	ocp.writeOctaveOutput("cartpole_multiple_shooting_out");
	
	// Print the optimal cost
	cout << "optimal objFcn: " << si.F[0] << endl;
	
	
	/************** run lqr ****************/
	cout << "\nrunning ILQR\n";
	double t0 = 0;
	double tf = ocp.getParamSolution("tEnd");
 	Lqr lqr(ode, t0, tf, ms.N, &cost);

	// regularization
	lqr.stateRegularization[0,0]  = 0.01;
	lqr.stateRegularization[1,1]  = 0.01;
	lqr.stateRegularization[2,2]  = 0.01;
	lqr.stateRegularization[3,3]  = 0.01;
	lqr.actionRegularization[0,0] = 0.1;

	// action bounding
	vector<double> ub(1);
	vector<double> lb(1);
	ub.at(0) = 20;
	lb.at(0) = -30;
	lqr.boundAction( lb, ub );

	// load ms solution into lqr
	for (int k=0; k<ms.N; k++)
		lqr.xTrajectory.at(k) = ocp.getStateSolution(k);
	for (int k=0; k<ms.N-1; k++)
		lqr.uTrajectory.at(k) = ocp.getActionSolution(k);

	// run lqr
	for (int k=0; k<1000; k++){
		lqr.runBackwardSweep();
		lqr.runForwardSweep();

		// manually calculate value function
		double value = 0;
		for (int j=0; j<lqr.N - 1; j++)
			value += lqr.cost_0.at(j).at(0);
		value += lqr.V_0.at(lqr.N-1).at(0);
		cout << "ilqr iter: " << k << ", value: " << lqr.V_0.at(0) << ", my value: " << value << endl;
	}

	ocp.setStates(  lqr.xTrajectory );
	ocp.setActions( lqr.uTrajectory );

	ocp.writeOctaveOutput("cartpole_lqr_out");
	
	cout << "successful finish\n";
	return 0;
}
