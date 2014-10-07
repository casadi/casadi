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

// ddpSpring.cpp

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
	// constants
	double mass = 1;
	double k    = 10;
	double b    = 0.2;

	SX x = state["x"];
	SX v = state["v"];
	SX u = action["u"];

	xDot["x"] = v;
	xDot["v"] = ( - k*x - b*v + u )/mass;
}

Ode
getOde()
{
	Ode ode("spring");

	ode.addState("x");
	ode.addState("v");

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
	SX cost;
	if (timestep == N-1){
		cost = SQR(state["x"]) + SQR(state["v"]);
		return cost;
	}

	cost = SQR(state["x"]) + SQR(state["v"]) + 0.2*SQR(action["u"]);
	
	return cost;
}


int
main()
{
	Ode ode = getOde();
	Ocp ocp;
	double tEnd = 6;
  
	MultipleShooting & ms = ocp.addMultipleShooting( "spring", ode, 0.0, tEnd, 100 );

	// cost function
	// SX x       = ms.getState(      "x", ms.N-1 );
	// SX thetaf  = ms.getState(      "v", ms.N-1 );
	// SX vthetaf = ms.getState( "vtheta", ms.N-1 );

	ocp.objFun = 0;
	for (int k=0; k < ms.N; k++){
		SX xk = ms.getState(  "x", k );
		SX vk = ms.getState(  "v", k );
		SX uk = ms.getAction( "u", k );

		ocp.objFun += SQR( xk ) + SQR( vk );
		if ( k < ms.N-1 )
			ocp.objFun += 0.2*SQR( uk );
	}


	// bounds
	ms.boundStateAction( "v", -22, 22 );

	/// initial conditions
	ms.boundStateAction( "x", 1.0, 1.0, 0 );
	ms.boundStateAction( "v", 0.0, 0.0, 0 );

	// // final conditions
	// ms.boundStateAction( "x", 0.0, 0.0, ms.N-1 );
	// ms.boundStateAction( "v", 0.0, 0.0, ms.N-1 );
	// ms.boundStateAction( "u", 0.0, 0.0, ms.N-1 );

	// // initial guess
	// for (int k=0; k<ms.N; k++){
	// 	double zero_to_one = k/(ms.N - 1.0);
	// 	ms.setStateActionGuess( "u", sin(2*M_PI*zero_to_one), k);
	// }

	// solve
	SnoptInterface si(ocp);
	
	si.run();
	ocp.writeOctaveOutput("spring_multiple_shooting_out");
	
	// Print the optimal cost
	cout << "optimal objFcn: " << si.F[0] << endl;
	
	
	/************** run ddp ****************/
	cout << "\nrunning DDP\n";
	double t0 = 0;
 	Ddp ddp(ode, t0, tEnd, ms.N, &cost);

	ddp.xTrajectory.at(0) = ocp.getStateSolution(0);

	// regularization
	// ddp.stateRegularization[0,0]  = 0.001;
	// ddp.stateRegularization[1,1]  = 0.001;
	// ddp.actionRegularization[0,0] = 0.001;

	// for (int k=0; k<ms.N; k++)
	// 	ddp.xTrajectory.at(k) = ocp.getStateSolution(k);
	// for (int k=0; k<ms.N-1; k++)
	// 	ddp.uTrajectory.at(k) = ocp.getActionSolution(k);

	for (int k=0; k<2; k++){
		cout << "ddp iter: " << k << endl;
		ddp.runBackwardSweep();
		ddp.runForwardSweep();
	}

	ocp.setStates(  ddp.xTrajectory );
	ocp.setActions( ddp.uTrajectory );
	
	ocp.writeOctaveOutput("spring_ddp_out");
	
	cout << "successful finish\n";
	return 0;
}
