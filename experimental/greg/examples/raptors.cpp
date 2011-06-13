/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include <iostream>
#include <cstdlib>

#include <casadi/stl_vector_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/fx/sx_function.hpp>
//#include <casadi/fx/jacobian.hpp>


#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace CasADi;
using namespace std;

#define EPS 1e-12

void
dxdt(map<string,SX> &xDot, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
    SX xh = state["xh"];
	SX yh = state["yh"];
	SX x1 = state["x1"];
	SX y1 = state["y1"];
	SX x2 = state["x2"];
	SX y2 = state["y2"];
	SX x3 = state["x3"];
	SX y3 = state["y3"];
	SX s1 = state["s1"];
	SX s2 = state["s2"];
	SX s3 = state["s3"];
	SX theta = state["theta"];
 
	#define VHUMAN 6.0

	// distance from human to each raptor
	SX d1 = sqrt( (x1-xh)*(x1-xh) + (y1-yh)*(y1-yh) ) + EPS;
	SX d2 = sqrt( (x2-xh)*(x2-xh) + (y2-yh)*(y2-yh) ) + EPS;
	SX d3 = sqrt( (x3-xh)*(x3-xh) + (y3-yh)*(y3-yh) ) + EPS;

	// velocity of the human
	SX vxh = VHUMAN*cos(theta);
	SX vyh = VHUMAN*sin(theta);
	
	// velocity of each raptor
	SX vx1 = s1 * (xh - x1) / d1;
	SX vy1 = s1 * (yh - y1) / d1;
	
	SX vx2 = s2 * (xh - x2) / d2;
	SX vy2 = s2 * (yh - y2) / d2;
	
	SX vx3 = s3 * (xh - x3) / d3;
	SX vy3 = s3 * (yh - y3) / d3;

	// outputs
	xDot["xh"] = vxh;
	xDot["yh"] = vyh;

	xDot["x1"] = vx1;
	xDot["y1"] = vy1;

	xDot["x2"] = vx2;
	xDot["y2"] = vy2;

	xDot["x3"] = vx3;
	xDot["y3"] = vy3;

	xDot["s1"] = (10 - s1)/1.0;
	xDot["s2"] = (25 - s2)/2.5;
	xDot["s3"] = (25 - s3)/2.5;

	xDot["theta"] = action["dtheta"];
}

Ode
getRaptorOde()
{
	Ode ode("raptor");
	ode.addState("xh");
	ode.addState("yh");

	ode.addState("x1");
	ode.addState("y1");

	ode.addState("y2");
	ode.addState("x2");

	ode.addState("y3");
	ode.addState("x3");

	ode.addState("s1");
	ode.addState("s2");
	ode.addState("s3");

	ode.addState("theta");

	ode.addAction("dtheta");

	ode.dxdt = &dxdt;

	return ode;
}


int
main()
{
	Ode ode = getRaptorOde();
	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");
	MultipleShooting & r0 = ocp.addMultipleShooting("raptor0", ode, 0, tEnd, 60);

	ocp.objFun = -tEnd; // maximum time

	// Bounds
	ocp.boundParam("tEnd", 0.1, 20);
	for (int k=0; k<r0.N; k++){
		r0.boundStateAction( "xh", -100, 100, k );
		r0.boundStateAction( "yh", -100, 100, k );
		r0.boundStateAction( "x1", -100, 100, k );
		r0.boundStateAction( "y1", -100, 100, k );
		r0.boundStateAction( "x2", -100, 100, k );
		r0.boundStateAction( "y2", -100, 100, k );
		r0.boundStateAction( "x3", -100, 100, k );
		r0.boundStateAction( "y3", -100, 100, k );

		r0.boundStateAction( "s1", 0.0, 30.0, k);
		r0.boundStateAction( "s2", 0.0, 30.0, k);
		r0.boundStateAction( "s3", 0.0, 30.0, k);

		r0.boundStateAction( "theta", -10, 10.0, k);
		r0.boundStateAction("dtheta", -30*M_PI/180.0, 30*M_PI/180.0, k);

		SX x1 = r0.getState("x1", k);
		SX y1 = r0.getState("y1", k);
		SX x2 = r0.getState("x2", k);
		SX y2 = r0.getState("y2", k);
		SX x3 = r0.getState("x3", k);
		SX y3 = r0.getState("y3", k);
		SX xh = r0.getState("xh", k);
		SX yh = r0.getState("yh", k);
		SX d1 = sqrt( (x1-xh)*(x1-xh) + (y1-yh)*(y1-yh) ) + EPS;
		SX d2 = sqrt( (x2-xh)*(x2-xh) + (y2-yh)*(y2-yh) ) + EPS;
		SX d3 = sqrt( (x3-xh)*(x3-xh) + (y3-yh)*(y3-yh) ) + EPS;

		ocp.addNonlconIneq( 0.1 - d1 );
		ocp.addNonlconIneq( 0.1 - d2 );
		ocp.addNonlconIneq( 0.1 - d3 );
	}

	//	exit(1);

	// initial conditions
	r0.boundStateAction("xh", 0, 0, 0);
	r0.boundStateAction("yh", 10.0*sqrt(3.0)/3.0, 10.0*sqrt(3.0)/3.0, 0);

	r0.boundStateAction("x1", 0.0, 0.0, 0 );
	r0.boundStateAction("y1", 20.0*sqrt(3.0)/2.0, 20.0*sqrt(3.0)/2.0, 0 );

	r0.boundStateAction("x2", -10.0, -10.0, 0 );
	r0.boundStateAction("y2", 0.0, 0.0, 0 );

	r0.boundStateAction("x3", 10.0, 10.0, 0 );
	r0.boundStateAction("y3", 0.0, 0.0, 0 );

	r0.boundStateAction("s1", 0.0, 0.0, 0 );
	r0.boundStateAction("s2", 0.0, 0.0, 0 );
	r0.boundStateAction("s3", 0.0, 0.0, 0 );


	// initial guesses
	ocp.setParamGuess("tEnd", 2);

	for (int k=0; k<r0.N; k++){
		r0.setStateActionGuess("xh", 0, k);
		r0.setStateActionGuess("yh", 10.0*sqrt(3.0)/3.0, k);

		r0.setStateActionGuess("x1", 0.0, k );
		r0.setStateActionGuess("y1", 20.0*sqrt(3.0)/2.0, k );

		r0.setStateActionGuess("x2", -10.0, k );
		r0.setStateActionGuess("y2", 0.0, k );

		r0.setStateActionGuess("x3", 10.0, k );
		r0.setStateActionGuess("y3", 0.0, k );

		r0.setStateActionGuess("s1", 1.0, k );
		r0.setStateActionGuess("s2", 1.0, k );
		r0.setStateActionGuess("s3", 1.0, k );
	// 	r0.setStateActionGuess("x", double(k)/double(r0.N)*20, k);
	// 	r0.setStateActionGuess("z", double(k)/double(r0.N)*0, k);
	// 	r0.setStateActionGuess("vx", 4, k);
	}

	// Create the NLP solver
	SnoptInterface si(ocp);

	si.run();

	// Print the optimal cost
	cout << "optimal time: " << -si.F[0] << endl;
	
	// Print the optimal solution
	// vector<double>xopt(ocp.getBigN());
	// solver.getOutput(xopt,NLP_X_OPT);
	//cout << "optimal solution: " << xopt << endl;

	ocp.writeMatlabOutput( "raptor_param_out", si.x );
	r0.writeMatlabOutput( "r0_out", si.x );

	return 0;
}
