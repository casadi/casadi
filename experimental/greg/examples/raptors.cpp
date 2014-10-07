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


#include <iostream>
#include <cstdlib>

#include <core/std_vector_tools.hpp>
#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>
//#include <core/function/jacobian.hpp>


#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace casadi;
using namespace std;

#define EPS 1e-12

void
dxdt0(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
    SX xh = state["xh"];
	SX yh = state["yh"];
	SX x1 = state["x1"];
	SX y1 = state["y1"];
	SX x2 = state["x2"];
	SX y2 = state["y2"];
	SX x3 = state["x3"];
	SX y3 = state["y3"];
	SX s1 = 4*t;
	SX s2 = 4*t;
	SX s3 = 4*t;
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

	xDot["theta"] = action["dtheta"];

	outputs["d1"] = d1;
	outputs["d2"] = d2;
	outputs["d3"] = d3;
}

void
dxdt1(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
    SX xh = state["xh"];
	SX yh = state["yh"];
	SX x1 = state["x1"];
	SX y1 = state["y1"];
	SX x2 = state["x2"];
	SX y2 = state["y2"];
	SX x3 = state["x3"];
	SX y3 = state["y3"];
	SX s1 = 10;
	SX s2 = 4*t;
	SX s3 = 4*t;
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

	xDot["theta"] = action["dtheta"];

	outputs["d1"] = d1;
	outputs["d2"] = d2;
	outputs["d3"] = d3;
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

	ode.addState("theta");

	ode.addAction("dtheta");

	ode.addOutput("d1");
	ode.addOutput("d2");
	ode.addOutput("d3");

	return ode;
}

#define DTHETA_MAX_DEG 80.0

int
main()
{
	Ode ode0 = getRaptorOde();
	Ode ode1 = getRaptorOde();

	ode0.dxdt = &dxdt0;
	ode1.dxdt = &dxdt1;

	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");

	MultipleShooting & r0 = ocp.addMultipleShooting("r0", ode0,    0,  2.5, 60);
	MultipleShooting & r1 = ocp.addMultipleShooting("r1", ode1,  2.5, tEnd, 20);

	ocp.objFun = -tEnd; // maximum time

	// Bounds
	ocp.boundParam("tEnd", 2.6, 6.25);
	for (int k=0; k<r0.N; k++){
		r0.boundStateAction( "xh", -100, 100, k );
		r0.boundStateAction( "yh", -100, 100, k );
		r0.boundStateAction( "x1", -100, 100, k );
		r0.boundStateAction( "y1", -100, 100, k );
		r0.boundStateAction( "x2", -100, 100, k );
		r0.boundStateAction( "y2", -100, 100, k );
		r0.boundStateAction( "x3", -100, 100, k );
		r0.boundStateAction( "y3", -100, 100, k );

		r0.boundStateAction( "theta", -10, 10.0, k);
		r0.boundStateAction("dtheta", -DTHETA_MAX_DEG*M_PI/180.0, DTHETA_MAX_DEG*M_PI/180.0, k);

		SX d1 = r0.getOutput("d1", k);
		SX d2 = r0.getOutput("d2", k);
		SX d3 = r0.getOutput("d3", k);

		ocp.addNonlconIneq( 0.1 - d1 );
		ocp.addNonlconIneq( 0.1 - d2 );
		ocp.addNonlconIneq( 0.1 - d3 );
	}


	for (int k=0; k<r1.N; k++){
		r1.boundStateAction( "xh", -100, 100, k );
		r1.boundStateAction( "yh", -100, 100, k );
		r1.boundStateAction( "x1", -100, 100, k );
		r1.boundStateAction( "y1", -100, 100, k );
		r1.boundStateAction( "x2", -100, 100, k );
		r1.boundStateAction( "y2", -100, 100, k );
		r1.boundStateAction( "x3", -100, 100, k );
		r1.boundStateAction( "y3", -100, 100, k );

		r1.boundStateAction( "theta", -10, 10.0, k);
		r1.boundStateAction("dtheta", -DTHETA_MAX_DEG*M_PI/180.0, DTHETA_MAX_DEG*M_PI/180.0, k);

		SX d1 = r1.getOutput("d1", k);
		SX d2 = r1.getOutput("d2", k);
		SX d3 = r1.getOutput("d3", k);

		ocp.addNonlconIneq( 0.1 - d1 );
		ocp.addNonlconIneq( 0.1 - d2 );
		ocp.addNonlconIneq( 0.1 - d3 );
	}


	// join stages
	SX x0finish = r0.getStateMat ( r0.N - 1 );
	SX u0finish = r0.getActionMat( r0.N - 1 );
	SX x1start = r1.getStateMat ( 0 );
	SX u1start = r1.getActionMat( 0 );

	ocp.addNonlconEq( x0finish - x1start );
	ocp.addNonlconEq( u0finish - u1start );


	// initial conditions
	r0.boundStateAction("xh", 0, 0, 0);
	r0.boundStateAction("yh", 10.0*sqrt(3.0)/3.0, 10.0*sqrt(3.0)/3.0, 0);

	r0.boundStateAction("x1", 0.0, 0.0, 0 );
	r0.boundStateAction("y1", 20.0*sqrt(3.0)/2.0, 20.0*sqrt(3.0)/2.0, 0 );

	r0.boundStateAction("x2", -10.0, -10.0, 0 );
	r0.boundStateAction("y2", 0.0, 0.0, 0 );

	r0.boundStateAction("x3", 10.0, 10.0, 0 );
	r0.boundStateAction("y3", 0.0, 0.0, 0 );

	// r0.boundStateAction("s1", 0.0, 0.0, 0 );
	// r0.boundStateAction("s2", 0.0, 0.0, 0 );
	// r0.boundStateAction("s3", 0.0, 0.0, 0 );


	// initial guesses
	ocp.setParamGuess("tEnd", 3.0);

	r0.setStateActionGuess("theta", 45*M_PI/180.0);
	r1.setStateActionGuess("theta", 45*M_PI/180.0);

	for (int k=0; k<r0.N; k++){
		r0.setStateActionGuess("xh", 0, k);
		r0.setStateActionGuess("yh", 10.0*sqrt(3.0)/3.0, k);

		r0.setStateActionGuess("x1", 0.0, k );
		r0.setStateActionGuess("y1", 20.0*sqrt(3.0)/2.0, k );

		r0.setStateActionGuess("x2", -10.0, k );
		r0.setStateActionGuess("y2", 0.0, k );

		r0.setStateActionGuess("x3", 10.0, k );
		r0.setStateActionGuess("y3", 0.0, k );

		// r0.setStateActionGuess("s1", 5.0, k );
		// r0.setStateActionGuess("s2", 5.0, k );
		// r0.setStateActionGuess("s3", 5.0, k );
	}

	for (int k=0; k<r1.N; k++){
		r1.setStateActionGuess("xh", 0, k);
		r1.setStateActionGuess("yh", 10.0*sqrt(3.0)/3.0, k);

		r1.setStateActionGuess("x1", 0.0, k );
		r1.setStateActionGuess("y1", 20.0*sqrt(3.0)/2.0, k );

		r1.setStateActionGuess("x2", -10.0, k );
		r1.setStateActionGuess("y2", 0.0, k );

		r1.setStateActionGuess("x3", 10.0, k );
		r1.setStateActionGuess("y3", 0.0, k );

		// r1.setStateActionGuess("s1", 10.0, k );
		// r1.setStateActionGuess("s2", 12.5, k );
		// r1.setStateActionGuess("s3", 12.5, k );
	}

	// Create the NLP solver
	SnoptInterface si(ocp);

	si.run();

	// Print the optimal cost
	cout << "optimal time: " << -si.F[0] << endl;
	
	ocp.writeOctaveOutput( "raptor_out" );

	return 0;
}
