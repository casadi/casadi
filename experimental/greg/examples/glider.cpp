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

#include <interfaces/ipopt/ipopt_solver.hpp>

#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"

// #include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace casadi;
using namespace std;

void
dxdt(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
	// constants
	double AR = 6;     // aspect ration
	double Cd0 = 0.03; // parasitic drag
	double m = 2.0;    // mass
	double rho = 1.22; // air density
	double A = 1.0;    // area
	double g = 9.8;    // acceleration due to gravity
	
	// eom
	SX alpha = action["alphaDeg"]*3.14159/180.0;
	SX CL = 2.0*3.14159*alpha;
	SX Cd = Cd0 + 1.1*CL*CL/(3.14159*AR);

	// airspeed
	SX vx = state["vx"];
	SX vz = state["vz"];
	SX norm_v = sqrt( vx*vx + vz*vz );
        
//        Lift = 0.5*rho*A*norm_v*CL*[ vz, -vx];
//        Drag = 0.5*rho*A*norm_v*Cd*[-vx, -vz];
	SX cAero = 0.5*rho*A*norm_v;
	SX Fx = cAero*( CL*vz - Cd*vx);
	SX Fz = cAero*(-CL*vx - Cd*vz);
	SX ax = Fx/m;
	SX az = Fz/m + g;

	xDot["x"] = vx;
	xDot["z"] = vz;
	xDot["vx"] = ax;
	xDot["vz"] = az;

	//	intermediate["airspeed"] = norm_v;
}

Ode
getGliderOde()
{
	Ode ode("glider");
	ode.addState("x");
	ode.addState("z");
	ode.addState("vx");
	ode.addState("vz");

	ode.addAction("alphaDeg");

	ode.dxdt = &dxdt;

	return ode;
}


int
main()
{
	Ode ode = getGliderOde();
	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");

	MultipleShooting & ms = ocp.addMultipleShooting("glider", ode, 0, tEnd, 100);
	
	ocp.objFun = -tEnd; // maximum time

	// Bounds/initial condition
	ocp.boundParam("tEnd", 2, 200);

	for (int k=0; k<ms.N; k++){
		ms.boundStateAction("x", 0, 1e3, k);
		ms.boundStateAction("z", -100, 0, k);
		ms.boundStateAction("vx", 0, 100, k);
		ms.boundStateAction("vz", -100, 100, k);

		ms.boundStateAction("alphaDeg", -15, 15, k);
	}

	#define INCLINATION0_DEG 0.0
	#define V0 20.0
	ms.boundStateAction("x", 0, 0, 0);
	ms.boundStateAction("z", 0, 0, 0);
	ms.boundStateAction("vx",  V0*cos(INCLINATION0_DEG*M_PI/180.0),  V0*cos(INCLINATION0_DEG*M_PI/180.0),  0);
	ms.boundStateAction("vz", -V0*sin(INCLINATION0_DEG*M_PI/180.0), -V0*sin(INCLINATION0_DEG*M_PI/180.0),  0);

	// initial guesses
	ocp.setParamGuess("tEnd", 80);
	for (int k=0; k<ms.N; k++){
		ms.setStateActionGuess("x", double(k)/double(ms.N)*20, k);
		ms.setStateActionGuess("z", double(k)/double(ms.N)*0, k);
		ms.setStateActionGuess("vx", 4, k);
	}

	// Create the NLP solver
	SXFunction ffcn(ocp.designVariables, ocp.objFun); // objective function
	SXFunction gfcn(ocp.designVariables, ocp.g); // constraint
	gfcn.setOption("ad_mode","reverse");
	gfcn.setOption("symbolic_jacobian",false);

	IpoptSolver solver(ffcn,gfcn);
	//IpoptSolver solver(ffcn,gfcn,Function(),Jacobian(gfcn));

	// Set options
	solver.setOption("tol",1e-10);
	solver.setOption("hessian_approximation","limited-memory");

	// initialize the solver
	solver.init();

	solver.setInput(    ocp.lb, "lbx");
	solver.setInput(    ocp.ub, "ubx");
	solver.setInput( ocp.guess, "x0");

	// Bounds on g
	solver.setInput( ocp.gMin, "lbg");
	solver.setInput( ocp.gMax, "ubg");

	// Solve the problem
	solver.solve();

	// Print the optimal cost
	double cost;
	solver.getOutput(cost,"f");
	cout << "optimal time: " << cost << endl;

	// Print the optimal solution
	vector<double>xopt(ocp.designVariables.size1());
	solver.getOutput(xopt,"x");
	cout << "optimal solution: " << xopt << endl;

	ocp.writeOctaveOutput( "glider_out" );

	return 0;
}
