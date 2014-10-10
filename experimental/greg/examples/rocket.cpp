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

#include <core/std_vector_tools.hpp>
#include <core/sx/sx_tools.hpp>
#include <core/function/sx_function.hpp>
#include <core/function/jacobian.hpp>

#include <interfaces/ipopt/ipopt_solver.hpp>


#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"

#include <string>
#include <map>

using namespace casadi;
using namespace std;

void
dxdt(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
	SX x = state["x"];
	SX v = state["v"];
	SX mass = state["mass"];
	SX thrust = action["thrust"];

	double massBurnedPerThrustSquared = 1e-6;

	xDot["x"] = v;
	xDot["v"] = thrust/mass;
	xDot["mass"] = -thrust*thrust*massBurnedPerThrustSquared;
}

int
main()
{
	Ode ode("rocket");
	ode.addState("x");
	ode.addState("v");
	ode.addState("mass");
	ode.addAction("thrust");

	ode.dxdt = &dxdt;

	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");
	MultipleShooting & ms = ocp.addMultipleShooting("rocket", ode, 0, tEnd, 80);

	ocp.objFun = tEnd;

	// Bounds/initial condition
	double x0 = 0;
	double xf = 1000;
	double massFuel0 = 50;
	double massShip = 100;
	ocp.boundParam("tEnd", 1, 1000);
	for (int k=0; k<ms.N; k++){
		ms.boundStateAction("x", x0, xf, k);
		ms.boundStateAction("v", -1000, 1000, k);
		ms.boundStateAction("mass", massShip, massShip + massFuel0, k);
		ms.boundStateAction("thrust", -3000, 3000, k);
	}

	// initial mass
	ms.boundStateAction("mass", massShip + massFuel0, massShip + massFuel0, 0);

	// initial/final position
	ms.boundStateAction("x", x0, x0, 0);
	ms.boundStateAction("x", xf, xf, ms.N-1);

	// velocity == 0 at start/finish
	ms.boundStateAction("v", 0, 0, 0);
	ms.boundStateAction("v", 0, 0, ms.N-1);

	// hessian
	SX sigma = ssym("sigma");
	SX lambda = ssym("lambda", ocp.g.size1());
	SX lagrangian = sigma.at(0)*ocp.objFun;
	for (int k=0; k<ocp.g.size1(); k++)
		lagrangian += lambda.at(k)*ocp.g.at(k);

	SX h = hessian(lagrangian, ocp.designVariables);

	vector<SX> inputs(3);
	inputs[0] = ocp.designVariables;
	inputs[1] = lambda;
	inputs[2] = sigma;

	// Create the NLP solver
	SXFunction ffcn(ocp.designVariables, ocp.objFun); // objective function
	SXFunction gfcn(ocp.designVariables, ocp.g); // constraint
	SXFunction hfcn(inputs, h); // hessian

	gfcn.setOption("ad_mode","reverse");
	gfcn.setOption("symbolic_jacobian",false);

	//IpoptSolver solver(ffcn,gfcn);
	//IpoptSolver solver(ffcn,gfcn,Function(),Jacobian(gfcn));
	IpoptSolver solver( ffcn, gfcn, hfcn,  Function());

	// Set options
	solver.setOption("tol",1e-8);
	//solver.setOption("hessian_approximation","limited-memory");
	solver.setOption("hessian_approximation","exact");

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

	// // Print the optimal solution
	vector<double>xopt(ocp.designVariables.size1());
	solver.getOutput(xopt,"x");
	// cout << "optimal solution: " << xopt << endl;

	ocp.writeOctaveOutput( "rocket_out" );
	return 0;
}
