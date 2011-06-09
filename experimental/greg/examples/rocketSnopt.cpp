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

#include <casadi/stl_vector_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/fx/sx_function.hpp>


#include <Ode.hpp>
#include <Ocp.hpp>
#include <OcpMultipleShooting.hpp>
#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace CasADi;
using namespace std;

void
dxdt(map<string,SX> &xDot, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t __attribute__((unused)))
{
	SX x = state["x"];
	SX v = state["v"];
	SX mass = state["mass"];
	SX thrust = action["thrust"];
	SX tEnd = param["tEnd"];

	double massBurnedPerThrustSquared = 1.0;

	xDot["x"] = v;
	xDot["v"] = thrust/mass;
	xDot["mass"] = -thrust*thrust*massBurnedPerThrustSquared;
}

OcpMultipleShooting
getRocketOcp(void)
{
	Ode ode("rocket");
	ode.addState("x");
	ode.addState("v");
	ode.addState("mass");
	ode.addAction("thrust");
	ode.addParam("tEnd");

	ode.dxdt = &dxdt;

	OcpMultipleShooting ocp(&ode);

	//ocp.discretize(150);
	ocp.discretize(80);
	//	ocp.discretize(10);

	SX tEnd = ocp.getParam("tEnd");
	ocp.setTimeInterval(0.0, tEnd);
	ocp.objFun = tEnd;

	// Bounds/initial condition
	double x0 = 0;
	double xf = 15;
	double mass0 = 150;
	double massShip = 100;
	ocp.boundParam("tEnd", 1, 1000);
	for (int k=0; k<ocp.N; k++){
		ocp.boundStateAction("x", x0, xf, k);
		ocp.boundStateAction("v", -100, 100, k);
		ocp.boundStateAction("mass", massShip, mass0, k);
		ocp.boundStateAction("thrust", -100, 100, k);
	}

	// initial mass
	ocp.boundStateAction("mass", mass0, mass0, 0);

	// initial/final position
	ocp.boundStateAction("x", x0, x0, 0);
	ocp.boundStateAction("x", xf, xf, ocp.N-1);

	// velocity == 0 at start/finish
	ocp.boundStateAction("v", 0, 0, 0);
	ocp.boundStateAction("v", 0, 0, ocp.N-1);

	return ocp;
}

int
main()
{
	OcpMultipleShooting ocp = getRocketOcp();

	SnoptInterface si(ocp);

	si.run();

	// cout << endl;
	// for (int k=0; k<si.n; k++)
	// 	printf("x[%d]: %f\t", k, si.x[k]);
	// cout << endl;
}
