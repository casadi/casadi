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

void
dxdt(map<string,SX> &xDot, map<string,SX> &outputs, map<string,SX> state, map<string,SX> action, map<string,SX> param, SX t)
{
	// constants
	double AR  =     6; // aspect ratio
	double Cd0 =  0.03; // parasitic drag
	double m   =   0.4; // mass
	double rho =  1.22; // air density
	double A   = 0.625; // area
	double g   =   9.8; // acceleration due to gravity
	
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

	outputs["airspeed"] = norm_v;
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
	
	ode.addOutput("airspeed");

	ode.dxdt = &dxdt;

	return ode;
}

int
main()
{
	Ode ode = getGliderOde();
	Ocp ocp;
	SX tEnd = ocp.addParam("tEnd");

	MultipleShooting & msTakeoff = ocp.addMultipleShooting("takeoff", ode, 0.0, 6, 100);
	MultipleShooting & msGlide = ocp.addMultipleShooting("glide", ode, 6, tEnd, 100);
	
	ocp.objFun = -tEnd; // maximum time

	// tie stages together
	SX xTakeoffEnd = msTakeoff.getStateMat(msTakeoff.N-1);
	SX uTakeoffEnd = msTakeoff.getActionMat(msTakeoff.N-1);
	SX xGlideBegin = msGlide.getStateMat(0);
	SX uGlideBegin = msGlide.getActionMat(0);

	ocp.addNonlconEq( xTakeoffEnd - xGlideBegin );
	ocp.addNonlconEq( uTakeoffEnd - uGlideBegin );

	// bounds
	ocp.boundParam("tEnd", 2, 200);

	msTakeoff.boundStateAction("x", 0, 1e3);
	msTakeoff.boundStateAction("z", -300, 0);
	msTakeoff.boundStateAction("vx", 0, 100);
	msTakeoff.boundStateAction("vz", -100, 100);

	msGlide.boundStateAction("x", 0, 1e3);
	msGlide.boundStateAction("z", -300, 0);
	msGlide.boundStateAction("vx", 0, 100);
	msGlide.boundStateAction("vz", -100, 100);

	msTakeoff.boundStateAction("alphaDeg", -15, 15);
	msGlide.boundStateAction("alphaDeg", -15, 15);

	// initial condition
	#define INCLINATION0_DEG 0.0
	#define V0 20.0
	msTakeoff.boundStateAction("x", 0.0, 0.0, 0);
	msTakeoff.boundStateAction("z", 0.0, 0.0, 0);
	msTakeoff.boundStateAction("vx",  V0*cos(INCLINATION0_DEG*M_PI/180.0),  V0*cos(INCLINATION0_DEG*M_PI/180.0),  0);
	msTakeoff.boundStateAction("vz", -V0*sin(INCLINATION0_DEG*M_PI/180.0), -V0*sin(INCLINATION0_DEG*M_PI/180.0),  0);

	// // initial guesses
	ocp.setParamGuess("tEnd", 80.0);

	for (int k=0; k< msTakeoff.N; k++){
		msTakeoff.setStateActionGuess("x", double(k)/double(msTakeoff.N)*20, k);
		msTakeoff.setStateActionGuess("z", -double(k)/double(msTakeoff.N)*15, k);
		msTakeoff.setStateActionGuess("vx", 4.0, k);
	}

	for (int k=0; k< msTakeoff.N; k++){
		msGlide.setStateActionGuess("x", 20 + double(k)/double(msTakeoff.N)*20, k);
		msGlide.setStateActionGuess("z", -double(msTakeoff.N - k - 1)/double(msTakeoff.N)*15, k);
		msGlide.setStateActionGuess("vx", 4.0, k);
	}

	// solve
	SnoptInterface si(ocp);

	si.run();

	// Print the optimal cost
	cout << "optimal time: " << -si.F[0] << endl;

	// write output
	ocp.writeOctaveOutput( "glider_out" );
	
	return 0;
}
