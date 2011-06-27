// lqrTest.cpp

#include <iostream>
#include <cstdlib>

#include <casadi/stl_vector_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/fx/sx_function.hpp>
//#include <casadi/fx/jacobian.hpp>

#include "Ode.hpp"
#include "Ocp.hpp"
#include "MultipleShooting.hpp"
#include "Lqr.hpp"

#include <SnoptInterface.hpp>

#include <string>
#include <map>

using namespace CasADi;
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
	outputs["bob_y"]  = x - l*cos(theta);
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
cost(map<string,SX> state, map<string,SX> action, int timestep, int N)
{

	SX cost;
	if (timestep == N-1){
		cost = SQR(state["x"]) + SQR(state["theta"]) + SQR(state["vx"]) + SQR(state["vtheta"]);
		return cost;
	}

	//	cost = 2*SQR(state["x"]) + 3*SQR(state["theta"]) + 4*SQR(state["vx"]) + 5*SQR(state["vtheta"]) + 6*SQR(action["u"]);
	cost = 0.1*SQR(action["u"]);
	
	return cost;
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
	ocp.objFun = 50*cos(thetaf) + 5*vthetaf*vthetaf;
	for (int k=0; k<ms.N-1; k++){
		SX u_k   = ms.getAction( "u", k );
		SX u_kp1 = ms.getAction( "u", k+1 );
		// ocp.objFun += 3*SQR( u_k - u_kp1 );
		ocp.objFun += 3*SQR( u_k );
	}

	// bounds
	ocp.boundParam("tEnd", 6, 6);
	
	ms.boundStateAction("x", -trackLength/2, trackLength/2);
	ms.boundStateAction("vx", -22, 22);
	ms.boundStateAction("theta", -50, 50);
	ms.boundStateAction("vtheta", -50, 50);

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

	// initial gues
	for (int k=0; k<ms.N; k++){
		double zero_to_one = k/(ms.N - 1.0);
		ms.setStateActionGuess( "u", sin(2*M_PI*zero_to_one), k);
	}

	// solve
	SnoptInterface si(ocp);
	
	si.run();
	ocp.writeOctaveOutput("cartpole_out");
	
	// Print the optimal cost
	cout << "optimal objFcn: " << si.F[0] << endl;
	
	
	/************** run lqr ****************/
	double t0 = 0;
	double tf = ocp.getParamSolution("tEnd");
 	Lqr lqr(ode, t0, tf, ms.N, &cost);

	for (int k=0; k<ms.N; k++)
		lqr.x_trajectory.at(k) = ocp.getStateSolution(k);
	for (int k=0; k<ms.N-1; k++)
		lqr.u_trajectory.at(k) = ocp.getActionSolution(k);


	// for (int k=0; k<N; k++)
	// 	cout << "lqr.x_trajectory.at(" << k << "):\n" << lqr.x_trajectory.at(k) << endl;
	// for (int k=0; k<N-1; k++)
	// 	cout << "lqr.u_trajectory.at(" << k << "):\n" << lqr.u_trajectory.at(k) << endl;


	lqr.runBackwardSweep();
	
	cout << "successful finish\n";
	return 0;
}
