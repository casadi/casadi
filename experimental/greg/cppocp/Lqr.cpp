// Lqr.cpp
// Greg Horn
// Casadi 2011

#include "Lqr.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace CasADi;
using namespace std;

Lqr::Lqr(Ode & _ode, double t0_, double tf_, int N_, SX (*cost_)(map<string,SX> state, map<string,SX> action)) : ode(_ode)
{
	ode.locked = 1;
	t0 = t0_;
	tf = tf_;
	N = N_;

	costFcnExt = cost_;
	
	setupFunctions();
}

Lqr::~Lqr() {}

void Lqr::setupFunctions()
{
	// inputs
	SXMatrix xk = create_symbolic("xk", ode.nx());
	SXMatrix uk = create_symbolic("uk", ode.nu());

	vector<SXMatrix> inputs_xu(2);
	inputs_xu[0] = xk;
	inputs_xu[1] = uk;

	map<string,SX> xkMap = ode.getStateMap(xk);
	map<string,SX> ukMap = ode.getActionMap(uk);


	/**************** cost function *********************/
	SXMatrix c(costFcnExt( xkMap, ukMap ));

	// function
	cost = SXFunction( inputs_xu, c );
	cost.init();

	// jacobian
	cost_x = cost.jacobian(0);
	cost_u = cost.jacobian(1);
	cost_x.init();
	cost_u.init();

	// hessian
	cost_xx = cost_x.jacobian(0);
	cost_xu = cost_x.jacobian(1); // == cost_u.jacobian(0).trans()
	cost_uu = cost_u.jacobian(1);
	cost_xx.init();
	cost_xu.init();
	cost_uu.init();

	// cout << "cost.outputSX():\n" << cost.outputSX() << endl;
	// cout << "cost_x.outputSX():\n" << cost_x.outputSX() << endl;
	// cout << "cost_u.outputSX():\n" << cost_u.outputSX() << endl;
	// cout << "cost_xx.outputSX():\n" << cost_xx.outputSX() << endl;
	// cout << "cost_xu.outputSX():\n" << cost_xu.outputSX() << endl;
	// cout << "cost_uu.outputSX():\n" << cost_uu.outputSX() << endl;


	/**************** dynamics function *********************/
	// dummy params for now
	map<string,SX> dummyParams;
	double dt = (tf - t0)/(N - 1);
	//SXMatrix xkp1 = ode.rk4Step( xk, uk, uk, dummyParams, t0, dt); // different one for each time step for f(x,u,__t__) - same w cost
	SXMatrix xkp1 = ode.eulerStep( xk, uk, dummyParams, SX(t0), SX(dt)); // different one for each time step for f(x,u,__t__) - same w cost

	f = SXFunction( inputs_xu, xkp1 );
	f.init();

	f_x = f.jacobian(0);
	f_u = f.jacobian(1);
	f_x.init();
	f_u.init();
}

SXMatrix Lqr::takeBackwardsStep(SXMatrix Vxx_kp1, SXMatrix Vx_kp1, SXMatrix V0_kp1)
{
	
	
}



// SXMatrix Lqr::getOutput(string o)
// {
// 	SXMatrix ret = create_symbolic(o, N);
// 	for (int k=0; k<N; k++)
// 		ret.at(k) = getOutput(o, k);
	
// 	return ret;
// }

// SX Lqr::getOutput(string o, int timeStep)
// {
// 	if (timeStep > N-1){
// 		cerr << "In SX Lqr::getOutput(string o, int timeStep),\n";
// 		cerr << "timeStep " << timeStep << " > N-1\n";
// 		exit(1);
// 	}

// 	// dynamics constraint
// 	SX dt = (tf - t0)/(N - 1);

// 	SXMatrix xk = getStateMat(timeStep);
// 	SXMatrix uk = getActionMat(timeStep);
	
// 	SX tk = t0 + timeStep*dt;
	
// 	map<string,SX> output = ode.getOutputFromDxdt( xk, uk, params, tk );

// 	// make sure output exists
// 	map<string, SX>::const_iterator oIter;
// 	oIter = output.find(o);
// 	if (oIter == output.end()){
// 		cerr << "Error - SX Lqr::getOutput(string o, int timeStep) could not find output \"" << o << "\"\n";
// 		throw 1;
// 	}
	
// 	return output[o];
// }

// SXMatrix Lqr::getDynamicsConstraintError(int timeStep)
// {
// 	if (timeStep > N-2){
// 		cerr << "In SXMatrix Lqr::getDynamicsConstraintError(int timeStep),\n";
// 		cerr << "timeStep: " << timeStep << " > N-2   (N==" << N << ")\n";
// 		exit(1);
// 	}

// 	// dynamics constraint
// 	SX dt = (tf - t0)/(N - 1);

// 	SXMatrix x0 = getStateMat(timeStep);
// 	SXMatrix x1 = getStateMat(timeStep + 1);
	
// 	SXMatrix u0 = getActionMat(timeStep);
// 	SXMatrix u1 = getActionMat(timeStep + 1);
	
// 	SX tk0 = t0 + timeStep*dt;
	
// 	//SXMatrix xErr = x1 - ode.rk4Step( x0, u0, u1, params, tk0, dt);
// 	//SXMatrix xErr = x1 - ode.eulerStep( x0, u0, params, tk0, dt);
// 	SXMatrix xErr = ode.simpsonsRuleError( x0, x1, u0, u1, params, tk0, dt);

// 	return xErr;
// }

// // total number of discretized states/actions/params
// int Lqr::getBigN()
// {
// 	return N*ode.nxu();
// }

// // get state at proper timestep
// SX Lqr::getState(string x, int timeStep)
// {
// 	return dv.at(getIdx(x, timeStep));
// }

// // get action at proper timestep
// SX Lqr::getAction(string u, int timeStep)
// {
// 	return dv.at(getIdx(u, timeStep));
// }

// // calculate the index of states/actions at proper timestep
// int Lqr::getIdx(string xu, int timeStep)
// {
// 	map<string, int>::const_iterator xIter, uIter;

// 	xIter = ode.states.find(xu);
// 	uIter = ode.actions.find(xu);

// 	if (xIter != ode.states.end()) // state
// 		return idx0 + timeStep*ode.nxu() + xIter->second;
// 	else if (uIter != ode.actions.end()) // action
// 		return idx0 + timeStep*ode.nxu() + ode.nx() + uIter->second;
// 	else {
// 		cerr << "Error - \"" << xu << "\" not found in states/actions" << endl;
// 		throw -1;
// 	}
// 	return -1;
// }

// void Lqr::setStateActionGuess(string xu, double guess_)
// {
// 	for (int k=0; k<N; k++)
// 		setStateActionGuess(xu, guess_, k);
// }

// void Lqr::setStateActionGuess(string xu, double guess_, int timeStep)
// {
// 	int idx = getIdx(xu, timeStep);
// 	guess[idx] = guess_;
// 	if (guess[idx] < lb[idx] ){
// 		cerr << "Requested initial guess " << xu << "[" << timeStep << "] == " << guess[idx];
// 		cerr << " is less than lb[" << idx << "] == " << lb[idx] << endl;
// 		cerr << "Setting guess " << xu << "[" << timeStep << "] = " << lb[idx] << endl;
// 		guess[idx] = lb[idx];
// 	}
// 	if (guess[idx] > ub[idx] ){
// 		cerr << "Requested initial guess " << xu << "[" << timeStep << "] == " << guess[idx];
// 		cerr << " is greater than ub[" << idx << "] == " << ub[idx] << endl;
// 		cerr << "Setting guess " << xu << "[" << timeStep << "] = " << ub[idx] << endl;
// 		guess[idx] = ub[idx];
// 	}
// }

// void Lqr::boundStateAction(string xu, double lb_, double ub_, int timeStep)
// {
// 	int idx = getIdx(xu, timeStep);
// 	lb[idx] = lb_;
// 	ub[idx] = ub_;
// 	if (guess[idx] < lb[idx]) guess[idx] = lb[idx];
// 	if (guess[idx] > ub[idx]) guess[idx] = ub[idx];
// }

// void Lqr::boundStateAction(string xu, double lb_, double ub_)
// {
// 	for (int k=0; k<N; k++)
// 		boundStateAction(xu, lb_, ub_, k);
// }


// SXMatrix Lqr::getStateMat(int timeStep)
// {
// 	SXMatrix ret = create_symbolic("a_state", ode.nx());
// 	// SXMatrix ret(ode.nx(), 1);
// 	map<string,int>::const_iterator xIter;
// 	for (xIter = ode.states.begin(); xIter != ode.states.end(); xIter++)
// 		ret[xIter->second] = getState(xIter->first, timeStep);

// 	return ret;
// }

// SXMatrix Lqr::getActionMat(int timeStep)
// {
// 	SXMatrix ret = create_symbolic("an_action", ode.nu());
// 	//	SXMatrix ret(ode.nu(),1);
// 	map<string,int>::const_iterator uIter;
// 	for (uIter = ode.actions.begin(); uIter != ode.actions.end(); uIter++)
// 		ret[uIter->second] = getAction(uIter->first, timeStep);

// 	return ret;
// }

// void Lqr::writeOctaveOutput( ofstream & f, double * xOpt )
// {
// 	f.precision(10);

// 	f << "function multipleShooting = ms_stage_" << name << "_out()" << endl;

// 	map<string, int>::const_iterator iter;

// 	// states
// 	f << "% states\n";
// 	f << "multipleShooting.states = struct();\n";
// 	for (iter = ode.states.begin(); iter != ode.states.end(); iter++){
// 		f << "multipleShooting.states." << iter->first << " = [" ;
// 		for (int k=0; k<N; k++){
// 			int idx = getIdx( iter->first, k );
// 			if (k < N-1)
// 				f << xOpt[idx] << ", ";
// 			else
// 				f << xOpt[idx] << "];" << endl;
// 		}
// 	}
// 	f << endl;


// 	// actions
// 	f << "% actions\n";
// 	f << "multipleShooting.actions = struct();\n";
// 	for (iter = ode.actions.begin(); iter != ode.actions.end(); iter++){
// 		f << "multipleShooting.actions." << iter->first << " = [" ;
// 		for (int k=0; k<N; k++){
// 			int idx = getIdx( iter->first, k );
// 			if (k < N-1)
// 				f << xOpt[idx] << ", ";
// 			else
// 				f << xOpt[idx] << "];" << endl;
// 		}
// 	}
// 	f << endl;


// 	// outputs
// 	f << "% outputs\n";
// 	f << "multipleShooting.outputs = struct();\n";
// 	for (iter = ode.outputs.begin(); iter != ode.outputs.end(); iter++){

// 		SXMatrix outSXMatrix = getOutput(iter->first);
// 		SXFunction outputFcn(dv, outSXMatrix);
// 		outputFcn.init();
// 		outputFcn.setInput(xOpt);
// 		outputFcn.evaluate();
// 		vector<double>outDouble(N);
// 		outputFcn.getOutput(outDouble);

// 		f << "multipleShooting.outputs." << iter->first << " = [";
// 		for (int k=0; k<N; k++){
// 			if (k < N - 1)
// 				f << outDouble[k] << ", ";
// 			else
// 				f << outDouble[k] << "];" << endl;
// 		}
// 	}
// 	f << endl;


// 	// start/end times
// 	f << "% time\n";
// 	SXFunction timeFcn( dv, vertcat( SXMatrix(t0), SXMatrix(tf) ) );
// 	timeFcn.init();
// 	timeFcn.setInput( xOpt );
// 	timeFcn.evaluate();
// 	double timeNum[2];
// 	timeFcn.getOutput( timeNum );

// 	f << "multipleShooting.time = linspace(" << timeNum[0] << ", " << timeNum[1] << ", " << N << ");\n\n";
// }
