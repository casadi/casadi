// MultipleShooting.hpp
// Greg Horn
// Casadi 2011

#include "MultipleShooting.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace CasADi;
using namespace std;

MultipleShooting::MultipleShooting(string name_,
								   Ode & _ode,
								   SX t0_,
								   SX tf_,
								   int N_,
								   vector<double>&lb_,
								   vector<double>&ub_,
								   vector<double>&guess_,
								   SXMatrix & dv_,
								   int idx0_,
								   map<string,SX> & params_) : ode(_ode), dv(dv_), lb(lb_), ub(ub_),
															   guess(guess_), params(params_)
{
	name = name_;
	t0 = t0_;
	tf = tf_;
	N = N_;
	idx0 = idx0_;

	ode.locked = 1;
}

MultipleShooting::~MultipleShooting() {}


SX MultipleShooting::getOutput(string o, int timeStep)
{
	if (timeStep > N-1){
		cerr << "In SX MultipleShooting::getOutput(string o, int timeStep),\n";
		cerr << "timeStep " << timeStep << " > N-1\n";
		exit(1);
	}

	// dynamics constraint
	SX dt = (tf - t0)/(N - 1);

	SXMatrix xk = getStateMat(timeStep);
	SXMatrix uk = getActionMat(timeStep);
	
	SX tk = t0 + timeStep*dt;
	
	map<string,SX> output = ode.getOutputFromDxdt( xk, uk, params, tk );
	
	return output[o];
}



SXMatrix MultipleShooting::getDynamicsConstraintError(int timeStep, map<string,SX> params)
{
	if (timeStep >= N-1){
		cerr << "In SXMatrix MultipleShooting::getDynamicsConstraintError(int timeStep),\n";
		cerr << "timeStep " << timeStep << " >= N-1\n";
		exit(1);
	}

	// dynamics constraint
	SX dt = (tf - t0)/(N - 1);

	SXMatrix x0 = getStateMat(timeStep);
	SXMatrix x1 = getStateMat(timeStep+1);
	
	SXMatrix u0 = getActionMat(timeStep);
	SXMatrix u1 = getActionMat(timeStep+1);
	
	SX tk   = t0 + timeStep*dt;
	SX tkp1 = t0 + (timeStep+1)*dt;
	
	SXMatrix xErr = x1 - ode.rk4Step( x0, u0, u1, params, tk, tkp1);
	return xErr;
}

// total number of discretized states/actions/params
int MultipleShooting::getBigN()
{
	return N*ode.nxu();
}

// get state at proper timestep
SX MultipleShooting::getState(string x, int timeStep)
{
	return dv.at(getIdx(x, timeStep));
}

// get action at proper timestep
SX MultipleShooting::getAction(string u, int timeStep)
{
	return dv.at(getIdx(u, timeStep));
}

// calculate the index of states/actions at proper timestep
int MultipleShooting::getIdx(string xu, int timeStep)
{
	map<string, int>::const_iterator xIter, uIter;

	xIter = ode.states.find(xu);
	uIter = ode.actions.find(xu);

	if (xIter != ode.states.end()) // state
		return idx0 + timeStep*ode.nxu() + xIter->second;
	else if (uIter != ode.actions.end()) // action
		return idx0 + timeStep*ode.nxu() + ode.nx() + uIter->second;
	else {
		cerr << "Error - \"" << xu << "\" not found in states/actions" << endl;
		throw -1;
	}
	return -1;
}

void MultipleShooting::setStateActionGuess(string xu, double guess_)
{
	for (int k=0; k<N; k++)
		setStateActionGuess(xu, guess_, k);
}

void MultipleShooting::setStateActionGuess(string xu, double guess_, int timeStep)
{
	int idx = getIdx(xu, timeStep);
	guess[idx] = guess_;
	if (guess[idx] < lb[idx] ){
		cerr << "Requested initial guess " << xu << "[" << timeStep << "] == " << guess[idx];
		cerr << " is less than lb[" << idx << "] == " << lb[idx] << endl;
		cerr << "Setting guess " << xu << "[" << timeStep << "] = " << lb[idx] << endl;
		guess[idx] = lb[idx];
	}
	if (guess[idx] > ub[idx] ){
		cerr << "Requested initial guess " << xu << "[" << timeStep << "] == " << guess[idx];
		cerr << " is greater than ub[" << idx << "] == " << ub[idx] << endl;
		cerr << "Setting guess " << xu << "[" << timeStep << "] = " << ub[idx] << endl;
		guess[idx] = ub[idx];
	}
}

void MultipleShooting::boundStateAction(string xu, double lb_, double ub_, int timeStep)
{
	int idx = getIdx(xu, timeStep);
	lb[idx] = lb_;
	ub[idx] = ub_;
	if (guess[idx] < lb[idx]) guess[idx] = lb[idx];
	if (guess[idx] > ub[idx]) guess[idx] = ub[idx];
}

void MultipleShooting::boundStateAction(string xu, double lb_, double ub_)
{
	for (int k=0; k<N; k++)
		boundStateAction(xu, lb_, ub_, k);
}


SXMatrix MultipleShooting::getStateMat(int timeStep)
{
	SXMatrix ret = create_symbolic("a_state", ode.nx());
	// SXMatrix ret(ode.nx(), 1);
	map<string,int>::const_iterator xIter;
	for (xIter = ode.states.begin(); xIter != ode.states.end(); xIter++)
		ret[xIter->second] = getState(xIter->first, timeStep);

	return ret;
}

SXMatrix MultipleShooting::getActionMat(int timeStep)
{
	SXMatrix ret = create_symbolic("an_action", ode.nu());
	//	SXMatrix ret(ode.nu(),1);
	map<string,int>::const_iterator uIter;
	for (uIter = ode.actions.begin(); uIter != ode.actions.end(); uIter++)
		ret[uIter->second] = getAction(uIter->first, timeStep);

	return ret;
}

void MultipleShooting::writeMatlabOutput( const char * filename, double * xOpt)
{
	char filename2[200];
	sprintf(filename2, "%s.m", filename);
	ofstream f(filename2);

	if (!f){
		cerr << "error opening " << filename2 << endl;
		exit(1);
	}
	f.precision(10);

	f << "function [x,u] = " << filename << "()" << endl;

	map<string, int>::const_iterator iter;

	// states
	for (iter = ode.states.begin(); iter != ode.states.end(); iter++){
		f << "x." << iter->first << " = [" ;
		for (int k=0; k<N; k++){
			int idx = getIdx( iter->first, k );
			if (k < N-1)
				f << xOpt[idx] << ", ";
			else
				f << xOpt[idx] << "];" << endl;
		}
	}
	f << endl;

	// actions
	for (iter = ode.actions.begin(); iter != ode.actions.end(); iter++){
		f << "u." << iter->first << " = [" ;
		for (int k=0; k<N; k++){
			int idx = getIdx( iter->first, k );
			if (k < N-1)
				f << xOpt[idx] << ", ";
			else
				f << xOpt[idx] << "];" << endl;
		}
	}
	f << endl;

	// // params
	// for (iter = ode.params.begin(); iter != ode.params.end(); iter++){
	// 	int idx = getParamIdx( iter->first );
	// 	f << "p." << iter->first << " = " << xOpt[idx] << ";" << endl;
	// }
	// f << endl;

	f.close();
}
