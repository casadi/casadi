// OcpMultipleShooting.hpp
// Greg Horn
// Casadi 2011

#include "OcpMultipleShooting.hpp"

using namespace CasADi;
using namespace std;

OcpMultipleShooting::OcpMultipleShooting(Ode * _ode) : Ocp(_ode) {}

OcpMultipleShooting::~OcpMultipleShooting() {}

void OcpMultipleShooting::setTimeInterval(CasADi::SX _t0, CasADi::SX _tf)
{
	t0 = _t0;
	tf = _tf;

	// dynamics constraint
	SX dt = (tf - t0)/(N - 1);

	vector<SX>::const_iterator iter;

	for (int k=0; k<N-1; k++){
		SXMatrix x0 = getStateMat(k);
		SXMatrix x1 = getStateMat(k+1);
		
		SXMatrix u0 = getActionMat(k);
		SXMatrix u1 = getActionMat(k+1);

		SXMatrix p = getParamMat();

		SX tk   = k*dt;
		SX tkp1 = (k+1)*dt;
		
		SXMatrix xErr = x1 - ode->rk4Step( x0, u0, u1, p, tk, tkp1);

		addNonlconEq( xErr );
	}
}

void OcpMultipleShooting::discretize(int _N)
{
	ode->assertUnlocked();
	ode->locked = 1;

	N = _N;

	designVariables = create_symbolic("designVariables", getBigN());

	for (int k=0; k<getBigN(); k++){
		lb.push_back(-1e50);
		ub.push_back( 1e50);
		//	-numeric_limits<double>::infinity()

		guess.push_back(0);
	}
}


// total number of discretized states/actions/params
int OcpMultipleShooting::getBigN()
{
	return N*ode->nxu() + ode->np();
}

// get state at proper timestep
SX OcpMultipleShooting::getState(string x, int timeStep)
{
	return designVariables[getStateActionIdx(x, timeStep)];
}

// get action at proper timestep
SX OcpMultipleShooting::getAction(string u, int timeStep)
{
	return designVariables[getStateActionIdx(u, timeStep)];
}

// get param at proper timestep
SX OcpMultipleShooting::getParam(string p)
{
	return designVariables[getParamIdx(p)];
}

// calculate the index of states/actions at proper timestep
int OcpMultipleShooting::getStateActionIdx(string xu, int timeStep)
{
	map<string, int>::const_iterator xIter, uIter;

	xIter = ode->states.find(xu);
	uIter = ode->actions.find(xu);

	if (xIter != ode->states.end()) // state
		return timeStep*ode->nxu() + xIter->second;
	else if (uIter != ode->actions.end()) // action
		return timeStep*ode->nxu() + ode->nx() + uIter->second;
	else {
		cerr << "Error - \"" << xu << "\" not found in states/actions" << endl;
		throw -1;
	}
	return -1;
}

// calculate the index of params
int OcpMultipleShooting::getParamIdx(string p)
{
	map<string, int>::const_iterator pIter;

	pIter = ode->params.find(p);

	if (pIter != ode->params.end())
		return N*ode->nxu() + pIter->second;
	else {
		cerr << "Error - \"" << p << "\" not found in params" << endl;
		throw -1;
	}
	return -1;
}


void OcpMultipleShooting::boundStateAction(string xu, double _lb, double _ub, int timeStep)
{
	int idx = getStateActionIdx(xu, timeStep);
	lb[idx] = _lb;
	ub[idx] = _ub;
}


void OcpMultipleShooting::boundParam(string p, double _lb, double _ub)
{
	int idx = getParamIdx(p);
	lb[idx] = _lb;
	ub[idx] = _ub;
}

SXMatrix OcpMultipleShooting::getStateMat(int timeStep)
{
	SXMatrix ret = create_symbolic("a_state", ode->nx());
	// SXMatrix ret(ode->nx(), 1);
	map<string,int>::const_iterator xIter;
	for (xIter = ode->states.begin(); xIter != ode->states.end(); xIter++)
		ret[xIter->second] = getState(xIter->first, timeStep);

	return ret;
}

SXMatrix OcpMultipleShooting::getActionMat(int timeStep)
{
	SXMatrix ret = create_symbolic("an_action", ode->nu());
	//	SXMatrix ret(ode->nu(),1);
	map<string,int>::const_iterator uIter;
	for (uIter = ode->actions.begin(); uIter != ode->actions.end(); uIter++)
		ret[uIter->second] = getAction(uIter->first, timeStep);

	return ret;
}

SXMatrix OcpMultipleShooting::getParamMat()
{
	SXMatrix ret = create_symbolic("a_param", ode->np());
	// SXMatrix ret(ude->np(), 1);
	map<string,int>::const_iterator pIter;
	for (pIter = ode->params.begin(); pIter != ode->params.end(); pIter++)
		ret[pIter->second] = getParam(pIter->first);

	return ret;
}

