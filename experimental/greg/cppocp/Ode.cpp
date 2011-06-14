// Ode.cpp
// Greg Horn
// Casadi 2011

#include <iostream>

#include "Ode.hpp"

using namespace std;
using namespace CasADi;

Ode::Ode(string _name)
{
	name = _name;
	locked = 0;
	dxdt = NULL;
}

Ode::~Ode(){}

// get state/action/param dimensions 
int Ode::nx(){ return  states.size(); }
int Ode::nu(){ return actions.size(); }
int Ode::no(){ return outputs.size(); }
int Ode::nxu(){ return nx()+nu();}

// throw error if ode has already been discretized
void Ode::assertUnlocked()
{
	if (locked){
		cerr << "Error - Ode \"" << name << "\" has already been discretized" << endl;
		throw 1;
	}
}

// throw error if name is not unique
void Ode::assertUniqueName(string newName)
{
	if (isState(newName) || isAction(newName) || isOutput(newName)) {
		cerr << "Error - new state/action/output \"" << newName << "\" is not unique" << endl;
		throw "1";
	}
}

int Ode::isState(string stateName)
{
	map<string, int>::const_iterator xIter;

	xIter = states.find(stateName);

	if (xIter != states.end())
		return 1;
	return 0;
}

int Ode::isAction(string actionName)
{
	map<string, int>::const_iterator uIter;

	uIter = actions.find(actionName);

	if (uIter != actions.end())
		return 1;
	return 0;
}

int Ode::isOutput(string outputName)
{
	map<string, int>::const_iterator oIter;

	oIter = outputs.find(outputName);

	if (oIter != outputs.end())
		return 1;
	return 0;
}


// add new states/actions/params
void Ode::addState(string _newState)
{
	assertUnlocked();
	assertUniqueName(_newState);
	states[_newState] = nx() - 1;
}

void Ode::addAction(string _newAction)
{
	assertUnlocked();
	assertUniqueName(_newAction);
	actions[_newAction] = nu() - 1;
}

void Ode::addOutput(string _newOutput)
{
	assertUnlocked();
	assertUniqueName(_newOutput);
	outputs[_newOutput] = no() - 1;
}

map<string,SX> Ode::getStateMap( SXMatrix & x )
{
	map<string,int>::const_iterator iter;
	map<string,SX> xMap;
	for (iter = states.begin(); iter != states.end(); iter++)
		xMap[iter->first] = x.at(iter->second);

	return xMap;
}

map<string,SX> Ode::getActionMap( SXMatrix & u )
{
	map<string,int>::const_iterator iter;
	map<string,SX> uMap;
	for (iter = actions.begin(); iter != actions.end(); iter++)
		uMap[iter->first] = u.at(iter->second);

	return uMap;
}

map<string,SX> Ode::getOutputFromDxdt( SXMatrix x, SXMatrix u, map<string,SX> & p, SX t )
{
	// state/action maps to wrap input SXMatrices
	map<string,SX> xMap = getStateMap(x);
	map<string,SX> uMap = getActionMap(u);

	// dummy xDot to be discarded
	map<string,SX> dummyXDotMap;

	// output map
	map<string,SX> outputMap;

	// call dxdt
	dxdt( dummyXDotMap, outputMap, xMap, uMap, p, t );

	return outputMap;
}


SXMatrix Ode::dxVectorDt( SXMatrix x, SXMatrix u, map<string,SX> & p, SX t )
{
	// state/action maps to wrap input SXMatrices
	map<string,SX> xMap = getStateMap(x);
	map<string,SX> uMap = getActionMap(u);

	// xDot will be populated by dxdt()
	map<string,SX> xDotMap;

	// dummy outputs to be discarded
	map<string,SX> dummyOutputMap;

	// call dxdt
	dxdt( xDotMap, dummyOutputMap, xMap, uMap, p, t );

	// make output SXMatrix
	SXMatrix xDotMat = create_symbolic("an_xDotMat", nx());
	map<string,int>::const_iterator iter;
	for (iter = states.begin(); iter != states.end(); iter++)
		xDotMat[iter->second] = xDotMap[iter->first];

	return xDotMat;
}

SXMatrix Ode::rk4Step( SXMatrix x0Vec, SXMatrix u0Vec, SXMatrix u1Vec, map<string,SX> & p, SX t0, SX t1)
{
	SX dt = t1-t0;
	SXMatrix k1 = dxVectorDt( x0Vec            ,             u0Vec, p, t0          );
	SXMatrix k2 = dxVectorDt( x0Vec + 0.5*dt*k1, 0.5*(u0Vec+u1Vec), p, t0 + 0.5*dt );
	SXMatrix k3 = dxVectorDt( x0Vec + 0.5*dt*k2, 0.5*(u0Vec+u1Vec), p, t0 + 0.5*dt );
	SXMatrix k4 = dxVectorDt( x0Vec +     dt*k3,             u1Vec, p, t0 +     dt );
	
	//	return x0Vec + dt*k1; // euler
	return x0Vec + dt*(k1 + 2*k2 + 2*k3 + k4)/6; // rk4
}
