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
int Ode::np(){ return  params.size(); }
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
	map<string, int>::const_iterator xIter, uIter, pIter;

	xIter = states.find(newName);
	uIter = actions.find(newName);
	pIter = params.find(newName);

	if (xIter != states.end() || uIter != actions.end() || pIter != params.end() ){
		cerr << "Error - new state/action/param \"" << newName << "\" is not unique" << endl;
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

int Ode::isParam(string paramName)
{
	map<string, int>::const_iterator pIter;

	pIter = params.find(paramName);

	if (pIter != params.end())
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

void Ode::addParam(string _newParam)
{
	assertUnlocked();
	assertUniqueName(_newParam);
	params[_newParam] = np() - 1;
}

SXMatrix Ode::dxVectorDt( SXMatrix x, SXMatrix u, SXMatrix p, SX t )
{
	map<string,int>::const_iterator iter;

	// create maps to wrap SXMatrices
	map<string,SX> xMap, uMap, pMap;
	for (iter = states.begin(); iter != states.end(); iter++)
		xMap[iter->first] = x.at(iter->second);
	for (iter = actions.begin(); iter != actions.end(); iter++)
		uMap[iter->first] = u.at(iter->second);
	for (iter = params.begin(); iter != params.end(); iter++)
		pMap[iter->first] = p.at(iter->second);

	// call dxdt
	map<string,SX> xDotMap;
	dxdt( xDotMap, xMap, uMap, pMap, t );
	
	// make output SXMatrix
	SXMatrix xDotMat = create_symbolic("an_xDotMat", nx());

	//	SXMatrix xDotMat(nx(), 1);
	for (iter = states.begin(); iter != states.end(); iter++)
		xDotMat[iter->second] = xDotMap[iter->first];

	return xDotMat;
}

SXMatrix Ode::rk4Step( SXMatrix x0Vec, SXMatrix u0Vec, SXMatrix u1Vec, SXMatrix pVec, SX t0, SX t1)
{
	SX dt = t1-t0;
	
	SXMatrix k1 = dxVectorDt( x0Vec            ,             u0Vec, pVec, t0          );
	SXMatrix k2 = dxVectorDt( x0Vec + 0.5*dt*k1, 0.5*(u0Vec+u1Vec), pVec, t0 + 0.5*dt );
	SXMatrix k3 = dxVectorDt( x0Vec + 0.5*dt*k2, 0.5*(u0Vec+u1Vec), pVec, t0 + 0.5*dt );
	SXMatrix k4 = dxVectorDt( x0Vec +     dt*k3,             u1Vec, pVec, t0 +     dt );
    
//	return x0Vec + dt*k1; // euler
	return x0Vec + dt*(k1 + 2*k2 + 2*k3 + k4)/6; // rk4
}


//    def setDxdt(self, dxdt):
//        if self.locked:
//            errStr = "Ode "+self.name+" has already been assigned to an Ocp and is in read-only mode"
//            raise ValueError(errStr)
//        self.dxdt = dxdt
//
//    # multiple rk4 steps on a single shooting interval
//    def rk4Steps(self, x0Vec, u0Vec, u1Vec, pVec, t0, t1):
//        N = 1
//        
//        x = x0Vec
//
//        for k in range(N):
//            dt = (t1 - t0)/np.double(N)
//            t0_ = t0+k*dt
//            t1_ = t0+(k+1)*dt
//
//            slider0 = (t0_-t0)/(t1 - t0)
//            slider1 = (t1_-t0)/(t1 - t0)
//
//            u0_ = u0Vec*(1.0 - slider0) + u1Vec*slider0
//            u1_ = u0Vec*(1.0 - slider1) + u1Vec*slider1
//
//            x = self.rk4Step(x, u0_, u1_, pVec, t0_, t1_)
//
//        return x
//
//
//    def runSim(self, time, x0, u, p):
//        N = len(time)
//        X = np.matrix(np.zeros([self._Nx(), N]))
//        X[:,0] = x0
//
//        for k in range(N-1):
//            u0 = u[:,k]
//            u1 = u[:,k+1]
//            dt = time[k+1] - time[k]
//            X[:,k+1] = self.rk4Step(X[:,k], u0, u1, p, time[k], time[k+1])
//
//        return X
