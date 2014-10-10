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

// Ode.cpp
// Greg Horn
// Casadi 2011

#include <iostream>

#include "Ode.hpp"

using namespace std;
using namespace casadi;

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
	  cerr << "Error - Ode \"" << name << "\" has already been initialized" << endl;
	  throw 1;
     }
}

void Ode::init()
{
     if (locked)
	  return;
     locked = 1;
}

void Ode::setupIntegrators()
{
     map<string,double> dummyParams;
     setupIntegrators(dummyParams);
}

void Ode::setupIntegrators(map<string,double> & params)
{
     /*********** set up rk4Step and eulerStep SXFunctions **********/
     // inputs
     SX xk = ssym("xk", nx());
     SX uk = ssym("uk", nu());
     SX t0 = ssym("t0",    1);
     SX dt = ssym("dt",    1);

     vector<SX> stepInputs(NUM_ODE_STEP_INPUTS);
     stepInputs.at(IDX_ODE_STEP_STATE)  = xk;
     stepInputs.at(IDX_ODE_STEP_ACTION) = uk;
     stepInputs.at(IDX_ODE_STEP_T0)     = t0;
     stepInputs.at(IDX_ODE_STEP_DT)     = dt;

     // "cast" map<string,double> to map<string,SX>
     map<string,SX> sxParams;
     map<string,double>::const_iterator iter;
     for (iter = params.begin(); iter != params.end(); iter++){
       sxParams[iter->first] = SX(iter->second);
     }

     // call fcns
     SX xNextRK4   = rk4Step(   xk, uk, uk, sxParams, t0.at(0), dt.at(0) );
     SX xNextEuler = eulerStep( xk, uk,     sxParams, t0.at(0), dt.at(0) );

     // outputs
     rk4StepFcn   = SXFunction( stepInputs, xNextRK4 );
     eulerStepFcn = SXFunction( stepInputs, xNextEuler );
     rk4StepFcn.init();
     eulerStepFcn.init();
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

map<string,SX> Ode::getStateMap( SX & x )
{
     map<string,int>::const_iterator iter;
     map<string,SX> xMap;
     for (iter = states.begin(); iter != states.end(); iter++)
	  xMap[iter->first] = x.at(iter->second);

     return xMap;
}

map<string,SX> Ode::getActionMap( SX & u )
{
     map<string,int>::const_iterator iter;
     map<string,SX> uMap;
     for (iter = actions.begin(); iter != actions.end(); iter++)
	  uMap[iter->first] = u.at(iter->second);

     return uMap;
}

map<string,SX> Ode::getOutputFromDxdt( SX x, SX u, map<string,SX> & p, SX t )
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


SX Ode::dxVectorDt( SX x, SX u, map<string,SX> & p, SX t )
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

     // make output SX
     SX xDotMat = ssym("an_xDotMat", nx());
     map<string,int>::const_iterator iter;
     for (iter = states.begin(); iter != states.end(); iter++)
	  xDotMat[iter->second] = xDotMap[iter->first];

     return xDotMat;
}

SX Ode::rk4Step( SX x0Vec, SX u0Vec, SX u1Vec, map<string,SX> & p, SX t0, SX dt)
{
     SX k1 = dxVectorDt( x0Vec            ,             u0Vec, p, t0          );
     SX k2 = dxVectorDt( x0Vec + 0.5*dt*k1, 0.5*(u0Vec+u1Vec), p, t0 + 0.5*dt );
     SX k3 = dxVectorDt( x0Vec + 0.5*dt*k2, 0.5*(u0Vec+u1Vec), p, t0 + 0.5*dt );
     SX k4 = dxVectorDt( x0Vec +     dt*k3,             u1Vec, p, t0 +     dt );
	
     return x0Vec + dt*(k1 + 2*k2 + 2*k3 + k4)/6; // rk4
}

SX Ode::eulerStep( SX x0Vec, SX u0Vec, map<string,SX> & p, SX t0, SX dt)
{
     return x0Vec + dt*dxVectorDt( x0Vec, u0Vec, p, t0 );
}

SX Ode::simpsonsRuleError( SX x0Vec, SX x1Vec, SX u0Vec, SX u1Vec, map<string,SX> & p, SX t0, SX dt)
{
     SX f0 = dxVectorDt( x0Vec, u0Vec, p, t0          );
     SX f1 = dxVectorDt( x1Vec, u1Vec, p, t0 +     dt );

     SX um = 0.5*( u0Vec + u1Vec );
     SX xm = 0.5*( x0Vec + x1Vec ) - 0.125*(f1-f0)*dt;

     SX fm = dxVectorDt( xm, um, p, t0 + 0.5*dt );

     return x1Vec - x0Vec - dt/6.0*(f0 + 4*fm + f1);
}


DMatrix Ode::rk4Step( DMatrix & xk, DMatrix & uk, double t0_, double dt_)
{
     if (!locked)
	  init();

     DMatrix t0(t0_);
     DMatrix dt(dt_);

     // inputs
     rk4StepFcn.setInput( xk, IDX_ODE_STEP_STATE  );
     rk4StepFcn.setInput( uk, IDX_ODE_STEP_ACTION );
     rk4StepFcn.setInput( t0, IDX_ODE_STEP_T0     );
     rk4StepFcn.setInput( dt, IDX_ODE_STEP_DT     );

     // evaluate
     rk4StepFcn.evaluate();
	
     // output
     DMatrix out( nx(), 1, 0.0 );
     rk4StepFcn.getOutput(out);
     return out;
}

DMatrix Ode::eulerStep( DMatrix & xk, DMatrix & uk, double t0_, double dt_)
{
     if (!locked)
	  init();

     DMatrix t0(t0_);
     DMatrix dt(dt_);

     // inputs
     eulerStepFcn.setInput( xk, IDX_ODE_STEP_STATE  );
     eulerStepFcn.setInput( uk, IDX_ODE_STEP_ACTION );
     eulerStepFcn.setInput( t0, IDX_ODE_STEP_T0     );
     eulerStepFcn.setInput( dt, IDX_ODE_STEP_DT     );

     // evaluate
     eulerStepFcn.evaluate();

     // output
     DMatrix out( nx(), 1, 0.0 );
     eulerStepFcn.getOutput(out);
     return out;
}
