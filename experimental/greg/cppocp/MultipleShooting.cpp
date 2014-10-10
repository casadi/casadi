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

// MultipleShooting.cpp
// Greg Horn
// Casadi 2011

#include "MultipleShooting.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <core/function/sx_function.hpp>

using namespace casadi;
using namespace std;

MultipleShooting::MultipleShooting(string name_,
				   Ode & _ode,
				   SX t0_,
				   SX tf_,
				   int N_,
				   vector<double>&lb_,
				   vector<double>&ub_,
				   vector<double>&guess_,
				   SX & dv_,
				   int idx0_,
				   map<string,SX> & params_) : ode(_ode), dv(dv_), lb(lb_), ub(ub_),
							       guess(guess_), params(params_)
{
     name = name_;
     t0 = t0_;
     tf = tf_;
     N = N_;
     idx0 = idx0_;

     ode.init();

     // make sure no ocp parameter names conflict with ode states/actions
     map<string, SX>::const_iterator iter;
     for (iter = params.begin(); iter != params.end(); iter++)
	  ode.assertUniqueName(iter->first);
}

MultipleShooting::~MultipleShooting() {}

SX MultipleShooting::getOutput(string o)
{
     SX ret = ssym(o, N);
     for (int k=0; k<N; k++)
	  ret.at(k) = getOutput(o, k);
	
     return ret;
}

SX MultipleShooting::getOutput(string o, int timeStep)
{
     if (timeStep > N-1){
	  cerr << "In SX MultipleShooting::getOutput(string o, int timeStep),\n";
	  cerr << "timeStep " << timeStep << " > N-1\n";
	  exit(1);
     }

     // dynamics constraint
     SX dt = (tf - t0)/(N - 1);

     SX xk = getStateMat(timeStep);
     SX uk = getActionMat(timeStep);
	
     SX tk = t0 + timeStep*dt;
	
     map<string,SX> output = ode.getOutputFromDxdt( xk, uk, params, tk );

     // make sure output exists
     map<string, SX>::const_iterator oIter;
     oIter = output.find(o);
     if (oIter == output.end()){
	  cerr << "Error - SX MultipleShooting::getOutput(string o, int timeStep) could not find output \"" << o << "\"\n";
	  throw 1;
     }
	
     return output[o];
}

SX MultipleShooting::getDynamicsConstraintError(int timeStep)
{
     if (timeStep > N-2){
	  cerr << "In SX MultipleShooting::getDynamicsConstraintError(int timeStep),\n";
	  cerr << "timeStep: " << timeStep << " > N-2   (N==" << N << ")\n";
	  exit(1);
     }

     // dynamics constraint
     SX dt = (tf - t0)/(N - 1);

     SX x0 = getStateMat(timeStep);
     SX x1 = getStateMat(timeStep + 1);
	
     SX u0 = getActionMat(timeStep);
     SX u1 = getActionMat(timeStep + 1);
	
     SX tk0 = t0 + timeStep*dt;
	
     //SX xErr = x1 - ode.rk4Step( x0, u0, u1, params, tk0, dt);
     //SX xErr = x1 - ode.eulerStep( x0, u0, params, tk0, dt);
     SX xErr = ode.simpsonsRuleError( x0, x1, u0, u1, params, tk0, dt);

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


SX MultipleShooting::getStateMat(int timeStep)
{
     SX ret = ssym("a_state", ode.nx());
     // SX ret(ode.nx(), 1);
     map<string,int>::const_iterator xIter;
     for (xIter = ode.states.begin(); xIter != ode.states.end(); xIter++)
	  ret[xIter->second] = getState(xIter->first, timeStep);

     return ret;
}

DMatrix MultipleShooting::getState(int timeStep, vector<double> & xopt)
{
     DMatrix ret(ode.nx(), 1, 0.0);
     map<string,int>::const_iterator xIter;
     for (xIter = ode.states.begin(); xIter != ode.states.end(); xIter++)
	  ret.at(xIter->second) = xopt.at(getIdx( xIter->first, timeStep ));

     return ret;
}

DMatrix MultipleShooting::getAction(int timeStep, vector<double> & xopt)
{
     DMatrix ret(ode.nu(), 1, 0.0);
     map<string,int>::const_iterator uIter;
     for (uIter = ode.actions.begin(); uIter != ode.actions.end(); uIter++)
	  ret.at(uIter->second) = xopt.at(getIdx( uIter->first, timeStep ));

     return ret;
}

SX MultipleShooting::getActionMat(int timeStep)
{
     SX ret = ssym("an_action", ode.nu());
     //	SX ret(ode.nu(),1);
     map<string,int>::const_iterator uIter;
     for (uIter = ode.actions.begin(); uIter != ode.actions.end(); uIter++)
	  ret[uIter->second] = getAction(uIter->first, timeStep);

     return ret;
}

void MultipleShooting::writeOctaveOutput( ostream & f, vector<double> & xopt )
{
     f.precision(12);

     f << "function multipleShooting = ms_stage_" << name << "_out()" << endl;

     map<string, int>::const_iterator iter;

     // states
     f << "% states\n";
     f << "multipleShooting = struct();\n";
     f << "multipleShooting.states = struct();\n";
     for (iter = ode.states.begin(); iter != ode.states.end(); iter++){
	  f << "multipleShooting.states." << iter->first << " = [" ;
	  for (int k=0; k<N; k++){
	       int idx = getIdx( iter->first, k );
	       if (k < N-1)
		    f << xopt.at(idx) << ", ";
	       else
		    f << xopt.at(idx) << "];" << endl;
	  }
     }
     f << endl;


     // actions
     f << "% actions\n";
     f << "multipleShooting.actions = struct();\n";
     for (iter = ode.actions.begin(); iter != ode.actions.end(); iter++){
	  f << "multipleShooting.actions." << iter->first << " = [" ;
	  for (int k=0; k<N; k++){
	       int idx = getIdx( iter->first, k );
	       if (k < N-1)
		    f << xopt.at(idx) << ", ";
	       else
		    f << xopt.at(idx) << "];" << endl;
	  }
     }
     f << endl;


     // outputs
     f << "% outputs\n";
     f << "multipleShooting.outputs = struct();\n";
     for (iter = ode.outputs.begin(); iter != ode.outputs.end(); iter++){

	  SX outSX = getOutput(iter->first);
	  SXFunction outputFcn(dv, outSX);
	  outputFcn.init();
	  outputFcn.setInput(xopt);
	  outputFcn.evaluate();
	  vector<double>outDouble(N);
	  outputFcn.getOutput(outDouble);

	  f << "multipleShooting.outputs." << iter->first << " = [";
	  for (int k=0; k<N; k++){
	       if (k < N - 1)
		    f << outDouble[k] << ", ";
	       else
		    f << outDouble[k] << "];" << endl;
	  }
     }
     f << endl;


     // start/end times
     f << "% time\n";
     SXFunction timeFcn( dv, vertcat( SX(t0), SX(tf) ) );
     timeFcn.init();
     timeFcn.setInput( xopt );
     timeFcn.evaluate();
     double timeNum[2];
     timeFcn.getOutput( timeNum );

     f << "multipleShooting.time = linspace(" << timeNum[0] << ", " << timeNum[1] << ", " << N << ");\n\n";
}
