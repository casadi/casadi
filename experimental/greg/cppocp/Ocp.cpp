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

// Ocp.cpp
// Greg Horn
// Casadi 2011

#include "Ocp.hpp"
#include <cstdlib>

#include <core/function/sx_function.hpp>

using namespace std;
using namespace casadi;

Ocp::Ocp(){}

Ocp::~Ocp()
{
     map<string,MultipleShooting*>::const_iterator msIter;
     for (msIter = ms.begin(); msIter != ms.end(); msIter++){
          delete msIter->second;
     }
}

void Ocp::addNonlconIneq( SX gNew, string name)
{
     if (gNew.size2() != 1){
          cerr << "gNew.size2() != 1" << endl;
          throw 1;
     }

     g = vertcat(g, gNew);
     gSizes.push_back(gNew.size1());
     gLabels.push_back(name);

     for (int k=0; k<gNew.size1(); k++){
          gMin.push_back(-1.0e50);
          gMax.push_back(0.0);
          //    -numeric_limits<double>::infinity()
     }
}
void Ocp::addNonlconIneq( SX gNew )
{
     addNonlconIneq(gNew, "");
}

void Ocp::addNonlconEq( SX gNew, string name)
{
     if (gNew.size2() != 1){
          cerr << "gNew.size2() != 1" << endl;
          throw 1;
     }

     g = vertcat(g, gNew);
     gSizes.push_back(gNew.size1());
     gLabels.push_back(name);

     for (int k=0; k<gNew.size1(); k++){
          gMin.push_back(0);
          gMax.push_back(0);
     }
}
void Ocp::addNonlconEq( SX gNew )
{
     addNonlconEq(gNew, "");
}

void Ocp::assertUniqueName(string s)
{
     if (isMultipleShooting(s) || isParam(s)) {
          cerr << "Error - new MultipleShooting or param \"" << s << "\" is not unique" << endl;
          throw "1";
     }

     // make sure name is not state/action in an ode
     map<string, MultipleShooting*>::const_iterator msIter;
     for (msIter = ms.begin(); msIter != ms.end(); msIter++)
          (msIter->second)->ode.assertUniqueName(s);
}


/****************** multiple shooting stuff ****************/
int Ocp::isMultipleShooting(string msName)
{
     map<string, MultipleShooting*>::const_iterator msIter;

     msIter = ms.find(msName);

     if (msIter != ms.end())
          return 1;
     return 0;
}

MultipleShooting & Ocp::addMultipleShooting(string name, Ode & ode, SX t0, SX tf, int N)
{
     assertUniqueName(name);

     int numNew = ode.nxu()*N;
     if (designVariables.size1() == 0)
          designVariables = ssym(name, numNew);
     else
          designVariables = vertcat( designVariables, ssym(name, numNew) );

     for (int k=0; k<numNew; k++){
          lb.push_back(-1e50);
          ub.push_back(1e50);
          guess.push_back(0);
     }

     ms[name] = new MultipleShooting(name, ode, t0, tf, N, lb, ub, guess, designVariables, designVariables.size1() - ode.nxu()*N, params);

     // dynamics constraint
     for (int k=0; k<N-1; k++)
          addNonlconEq( ms[name]->getDynamicsConstraintError(k) );

     return *(ms[name]);
}


/****************** param stuff ******************/
// get param
SX & Ocp::getParam(string p)
{
     if (!isParam(p)){
          cerr << "Error in Ocp::getParam - \"" << p << "\" is not a known parameter\n";
          throw 1;
     }

     return params[p];
}

double Ocp::getParamSolution(string p)
{
     if (!isParam(p)){
          cerr << "Error in Ocp::getParamSolution - \"" << p << "\" is not a known parameter\n";
          throw 1;
     }
     return xopt.at(paramsIdx[p]);
}

// add a param
SX & Ocp::addParam(string _newParam)
{
     assertUniqueName(_newParam);

     int idx = designVariables.size1();

     SX newDv = ssym(_newParam, 1);
//      designVariables = vertcat( designVariables, ssym(_newParam, 1) );
     if (designVariables.size1() == 0)
          designVariables = newDv;
     else
          designVariables = vertcat( designVariables, newDv );

     guess.push_back(0);
     lb.push_back(-1e50);
     ub.push_back(1e50);

     paramsIdx[_newParam] = idx;
     params[_newParam] = designVariables.at(idx);

     return getParam(_newParam);

}

void Ocp::boundParam(string p, double lb_, double ub_)
{
     if (!isParam(p)){
          cerr << "Error in Ocp::boundParam - \"" << p << "\" is not a known parameter\n";
          throw 1;
     }
     int idx = paramsIdx[p];
     lb[idx] = lb_;
     ub[idx] = ub_;
     if (guess[idx] < lb[idx]) guess[idx] = lb[idx];
     if (guess[idx] > ub[idx]) guess[idx] = ub[idx];
}



void Ocp::setParamGuess(string p, double guess_)
{
     int idx = paramsIdx[p];
     guess[idx] = guess_;
     if (guess[idx] < lb[idx] ){
          cerr << "Requested initial guess " << p << " == " << guess[idx];
          cerr << " is less than lb[" << idx << "] == " << lb[idx] << endl;
          cerr << "Setting guess " << p << " = " << lb[idx] << endl;
          guess[idx] = lb[idx];
     }
     if (guess[idx] > ub[idx] ){
          cerr << "Requested initial guess " << p << " == " << guess[idx];
          cerr << " is greater than ub[" << idx << "] == " << ub[idx] << endl;
          cerr << "Setting guess " << p << " = " << ub[idx] << endl;
          guess[idx] = ub[idx];
     }
}

int Ocp::isParam(string paramName)
{
     map<string, int>::const_iterator pIter;

     pIter = paramsIdx.find(paramName);

     if (pIter != paramsIdx.end())
          return 1;
     return 0;
}



void Ocp::setStates( vector<DMatrix> & x)
{
     for (int timestep=0; timestep < x.size(); timestep++)
          setState( x.at(timestep), timestep);
}

void Ocp::setActions( vector<DMatrix> & u)
{
     for (int timestep=0; timestep < u.size(); timestep++)
          setAction( u.at(timestep), timestep);
}


void Ocp::setState( DMatrix x, int timestep )
{
     if (ms.size() != 1){
          cerr << "Error in Ocp::setState - ms.size != 1\n";
          throw 1;
     }

     MultipleShooting * firstMs = ms.begin()->second;

     map<string, int>::const_iterator iter;

     for (iter = firstMs->ode.states.begin(); iter != firstMs->ode.states.end(); iter++)
          xopt.at( firstMs->getIdx( iter->first, timestep ) ) = x.at( iter->second );
}


void Ocp::setAction( DMatrix u, int timestep )
{
     if (ms.size() != 1){
          cerr << "Error in Ocp::setAction - ms.size != 1\n";
          throw 1;
     }

     MultipleShooting * firstMs = ms.begin()->second;

     map<string, int>::const_iterator iter;

     for (iter = firstMs->ode.actions.begin(); iter != firstMs->ode.actions.end(); iter++)
          xopt.at( firstMs->getIdx( iter->first, timestep ) ) = u.at( iter->second );
}





// get states/actions (only for single-stage ocp now)
DMatrix Ocp::getStateSolution(int timestep)
{
     if (ms.size() != 1){
          cerr << "Error in Ocp::getStateSolution - ms.size != 1\n";
          throw 1;
     }

     MultipleShooting * firstMs = ms.begin()->second;

     return firstMs->getState( timestep, xopt );
}

DMatrix Ocp::getActionSolution(int timestep)
{
     if (ms.size() != 1){
          cerr << "Error in Ocp::getActionSolution - ms.size != 1\n";
          throw 1;
     }

     MultipleShooting * firstMs = ms.begin()->second;

     return firstMs->getAction( timestep, xopt );
}


void Ocp::writeSolution( const char * filename )
{
     ofstream f(filename);
     if (!f){
          cerr << "error opening " << filename << endl;
          exit(1);
     }

     writeSolution(f);

     f.close();
}

void Ocp::writeSolution( ostream & f )
{
     f.precision(15);
     for (int k=0; k<guess.size(); k++)
          f << xopt.at(k) << endl;
}

void Ocp::loadGuess( const char * filename )
{
     ifstream f(filename);

     if (!f){
          cerr << "error opening " << filename << endl;
          exit(1);
     }
     
     f.close();
}

void Ocp::loadGuess( istream & f )
{
#define MAXLINE 200
     char oneline[MAXLINE];
     int k=0;
     while(f.getline(oneline, MAXLINE)){

          double value;
          if (1 == sscanf( oneline, "%lf", &value)){
               if (k >= guess.size()){
		    cerr << "Too many design variables in guess file\n";
                    throw -1;
               }

               guess[k] = value;
               //cout << "file[" << k << "]: " << value << ", guess.size: " << guess.size() << endl;
               k++;
          }
     }
}

void Ocp::writeOctaveOutput( string name )
{
     writeOctaveOutput(name, "");
}

void Ocp::writeOctaveOutput( string name, string add )
{
     string filename(name + ".m");
     writeOctaveOutput(name, add, filename);
}

void Ocp::writeOctaveOutput( string name, string add, string filename )
{
     ofstream f((char*)filename.c_str());

     if (!f){
          cerr << "error opening " << filename << endl;
          exit(1);
     }
     
     writeOctaveOutput(name, add, f);

     f.close();
}

void Ocp::writeOctaveOutput( string name, string add, ostream & f )
{
     f.precision(12);

     f << "function opt = " << name << "()" << endl;
     f << "% autogenerated by Ocp::writeOctaveOutput\n\n";

     f << "opt = struct();\n\n";

     // print addition stuff first
     f << add;

     map<string, int>::const_iterator siIter;

     // params
     f << "% parameters\n";
     f << "opt.params = struct();\n";
     for (siIter = paramsIdx.begin(); siIter != paramsIdx.end(); siIter++){
          int idx = paramsIdx[ siIter->first ];
          f << "opt.params." << siIter->first << " = " << xopt.at(idx) << ";" << endl;
     }
     f << endl;

     // multiple shootings
     map<string, MultipleShooting*>::const_iterator msIter;
     f << "% multiple shooting stages\n";
     f << "opt.multipleShootingStages = struct();\n";
     for (msIter = ms.begin(); msIter != ms.end(); msIter++){
          f << "opt.multipleShootingStages." << (msIter->second)->name << " = ms_stage_" << (msIter->second)->name << "_out();\n";
     }
     f << endl << endl;

     /************** tie multiple shootings stages together **************/
     f << "%%%%%%%%%% tie stages together in concatonated output %%%%%%%%%%%\n";

     // time
     f << "opt.time = [";
     for (msIter = ms.begin(); msIter != ms.end(); msIter++)
          f << "opt.multipleShootingStages." << (msIter->second)->name << ".time,";
     f << "];\n\n";
        
     map<string, MultipleShooting*>::const_iterator firstMs;
     firstMs = ms.begin();

     // states
     f << "% states\n";
     f << "opt.states = struct();\n";
     if (ms.size() > 0){
          // loop through states in first multiple shooting
          for (siIter = firstMs->second->ode.states.begin(); siIter != firstMs->second->ode.states.end(); siIter++){

               // for each state, concatonate all multiple shootings
               f << "try\n";
               f << "    opt.states." << siIter->first << " = [";
               for (msIter = ms.begin(); msIter != ms.end(); msIter++){
                    f << "opt.multipleShootingStages." << msIter->first << ".states." << siIter->first << ",";
               }
               f << "];\n";
               f << "end\n";
          }
          f << endl << endl;
     }


     // actions
     f << "% actions\n";
     f << "opt.actions = struct();\n";
     if (ms.size() > 0){

          // loop through actions in first multiple shooting
          for (siIter = firstMs->second->ode.actions.begin(); siIter != firstMs->second->ode.actions.end(); siIter++){

               // for each state, concatonate all multiple shootings
               f << "try\n";
               f << "    opt.actions." << siIter->first << " = [";
               for (msIter = ms.begin(); msIter != ms.end(); msIter++){
                    f << "opt.multipleShootingStages." << msIter->first << ".actions." << siIter->first << ",";
               }
               f << "];\n";
               f << "end\n";
          }
          f << endl << endl;
     }

     // outputs
     f << "% outputs\n";
     f << "opt.outputs = struct();\n";
     if (ms.size() > 0){

          // loop through outputs in first multiple shooting
          for (siIter = firstMs->second->ode.outputs.begin(); siIter != firstMs->second->ode.outputs.end(); siIter++){

               // for each state, concatonate all multiple shootings
               f << "try\n";
               f << "    opt.outputs." << siIter->first << " = [";
               for (msIter = ms.begin(); msIter != ms.end(); msIter++){
                    f << "opt.multipleShootingStages." << msIter->first << ".outputs." << siIter->first << ",";
               }
               f << "];\n";
               f << "end\n";
          }
          f << endl << endl;
     }

        
     /************** write all the multiple shooting outputs *********/
     f << "%%%%%%%%%% each stage's individual output %%%%%%%%%%%\n";
     for (msIter = ms.begin(); msIter != ms.end(); msIter++)
          (msIter->second)->writeOctaveOutput( f, xopt );
}
