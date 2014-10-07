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

// Ocp.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include <iostream>
#include <string>
#include <map>

#include "Ode.hpp"
#include "MultipleShooting.hpp"

#include <core/sx/sx_tools.hpp>
#include <core/std_vector_tools.hpp>

class Ocp
{
public:
     Ocp(void);
     ~Ocp(void);

     void addNonlconIneq(casadi::SX gNew);
     void addNonlconIneq(casadi::SX gNew, std::string name);      
     void addNonlconEq(casadi::SX gNew);
     void addNonlconEq(casadi::SX gNew, std::string name);

     casadi::SX designVariables;
     std::vector<double>lb;
     std::vector<double>ub;
     std::vector<double>guess;

     casadi::SX objFun;
     casadi::SX g;
     std::vector<int> gSizes;
     std::vector<std::string> gLabels;
     std::vector<double> gMin;
     std::vector<double> gMax;

     MultipleShooting & addMultipleShooting(std::string name, Ode & ode, casadi::SX t0, casadi::SX tf, int N);

     casadi::SX & addParam(std::string _newParam);
     casadi::SX & getParam(std::string p);
     double getParamSolution(std::string p);
     void boundParam(std::string p, double lb_, double ub_);
     void setParamGuess(std::string p, double guess_);
        
     casadi::DMatrix getStateSolution(int timestep);
     casadi::DMatrix getActionSolution(int timestep);

     void setState(  casadi::DMatrix x, int timestep );
     void setAction( casadi::DMatrix u, int timestep );
     void setStates(  std::vector<casadi::DMatrix> & x);
     void setActions( std::vector<casadi::DMatrix> & u);

     void writeSolution( std::ostream & f );
     void writeSolution( const char * filename );
     void writeOctaveOutput( std::string name );
     void writeOctaveOutput( std::string name, std::string add );
     void writeOctaveOutput( std::string name, std::string add, std::string filename );
     void writeOctaveOutput( std::string name, std::string add, std::ostream & f );
     void loadGuess( std::istream & f );
     void loadGuess( const char * filename );

     // multiple shooting instances
     std::map<std::string,MultipleShooting*> ms;

     std::vector<double> xopt;

private:

     void assertUniqueName(std::string s);
     int isMultipleShooting(std::string msName);



     // params
     std::map<std::string,int> paramsIdx;
     std::map<std::string,casadi::SX> params;
     int np(void);
     int isParam(std::string paramName);


     casadi::SX getParamMat(void);

     int getParamIdx(std::string p);



protected:
};
