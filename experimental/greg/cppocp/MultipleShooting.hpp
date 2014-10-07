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

// MultipleShooting.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <core/sx/sx_tools.hpp>

class MultipleShooting
{
public:
     MultipleShooting(std::string name_,
		      Ode & _ode,
		      casadi::SX t0_,
		      casadi::SX tf_,
		      int N_,
		      std::vector<double>&lb_,
		      std::vector<double>&ub_,
		      std::vector<double>&guess_,
		      casadi::SX & dv_,
		      int idx0_,
		      std::map<std::string,casadi::SX>&params_);
     ~MultipleShooting();

     int getBigN(void);

     void boundStateAction(std::string xu, double lb_, double ub_, int timeStep);
     void boundStateAction(std::string xu, double lb_, double ub_);
     void setStateActionGuess(std::string xu, double guess_, int timeStep);
     void setStateActionGuess(std::string xu, double guess_);

     void writeOctaveOutput( std::ostream & f, std::vector<double> & xopt);

     casadi::SX getState(std::string x, int timeStep);
     casadi::SX getAction(std::string u, int timeStep);
     casadi::SX getOutput(std::string o, int timeStep);
     casadi::SX getOutput(std::string o);

     casadi::SX getStateMat(int timeStep);
     casadi::SX getActionMat(int timeStep);

     casadi::DMatrix getState(int timeStep, std::vector<double> & xopt);
     casadi::DMatrix getAction(int timeStep, std::vector<double> & xopt);


     int N;
	
     casadi::SX getDynamicsConstraintError(int timeStep);

     Ode & ode;
     int getIdx(std::string xu, int timeStep);
  
     std::string name;

private:
     casadi::SX & dv;
     std::vector<double>&lb;
     std::vector<double>&ub;
     std::vector<double>&guess;
     std::map<std::string,casadi::SX> & params;
     int idx0;

     casadi::SX t0;
     casadi::SX tf;

};
