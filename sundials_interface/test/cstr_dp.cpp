/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include <iostream>
#include "casadi/sx/sx_matrix.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/mx/mx_constant.hpp"
#include "casadi/fx/external_function.hpp"
#include <limits>

using namespace std;
using namespace CasADi;

  // parameters
  double F0=100.0/1000.0/60.0;
  double c0=1000;
  double F=100.0/1000.0/60.0; 
  double T0 = 350;
  double r = 0.219;
  double k0 = 7.2e10/60.0;
  double EdivR = 8750;
  double U = 915.6;
  double rho = 1000;
  double Cp = 0.239*1000;
  double dH = -5e4;
  double V = 100;
  double c_init = 1000;
  double T_init = 350;

  double startTime=0.0;
  double finalTime=150;

  double c_ref = 500;
  double T_ref = 320;
  double Tc_ref = 300;
  double q_c = 1;
  double q_T = 1;
  double q_Tc = 1;      

  double lbu = 230;
  double ubu = 370;

// square
double sq(double x){
  return x*x;
}


// cost function
double L(const double* x, const double* u, const double* p){
/*  der(cost) = q_c*(c_ref-cstr.c)^2 + q_T*(T_ref-cstr.T)^2 + 
                  q_Tc*(Tc_ref-cstr.Tc)^2;*/
  double c = x[0];
  double T = x[1];
  double Tc = u[0];

  if(T>=350) return numeric_limits<double>::infinity();

  return q_c*sq(c_ref-c) + q_T*sq(T_ref-T) + 
                  q_Tc*sq(Tc_ref-Tc);
}

// ODE rhs
void f(const double* x, const double* u, const double* p, double* res){
//   der(c) = F0*(c0-c)/V-k0*c*exp(-EdivR/T);
//   der(T) = F0*(T0-T)/V-dH/(rho*Cp)*k0*c*exp(-EdivR/T)+2*U/(r*rho*Cp)*(Tc-T);
  double c = x[0];
  double T = x[1];
  double Tc = u[0];
  res[0] = F0*(c0-c)/V-k0*c*exp(-EdivR/T);
  res[1] = F0*(T0-T)/V-dH/(rho*Cp)*k0*c*exp(-EdivR/T)+2*U/(r*rho*Cp)*(Tc-T);
}

int main(){
  try{


  return 0;

  } catch (const char * str){
  cerr << str << endl;
  return 1;
}

  
}

