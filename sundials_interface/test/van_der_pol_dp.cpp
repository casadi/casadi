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
#include <casadi/stl_vector_tools.hpp>
#include <casadi/mx/mx_constant.hpp>
#include <casadi/fx/external_function.hpp>

using namespace std;
using namespace CasADi;

// "infinity", numeric_limits<double>::infinity() gives not-a-number when multiplied with 0!
double infty = numeric_limits<double>::infinity();
//double infty = 1e100;


        // Parameters
    double t0 = 0;
    double tf = 10;

    double p1 = 1;             // Parameter 1
    double p2 = 1;             // Parameter 2
    double p3 = 2;             // Parameter 3

// state space discretization
const int NX = 2;
double xmin[2] = {-1,-1};
double xmax[2] = {1,1};
const int nx[2] = {50, 50};
double delta_x[2] = {2.0/49,2.0/49};

// control space discretization
const int NU = 1;
double umin[1] = {-1.0};
double umax[1] = {0.75};
const int nu[1] = {100};
double delta_u[1] = {1.75/99};

// time discretization
const int NT = 100;
const double DT = tf/NT;

// the integer corresponding to a value of the state
int xToInt(double xval, int ind){
  return round((xval-xmin[ind])/delta_x[ind]);
}

int xError(double xval, int ind){
  double a = (xval-xmin[ind])/delta_x[ind];
  return a - round(a);
}
 

// the value of the state corresponding to an integer
double xFromInt(int i, int ind){
  return xmin[ind] + i*delta_x[ind];
}

// the value of a control corresponding to an integer
double uFromInt(int i, int ind){
  return umin[ind] + i*delta_u[ind];
}


// cost function
double L(const double* x, const double* u, const double* p){
//     der(cost) = exp(p3 * 1/*time*/) * (x1^2 + x2^2 + u^2);
  return exp(p3 * 1/*time*/) * (x[0]*x[0] + x[1]*x[1] + u[0]*u[0]);
}

// stage transition
void F(const double* x, const double* u, const double* p, double* res){
//     der(x1) = (1 - x2^2) * x1 - x2 + u;
//     der(x2) = p1 * x1;

  res[0] = x[0] + DT*((1 - x[1]*x[1]) * x[0] - x[1] + u[0]);
  res[1] = x[1] + DT*(p1 * x[0]);
}

int main(){
  try{

    
    // Initial conditions
    double x1_0 = 0;
    double x2_0 = 1;

  // discretize the cost to arrive
  double J[nx[0]][nx[1]];
  double Jnew[nx[0]][nx[1]];
  
  // Discretize the optimal control
  int u_stack[NT][nx[0]][nx[1]];
  int u_stack0[NT][nx[0]][nx[1]];
  int u_stack1[NT][nx[0]][nx[1]];
  
  // Set to infinity everywhere
  for(int i0=0; i0<nx[0]; ++i0)
    for(int i1=0; i1<nx[1]; ++i1){
      J[i0][i1] = numeric_limits<double>::infinity();
  }

  // Set cost at initial value to zero
  J[xToInt(x1_0,0)][xToInt(x2_0,1)] = 0;


  // state vector
  double x[2];
  
  // state derivative
  double xnext[2];

  // control vector
  double u[1];

  // Loop over all time points
  for(int k=0; k<NT; ++k){

  // Set the new cost to arrive to infinity
  for(int i0=0; i0<nx[0]; ++i0)
    for(int i1=0; i1<nx[1]; ++i1){
      Jnew[i0][i1] = numeric_limits<double>::infinity();
  }

  // Loop over all spatial points
  for(int i0=0; i0<nx[0]; ++i0){
    x[0] = xFromInt(i0,0);
    for(int i1=0; i1<nx[1]; ++i1){
      x[1] = xFromInt(i1,1);
      
      // Get current cost
      double Ji = J[i0][i1];

     // Loop over all possible controls
      for(int j0=0; j0<nu[0]; ++j0){
        u[0] = uFromInt(j0,0);

        // Evaluate the cost function
        double Lij = L(x, u, 0);
        
        // Take a step
        F(x,u,0,xnext);

        // Get the integer corresponding to the new x
        int r0 = xToInt(xnext[0],0);
        if(r0<0 || r0>=nx[0]) continue;
        int r1 = xToInt(xnext[1],1);
        if(r1<0 || r1>=nx[1]) continue;

        // Add a regularization term
        double err0 = xError(r0,0);
        double err1 = xError(r1,1);
        Lij += 10*(err0*err0 + err1*err1);

        // save if cost is smaller
        if(Jnew[r0][r1] > Lij + Ji){
          Jnew[r0][r1] = Lij + Ji;
          u_stack[k][r0][r1] = j0;
          u_stack0[k][r0][r1] = i0;
          u_stack1[k][r0][r1] = i1;
        }
      } // for j0 = ...
    } // for i1 = ...
  } // for i0 = ...

  // Copy the cost to arrive
  for(int i0=0; i0<nx[0]; ++i0)
    for(int i1=0; i1<nx[1]; ++i1){
      J[i0][i1] = Jnew[i0][i1];
    }
  } // for k = ...

  // print J
  cout.precision(1);
  cout.setf(ios::fixed);
  for(int i0=0; i0<nx[0]; ++i0){
    for(int i1=0; i1<nx[1]; ++i1){
      cout << J[i0][i1] << ",";
    }
    cout << endl;
  }

  // Find the optimal cost
  double Jmin = infty;
  int i0_opt, i1_opt;
  for(int i0=0; i0<nx[0]; ++i0)
    for(int i1=0; i1<nx[1]; ++i1){
      if(Jmin > J[i0][i1]){
        Jmin = J[i0][i1];
        i0_opt = i0;
        i1_opt = i1;
      }
  }
  cout << "optimal cost: " << Jmin*DT << endl;
  cout << "terminal value of x: " << xFromInt(i0_opt,0)<< "," << xFromInt(i1_opt,1) << endl;

  // get the optimal control
  int u_opt[NT];
  for(int k=NT-1; k>=0; --k){
    u_opt[k] = u_stack[NT-1][i0_opt][i1_opt];
    int i0_opt_old = i0_opt;
    i0_opt = u_stack0[k][i0_opt][i1_opt];
    i1_opt = u_stack1[k][i0_opt_old][i1_opt];
  }

  cout << "optimal control: {" << uFromInt(u_opt[0],0);
  for(int k=1; k<NT; ++k)
    cout << "," << uFromInt(u_opt[k],0);
  cout << "}" << endl;

  return 0;

  } catch (const char * str){
  cerr << str << endl;
  return 1;
}

  
}

