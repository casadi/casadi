/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
#include <ctime>
#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include <casadi/fx/jacobian.hpp>


using namespace CasADi;
using namespace std;

// Define the source code of the function which will be used both for code generation and evaluation
#define SOURCE_CODE(NAME)                        \
                                                 \
void NAME(){                                     \
  /* This is the hello world program! */         \
  cout << "hello world!" << endl;                \
}                                                \

// Create the function
SOURCE_CODE(hello)

// Convert a macro to a string using the double expansion trick
#define TOSTRING1(x)  #x
#define TOSTRING(x)  TOSTRING1(x)

void print_code(const std::string& code){
  for(string::const_iterator it=code.begin(); it!=code.end(); ++it){
    cout << *it;
    
    // insert line breaks after {, } and ; for readability
    if(*it=='{' || *it=='}' || *it==';'){
      cout << endl;
      
      // Skip space in the beginning of the next line
      if(it+1 != code.end() && *(it+1)==' ')
        it++;
    }
  }
}

int main(){
  
  // Run the program
  hello();
  
  // Print the code
  cout << "The code was : " << endl;
  print_code(TOSTRING(SOURCE_CODE(hello_modified)));
    
  cout << "program started" << endl;
      
  // Dimensions
  int nu = 1000;  // Number of control segments
  int nj = 1000; // 10000;  // // Number of integration steps per control segment

  // optimization variable
  vector<SX> u = create_symbolic("u",nu); // control

  SX s_0 = 0; // initial position
  SX v_0 = 0; // initial speed
  SX m_0 = 1; // initial mass
  
  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Trajectory
//  vector<SX> v_traj(nu);

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for(int k=0; k<nu; ++k){
    for(int j=0; j<nj; ++j){
      s += dt*v;
      v += dt / m * (u[k]- alpha * v*v);
      m += -dt * beta*u[k]*u[k];
    }
  //  v_traj[k] = v;
  }

  // Objective function
  SX f = 0;
  for(int i=0; i<u.size(); ++i)
    f += u[i]*u[i];
    
  // Terminal constraints
  vector<SX> g(2);
  g[0] = s;
  g[1] = v;
/*  for(vector<SX>::const_iterator it=v_traj.begin(); it!=v_traj.end(); ++it)
    g << *it;*/
  
  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  SXFunction gfcn(u,g); // constraint
  gfcn.setOption("ad_mode","reverse");
  gfcn.setOption("symbolic_jacobian",false);
  
  // Allocate an NLP solver
  IpoptSolver solver(ffcn,gfcn);
//  IpoptSolver solver(ffcn,gfcn,FX(),Jacobian(gfcn));

  // Set options
  solver.setOption("tol",1e-6);
  solver.setOption("hessian_approximation","limited-memory");

  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> umin(nu), umax(nu), usol(nu);
  for(int i=0; i<nu; ++i){
    umin[i] = -10;
    umax[i] =  10;
    usol[i] = 0.4;
  }
  solver.setInput(umin,NLP_LBX);
  solver.setInput(umax,NLP_UBX);
  solver.setInput(usol,NLP_X_INIT);
  
  // Bounds on g
  vector<double> gmin(2,-numeric_limits<double>::infinity()), gmax(2,1.1);
  gmin[0] = gmax[0] = 10;
  gmin[1] = gmax[1] =  0;
  solver.setInput(gmin,NLP_LBG);
  solver.setInput(gmax,NLP_UBG);

  // Solve the problem
  solver.solve();

  // Print the optimal cost
  double cost;
  solver.getOutput(cost,NLP_COST);
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  solver.getOutput(usol,NLP_X_OPT);
  cout << "optimal solution: " << usol << endl;

  return 0;
}

