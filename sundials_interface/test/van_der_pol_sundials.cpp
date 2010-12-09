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

#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>

#include <fstream>

using namespace std;
using namespace OPTICON;

int main(){
  try{

  // Create a single-stage optimal control problem
  OCP_old ocp;

  // Time variable
  Matrix t = ocp.getTime();

  // Specify the time horizon
  Matrix t0 = 0,  tf = 10;
  ocp.setInitialTime(t0);
  ocp.setFinalTime(tf);

  // Differential states
  Matrix x1("x1"), x1dot("x1dot");
  Matrix x2("x2"), x2dot("x2dot");
  ocp.addState(x1,x1dot);
  ocp.addState(x2,x2dot);
  
  // Control
  Matrix u("u");
  ocp.addControl(u);

  // load the optimal solution (generated with ACADO)
  double u_opt_acado[] = {
    6.3641264739528092e-01,7.4999999999998601e-01,
    7.5000000000000000e-01,4.1252273198200651e-01,
    -1.8398623108660558e-01,-1.7948787122670185e-01,
    -6.5156242641681933e-02,-2.8194925877886162e-03,
    1.0489476754129789e-02,6.3093344776444004e-03,
    1.5442153523679725e-03,-3.1367783059813460e-04,
    -4.5769613236962433e-04,-1.9472485878029131e-04,
    -2.2975976127258114e-05,2.3294940136911033e-05,
    1.7055777534618476e-05,5.0469128292661197e-06,
    -3.2215265127224607e-07,-6.5237699317386640e-07};
  // round off:
  double u_min = -0.75;
  double u_max = 1;
  int n_disc = 1000;
  double du = (u_max-u_min)/(n_disc-1);
  for(int i=0; i<20; ++i){
    int i_disc = round((u_opt_acado[i]-u_min)/du);
    u_opt_acado[i] = u_min + du*i_disc;
  }

  Matrix t_opt = linspace(0.1,19.9,19);
  Matrix u_opt = Matrix(&u_opt_acado[0],20,1);
  Matrix u_guess = pw_const(t,t_opt,u_opt);
  ocp.guessSolution(u,u_guess);

  // Differential equation
  ocp.addEquation(x1dot, (1.0-x2*x2)*x1 - x2 + u);
  ocp.addEquation(x2dot, x1);

  // Initial conditions
  ocp.addInitialCondition(x1, 0);
  ocp.addInitialCondition(x2, 1);

  // Parametrize controls into 20 uniform intervals
  Matrix u_disc = parametrizeControls(ocp,u,10);

  // The desired output times
  int nt_out = 100; // number of outputs
  Matrix t_out = linspace(0,tf,nt_out);

  // LAGRANGE TERM
  int l_ind = ocp.addOutput(tf, Matrix(), x1*x1 + x2*x2 + u*u); 

  // Set output function
  int oind_x1 = ocp.addOutput(t_out, x1);
  int oind_x2 = ocp.addOutput(t_out, x2);
  int oind_u = ocp.addOutput(t_out, u);

  // Eliminate dependent variables from the functions
  eliminateDependent(ocp);

  // Print to screen
  cout << ocp;

 // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Set the linear solver
   integrator->setOption("linear_solver","dense");
//  integrator->setOption("linear_solver","band");
//  integrator->setOption("linear_solver","sparse");

  // Use exact jacobian
  integrator->setOption("exact_jacobian",false);

  // Upper and lower band-widths (only relevant when using a band linear solver)
//   integrator->setOption("mupper",1);
//   integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-8);

  // Initialize the integrator
  integrator->init();
  cout << "initialized" << endl;

  // Integrate once for visualization
  integrator->evaluate();
  cout << "solved" << endl;

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_vdp.txt");

  // Save results to file
  resfile << "t_out " << t_out << endl;
  resfile << "x1 " << integrator->output[oind_x1].data() << endl;
  resfile << "x2 " << integrator->output[oind_x2].data() << endl;
  resfile << "u " << integrator->output[oind_u].data() << endl;
  resfile << "lfun " << integrator->output[l_ind].data() << endl;
  cout << "lfun " << integrator->output[l_ind].data() << endl;

  // Close the results file
  resfile.close();


  return 0;
} catch (const char * str){
  cerr << str << endl;
  return 1;
}

}
