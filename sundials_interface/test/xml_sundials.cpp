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
#include <fstream>
#include <string>

#include <integrator/sundials_interface/cvodes_integrator.hpp>
#include <ocp/xml_interface/fmi_parser.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <casadi/expression_tools.hpp>

using namespace std;
using namespace OPTICON;

int main(int argc, char *argv[]){

try{

  // Read the file name
  if(argc!=2) throw "The function must have exactly one argument";
  char *modelfile = argv[1];

  // Allocate a parser
  FMIParser parser;
  
  // Load the xml model
  parser.loadFile(modelfile);

  // Dump representation to screen
  parser.dump();

  // Create an optimal control problem
  OCP ocp;

  // Read from the file
  parser >> ocp;
  cout << ocp;

  // Make the ocp explicit
  makeExplicit(ocp);

  // Sort the initial conditions and make then explicit
  sortInitialConditions(ocp,true);

  // Print the ocp to screen
  cout << ocp;

  // Get the integration interval
  Matrix t0 = ocp.getInitialTime();
  Matrix tf = ocp.getFinalTime();

  // The desired output times
  int nt_out = 100; // number of outputs
  Matrix t_out = linspace(t0,tf,nt_out);

  // Get the variables
  Matrix t = ocp.getTime();
  Matrix u = ocp.variable("u");
  Matrix x1 = ocp.variable("x1");
  Matrix x2 = ocp.variable("x2");
  Matrix cost = ocp.variable("cost");
  Matrix p1 = ocp.variable("p1");
  Matrix p2 = ocp.variable("p2");
  Matrix p3 = ocp.variable("p3");

  // Guess the parameters (should read from file)
  ocp.guessSolution(p1, 1);
  ocp.guessSolution(p2, 1);
  ocp.guessSolution(p3, 2);

  // Set output function
  int outind_u = ocp.addOutput(u, t_out);
  int outind_x1 = ocp.addOutput(x1, t_out);
  int outind_x2 = ocp.addOutput(x2, t_out);

  // Give a guess for the control trajectory
  Matrix u0_t, u0_u; // piecewise linear approximation through these points
  u0_t <<  0.0 <<  0.2 <<    2 <<  4   << 6  <<  20;
  u0_u << -0.5 <<  0.7 <<  0.7 << -0.1 << 0  <<   0; 
  Matrix u0 = pw_lin(t, u0_t, u0_u);
//  ocp.guessSolution(u, u0);
  ocp.guessSolution(u, sin(t));
 
  // Discretize controls into 20 uniform intervals
  Matrix u_disc = parametrizeControls(ocp,u,100);
  cout << u_disc << endl;

  // Eliminate dependent variables from the functions (WILL BE REMOVED)
  eliminateDependent(ocp);

   // Allocate an integrator
  CvodesIntegrator integrator(ocp);

  // Set the linear solver
   integrator->setOption("linear_solver","dense");
  //  integrator->setOption("linear_solver","band");
  //  integrator->setOption("linear_solver","sparse");

  // Use the exact Jacobian (or an approximation using finite differences)
  integrator->setOption("exact_jacobian",true);

  // Upper and lower band-widths (only relevant when using a band linear solver)
  //   integrator->setOption("mupper",1);
  //   integrator->setOption("mlower",1);

  // set tolerances 
  integrator->setOption("reltol",1e-6);
  integrator->setOption("abstol",1e-8);

  // Initialize the integrator
  integrator->init();

  // Integrate once for visualization
  integrator->setOption("use_events",false);
  integrator->setOption("use_output",true);
  integrator->setOption("use_fsens",false);
  integrator->setOption("use_asens",false);
  integrator->setOption("use_quad",false);
  integrator->evaluateFwd();

  // Create a file for saving the results
  ofstream resfile;
  resfile.open ("results_vdp.txt");

  // Save results to file
  resfile << "t_out " << t_out << endl;
  resfile << "u " << integrator->getOutput(outind_u) << endl;
  resfile << "x1 " << integrator->getOutput(outind_x1) << endl;
  resfile << "x2 " << integrator->getOutput(outind_x2) << endl;

  
  // Close the results file
  cout << "wrote results to results_vdp.txt" << endl;
  resfile.close();

  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
