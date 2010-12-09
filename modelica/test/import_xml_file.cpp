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

#include <acado_modelling/modelica/fmi_parser.hpp>
#include <acado_modelling/ocp_tools.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <acado_interface/acado_interface.hpp>

using namespace std;
using namespace OPTICON;


int main(int argc, char *argv[]){

try{

  // Read the file name
  if(argc!=2) throw "The function must have exactly one argument";
  char *modelfile = argv[1];

  // Allocate a parser and load the xml model
  FMIParser parser(modelfile);
  
  // Dump representation to screen
  cout << "**********************" << endl;
  cout << "    XML representation of \"" << modelfile << "\":" << endl;
  cout << "**********************" << endl;
  parser.print();
  
  // Print the ocp to screen
  cout << endl;
  cout << "**********************"     << endl;
  cout << "    OCP: "                  << endl;
  cout << "**********************"     << endl;

  // Obtain the symbolic representation of the OCP
  AcadoOCP ocp = parser.parse();

  // Sort the variables according to type
  ocp.sortVariables();
 
  // Print to screen
  ocp.print();

  cout << "**********************"     << endl;
  cout << "    Play around with the variables: " << endl;
  cout << "**********************"     << endl;
  SX t = ocp.vars("time");
  cout << "t = " << t << endl;
  
  Variable x1 = ocp.vars("x1");
  cout << "x1 = " << x1 << endl;
  cout << "der(x1) = " << der(x1) << endl;

//  DifferentialState x1 = shared_cast<DifferentialState>(ocp.vars("x1"));
  
  // Get all variables
  vector<Variable> v = ocp.vars;
  cout << "v = " << v << endl;
  
  // Make explicit
  ocp.makeExplicit();
  
  // Print again
  ocp.print();
  
  
  cout << "**********************"     << endl;
  cout << "    Pass to ACADO: " << endl;
  cout << "**********************"     << endl;
  
    // DAE input
  vector<vector<SX> > ffcn_in(ACADO_FCN_NUM_IN);
  
  // Time
 ffcn_in[ACADO_FCN_T] = ocp.t;

  // Differential state
  ffcn_in[ACADO_FCN_XD] = ocp.xd;

  // Algebraic state
  ffcn_in[ACADO_FCN_XA] = ocp.xa;

  // Control
  ffcn_in[ACADO_FCN_U] = ocp.u;
  
  // Parameter
  ffcn_in[ACADO_FCN_P] = ocp.p;
  
  // DAE outputs
  vector<vector<SX> > ffcn_out(1);
  ffcn_out[0] = ocp.diffeq;
  ffcn_out[0].insert(ffcn_out[0].end(),ocp.algeq.begin(), ocp.algeq.end());
  
  // Create DAE function
  SXFunction ffcn(ffcn_in,ffcn_out);
  ffcn.setOption("ad_order",1);
  
  // Objective function
  SXFunction mfcn(ffcn_in,vector<vector<SX> >(1,ocp.mterm));
  mfcn.setOption("ad_order",1);
  
  // Path constraint function
  SXFunction cfcn(ffcn_in,vector<vector<SX> >(1,ocp.cfcn));
  cfcn.setOption("ad_order",1);
  
  // Initial constraint function
  SXFunction rfcn(ffcn_in,vector<vector<SX> >(1,ocp.initeq));
  rfcn.setOption("ad_order",1);

  // Create ACADO solver
  AcadoInterface ocp_solver(ffcn,mfcn,cfcn,rfcn);

  // Set options
  ocp_solver.setOption("start_time",ocp.t0);
  ocp_solver.setOption("final_time",ocp.tf);
  ocp_solver.setOption("number_of_shooting_nodes",30);
/*  ocp_solver.setOption("max_num_iterations",2);
  ocp_solver.setOption("kkt_tolerance",1e-2);*/
  
  // Initialize
  ocp_solver.init();

  // Set bounds on states
  ocp_solver.setInput(ocp.cfcn_lb,ACADO_LBC);
  ocp_solver.setInput(ocp.cfcn_ub,ACADO_UBC);
  
  // Give an initial guess
/*  double x_init[] = {0, 0.0410264, 0.0660755, 0.00393984, 0.00556818};
  vector<double>& x_init_all = ocp_solver.input(ACADO_X_GUESS).data();
  for(int i=0; i<5; ++i){
    for(int j=0; j<31; ++j){
      x_init_all[i+5*j] = x_init[i];
    }
  }*/
  ocp_solver.solve();
  
  double cost;
  ocp_solver.getOutput(cost,ACADO_COST);
  cout << "optimal cost = " << cost << endl;

  const vector<double> &xopt = ocp_solver.output(ACADO_X_OPT).data();
  cout << "xopt = " << xopt << endl;

  const vector<double> &uopt = ocp_solver.output(ACADO_U_OPT).data();
  cout << "uopt = " << uopt << endl;

  const vector<double> &popt = ocp_solver.output(ACADO_P_OPT).data();
  cout << "popt = " << popt << endl;
  
  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
