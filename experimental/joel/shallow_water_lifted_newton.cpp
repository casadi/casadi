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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSefcn.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <casadi/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>
#include <nonlinear_programming/lifted_sqp.hpp>

// Initialize the NLP at the optimal solution
bool inititialize_at_solution = false;

// Initialize the NLP at the optimal solution
bool with_ipopt = false;

// Use Gauss-Newton
bool gauss_newton = true;

// Lifted
bool lifted = true;

using namespace CasADi;
using namespace std;

int main(){

  // Physical parameters
  double g = 9.81; // gravity
  double poolwidth = 0.2;
  double drag0 = 2.0; // => u(0)
  double depth0 = 0.01; // => u(1)
  double sprad = 0.03;
  double spheight = 0.01;
  double endtime = 1.0;

  // Discretization
  int numboxes = 20;
  int num_eulersteps = 25;
  int num_measurements = 20;

  // Plotting
  bool plot_progress = false;
  int numboxes_per_plot = 1;

  // Discretization
  int ntimesteps = num_eulersteps*num_measurements;
  double dt = endtime/ntimesteps;
  double dx = poolwidth/numboxes;
  double dy = poolwidth/numboxes;
  vector<double> x(numboxes), y(numboxes);
  for(int i=0; i<numboxes; ++i){
    x[i] = (i+0.5)*dx;
    y[i] = (i+0.5)*dy;
  }

  // Initial conditions
  DMatrix u0 = DMatrix::zeros(numboxes+1,numboxes  );
  DMatrix v0 = DMatrix::zeros(numboxes  ,numboxes+1);
  DMatrix h0 = DMatrix::zeros(numboxes  ,numboxes  ); 
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes; ++j){
      double spdist = sqrt(pow((x[i]-0.04),2) + pow((y[j]-0.04),2));
      if(spdist<sprad/3.0){
	h0.elem(i,j) = spheight * cos(3.0*M_PI*spdist/(2.0*sprad));
      }
    }
  }
  //   h0.printDense();
  
  // Free parameters
  SX drag("b");
  SX depth("H");
  vector<SX> p(2); p[0]=drag; p[1]=depth;
  vector<double> p0(2); p0[0]=drag0; p0[1]=depth0;
  
  // The state at a measurement
  SXMatrix uk = ssym("uk",numboxes+1, numboxes);
  SXMatrix vk = ssym("vk",numboxes  , numboxes+1);
  SXMatrix hk = ssym("hk",numboxes  , numboxes);

  // Mayer term
  SXMatrix hmeas = ssym("h_meas",numboxes  , numboxes);
  vector<SXMatrix> m_in(2);
  m_in[0] = hk;
  m_in[1] = hmeas;
  SX mterm = 0;
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes; ++j){
      SX dev = hk.elem(i,j) - hmeas.elem(i,j);
      mterm += dev*dev;
    }
  }
  SXFunction m(m_in,mterm);
  m.init();
  
  // Take one step of the integrator
  SXMatrix u = uk;
  SXMatrix v = vk;
  SXMatrix h = hk;
  
  // Temporaries
  SX d1 = -dt*g/dx;
  SX d2 = dt*drag;
  
  // Update u
  for(int i=0; i<numboxes-1; ++i){
    for(int j=0; j<numboxes; ++j){
      u.elem(1+i,j) += d1*(h.elem(1+i,j)-h.elem(i,j))- d2*u.elem(1+i,j);
    }
  }
  
  // Update v
  d1 = -dt*g/dy;
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes-1; ++j){
      v.elem(i,j+1) += d1*(h.elem(i,j+1)-h.elem(i,j))- d2*v.elem(i,j+1);
    }
  }
    
  // Update h
  d1 = (-depth*dt)*(1.0/dx);
  d2 = (-depth*dt)*(1.0/dy);
  for(int i=0; i<numboxes; ++i){
    for(int j=0; j<numboxes; ++j){
      h.elem(i,j) += d1*(u.elem(1+i,j)-u.elem(i,j)) + d2*(v.elem(i,j+1)-v.elem(i,j));
    }
  }
  
  // Create an integrator function
  vector<SXMatrix> f_step_in(4);
  f_step_in[0] = p;
  f_step_in[1] = uk;
  f_step_in[2] = vk;
  f_step_in[3] = hk;
  vector<SXMatrix> f_step_out(3);
  f_step_out[0] = u;
  f_step_out[1] = v;
  f_step_out[2] = h;
  SXFunction f_step(f_step_in,f_step_out);
  f_step.setOption("live_variables",true);
  f_step.init();
  cout << "generated single step dynamics (" << f_step.getAlgorithmSize() << " nodes)" << endl;

  // Integrate over one interval
  vector<MX> f_in(4);
  MX P = msym("P",2);
  MX Uk = msym("Uk",numboxes+1, numboxes);
  MX Vk = msym("Vk",numboxes  , numboxes+1);
  MX Hk = msym("Hk",numboxes  , numboxes);
  f_in[0] = P;
  f_in[1] = Uk;
  f_in[2] = Vk;
  f_in[3] = Hk;
  vector<MX> f_inter = f_in;
  vector<MX> f_out;
  for(int j=0; j<num_eulersteps; ++j){
    // Create a call node
    f_out = f_step.call(f_inter);
    
    // Save intermediate state
    f_inter[1] = f_out[0];
    f_inter[2] = f_out[1];
    f_inter[3] = f_out[2];
  }

  // Create an integrator function
  FX f = MXFunction(f_in,f_out);
  f.init();
  cout << "generated discrete dynamics, MX (" << shared_cast<MXFunction>(f).countNodes() << " nodes)" << endl;

  // Expand the discrete dynamics
  if(false){
    f = SXFunction(shared_cast<MXFunction>(f));
    f.init();
    cout << "generated discrete dynamics, SX (" << shared_cast<SXFunction>(f).getAlgorithmSize() << " nodes)" << endl;
  }

  // Measurements
  vector<DMatrix> H_meas;
  H_meas.reserve(num_measurements);

  // Simulate once to generate "measurements"
  f.setInput(p0,0);
  f.setInput(u0,1);
  f.setInput(v0,2);
  f.setInput(h0,3);
  clock_t time1 = clock();
  for(int k=0; k<num_measurements; ++k){
    f.evaluate();
    const DMatrix& u = f.output(0);
    const DMatrix& v = f.output(1);
    const DMatrix& h = f.output(2);
    f.setInput(u,1);
    f.setInput(v,2);
    f.setInput(h,3);
    
    // Save a copy of h
    H_meas.push_back(h);
  }
  clock_t time2 = clock();
  double t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "measurements generated in " << t_elapsed << " seconds." << endl;
  
  // Lifted quantitites
  vector<MX> H_lifted = msym("h_lifted",numboxes,numboxes,num_measurements);

  // Generate full-space NLP
  vector<MX> H_eq;
  MX U = u0;
  MX V = v0;
  MX H = h0;
  int offset = 0;
  for(int k=0; k<num_measurements; ++k){
    // Take a step
    MX f_arg[4] = {P,U,V,H};
    vector<MX> f_res = f.call(vector<MX>(f_arg,f_arg+4));
    U = f_res[0];
    V = f_res[1];
    H = f_res[2];
    
    // Lift the variable
    MX H_def = H;
    H = H_lifted[k];
    H_eq.push_back(H_def-H);
  }
  
  // Lifted NLP function
  vector<MX> fff_in;
  fff_in.push_back(P);
  fff_in.insert(fff_in.end(),H_lifted.begin(),H_lifted.end());
  vector<MX> fff_out = H_eq;
  MXFunction fff(fff_in,fff_out);
  fff.init();
  cout << "Generated lifted problem" << endl;
  SXFunction fff_sx(fff);
  fff_sx.init();
  cout << "generated lifted problem, SX (" << fff_sx.getAlgorithmSize() << " nodes)" << endl;
  
  // NLP variables
  SXMatrix nlp_x;
  for(int k=0; k<fff_sx.getNumInputs(); ++k){
    nlp_x.append(flatten(fff_sx.inputSX(k)));
  }
  
  // NLP objective function terms
  SXMatrix nlp_f;
  for(int k=0; k<num_measurements; ++k){
    nlp_f.append(flatten(fff_sx.inputSX(1+k)-SXMatrix(H_meas[k])));
  }
  
  // NLP constraint function terms
  SXMatrix nlp_g;
  for(int k=0; k<num_measurements; ++k){
    nlp_g.append(flatten(fff_sx.outputSX(k)));
  }
  
  // Objective term
  if(!gauss_newton){
    nlp_f = inner_prod(nlp_f,nlp_f);
  }
  
  // Formulate the NLP
  SXFunction ffcn(nlp_x,nlp_f);
  SXFunction gfcn(nlp_x,nlp_g);
  
  // Gauss-Newton Hessian approximation for Ipopt
//   SXMatrix JF = jacobian(nlp_f,nlp_x);
//   SXFunction hfcn(nlp_x,mul(trans(JF),JF));
  
  // Solve with IPOPT
  NLPSolver nlp_solver;
  if(with_ipopt){
    casadi_assert(!gauss_newton);
    nlp_solver = IpoptSolver(ffcn,gfcn);
    nlp_solver.setOption("generate_hessian",true);
  } else {
    nlp_solver = LiftedSQP(ffcn,gfcn);
    if(gauss_newton) nlp_solver.setOption("gauss_newton",true);
    nlp_solver.setOption("qp_solver",Interfaces::QPOasesSolver::creator);
    Dictionary qp_solver_options;
    qp_solver_options["printLevel"] = "none";
    //qp_solver_options["verbose"] = true;
    nlp_solver.setOption("qp_solver_options",qp_solver_options);
    if(lifted) nlp_solver.setOption("num_lifted",nlp_g.size());
    nlp_solver.setOption("toldx",1e-9);
    nlp_solver.setOption("verbose",true);
  }
  nlp_solver.init();
  
  if(inititialize_at_solution){
    nlp_solver.input(NLP_X_INIT)[0] = drag0;
    nlp_solver.input(NLP_X_INIT)[1] = depth0;
    vector<double>::iterator it=nlp_solver.input(NLP_X_INIT).begin()+2;
    for(int k=0; k<num_measurements; ++k){
      copy(H_meas[k].begin(),H_meas[k].end(),it);
      it += numboxes*numboxes;
    }
  } else {
    nlp_solver.input(NLP_X_INIT)[0] = 0.5;
    nlp_solver.input(NLP_X_INIT)[1] = 0.01;
  }
  nlp_solver.input(NLP_LBG).setZero();
  nlp_solver.input(NLP_UBG).setZero();
//   nlp_solver.input(NLP_LBX).setAll(-1000.);
//   nlp_solver.input(NLP_UBX).setAll( 1000.);
//   nlp_solver.input(NLP_LBX)[0] = 0;
//   nlp_solver.input(NLP_LBX)[1] = 0;
  nlp_solver.solve();

//  cout << "x_opt = " << nlp_solver.output(NLP_X_OPT) << endl;
  cout << "f_opt = " << nlp_solver.output(NLP_COST) << endl;
  return 0;
}
