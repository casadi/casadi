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
  int numboxes = 30;
  int num_eulersteps = 100;
  int num_measurements = 100;

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

  // Measurement
  vector<double> h_meas;
  h_meas.reserve(numboxes*numboxes*num_measurements);

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
    h_meas.insert(h_meas.end(),h.begin(),h.end());
  }
  clock_t time2 = clock();
  double t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "measurements generated in " << t_elapsed << " seconds." << endl;
  
  // Difference vector d
  MX d = msym("d",numboxes*numboxes*num_measurements);

  // Generate reconstruction function
  vector<MX> H_sim;
  MX U = u0;
  MX V = v0;
  MX H = h0;
  int offset = 0;
  for(int k=0; k<num_measurements; ++k){
    // Get the local difference vector
    int offset_next = offset+numboxes*numboxes;
    MX dk = reshape(d[range(offset,offset_next)],H.sparsity());
    offset = offset_next;
    
    // Take a step
    MX f_arg[4] = {P,U,V,H};
    vector<MX> f_res = f.call(vector<MX>(f_arg,f_arg+4));
    U = f_res[0];
    V = f_res[1];
    H = f_res[2];
    H -= dk;
    H_sim.push_back(H);
  }
  MX zfcn_in[2] = {P,d};
  MX zfcn_out[1] = {vertcat(H_sim)};
  MXFunction zfcn(vector<MX>(zfcn_in,zfcn_in+2),vector<MX>(zfcn_out,zfcn_out+1));
  zfcn.init();
  cout << "generated zfcn, MX (" << shared_cast<MXFunction>(zfcn).countNodes() << " nodes)" << endl;

  time1 = clock();
  zfcn.evaluate(0,1);
  time2 = clock();
  t_elapsed = double(time2-time1)/CLOCKS_PER_SEC;
  cout << "zfcn evaluated in " << t_elapsed << " seconds." << endl;
  
  
//   vector<MX> f_in(4);
// 
//   
//   MX 
//   f.setInput(p0,0);
//   f.setInput(u0,1);
//   f.setInput(v0,2);
//   f.setInput(h0,3);
//   clock_t time1 = clock();
//     f.evaluate();
//     const DMatrix& u = f.output(0);
//     const DMatrix& v = f.output(1);
//     const DMatrix& h = f.output(2);
//     f.setInput(u,1);
//     f.setInput(v,2);
//     f.setInput(h,3);
//     
//     // Save a copy of h
//     h_meas.insert(h_meas.end(),h.begin(),h.end());
//   }
//   
//   
//   
//   
  
  
  
  cout << "numboxes*numboxes*num_measurements = " << numboxes*numboxes*num_measurements << endl;
  cout << "numboxes*numboxes*num_measurements = " << numboxes*numboxes*num_measurements << endl;
  
  
  // Now construct the

  
  return 0;

}
