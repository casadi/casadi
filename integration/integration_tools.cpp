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

#include "integration_tools.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/fx/implicit_function.hpp"
#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/integrator.hpp"

#include <vector>

using namespace std;

namespace CasADi{
  
  std::vector<double> collocationPoints(int order, const std::string& scheme) {
    if (scheme=="radau") {
      casadi_assert_message(order>0 && order<10,"Error in collocationPoints(order, scheme): only order up to 9 supported for scheme 'radau', but got " << order << ".");
      return std::vector<double>(radau_points[order],radau_points[order]+order+1);
    } else if (scheme=="legendre") {
      casadi_assert_message(order>0 && order<10,"Error in collocationPoints(order, scheme): only order up to 9 supported for scheme 'legendre', but got " << order << ".");
      return std::vector<double>(legendre_points[order],legendre_points[order]+order+1);
    } else {
      casadi_error("Error in collocationPoints(order, scheme): unknown scheme '" << scheme << "'. Select one of 'radau', 'legendre'.");
    }
  }
  
  FX explicitRK(FX& f, const MX& tf, int order, int ne) {
    casadi_assert_message(ne>=1,"Parameter ne (number of elements must be at least 1), but got " << ne << ".");
    casadi_assert_message(order==4,"Only RK order 4 is supported now.");
    casadi_assert_message(f.getNumInputs()==DAE_NUM_IN && f.getNumOutputs()==DAE_NUM_OUT,"Supplied function must adhere to dae scheme.");
    casadi_assert_message(f.output(DAE_ALG).empty() && f.input(DAE_Z).empty(),"Supplied function cannot have algebraic states.");
    casadi_assert_message(f.output(DAE_QUAD).empty(),"Supplied function cannot have quadrature states.");

    MX X = msym("X",f.input(DAE_X).sparsity());
    MX P = msym("P",f.input(DAE_P).sparsity());
    MX X0 = X;
    MX t = 0;
    MX dt = tf/ne;
    
    std::vector<double> b(order);
    b[0]=1.0/6;b[1]=1.0/3;b[2]=1.0/3;b[3]=1.0/6;

    std::vector<double> c(order);
    c[0]=0;c[1]=1.0/2;c[2]=1.0/2;c[3]=1;
    
    std::vector< std::vector<double> > A(order-1);
    A[0].resize(1);A[0][0]=1.0/2;
    A[1].resize(2);A[1][0]=0;A[1][1]=1.0/2;
    A[2].resize(3);A[2][0]=0;A[2][1]=0;A[2][2]=1;
    
    std::vector<MX> k(order);
    
    for (int i=0;i<ne;++i) {
      for (int j=0;j<order;++j) {
        MX XL = 0;
        for (int jj=0;jj<j;++jj) {
          XL+=k.at(jj)*A.at(j-1).at(jj);
        }
        //std::cout << "help: " << A.at(j-1) << "," << c.at(j) << std::endl;
        k[j] = dt*f.call(daeIn("x",X+XL,"p",P,"t",t+dt*c.at(j)))[DAE_ODE];
      }
      
      for (int j=0;j<order;++j) {
        X += b.at(j)*k.at(j);
 
      }
      t+= dt;
   }

   MXFunction ret(integratorIn("x0",X0,"p",P),integratorOut("xf",X));
    
   return ret;
  }
  
  void collocationInterpolators(const std::vector<double> & tau_root, std::vector< std::vector<double> > &C, std::vector< double > &D) {
  
    // Find the degree of the interpolation
    int deg = tau_root.size()-1;
    
    
    // Allocate storage space for resulting coefficients
    C.resize(deg+1);
    for (int i=0;i<deg+1;++i) {
      C[i].resize(deg+1);
    }
    D.resize(deg+1);

    // Collocation point
    SXMatrix tau = ssym("tau");

    // For all collocation points
    for(int j=0; j<deg+1; ++j){
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      SXMatrix L = 1;
      for(int j2=0; j2<deg+1; ++j2){
        if(j2 != j){
          L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2]);
        }
      }
    
      SXFunction lfcn(tau,L);
      lfcn.init();
  
      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      lfcn.setInput(1.0);
      lfcn.evaluate();
      D[j] = lfcn.output().at(0);

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      for(int j2=0; j2<deg+1; ++j2){
        
        lfcn.setInput(tau_root[j2]);
        lfcn.setFwdSeed(1.0);
        lfcn.evaluate(1,0);
        C[j2][j] = lfcn.fwdSens().at(0);
      }
    }

  }
  
  FX implicitRK(FX& f, implicitFunctionCreator impl, const Dictionary& impl_options, const MX& tf, int order, const std::string& scheme, int ne) {
    casadi_assert_message(ne>=1,"Parameter ne (number of elements must be at least 1), but got " << ne << ".");
    casadi_assert_message(order==4,"Only RK order 4 is supported now.");
    casadi_assert_message(f.getNumInputs()==DAE_NUM_IN && f.getNumOutputs()==DAE_NUM_OUT,"Supplied function must adhere to dae scheme.");
    casadi_assert_message(f.output(DAE_QUAD).empty(),"Supplied function cannot have quadrature states.");
    
    // Obtain collocation points
    std::vector<double> tau_root = collocationPoints(order,"legendre");
    
    // Retrieve collocation interpolating matrices
    std::vector < std::vector <double> > C;
    std::vector < double > D;
    collocationInterpolators(tau_root,C,D);
    
    // Retrieve problem dimensions
    int nx = f.input(DAE_X).size();
    int nz = f.input(DAE_Z).size();
    int np = f.input(DAE_P).size();
    
    //Variables for one finite element
    MX X = msym("X",nx);
    MX P = msym("P",np);
    MX V = msym("V",order*(nx+nz)); // Unknowns
    
    MX X0 = X;
    
    // Components of the unknowns that correspond to states at collocation points
    std::vector<MX> Xc;Xc.reserve(order);
    Xc.push_back(X0);
    
    // Components of the unknowns that correspond to algebraic states at collocation points
    std::vector<MX> Zc;Zc.reserve(order);
    
    // Splitting the unknowns
    std::vector<int> splitPositions = range(0,order*nx,nx);
    if (nz>0) {
      std::vector<int> Zc_pos = range(order*nx,order*nx+(order+1)*nz,nz);
      splitPositions.insert( splitPositions.end(), Zc_pos.begin(), Zc_pos.end() );
    } else {
      splitPositions.push_back(order*nx);
    }
    std::vector<MX> Vs = vertsplit(V,splitPositions);
    
    // Extracting unknowns from Z
    for (int i=0;i<order;++i) {
      Xc.push_back(X0+Vs[i]);
    }
    if (nz>0) {
      for (int i=0;i<order;++i) {
        Zc.push_back(Vs[order+i]);
      }
    }
    
    // Get the collocation Equations (that define V)
    std::vector<MX> V_eq;
    
    // Local start time
    MX t0_l=msym("t0");
    MX h = msym("h");

    for (int j=1;j<order+1;++j) {
      // Expression for the state derivative at the collocation point
      MX xp_j = 0;
      for (int r=0;r<order+1;++r) {
        xp_j+= C[j][r]*Xc[r];
      }
      // Append collocation equations & algebraic constraints
      std::vector<MX> f_out;
      MX t_l = t0_l+tau_root[j]*h;
      if (nz>0) {
        f_out = f.call(daeIn("t",t_l,"x",Xc[j],"p",P,"z",Zc[j-1]));
      } else {
        f_out = f.call(daeIn("t",t_l,"x",Xc[j],"p",P));
      }
      V_eq.push_back(h*f_out[DAE_ODE]-xp_j);
      V_eq.push_back(f_out[DAE_ALG]);
      
    }

    // Root-finding function, implicitly defines V as a function of X0 and P
    std::vector<MX> vfcn_inputs;
    vfcn_inputs.push_back(V);
    vfcn_inputs.push_back(X);
    vfcn_inputs.push_back(P);
    vfcn_inputs.push_back(t0_l);
    vfcn_inputs.push_back(h);
    
    FX vfcn = MXFunction(vfcn_inputs,vertcat(V_eq));
    vfcn.init();
    
    try {
      // Attempt to convert to SXFunction to decrease overhead
      vfcn = SXFunction(vfcn);
      vfcn.init();
    } catch (CasadiException & e) {
      //
    }
    
    // Create a implicit function instance to solve the system of equations
    ImplicitFunction ifcn = impl(vfcn,FX(),LinearSolver());
    ifcn.setOption(impl_options);
    ifcn.init();
    
    // Get an expression for the state at the end of the finite element
    std::vector<MX> ifcn_call_in(4);
    std::copy(vfcn_inputs.begin()+1,vfcn_inputs.end(),ifcn_call_in.begin());
    std::vector<MX> ifcn_call_out = ifcn.eval(ifcn_call_in);
    Vs = vertsplit(ifcn_call_out[0],splitPositions);
    
    MX XF = 0;
    for (int i=0;i<order+1;++i) {
      XF += D[i]*(i==0? X : X + Vs[i-1]);
    }
    
    
    // Get the discrete time dynamics
    MXFunction F = MXFunction(ifcn_call_in,XF);
    F.init();
    
    // Loop over all finite elements
    MX h_ = tf/ne;
    MX t0_ = 0;
    
    for (int i=0;i<ne;++i) {
      std::vector<MX> F_in;
      F_in.push_back(X);
      F_in.push_back(P);
      F_in.push_back(t0_);
      F_in.push_back(h_);
      t0_+= h_;
      std::vector<MX> F_out = F.call(F_in);
      X = F_out[0];
    }
    
    // Create a ruturn function with Integrator signature
    MXFunction ret = MXFunction(integratorIn("x0",X0,"p",P),integratorOut("xf",X));
    ret.init();
    
    return ret;
  
  }

} // namespace CasADi

