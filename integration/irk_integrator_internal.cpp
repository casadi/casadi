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

#include "irk_integrator_internal.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

  IRKIntegratorInternal::IRKIntegratorInternal(const FX& f, const FX& g) : RKBaseInternal(f,g){
    addOption("interpolation_order",           OT_INTEGER,  3,  "Order of the interpolating polynomials");
    addOption("collocation_scheme",            OT_STRING,  "radau",  "Collocation scheme","radau|legendre");
    setOption("name","unnamed_irk_integrator");
  }

  void IRKIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    RKBaseInternal::deepCopyMembers(already_copied);
  }

  IRKIntegratorInternal::~IRKIntegratorInternal(){
  }

  void IRKIntegratorInternal::init(){
  
    // Call the base class init
    RKBaseInternal::init();
  
  }

  void IRKIntegratorInternal::setupFG(){

    // Interpolation order
    int deg = getOption("interpolation_order");

    // All collocation time points
    std::vector<long double> tau_root = collocationPointsL(deg,getOption("collocation_scheme"));

    // Coefficients of the collocation equation
    vector<vector<double> > C(deg+1,vector<double>(deg+1,0));
      
    // Coefficients of the continuity equation
    vector<double> D(deg+1,0);
      
    // Coefficients of the quadratures
    vector<double> B(deg+1,0);

    // For all collocation points
    for(int j=0; j<deg+1; ++j){

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for(int r=0; r<deg+1; ++r){
        if(r!=j){
          p *= Polynomial(-tau_root[r],1)/(tau_root[j]-tau_root[r]);
        }
      }
    
      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D[j] = zeroIfSmall(p(1.0L));
    
      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for(int r=0; r<deg+1; ++r){
        C[j][r] = zeroIfSmall(dp(tau_root[r]));
      }
        
      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      B[j] = zeroIfSmall(ip(1.0L));
    }

    // Symbolic inputs
    MX x0 = msym("x0",f_.input(DAE_X).sparsity());
    MX p = msym("p",f_.input(DAE_P).sparsity());
    MX t = msym("t",f_.input(DAE_T).sparsity());

    // Implicitly defined variables (z and x)
    MX v = msym("v",deg*(nx_+nz_));
    vector<int> v_offset(1,0);
    for(int d=0; d<deg; ++d){
      v_offset.push_back(v_offset.back()+nx_);
      v_offset.push_back(v_offset.back()+nz_);
    }
    vector<MX> vv = vertsplit(v,v_offset);
    vector<MX>::const_iterator vv_it = vv.begin();

    // Collocated states
    vector<MX> x(deg+1), z(deg+1);
    x[0] = x0;
    for(int d=1; d<=deg; ++d){
      x[d] = *vv_it++;
      z[d] = *vv_it++;
    }

    // Collocation time points
    vector<MX> tt(deg+1);
    for(int d=0; d<=deg; ++d){
      tt[d] = t + h_*tau_root[d];
    }

    // Equations that implicitly define v
    vector<MX> eq;

    // Quadratures
    MX qf = MX::zeros(f_.output(DAE_QUAD).sparsity());

    // End state
    MX xf = D[0]*x0;

    // For all collocation points
    for(int j=1; j<deg+1; ++j){

      // Get an expression for the state derivative at the collocation point
      MX xp_j = 0;
      for(int r=0; r<deg+1; ++r){
        xp_j += (C[r][j]/h_) * x[r];
      }

      // Evaluate the DAE
      vector<MX> f_in(DAE_NUM_IN);
      f_in[DAE_T] = tt[j];
      f_in[DAE_P] = p;
      f_in[DAE_X] = x[j];
      f_in[DAE_Z] = z[j];
      vector<MX> f_out = f_.call(f_in);
      
      // Add collocation equation
      eq.push_back(f_out[DAE_ODE] - xp_j);
        
      // Add the algebraic conditions
      eq.push_back(f_out[DAE_ALG]);

      // Add contribution to the final state
      xf += D[j]*x[j];
        
      // Add contribution to quadratures
      qf += (B[j]*h_)*f_out[DAE_QUAD];
    }

    // Form forward discrete time dynamics
    vector<MX> F_in(DAE_NUM_IN);
    F_in[DAE_T] = t;
    F_in[DAE_X] = x0;
    F_in[DAE_P] = p;
    F_in[DAE_Z] = v;
    vector<MX> F_out(DAE_NUM_OUT);
    F_out[DAE_ODE] = xf;
    F_out[DAE_ALG] = vertcat(eq);
    F_out[DAE_QUAD] = qf;
    F_ = MXFunction(F_in,F_out);
    F_.init();
  }
  

  double IRKIntegratorInternal::zeroIfSmall(double x){
    return fabs(x) < numeric_limits<double>::epsilon() ? 0 : x;
  }

} // namespace CasADi
