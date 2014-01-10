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

#include "rk_integrator_internal.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

  RKIntegratorInternal::RKIntegratorInternal(const FX& f, const FX& g) : RKBaseInternal(f,g){
  }

  void RKIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    RKBaseInternal::deepCopyMembers(already_copied);
  }

  RKIntegratorInternal::~RKIntegratorInternal(){
  }

  void RKIntegratorInternal::init(){
    // Call the base class init
    RKBaseInternal::init();

    // Algebraic variables not (yet?) supported
    casadi_assert_message(nz_==0 && nrz_==0, "Explicit Runge-Kutta integrators do not support algebraic variables");

    // Symbolic inputs
    MX x0 = msym("x0",f_.input(DAE_X).sparsity());
    MX z0 = msym("z0",f_.input(DAE_Z).sparsity());
    MX p = msym("p",f_.input(DAE_P).sparsity());
    MX t = msym("t",f_.input(DAE_T).sparsity());
    
    // Arguments when calling f
    vector<MX> f_arg(DAE_NUM_IN);
    vector<MX> f_res;
    f_arg[DAE_P] = p;
    f_arg[DAE_Z] = z0;

    // k1
    f_arg[DAE_T] = t;
    f_arg[DAE_X] = x0;
    f_res = f_.call(f_arg);
    MX k1 = f_res[DAE_ODE];
    MX k1_quad = f_res[DAE_QUAD];

    // k2
    f_arg[DAE_T] = t + h_/2;
    f_arg[DAE_X] = x0 + (h_/2) * k1;
    f_res = f_.call(f_arg);
    MX k2 = f_res[DAE_ODE];
    MX k2_quad = f_res[DAE_QUAD];
    
    // k3
    f_arg[DAE_X] = x0 + (h_/2) * k2;
    f_res = f_.call(f_arg);
    MX k3 = f_res[DAE_ODE];
    MX k3_quad = f_res[DAE_QUAD];
    
    // k4
    f_arg[DAE_T] = t + h_;
    f_arg[DAE_X] = x0 + h_ * k3;
    f_res = f_.call(f_arg);
    MX k4 = f_res[DAE_ODE];
    MX k4_quad = f_res[DAE_QUAD];
    
    // Take step
    MX xf = x0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
    MX qf = (h_/6)*(k1_quad + 2*k2_quad + 2*k3_quad + k4_quad);

    // Define discrete time dynamics
    f_arg[DAE_T] = t;
    f_arg[DAE_X] = x0;
    f_arg[DAE_P] = p;
    f_arg[DAE_Z] = z0;;
    f_res[DAE_ODE] = xf;
    f_res[DAE_QUAD] = qf;
    f_res[DAE_ALG] = MX();
    F_ = MXFunction(f_arg,f_res);
    F_.init();

    // Algebraic variables not (yet?) supported
    casadi_assert_message(g_.isNull(), "Not implemented");
  }
  
} // namespace CasADi
