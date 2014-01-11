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
    MX t_k1 = f_arg[DAE_T] = t;
    MX x_k1 = f_arg[DAE_X] = x0;
    f_res = f_.call(f_arg);
    MX k1 = f_res[DAE_ODE];
    MX k1q = f_res[DAE_QUAD];

    // k2
    MX t_k2 = f_arg[DAE_T] = t + h_/2;
    MX x_k2 = f_arg[DAE_X] = x0 + (h_/2) * k1;
    f_res = f_.call(f_arg);
    MX k2 = f_res[DAE_ODE];
    MX k2q = f_res[DAE_QUAD];
    
    // k3
    MX t_k3 = t_k2;
    MX x_k3 = f_arg[DAE_X] = x0 + (h_/2) * k2;
    f_res = f_.call(f_arg);
    MX k3 = f_res[DAE_ODE];
    MX k3q = f_res[DAE_QUAD];
    
    // k4
    MX t_k4 = f_arg[DAE_T] = t + h_;
    MX x_k4 = f_arg[DAE_X] = x0 + h_ * k3;
    f_res = f_.call(f_arg);
    MX k4 = f_res[DAE_ODE];
    MX k4q = f_res[DAE_QUAD];
    
    // Take step
    MX xf = x0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
    MX qf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

    // Define discrete time dynamics
    f_arg[DAE_T] = t;
    f_arg[DAE_X] = x0;
    f_arg[DAE_P] = p;
    f_arg[DAE_Z] = z0;
    f_res[DAE_ODE] = xf;
    f_res[DAE_QUAD] = qf;
    f_res[DAE_ALG] = MX();
    F_ = MXFunction(f_arg,f_res);
    F_.init();

    // Backward integration
    if(!g_.isNull()){
      // Symbolic inputs
      MX rx0 = msym("x0",g_.input(RDAE_RX).sparsity());
      MX rz0 = msym("z0",g_.input(RDAE_RZ).sparsity());
      MX rp = msym("p",g_.input(RDAE_RP).sparsity());

      // Arguments when calling g
      vector<MX> g_arg(RDAE_NUM_IN);
      vector<MX> g_res;
      g_arg[RDAE_P] = p;
      g_arg[RDAE_Z] = z0;
      g_arg[RDAE_RZ] = rz0;
      g_arg[RDAE_RP] = rp;

      // k1
      g_arg[RDAE_T] = t_k4;
      g_arg[RDAE_X] = x_k4;
      g_arg[RDAE_RX] = rx0;
      g_res = g_.call(g_arg);
      k1 = g_res[RDAE_ODE];
      k1q = g_res[RDAE_QUAD];

      // k2
      g_arg[RDAE_T] = t_k3;
      g_arg[RDAE_X] = x_k3;
      g_arg[RDAE_RX] = rx0 + (h_/2) * k1;
      g_res = g_.call(g_arg);
      k2 = g_res[RDAE_ODE];
      k2q = g_res[RDAE_QUAD];

      // k3
      g_arg[RDAE_T] = t_k2;
      g_arg[RDAE_X] = x_k2;
      g_arg[RDAE_RX] = rx0 + (h_/2) * k2;
      g_res = g_.call(g_arg);
      k3 = g_res[RDAE_ODE];
      k3q = g_res[RDAE_QUAD];

      // k4
      g_arg[RDAE_T] = t_k1;
      g_arg[RDAE_X] = x_k1;
      g_arg[RDAE_RX] = rx0 + h_ * k3;
      g_res = g_.call(g_arg);
      k4 = g_res[RDAE_ODE];
      k4q = g_res[RDAE_QUAD];

      // Take step
      MX rxf = rx0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
      MX rqf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      g_arg[RDAE_T] = t;
      g_arg[RDAE_X] = x0;
      g_arg[RDAE_P] = p;
      g_arg[RDAE_Z] = z0;
      g_arg[RDAE_RX] = rx0;
      g_arg[RDAE_RP] = rp;
      g_arg[RDAE_RZ] = rz0;
      g_res[RDAE_ODE] = rxf;
      g_res[RDAE_QUAD] = rqf;
      g_res[RDAE_ALG] = MX();
      G_ = MXFunction(g_arg,g_res);
      G_.init();
    }
  }
  
} // namespace CasADi
