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
#include "symbolic/mx/mx_tools.hpp"
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
  
  FX explicitRKIntegrator(FX& f, const MX& tf, int order, int ne) {
    casadi_assert_message(ne>=1,"Parameter ne (number of elements must be at least 1), but got " << ne << ".");
    casadi_assert_message(order==4,"Only RK order 4 is supported now.");
    casadi_assert_message(f.getNumInputs()==DAE_NUM_IN && f.getNumOutputs()==DAE_NUM_OUT,"Supplied function must adhere to dae scheme.");
    casadi_assert_message(f.input(DAE_ALG).empty(),"Supplied function cannot have algebraic states.");

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
   //return FX();
  }

} // namespace CasADi

