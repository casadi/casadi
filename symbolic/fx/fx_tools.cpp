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

#include "fx_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "integrator.hpp"
#include <iostream>

using namespace std;

namespace CasADi{
    
FX parameterizeTime(FX dae) {

   // dimensionless time
   MX tau = MX::sym("tau");
   
   // The dae parameters augmented with t0 and tf
   MX P = MX::sym("P",1,2+dae.input(DAE_P).size());
   MX t0 = P[0];
   MX tf = P[1];
   
   std::vector<MX> dae_in(DAE_NUM_IN);
   std::vector<MX> dae_input = dae.symbolicInput();
   
   if (dae.input(DAE_T).size()==1) {
     dae_in[DAE_T]    = t0 + (tf-t0)*tau;
   }
   
   dae_in[DAE_P]    = reshape(trans(P[range(2,2+dae.input(DAE_P).size())]),dae.input(DAE_P).sparsity());
   dae_in[DAE_X]    = dae_input[DAE_X];

   std::vector<MX> ret_in(DAE_NUM_IN);
   ret_in[DAE_T]    = tau;
   ret_in[DAE_P]    = P;
   ret_in[DAE_X]    = dae_input[DAE_X];

   std::vector<MX> ret_out(DAE_NUM_OUT);
   ret_out[DAE_ODE] = (tf-t0)*dae.call(dae_in)[0];
   
   MXFunction ret(ret_in,ret_out);
   if (dae.isInit()) ret.init();
   
   // Expand if dae was an SXFunction
   if(is_a<SXFunction>(dae)){
     if(!ret.isInit()) ret.init();
     SXFunction ret_sx(ret);
     if(ret.isInit()) ret_sx.init();
     return ret_sx;
   } else {
     return ret;
   }
 }
 
FX parameterizeTimeOutput(FX f) {
  // dimensionless time
   MX tau = MX::sym("tau");
   
   // The f parameters augmented with t0 and tf
   MX P = MX::sym("P",1,2+f.input(DAE_P).size());
   MX t0 = P[0];
   MX tf = P[1];
   
   std::vector<MX> f_in(DAE_NUM_IN);
   std::vector<MX> f_input = f.symbolicInput();
   
   if (f.input(DAE_T).size()==1) {
     f_in[DAE_T]    = t0 + (tf-t0)*tau;
   }
   
   f_in[DAE_P]    = reshape(trans(P[range(2,2+f.input(DAE_P).size())]),f.input(DAE_P).sparsity());
   f_in[DAE_X]    = f_input[DAE_X];

   std::vector<MX> ret_in(DAE_NUM_IN);
   ret_in[DAE_T]    = tau;
   ret_in[DAE_P]    = P;
   ret_in[DAE_X]    = f_input[DAE_X];

   MXFunction ret(ret_in,f.call(f_in));
   
   if (f.isInit()) ret.init();
   
   // Expand if f was an SXFunction
   if(is_a<SXFunction>(f)){
     if(!ret.isInit()) ret.init();
     SXFunction ret_sx(ret);
     if(ret.isInit()) ret_sx.init();
     return ret_sx;
   } else {
     return ret;
   }
}



} // namespace CasADi

