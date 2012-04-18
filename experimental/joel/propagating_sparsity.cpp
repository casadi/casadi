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

#include <casadi/stl_vector_tools.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/fx/fx_internal.hpp>
#include <casadi/fx/sx_function.hpp>

using namespace CasADi;
using namespace std;

int main(){
  SXMatrix x = ssym("x",3);
  SXMatrix z = x[0]*x[0]+x[2] + 3;
  FX f = SXFunction(x,z);
  f.init();
  
  bvec_t* f_in = get_bvec_t(f.input().data());
  f_in[0] = bvec_t(1) << 0;
  f_in[1] = bvec_t(1) << 2;
  f_in[2] = (bvec_t(1) << 4) | (bvec_t(1) << 63);
  
  bvec_t* f_out = get_bvec_t(f.output().data());
  f_out[0] = 0;
  
  f->spInit(true);
  f->spInit(false);
  f->spEvaluate(true);
  
  for(int k=0; k<bvec_size; ++k){
    if(k%4==0) cout << " ";
    if(f_out[0] & (bvec_t(1)<<k)){
      cout << 1;
    } else {
      cout << 0;
    }
  }
  cout << endl;
  return 0;
}
