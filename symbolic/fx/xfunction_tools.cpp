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

#include "xfunction_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/sparsity_tools.hpp"
#include "../stl_vector_tools.hpp"

namespace CasADi{
    
  SXFunction vec(SXFunction f) {
    // Pass null if input is null
    if (f.isNull()) return SXFunction();
  
    /// Get the SXElement input and output vectors
    std::vector<SX> f_in = f.inputExpr();
    std::vector<SX> f_out = f.outputExpr();
  
    // Apply vec to them
    for(std::vector<SX>::iterator it=f_in.begin(); it!=f_in.end(); ++it) *it = vec(*it);
    for(std::vector<SX>::iterator it=f_out.begin(); it!=f_out.end(); ++it) *it = vec(*it);
    
    // Make a new function with the vectorized input/outputs
    SXFunction ret(f_in,f_out);
  
    // Initialize it if it was initialized
    if (f.isInit()) ret.init();
    return ret;
  }

  MXFunction vec(FX f) {

    // Pass null if input is null
    if (f.isNull()) return MXFunction();
  
    // Get the MX inputs, only used for shape
    const std::vector<MX> &f_in = f.symbolicInput();

    // Have a vector with MX that have the shape of vec(symbolicInput())
    std::vector<MX> f_in_vec(f_in.size());

    // Make vector valued MX's out of them
    std::vector<MX> f_in_vec_reshaped(f_in.size());

    // Apply the vec-transformation to the inputs
    for(int i=0; i<f_in.size(); ++i){
      std::stringstream s;
      s << "X_flat_" << i;
      f_in_vec[i] = MX::sym(s.str(),vec(f_in[i].sparsity()));
      f_in_vec_reshaped[i] = reshape(f_in_vec[i],f_in[i].sparsity());
    }
  
    // Call the original function with the vectorized inputs
    std::vector<MX> f_out = f.call(f_in_vec_reshaped);
  
    // Apply the vec-transformation to the outputs
    for(std::vector<MX>::iterator it=f_out.begin(); it!=f_out.end(); ++it) *it = vec(*it);
    
    // Make a new function with the vectorized input/outputs
    MXFunction ret(f_in_vec,f_out);
  
    // Initialize it if it was initialized
    if (f.isInit()) ret.init();
    return ret;
  }  
    
} // namespace CasADi

