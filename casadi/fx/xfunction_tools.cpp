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
    
SXFunction vec (const SXFunction &a) {
  // Pass null if input is null
  if (a.isNull()) return SXFunction();
  
  /// Get the SX input and output vectors
  std::vector<SXMatrix> symbolicInputSX = a.inputsSX();
  std::vector<SXMatrix> symbolicOutputSX = a.outputsSX();
  
  // Apply vec to them
  for (int i=0;i<symbolicInputSX.size();++i)
    symbolicInputSX[i] = vec(symbolicInputSX[i]);
  for (int i=0;i<symbolicOutputSX.size();++i)
    symbolicOutputSX[i] = vec(symbolicOutputSX[i]);
    
  // Make a new function with the vecced input/outputs
  SXFunction ret(symbolicInputSX,symbolicOutputSX);
  
  // Initialize it if a was
  if (a.isInit()) ret.init();
  return ret;
}


MXFunction vec (const FX &a_) {
  FX a = a_;

  // Pass null if input is null
  if (a.isNull()) return MXFunction();
  
  // Get the MX inputs, only used for shape
  const std::vector<MX> &symbolicInputMX = a.symbolicInput();
  // Have a vector with MX that have the shape of vec(symbolicInputMX )
  std::vector<MX> symbolicInputMX_vec(a.getNumInputs());
  // Make vector valued MX's out of them
  std::vector<MX> symbolicInputMX_vec_reshape(a.getNumInputs());

  // Apply the vec-transformation to the inputs
  for (int i=0;i<symbolicInputMX.size();++i) {
    std::stringstream s;
    s << "X_flat_" << i;
    symbolicInputMX_vec[i] = MX(s.str(),vec(symbolicInputMX[i].sparsity()));
    symbolicInputMX_vec_reshape[i] = reshape(symbolicInputMX_vec[i],symbolicInputMX[i].sparsity());
  }
  
  // Call the original function with the vecced inputs
  std::vector<MX> symbolicOutputMX = a.call(symbolicInputMX_vec_reshape);
  
  // Apply the vec-transformation to the outputs
  for (int i=0;i<symbolicOutputMX.size();++i)
    symbolicOutputMX[i] = vec(symbolicOutputMX[i]);
    
  // Make a new function with the vecced input/outputs
  MXFunction ret(symbolicInputMX_vec,symbolicOutputMX);
  
  // Initialize it if a was
  if (a.isInit()) ret.init();
  return ret;

}  
    
    
} // namespace CasADi

