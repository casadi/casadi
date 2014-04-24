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

%{
#include "casadi/core/mx/mx.hpp"
#include "casadi/core/mx/mx_tools.hpp"
%}

%include "casadi/core/mx/mx.hpp"




%template(SparsityVector) std::vector<casadi::Sparsity>;



%extend casadi::MX{
  
  %matrix_helpers(casadi::MX)
  
  #ifdef SWIGPYTHON
  %python_array_wrappers(1002.0)
  
  %pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    import numpy as np
    if np.__version__=="1.8.1": #1083
      return np.array(np.nan)
    raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")
    
  def __iter__(self):
    for i in range(self.size()):
      yield self[i]
    
  %}
  #endif //SWIGPYTHON
  
  binopsrFull(casadi::MX)
};



VECTOR_REPR(casadi::MX)
VECTOR_REPR(std::vector<casadi::MX>)

