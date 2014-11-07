/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
#ifndef CASADI_MX_I
#define CASADI_MX_I

%include <casadi/core/shared_object.i>
%include <casadi/core/matrix/matrix.i>
%include <casadi/core/matrix/generic_expression.i>

%include <casadi/core/mx/mx.hpp>

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
    return self.nz.__iter__()
    
  %}
  #endif //SWIGPYTHON
  
  binopsrFull(casadi::MX)
};

VECTOR_REPR(casadi::MX)
VECTOR_REPR(std::vector<casadi::MX>)

#endif // CASADI_MX_I
