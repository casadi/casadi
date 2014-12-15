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
#ifndef CASADI_SX_ELEMENT_I
#define CASADI_SX_ELEMENT_I

%include <casadi/core/sx/sx_element.hpp>

#ifdef SWIGPYTHON
%extend casadi::Sparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
            
        @property
        def T(self):
            return _casadi.transpose(self)
            
        def __array__(self,*args,**kwargs):
            return DMatrix(self,1).toArray()
    %}
};

#endif // SWIGPYTHON

VECTOR_REPR(casadi::Matrix<casadi::SXElement>)

#ifdef SWIGPYTHON
%pythoncode %{
try:
  import numpy
  def constpow(x,y):
    pass

  constpow=numpy.frompyfunc(constpow,2,1)
  
  def fmin(x,y):
    pass
    
  def fmax(x,y):
    pass
  
  _min_ufunc = numpy.frompyfunc(fmin,2,1)
  _max_ufunc = numpy.frompyfunc(fmax,2,1)
  
  _defaultmin = min
  def min(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmin') or hasattr(args[1],'fmin')):
      return _min_ufunc(*args)
    else:
      return _defaultmin(*args,**kwargs)
      
  _defaultmax = max
  def max(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmax') or hasattr(args[1],'fmax')):
      return _max_ufunc(*args)
    else:
      return _defaultmax(*args,**kwargs)
except:
  pass
%}
#endif // SWIGPYTHON

namespace casadi {
%extend Matrix<SXElement>{
    
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<casadi::SXElement>)
       
    #ifdef SWIGPYTHON
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for j in range(self.size2()):
        for el in range(self.colind(j),self.colind(j+1)):
          i=self.row(el)
          r[i,j] = self.nz[el]
      return r
    %}
    
  %python_array_wrappers(1001.0)
  #endif // SWIGPYTHON 
  
#ifdef SWIGPYTHON
  binopsrFull(casadi::Matrix<casadi::SXElement>)  
#endif // SWIGPYTHON

};

} // namespace casadi

#ifdef SWIGPYTHON
#include <arrayobject.h>

// Template instantiations
%template()    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SX)             casadi::Matrix<casadi::SXElement>;

%extend casadi::Matrix<casadi::SXElement> {
   %template(SX) Matrix<int>;
   %template(SX) Matrix<double>;
};

#endif // CASADI_SX_ELEMENT_I
