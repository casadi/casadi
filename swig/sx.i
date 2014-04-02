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
#include "symbolic/matrix/sparsity.hpp"
#include "symbolic/matrix/slice.hpp"
#include "symbolic/matrix/generic_expression.hpp"
#include "symbolic/matrix/generic_matrix.hpp"
#include "symbolic/matrix/matrix.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_element.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/mx/mx.hpp"

%}

// http://www.gnu.org/software/octave/doc/interpreter/Operator-Overloading.html#Operator-Overloading
// http://mentat.za.net/DaCodaAlFine.pdf
// "dispatch binary operator" can be found in octrun.swg: look for dispatch_binary_op in swig generated


#ifndef SWIGXML
%include "typemaps.i"
#endif

%include "symbolic/matrix/sparsity.hpp"
%include "symbolic/matrix/slice.hpp"

%include "symbolic/matrix/generic_expression.hpp"
%template(ExpIMatrix)        CasADi::GenericExpression<CasADi::Matrix<int> >;
%template(ExpDMatrix)        CasADi::GenericExpression<CasADi::Matrix<double> >;
%template(ExpSX)       CasADi::GenericExpression<CasADi::Matrix<CasADi::SXElement> >;
%template(ExpMX)             CasADi::GenericExpression<CasADi::MX>;
%template(ExpSXElement)             CasADi::GenericExpression<CasADi::SXElement>;

%include "symbolic/matrix/generic_matrix.hpp"
%template(GenIMatrix)        CasADi::GenericMatrix<CasADi::Matrix<int> >;
%template(GenDMatrix)        CasADi::GenericMatrix<CasADi::Matrix<double> >;
%template(GenSX)       CasADi::GenericMatrix<CasADi::Matrix<CasADi::SXElement> >;
%template(GenMX)             CasADi::GenericMatrix<CasADi::MX>;

%include "symbolic/matrix/matrix.hpp"
%template(IMatrix)           CasADi::Matrix<int>;
%template(DMatrix)           CasADi::Matrix<double>;

%extend CasADi::Matrix<double> {
   %template(DMatrix) Matrix<int>;
};

%include "symbolic/sx/sx_element.hpp"



#ifdef SWIGPYTHON
%extend CasADi::Sparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
            
        @property
        def T(self):
            return self.transpose()
            
        def __array__(self,*args,**kwargs):
            return DMatrix(self,1).toArray()
    %}
};

#endif // SWIGPYTHON

VECTOR_REPR(CasADi::SXElement)
VECTOR_REPR(std::vector<CasADi::SXElement>)
VECTOR_REPR(CasADi::Matrix<CasADi::SXElement>)

#ifdef SWIGPYTHON
%pythoncode %{
try:
  # This following code attempts to make constpow to work.
  # We use numpy's ufunc functionality for that.
  #
  #   myfunc = numpy.frompyfunc(myfunction,2,1)    creates a ufunc.
  #
  #   What happens when we do myfunc(x,y)?
  #
  #     The ufunc will call myfunction(x.__array__(),y.__array__())
  #        It is an error if __array__ does not return a numpy.array/numpy.matrix
  #        This is problematic since e.g. MX is not expandable in a numpy matrix
  #        We let __array__() return a dummy 1x1 numpy.array.
  #        Since the arguments to myfunction are dummy, there is no reason to implement myfunction.
  #     
  #     Great, so we have a ufunctor that takes its arguments, makes dummy values of it
  #     and feeds those dummies to a function that does nothing whatsoever with them?
  #     Hold on, this is not the end of the story.
  #
  #     The ufunc functionality has a provision to post-process the result that was made from myfunction.
  #       It first checks which of the arguments x, y has the largest __array_priority__.
  #       We made it such that casadi types always have the highest priority.
  #       From the argument with the largest priority, the __array_wrap__ member is called.
  #       The __array__wrap__ method can have a 'context' argument, which gives information about the ufunctor.
  #       What we do is check the name of the ufunctor and try to call x.name(y) x.__name__(y) y.__rname__(x)
  #       It's the result of this call that is finally returned to the user.
  #
  #  Why don't you just let __array__() return a 1x1 object numpy.array with the casadi
  #  object as its only element, and implement myfunction so that it unwraps from the 1x1
  #  and returns a 1x1 wrapped result?
  #
  #  Good question. The problem is broadcasting.
  #  Consider myfunc(array([[4,5],[6,7]]),MX("x",2,2)).
  #   Because the y argument is wrapped in a 1x1, numpy will expand it to 2x2 to match the dimensions of argument x.
  #   The result will be 4 calls to myfunction: myfunction(4,MX("x,2,2")),  myfunction(4,MX("x,2,2"))
  #
  #  Why do you use ufunc in the first places? What's wrong with just global functions?       
  #
  #  Good point. The problem is user-friendliness.
  #  We don't want a mere 'from numpy import *' statement to destroy functionality of x**y for example.
  #  Come to think of it, this is not really a terrific argument.
  #  We might as well have an __array__() that warns the user if he should have overloaded the global functions.
  #
  #     In fact, we should refactor and forget about this awful mess.
  
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

namespace CasADi {


%extend SXElement {
#ifdef SWIGPYTHON

  %python_array_wrappers(1000.0)
  
  %pythoncode %{

  def toArray(self):
    import numpy as n
    r = n.array((),dtype=object)
    r.resize(1,1)
    r[0,0] = self
    return r
  %}
  
  

  #endif // SWIGPYTHON
  
  #ifdef SWIGOCTAVE
  std::vector<int> __dims__() const {
    std::vector<int> ret(2);
    ret[0] = 1;
    ret[1] = 1;
    return ret;
  }
  
  #endif // SWIGOCTAVE
  
  binopsrFull(CasADi::SXElement)
  // a+b when a is SXElement, b is numpy.array. __array_priority works, but does not suffice to yield implicit casting
  binopsFull(const CasADi::Matrix<CasADi::SXElement> & b,,CasADi::Matrix<CasADi::SXElement>,CasADi::Matrix<CasADi::SXElement>)

};



%extend Matrix<SXElement>{
    
    %matrix_convertors
    %matrix_helpers(CasADi::Matrix<CasADi::SXElement>)
       
    #ifdef SWIGPYTHON
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for j in range(self.size2()):  # loop over columns
        for el in range(self.colind(j),self.colind(j+1)): # loop over the non-zero elements
          i=self.row(el)  # column
          r[i,j] = self.at(el) # add the non-zero element

      return r
    %}
    
  %python_array_wrappers(1001.0)
  #endif // SWIGPYTHON 
  
  binopsrFull(CasADi::Matrix<CasADi::SXElement>)  
};
 


} // namespace CasADi

#ifdef SWIGPYTHON
#ifdef WITH_NUMPY
#include <arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template()    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SX)             CasADi::Matrix<CasADi::SXElement>;

%extend CasADi::Matrix<CasADi::SXElement> {
   %template(SX) Matrix<int>;
   %template(SX) Matrix<double>;
};

