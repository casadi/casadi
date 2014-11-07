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

#ifndef SWIGMATLAB
%extend SXElement {
  double __float__() { return $self->getValue();}
  int __int__() { return $self->getIntValue();}
}
#endif

#ifdef SWIGPYTHON
%extend casadi::Sparsity{
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

VECTOR_REPR(casadi::SXElement)
VECTOR_REPR(std::vector<casadi::SXElement>)
VECTOR_REPR(casadi::Matrix<casadi::SXElement>)

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

namespace casadi {


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
  
  binopsrFull(casadi::SXElement)
  // a+b when a is SXElement, b is numpy.array. __array_priority works, but does not suffice to yield implicit casting
  binopsFull(const casadi::Matrix<casadi::SXElement> & b,,casadi::Matrix<casadi::SXElement>,casadi::Matrix<casadi::SXElement>)

};



%extend Matrix<SXElement>{
    
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<casadi::SXElement>)
       
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
  
  binopsrFull(casadi::Matrix<casadi::SXElement>)  
};
 


} // namespace casadi

#ifdef SWIGPYTHON
#ifdef WITH_NUMPY
#include <arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template()    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SX)             casadi::Matrix<casadi::SXElement>;

%extend casadi::Matrix<casadi::SXElement> {
   %template(SX) Matrix<int>;
   %template(SX) Matrix<double>;
};

#endif // CASADI_SX_ELEMENT_I
