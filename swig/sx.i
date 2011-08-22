%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/slice.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"


%}

// http://www.gnu.org/software/octave/doc/interpreter/Operator-Overloading.html#Operator-Overloading
// http://mentat.za.net/DaCodaAlFine.pdf
// "dispatch binary operator" can be found in octrun.swg: look for dispatch_binary_op in swig generated

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/slice.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"


#ifdef SWIGPYTHON
%extend CasADi::CRSSparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
            
        def __array__(self,*args,**kwargs):
            return DMatrix(self,1).toArray()
    %}
};

#endif // SWIGPYTHON

%extend std::vector<CasADi::SX>{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<CasADi::Matrix< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<std::vector< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

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
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmin') or hasattr(args[1],'fmin')):
      return _max_ufunc(*args)
    else:
      return _defaultmax(*args,**kwargs)
except:
  pass
%}
#endif // SWIGPYTHON

namespace CasADi {


%extend SX {
#ifdef SWIGPYTHON
  %pythoncode %{
    def __lt__(self,other):
      return _casadi.__lt__(self,other)
    def __le__(self,other):
      return _casadi.__le__(self,other)
    def __eq__(self,other):
      return _casadi.__eq__(self,other)
    def __ne__(self,other):
      return _casadi.__ne__(self,other)
    def __gt__(self,other):
      return _casadi.__gt__(self,other)
    def __ge__(self,other):
      return _casadi.__ge__(self,other)
  %}
  
  
  %pythoncode %{
  __array_priority__ = 1000.0
  
  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    selfM = SXMatrix(self)
    if "vectorized" in name:
      name = name[:-len(" (vectorized)")]

    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
      args.reverse()
    if not(hasattr(selfM,name)):
      name = '__' + name + '__'
    fun=getattr(selfM, name)
    return fun(*args[1:])

  def toArray(self):
    import numpy as n
    r = n.array((),dtype=object)
    r.resize(1,1)
    r[0,0] = self
    return r
      
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      return n.array([1])
    else:
      return self.toArray()
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
  
  binopsFull(const Matrix<double>& b,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
  binopsFull(const CasADi::Matrix<CasADi::SX> & b,,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
    
};
  
%extend Matrix<SX>{
    // The constructor has to be added since SX::operator Matrix<SX does not work
    // Matrix<SX>(const SX&){ *$self
    #ifdef SWIGPYTHON
    %pythoncode %{
      def __lt__(self,other):
        return _casadi.__lt__(self,other)
      def __le__(self,other):
        return _casadi.__le__(self,other)
      def __eq__(self,other):
        return _casadi.__eq__(self,other)
      def __ne__(self,other):
        return _casadi.__ne__(self,other)
      def __gt__(self,other):
        return _casadi.__gt__(self,other)
      def __ge__(self,other):
        return _casadi.__ge__(self,other)
    %}
    #endif

    binopsFull(double b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    binopsFull(const CasADi::Matrix<double>& b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    #ifdef SWIGOCTAVE
    binopsFull(const CasADi::SX& b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    #endif // SWIGOCTAVE
  
    %python_matrix_convertors
    %python_matrix_helpers(CasADi::Matrix<CasADi::SX>)
       
    #ifdef SWIGPYTHON

    
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for i in range(self.size1()):  # loop over rows
        for el in range(self.rowind(i),self.rowind(i+1)): # loop over the non-zero elements
          j=self.col(el)  # column
          r[i,j] = self[el] # add the non-zero element

      return r
    %}
    
  %pythoncode %{
  __array_priority__ = 1001.0
  %}
    
  %pythoncode %{
  def __array_wrap__(self,out_arr,context=None):
    name = context[0].__name__
    args = list(context[1])
    if "vectorized" in name:
      name = name[:-len(" (vectorized)")]
    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
  %}

  %pythoncode %{
    def __array__(self,*args,**kwargs):
      import numpy as n
      if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
        return n.array([1])
      else:
        return self.toArray()
  %}
  #endif // SWIGPYTHON    
};



} // namespace CasADi

#ifdef SWIGPYTHON
#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SXMatrixVector)       std::vector<CasADi::Matrix<CasADi::SX> > ;
%template(SXMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;




%my_generic_const_typemap(PRECEDENCE_SX,CasADi::SX);
%my_generic_const_typemap(PRECEDENCE_SXMatrix,CasADi::Matrix<CasADi::SX>);
%my_generic_const_typemap(PRECEDENCE_SXMatrixVector,std::vector< CasADi::Matrix<CasADi::SX> >);

#ifdef SWIGOCTAVE
%meta_vector(CasADi::Matrix<CasADi::SX>);
#endif //SWIGOCTAVE
%meta_vector(std::vector<CasADi::SX>);
%meta_vector(CasADi::SX);

%template(SXVector)             std::vector<CasADi::SX>;
%template(SXVectorVector)       std::vector<std::vector<CasADi::SX> > ;
%template(SXVectorVectorVector) std::vector< std::vector<std::vector<CasADi::SX> > > ;
%template(SXMatrix)             CasADi::Matrix<CasADi::SX>;

