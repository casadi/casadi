%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%my_generic_const_typemap(PRECEDENCE_MX,CasADi::MX);

%include "casadi/mx/mx.hpp"




%template(sparsity_vector) std::vector<CasADi::CRSSparsity>;

#ifdef SWIGPYTHON
%meta_pair(CasADi::MX, std::vector< CasADi::MX >)
%typemap(out) std::pair< CasADi::MX, std::vector< CasADi::MX >  > {
    bool ret = meta< std::pair< CasADi::MX, std::vector< CasADi::MX >  > >::toPython($1,$result);
    if (!ret) SWIG_exception_fail(SWIG_TypeError,"Could not convert to (MX,std::vector<MX>)");
}
#endif //SWIGPYTHON

%extend CasADi::MX{
  %python_matrix_helpers(CasADi::MX)
  #ifdef SWIGPYTHON
  %pythoncode %{
  __array_priority__ = 1002.0
  
  
  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
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
      
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      return n.array([1])
    else:
      raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")
  %}
  #endif //SWIGPYTHON
};

#ifdef SWIGPYTHON
%extend CasADi::MX{ 
  #binopsFull(double b,CasADi::MX,,CasADi::MX)
  #binopsFull(const CasADi::Matrix<double>& b,CasADi::MX,,CasADi::MX)
};
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
#endif // SWIGPYTHON

%include "casadi/mx/mx_tools.hpp"

