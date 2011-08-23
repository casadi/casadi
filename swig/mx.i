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
  
  %python_array_wrappers(1002.0)
  
  %pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")
  %}
  #endif //SWIGPYTHON
  
  binopsrFull(CasADi::MX)
};

%my_generic_const_typemap(PRECEDENCE_MXVector,std::vector< CasADi::MX >);
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DMatrixVector,std::vector< CasADi::Matrix<double> >);
#endif // SWIGPYTHON

%include "casadi/mx/mx_tools.hpp"

