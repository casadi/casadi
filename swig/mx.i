%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%include "casadi/mx/mx.hpp"




%template(sparsity_vector) std::vector<CasADi::CRSSparsity>;

#ifdef SWIGPYTHON
%{
template<> swig_type_info** meta< std::pair< CasADi::MX, std::vector< CasADi::MX> > >::name = &SWIGTYPE_p_std__pairT_CasADi__MX_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t_t;
%}

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



%include "casadi/mx/mx_tools.hpp"

