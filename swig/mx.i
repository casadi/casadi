%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%include "casadi/mx/mx.hpp"




%template(sparsity_vector) std::vector<CasADi::CRSSparsity>;

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


