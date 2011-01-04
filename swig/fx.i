%{
#include "casadi/fx/function_io.hpp"
#include "casadi/fx/fx.hpp"
#include "casadi/fx/jacobian.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/fx/integrator.hpp"
#include "casadi/fx/simulator.hpp"
%}

%include "casadi/fx/function_io.hpp"
%include "casadi/fx/fx.hpp"
%include "casadi/fx/jacobian.hpp"
%include "casadi/fx/sx_function.hpp"
%include "casadi/fx/mx_function.hpp"
%include "casadi/fx/linear_solver.hpp"
%include "casadi/fx/integrator_jacobian.hpp"
%include "casadi/fx/integrator.hpp"
%include "casadi/fx/simulator.hpp"

%template(vector_integrator) std::vector<CasADi::Integrator>;

/*#ifdef WITH_NUMPY*/
namespace CasADi {
%extend FunctionIO {
  %pythoncode %{
  def getArray(self,dir=0):
    import numpy as n
    v = self.data(dir)
    r = n.array((),dtype=float)
    r.resize(self.size1(),self.size2())
    for i in range(self.size1()):  # loop over rows
      for el in range(self.rowind_[i],self.rowind_[i+1]): # loop over the non-zero elements
        j=self.col_[el]  # column
        r[i,j] = v[el] # add the non-zero element

    return r
  %}
      
  %pythoncode %{
  def getMatrix(self,dir=0):
    import numpy as n
    return n.matrix(self.getArray(dir))
  %}
} // %extend FX
} // namespace CasADi
// #endif // WITH_NUMPY
