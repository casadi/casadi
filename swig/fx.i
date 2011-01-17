%{
#include "casadi/fx/function_io.hpp"
#include "casadi/fx/fx.hpp"
#include "casadi/fx/jacobian.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/fx/integrator.hpp"
#include "casadi/fx/simulator.hpp"
#include "casadi/fx/nlp_solver.hpp"
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
%include "casadi/fx/nlp_solver.hpp"

%template(IntegratorVector) std::vector<CasADi::Integrator>;

/*#ifdef WITH_NUMPY*/
namespace CasADi {
%extend FunctionIO {
  %pythoncode %{
  def getArray(self):
    import numpy as n
    v = self.get()
    r = n.array((),dtype=float)
    r.resize(v.size1(),v.size2())
    for i in range(v.size1()):  # loop over rows
      for el in range(v.rowind(i),v.rowind(i+1)): # loop over the non-zero elements
        j=v.col(el)  # column
        r[i,j] = v[el] # add the non-zero element

    return r
  %}
      
  %pythoncode %{
  def getMatrix(self):
    import numpy as n
    return n.matrix(self.getArray())
  %}
} // %extend FX

%extend FX {


  %pythoncode %{
  def setInput(self,num,ind=0):
    """ A wrapper around setInput, that allows 2D numpy.arrays.
    Should be replaced by a typemap later on. """
    import numpy as n
    if (isinstance(num,n.ndarray)):
     temp=n.array(num)
     if len(temp.shape)>1 and not(temp.shape[0]==self.argument(ind).size1() and temp.shape[1]==self.argument(ind).size2()):
       raise Exception("setInput dimension mismatch. You provided a non-vector matrix (%d,%d), but the dimensions don't match with (%d,%d). " % (temp.shape[0],temp.shape[1],self.argument(ind).size1(),self.argument(ind).size2()))
     self.setInputOriginal(temp.flatten().tolist(),ind)
    else:
     self.setInputOriginal(num,ind)
  %}
} // %extend FX

} // namespace CasADi
// #endif // WITH_NUMPY
