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
        r[i*self.size2()+j] = v[el] # add the non-zero element

    return r
  %}
      
  %pythoncode %{
  def getMatrix(self,dir=0):
    import numpy as n
    return n.matrix(getArray(self,dir))
  %}
} // %extend FX
} // namespace CasADi
// #endif // WITH_NUMPY


#if 0
namespace CasADi {
%extend FX {
  %pythoncode %{
  def getOutput(self,ind=0):
   """ A wrapper around getOutputData, that outputs as a numpy.ndarray """
   import numpy as n
   return n.reshape(n.array(self.getOutputData(ind)),(self.output(ind).size1(),self.output(ind).size2()))
  %}
  
  %pythoncode %{
  def getInput(self,ind=0):
   """ A wrapper around getInputData, that outputs as a numpy.ndarray """
   import numpy as n
   return n.reshape(n.array(self.getInputData(ind)),(self.input(ind).size1(),self.input(ind).size2()))
  %}
  
  %pythoncode %{
  def setInput(self,num,ind=0):
   """ A wrapper around setInput, that allows list,tuples and numpy.arrays """
   import numpy as n
   if type(num)==type([]) or type(num)==type((1,)) :
    self.setInputData(num,ind)
   else:
    temp=n.array(num)
    if len(temp.shape)>1 and not(temp.shape[0]==self.input(ind).size1() and temp.shape[1]==self.input(ind).size2()):
      raise Exception("setInput dimension mismatch. You provided a non-vector matrix (%d,%d), but the dimensions don't match with (%d,%d). " % (temp.shape[0],temp.shape[1],self.input(ind).size1(),self.input(ind).size2()))
    self.setInputData(temp.flatten().tolist(),ind)
  %}
  
  %pythoncode %{
  def __call__(self,x,iind=0,oind=0):
   """ setInput, evaluate and getOutput combined.
   
     def __call__(self,x,iind=0,oind=0)
     
     When x is numeric, does the same as:
      if (self.isInit()):
       self.init()
      self.setInput(x,iind)
      self.evaluate()
      return n.array(self.getOutput(oind))
   
   """
   if isinstance(x[0],MX):
    return self.call(x,iind)
   else:
    import numpy as n
    if not(self.isInit()):
      self.init()
    self.setInput(x,iind)
    self.evaluate()
    return n.array(self.getOutput(oind))
  %}
} // %extend FX
} // namespace CasADi

#endif
