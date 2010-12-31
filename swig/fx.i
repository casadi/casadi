%{
#include "casadi/fx/fx.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/fx/integrator.hpp"
#include "casadi/fx/simulator.hpp"
%}

namespace CasADi {

class FunctionIO {
  public:
    int size1() const;
    int size2() const;
    int numel() const;
};

class FX : public OptionsFunctionality{
  public:
  void init();
  void evaluate(int fsens_order=0, int asens_order=0);
  void solve();
  int getNumInputs() const;
  int getNumOutputs() const;
  FX jacobian(int iind=0, int oind=0);
  FX hessian(int iind=0, int oind=0);
  
  const FunctionIO & 	input (int i=0) const;
  const FunctionIO & 	output (int i=0) const;

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
     
      self.setInput(x,iind)
      self.evaluate()
      return n.array(self.getOutput(oind))
   
   """
   if isinstance(x[0],MX):
    return self.call(x,iind)
   else:
    import numpy as n
    self.setInput(x,iind)
    self.evaluate()
    return n.array(self.getOutput(oind))
  %}
  
  
  %rename(setInputData) setInput;
  
  // Forward renaming declarations
  %rename(getFwdSeed) getFwdSeedData;
  %rename(getFwdSens) getFwdSensData;
  %rename(getAdjSeed) getAdjSeedData;
  %rename(getAdjSens) getAdjSensData;

  // Get the data structure
  const std::vector<double>& getInputData(int ind=0) const;
  const std::vector<double>& getOutputData(int ind=0) const;
  const std::vector<double>& getFwdSeedData(int ind=0, int dir=0) const;
  const std::vector<double>& getFwdSensData(int ind=0, int dir=0) const;
  const std::vector<double>& getAdjSeedData(int ind=0, int dir=0) const;
  const std::vector<double>& getAdjSensData(int ind=0, int dir=0) const;

  // Set a value
#define SETTERS(T)\
  void setInput(T val, int ind=0); \
  void setOutput(T val, int ind=0); \
  void setFwdSeed(T val, int ind=0, int dir=0); \
  void setFwdSens(T val, int ind=0, int dir=0); \
  void setAdjSeed(T val, int ind=0, int dir=0); \
  void setAdjSens(T val, int ind=0, int dir=0);
SETTERS(double)
SETTERS(const std::vector<double>&)
#undef SETTERS

};



  
  
%extend FX{
  MX call(const std::vector<MX> &x, int ind=0) const{
    return $self->operator()(x,ind);
}


};

class SXFunction : public FX{
public:
  /** \brief  Multiple (vector valued) input, multiple (vector valued) output */
  SXFunction(const std::vector< std::vector<SX> >& arg, const std::vector< std::vector<SX> >& res);

  /** \brief  Multiple (matrix valued) input, multiple (matrix valued) output */
  SXFunction(const std::vector< SXMatrix>& arg, const std::vector<SXMatrix>& res);

  /** \brief get the input arguments symbolically */
  SXMatrix getArgumentIn(int iind=0) const;

  /** \brief get the output arguments symbolically */
  SXMatrix getArgumentOut(int iind=0) const;

  /** \brief  evaluate symbolically */
  std::vector<SXMatrix> eval(const std::vector<SXMatrix>& arg);

  /** \brief  evaluate symbolically, non-zero entries */
  std::vector< std::vector<SX> > eval(const std::vector< std::vector<SX> >& arg);
  
  /** \brief Jacobian of output oind with respect to input iind */
  SXFunction jacobian(int iind=0, int oind=0);
  
  /** \brief Hessian of output oind with respect to input iind */
  SXFunction hessian(int iind=0, int oind=0);

  /// Jacobian via source code transformation
  SXMatrix jac(int iind=0, int oind=0);

  /// Gradient via source code transformation
  SXMatrix grad(int iind=0, int oind=0);
  
  /// Hessian (forward over adjoint) via source code transformation
  SXMatrix hess(int iind=0, int oind=0);

};

class MXFunction : public FX{
public:
/** \brief  Single input, single output */
  MXFunction(const MX& input, const MX& output);

/** \brief  Single input, multiple output */
  MXFunction(const MX& input, const std::vector<MX>& output);

/** \brief  Multiple input, single output */
  MXFunction(const std::vector<MX>& input, const MX& output);

/** \brief  Multiple input, multiple output*/
  MXFunction(const std::vector<MX>& input, const std::vector<MX>& output);
};

} // namespace CasADi

%include "casadi/fx/linear_solver.hpp"
%include "casadi/fx/integrator_jacobian.hpp"
%include "casadi/fx/integrator.hpp"
%include "casadi/fx/simulator.hpp"
