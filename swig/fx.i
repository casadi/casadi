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
   import numpy as n
   return n.reshape(n.array(self.getOutputData(ind)),(self.output(ind).size1(),self.output(ind).size2()))
  %}
  
  %pythoncode %{
  def getInput(self,ind=0):
   import numpy as n
   return n.reshape(n.array(self.getInputData(ind)),(self.input(ind).size1(),self.input(ind).size2()))
  %}
  
  
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
  MX __call__(const std::vector<MX> &x, int ind=0) const{
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
