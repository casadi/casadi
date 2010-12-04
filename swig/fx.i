%{
#include "casadi/fx/fx.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/fx/integrator.hpp"
#include "casadi/fx/simulator.hpp"
%}

namespace CasADi {

class FX : public OptionsFunctionality{
  public:
  void init();
  void evaluate(int fsens_order=0, int asens_order=0);
  void solve();
  int getNumInputs() const;
  int getNumOutputs() const;
  FX jacobian(int iind=0, int oind=0);
  FX hessian(int iind=0, int oind=0);

  // Forward renaming declarations
  %rename(getInput) getInputData;
  %rename(getOutput) getOutputData;
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

/// Public class
class LinearSolver : public FX{
public:
  /// Factorize the matrix / prepare the solution
  void prepare();
  
  /// Solve the system of equations
  void solve();
};

// Input arguments of an integrator 
enum IntegratorInput{INTEGRATOR_T0, INTEGRATOR_TF, INTEGRATOR_X0, INTEGRATOR_P, INTEGRATOR_XP0, INTEGRATOR_NUM_IN};

// Output arguments of an integrator 
enum IntegratorOutput{INTEGRATOR_XF, INTEGRATOR_XPF, INTEGRATOR_NUM_OUT};

// Input arguments of an explicit ODE right hand side 
enum ODEInput{ODE_T, ODE_Y, ODE_P, ODE_NUM_IN};

// Output arguments of an explicit ODE right hand side 
enum ODEOutput{ODE_RHS, ODE_NUM_OUT};

// Input arguments of an DAE residual function 
enum DAEInput{DAE_T, DAE_Y, DAE_YDOT, DAE_P, DAE_NUM_IN};

// Output arguments of an DAE residual function
enum DAEOutput{DAE_RES, DAE_NUM_OUT};

//Public class 
class Integrator : public FX{
public:
  // Print solver statistics 
  void printStats(std::ostream &stream=std::cout) const;
  
  // Reset the solver and bring the time back to t0 
  void reset(int fsens_order, int asens_order=0);

  //Integrate until a specified time point 
  void integrate(double t_out);

  // Set a stop time for the forward integration 
  void setStopTime(double tf);  
};

// Indices of the inputs of the output function
enum OutputInput{OUTPUT_T, OUTPUT_X, OUTPUT_P, OUTPUT_NUM_IN};

// Indices of the inputs of the function
enum SimulatorInput{SIMULATOR_X0, SIMULATOR_P, SIMULATOR_NUM_IN};

class Simulator : public FX{
public:
  
  // Constructor
  Simulator(const Integrator& integrator, const FX& output_fcn, const std::vector<double>& grid);
  
  // Output function equal to the state
  Simulator(const Integrator& integrator, const std::vector<double>& grid);
};








} // namespace CasADi

