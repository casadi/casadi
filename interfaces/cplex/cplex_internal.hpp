#ifndef CPLEX_INTERNAL_HPP
#define CPLEX_INTERNAL_HPP

#include "cplex/cplex_interface.h"
#include "ilcplex/cplex.h"
#include "casadi/fx/nlp_solver_internal.hpp"

namespace CasADi{
    
class CplexInternal : public NLPSolverInternal{

public:
  explicit CplexInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF);
  virtual ~CplexInternal();

  virtual void init();
  virtual void evaluate(int fsens_order, int asens_order);
  
  /// objective function
  FX F_;
  /// constraint function
  FX G_; 
  /// Hessian of the Lagrangian function
  FX H_;
  /// Jacobian of the constraint function
  FX J_; 
  /// Gradient of the objective function
  FX GF_; 
  
  // CPLEX environment pointer
  CPXENVptr env_;
  // CPLEX lp pointer
  CPXLPptr lp_;
  
  // CPLEX double parameter
  std::map<std::string, double> double_param_;
  
  // CPLEX int parameter
  std::map<std::string, int> int_param_;
  
  // sense of the optimization (min or max)
  int sense_;
};

} // namespace CasADi

#endif //CPLEX_INTERNAL_HPP
