#ifndef CPLEX_INTERNAL_HPP
#define CPLEX_INTERNAL_HPP

#include "cplex/cplex_interface.h"
#include "ilcplex/cplex.h"
#include "casadi/fx/nlp_solver_internal.hpp"

namespace CasADi{

/** \brief CplexMatrix is a class used to convert CasADi matrices to CPLEX format (similar to CSC).
  The class definition can be found in cplex_internal.cpp.
  \author Carlo Savorgnan
  \date 2011
*/
class CplexMatrix;

class CplexInternal : public NLPSolverInternal{
  // TODO comment me!!!!
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
  /// Hessian of the Lagrangian function (used for format conversion)
  CplexMatrix H_mat_;
  /// Jacobian of the constraint function (used for format conversion)
  CplexMatrix J_mat_; 
  
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
