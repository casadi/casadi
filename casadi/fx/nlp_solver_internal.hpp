#ifndef NLP_SOLVER_INTERNAL_HPP
#define NLP_SOLVER_INTERNAL_HPP

#include "nlp_solver.hpp"

namespace CasADi{
    
/** \brief NLP solver storage class
  \author Joel Andersson 
  \date 2010
*/
class NLPSolverInternal : public FXNode{

public:
  explicit NLPSolverInternal(const FX& F, const FX& G, const FX& H, const FX& J);
  virtual ~NLPSolverInternal() = 0;

  virtual void init();

protected:
  int n_,m_;

  // The NLP functions
  /// objective function
  FX F_;
  /// constraint function
  FX G_; 
  /// Hessian of the Lagrangian function
  FX H_;
  /// Jacobian of the constraint function
  FX J_; 

};

} // namespace CasADi

#endif //NLP_SOLVER_INTERNAL_HPP
