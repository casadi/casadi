#ifndef NLP_SOLVER_HPP
#define NLP_SOLVER_HPP

#include "../expression_tools.hpp"
#include "fx.hpp"

namespace CasADi{

  /// Inputs of an NLP Solver
  enum NLPInput{NLP_X_INIT,NLP_LBX,NLP_UBX,NLP_LBG,NLP_UBG,NLP_LAMBDA_INIT,NLP_NUM_IN};

  /// Outputs of an NLP Solver
  enum NLPOutput{NLP_X_OPT,NLP_COST,NLP_LAMBDA_OPT,NLP_LAMBDA_LBX,NLP_LAMBDA_UBX,NLP_NUM_OUT};

class NLPSolverInternal;

/** \brief NLPSolver
  \author Joel Andersson 
  \date 2010
*/
class NLPSolver : public FX{
  public:

  /// Default constructor
  NLPSolver();

  /// Access functions of the node
  NLPSolverInternal* operator->();
  const NLPSolverInternal* operator->() const;

  /// Assert that the node is pointing to the right type of object
  void assertNode() const;
};

} // namespace CasADi

#endif // NLP_SOLVER_HPP

