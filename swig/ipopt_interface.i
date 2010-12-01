%{
#include "ipopt_interface/ipopt_solver.hpp"
%}

namespace CasADi{

// Inputs of an NLP Solver
enum NLPInput{NLP_X_INIT,NLP_LBX,NLP_UBX,NLP_LBG,NLP_UBG,NLP_LAMBDA_INIT,NLP_NUM_IN};

// Outputs of an NLP Solver
enum NLPOutput{NLP_X_OPT,NLP_COST,NLP_LAMBDA_OPT,NLP_LAMBDA_LBX,NLP_LAMBDA_UBX,NLP_NUM_OUT};

class IpoptSolver : public FX{
  public:
    explicit IpoptSolver(const FX& F,         /**< F objective function */
                         const FX& G = FX(),  /**< constraint function (default only bound constraints) */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory) */
                         const FX& J = FX()   /**< Jacobian of G (default -> differentiate) */
    );
    
};

} // namespace CasADi


