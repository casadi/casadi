#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class CplexInternal;
  
/** \brief Interface to CPLEX solver.
  @copydoc NLPSolver_doc
  Attention! The interface is not complete yet.
  Also if a quadratic term can be set with this interface, it is ignored!
  \author Carlo Savorgnan
  \date 2011
*/
class CplexSolver : public NLPSolver {
  // TODO comment me!!!!
  public:
    /// Default constructor
    CplexSolver();
    
    /// Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit CplexSolver(const FX& F,         /**< F objective function */
                         const FX& G = FX(),  /**< constraint function (default only bound constraints) */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory). NOT USED*/
                         const FX& J = FX(),  /**< Jacobian of G (default -> differentiate) */
                         const FX& GF = FX()  /**< Gradient of the objective function (default: adjoint mode AD on F) */
                        );

    /// Access functions of the node
    CplexInternal* operator->();
    const CplexInternal* operator->() const;
    
    /// Set CPLEX integer parameters
    void setIntParam(const std::string& name, int val);
    
    /// Set CPLEX double parameters
    void setDoubleParam(const std::string& name, double val);

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    
};

} // namespace CasADi

#endif //CPLEX_SOLVER_HPP
