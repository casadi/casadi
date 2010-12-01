#ifndef IPOPT_SOLVER_HPP
#define IPOPT_SOLVER_HPP

#include "casadi/fx/nlp_solver.hpp"

namespace CasADi{
  
class IpoptInternal;
  
class IpoptSolver : public NLPSolver {
  public:
    /// Default constructor
    IpoptSolver();
    
    /// Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit IpoptSolver(const FX& F,         /**< F objective function */
                         const FX& G = FX(),  /**< constraint function (default only bound constraints) */
                         const FX& H = FX(),  /**< Hessian of the lagrangian function (default: limited memory) */
                         const FX& J = FX()   /**< Jacobian of G (default -> differentiate) */
                        );

    /// Access functions of the node
    IpoptInternal* operator->();
    const IpoptInternal* operator->() const;

    /// Assert that the node is pointing to the right type of object
    void assertNode() const;

    
};

} // namespace CasADi

#endif //IPOPT_SOLVER_HPP
