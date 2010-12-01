#ifndef ACADO_INTERNAL_HPP
#define ACADO_INTERNAL_HPP

#include "acado_interface.hpp"
#include "acado_function.hpp"

/** \brief  forward declarations */
namespace ACADO{
  class OCP;
  class Time;
  class DifferentialState;
  class DifferentialStateDerivative;
  class AlgebraicState;
  class Control;
  class Parameter;
  class DifferentialEquation;
  class IntermediateState;
  class OCP;
  class OptimizationAlgorithm;
}
namespace CasADi{

  
class AcadoInternal : public FXNode{
  friend class AcadoInterface;
  
  /** \brief  Constructor only accessable from the AcadoOCPSolver pointer class */
  explicit AcadoInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn);
  
  public:
    
    /** \brief  Destructor */
    virtual ~AcadoInternal();
    
    /** \brief  Initialize the solver */
    virtual void init();
    
    /** \brief  Solve the problem */
    virtual void evaluate(int fsens_order=0, int asens_order=0);
    
    // Dimensions
    int nt_, nxd_, nxa_, nu_, np_, nxdot_;
    int nx_;
    
    // Number of shooting nodes
    int n_nodes_;
    
    // Number of non-linear path constraints
    int nc_;
    
    // Number of initial value constraints
    int nr_;
        
    /** \brief  Public pointers to ACADO classes  */
    ACADO::TIME                   *t_;
    ACADO::DifferentialState     *xd_;
    ACADO::AlgebraicState        *xa_;
    ACADO::Control                *u_;
    ACADO::Parameter              *p_;
    ACADO::DifferentialStateDerivative *xdot_;
    ACADO::IntermediateState     *arg_;

    ACADO::DifferentialEquation   *f_;
    ACADO::OCP                   *ocp_;
    ACADO::OptimizationAlgorithm *algorithm_;

    // DAE rhs
    AcadoFunction ffcn_;

    // Meyer term
    AcadoFunction mfcn_;

    // Constraint function
    AcadoFunction cfcn_;
    
    // Initial equation
    AcadoFunction rfcn_;
    
    
};

} // namespace CasADi

#endif //ACADO_INTERNAL_HPP
