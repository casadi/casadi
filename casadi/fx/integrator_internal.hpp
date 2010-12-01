#ifndef INTEGRATOR_INTERNAL_HPP
#define INTEGRATOR_INTERNAL_HPP

#include "integrator.hpp"

namespace CasADi{

/** \brief Internal storage for integrator related data
  \author Joel Andersson 
  \date 2010
*/
class IntegratorInternal : public FXNode{
public:
  /** \brief  Constructor */
  IntegratorInternal(int nx, int np);

  /** \brief  Destructor */
  virtual ~IntegratorInternal()=0;
  
  /** \brief  Print solver statistics */
  virtual void printStats(std::ostream &stream) const = 0;

    /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int fsens_order, int asens_order) = 0;

    /** \brief  Reset the solver of the adjoint problem and take time to tf */
  virtual void resetAdj() = 0;

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out) = 0;

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out) = 0;

  /** \brief  Set stop time for the integration */
  virtual void setStopTime(double tf) = 0;
  
  /** \brief  evaluate */
  virtual void evaluate(int fsens_order, int asens_order);

  /** \brief  Initialize */
  virtual void init();

  
  /** Lenght of the state vector 
   (also includes "states" evaluated from quadrature formulas)
  */
  int nx_;
  
  /// Number of parameters
  int np_;
  
  /// Current time
  double t_;

  //@{
  /// options
  bool exact_jacobian_;
  double abstol_, reltol_;
  double fsens_abstol_, fsens_reltol_;
  double asens_abstol_, asens_reltol_;
  int max_num_steps_;
  bool finite_difference_fsens_;  
  //@}
};
  
} // namespace CasADi

#endif // INTEGRATOR_INTERNAL_HPP
