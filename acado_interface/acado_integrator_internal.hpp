#ifndef ACADO_INTEGRATOR_INTERNAL_HPP
#define ACADO_INTEGRATOR_INTERNAL_HPP

#include <acado_integrators.hpp>
#include "acado_integrator.hpp"
#include "casadi/fx/integrator_internal.hpp"

namespace CasADi{
  
class ACADOIntegratorInternal : public IntegratorInternal{
public:
  /** \brief  Constructor */
  explicit ACADOIntegratorInternal(const FX& f);

  /** \brief  Destructor */
  virtual ~ACADOIntegratorInternal();

  /** \brief  Initialize stage */
  virtual void init();

  /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int fsens_order, int asens_order);

  /** \brief  Reset the solver of the adjoint problem at the current time */
  virtual void resetAdj();

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out);

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out);

  /** \brief  Print integrator statistics */
  virtual void printStats(std::ostream &stream) const;

  /** \brief  ACADO integrators redefines the following function (should not be necessary) */
  virtual void evaluate(int fsens_order, int asens_order);

  virtual void setStopTime(double tf){}

  
  /// Function called by ACADO
  void ffcn( double *x, double *f);
  void ffcn_fwd(int number, double *x, double *seed, double *f, double *df);
  void ffcn_adj(int number, double *x, double *seed, double *f, double *df);

  /// Wrapper
  static void ffcn_wrapper( double *x, double *f, void *user_data );
  static void ffcn_fwd_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);
  static void ffcn_adj_wrapper(int number, double *x, double *seed, double *f, double *df, void *user_data);
  
  protected:
    
  // Auxiliary
  static int getNX(const FX& f); // count the total number of states
  static int getNP(const FX& f); // count the number of parameters

  // DAE rhs
  FX f_;

  ACADO::DifferentialState *x_;
  ACADO::Parameter *p_;
  ACADO::TIME *t_;
  ACADO::IntermediateState *augx_;
  ACADO::CFunction *model_;
  ACADO::DifferentialEquation *dae_;
  ACADO::Integrator *integrator_;

  // Seed vectors in acado format
  ACADO::Vector *x_temp_;
  ACADO::Vector *p_temp_;

  bool is_init_;
};


} // namespace CasADi

#endif //ACADO_INTEGRATOR_INTERNAL_HPP

