#ifndef IDAS_INTEGRATOR_HPP
#define IDAS_INTEGRATOR_HPP

#include "casadi/fx/integrator.hpp"

/** Integrator where the differential equation is assumed to be in the form:

  index-1 DAE:
  f(t,y,der(y),p) == 0

  pure quadratures:
  der(q) = g(t,y,p)
  
  The state vector is x = [x,q]
  
*/

namespace CasADi{
namespace Sundials{

/** \brief  Forward declaration of internal class */
class IdasInternal;

/** \brief  Input arguments of a jacobian function: J = df/dy + cj*df/dydot */
enum JACInput{JAC_T, JAC_Y, JAC_YDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/** \brief  Output arguments of an DAE residual function */
enum JACOutput{JAC_J, JAC_NUM_OUT};

/** \brief  Public class */
class IdasIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  IdasIntegrator();
  
  /** \brief  Create an integrator for explicit ODEs */
  explicit IdasIntegrator(const FX& f, const FX& q=FX(), const FX& jacx=FX(), const FX& jacp=FX());

  /** \brief  Access functions of the node */
  IdasInternal* operator->();
  const IdasInternal* operator->() const;
  
};


} // namespace Sundials
} // namespace CasADi

#endif //IDAS_INTEGRATOR_HPP

