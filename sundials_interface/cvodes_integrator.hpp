#ifndef CVODES_INTEGRATOR_HPP
#define CVODES_INTEGRATOR_HPP

#include "casadi/fx/integrator.hpp"

/** Function that integrates the ODE:

  ydot == f(t,y,p)
  from t0 to tf
  
  given the initial condition
  y(t0) == y0;
  
*/


namespace CasADi{
namespace Sundials{
  
/** \brief  Forward declaration of internal class */
class CVodesInternal;

/** \brief  Public class */
class CVodesIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  CVodesIntegrator();
  
  /** \brief  Create an integrator for explicit ODEs */
  explicit CVodesIntegrator(const FX& f, const FX& q=FX());
  
  /** \brief  Access functions of the node */
  CVodesInternal* operator->();
  const CVodesInternal* operator->() const;
};


} // namespace Sundials
} // namespace CasADi

#endif //CVODES_INTEGRATOR_HPP

