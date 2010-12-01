#ifndef ACADO_INTEGRATOR_HPP
#define ACADO_INTEGRATOR_HPP

#include "casadi/fx/integrator.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class ACADOIntegratorInternal;

/** \brief  Public class */
class ACADOIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  ACADOIntegrator();
  
  /** \brief  Create an integrator for semi-explicit ODEs/DAEs */
  explicit ACADOIntegrator(const FX& f);
  
  /** \brief  Access functions of the node */
  ACADOIntegratorInternal* operator->();
  const ACADOIntegratorInternal* operator->() const;
};


} // namespace CasADi

#endif //ACADO_INTEGRATOR_HPP

