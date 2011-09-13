/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADO_INTEGRATOR_INTERNAL_HPP
#define ACADO_INTEGRATOR_INTERNAL_HPP

#include "acado_integrator.hpp"
#include "casadi/fx/integrator_internal.hpp"
#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  
/**
@copydoc AcadoIntegrator_doc
*/
class AcadoIntegratorInternal : public IntegratorInternal{
  friend class AcadoIntegrator;

  public:
  
  /** \brief  Constructor */
  explicit AcadoIntegratorInternal(const FX& f, const FX& q);

  /** \brief  Clone */
  virtual AcadoIntegratorInternal* clone() const;
  
  /** \brief  Create a new integrator */
  virtual AcadoIntegratorInternal* create(const FX& f, const FX& q) const{ return new AcadoIntegratorInternal(f,q);}

  /** \brief  Destructor */
  virtual ~AcadoIntegratorInternal();

  /** \brief  Initialize */
  virtual void init();

  /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int nfsens=0, int nasens=0){}

  /** \brief  Reset the solver of the adjoint problem and take time to tf */
  virtual void resetAdj(){}

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out){}

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out){}

  /** \brief  Set the stop time of the forward integration */
  virtual void setStopTime(double tf){}
  
  /** \brief  Print solver statistics */  
  virtual void printStats(std::ostream &stream) const{}

  /** \brief Set linear solver */
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac){}
  
  /** \brief Get the Jacobian */
  virtual FX getJacobian(){}

  /** \brief Get the Linear solver */
  virtual LinearSolver getLinearSolver(){}

};


} // namespace CasADi

#endif //ACADO_INTEGRATOR_INTERNAL_HPP

