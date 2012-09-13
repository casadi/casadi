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
#include "acado_forward_declarations.hpp"
#include "acado_function.hpp"
#include "casadi/fx/integrator_internal.hpp"
#include "casadi/fx/linear_solver.hpp"
#include "casadi/fx/fx_internal.hpp"

namespace CasADi{
  
/**
@copydoc AcadoIntegrator_doc
*/
class AcadoIntegratorInternal : public IntegratorInternal{
  friend class AcadoIntegrator;

  public:
  
  /** \brief  Constructor */
  explicit AcadoIntegratorInternal(const FX& f, const FX& g);

  /** \brief  Clone */
  virtual AcadoIntegratorInternal* clone() const;
  
  /** \brief  Create a new integrator */
  virtual AcadoIntegratorInternal* create(const FX& f, const FX& g) const{ return new AcadoIntegratorInternal(f,g);}

  /** \brief  Destructor */
  virtual ~AcadoIntegratorInternal();

  /** \brief  Set all pointers to NULL */
  void setNull();
  
  /** \brief  Deallocate memory */
  void freeMem();
  
  /** \brief  Initialize */
  virtual void init();

  /** \brief  Reset the forward problem and bring the time back to t0 */
  virtual void reset();

  /** \brief  Reset the backward problem and take time to tf */
  virtual void resetB(){}

  /** \brief  Integrate forward until a specified time point */
  virtual void integrate(double t_out);

  /** \brief  Integrate backward until a specified time point */
  virtual void integrateB(double t_out){}

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

    // Dimensions
    int nt_, nxd_, nxa_;
        
    /** \brief  Public pointers to ACADO classes  */
    ACADO::TIME                 *t_;
    ACADO::DifferentialState    *xd_;
    ACADO::AlgebraicState       *xa_;
    ACADO::Parameter            *p_;
    ACADO::IntermediateState    *arg_;

    ACADO::DifferentialEquation *diff_eq_;
    ACADO::IntegratorBDF        *integrator_;

    // DAE rhs
    AcadoFunction rhs_;

    // Integration grid
    ACADO::Grid                 *interval_;
    ACADO::VariablesGrid        *differentialStates_;
    ACADO::VariablesGrid        *algebraicStates_;
    ACADO::Vector               *tmp_;
    
    // Seed/sensitivity vectors
    ACADO::Vector               *x_tmp_, *p_tmp_;
    
    // Number of interpolation time points
    int num_grid_points_;
    
    // Has the system been integrated once after the last initialization?
    bool has_been_integrated_;
    
    // Derivatives
    int nfsens_, nasens_;
    
    // Grid frozen
    bool frozen_grid_;
};


} // namespace CasADi

#endif //ACADO_INTEGRATOR_INTERNAL_HPP

