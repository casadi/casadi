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

#ifndef MULTIPLE_SHOOTING_HPP
#define MULTIPLE_SHOOTING_HPP

#include "../casadi/fx/ocp_solver.hpp"
#include "../casadi/fx/nlp_solver.hpp"

namespace CasADi{
  namespace OptimalControl{
  
    class MultipleShootingInternal;
    
  /** \brief Multiple Shooting
   *
   *   ns: Number of shooting nodes: from option number_of_grid_points\n
   *   nx: Number of differential states: from ffcn.input(INTEGRATOR_X0).size() \n
   *   nu: Number of controls: from ffcn.input(INTEGRATOR_P).size() - np \n
   *   np: Number of parameters: from option number_of_parameters\n
   *   nh: Number of point constraints: from cfcn.input(0).size()
   *
   * MultipleShooting is an CasADi::FX mapping from CasADi::OCPInput to CasADi::OCPOutput
   *
   *   \author Joel Andersson
   *   \date 2011
  */ 
class MultipleShooting : public OCPSolver{
  public:
    /// Default constructor
    MultipleShooting();
  
    /// Create a multiple shooting OCP solver
    /**
    * \param ffcn Discrete time dynamics, CasADi::FX mapping from CasADi::IntegratorInput to CasADi::IntegratorOutput
    * \param mfcn Mayer term, mapping endstate (nx x 1) to cost (1 x 1)
    * \param cfcn Path constraints, CasADi::FX mapping from CasADi::DAEInput to (nh x 1)
    * \param rfcn Initial value constraints
    */
    explicit MultipleShooting(const FX& ffcn, const FX& mfcn, const FX& cfcn=FX(), const FX& rfcn=FX());

    /// Access functions of the node
    MultipleShootingInternal* operator->();

    /// Const access functions of the node
    const MultipleShootingInternal* operator->() const;

    /// Get the NLP cost function
    FX getF() const;
    
    /// Get the NLP constraint function
    FX getG() const;
    
    /// Get the NLP Jacobian function
    FX getJ() const;

    /// Get the variables
    void getGuess(std::vector<double>& V_init) const;
    
    /// Get the variables
    void getVariableBounds(std::vector<double>& V_min, std::vector<double>& V_max) const;
    
    /// Get the constraints
    void getConstraintBounds(std::vector<double>& G_min, std::vector<double>& G_max) const;

    /// Set the optimal solution
    void setOptimalSolution( const std::vector<double> &V_opt );
    
    /// Set NLP solver
    void setNLPSolver(const NLPSolver& nlp_solver);
};
                        
                        
  } // namespace OptimalControl
} // namespace CasADi


#endif // MULTIPLE_SHOOTING_HPP
