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
   * Available options, from OCPSolver:
   * "number_of_parameters", OT_INTEGER,  0
   * "number_of_grid_points", OT_INTEGER,  20
   * "final_time",OT_REAL, 1.0
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
    * \param ffcn Discrete time dynamics, should have I/O of an integrator
    * \param mfcn Mayer term, mappping endstate -> cost
    * \param cfcn Path constraints
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
