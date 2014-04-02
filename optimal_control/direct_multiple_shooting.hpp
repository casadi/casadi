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

#ifndef DIRECT_MULTIPLE_SHOOTING_HPP
#define DIRECT_MULTIPLE_SHOOTING_HPP

#include "../symbolic/function/ocp_solver.hpp"
#include "../symbolic/function/nlp_solver.hpp"

namespace CasADi{
  class DirectMultipleShootingInternal;
    
    
  /** \brief Direct Multiple Shooting
   *
   *   ns: Number of shooting nodes: from option number_of_grid_points\n
   *   nx: Number of differential states: from ffcn.input(INTEGRATOR_X0).size() \n
   *   nc: Number of constants during intergation: ffcn.input(INTEGRATOR_P).size()
   *   nu: Number of controls: from nc - np \n
   *   np: Number of parameters: from option number_of_parameters\n
   *   nh: Number of point constraints: from cfcn.input(0).size()
   *
   *
   *   \author Joel Andersson
   *   \date 2011
  */ 
class DirectMultipleShooting : public OCPSolver{
  public:
    /// Default constructor
    DirectMultipleShooting();
  
    /** \brief Create a multiple shooting OCP solver
    * \param ffcn Continuous time dynamics, an CasADi::Function with the folowing mapping:
    * \copydoc scheme_DAEInput
    * \copydoc scheme_DAEOutput
    * Important notes:
    *  - In the above table, INTEGRATOR_P input is not really of shape (np x 1), but rather ( (np+nu) x 1 ).
    *  - The first np entries of the INTEGRATOR_P input are interpreted as parameters to be optimized but constant over the whole domain. The remainder are interpreted as controls. 
    *  - BEWARE: if the right hand side of ffcn is dependent on time, the results will be incorrect.
    *
    * \param mfcn Mayer term, CasADi::Function mapping to cost (1 x 1)
    * @copydoc scheme_MayerInput
    * \param cfcn Path constraints, CasADi::Function mapping to (nh x 1)
    * @copydoc scheme_DAEInput
    * \param rfcn Initial value constraints
    */
    explicit DirectMultipleShooting(const Function& ffcn, const Function& mfcn, const Function& cfcn=Function(), const Function& rfcn=Function());

    /// Access functions of the node
    DirectMultipleShootingInternal* operator->();

    /// Const access functions of the node
    const DirectMultipleShootingInternal* operator->() const;
    
    /// Get the variables
    void getGuess(std::vector<double>& V_init) const;
    
    /// Get the variables
    void getVariableBounds(std::vector<double>& V_min, std::vector<double>& V_max) const;
    
    /// Get the constraints
    void getConstraintBounds(std::vector<double>& G_min, std::vector<double>& G_max) const;

    /// Set the optimal solution
    void setOptimalSolution( const std::vector<double> &V_opt );
    
    // Access the underlying NLPSolver object
    NLPSolver getNLPSolver() const;

    // Prints out a human readable report about possible constraint violations, after solving 
    void reportConstraints(std::ostream &stream=std::cout);

    std::string getReportConstraints() { std::stringstream s; reportConstraints(s); return s.str(); }

    
    
    
};
                        
} // namespace CasADi


#endif // DIRECT_MULTIPLE_SHOOTING_HPP
