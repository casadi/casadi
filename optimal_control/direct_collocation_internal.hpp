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

#ifndef DIRECT_COLLOCATION_INTERNAL_HPP
#define DIRECT_COLLOCATION_INTERNAL_HPP

#include "direct_collocation.hpp"
#include "../symbolic/fx/ocp_solver_internal.hpp"
#include "../symbolic/fx/parallelizer.hpp"
#include "../symbolic/fx/mx_function.hpp"
#include "../symbolic/fx/sx_function.hpp"
#include "integration/integration_tools.hpp"

namespace CasADi{
  
class DirectCollocationInternal : public OCPSolverInternal{
  friend class DirectCollocation;
  
  public:
    // Constructor
    DirectCollocationInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn);

    // clone
    virtual DirectCollocationInternal* clone() const{ return new DirectCollocationInternal(*this);}

    // Destructor
    virtual ~DirectCollocationInternal();
    
    // Initialize
    virtual void init();

    // Solve the OCP
    virtual void evaluate(int nfdir, int nadir);
   
    // Get the variables
    void getGuess(std::vector<double>& V_init) const;
    
    // Get the variables
    void getVariableBounds(std::vector<double>& V_min, std::vector<double>& V_max) const;
    
    // Get the constraints
    void getConstraintBounds(std::vector<double>& G_min, std::vector<double>& G_max) const;

    // Set the optimal solution
    void setOptimalSolution( const std::vector<double> &V_opt );
    
    // Prints out a human readable report about possible constraint violations - all constraints
    void reportConstraints(std::ostream &stream=std::cout);
    
    // NLP
    MXFunction nlp_;

    // NLP solver
    NLPSolver nlp_solver_;

    // Interpolation order
    int deg_;

};
                        
} // namespace CasADi


#endif // DIRECT_COLLOCATION_INTERNAL_HPP
