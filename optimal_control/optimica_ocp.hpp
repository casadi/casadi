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

#ifndef OPTIMICA_OCP_HPP
#define OPTIMICA_OCP_HPP

#include "casadi/printable_object.hpp"
#include "variable.hpp"

namespace CasADi{
  namespace OptimalControl{
    
/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCP : public PrintableObject{
  public:
    /// OCP
    OCP();

#ifndef SWIG
    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;
    
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;
#endif
    /// Sort variables according to type
    void sortType();

    /// Sort the variables and equations according to BLT
    void sortBLT();
    
    /// Try to make explicit by symbolically solving for xdot (experimental, only small systems)
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();
    
    /// Create a new, scaled OCP
    void scale();

    
    
    
    
    /// Access the variables in a class hierarchy -- public data member
    Variable variables;

    /// Time
    Variable t_;
    
    /// Differential states
    std::vector<Variable> x_;

    /// Algebraic states
    std::vector<Variable> z_;
    
    /// Controls
    std::vector<Variable> u_;
    
    /// Free parameters
    std::vector<Variable> p_;

    /// Constants
    std::vector<Variable> c_;

    /// Dependent
    std::vector<Variable> d_;
        
    /// Differential algebraic equations
    std::vector<SX> dae;
    
    /// Initial equations
    std::vector<SX> initeq;

    /// Path constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    /// Mayer objective terms
    std::vector<SX> mterm;

    /// Mayer time time points TODO: Remove this when WITH_TIMEDVARIABLE is default
    std::vector<double> mtp;
    
    /// Lagrange objective terms
    std::vector<SX> lterm;
    
    /// Initial time
    double t0;
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    
    /// Final time is free
    bool tf_free;
    
    /// Is scaled?
    bool is_scaled_;

};

#ifdef SWIG
%extend OCP{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif

  } // namespace OptimalControl
} // namespace CasADi

#endif // OPTIMICA_OCP_HPP


