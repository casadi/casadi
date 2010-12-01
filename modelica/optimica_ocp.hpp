/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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
  namespace Modelica{

/** Symbolic, object oriented representation of an optimal control problem (OCP) */
class OCP : public PrintableObject{
  public:    
    /** \brief OCP */
    OCP();

    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Sort variables
    void sortVariables();

    /// Try to make explicit by symbolically solving for xdot (experimental, only small systems)
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();

    /// Access the variables in a class hierarchy -- public data member
    Variable variables;
    
  public:

    /// Time variable(s) TODO: Should not be a vector
    std::vector<SX> t;
    
    /// Differential states appearing implicitly
    std::vector<SX> x;

    /// Time derivative of the differential states appearing implicitly
    std::vector<SX> xdot;
 
    /// Differential states
    std::vector<SX> xd;
 
    /// Algebraic states
    std::vector<SX> xa;
    
    /// Controls
    std::vector<SX> u;
    
    /// Parameters
    std::vector<SX> p;

    /// Dependent variables and constants
    std::vector<Variable> d;
    
    // EQUATIONS

    /// Fully implicit equations
    std::vector<SX> dyneq;
    
    /// Explicit differential equations
    std::vector<SX> diffeq;

    /// Algebraic equations
    std::vector<SX> algeq;
    
    /// Initial equations
    std::vector<SX> initeq;
    
    /// Definition of dependent variables
    std::vector<SX> depdef;

    // OBJECTIVE
    /// Mayer terms
    std::vector<SX> mterm;
    
    /// Mayer time time point
    std::vector<double> mtp;
    
    /// Constraint function with upper and lower bounds
    std::vector<SX> cfcn, cfcn_lb, cfcn_ub;

    /// Lagrange terms (symbolic/numeric)
    // std::vector<LagrangeTerm> lterm;
    
    /// Least squares terms (symbolic/numeric)
    // std::vector<LagrangeTerm> lsqterm;
    
    /// Initial time
    double t0;
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    
    /// Final time is free
    bool tf_free;

};

  } // namespace Modelica
} // namespace CasADi

#endif // OPTIMICA_OCP_HPP


