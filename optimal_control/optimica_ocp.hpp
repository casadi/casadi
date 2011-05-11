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

    /// Tree structure for storing variables
    class VariableTree{
      public:
        /// Access a sub-collection by name
        VariableTree& subByName(const std::string& name, bool allocate=false);

        /// Access a sub-collection by index
        VariableTree& subByIndex(int ind, bool allocate=false);
        
        /// Get all variables
        void getAll(std::vector<Variable>& v) const;
    
        /// Get all names
        std::vector<std::string> getNames() const;
        
        /// Print node
        #ifndef SWIG
        void print(std::ostream &stream, int indent=0) const;
        #endif // SWIG

        /// Variable
        Variable var_;
        
        /// Children nodes
        std::vector<VariableTree> children_;
        
        /// Names of children
        std::map<std::string,int> name_part_;
    };

#ifdef SWIG
// Make sure that a copy constructor is created
%copyctor OCP;
#endif // SWIG
    
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

    /// Eliminate dependent equations
    void eliminateDependent();
    
    /// Sort the variables and equations according to BLT, with or without including the differentiated states in the dependency graph
    void sortBLT(bool with_x=false);
    
    /// Symbolically solve for xdot and z
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();
    
    /// Add a binding equation
    void addExplicitEquation(const Matrix<SX>& var, const Matrix<SX>& bind_eq);
    
    /// Scale the variables
    void scaleVariables();
    
    /// Scale the implicit equations
    void scaleEquations();
    
    /// Access the variables in a class hierarchy -- public data member
    VariableTree variables_;
    
    /// Create the implicit/explict ODE functions and quadrature state functions
    void createFunctions(bool create_dae=true, bool create_ode=true, bool create_quad=true);
    
    /// Time
    SX t_;
    
    /// Differential states
    std::vector<Variable> x_;

    /// Algebraic states
    std::vector<Variable> z_;
    
    /// Controls
    std::vector<Variable> u_;
        
    /// Free parameters
    std::vector<Variable> p_;
    
    /// Dependent variables
    std::vector<Variable> d_;

    /// Explicit equations
    std::vector<SX> explicit_var_, explicit_fcn_;
    
    /// Implicit equations
    std::vector<SX> implicit_var_, implicit_fcn_;
    
    /// Initial equations
    std::vector<SX> initial_eq_;

    /// Constraint function with upper and lower bounds
    std::vector<SX> path_fcn_;
    std::vector<double> path_min_, path_max_;
    FX pathfcn_;

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
    bool scaled_variables_, scaled_equations_;

    /// Has the dependents been eliminated
    bool eliminated_dependents_;
    
    /// BLT blocks
    std::vector<int> rowblock_;  // block k is rows r[k] to r[k+1]-1
    std::vector<int> colblock_;  // block k is cols s[k] to s[k+1]-1
    int nb_;

    /// BLT sorted?
    bool blt_sorted_;
    
    /// ODE right hand side function
    FX oderhs_;
    
    /// DAE residual function
    FX daeres_;
    
    /// Quadrature right hand side
    FX quadrhs_;

    /// Costs function
    FX costfcn_;
};

#ifdef SWIG
%extend OCP{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }

  %pythoncode %{
     def __deepcopy__(self,memo):
        return OCP(self)
  %}
}
#endif

  } // namespace OptimalControl
} // namespace CasADi

#endif // OPTIMICA_OCP_HPP


