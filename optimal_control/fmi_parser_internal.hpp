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

#ifndef FMI_PARSER_INTERNAL_HPP
#define FMI_PARSER_INTERNAL_HPP

#include "fmi_parser.hpp"
#include "optimica_ocp.hpp"
#include "xml_node.hpp"

/** \brief  Forward declarations */
class TiXmlElement;
class TiXmlNode;

namespace CasADi{
namespace OptimalControl{

class FMIParserInternal : public SharedObjectNode{
  public:
    /// Constructor
    explicit FMIParserInternal(const std::string& filename);
  
    /// destructor
    virtual ~FMIParserInternal(); 

    /// clone
    virtual FMIParserInternal* clone() const{ return new FMIParserInternal(*this);}

    /// Parse from XML to C++ format
    void parse();

    ///  Print representation
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print description 
    virtual void print(std::ostream &stream=std::cout) const;

    /// Add model variables */
    void addModelVariables();

    /// Add binding equations */
    void addBindingEquations();

    /// Add dynamic equations
    void addDynamicEquations();

    /// Read an equation
    SX readExpr(const XMLNode& odenode);

    /// Read a variable
    Variable& readVariable(const XMLNode& node);

    /// Add initial equations
    void addInitialEquations();

    /// Add optimization */
    void addOptimization();
    void addObjectiveFunction(const XMLNode& onode);
    void addIntegrandObjectiveFunction(const XMLNode& onode);
    void addConstraints(const XMLNode& onode);
    void addIntervalStartTime(const XMLNode& onode);
    void addIntervalFinalTime(const XMLNode& onode);

    // NOTE 1: Joel: The FMIParserInternal class will later have to be changed to work with the MX class instead of SX, 
    //               therefore I had to change the implementation so that it is more generic

    // NOTE 2: Joel: Will there really ever be so many functions that it will motivate a binary search of the functions rather than a simple linear search?

    /// Look-up table mapping XML names to SX unary functions
    std::map<std::string,SX (*)(const SX&)> unary_;

    /// Look-up table mapping XML names to SX binary functions
    std::map<std::string,SX (*)(const SX&,const SX&)> binary_;

    /// The optimal control problem representation -- keep synchronized with the XML representation!
    OCP ocp_;

    /// Parsed XML document
    XMLNode document_;

    /// Filename
    std::string filename_;
    
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
    void addExplicitEquation(const Matrix<SX>& var, const Matrix<SX>& bind_eq, bool to_front=false);
    
    /// Scale the variables
    void scaleVariables();
    
    /// Scale the implicit equations
    void scaleEquations();
    
    /// Create the implicit/explict ODE functions and quadrature state functions
    void createFunctions(bool create_dae=true, bool create_ode=true, bool create_quad=true);
    
    /// Make a differential state algebraic by replacing its time derivative by 0
    void makeAlgebraic(const Variable& v);
    
    /// Update the initial values for the dependent variables
    void findConsistentIC();

    /// Access the variables in a class hierarchy -- public data member
    VariableTree variables_;
    
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

    /// Lagrange objective terms
    std::vector<SX> lterm;
    
    /// Initial time
    double t0;
    double getStartTime() const{ return t0;}
    
    /// Initial time is free
    bool t0_free;
    
    /// Final time
    double tf;
    double getFinalTime() const{ return tf;}
    
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
    
    /// ODE right hand side function with lagrange term
    FX oderhs_with_lterm_;
    
    /// Output function
    FX output_fcn_;
    
    /// DAE residual function
    FX daeres_;
    
    /// Quadrature right hand side
    FX quadrhs_;

    /// Costs function
    FX costfcn_;
    
    /// Get the explicit expression for a variable
    SX getExplicit(const SX& v) const;
};

} // namespace OptimalControl
} // namespace CasADi

#endif //FMI_PARSER_INTERNAL_HPP
