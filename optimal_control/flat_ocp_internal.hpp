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

#ifndef FLAT_OCP_INTERNAL_HPP
#define FLAT_OCP_INTERNAL_HPP

#include "flat_ocp.hpp"
#include "xml_node.hpp"

/** \brief  Forward declarations */
class TiXmlElement;
class TiXmlNode;

namespace CasADi{
namespace OptimalControl{

class FlatOCPInternal : public OptionsFunctionalityNode{
  public:
    /// Constructor
    explicit FlatOCPInternal(const std::string& filename);
  
    /// destructor
    virtual ~FlatOCPInternal(); 

    /// Initialize
    virtual void init();
    
    /// clone
    virtual FlatOCPInternal* clone() const{ return new FlatOCPInternal(*this);}

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

    //@{
    /// Add optimization
    void addOptimization();
    void addObjectiveFunction(const XMLNode& onode);
    void addIntegrandObjectiveFunction(const XMLNode& onode);
    void addConstraints(const XMLNode& onode);
    void addIntervalStartTime(const XMLNode& onode);
    void addIntervalFinalTime(const XMLNode& onode);
    //@}

    /// Parsed XML document
    XMLNode document_;

    /// Time
    SX t_;

    /// States (implicitly defined)
    std::vector<Variable> x_;

    /// Differential states
    std::vector<Variable> xd_;

    /// Algebraic states
    std::vector<Variable> xa_;

    /// Quadrature states
    std::vector<Variable> q_;

    /// Controls
    std::vector<Variable> u_;
        
    /// Free parameters
    std::vector<Variable> p_;

    /// Dependent variables
    std::vector<Variable> y_;

    /// fully implicit DAE
    std::vector<SX> dae_;
    
    /// explicit OD
    std::vector<SX> ode_;
    
    /// quadratures
    std::vector<SX> quad_;
    
    /// algebraic equations
    std::vector<SX> alg_;
    
    /// dependent equations
    std::vector<SX> dep_;
    
    /// initial equations
    std::vector<SX> initial_;
    
    /// Filename
    std::string filename_;
    
    /// Sort variables according to type
    void sortType();
    
    /// Eliminate interdependencies in the dependent equations
    void eliminateInterdependencies();

    /// Eliminate dependent equations
    void eliminateDependent(bool eliminate_dependents_with_bounds=false);
    
    /// Sort the variables and equations according to BLT, with or without including the differentiated states in the dependency graph
    void sortBLT(bool with_x=false);
    
    /// Symbolically solve for xdot and z
    void makeExplicit();

    /// Replace all state derivatives by algebraic variables with the same name
    void makeSemiExplicit();
    
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

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);
    
    /// Access a variable by name
    Variable& variable(const std::string& name);

    /// Get the qualified name
    static std::string qualifiedName(const XMLNode& nn);
    
    /// Find of variable by name
    std::map<std::string,Variable> varmap_;
    
    /// Constraint function
    FX pathfcn_;

    /// Mayer objective terms
    std::vector<SX> mterm_;

    /// Lagrange objective terms
    std::vector<SX> lterm_;

    /// Path constraints
    std::vector<SX> path_;

    /// Path constraints upper and lower bounds
    std::vector<double> path_min_, path_max_;
    
    /// Initial time
    double t0_;
    
    /// Initial time is free
    bool t0_free_;
    
    /// Final time
    double tf_;
    
    /// Final time is free
    bool tf_free_;

    /// Verbose parsing
    bool verbose_;
    
    /// Have the variables been scaled
    bool scaled_variables_;

    /// Have the equations been scaled
    bool scaled_equations_;
    
};

} // namespace OptimalControl
} // namespace CasADi

#endif //FLAT_OCP_INTERNAL_HPP
