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

    ///  Print representation
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print description 
    virtual void print(std::ostream &stream=std::cout) const;

    /// Parse from XML to C++ format
    void parse();

    /// Add model variables
    void addModelVariables();

    /// Add binding equations
    void addBindingEquations();

    /// Add dynamic equations
    void addDynamicEquations();

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

    /// Read an equation
    SX readExpr(const XMLNode& odenode);

    /// Read a variable
    Variable& readVariable(const XMLNode& node);

    /// Sort variables according to type
    void sortType();
    
    /// Eliminate interdependencies in the dependent equations
    void eliminateInterdependencies();

    /// Eliminate dependent equations
    void eliminateDependent();

    /// Eliminate Lagrange terms from the objective function and make them quadrature states
    void eliminateLagrangeTerms();

    /// Eliminate quadrature states and turn them into ODE states
    void eliminateQuadratureStates();

    /// Sort the DAE equations and variables
    void sortDAE();
    
    /// Transform the fully implicit DAE to a explicit or semi-explicit form
    void makeExplicit();

    /// All states, differential and algebraic
    std::vector<Variable> x_all() const;

    /// Get the DAE input arguments
    std::vector<SXMatrix> daeArg() const;

    /// Substitute the dependents from a set of expressions
    std::vector<SXMatrix> substituteDependents(const std::vector<SXMatrix>& x) const;
    
    /// Get the ODE/DAE right hand side function
    FX daeFcn() const;

    /// Generate a MUSCOD-II compatible DAT file
    void generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops) const;

    /// Scale the variables
    void scaleVariables();
    
    /// Scale the implicit equations
    void scaleEquations();
    
    /// Make a differential state algebraic by replacing its time derivative by 0
    void makeAlgebraic(const Variable& v);
    
    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);
    
    /// Access a variable by name
    Variable& variable(const std::string& name);

    /// Get the qualified name
    static std::string qualifiedName(const XMLNode& nn);
    
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
    std::vector<Variable> xq_;

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
    
};

} // namespace OptimalControl
} // namespace CasADi

#endif //FLAT_OCP_INTERNAL_HPP
