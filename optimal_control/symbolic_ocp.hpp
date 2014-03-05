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

#ifndef SYMBOLIC_OCP_HPP
#define SYMBOLIC_OCP_HPP

#include "variable.hpp"

namespace CasADi{
  
  // Forward declarations
  class XMLNode;
    
  // All the types of variables in the SymbolicOCP class
  enum SymbolicOCPVariables{VAR_X, VAR_Z, VAR_Q, VAR_CI, VAR_CD, VAR_PI, VAR_PD, VAR_PF, VAR_Y, VAR_U, NUM_VAR};
  
/** \brief A flat OCP representation coupled to an XML file

 <H3>Variables:  </H3>
  \verbatim
   x:      differential states
   z:      algebraic states
   p :     independent parameters
   t :     time
   u :     control signals
   q :     quadrature states
   y :     dependent variables
  \endverbatim 

  <H3>Equations:  </H3>
  \verbatim
  explicit or implicit ODE: \dot{x} = ode(t,x,z,u,p_free,pi,pd)
     or                           0 = ode(t,x,z,\dot{x},u,p_free,pi,pd)
  algebraic equations:            0 = alg(t,x,z,u,p_free,pi,pd)
  quadratures:              \dot{q} = quad(t,x,z,u,p_free,pi,pd)
  dependent equations:            y = dep(t,x,z,u,p_free,pi,pd)
  initial equations:              0 = initial(t,x,z,u,p_free,pi,pd)
  \endverbatim 

  <H3>Objective function terms:  </H3>
  \verbatim
  Mayer terms:          \sum{mterm_k}
  Lagrange terms:       \sum{\integral{mterm}}
  \endverbatim

  Note that when parsed, all dynamic equations end up in the implicit category "dae". 
  At a later state, the DAE can be reformulated, for example in semi-explicit form, 
  possibly in addition to a set of quadrature states.
 
  The functions for reformulation is are provided as member functions to this class or as independent
  functions located in the header file "ocp_tools.hpp".

  <H3>Usage skeleton:</H3>
  
  1. Call default constructor 
  > SymbolicOCP ocp;
  
  2. Parse an FMI conformant XML file <BR>
  > ocp.parseFMI(xml_file_name)
  
  3. Modify/add variables, equations, optimization <BR>
  > ...
  
  When the optimal control problem is in a suitable form, it is possible to either generate functions
  for numeric/symbolic evaluation or exporting the OCP formulation into a new FMI conformant XML file.
  The latter functionality is not yet available.

  \date 2012
  \author Joel Andersson
*/
class SymbolicOCP : public PrintableObject{
  public:

    /// Default constructor
    SymbolicOCP();
    
    /** @name Variables categories
    *  Public data members
    */
    //@{
    /** \brief Time */
    SX t;
    
    /** \brief Differential states */
    std::vector<Variable> x;
    
    /** \brief Algebraic states */
    std::vector<Variable> z;
    
    /** \brief Quadrature states (length == quad().size()) */
    std::vector<Variable> q;

    /** \brief Independent constants */
    std::vector<Variable> ci;

    /** \brief Dependent constants */
    std::vector<Variable> cd;

    /** \brief Independent parameters 
     An independent parameter is a parameter whose value is determined by an expression that contains only literals: "parameter Real p1=2" or "parameter Boolean b(start=true)". In the latter case, the value of the parameter becomes true, and the Modelica compiler will generate a warning since there is no binding expression for the parameter. An independent parameter is fixed after the DAE has been initialized. */
    std::vector<Variable> pi;

    /** \brief Dependent parameters 
     A dependent parameter is a parameter whose value is determined by an expression which contains references to other parameters: "parameter Real p2=2*p1". A dependent parameter is fixed after the DAE has been initialized. */
    std::vector<Variable> pd;

    /** \brief Free parameters 
     A free parameter (which is Optimica specific without correspondance in Modelica) is a parameter that the optimization algorithm can change in order to minimize the cost function: "parameter Real x(free=true)". Note that these parameters in contrast to dependent/independent parameters may change after the DAE has been initialized. A free parameter should not have any binding expression since it would then no longer be free. The compiler will transform non-free parameters to free parameters if they depend on a free parameters. The "free" attribute thus propagage through the parameter binding equations. */
    std::vector<Variable> pf;
    
    /** \brief Dependent variables (length == dep().size()) */
    std::vector<Variable> y;
    
    /** \brief Control signals */
    std::vector<Variable> u;

    /** \brief Get all variables of a certain type */
    std::vector<Variable>& variableByType(SymbolicOCPVariables type);
    
    /** \brief Get all variables of a certain type */
    const std::vector<Variable>& variableByType(SymbolicOCPVariables type) const;
    //@}
    
    /** @name Equations
    *  Get all equations of a particular type 
    */
    //@{
      
    /// Explicit or implicit ODE
    SX ode;
    
    /// Algebraic equations
    SX alg;
    
    /// Quadrature equations
    SX quad;
    
    /// Dependent equations
    SX dep;
    
    /// Initial equations (remove?)
    SX initial;
    //@}
    
    /** @name Time points
    */
    //@{

    /// Interval start time
    double t0;
    
    /// Interval final time
    double tf;
    
    /// Interval start time is free
    bool t0_free;
    
    /// Interval final time is free
    bool tf_free;
    
    /// Interval start time initial guess
    double t0_guess;
    
    /// Interval final time initial guess
    double tf_guess;
    
    /// Time points
    std::vector<double> tp;
    
    //@}

    /** @name Objective function terms
    *  Terms in the objective function.
    */
    //@{
      
    /// Mayer terms in the objective (point terms)
    SX mterm;
    
    /// Lagrange terms in the objective (integral terms)
    SX lterm;
    //@}

    /** @name Path constraints of the optimal control problem
    */
    //@{

    /// Path constraint functions
    SX path;
    
    /// Path constraint functions bounds
    DMatrix path_min, path_max;
    //@}

    /** @name Point constraints of the optimal control problem
    */
    //@{

    /// Point constraint functions
    SX point;
    
    /// Path constraint functions bounds
    DMatrix point_min, point_max;
    //@}

    /// Parse from XML to C++ format
    void parseFMI(const std::string& filename, const Dictionary& options = Dictionary());

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);
    
    /// Access a variable by name
    Variable& variable(const std::string& name);

    /// Make a differential state algebraic by replacing its time derivative by 0
    void makeAlgebraic(const Variable& v);

    /// Make a differential state algebraic by replacing its time derivative by 0
    void makeAlgebraic(const std::string& name);
    
    /** @name Manipulation
    *  Reformulate the dynamic optimization problem.
    */
    //@{
    /// Eliminate interdependencies in the dependent equations
    void eliminateInterdependencies();
    
    /// Eliminate dependent equations, by default sparing the dependent variables with upper or lower bounds
    void eliminateDependent(bool eliminate_dependents_with_bounds=true);

    /// Eliminate Lagrange terms from the objective function and make them quadrature states
    void eliminateLagrangeTerms();
    
    /// Eliminate quadrature states and turn them into ODE states
    void eliminateQuadratureStates();
    
    /// Sort the ODE and differential states
    void sortODE();

    /// Sort the algebraic equations and algebraic states
    void sortALG();

    /// Sort the dependent parameters
    void sortDependentParameters();
    
    /// Transform the implicit ODE to an explicit ODE
    void makeExplicit();
    
    /// Eliminate algebraic states, transforming them into outputs
    void eliminateAlgebraic();
    
    /// Substitute the dependents from a set of expressions
    std::vector<SX> substituteDependents(const std::vector<SX>& x) const;
    
    /// Generate a MUSCOD-II compatible DAT file
    void generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops=Dictionary()) const;
    
    //@}
    
    /// Get the qualified name
    static std::string qualifiedName(const XMLNode& nn);
    
    /// Find of variable by name
    std::map<std::string,Variable> varmap_;

    /// Read an equation
    SXElement readExpr(const XMLNode& odenode, bool& has_der, bool elim_binding);

    /// Read a variable
    Variable& readVariable(const XMLNode& node);

    /// Scale the variables
    void scaleVariables();
    
    /// Scale the implicit equations
    void scaleEquations();
    
    #ifndef SWIG
    ///  Print representation
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print description 
    virtual void print(std::ostream &stream=std::cout) const;
    #endif // SWIG

};

} // namespace CasADi

#endif //SYMBOLIC_OCP_HPP
