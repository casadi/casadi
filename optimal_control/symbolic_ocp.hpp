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
    
    /** @name Variables and equations
     *  Public data members
     */
    //@{
    /** \brief Independent variable (usually time) */
    SX t;
    
    /** \brief Differential-algebraic equation (DAE) with corresponding state vector and initial conditions
     * DAE in fully-implicit form and corresponding states and algebraic variables. 
     * dae and s have matching dimensions and 0 == dae(der(s),s,...) implicitly defines der(s).
     * At t==0, 0 == initial(der(s),s,...) holds in addition to the dae.
     */
    SX s, dae, initial;

    /** \brief Ordinary differential equation (ODE) and corresponding state vector
     * ODE in explicit form and corresponding state vector.
     * ode and x have matching dimensions and der(x) == ode(x,...).
     */
    SX x, ode;

    /** \brief Algebraic equations and corresponding algebraic variables
     * alg and z have matching dimensions and 0 == alg(z,...) implicitly defines z.
     */
    SX z, alg;

    /** \brief Quadrature equations and corresponding quadrature states
     * Quadrature equation, e.g. an ODE whose state does not enter in the right-hand-side.
     * quad and q have matching dimensions and der(q) == quad(...)
     */
    SX q, quad;

    /** \brief Output variables and corresponding definitions
     * Interdependencies are allowed but must be non-cyclic.
     * y and def_y have matching dimensions and y == y_def(y,...)
     */
    SX y, y_def;

    /** \brief Free controls 
     * The trajectories of the free controls are decision variables of the optimal control problem. They are chosen by
     * the optimization algorithm in order to minimize the cost functional.
     */
    SX u;
    
    /** \brief Free parameters 
     * A free parameter is variables which is constant over time, but whose value is chosen by the optimization algorithm
     * in order to minimize the cost functional.
     */
    SX p;

    /** \brief Independent parameters
     * An independent parameter is a parameter whose value is determined by an expression that contains only literals.
     * An independent parameter is fixed after the DAE has been initialized.
     * Its value is located in the "value" attribute.
     */
    SX pi;

    /** \brief Dependent parameters and corresponding definitions
        A dependent parameter is a parameter whose value is determined by an expression which contains references to other parameters.
        A dependent parameter is fixed after the DAE has been initialized.        
    */
    SX pd, pd_def;

    /** \brief Independent constant
     * An independent constant is a constant whose value is determined by an expression that contains only literals.
     * Its value is located in the "value" attribute.
     */
    SX ci;

    /** \brief Dependent constants and correspinding definitions
     * A dependent constant is a constant whose value is determined by an expression which contains references to other constants.
    */
    SX cd, cd_def;
    //@}

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
    void parseFMI(const std::string& filename);

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);
    
    //@{
    /// Access a variable by name
    Variable& variable(const std::string& name);
    const Variable& variable(const std::string& name) const;
    //@}
    
    /** @name Manipulation
     *  Reformulate the dynamic optimization problem.
     */
    //@{

    /// Identify and separate the algebraic variables and equations in the DAE
    void separateAlgebraic();

    /// Eliminate algebraic variables, transforming them into outputs
    void eliminateAlgebraic();

    /// Transform the implicit DAE to a semi-explicit DAE
    void makeSemiExplicit();

    /// Transform the implicit DAE or semi-explicit DAE into an explicit ODE
    void makeExplicit();

    /// Eliminate independent parameters
    void eliminateIndependentParameters();
    
    /// Sort the dependent parameters
    void sortDependentParameters();
    
    /// Eliminate interdependencies amongst the dependent parameters
    void eliminateDependentParameterInterdependencies();

    /// Eliminate dependent parameters
    void eliminateDependentParameters();

    /// Sort the outputs
    void sortOutputs();

    /// Eliminate interdependencies amongst the outputs
    void eliminateOutputInterdependencies();

    /// Eliminate outputs
    void eliminateOutputs();

    /// Eliminate Lagrange terms from the objective function and make them quadrature states
    void eliminateLagrangeTerms();
    
    /// Eliminate quadrature states and turn them into ODE states
    void eliminateQuadratureStates();
    
    /// Sort the DAE and implictly defined states
    void sortDAE();

    /// Sort the algebraic equations and algebraic states
    void sortALG();
        
    /// Generate a MUSCOD-II compatible DAT file
    void generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops=Dictionary()) const;
    
    //@}

    /// Scale the variables
    void scaleVariables();
    
    /// Scale the implicit equations
    void scaleEquations();

    /// Find an expression by name
    SX operator()(const std::string& name) const;

    /// Find an derivative expression by name
    SX der(const std::string& name) const;

    /// Find an derivative expression by non-differentiated expression
    SX der(const SX& var) const;

    /// Find a binding expression by name
    SX binding(const std::string& name) const;

    /// Find an binding expression by non-differentiated expression
    SX binding(const SX& var) const;
    
    /// Get the nominal value by name
    double nominal(const std::string& name) const;
    
    /// Get the nominal value(s) by expression
    std::vector<double> nominal(const SX& var) const;

    /// Set the nominal value by name
    void setNominal(const std::string& name, double val);

    /// Set the nominal value(s) by expression
    void setNominal(const SX& var, const std::vector<double>& val);

    /// Get the (optionally normalized) current value by name
    double value(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) current value(s) by expression
    std::vector<double> value(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) current value by name
    void setValue(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) current value(s) by expression
    void setValue(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) lower bound by name
    double min(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) lower bound(s) by expression
    std::vector<double> min(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) lower bound by name
    void setMin(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) lower bound(s) by expression
    void setMin(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) upper bound by name
    double max(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) upper bound(s) by expression
    std::vector<double> max(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) upper bound by name
    void setMax(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) upper bound(s) by expression
    void setMax(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) value at time 0 by name
    double start(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) value(s) at time 0 by expression
    std::vector<double> start(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) value at time 0 by name
    void setStart(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) value(s) at time 0 by expression
    void setStart(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) initial guess by name
    double initialGuess(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) initial guess(es) by expression
    std::vector<double> initialGuess(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) initial guess by name
    void setInitialGuess(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) initial guess(es) by expression
    void setInitialGuess(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the (optionally normalized) derivative value at time 0 by name
    double derivativeStart(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) derivative value(s) at time 0 by expression
    std::vector<double> derivativeStart(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) derivative value at time 0 by name
    void setDerivativeStart(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) derivative value(s) at time 0 by expression
    void setDerivativeStart(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the unit for a component
    std::string unit(const std::string& name) const;

    /// Get the unit given a vector of symbolic variables (all units must be identical)
    std::string unit(const SX& var) const;

    /// Set the unit for a component
    void setUnit(const std::string& name, const std::string& val);

    /// Timed variable (never allocate)
    SX atTime(const std::string& name, double t, bool allocate=false) const;

    /// Timed variable (allocate if necessary)
    SX atTime(const std::string& name, double t, bool allocate=false);


#ifndef SWIG
    ///  Print representation
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print description 
    virtual void print(std::ostream &stream=std::cout) const;

    // Internal methods
  protected:

    /// Get the qualified name
    static std::string qualifiedName(const XMLNode& nn);
    
    /// Find of variable by name
    std::map<std::string,Variable> varmap_;

    /// Read an equation
    SX readExpr(const XMLNode& odenode);

    /// Read a variable
    Variable& readVariable(const XMLNode& node);

    /// Get an attribute by expression
    typedef double (SymbolicOCP::*getAtt)(const std::string& name, bool normalized) const;
    std::vector<double> attribute(getAtt f, const SX& var, bool normalized) const;

    /// Set an attribute by expression
    typedef void (SymbolicOCP::*setAtt)(const std::string& name, double val, bool normalized);  
    void setAttribute(setAtt f, const SX& var, const std::vector<double>& val, bool normalized);

#endif // SWIG

  };

} // namespace CasADi

#endif //SYMBOLIC_OCP_HPP
