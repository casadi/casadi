/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SYMBOLIC_OCP_HPP
#define CASADI_SYMBOLIC_OCP_HPP

#include "variable.hpp"

namespace casadi {

  // Forward declarations
  class XmlNode;

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
      explicit or implicit ODE: \dot {x} = ode(t, x, z, u, p_free, pi, pd)
      or                           0 = ode(t, x, z,\dot {x}, u, p_free, pi, pd)
      algebraic equations:            0 = alg(t, x, z, u, p_free, pi, pd)
      quadratures:              \dot {q} = quad(t, x, z, u, p_free, pi, pd)
      dependent equations:            y = dep(t, x, z, u, p_free, pi, pd)
      initial equations:              0 = initial(t, x, z, u, p_free, pi, pd)
      \endverbatim

      <H3>Objective function terms:  </H3>
      \verbatim
      Mayer terms:          \sum {mterm_k}
      Lagrange terms:       \sum {\integral{mterm}}
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

      When the optimal control problem is in a suitable form, it is possible to either
      generate functions for numeric/symbolic evaluation or exporting the OCP formulation
      into a new FMI conformant XML file. The latter functionality is not yet available.

      \date 2012
      \author Joel Andersson
  */
  class CASADI_CORE_EXPORT SymbolicOCP : public PrintableObject<SymbolicOCP> {
  public:

    /// Default constructor
    SymbolicOCP(bool ignore_timed_variables=true);

    /** @name Variables and equations
     *  Public data members
     */
    ///@{
    /** \brief Independent variable (usually time) */
    SX t;

    /** \brief Differential-algebraic equation (DAE) with corresponding state vector and initial
     * conditions
     * DAE in fully-implicit form and corresponding states and algebraic variables.
     * dae and s have matching dimensions and <tt>0 == dae(der(s), s, ...)</tt>
     * implicitly defines <tt>der(s)</tt>.
     * At <tt>t==0</tt>, <tt>0 == initial(der(s), s, ...)</tt> holds in addition to the dae.
     */
    SX s, dae, initial;

    /** \brief Differential states defined by ordinary differential equations (ODE)
     * The ODE can be retrieved by calling the method #ode with x as argument.
     */
    SX x;

    /** \brief Algebraic equations and corresponding algebraic variables
     * \a alg and \a z have matching dimensions and
     * <tt>0 == alg(z, ...)</tt> implicitly defines \a z.
     */
    SX z, alg;

    /** \brief Quadrature states
     * Quadrature states are defined by ODEs whose state does not enter in the right-hand-side.
     * The ODE can be retrieved by calling the method #ode with q as argument.
     */
    SX q;

    /** \brief Output variables
     * Interdependencies are allowed but must be non-cyclic.
     * The definitions can be retrieved by calling the method #beq with y as argument.
     */
    SX y;

    /** \brief Free controls
     * The trajectories of the free controls are decision variables of the optimal control problem.
     * They are chosen by the optimization algorithm in order to minimize the cost functional.
     */
    SX u;

    /** \brief Free parameters
     * A free parameter is variables which is constant over time, but whose value is chosen by the
     * optimization algorithm in order to minimize the cost functional.
     */
    SX p;

    /** \brief Independent parameters
     * An independent parameter is a parameter whose value is determined by an expression that
     * contains only literals.
     * An independent parameter is fixed after the DAE has been initialized.
     * The definitions can be retrieved by calling the method #beq with pi as argument.
     */
    SX pi;

    /** \brief Dependent parameters and corresponding definitions
     * A dependent parameter is a parameter whose value is determined by an expression which
     * contains references to other parameters.
     * A dependent parameter is fixed after the DAE has been initialized.
     * Interdependencies are allowed but must be non-cyclic.
     * The definitions can be retrieved by calling the method #beq with pd as argument.
    */
    SX pd;

    /** \brief Independent constant
     * An independent constant is a constant whose value is determined by an expression that
     * contains only literals.
     * The definitions can be retrieved by calling the method #beq with ci as argument.
     */
    SX ci;

    /** \brief Dependent constants and corresponding definitions
     * A dependent constant is a constant whose value is determined by an expression which
     * contains references to other constants.
     * Interdependencies are allowed but must be non-cyclic.
     * The definitions can be retrieved by calling the method #beq with cd as argument.
    */
    SX cd;
    ///@}

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

    ///@}

    /** @name Objective function terms
     *  Terms in the objective function.
     */
    ///@{

    /// Mayer terms in the objective (point terms)
    SX mterm;

    /// Lagrange terms in the objective (integral terms)
    SX lterm;
    ///@}

    /** \brief Path constraints of the optimal control problem
     */
    SX path;

    /** \brief Point constraints of the optimal control problem
     */
    SX point;

    /// Parse from XML to C++ format
    void parseFMI(const std::string& filename);

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);

    ///@{
    /// Access a variable by name
    Variable& variable(const std::string& name);
    const Variable& variable(const std::string& name) const;
    ///@}

    /** @name Manipulation
     *  Reformulate the dynamic optimization problem.
     */
    ///@{

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

    /// Sort the DAE and implicitly defined states
    void sortDAE();

    /// Sort the algebraic equations and algebraic states
    void sortALG();

    /// Generate a <tt>MUSCOD-II</tt> compatible DAT file
    void generateMuscodDatFile(const std::string& filename,
                               const Dictionary& mc2_ops=Dictionary()) const;

    ///@}

    /// Scale the variables
    void scaleVariables();

    /// Scale the implicit equations
    void scaleEquations();

    /// Get variable expression by name
    SX operator()(const std::string& name) const;

    /// Get a derivative expression by name
    SX der(const std::string& name) const;

    /// Get a derivative expression by non-differentiated expression
    SX der(const SX& var) const;

    /// Get a binding equation by name
    SX beq(const std::string& name) const;

    /// Get a binding equation by non-differentiated expression
    SX beq(const SX& var) const;

    /// Set a binding equation by name
    void setBeq(const std::string& name, const SX& val);

    /// Set an binding expression by non-differentiated expression
    void setBeq(const SX& var, const SX& val);

    /** \brief Get a derivative binding equation (i.e. ordinary differential equation, ODE)
     * by name.
     *
     * Returns variable expression if unknown.
     */
    SX ode(const std::string& name) const;

    /** \brief Get a derivative binding expression (i.e. ordinary differential equation, ODE)
     * by non-differentiated expression.
     *
     * Returns derivative expression if unknown.
     */
    SX ode(const SX& var) const;

    /// Set a derivative binding equation by name
    void setOde(const std::string& name, const SX& val);

    /// Set an derivative binding expression by non-differentiated expression
    void setOde(const SX& var, const SX& val);

    /// Get the nominal value by name
    double nominal(const std::string& name) const;

    /// Get the nominal value(s) by expression
    std::vector<double> nominal(const SX& var) const;

    /// Set the nominal value by name
    void setNominal(const std::string& name, double val);

    /// Set the nominal value(s) by expression
    void setNominal(const SX& var, const std::vector<double>& val);

    /// Get the lower bound by name
    SX min(const std::string& name) const;

    /// Get the lower bound(s) by expression
    SX min(const SX& var) const;

    /// Set the lower bound by name
    void setMin(const std::string& name, const SX& val);

    /// Set the lower bound(s) by expression
    void setMin(const SX& var, const SX& val);

    /// Get the upper bound by name
    SX max(const std::string& name) const;

    /// Get the upper bound(s) by expression
    SX max(const SX& var) const;

    /// Set the upper bound by name
    void setMax(const std::string& name, const SX& val);

    /// Set the upper bound(s) by expression
    void setMax(const SX& var, const SX& val);

    /// Get the initial guess by name
    SX initialGuess(const std::string& name) const;

    /// Get the initial guess(es) by expression
    SX initialGuess(const SX& var) const;

    /// Set the initial guess by name
    void setInitialGuess(const std::string& name, const SX& val);

    /// Set the initial guess(es) by expression
    void setInitialGuess(const SX& var, const SX& val);

    /// Get the (optionally normalized) value at time 0 by name
    double start(const std::string& name, bool normalized=false) const;

    /// Get the (optionally normalized) value(s) at time 0 by expression
    std::vector<double> start(const SX& var, bool normalized=false) const;

    /// Set the (optionally normalized) value at time 0 by name
    void setStart(const std::string& name, double val, bool normalized=false);

    /// Set the (optionally normalized) value(s) at time 0 by expression
    void setStart(const SX& var, const std::vector<double>& val, bool normalized=false);

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

    ///  Print representation
    void repr(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print description
    void print(std::ostream &stream=std::cout, bool trailing_newline=true) const;

#ifndef SWIG
    // Internal methods
  protected:

    /// Get the qualified name
    static std::string qualifiedName(const XmlNode& nn);

    /// Find of variable by name
    typedef std::map<std::string, Variable> VarMap;
    VarMap varmap_;

    /// Allow timed variables?
    bool ignore_timed_variables_;

    /// Read an equation
    SX readExpr(const XmlNode& odenode);

    /// Read a variable
    Variable& readVariable(const XmlNode& node);

    /// Get an attribute by expression
    typedef double (SymbolicOCP::*getAtt)(const std::string& name, bool normalized) const;
    std::vector<double> attribute(getAtt f, const SX& var, bool normalized) const;

    /// Get a symbolic attribute by expression
    typedef SX (SymbolicOCP::*getAttS)(const std::string& name) const;
    SX attribute(getAttS f, const SX& var) const;

    /// Set an attribute by expression
    typedef void (SymbolicOCP::*setAtt)(const std::string& name, double val, bool normalized);
    void setAttribute(setAtt f, const SX& var, const std::vector<double>& val, bool normalized);

    /// Set a symbolic attribute by expression
    typedef void (SymbolicOCP::*setAttS)(const std::string& name, const SX& val);
    void setAttribute(setAttS f, const SX& var, const SX& val);

#endif // SWIG

  };

} // namespace casadi

#endif // CASADI_SYMBOLIC_OCP_HPP
