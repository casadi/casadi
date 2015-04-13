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

  /** \brief A flat OCP representation
      <H3>Independent variables:  </H3>
      \verbatim
      t:      time
      \endverbatim

      <H3>Time-continuous variables:  </H3>
      \verbatim
      x:      states defined by ODE
      s:      implicitly defined states
      z:      algebraic variables
      u:      control signals
      q:      quadrature states
      i:      intermediate variables
      y:      outputs
      \endverbatim

      <H3>Time-constant variables:  </H3>
      \verbatim
      p:      free parameters
      \endverbatim

      <H3>Dynamic constraints (imposed everywhere):  </H3>
      \verbatim
      ODE                    \dot{x} ==  ode(t, x, s, z, u, p, i)
      DAE or implicit ODE:         0 ==  dae(t, x, s, z, u, p, i, sdot)
      algebraic equations:         0 ==  alg(t, x, s, z, u, p, i)
      quadrature equations:  \dot{q} == quad(t, x, s, z, u, p, i)
      intermediate equations:      i == idef(t, x, s, z, u, p, i)
      output equations:            y == ydef(t, x, s, z, u, p, i)
      \endverbatim

      <H3>Point constraints (imposed pointwise):  </H3>
      \verbatim
      Initial equations:           0 == init(t, x, s, z, u, p, i, sdot)
      \endverbatim

      <H3>Contributions to the objective function:  </H3>
      \verbatim
      Mayer terms:          \sum {mterm_k}
      Lagrange terms:       \sum {\integral{mterm}}
      \endverbatim

      \date 2012-2015
      \author Joel Andersson
  */
  class CASADI_EXPORT SymbolicOCP : public PrintableObject<SymbolicOCP> {
  public:

    /// Default constructor
    SymbolicOCP();

    /** @name Variables and equations
     *  Public data members
     */
    ///@{
    /** \brief Independent variable (usually time) */
    SX t;

    /** \brief Differential states defined by ordinary differential equations (ODE)
     */
    std::vector<SX> x, ode;

    /** \brief Differential-algebraic equation (DAE) with corresponding state vector,
     * state derivatives.
     */
    std::vector<SX> s, sdot, dae;

    /** \brief Algebraic equations and corresponding algebraic variables
     * \a alg and \a z have matching dimensions and
     * <tt>0 == alg(z, ...)</tt> implicitly defines \a z.
     */
    std::vector<SX> z, alg;

    /** \brief Quadrature states
     * Quadrature states are defined by ODEs whose state does not enter in the right-hand-side.
     */
    std::vector<SX> q, quad;

    /** \brief Intermediate variables and definitions definitions
     * Interdependencies are allowed but must be non-cyclic.
     */
    std::vector<SX> i, idef;

    /** \brief Output variables and corresponding definitions
     */
    std::vector<SX> y, ydef;

    /** \brief Free controls
     * The trajectories of the free controls are decision variables of the optimal control problem.
     * They are chosen by the optimization algorithm in order to minimize the cost functional.
     */
    std::vector<SX> u;

    /** \brief Free parameters
     * A free parameter is variables which is constant over time, but whose value is chosen by the
     * optimization algorithm in order to minimize the cost functional.
     */
    std::vector<SX> p;
    ///@}

    /** \brief Initial conditions
     * At <tt>t==0</tt>, <tt>0 == init(sdot, s, ...)</tt> holds in addition to
     * the ode and/or dae.
     */
    std::vector<SX> init;

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
    std::vector<SX> mterm;

    /// Lagrange terms in the objective (integral terms)
    std::vector<SX> lterm;
    ///@}

    /** @name Symbolic modeling
     *  Formulate an optimal control problem
     */
    ///@{
    /// Add a new differential state
    SX add_x(const std::string& name);

    /// Add a implicit state
    std::pair<SX, SX> add_s(const std::string& name);

    /// Add a new algebraic variable
    SX add_z(const std::string& name);

    /// Add a new parameter
    SX add_p(const std::string& name);

    /// Add a new control
    SX add_u(const std::string& name);

    /// Add an ordinary differential equation
    void add_ode(const SX& new_ode);

    /// Add a differential-algebraic equation
    void add_dae(const SX& new_dae);

    /// Add an algebraic equation
    void add_alg(const SX& new_alg);

    /// Check if dimensions match
    void sanityCheck() const;
    ///@}

    /** @name Manipulation
     *  Reformulate the dynamic optimization problem.
     */
    ///@{

    /// Identify and separate the algebraic variables and equations in the DAE
    void split_dae();

    /// Eliminate algebraic variables and equations transforming them into outputs
    void eliminate_alg();

    /// Transform the implicit DAE to a semi-explicit DAE
    void makeSemiExplicit();

    /// Transform the implicit DAE or semi-explicit DAE into an explicit ODE
    void makeExplicit();

    /// Sort intermediate variables
    void sort_i();

    /// Eliminate interdependencies amongst intermediate variables
    void split_i();

    /// Eliminate intermediate variables
    void eliminate_i();

    /// Eliminate Lagrange terms from the objective function and make them quadrature states
    void eliminate_lterm();

    /// Eliminate quadrature states and turn them into ODE states
    void eliminate_quad();

    /// Sort the DAE and implicitly defined states
    void sort_dae();

    /// Sort the algebraic equations and algebraic states
    void sort_alg();

    /// Scale the variables
    void scaleVariables();

    /// Scale the implicit equations
    void scaleEquations();
    ///@}

    /** @name Import and export
     */
    ///@{
    /// Import existing problem from FMI/XML
    void parseFMI(const std::string& filename);

    /// Generate a <tt>MUSCOD-II</tt> compatible DAT file
    void generateMuscodDatFile(const std::string& filename,
                               const Dictionary& mc2_ops=Dictionary()) const;

    /// Generate a file for numerical evaluation
    void generateCode(const std::string& filename,
                      const Dictionary& options=Dictionary());

    /// Generate a header file for generateCode
    static void generateHeader(const std::string& filename, const std::string& prefix="");

    /// Generate code for a particular function
    static void generateFunction(std::ostream &stream, const std::string& fname,
                                 const std::vector<SX>& f_in, const std::vector<SX>& f_out,
                                 CodeGenerator& gen,
                                 bool fwd=false, bool adj=false, bool foa=false);

    /// Corresponding header
    static void generateFunctionHeader(std::ostream &stream, const std::string& fname,
                                       bool fwd=false, bool adj=false, bool foa=false);
    ///@}

    /// Get variable expression by name
    SX operator()(const std::string& name) const;

    /// Get a derivative expression by name
    SX der(const std::string& name) const;

    /// Get a derivative expression by non-differentiated expression
    SX der(const SX& var) const;

    /// Get the nominal value by name
    double nominal(const std::string& name) const;

    /// Get the nominal value(s) by expression
    std::vector<double> nominal(const SX& var) const;

    /// Set the nominal value by name
    void setNominal(const std::string& name, double val);

    /// Set the nominal value(s) by expression
    void setNominal(const SX& var, const std::vector<double>& val);

    /// Get the lower bound by name
    double min(const std::string& name, bool normalized=false) const;

    /// Get the lower bound(s) by expression
    std::vector<double> min(const SX& var, bool normalized=false) const;

    /// Set the lower bound by name
    void setMin(const std::string& name, double val, bool normalized=false);

    /// Set the lower bound(s) by expression
    void setMin(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the upper bound by name
    double max(const std::string& name, bool normalized=false) const;

    /// Get the upper bound(s) by expression
    std::vector<double> max(const SX& var, bool normalized=false) const;

    /// Set the upper bound by name
    void setMax(const std::string& name, double val, bool normalized=false);

    /// Set the upper bound(s) by expression
    void setMax(const SX& var, const std::vector<double>& val, bool normalized=false);

    /// Get the initial guess by name
    double initialGuess(const std::string& name, bool normalized=false) const;

    /// Get the initial guess(es) by expression
    std::vector<double> initialGuess(const SX& var, bool normalized=false) const;

    /// Set the initial guess by name
    void setInitialGuess(const std::string& name, double val, bool normalized=false);

    /// Set the initial guess(es) by expression
    void setInitialGuess(const SX& var, const std::vector<double>& val, bool normalized=false);

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

    ///  Print representation
    void repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print description
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Add a variable
    void addVariable(const std::string& name, const Variable& var);

    /// Add a new variable: returns corresponding symbolic expression
    SX addVariable(const std::string& name);

    ///@{
    /// Access a variable by name
    Variable& variable(const std::string& name);
    const Variable& variable(const std::string& name) const;
    ///@}

#ifndef SWIG
    // Internal methods
  protected:

    /// Get the qualified name
    static std::string qualifiedName(const XmlNode& nn);

    /// Find of variable by name
    typedef std::map<std::string, Variable> VarMap;
    VarMap varmap_;

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
