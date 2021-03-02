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


#ifndef CASADI_DAE_BUILDER_HPP
#define CASADI_DAE_BUILDER_HPP

#include <unordered_map>

#include "function.hpp"

namespace casadi {

// Forward declarations
class XmlNode;
class DaeBuilder;

#ifndef SWIG
/** \brief Holds expressions and meta-data corresponding to a physical quantity evolving in time
    \date 2012-2021
    \author Joel Andersson
 */
struct CASADI_EXPORT Variable : public Printable<Variable> {
  /// Variable type
  enum Type {REAL, INTEGER, BOOLEAN, STRING, ENUM, N_TYPE};

  /// Causality: FMI 2.0 specification, section 2.2.7
  enum Causality {PARAMETER, CALCULATED_PARAMETER, INPUT, OUTPUT, LOCAL, INDEPENDENT, N_CAUSALITY};

  /// Variability: FMI 2.0 specification, section 2.2.7
  enum Variability {CONSTANT, FIXED, TUNABLE, DISCRETE, CONTINUOUS, N_VARIABILITY};

  /// Initial: FMI 2.0 specification, section 2.2.7
  enum Initial {EXACT, APPROX, CALCULATED, INITIAL_NA, N_INITIAL};

  /// Constructor
  Variable(const std::string& name = "");

  /** Attributes common to all types of variables, cf. FMI specification */
  ///@{
  std::string name;
  casadi_int value_reference;
  std::string description;
  Type type;
  Causality causality;
  Variability variability;
  Initial initial;
  ///@}

  /** Attributes specific to Real, cf. FMI specification */
  ///@{
  // std::string declared_type;
  // std::string quantity;
  std::string unit;
  std::string display_unit;
  // bool relative_quantity;
  double min;
  double max;
  double nominal;
  // bool unbounded;
  double start;
  casadi_int derivative;
  casadi_int antiderivative;
  // bool reinit;
  ///@}

  /// Do other expressions depend on this variable
  bool dependency;

  /// Variable expression
  MX v;

  /// Binding equation
  MX beq;

  /// Readable name of the class
  std::string type_name() const {return "Variable";}

  /// Print a description of the object
  void disp(std::ostream& stream, bool more=false) const;

  /// Get string representation
  std::string get_str(bool more=false) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  // Default initial attribute, per specification
  static Initial default_initial(Causality causality, Variability variability);
};

#endif  // SWIG

/** \brief An initial-value problem in differential-algebraic equations
    <H3>Independent variables:  </H3>
    \verbatim
    t:      time
    \endverbatim

    <H3>Variables:  </H3>
    \verbatim
    x:      differential states
    z:      algebraic variables
    u:      control signals
    q:      quadrature states
    p:      free parameters
    d:      dependent parameters
    w:      dependent variables
    y:      outputs
    \endverbatim

    <H3>Dynamic constraints (imposed everywhere):  </H3>
    \verbatim
    ODE                    \dot{x} ==  ode(t, x, z, u, p, w, d)
    algebraic equations:         0 ==  alg(t, x, z, u, p, w, d)
    quadrature equations:  \dot{q} == quad(t, x, z, u, p, w, d)
    dependent parameters:        d == ddef(t, x, z, u, p, w, d)
    dependent parameters:        w == wdef(t, x, z, u, p, w, d)
    output equations:            y == ydef(t, x, z, u, p, w, d)
    \endverbatim

    <H3>Point constraints (imposed pointwise):  </H3>
    \verbatim
    Initial equations:           0 == init(t, x, z, u, p, v, sdot)
    \endverbatim

    \date 2012-2021
    \author Joel Andersson
*/
class CASADI_EXPORT DaeBuilder
  : public SWIG_IF_ELSE(PrintableCommon, Printable<DaeBuilder>) {
public:

  /// Default constructor
  DaeBuilder();

  /** @name Variables and equations
   *  Public data members
   */
  ///@{
  /** \brief Independent variable (usually time) */
  MX t;

  /** \brief Differential states defined by ordinary differential equations (ODE)
   */
  std::vector<MX> x, ode, lam_ode;

  /** \brief Algebraic equations and corresponding algebraic variables
   * \a alg and \a z have matching dimensions and
   * <tt>0 == alg(z, ...)</tt> implicitly defines \a z.
   */
  std::vector<MX> z, alg, lam_alg;

  /** \brief Quadrature states
   * Quadrature states are defined by ODEs whose state does not enter in the right-hand-side.
   */
  std::vector<MX> q, quad, lam_quad;

  /** \brief Output variables and corresponding definitions
   */
  std::vector<MX> y, ydef, lam_ydef;

  /** \brief Free controls
   * The trajectories of the free controls are decision variables of the optimal control problem.
   * They are chosen by the optimization algorithm in order to minimize the cost functional.
   */
  std::vector<MX> u;

  /** \brief Parameters
   * A parameter is constant over time, but whose value is chosen by e.g. an
   * optimization algorithm.
   */
  std::vector<MX> p;

  /** \brief Named constants */
  std::vector<MX> c, cdef;

  /** \brief Dependent parameters and corresponding definitions
   * Interdependencies are allowed but must be non-cyclic.
   */
  std::vector<MX> d, ddef, lam_ddef;
  ///@}

  /** \brief Dependent variables and corresponding definitions
   * Interdependencies are allowed but must be non-cyclic.
   */
  std::vector<MX> w, wdef, lam_wdef;
  ///@}

  /** \brief Auxiliary variables: Used e.g. to define functions */
  std::vector<MX> aux;

  /** \brief Initial conditions
   * At <tt>t==0</tt>, <tt>0 == init(sdot, s, ...)</tt> holds in addition to
   * the ode and/or dae.
   */
  std::vector<MX> init;
  ///@}

  /** @name Symbolic modeling
   *  Formulate an optimal control problem
   */
  ///@{
  /// Add a new parameter
  MX add_p(const std::string& name=std::string(), casadi_int n=1);

  /// Add a new control
  MX add_u(const std::string& name=std::string(), casadi_int n=1);

  /// Add a new differential state
  MX add_x(const std::string& name=std::string(), casadi_int n=1);

  /// Add a new algebraic variable
  MX add_z(const std::string& name=std::string(), casadi_int n=1);

  /// Add a new quadrature state
  MX add_q(const std::string& name=std::string(), casadi_int n=1);

  /// Add a new dependent parameter
  MX add_d(const std::string& name, const MX& new_ddef);

  /// Add a new dependent variable
  MX add_w(const std::string& name, const MX& new_wdef);

  /// Add a new output
  MX add_y(const std::string& name, const MX& new_ydef);

  /// Add an ordinary differential equation
  void add_ode(const std::string& name, const MX& new_ode);

  /// Add an algebraic equation
  void add_alg(const std::string& name, const MX& new_alg);

  /// Add a quadrature equation
  void add_quad(const std::string& name, const MX& new_quad);

  /// Add an auxiliary variable
  MX add_aux(const std::string& name=std::string(), casadi_int n=1);

  /// Check if dimensions match
  void sanity_check() const;
  ///@}

  /** @name Register an existing variable */
  ///@{
  /// Register time variable
  void register_t(const MX& new_t);

  /// Register constant
  void register_c(const MX& new_c, const MX& new_cdef);

  /// Register dependent parameter
  void register_d(const MX& new_d, const MX& new_ddef);

  /// Register dependent variable
  void register_w(const MX& new_w, const MX& new_wdef);

  /// Register differential state
  void register_x(const MX& new_x);

  /// Register algebraic variable
  void register_z(const MX& new_z);

  /// Register input
  void register_u(const MX& new_u);

  /// Register free parameter
  void register_p(const MX& new_p);

  /// Register output variable
  void register_y(const MX& new_y, const MX& new_ydef);
  ///@}

  /** @name Manipulation
   *  Reformulate the dynamic optimization problem.
   */
  ///@{

  /// Lift problem formulation by extracting shared subexpressions
  void lift(bool lift_shared = true, bool lift_calls = true);

  /// Eliminate quadrature states and turn them into ODE states
  void eliminate_quad();
  ///@}

  /** @name Functions
   *  Add or load auxiliary functions
   */
  ///@{

  /// Add a function from loaded expressions
  Function add_fun(const std::string& name,
                   const std::vector<std::string>& arg,
                   const std::vector<std::string>& res, const Dict& opts=Dict());

  /// Add an already existing function
  Function add_fun(const Function& f);

  /// Add an external function
  Function add_fun(const std::string& name, const Importer& compiler,
                   const Dict& opts=Dict());

  /// Does a particular function already exist?
  bool has_fun(const std::string& name) const;

  /// Get function by name
  Function fun(const std::string& name) const;
///@}

  /** @name Import and export
   */
  ///@{
  /// Import existing problem from FMI/XML
  void parse_fmi(const std::string& filename);

#ifndef SWIG
  // Input convension in codegen
  enum DaeBuilderIn {
    DAE_BUILDER_T,
    DAE_BUILDER_C,
    DAE_BUILDER_P,
    DAE_BUILDER_D,
    DAE_BUILDER_W,
    DAE_BUILDER_U,
    DAE_BUILDER_X,
    DAE_BUILDER_Z,
    DAE_BUILDER_Q,
    DAE_BUILDER_Y,
    DAE_BUILDER_NUM_IN
  };

  // Output convension in codegen
  enum DaeBuilderOut {
    DAE_BUILDER_DDEF,
    DAE_BUILDER_WDEF,
    DAE_BUILDER_ODE,
    DAE_BUILDER_ALG,
    DAE_BUILDER_QUAD,
    DAE_BUILDER_YDEF,
    DAE_BUILDER_NUM_OUT
  };

  // Get input expression, given enum
  std::vector<MX> input(DaeBuilderIn ind) const;

  // Get output expression, given enum
  std::vector<MX> output(DaeBuilderOut ind) const;

  // Get input expression, given enum
  std::vector<MX> input(const std::vector<DaeBuilderIn>& ind) const;

  // Get output expression, given enum
  std::vector<MX> output(const std::vector<DaeBuilderOut>& ind) const;
#endif // SWIG

  /// Add a named linear combination of output expressions
  void add_lc(const std::string& name, const std::vector<std::string>& f_out);

  /// Construct a function object
  Function create(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out, bool sx = false, bool lifted_calls = false) const;
  ///@}

  /// Get variable expression by name
  MX var(const std::string& name) const;

  /// Get variable expression by name
  MX operator()(const std::string& name) const {return var(name);}

  /// Get a derivative expression by name
  MX der(const std::string& name) const;

  /// Get a derivative expression by non-differentiated expression
  MX der(const MX& var) const;

  ///@{
  /// Get/set description
  std::string description(const std::string& name) const;
  void set_description(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the type
  std::string type(const std::string& name) const;
  void set_type(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the causality
  std::string causality(const std::string& name) const;
  void set_causality(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the variability
  std::string variability(const std::string& name) const;
  void set_variability(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the initial property
  std::string initial(const std::string& name) const;
  void set_initial(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the unit
  std::string unit(const std::string& name) const;
  std::string unit(const MX& var) const;
  void set_unit(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the display unit
  std::string display_unit(const std::string& name) const;
  void set_display_unit(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the lower bound
  double min(const std::string& name) const;
  std::vector<double> min(const MX& var) const;
  void set_min(const std::string& name, double val);
  void set_min(const MX& var, const std::vector<double>& val);
  ///@}

  ///@{
  /// Get/set the upper bound
  double max(const std::string& name) const;
  std::vector<double> max(const MX& var) const;
  void set_max(const std::string& name, double val);
  void set_max(const MX& var, const std::vector<double>& val);
  ///@}

  ///@{
  /// Get/set the nominal value
  double nominal(const std::string& name) const;
  std::vector<double> nominal(const MX& var) const;
  void set_nominal(const std::string& name, double val);
  void set_nominal(const MX& var, const std::vector<double>& val);
  ///@}

  ///@{
  /// Get/set the value at time 0
  double start(const std::string& name) const;
  std::vector<double> start(const MX& var) const;
  void set_start(const std::string& name, double val);
  void set_start(const MX& var, const std::vector<double>& val);
  ///@}

  ///@{
  /// Get/set the binding equation
  const casadi::MX& binding_equation(const std::string& name) const;
  void set_binding_equation(const std::string& name, const MX& val);
  ///@}

  /// Readable name of the class
  std::string type_name() const {return "DaeBuilder";}

  ///  Print representation
  void disp(std::ostream& stream, bool more=false) const;

  /// Get string representation
  std::string get_str(bool more=false) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  /// Add a variable
  void add_variable(const std::string& name, const Variable& var);

  /// Add a new variable: returns corresponding symbolic expression
  MX add_variable(const std::string& name, casadi_int n=1);

  /// Add a new variable: returns corresponding symbolic expression
  MX add_variable(const std::string& name, const Sparsity& sp);

  /// Add a new variable from symbolic expressions
  void add_variable(const MX& new_v);

  /// Check if a particular variable exists
  bool has_variable(const std::string& name) const;

  ///@{
  /// Access a variable by name
  Variable& variable(const std::string& name);
  const Variable& variable(const std::string& name) const;
  ///@}

    /// Get the (cached) oracle, SX or MX
  const Function& oracle(bool sx = false, bool elim_v = false, bool lifted_calls = false) const;

#ifndef SWIG
  // Internal methods
protected:

  /// Get the qualified name
  static std::string qualified_name(const XmlNode& nn);

  /// All variables
  std::vector<Variable> variables_;

  /// Find of variable by name
  std::unordered_map<std::string, size_t> varind_;

  /// Linear combinations of output expressions
  Function::AuxOut lc_;

  /** \brief Functions */
  std::vector<Function> fun_;

  /** \brief Function oracles (cached) */
  mutable Function oracle_[2][2][2];

  /// Should the cache be cleared?
  mutable bool clear_cache_;

  /// Read an equation
  MX read_expr(const XmlNode& node);

  /// Read a variable
  Variable& read_variable(const XmlNode& node);

  /// Get an attribute by expression
  typedef double (DaeBuilder::*getAtt)(const std::string& name) const;
  std::vector<double> attribute(getAtt f, const MX& var) const;

  /// Get a symbolic attribute by expression
  typedef MX (DaeBuilder::*getAttS)(const std::string& name) const;
  MX attribute(getAttS f, const MX& var) const;

  /// Set an attribute by expression
  typedef void (DaeBuilder::*setAtt)(const std::string& name, double val);
  void set_attribute(setAtt f, const MX& var, const std::vector<double>& val);

  /// Set a symbolic attribute by expression
  typedef void (DaeBuilder::*setAttS)(const std::string& name, const MX& val);
  void set_attribute(setAttS f, const MX& var, const MX& val);

  /// Problem structure has changed: Clear cache
  void clear_cache() const;

  /// Helper class, represents inputs and outputs for a function call node
  struct CallIO {
    // Function instances
    Function f, adj1_f, J, H;
    // Index in v and vdef
    std::vector<size_t> v, vdef;
    // Nondifferentiated inputs
    std::vector<MX> arg;
    // Nondifferentiated inputs
    std::vector<MX> res;
    // Jacobian outputs
    std::vector<MX> jac_res;
    // Adjoint seeds
    std::vector<MX> adj1_arg;
    // Adjoint sensitivities
    std::vector<MX> adj1_res;
    // Hessian outputs
    std::vector<MX> hess_res;
    // Calculate Jacobian blocks
    void calc_jac();
    // Calculate gradient of Lagrangian
    void calc_grad();
    // Calculate Hessian of Lagrangian
    void calc_hess();
    // Access a specific Jacobian block
    const MX& jac(casadi_int oind, casadi_int iind) const;
    // Access a specific Hessian block
    const MX& hess(casadi_int iind1, casadi_int iind2) const;
  };

  /// Calculate contribution to jac_vdef_v from lifted calls
  MX jac_vdef_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
    const std::vector<casadi_int>& h_offsets) const;

  /// Calculate contribution to hess_?_v_v from lifted calls
  MX hess_v_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
    const std::vector<casadi_int>& h_offsets) const;

#endif // SWIG
};

#ifndef SWIG

/// Helper class: Specify number of entries in an enum
template<typename T>
struct enum_traits {};

// Helper function: Convert string to enum
template<typename T>
T to_enum(const std::string& s) {
  // Linear search over permitted values
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
    if (s == to_string(static_cast<T>(i))) return static_cast<T>(i);
  }
  // Informative error message
  std::stringstream ss;
  ss << "No such enum: '" << s << "'. Permitted values: ";
  for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
    // Separate strings
    if (i > 0) ss << ", ";
    // Print enum name
    ss << "'" << to_string(static_cast<T>(i)) << "'";
  }
  casadi_error(ss.str());
  return enum_traits<T>::n_enum;  // never reached
}

///@{
/// Number of entries in enums
template<> struct enum_traits<Variable::Type> {
  static const Variable::Type n_enum = Variable::N_TYPE;
};
template<> struct enum_traits<Variable::Causality> {
  static const Variable::Causality n_enum = Variable::N_CAUSALITY;
};
template<> struct enum_traits<Variable::Variability> {
  static const Variable::Variability n_enum = Variable::N_VARIABILITY;
};
template<> struct enum_traits<Variable::Initial> {
  static const Variable::Initial n_enum = Variable::N_INITIAL;
};
template<> struct enum_traits<DaeBuilder::DaeBuilderIn> {
  static const DaeBuilder::DaeBuilderIn n_enum = DaeBuilder::DAE_BUILDER_NUM_IN;
};
template<> struct enum_traits<DaeBuilder::DaeBuilderOut> {
  static const DaeBuilder::DaeBuilderOut n_enum = DaeBuilder::DAE_BUILDER_NUM_OUT;
};
///@}

///@{
/// Convert to string
CASADI_EXPORT std::string to_string(Variable::Type v);
CASADI_EXPORT std::string to_string(Variable::Causality v);
CASADI_EXPORT std::string to_string(Variable::Variability v);
CASADI_EXPORT std::string to_string(Variable::Initial v);
CASADI_EXPORT std::string to_string(DaeBuilder::DaeBuilderIn v);
CASADI_EXPORT std::string to_string(DaeBuilder::DaeBuilderOut v);
///@}

#endif  // SWIG

} // namespace casadi

#endif // CASADI_DAE_BUILDER_HPP
