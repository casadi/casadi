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

  // Attributes
  enum Attribute {MIN, MAX, NOMINAL, START, N_ATTRIBUTE};

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
  MX min;
  MX max;
  MX nominal;
  // bool unbounded;
  MX start;
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

  // Get attribute by enum
  MX attribute(Attribute att) const;
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

  /// Constructor
  DaeBuilder(const std::string& name);

  /** @name Variables and equations */
  ///@{
  /** \brief Independent variable (usually time) */
  const MX t() const {return t_;}

  /** \brief Differential states */
  const std::vector<MX>& x() const {return x_;}

  /** \brief Ordinary differential equations (ODE) */
  const std::vector<MX>& ode() const {return ode_;}

  /** \brief Algebraic variables */
  const std::vector<MX>& z() const {return z_;}

  /** \brief Algebraic equations */
  const std::vector<MX>& alg() const {return alg_;}

  /** \brief Quadrature states */
  const std::vector<MX>& q() const {return q_;}

  /** \brief Quadrature equations */
  const std::vector<MX>& quad() const {return quad_;}

  /** \brief Output variables */
  const std::vector<MX>& y() const {return y_;}

  /** \brief Definitions of output variables */
  const std::vector<MX>& ydef() const {return ydef_;}

  /** \brief Free controls */
  const std::vector<MX>& u() const {return u_;}

  /** \brief Parameters */
  const std::vector<MX>& p() const {return p_;}

  /** \brief Named constants */
  const std::vector<MX>& c() const {return c_;}

  /** \brief Definitions of named constants */
  const std::vector<MX>& cdef() const {return cdef_;}

  /** \brief Dependent parameters */
  const std::vector<MX>& d() const {return d_;}

  /** \brief Definitions of dependent parameters
    * Interdependencies are allowed but must be non-cyclic.
    */
  const std::vector<MX>& ddef() const {return ddef_;}

  /** \brief Dependent variables */
  const std::vector<MX>& w() const {return w_;}

  /** \brief Dependent variables and corresponding definitions
   * Interdependencies are allowed but must be non-cyclic.
   */
  const std::vector<MX>& wdef() const {return wdef_;}

  /** \brief Auxiliary variables: Used e.g. to define functions */
  const std::vector<MX>& aux() const {return aux_;}

  /** \brief Initial conditions, left-hand-side */
  const std::vector<MX>& init_lhs() const {return init_lhs_;}

  /** \brief Initial conditions, right-hand-side */
  const std::vector<MX>& init_rhs() const {return init_rhs_;}
  ///@}

  /** @name Symbolic modeling
   *  Formulate a dynamic system model
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

  /// Add an initial equation
  void add_init(const MX& lhs, const MX& rhs);

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

  /// Eliminate all dependent variables
  void eliminate_w();

  /// Lift problem formulation by extracting shared subexpressions
  void lift(bool lift_shared = true, bool lift_calls = true);

  /// Eliminate quadrature states and turn them into ODE states
  void eliminate_quad();

  /// Sort dependent parameters
  void sort_d();

  /// Sort dependent variables
  void sort_w();

  /// Sort algebraic variables
  void sort_z(const std::vector<std::string>& z_order);

  /// Prune dependent parameters that cannot be evaluated
  void prune_d();

  /// Prune unused controls
  void prune(bool prune_p = true, bool prune_u = true);
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

  /// Construct a function object for evaluating attributes
  Function attribute_fun(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out) const;

  /// Construct a function for evaluating dependent parameters
  Function dependent_fun(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out) const;
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
  MX min(const std::string& name) const;
  void set_min(const std::string& name, const MX& val);
  ///@}

  ///@{
  /// Get/set the upper bound
  MX max(const std::string& name) const;
  void set_max(const std::string& name, const MX& val);
  ///@}

  ///@{
  /// Get/set the nominal value
  MX nominal(const std::string& name) const;
  void set_nominal(const std::string& name, const MX& val);
  ///@}

  ///@{
  /// Get/set the value at time 0
  MX start(const std::string& name) const;
  void set_start(const std::string& name, const MX& val);
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
  const Function& oracle(bool sx = false, bool elim_w = false, bool lifted_calls = false) const;

#ifndef SWIG
  // Internal methods
protected:

  /// Get the qualified name
  static std::string qualified_name(const XmlNode& nn);

  /// Name of instance
  std::string name_;

  /// All variables
  std::vector<Variable> variables_;

  /// Find of variable by name
  std::unordered_map<std::string, size_t> varind_;

  ///@{
  /// Ordered variables and equations
  MX t_;
  std::vector<MX> x_, ode_;
  std::vector<MX> z_, alg_;
  std::vector<MX> q_, quad_;
  std::vector<MX> y_, ydef_;
  std::vector<MX> u_;
  std::vector<MX> p_;
  std::vector<MX> c_, cdef_;
  std::vector<MX> d_, ddef_;
  std::vector<MX> w_, wdef_;
  std::vector<MX> aux_;
  std::vector<MX> init_lhs_, init_rhs_;
  ///@}

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

  // Sort dependent variables/parameters
  void sort_dependent(std::vector<MX>& v, std::vector<MX>& vdef);

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
template<> struct enum_traits<Variable::Attribute> {
  static const Variable::Attribute n_enum = Variable::N_ATTRIBUTE;
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
CASADI_EXPORT std::string to_string(Variable::Attribute v);
CASADI_EXPORT std::string to_string(DaeBuilder::DaeBuilderIn v);
CASADI_EXPORT std::string to_string(DaeBuilder::DaeBuilderOut v);
///@}

#endif  // SWIG

} // namespace casadi

#endif // CASADI_DAE_BUILDER_HPP
