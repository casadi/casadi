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


#ifndef CASADI_DAE_BUILDER_INTERNAL_HPP
#define CASADI_DAE_BUILDER_INTERNAL_HPP

#include <unordered_map>

#include "dae_builder.hpp"
#include "shared_object_internal.hpp"

namespace casadi {

// Forward declarations
class XmlNode;

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

/// \cond INTERNAL
/// Internal class for DaeBuilder, see comments on the public class.
class CASADI_EXPORT DaeBuilderInternal : public SharedObjectInternal {
  friend class DaeBuilder;

 public:

  /// Constructor
  explicit DaeBuilderInternal(const std::string& name);

  /// Destructor
  ~DaeBuilderInternal() override;

  /// Readable name of the internal class
  std::string class_name() const override {return "DaeBuilderInternal";}

  /// Check if dimensions match
  void sanity_check() const;

  /** @name Manipulation
   *  Reformulate the dynamic optimization problem.
   */
  ///@{

  /// Eliminate all dependent variables
  void eliminate_w();

  /// Lift problem formulation by extracting shared subexpressions
  void lift(bool lift_shared, bool lift_calls, bool inline_calls);

  /// Eliminate quadrature states and turn them into ODE states
  void eliminate_quad();

  /// Sort dependent parameters
  void sort_d();

  /// Sort dependent variables
  void sort_w();

  /// Sort algebraic variables
  void sort_z(const std::vector<std::string>& z_order);

  /// Clear input variable
  void clear_in(const std::string& v);

  /// Clear output variable
  void clear_out(const std::string& v);

  /// Prune dependent parameters that cannot be evaluated
  void prune_d();

  /// Prune unused controls
  void prune(bool prune_p, bool prune_u);
  ///@}

  /** @name Import and export
   */
  ///@{
  /// Import existing problem from FMI/XML
  void parse_fmi(const std::string& filename);

#ifndef SWIG
  // Input convension in codegen
  enum DaeBuilderInternalIn {
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
  enum DaeBuilderInternalOut {
    DAE_BUILDER_DDEF,
    DAE_BUILDER_WDEF,
    DAE_BUILDER_ODE,
    DAE_BUILDER_ALG,
    DAE_BUILDER_QUAD,
    DAE_BUILDER_YDEF,
    DAE_BUILDER_NUM_OUT
  };

  // Get input expression, given enum
  const std::vector<MX>& input(DaeBuilderInternalIn ind) const;

  // Get output expression, given enum
  std::vector<MX> output(DaeBuilderInternalOut ind) const;

  // Get input expression, given enum
  std::vector<MX> input(const std::vector<DaeBuilderInternalIn>& ind) const;

  // Get output expression, given enum
  std::vector<MX> output(const std::vector<DaeBuilderInternalOut>& ind) const;
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

  /// Function corresponding to all equations
  Function gather_eq() const;

  /// Get variable expression by name
  MX var(const std::string& name) const;

  /// Get variable expression by name
  MX operator()(const std::string& name) const {return var(name);}

  /// Get a derivative expression by name
  MX der(const std::string& name) const;

  /// Get a derivative expression by non-differentiated expression
  MX der(const MX& var) const;

  /// Readable name of the class
  std::string type_name() const {return "DaeBuilderInternal";}

  /// Print description
  void disp(std::ostream& stream, bool more) const override;

  /// Get string representation
  std::string get_str(bool more=false) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  /// Add a variable
  void add_variable(const std::string& name, const Variable& var);

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

  /// All input variables
  std::vector<std::vector<MX>> in_;

  /// Explicit time dependence?
  bool has_t() const { return !in_[DAE_BUILDER_T].empty();};

  /** \brief Shorthands for variable sets */
  ///@{
  const MX& t() const { return in_[DAE_BUILDER_T].at(0);}
  const std::vector<MX>& x() const { return in_[DAE_BUILDER_X];}
  const std::vector<MX>& z() const { return in_[DAE_BUILDER_Z];}
  const std::vector<MX>& q() const { return in_[DAE_BUILDER_Q];}
  const std::vector<MX>& y() const { return in_[DAE_BUILDER_Y];}
  const std::vector<MX>& u() const { return in_[DAE_BUILDER_U];}
  const std::vector<MX>& p() const { return in_[DAE_BUILDER_P];}
  const std::vector<MX>& c() const { return in_[DAE_BUILDER_C];}
  const std::vector<MX>& d() const { return in_[DAE_BUILDER_D];}
  const std::vector<MX>& w() const { return in_[DAE_BUILDER_W];}
  ///@}

  /** \brief Shorthands for variable sets: Non-const access */
  ///@{
  MX& t() { return in_[DAE_BUILDER_T].at(0);}
  std::vector<MX>& x() { return in_[DAE_BUILDER_X];}
  std::vector<MX>& z() { return in_[DAE_BUILDER_Z];}
  std::vector<MX>& q() { return in_[DAE_BUILDER_Q];}
  std::vector<MX>& y() { return in_[DAE_BUILDER_Y];}
  std::vector<MX>& u() { return in_[DAE_BUILDER_U];}
  std::vector<MX>& p() { return in_[DAE_BUILDER_P];}
  std::vector<MX>& c() { return in_[DAE_BUILDER_C];}
  std::vector<MX>& d() { return in_[DAE_BUILDER_D];}
  std::vector<MX>& w() { return in_[DAE_BUILDER_W];}
  ///@}

  ///@{
  /// Ordered variables and equations
  std::vector<MX> ode_;
  std::vector<MX> alg_;
  std::vector<MX> quad_;
  std::vector<MX> y_, ydef_;
  std::vector<MX> c_, cdef_;
  std::vector<MX> d_, ddef_;
  std::vector<MX> w_, wdef_;
  std::vector<MX> aux_;
  std::vector<MX> init_lhs_, init_rhs_;
  ///@}

  ///@{
  /// Add a new variable
  MX add_t(const std::string& name);
  MX add_p(const std::string& name, casadi_int n);
  MX add_u(const std::string& name, casadi_int n);
  MX add_x(const std::string& name, casadi_int n);
  MX add_z(const std::string& name, casadi_int n);
  MX add_q(const std::string& name, casadi_int n);
  MX add_c(const std::string& name, const MX& new_cdef);
  MX add_d(const std::string& name, const MX& new_ddef);
  MX add_w(const std::string& name, const MX& new_wdef);
  MX add_y(const std::string& name, const MX& new_ydef);
  ///@}

  ///@{
  /// Register an existing variable
  void register_t(const MX& new_t);
  void register_x(const MX& new_x);
  void register_z(const MX& new_z);
  void register_u(const MX& new_u);
  void register_p(const MX& new_p);
  void register_c(const MX& new_c, const MX& new_cdef);
  void register_d(const MX& new_d, const MX& new_ddef);
  void register_w(const MX& new_w, const MX& new_wdef);
  void register_y(const MX& new_y, const MX& new_ydef);
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
  static void sort_dependent(std::vector<MX>& v, std::vector<MX>& vdef);

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
template<> struct enum_traits<DaeBuilderInternal::DaeBuilderInternalIn> {
  static const DaeBuilderInternal::DaeBuilderInternalIn n_enum
    = DaeBuilderInternal::DAE_BUILDER_NUM_IN;
};
template<> struct enum_traits<DaeBuilderInternal::DaeBuilderInternalOut> {
  static const DaeBuilderInternal::DaeBuilderInternalOut n_enum
    = DaeBuilderInternal::DAE_BUILDER_NUM_OUT;
};
///@}

///@{
/// Convert to string
CASADI_EXPORT std::string to_string(Variable::Type v);
CASADI_EXPORT std::string to_string(Variable::Causality v);
CASADI_EXPORT std::string to_string(Variable::Variability v);
CASADI_EXPORT std::string to_string(Variable::Initial v);
CASADI_EXPORT std::string to_string(Variable::Attribute v);
CASADI_EXPORT std::string to_string(DaeBuilderInternal::DaeBuilderInternalIn v);
CASADI_EXPORT std::string to_string(DaeBuilderInternal::DaeBuilderInternalOut v);
///@}

#endif  // SWIG

} // namespace casadi

#endif // CASADI_DAE_BUILDER_INTERNAL_HPP
