/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "shared_object.hpp"
#include "casadi_enum.hpp"

namespace casadi {

// Forward declarations
struct XmlNode;

/// Variable type (FMI 2)
enum class TypeFmi2 {REAL, INTEGER, BOOLEAN, STRING, ENUM, NUMEL};

/// Variable type (FMI 3)
enum class Type {FLOAT32, FLOAT64, INT8, UINT8, INT16, UINT16, INT32, UINT32, INT64, UINT64,
  BOOLEAN, STRING, BINARY, ENUMERATION, CLOCK, NUMEL};

/// Causality: FMI 2.0 specification, section 2.2.7 or FMI 3.0 specification, section 2.4.7.4
enum class Causality {PARAMETER, CALCULATED_PARAMETER, INPUT, OUTPUT, LOCAL, INDEPENDENT, NUMEL};

/// Variability: FMI 2.0 specification, section 2.2.7 or FMI 3.0 specification, section 2.4.7.4
enum class Variability {CONSTANT, FIXED, TUNABLE, DISCRETE, CONTINUOUS, NUMEL};

// CasADi classification of model variables, cf. Table 17 in FMI specification, 3.0.2
//              PARAMETER  CALCULATED_PARAMETER  INPUT  OUTPUT  LOCAL  INDEPENDENT
// CONSTANT     -          -                     -      C       C      -
// FIXED        C          D                     -      -       D      -
// TUNABLE      P          D                     -      -       D      -
// DISCRETE     -          -                     U      Y       X      -
// CONTINUOUS   -          -                     U      Y/Q     X      T

// Input convention in codegen
enum class Category {T, C, P, D, W, U, X, Z, Q, Y, E, DER, RES, ASSIGN, REINIT, NUMEL};

/// Initial: FMI 2.0 specification, section 2.2.7 or FMI 3.0 specification, section 2.4.7.5
enum class Initial {EXACT, APPROX, CALCULATED, NA, NUMEL};

// Attributes
enum class Attribute {MIN, MAX, NOMINAL, START, VALUE, STRINGVALUE, NUMEL};

// Permitted dependenciesKind values
enum class DependenciesKind {DEPENDENT, CONSTANT, FIXED, TUNABLE, DISCRETE, NUMEL};

/** \brief Holds expressions and meta-data corresponding to a physical quantity evolving in time

    \date 2012-2021
    \author Joel Andersson

    \identifier{t} */
struct CASADI_EXPORT Variable {
  friend class DaeBuilderInternal;

 private:
  /// Constructor (only accessible via DaeBuilderInternal::new_variable)
  Variable(casadi_int index, const std::string& name,
    const std::vector<casadi_int>& dimension);

 public:
  /// @brief Location in variable vector
  casadi_int index;

  /// @brief  Name of the variable
  std::string name;

  /// Dimensions
  std::vector<casadi_int> dimension;

  /// Number of elements - product of all dimensions
  casadi_int numel;

  /** Attributes common to all types of variables, cf. Table 17 in FMI specification */
  ///@{
  unsigned int value_reference;
  std::string description;
  Type type;
  Causality causality;
  Variability variability;
  ///@}

  /// CasADi's classification of the variable
  Category category;

  /** Type specific attributes common to all types, cf. Table FMI 3.0 specification */
  ///@{
  // std::string declared_type;
  std::string unit;
  std::string display_unit;
  Initial initial;
  // std::string quantity;
  // bool relative_quantity;
  // bool unbounded;
  double min;
  double max;
  double nominal;
  std::vector<double> start;
  casadi_int der_of;  // 'derivative' in FMI specification
  casadi_int parent;
  ///@}

  // corresponding derivative variable
  casadi_int der;

  /// Numerical value (also for booleans, integers, enums)
  std::vector<double> value;

  /// String value (if string-valued)
  std::string stringvalue;

  /// Do other expressions depend on this variable
  bool dependency;

  /// Dependencies
  mutable std::vector<casadi_int> dependencies;

  /// Dependencies
  mutable std::vector<DependenciesKind> dependenciesKind;

  /// Variable expression (always a vector)
  MX v;

  /// Binding equation
  MX beq;

  /// Initial equation
  MX ieq;

  /// Total number of elements for a particular attribute
  casadi_int size(Attribute a) const;

  ///@{
  /// Get by attribute name
  void get_attribute(Attribute a, double* val) const;
  void get_attribute(Attribute a, std::vector<double>* val) const;
  void get_attribute(Attribute a, std::string* val) const;
  ///@}

  ///@{
  /// Set by attribute name
  void set_attribute(Attribute a, double val);
  void set_attribute(Attribute a, const std::vector<double>& val);
  void set_attribute(Attribute a, const std::string& val);
  ///@}

  // Default initial attribute, per specification
  static Initial default_initial(Causality causality, Variability variability);

  // Export as XML
  XmlNode export_xml(const DaeBuilderInternal& self) const;

  // Is the variable real?
  bool is_real() const {return type == Type::FLOAT32 || type == Type::FLOAT64;}

  // Does the variable need a start attribute?
  bool has_start() const;

  // Has the variable been set
  bool is_set() const {
    return !(type==Type::STRING ? stringvalue.empty() : std::isnan(value.front()));
  }

  // Does the variable have a binding equation
  bool has_beq() const {return !beq.is_empty();}

  // Does the variable need a derivative variable?
  bool needs_der() const;

  // Derivative of the expression
  MX dot(const DaeBuilderInternal& self) const;
};

/// \cond INTERNAL
/// Internal class for DaeBuilder, see comments on the public class.
class CASADI_EXPORT DaeBuilderInternal : public SharedObjectInternal {
  friend class DaeBuilder;
  friend class FmuInternal;
  friend class FmuFunction;

 public:

  /// Constructor
  explicit DaeBuilderInternal(const std::string& name, const std::string& path, const Dict& opts);

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

  /// Eliminate all dependent parameters
  void eliminate_d();

  /// Eliminate all dependent variables
  void eliminate_w();

  /// Lift problem formulation by extracting shared subexpressions
  void lift(bool lift_shared, bool lift_calls);

  /// Eliminate quadrature states and turn them into ODE states
  void eliminate_quad();

  /// Sort dependent parameters
  void sort_d();

  /// Sort dependent variables
  void sort_w();

  /// Sort algebraic variables
  void sort_z(const std::vector<std::string>& z_order);

  /// Classified variable indices (mutable)
  std::vector<size_t>& indices(Category cat);

  /// Classified variable indices (immutable)
  const std::vector<size_t>& indices(Category cat) const;

  /// Number of indices with a particular category
  size_t size(Category cat) const {return indices(cat).size();}

  /// Reorder variables in a category
  void reorder(Category cat, const std::vector<size_t>& v);

  /// Prune unused controls
  void prune(bool prune_p, bool prune_u);

  /// Identify free variables and residual equations
  void tearing_variables(std::vector<std::string>* res, std::vector<std::string>* iv,
    std::vector<std::string>* iv_on_hold) const;

  /// Identify free variables and residual equations
  void tear();
  ///@}

  /** @name Import and export
   */
  ///@{
  /// Import existing problem from FMI/XML
  void load_fmi_description(const std::string& filename);

  /// Get current date and time in the ISO 8601 format
  static std::string iso_8601_time();

  // Generate a random 32 digit hexadecimal number
  static std::string generate_guid();

  /// Export instance into an FMU (experimental)
  std::vector<std::string> export_fmu(const Dict& opts) const;

  /// Generate FMU wrapper file (fmi3Functions.c)
  std::string generate_wrapper(const std::string& guid, const CodeGenerator& gen) const;

  /// Generate buildDescription.xml
  std::string generate_build_description(const std::vector<std::string>& cfiles) const;

  /// Generate modelDescription.xml
  std::string generate_model_description(const std::string& guid) const;

  /// Generate FMU ModelVariables
  XmlNode generate_model_variables() const;

  /// Generate FMU ModelStructure
  XmlNode generate_model_structure() const;

  /// Update model variable dependencies
  void update_dependencies() const;

  ///@{
  /// Helper function: generate constants
  static std::string generate(const std::vector<size_t>& v);
  static std::string generate(const std::vector<double>& v);
  ///@}

  // Output convention in codegen
  enum DaeBuilderInternalOut {
    DAE_BUILDER_ODE,
    DAE_BUILDER_ALG,
    DAE_BUILDER_QUAD,
    DAE_BUILDER_ZERO,
    DAE_BUILDER_DDEF,
    DAE_BUILDER_WDEF,
    DAE_BUILDER_YDEF,
    DAE_BUILDER_NUM_OUT
  };

  // Get input expression, given enum
  std::vector<MX> input(Category ind) const;

  // Get output expression, given enum
  std::vector<MX> output(DaeBuilderInternalOut ind) const;

  // Get input expression, given enum
  std::vector<MX> input(const std::vector<Category>& ind) const;

  // Get output expression, given enum
  std::vector<MX> output(const std::vector<DaeBuilderInternalOut>& ind) const;

  /// Add a named linear combination of output expressions
  void add_lc(const std::string& name, const std::vector<std::string>& f_out);

  /// Construct a function object
  Function create(const std::string& fname,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts, bool sx, bool lifted_calls) const;

  /// Construct function from an FMU DLL
  Function fmu_fun(const std::string& fname,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts) const;

  /// Construct a function for evaluating dependent parameters
  Function dependent_fun(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out) const;

  /// Construct a function describing transition at a specific event
  Function transition(const std::string& fname, casadi_int index,
    bool dummy_index_input = false) const;

  /// Construct a function describing transition at all events
  Function transition(const std::string& fname) const;

  /// Function corresponding to all equations
  Function gather_eq() const;

  /// Get variable expression by name
  const MX& var(const std::string& name) const;

  /// Get a derivative expression by variable index
  MX der(size_t ind) const {return variable(ind).dot(*this);}

  /// Get a derivative expression by name
  MX der(const std::string& name) const {return der(find(name));}

  /// Get a derivative expression by non-differentiated expression
  MX der(const MX& var) const;

  /// Find a unique name, with a specific prefix
  std::string unique_name(const std::string& prefix) const;

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

  /// Create a new variable
  Variable& new_variable(const std::string& name,
    const std::vector<casadi_int>& dimension = {1});

  /// Check if a particular variable exists
  bool has(const std::string& name) const;

  /// Get a list of all variables
  std::vector<std::string> all() const;

  /// Get a list of all variables of a particular category
  std::vector<std::string> all(Category cat) const;

  /// Length of variables array
  size_t n_variables() const {return variables_.size();}

  /// Length of memory for all variables
  size_t n_mem() const;

  /// Start values for all variables
  std::vector<double> start_all() const;

  ///@{
  /// Access a variable by index
  Variable& variable(size_t ind) {return *variables_.at(ind);}
  const Variable& variable(size_t ind) const {return *variables_.at(ind);}
  ///@}

  ///@{
  /// Access a variable by Category and index
  Variable& variable(Category cat, size_t ind) {return variable(indices(cat).at(ind));}
  const Variable& variable(Category cat, size_t ind) const {
    return variable(indices(cat).at(ind));
  }
  ///@}

  ///@{
  /// Access a variable by name
  Variable& variable(const std::string& name) {return variable(find(name));}
  const Variable& variable(const std::string& name) const {return variable(find(name));}
  ///@}

  ///@{
  /// Access a variable by expression
  Variable& variable(const MX& v) {return variable(find(v));}
  const Variable& variable(const MX& v) const {return variable(find(v));}
  ///@}

  /// Get variable expression by index
  const MX& var(size_t ind) const;

  /// Get variable expression by category and index
  const MX& var(Category cat, size_t ind) const {return var(indices(cat).at(ind));}

  /// Get variable expressions by index
  std::vector<MX> var(const std::vector<size_t>& ind) const;

  /// Get variable expressions by category
  std::vector<MX> var(Category cat) const {return var(indices(cat));}

  /// Get index of variable, given name
  size_t find(const std::string& name) const;

  /// Get index of variable, given expression
  size_t find(const MX& v) const;

  /// Get indices of variable, given multiple names
  std::vector<size_t> find(const std::vector<std::string>& name) const;

  /// Get indices of variable, given multiple expressions
  std::vector<size_t> find(const std::vector<MX>& v) const;

  /** \brief Get variable name by index

      \identifier{2bv} */
  const std::string& name(size_t ind) const;

  /** \brief Get variable names by indices

      \identifier{2bw} */
  std::vector<std::string> name(const std::vector<size_t>& ind) const;

  /// Get the (cached) oracle, SX or MX
  const Function& oracle(bool sx = false, bool elim_w = false, bool lifted_calls = false) const;

  /// Get Jacobian sparsity
  Sparsity jac_sparsity(const std::vector<size_t>& oind, const std::vector<size_t>& iind) const;

  /// Get what is known of the Hessian sparsity
  Sparsity hess_sparsity(const std::vector<size_t>& oind, const std::vector<size_t>& iind) const;

  // Internal methods
protected:

  /// Get the qualified name
  static std::string qualified_name(const XmlNode& nn, Attribute* att = 0);

  // User-set options
  bool debug_;
  double fmutol_;
  bool ignore_time_;

  // FMI attributes
  std::string fmi_version_;
  std::string model_name_;
  std::string instantiation_token_;  // In FMI 2: guid
  std::string description_;
  std::string author_;
  std::string copyright_;
  std::string license_;
  std::string generation_tool_;
  std::string generation_date_and_time_;
  std::string variable_naming_convention_;
  casadi_int number_of_event_indicators_;

  // Default experiment
  double start_time_, stop_time_, tolerance_, step_size_;

  // Model Exchange
  std::string model_identifier_;
  bool provides_directional_derivatives_;
  bool provides_adjoint_derivatives_;
  bool can_be_instantiated_only_once_per_process_;
  std::vector<std::string> source_files_;

  /// Name of instance
  std::string name_;

  // Path to FMU, if any
  std::string path_;

  // Symbolic representation of the model equations?
  bool symbolic_;

  // Disallow quadrature states
  bool no_q_;

  // FMI major version
  casadi_int fmi_major_;

  /// All variables
  std::vector<Variable*> variables_;

  // Model structure
  std::vector<size_t> outputs_, derivatives_, initial_unknowns_;

  /// Find of variable by name
  std::unordered_map<std::string, size_t> varind_;

  /// Find of variable by value reference
  std::unordered_map<unsigned int, size_t> vrmap_;

  /// Ordered variables
  std::vector<std::vector<size_t>> indices_;

  // Initial equations
  std::vector<size_t> init_;

  // Event conditions and transition equations
  std::vector<MX> when_cond_, when_lhs_, when_rhs_;

 /** \brief Is there a time variable?

     \identifier{2bx} */
  bool has_t() const;

  // Time variable
  const MX& time() const;

  /** \brief Definitions of dependent constants

      \identifier{u} */
  std::vector<MX> cdef() const;

  /** \brief Definitions of dependent parameters

      \identifier{v} */
  std::vector<MX> ddef() const;

  /** \brief Definitions of dependent variables

      \identifier{w} */
  std::vector<MX> wdef() const;

  /** \brief Definitions of output variables

      \identifier{x} */
  std::vector<MX> ydef() const;

  /** \brief ODE right hand sides

      \identifier{y} */
  std::vector<MX> ode() const;

  /** \brief Algebraic right hand sides

      \identifier{z} */
  std::vector<MX> alg() const;

  /** \brief Quadrature right hand sides

      \identifier{10} */
  std::vector<MX> quad() const;

  /** \brief Zero crossing functions

      \identifier{2a6} */
  std::vector<MX> zero() const;

  /** \brief Initial conditions, left-hand-side

      \identifier{2a7} */
  std::vector<MX> init_lhs() const;

  /** \brief Initial conditions, right-hand-side

      \identifier{2a8} */
  std::vector<MX> init_rhs() const;

  /// Add a new variable
  Variable& add(const std::string& name, Causality causality, Variability variability,
    const Dict& opts);

  /// Add a new variable, default variability
  Variable& add(const std::string& name, Causality causality, const Dict& opts);

  /// Add a new variable, default variability and causality
  Variable& add(const std::string& name, const Dict& opts) {
    // Per FMI 3.0.2 specification, section 2.4.7.4: Default causality is LOCAL
    return add(name, Causality::LOCAL, opts);
  }

  /// Set or change the category for a variable
  void categorize(size_t ind, Category cat);

  /// Insert into list of variables, keeping it ordered
  void insert(std::vector<size_t>& v, size_t ind) const;

  /// Remove from list of variables
  void remove(std::vector<size_t>& v, size_t ind) const;

  /// Get causality
  Causality causality(size_t ind) const;

  /// Set causality
  void set_causality(size_t ind, Causality causality);

  /// Get variability
  Variability variability(size_t ind) const;

  /// Set variability
  void set_variability(size_t ind, Variability variability);

  /// Get category
  Category category(size_t ind) const;

  /// Set category
  void set_category(size_t ind, Category cat);

  /// Add a simple equation
  void eq(const MX& lhs, const MX& rhs, const Dict& opts);

  /// Add when equations
  void when(const MX& cond, const std::vector<std::string>& eqs, const Dict& opts);

  /// Assignment inside when-equations or if-else equations
  std::string assign(const std::string& name, const MX& val);

  /// Reinitialize a state inside when-equations
  std::string reinit(const std::string& name, const MX& val);

  /// Set a initial equation
  void set_init(const std::string& name, const MX& init_rhs);

  /// Linear combinations of output expressions
  Function::AuxOut lc_;

  /** \brief Functions

      \identifier{11} */
  std::vector<Function> fun_;

  /** \brief Function oracles (cached)

      \identifier{12} */
  mutable Function oracle_[2][2][2];

  /// Should the cache be cleared?
  mutable bool clear_cache_;

  /// Read an equation
  MX read_expr(const XmlNode& node);

  /// Read an identifier expression
  MX read_identifier(const XmlNode& node);

  /// Read a variable
  Variable& read_variable(const XmlNode& node, Attribute* att = 0);

  // Read DefaultExperiment
  void import_default_experiment(const XmlNode& n);

  // Read dependencies node
  std::vector<casadi_int> read_dependencies(const XmlNode& n);

  // Read dependenciesKind node
  std::vector<DependenciesKind> read_dependencies_kind(const XmlNode& n, size_t ndep);

  // Read ModelExchange
  void import_model_exchange(const XmlNode& n);

  // Read ModelVariables
  void import_model_variables(const XmlNode& modvars);

  // Read ModelStructure
  void import_model_structure(const XmlNode& n);

  // Read symbolic binding equations
  void import_binding_equations(const XmlNode& eqs);

  // Read symbolic dynamic equations
  void import_dynamic_equations(const XmlNode& eqs);

  // Read symbolic initial equations
  void import_initial_equations(const XmlNode& eqs);

  /// Problem structure has changed: Clear cache
  void clear_cache() const;

  /// Add a function from loaded expressions
  Function add_fun(const std::string& name,
                   const std::vector<std::string>& arg,
                   const std::vector<std::string>& res, const Dict& opts=Dict());

  /// Add an already existing function
  Function add_fun(const Function& f);

  /// Does a particular function already exist?
  bool has_fun(const std::string& name) const;

  /// Get function by name
  Function fun(const std::string& name) const;

  // Reset value attributes
  void reset();

  ///@{
  /// Get by attribute name
  double attribute(Attribute a, const std::string& name) const;
  std::vector<double> attribute(Attribute a, const std::vector<std::string>& name) const;
  ///@}

  ///@{
  /// Set by attribute name
  void set_attribute(Attribute a, const std::string& name, double val);
  void set_attribute(Attribute a, const std::vector<std::string>& name,
    const std::vector<double>& val);
  ///@}

  ///@{
  /// Get by attribute name (string-valued)
  std::string string_attribute(Attribute a, const std::string& name) const;
  std::vector<std::string> string_attribute(Attribute a,
    const std::vector<std::string>& name) const;
  ///@}

  ///@{
  /// Set by attribute name (string-valued)
  void set_string_attribute(Attribute a, const std::string& name, const std::string& val);
  void set_string_attribute(Attribute a, const std::vector<std::string>& name,
    const std::vector<std::string>& val);
  ///@}

  /// Total number of elements for a particular attribute
  casadi_int size(Attribute a, const std::vector<std::string>& name) const;

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
};

///@{
/// Number of entries in enums
template<> struct enum_traits<DaeBuilderInternal::DaeBuilderInternalOut> {
  static const size_t n_enum = DaeBuilderInternal::DAE_BUILDER_NUM_OUT;
};
///@}

///@{
/// Version mappings
CASADI_EXPORT Type from_fmi2(TypeFmi2 v);
CASADI_EXPORT TypeFmi2 to_fmi2(Type v);
///@}

///@{
/// Convert to string
CASADI_EXPORT std::string to_string(TypeFmi2 v);
CASADI_EXPORT std::string to_string(Type v);
CASADI_EXPORT std::string to_string(Causality v);
CASADI_EXPORT std::string to_string(Variability v);
CASADI_EXPORT std::string to_string(Initial v);
CASADI_EXPORT std::string to_string(Attribute v);
CASADI_EXPORT std::string to_string(DependenciesKind v);
CASADI_EXPORT std::string to_string(Category v);
CASADI_EXPORT std::string to_string(DaeBuilderInternal::DaeBuilderInternalOut v);
///@}

/// \endcond

} // namespace casadi

#endif // CASADI_DAE_BUILDER_INTERNAL_HPP
