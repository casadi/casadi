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


#ifndef CASADI_DAE_BUILDER_HPP
#define CASADI_DAE_BUILDER_HPP

#include "function.hpp"

namespace casadi {

// Forward declarations
class DaeBuilderInternal;

/** \brief A symbolic representation of a differential-algebraic equations model

    <H3>Variables:  </H3>
    \verbatim
    t:      independent variable (usually time)
    c:      constants
    p:      parameters
    d:      dependent parameters (time independent)
    u:      controls
    w:      dependent variables  (time dependent)
    x:      differential states
    z:      algebraic variables
    q:      quadrature states
    y:      outputs
    \endverbatim

    <H3>Equations:  </H3>
    \verbatim
    differential equations: \dot{x} ==  ode(...)
    algebraic equations:          0 ==  alg(...)
    quadrature equations:   \dot{q} == quad(...)
    dependent parameters:         d == ddef(d_prev,p)
    dependent variables:          w == wdef(w_prev,x,z,u,p,t)
    output equations:             y == ydef(...)
    initial equations:     init_lhs == init_rhs(...)
    events:      when when_cond < 0: when_lhs := when_rhs
    \endverbatim

    \date 2012-2021
    \author Joel Andersson

    \identifier{5c} */
class CASADI_EXPORT DaeBuilder
  : public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<DaeBuilder>) {
 public:

  /// Readable name of the class
  std::string type_name() const {return "DaeBuilder";}

  /// Default constructor
  DaeBuilder();

  /// Construct a DaeBuilder instance
  explicit DaeBuilder(const std::string& name, const std::string& path = "",
    const Dict& opts = Dict());

  /** \brief Name of instance

      \identifier{5d} */
  const std::string& name() const;

  /** @name Variables and equations */
  ///@{
  /** \brief Expression for independent variable (usually time)

      \identifier{2by} */
  const MX& time() const;

  /** \brief Independent variable (usually time)

      \identifier{2bz} */
  std::vector<std::string> t_new() const {return all("t");}

  /** \brief Differential states

      \identifier{5f} */
  std::vector<std::string> x() const {return all("x");}

  /// Outputs */
  std::vector<std::string> y() const;

  /** \brief Ordinary differential equations (ODE)

      \identifier{5g} */
  std::vector<MX> ode() const;

  /** \brief Algebraic variables

      \identifier{5h} */
  std::vector<std::string> z() const {return all("z");}

  /** \brief Algebraic equations

      \identifier{5i} */
  std::vector<MX> alg() const;

  /** \brief Quadrature states

      \identifier{5j} */
  std::vector<std::string> q() const {return all("q");}

  /** \brief Quadrature equations

      \identifier{5k} */
  std::vector<MX> quad() const;

  /** \brief Zero-crossing functions

      \identifier{2b0} */
  std::vector<MX> zero() const;

  /** \brief Definitions of output variables

      \identifier{5m} */
  std::vector<MX> ydef() const;

  /** \brief Set all output variables

      \identifier{2db} */
  void set_y(const std::vector<std::string>& name);

  /** \brief Get all rate variables

      \identifier{2dc} */
  std::vector<std::string> rate() const;

  /** \brief Set rate variables

      \identifier{2dd} */
  void set_rate(const std::vector<std::string>& name);

  /** \brief Free controls

      \identifier{5n} */
  std::vector<std::string> u() const {return all("u");}

  /** \brief Parameters

      \identifier{5o} */
  std::vector<std::string> p() const {return all("p");}

  /** \brief Named constants

      \identifier{5p} */
  std::vector<std::string> c() const {return all("c");}

  /** \brief Definitions of named constants

      \identifier{5q} */
  std::vector<MX> cdef() const;

  /** \brief Dependent parameters

      \identifier{5r} */
  std::vector<std::string> d() const {return all("d");}

  /** \brief Definitions of dependent parameters

    * Interdependencies are allowed but must be non-cyclic.

      \identifier{5s} */
  std::vector<MX> ddef() const;

  /** \brief Dependent variables

      \identifier{5t} */
  std::vector<std::string> w() const {return all("w");}

  /** \brief Dependent variables and corresponding definitions

   * Interdependencies are allowed but must be non-cyclic.

      \identifier{5u} */
  std::vector<MX> wdef() const;

  /** \brief Initial conditions, left-hand-side

      \identifier{2b1} */
  std::vector<MX> init_lhs() const;

  /** \brief Initial conditions, right-hand-side

      \identifier{2b2} */
  std::vector<MX> init_rhs() const;

  /** \brief Model structure: outputs

      \identifier{61} */
  std::vector<std::string> outputs() const;

  /** \brief Model structure: derivatives

      \identifier{62} */
  std::vector<std::string> derivatives() const;

  /** \brief Model structure: initial unknowns

      \identifier{63} */
  std::vector<std::string> initial_unknowns() const;

  /** @name Variables and equations */
  ///@{

  /** \brief Is there a time variable?

      \identifier{64} */
  bool has_t() const;

  /** \brief Is there a rate output?

      \identifier{2dg} */
  bool has_rate() const;

  /** \brief Differential states

      \identifier{65} */
  casadi_int nx() const;

  /** \brief Algebraic variables

      \identifier{66} */
  casadi_int nz() const;

  /** \brief Quadrature states

      \identifier{67} */
  casadi_int nq() const;

  /** \brief Zero-crossing functions

      \identifier{2cb} */
  casadi_int nzero() const;

  /** \brief Output variables

      \identifier{68} */
  casadi_int ny() const;

  /** \brief Free controls

      \identifier{69} */
  casadi_int nu() const;

  /** \brief Parameters

      \identifier{6a} */
  casadi_int np() const;

  /** \brief Named constants

      \identifier{6b} */
  casadi_int nc() const;

  /** \brief Dependent parameters

      \identifier{6c} */
  casadi_int nd() const;

  /** \brief Dependent variables

      \identifier{6d} */
  casadi_int nw() const;
  ///@}

  /** @name Symbolic modeling
   *  Formulate a dynamic system model
   */
  ///@{

  /// Add a new model variable
  MX add(const std::string& name,
    const std::string& causality,
    const std::string& variability,
    const Dict& opts=Dict());

  /// Add a new model variable, default variability
  MX add(const std::string& name,
    const std::string& causality,
    const Dict& opts=Dict());

  /// Add a new model variable, default variability and causality
  MX add(const std::string& name,
    const Dict& opts=Dict());

  /// Add a new model variable, symbolic expression already available
  void add(const std::string& name,
    const std::string& causality,
    const std::string& variability,
    const MX& expr,
    const Dict& opts=Dict());

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Renamed "time"
  const MX& t() const { return time();}

  /// [DEPRECATED] Replaced by add
  MX add_t(const std::string& name="t");

  /// [DEPRECATED] Replaced by add
  MX add_p(const std::string& name=std::string());

  /// [DEPRECATED] Replaced by add
  MX add_u(const std::string& name=std::string());

  /// [DEPRECATED] Replaced by add
  MX add_x(const std::string& name=std::string());

  /// [DEPRECATED] Replaced by add
  MX add_z(const std::string& name=std::string());

  /// [DEPRECATED] Replaced by add
  MX add_q(const std::string& name=std::string());

  /// [DEPRECATED] Replaced by add and eq
  MX add_c(const std::string& name, const MX& new_cdef);

  /// [DEPRECATED] Replaced by add and eq
  MX add_d(const std::string& name, const MX& new_ddef);

  /// [DEPRECATED] Replaced by add and eq
  MX add_w(const std::string& name, const MX& new_wdef);

  /// [DEPRECATED] Replaced by add and eq
  MX add_y(const std::string& name, const MX& new_ydef);

  /// [DEPRECATED] Replaced by eq
  void set_beq(const std::string& name, const MX& val);

  #endif  // WITH_DEPRECATED_FEATURES

  /// Add a simple equation
  void eq(const MX& lhs, const MX& rhs, const Dict& opts=Dict());

  /// Add when equations
  void when(const MX& cond, const std::vector<std::string>& eqs, const Dict& opts=Dict());

  /// Assignment inside a when-equation or if-else equation
  std::string assign(const std::string& name, const MX& val);

  /// Reinitialize a state inside when-equations
  std::string reinit(const std::string& name, const MX& val);

  /// Specify the initial equation for a variable
  void set_init(const std::string& name, const MX& init_rhs);

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Replaced by eq
  void set_ode(const std::string& name, const MX& ode_rhs) {
    eq(var(name), ode_rhs);
  }

  /// [DEPRECATED] Replaced by eq
  void set_alg(const std::string& name, const MX& alg_rhs) {
    (void)name;
    eq(0, alg_rhs);
  }

  /// [DEPRECATED] Replaced by set_init
  void add_init(const MX& lhs, const MX& rhs) {
    set_init(lhs.name(), rhs);
  }

  /// [DEPRECATED] Replaced by nzero()
  casadi_int ne() const {return nzero();}

  /// [DEPRECATED] Use all("zero") */
  std::vector<std::string> e() const {return all("zero");}

  #endif  // WITH_DEPRECATED_FEATURES

  /// Check if dimensions match
  void sanity_check() const;
  ///@}

  /// Reorder variables in a category
  void reorder(const std::string& cat, const std::vector<std::string>& v);

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Use set_variability, set_causality or set_category to change variable category
  void clear_all(const std::string& v);

  /// [DEPRECATED] Use set_variability, set_causality, set_category and/or reorder
  void set_all(const std::string& v, const std::vector<std::string>& name);

  /** @name [DEPRECATED] Register an existing variable */
  ///@{
  void register_t(const std::string& name);
  void register_p(const std::string& name);
  void register_u(const std::string& name);
  void register_x(const std::string& name);
  void register_z(const std::string& name);
  void register_q(const std::string& name);
  void register_c(const std::string& name);
  void register_d(const std::string& name);
  void register_w(const std::string& name);
  void register_y(const std::string& name);
  void register_e(const std::string& name);
  ///@}

  /// [DEPRECATED] Use eliminate("d")
  void eliminate_d();

  /// [DEPRECATED] Use eliminate("w")
  void eliminate_w();

  /// [DEPRECATED] Use eliminate("q")
  void eliminate_quad();

  /// [DEPRECATED] Use sort("d")
  void sort_d();

  /// [DEPRECATED] Use sort("w")
  void sort_w();

  /// [DEPRECATED] Use reorder("z", new_order)
  void sort_z(const std::vector<std::string>& z_order);

  #endif // WITH_DEPRECATED_FEATURES

  /** @name Manipulation
   *  Reformulate the dynamic optimization problem.
   */
  ///@{


  /// Eliminate all dependent parameters
  void eliminate(const std::string& cat);

  /// Sort dependent parameters
  void sort(const std::string& cat);

  /// Lift problem formulation by extracting shared subexpressions
  void lift(bool lift_shared = true, bool lift_calls = true);

  /// Prune unused controls
  void prune(bool prune_p = true, bool prune_u = true);

  /// Identify iteration variables and residual equations using naming convention
  void tear();
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

  /// Get all functions
  std::vector<Function> fun() const;

  /// Collect embedded functions from the expression graph
  void gather_fun(casadi_int max_depth = -1);
///@}

  /** @name Import and export
   */
  ///@{
  /// Import existing problem from FMI/XML
  void parse_fmi(const std::string& filename) {load_fmi_description(filename); }

  /// Does the FMU provide support for analytic derivatives
  bool provides_directional_derivatives() const;

  /// Does the FMU provide support for analytic derivatives (FMI 2 naming)
  bool provides_directional_derivative() const {return provides_directional_derivatives();}

  /// Import problem description from FMI or XML
  void load_fmi_description(const std::string& filename);

  /// Export instance into an FMU
  std::vector<std::string> export_fmu(const Dict& opts=Dict());

  /// Add a named linear combination of output expressions
  void add_lc(const std::string& name, const std::vector<std::string>& f_out);

  /// Construct a function object, legacy syntax
  Function create(const std::string& fname,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out, bool sx, bool lifted_calls = false) const;

  /** \brief  Construct a function object, names provided

    \param name    Name assigned to the resulting function object
    \param name_in   Names of all the inputs
    \param name_out  Names of all the outputs
    \param opts    Optional settings

      \identifier{6e} */
  Function create(const std::string& name,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const Dict& opts=Dict()) const;
  ///@}

  /** \brief Create a function with standard integrator DAE signature

    \param name    Name assigned to the resulting function object
    \param opts    Optional settings

      \identifier{2c0} */
  Function create(const std::string& fname, const Dict& opts=Dict()) const;

  /** \brief Create a function with standard integrator DAE signature, default naming

      \identifier{2c1} */
  Function create() const {return create(name() + "_dae");}

  /// Construct a function for evaluating dependent parameters
  Function dependent_fun(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out) const;

  /// Construct a function describing transition at a specific events
  Function transition(const std::string& fname, casadi_int index) const;

  /// Construct a function describing transition at any events
  Function transition(const std::string& fname) const;

  /// Construct an event transition function, default naming
  Function transition() const {return transition(name() + "_transition");}

  ///@{
  /// Get variable expression by name
  MX var(const std::string& name) const;
  MX operator()(const std::string& name) const {return var(name);}
  ///@}

  /// Get the time derivative of model variables
  std::vector<std::string> der(const std::vector<std::string>& name) const;

  ///@{
  /// Differentiate an expression with respect to time
  MX der(const MX& v) const;
  MX der(const MX& v);
  ///@}

  /// Get the pre-variables of model variables
  std::vector<std::string> pre(const std::vector<std::string>& name) const;

  /// Get the pre-expression given variable expression
  MX pre(const MX& v) const;

  /// Does a variable have a binding equation?
  bool has_beq(const std::string& name) const;

  /// Get the binding equation for a variable
  MX beq(const std::string& name) const;

  ///@{
  /// Get/set value reference
  casadi_int value_reference(const std::string& name) const;
  void set_value_reference(const std::string& name, casadi_int val);
  ///@}

  ///@{
  /// Get/set description
  std::string description(const std::string& name) const;
  void set_description(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the type
  std::string type(const std::string& name, casadi_int fmi_version = 3) const;
  void set_type(const std::string& name, const std::string& val);
  ///@}

  /// Get the causality
  std::string causality(const std::string& name) const;

  /** \brief Set the causality, if permitted

  The following changes are permitted:
    * For controls 'u' (variability 'continuous', causality 'input'), free parameters 'p'
    (variability 'tunable', causality 'parameter') and fixed parameters 'c' (variability
    'fixed', causality 'parameter'), causality can only be changed indirectly, by updating
    the variability
    * Add or remove an output 'y' by setting the causality to 'output' or 'local',
    respectively

    No other changes are permitted.

      \identifier{2c2} */
  void set_causality(const std::string& name, const std::string& val);

  /// Get the variability
  std::string variability(const std::string& name) const;

  /** \brief Set the variability, if permitted

  For controls 'u' (variability 'continuous', causality 'input'), free parameters 'p'
  (variability 'tunable', causality 'parameter') and fixed parameters 'c'
  (variability 'fixed', causality 'parameter'), update variability in
  order to change the category. Causality is updated accordingly.

  Other changes are not permitted

      \identifier{2c3} */
  void set_variability(const std::string& name, const std::string& val);

  /// Get the variable category
  std::string category(const std::string& name) const;

  /** \brief Set the variable category, if permitted

  The following changes are permitted:
    * Controls 'u' can be changed to/from tunable parameters 'p' or fixed parameters 'c'
    * Differential states that do not appear in the right-hand-sides can be changed between
    regular states 'x' and quadrature states 'q'

    Other changes are not permitted. Causality and variability is updated accordingly.

      \identifier{2c4} */
  void set_category(const std::string& name, const std::string& val);

  ///@{
  /// Get/set the initial property
  std::string initial(const std::string& name) const;
  void set_initial(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the unit
  std::string unit(const std::string& name) const;
  void set_unit(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /// Get/set the display unit
  std::string display_unit(const std::string& name) const;
  void set_display_unit(const std::string& name, const std::string& val);
  ///@}

  /// Get the number of elements of a variable
  casadi_int numel(const std::string& name) const;

  /// Get the dimensions of a variable
  std::vector<casadi_int> dimension(const std::string& name) const;

  /// Get the start time
  double start_time() const;

  /// Set the start time
  void set_start_time(double val);

  /// Get the stop time
  double stop_time() const;

  /// Set the stop time
  void set_stop_time(double val);

  /// Get the tolerance
  double tolerance() const;

  /// Set the tolerance
  void set_tolerance(double val);

  /// Get the step size
  double step_size() const;

  /// Set the step size
  void set_step_size(double val);

  // The following routines are not needed in MATLAB and would cause ambiguity
  // Note that a multirow strings can be interpreted as a vector of strings
#if !(defined(SWIG) && defined(SWIGMATLAB))
  /// Get the time derivative of model variables, single variable
  std::string der(const std::string& name) const;

  /// Get the pre-variables of model variables
  std::string pre(const std::string& name) const;

  /// Get an attribute, single variable
  double attribute(const std::string& a, const std::string& name) const;

  /// Set an attribute, single variable
  void set_attribute(const std::string& a, const std::string& name, double val);

  /// Get the lower bound, single variable
  double min(const std::string& name) const;

  /// Set the lower bound, single variable
  void set_min(const std::string& name, double val);

  /// Get the upper bound, single variable
  double max(const std::string& name) const;

  /// Set the upper bound, single variable
  void set_max(const std::string& name, double val);

  /// Get the nominal value, single variable
  double nominal(const std::string& name) const;

  /// Set the nominal value, single variable
  void set_nominal(const std::string& name, double val);

  /// Get the start attribute, single variable
  std::vector<double> start(const std::string& name) const;

  /// Set the start attribute, single variable
  void set_start(const std::string& name, double val);

  /// Set the start attribute, vector argument
  void set_start(const std::string& name, const std::vector<double>& val);

  // Clear all set values
  void reset();

  // Set the current value, single value
  void set(const std::string& name, double val);

  // Set the current value, single value (string)
  void set(const std::string& name, const std::string& val);

  /// Evaluate the values for a set of variables at the initial time, single value
  GenericType get(const std::string& name) const;

#endif  // !SWIGMATLAB

  /// Get an attribute
  std::vector<double> attribute(const std::string& a, const std::vector<std::string>& name) const;

  /// Set an attribute
  void set_attribute(const std::string& a, const std::vector<std::string>& name,
    const std::vector<double>& val);

  /// Get the lower bound
  std::vector<double> min(const std::vector<std::string>& name) const;

  /// Set the lower bound
  void set_min(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the upper bound
  std::vector<double> max(const std::vector<std::string>& name) const;

  /// Set the upper bound
  void set_max(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the nominal value
  std::vector<double> nominal(const std::vector<std::string>& name) const;

  /// Set the nominal value
  void set_nominal(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the start attribute
  std::vector<double> start(const std::vector<std::string>& name) const;

  /// Set the start attribute
  void set_start(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Set the current value
  void set(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Set the current value (string)
  void set(const std::vector<std::string>& name, const std::vector<std::string>& val);

  /// Evaluate the values for a set of variables at the initial time
  std::vector<GenericType> get(const std::vector<std::string>& name) const;

  /// Check if a particular variable exists
  bool has(const std::string& name) const;

  /// Get a list of all variables
  std::vector<std::string> all() const;

  /// Get a list of all variables of a particular category
  std::vector<std::string> all(const std::string& cat) const;

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Use add
  MX add_variable(const std::string& name, casadi_int n=1);

  /// [DEPRECATED] Use add
  MX add_variable(const std::string& name, const Sparsity& sp);

  /// Add a new variable from symbolic expressions
  void add_variable(const MX& new_v);

  /// [DEPRECATED] Use add
  size_t add_variable_new(const std::string& name, casadi_int n=1);

  /// [DEPRECATED] Use add
  size_t add_variable_new(const std::string& name, const Sparsity& sp);

  /// [DEPRECATED] Use add
  size_t add_variable_new(const MX& new_v);

  /// [DEPRECATED] Ranamed "has"
  bool has_variable(const std::string& name) const {return has(name);}

  /// Get a list of all variables
  std::vector<std::string> all_variables() const {return all();}
#endif // WITH_DEPRECATED_FEATURES

  /// Get the (cached) oracle, SX or MX
  Function oracle(bool sx = false, bool elim_w = false, bool lifted_calls = false) const;

  /** \brief Get Jacobian sparsity

      \identifier{6g} */
  Sparsity jac_sparsity(const std::vector<std::string>& onames,
    const std::vector<std::string>& inames) const;

#ifndef SWIG
#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Use add
  Variable& new_variable(const std::string& name, casadi_int numel = 1);

  ///@{
  /// [DEPRECATED] Access to internal class and corresponding indexing removed
  Variable& variable(const std::string& name);
  const Variable& variable(const std::string& name) const;
  Variable& variable(size_t ind);
  const Variable& variable(size_t ind) const;
  size_t find(const std::string& name) const;
  std::vector<size_t> find(const std::vector<std::string>& name) const;
  const std::string& name(size_t ind) const;
  std::vector<std::string> name(const std::vector<size_t>& ind) const;
  ///@}

  ///@{
  /// [DEPRECATED] Use string name, not internal index to access variables
  const MX& var(size_t ind) const;
  std::vector<MX> var(const std::vector<size_t>& ind) const;
  ///@}

#endif // WITH_DEPRECATED_FEATURES

  /// Access a member function or object
  const DaeBuilderInternal* operator->() const;

  /// Access a member function or object
  DaeBuilderInternal* operator->();

  /// Check if a particular cast is allowed
  static bool test_cast(const SharedObjectInternal* ptr);

#endif // SWIG
};

} // namespace casadi

#endif // CASADI_DAE_BUILDER_HPP
