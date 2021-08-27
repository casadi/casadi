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
*/
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

  /** \brief Name of instance */
  const std::string& name() const;

  /** @name Variables and equations */
  ///@{
  /** \brief Independent variable (usually time) */
  const MX& t() const;

  /** \brief Differential states */
  std::vector<std::string> x() const;

  /** \brief Ordinary differential equations (ODE) */
  std::vector<MX> ode() const;

  /** \brief Algebraic variables */
  std::vector<std::string> z() const;

  /** \brief Algebraic equations */
  std::vector<MX> alg() const;

  /** \brief Quadrature states */
  std::vector<std::string> q() const;

  /** \brief Quadrature equations */
  std::vector<MX> quad() const;

  /** \brief Output variables */
  std::vector<std::string> y() const;

  /** \brief Definitions of output variables */
  std::vector<MX> ydef() const;

  /** \brief Free controls */
  std::vector<std::string> u() const;

  /** \brief Parameters */
  std::vector<std::string> p() const;

  /** \brief Named constants */
  std::vector<std::string> c() const;

  /** \brief Definitions of named constants */
  std::vector<MX> cdef() const;

  /** \brief Dependent parameters */
  std::vector<std::string> d() const;

  /** \brief Definitions of dependent parameters
    * Interdependencies are allowed but must be non-cyclic.
    */
  std::vector<MX> ddef() const;

  /** \brief Dependent variables */
  std::vector<std::string> w() const;

  /** \brief Dependent variables and corresponding definitions
   * Interdependencies are allowed but must be non-cyclic.
   */
  std::vector<MX> wdef() const;

  /** \brief Auxiliary variables: Used e.g. to define functions */
  const std::vector<MX>& aux() const;

  /** \brief Initial conditions, left-hand-side */
  const std::vector<MX>& init_lhs() const;

  /** \brief Initial conditions, right-hand-side */
  const std::vector<MX>& init_rhs() const;

  /** \brief When statement: triggering condition */
  const std::vector<MX>& when_cond() const;

  /** \brief When statement: left-hand-side */
  const std::vector<MX>& when_lhs() const;

  /** \brief When statement: right-hand-side */
  const std::vector<MX>& when_rhs() const;
  ///@}

  /** @name Variables and equations */
  ///@{

  /** \brief Is there a time variable? */
  bool has_t() const;

  /** \brief Differential states */
  casadi_int nx() const;

  /** \brief Algebraic variables */
  casadi_int nz() const;

  /** \brief Quadrature states */
  casadi_int nq() const;

  /** \brief Output variables */
  casadi_int ny() const;

  /** \brief Free controls */
  casadi_int nu() const;

  /** \brief Parameters */
  casadi_int np() const;

  /** \brief Named constants */
  casadi_int nc() const;

  /** \brief Dependent parameters */
  casadi_int nd() const;

  /** \brief Dependent variables */
  casadi_int nw() const;
  ///@}

  /** @name Symbolic modeling
   *  Formulate a dynamic system model
   */
  ///@{
  /// Add an independent variable (time)
  MX add_t(const std::string& name="t");

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

  /// Add a new constant
  MX add_c(const std::string& name, const MX& new_cdef);

  /// Add a new dependent parameter
  MX add_d(const std::string& name, const MX& new_ddef);

  /// Add a new dependent variable
  MX add_w(const std::string& name, const MX& new_wdef);

  /// Add a new output
  MX add_y(const std::string& name, const MX& new_ydef);

  /// Add an ordinary differential equation
  MX add_ode(const std::string& name, const MX& new_ode);

  /// Add an algebraic equation
  MX add_alg(const std::string& name, const MX& new_alg);

  /// Add a quadrature equation
  MX add_quad(const std::string& name, const MX& new_quad);

  /// Add an auxiliary variable
  MX add_aux(const std::string& name=std::string(), casadi_int n=1);

  /// Add an initial equation
  void add_init(const MX& lhs, const MX& rhs);

  /// Add a when statement
  void add_when(const MX& cond, const MX& lhs, const MX& rhs);

  /// Check if dimensions match
  void sanity_check() const;
  ///@}

  /** @name Register an existing variable */
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
  ///@}

  /** @name Manipulation
   *  Reformulate the dynamic optimization problem.
   */
  ///@{

  /// Clear input variable
  void clear_in(const std::string& v);

  /// Clear output variable
  void clear_out(const std::string& v);

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

  /// Import problem description from FMI or XML
  void load_fmi_description(const std::string& filename);

  /// Import FMI functions from DLL
  void load_fmi_functions(const std::string& path);

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

  /** \brief  Load a function from an FMU DLL, indices provided
    \param name    Name assigned to the resulting function object
    \param comp_in   Names of all the input components
    \param comp_out  Names of all the output components
    \param name_in   Names of all the inputs
    \param name_out  Names of all the outputs
    \param opts    Optional settings
  */
  Function fmu_fun(const std::string& name,
      const std::vector<std::vector<std::string>>& comp_in,
      const std::vector<std::vector<std::string>>& comp_out,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts=Dict()) const;
  ///@}

  /** \brief  Load a function from an FMU DLL, names provided
    \param name    Name assigned to the resulting function object
    \param name_in   Names of all the inputs
    \param name_out  Names of all the outputs
    \param opts    Optional settings
  */
  Function fmu_fun(const std::string& name,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts=Dict()) const;
  ///@}

  /** \brief  Load a function from an FMU DLL, standard IO conforming with simulator
    \param name    Name assigned to the resulting function object
    \param opts    Optional settings
  */
  Function fmu_fun(const std::string& name, const Dict& opts=Dict()) const;

  /// Free all FMU instances
  void reset_fmu() const;

  ///@{
  /// Get variable expression by name
  MX var(const std::string& name) const;
  std::vector<MX> var(const std::vector<std::string>& name) const;
  ///@}

  ///@{
  /// Get binding equation by name
  MX beq(const std::string& name, bool derivative = false) const;
  std::vector<MX> beq(const std::vector<std::string>& name, bool derivative = false) const;
  ///@}

  /// Get variable expression by name
  ///@{
  MX operator()(const std::string& name) const {return var(name);}
  std::vector<MX> operator()(const std::vector<std::string>& name) const {return var(name);}
  ///@}

  /// Get a derivative expression by name
  MX der(const std::string& name) const;

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

  /// Add a new variable: returns corresponding symbolic expression
  MX add_variable(const std::string& name, casadi_int n=1);

  /// Add a new variable: returns corresponding symbolic expression
  MX add_variable(const std::string& name, const Sparsity& sp);

  /// Add a new variable from symbolic expressions
  void add_variable(const MX& new_v);

  /// Add a new variable: returns corresponding symbolic expression
  size_t add_variable_new(const std::string& name, casadi_int n=1);

  /// Add a new variable: returns corresponding symbolic expression
  size_t add_variable_new(const std::string& name, const Sparsity& sp);

  /// Add a new variable from symbolic expressions
  size_t add_variable_new(const MX& new_v);

  /// Check if a particular variable exists
  bool has_variable(const std::string& name) const;

    /// Get the (cached) oracle, SX or MX
  Function oracle(bool sx = false, bool elim_w = false, bool lifted_calls = false) const;

#ifndef SWIG
  /// Add a variable
  size_t add_variable(const std::string& name, const Variable& var);

  ///@{
  /// Access a variable by name
  Variable& variable(const std::string& name);
  const Variable& variable(const std::string& name) const;
  ///@}

  /// Access a member function or object
  const DaeBuilderInternal* operator->() const;

  /// Access a member function or object
  DaeBuilderInternal* operator->();

  /// Check if a particular cast is allowed
  static bool test_cast(const SharedObjectInternal* ptr);

  /// Get single variable expression by index
  const MX& var(size_t ind) const;

  /// Get variable expressions by index
  std::vector<MX> var(const std::vector<size_t>& ind) const;

  /// Get index of variable
  size_t find(const std::string& name) const;

  /// Get indices of variable
  std::vector<size_t> find(const std::vector<std::string>& name) const;

  /** \brief Get variable name by index */
  const std::string& name(size_t ind) const;

  /** \brief Get variable names by indices */
  std::vector<std::string> name(const std::vector<size_t>& ind) const;

#endif // SWIG
};

} // namespace casadi

#endif // CASADI_DAE_BUILDER_HPP
