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

#ifndef CASADI_OPTISTACK_HPP
#define CASADI_OPTISTACK_HPP

#include "function.hpp"
#include "callback.hpp"

namespace casadi {

typedef DM native_DM;

class OptiNode;
class OptiSol;
class OptiAdvanced;
class OptiCallback;
/** \brief A simplified interface for NLP modeling/solving

      This class offers a view with model description facilities
      The API is guaranteed to be stable.


      Example NLP:
      \verbatim
      opti = casadi.Opti();

      x = opti.variable();
      y = opti.variable();

      opti.minimize(  (y-x^2)^2   );
      opti.subject_to( x^2+y^2==1 );
      opti.subject_to(     x+y>=1 );

      opti.solver('ipopt');
      sol = opti.solve();

      sol.value(x)
      sol.value(y)
      \endverbatim

      Example parametric NLP:
      \verbatim
      opti = casadi.Opti();

      x = opti.variable(2,1);
      p = opti.parameter();

      opti.minimize(  (p*x(2)-x(1)^2)^2   );
      opti.subject_to( 1<=sum(x)<=2 );

      opti.solver('ipopt');

      opti.set_value(p, 3);
      sol = opti.solve();
      sol.value(x)

      opti.set_value(p, 5);
      sol = opti.solve();
      sol.value(x)
      \endverbatim

      \date 2017
      \author Joris Gillis, Erik Lambrechts, Joel Andersson

    \identifier{16} */
class CASADI_EXPORT Opti
  : public SWIG_IF_ELSE(PrintableCommon, Printable<Opti>),
  public SharedObject {
  friend class InternalOptiCallback;
public:

  /** \brief Create Opti Context
   * 
   * \param[in] problem_type of optimization 'nlp' or 'conic' (default nlp)

      \identifier{17} */
  Opti(const std::string& problem_type="nlp");

  /** \brief Create a decision variable (symbol)
  *
  * The order of creation matters.
  * The order will be reflected in the optimization problem.
  * It is not required for decision variables to actualy appear in the optimization problem.
  *
  * \param[in] n number of rows (default 1)
  * \param[in] m number of columnss (default 1)
  * \param[in] attribute: 'full' (default) or 'symmetric'

      \identifier{18} */
  MX variable(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");
  MX variable(const Sparsity& sp, const std::string& attribute="full");
  MX variable(const MX& symbol, const std::string& attribute="full");

  /** \brief Create a parameter (symbol); fixed during optimization
  *
  * The order of creation does not matter.
  * It is not required for parameter to actualy appear in the optimization problem.
  * Parameters that do appear, must be given a value before the problem can be solved.
  *
  * \param[in] n number of rows (default 1)
  * \param[in] m number of columnss (default 1)
  * \param[in] attribute: 'full' (default) or 'symmetric'

      \identifier{19} */
  MX parameter(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");
  MX parameter(const Sparsity& sp, const std::string& attribute="full");
  MX parameter(const MX& symbol, const std::string& attribute="full");

  /** \brief Set objective
  *
  * Objective must be a scalar. Default objective: 0
  * When method is called multiple times, the last call takes effect

      \identifier{1a} */
  void minimize(const MX& f, double linear_scale=1);

  /// @{
  /** \brief Add constraints
  *
  * Examples:
  * \verbatim
  * \begin{itemize}
  * opti.subject_to( sqrt(x+y) >= 1);
  * opti.subject_to( sqrt(x+y) > 1)}: same as above
  * opti.subject_to( 1<= sqrt(x+y) )}: same as above
  * opti.subject_to( 5*x+y==1 )}: equality
  *
  * Python
  * opti.subject_to([x*y>=1,x==3])
  * opti.subject_to(opti.bounded(0,x,1))
  *
  * MATLAB
  * opti.subject_to({x*y>=1,x==3})
  * opti.subject_to( 0<=x<=1 )
  * \endverbatim
  *
  *
  * Related functionalities:
  *  - opti.lbg,opti.g,opti.ubg represent the vector of flattened constraints
  *  - opti.debug.show_infeasibilities() may be used to inspect which constraints are violated
  *
      \identifier{1b} */
  void subject_to(const MX& g, const Dict& options=Dict());
  void subject_to(const std::vector<MX>& g, const Dict& options=Dict());
  void subject_to(const MX& g, const DM& linear_scale, const Dict& options=Dict());
  void subject_to(const std::vector<MX>& g, const DM& linear_scale, const Dict& options=Dict());
  /// @}

  /// Clear constraints
  void subject_to();

  /** \brief Set a solver
  *
  * \param[in] solver any of the nlpsol plugins can be used here
  *            In practice, not all nlpsol plugins may be supported yet
  * \param[in] options passed on to nlpsol plugin
  *            No stability can be guaranteed about this part of the API
  * \param[in] options to be passed to nlpsol solver
  *            No stability can be guaranteed about this part of the API

      \identifier{1c} */
  void solver(const std::string& solver,
              const Dict& plugin_options=Dict(),
              const Dict& solver_options=Dict());

  /// @{
  /** Set initial guess for decision variables
  * \verbatim
  * opti.set_initial(x, 2)
  * opti.set_initial(10*x(1), 2)
  * \endverbatim
  */
  void set_initial(const MX& x, const DM& v);
  void set_initial(const std::vector<MX>& assignments);
  /// @}

  /// @{
  /** \brief Set value of parameter
  *
  * Each parameter must be given a value before 'solve' can be called

      \identifier{1d} */
  void set_value(const MX& x, const DM& v);
  void set_value(const std::vector<MX>& assignments);
  /// @}

  /// @{
  /** \brief Set domain of a decision variable
  *
  * \param[in] x decision variable
  * \param[in] type 'real', 'integer' (default: real)
  *
  * \verbatim
  * opti.set_domain(x, "real")
  * opti.set_domain(x, "integer")
  * \endverbatim

      \identifier{27t} */
  void set_domain(const MX& x, const std::string& domain);

  /** \brief Set scale of a decision variable
  * 
  * (x-offset)/scale will be used in the optimization problem
  *
  * \param[in] x decision variable
  * \param[in] scale scaling value (default: 1)
  * \param[in] offset scaling value (default: 0)
  *
  * \verbatim
  * opti.set_linear_scale(x, 20)
  * opti.set_linear_scale(x, 20, 273.15)
  * \endverbatim

      \identifier{2bs} */
  void set_linear_scale(const MX& x, const DM& scale, const DM& offset=0);

  /// Crunch the numbers; solve the problem
  OptiSol solve();

  /** \brief Crunch the numbers; solve the problem
   * 
   * Allows the solver to return without error when
   * an iteration or time limit is reached

      \identifier{1e} */
  OptiSol solve_limited();

  /// @{
  /** Obtain value of expression at the current value
  *
  * In regular mode, teh current value is the converged solution
  * In debug mode, the value can be non-converged
  *
  * \param[in] values Optional assignment expressions (e.g. x==3)
  *            to overrule the current value
  */
  native_DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
  native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const;
  native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const;
  /// @}

  /** \brief Get statistics
  *
  * nlpsol stats are passed as-is.
  * No stability can be guaranteed about this part of the API

      \identifier{1f} */
  Dict stats() const;

  /** \brief Get return status of solver

  *          passed as-is from nlpsol
  * No stability can be guaranteed about this part of the API

      \identifier{1g} */
  std::string return_status() const;

  /** \brief get assignment expressions for initial values

      \identifier{266} */
  std::vector<MX> initial() const;

  /** \brief get assignment expressions for latest values

      \identifier{267} */
  std::vector<MX> value_variables() const;
  std::vector<MX> value_parameters() const;

  /** \brief Scale a helper function constructed via opti.x, opti.g, ...

      \identifier{2ci} */
  Function scale_helper(const Function& h) const;

  /** \brief get the dual variable
  *
  * m must be a constraint expression.
  * The returned value is still a symbolic expression.
  * Use `value` on it to obtain the numerical value.

      \identifier{1h} */
  MX dual(const MX& m) const;

  /** \brief Number of (scalarised) decision variables

      \identifier{268} */
  casadi_int nx() const;

  /** \brief Number of (scalarised) parameters

      \identifier{269} */
  casadi_int np() const;

  /** \brief Number of (scalarised) constraints

      \identifier{26a} */
  casadi_int ng() const;

  /** \brief Get all (scalarised) decision variables as a symbolic column vector

      \identifier{26b} */
  MX x() const;

  /** \brief Get all (scalarised) parameters as a symbolic column vector

      \identifier{26c} */
  MX p() const;

  /** \brief Get all (scalarised) constraint expressions as a column vector

      \identifier{26d} */
  MX g() const;

  /** \brief Get objective expression

      \identifier{26e} */
  MX f() const;

  /** \brief Get all (scalarised) bounds on constraints as a column vector

      \identifier{26f} */
  MX lbg() const;
  MX ubg() const;

  DM x_linear_scale() const;
  DM x_linear_scale_offset() const;
  DM g_linear_scale() const;
  double f_linear_scale() const;

  /** \brief Get all (scalarised) dual variables as a symbolic column vector
  *
  * Useful for obtaining the Lagrange Hessian:
  * \verbatim
  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * \endverbatim

      \identifier{1i} */
  MX lam_g() const;

  /// @{
  /** \brief Create a CasADi Function from the Opti solver
   * 
   * \param[in] name Name of the resulting CasADi Function
   * \param[in] args List of parameters and decision/dual variables
   *                (which can be given an initial guess) with the resulting Function
   * \param[in] res List of expressions that will get evaluated at the optimal solution
   * \param[in] opts Standard CasADi Funcion options

      \identifier{1j} */
  Function to_function(const std::string& name,
      const std::vector<MX>& args, const std::vector<MX>& res,
      const Dict& opts = Dict());

  Function to_function(const std::string& name,
      const std::vector<MX>& args, const std::vector<MX>& res,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts = Dict());

  Function to_function(const std::string& name,
      const std::map<std::string, MX>& dict,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts = Dict());
  /// @}

  #ifndef SWIGMATLAB
  /** \brief Construct a double inequality
  *
  * Constructs:  lb(p) <= g(x,p) <= ub(p)
  *
  * Python prohibits such syntax directly

      \identifier{1k} */
  static MX bounded(const MX& lb, const MX& expr, const MX& ub) { return (lb<=expr)<= ub; }
  #endif

  /** \brief Get a copy with advanced functionality
   *
   * You get access to more methods, but you have no guarantees about API stability
   *
   * The copy is effectively a deep copy:
   * Updating the state of the copy does not update the original.
   *
      \identifier{1l} */
  OptiAdvanced debug() const;

  /** \brief Get a copy with advanced functionality
   *
   * You get access to more methods, but you have no guarantees about API stability
   *
   * The copy is effectively a deep copy:
   * Updating the state of the copy does not update the original.
   *
      \identifier{1m} */
  OptiAdvanced advanced() const;

  /** \brief Get a copy of the
   *
   * The copy is effectively a deep copy:
   * Updating the state of the copy does not update the original.
   *
      \identifier{1n} */
  Opti copy() const;

  /** \brief add user data

  * Add arbitrary data in the form of a dictionary to symbols
  * or constraints

      \identifier{1o} */
  void update_user_dict(const MX& m, const Dict& meta);
  void update_user_dict(const std::vector<MX>& m, const Dict& meta);
  /// Get user data
  Dict user_dict(const MX& m) const;

  /// Readable name of the class
  std::string type_name() const { return "Opti"; }

  ///  Print representation
  void disp(std::ostream& stream, bool more=false) const;

  /// Get string representation
  std::string get_str(bool more=false) const;

  ///@{
  /** \brief Helper methods for callback()
   *
   * Do not use directly.

      \identifier{1p} */
  void callback_class(OptiCallback* callback);
  void callback_class();
  ///@}

#ifndef SWIG
  Opti(const Opti& x);

  /** \brief Destructor

      \identifier{1q} */
  ~Opti() {}

  static Opti create(OptiNode* node);
  /// \cond INTERNAL
  ///@{
  /** \brief  Access a member of the node

      \identifier{1r} */
  OptiNode* operator->();

  /** \brief  Const access a member of the node

      \identifier{1s} */
  const OptiNode* operator->() const;
  ///@}
  /// \endcond

  Opti(OptiNode* node);

#endif // SWIG

};

  enum ConstraintType {
    OPTI_GENERIC_EQUALITY,  // g1(x,p) == g2(x,p)
    OPTI_GENERIC_INEQUALITY,  // g1(x,p) <= g2(x,p)
    OPTI_EQUALITY, // g(x,p) == bound(p)
    OPTI_INEQUALITY,  // g(x,p) <= bound(p)
    OPTI_DOUBLE_INEQUALITY,  // lb(p) <= g(x,p) <= ub(p)
    OPTI_PSD, // A(x,p) >= b(p)
    OPTI_UNKNOWN};
  enum VariableType {
    OPTI_VAR, // variable
    OPTI_PAR,  // parameter
    OPTI_DUAL_G // dual
  };
  enum DomainType {
    OPTI_DOMAIN_REAL,
    OPTI_DOMAIN_INTEGER
  };


  struct IndexAbstraction {
    IndexAbstraction() : start(0), stop(0) {}
    casadi_int start;
    casadi_int stop;
  };
  struct MetaCon : IndexAbstraction {
    MetaCon() :  n(1), flipped(false), linear_scale(1) {}
    MX original;  // original expression
    MX canon; // Canonical expression
    ConstraintType type;
    MX lb;
    MX ub;
    casadi_int n;
    bool flipped;
    MX dual_canon;
    MX dual;
    Dict extra;
    DM linear_scale;
  };
  struct MetaVar : IndexAbstraction {
    std::string attribute;
    casadi_int n;
    casadi_int m;
    VariableType type;
    DomainType domain;
    casadi_int count;
    casadi_int i;
    casadi_int active_i;
    Dict extra;
  };


class OptiCallback {
public:
  OptiCallback() {
  }
  OptiCallback(const OptiCallback& obj) {
    casadi_error("Callback objects cannot be copied");
  }
  virtual void call(casadi_int i) {
    uout() << "This is a simple callback at iteration" << i << std::endl;
  }
  virtual ~OptiCallback() {}
};

class CASADI_EXPORT OptiAdvanced : public Opti {
  friend class InternalOptiCallback;
public:

  OptiAdvanced(const Opti& x);

  /** \brief Destructor

      \identifier{1t} */
  ~OptiAdvanced() {}


  /// Get the underlying CasADi solver of the Opti stack
  Function casadi_solver() const;

  /// return true if expression is only dependant on Opti parameters, not variables
  bool is_parametric(const MX& expr) const;

  /// @{
  /** \brief Get symbols present in expression
  *
  *  Returned vector is ordered according to the order of
  *  variable()/parameter() calls used to create the variables

      \identifier{1u} */
  std::vector<MX> symvar() const;
  std::vector<MX> symvar(const MX& expr) const;
  std::vector<MX> symvar(const MX& expr, VariableType type) const;
  /// @}

  /// Interpret an expression (for internal use only)
  MetaCon canon_expr(const MX& expr) const;

  /// Get meta-data of symbol (for internal use only)
  MetaVar get_meta(const MX& m) const;

  /// Get meta-data of symbol (for internal use only)
  MetaCon get_meta_con(const MX& m) const;

  /// Set meta-data of an expression
  void set_meta(const MX& m, const MetaVar& meta);

  /// Set meta-data of an expression
  void set_meta_con(const MX& m, const MetaCon& meta);

  void assert_active_symbol(const MX& m) const;

  std::vector<MX> active_symvar(VariableType type) const;
  std::vector<DM> active_values(VariableType type) const;

  MX x_lookup(casadi_index i) const;
  MX g_lookup(casadi_index i) const;

  casadi_index g_index_reduce_g(casadi_index i) const;
  casadi_index g_index_reduce_x(casadi_index i) const;
  casadi_index g_index_unreduce_g(casadi_index i) const;

  std::string x_describe(casadi_index i, const Dict& opts=Dict()) const;
  std::string g_describe(casadi_index i, const Dict& opts=Dict()) const;
  std::string describe(const MX& x, casadi_index indent=0, const Dict& opts=Dict()) const;

  void show_infeasibilities(double tol=0, const Dict& opts=Dict()) const;

  void solve_prepare();
  DMDict solve_actual(const DMDict& args);

  DMDict arg() const;
  void res(const DMDict& res);
  DMDict res() const;
  std::vector<MX> constraints() const;
  MX objective() const;

  OptiAdvanced baked_copy() const;

  void assert_empty() const;


  /// Fix the structure of the optimization problem
  void bake();

  bool problem_dirty_;
  void mark_problem_dirty(bool flag=true);
  bool problem_dirty() const;

  bool solver_dirty_;
  void mark_solver_dirty(bool flag=true);
  bool solver_dirty() const;

  bool solved_;
  void mark_solved(bool flag=true);
  bool solved() const;

  void assert_solved() const;
  void assert_baked() const;

  casadi_int instance_number() const;

protected:
  OptiAdvanced() {}
};

/** \brief A simplified interface for NLP modeling/solving

      This class offers a view with solution retrieval facilities
      The API is guaranteed to be stable.


      \date 2017
      \author Joris Gillis, Erik Lambrechts

    \identifier{1v} */
class CASADI_EXPORT OptiSol : public SWIG_IF_ELSE(PrintableCommon, Printable<OptiAdvanced>) {
  friend class OptiNode;
  public:
    std::string type_name() const {return "OptiSol";}
    void disp(std::ostream& stream, bool more=false) const;
    std::string get_str(bool more=false) const;
    /// @{
    /** Obtain value of expression at the current value
    *
    * In regular mode, teh current value is the converged solution
    * In debug mode, the value can be non-converged
    *
    * \param[in] values Optional assignment expressions (e.g. x==3)
    *            to overrule the current value
    */
    native_DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const;
    /// @}

    /// get assignment expressions for the optimal solution
    std::vector<MX> value_variables() const;
    std::vector<MX> value_parameters() const;

    /** \brief Get statistics
    *
    * nlpsol stats are passed as-is.
    * No stability can be guaranteed about this part of the API

        \identifier{1w} */
    Dict stats() const;

    Opti opti() const { return optistack_; } // NOLINT(cppcoreguidelines-slicing)

  protected:
    OptiSol(const Opti& opti);
    OptiAdvanced optistack_;
};


} // namespace casadi

#endif // CASADI_OPTI_HPP
