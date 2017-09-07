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

#ifndef CASADI_OPTISTACK_HPP
#define CASADI_OPTISTACK_HPP

#include "function.hpp"
#include "callback.hpp"

namespace casadi {

class InternalOptiCallback;

typedef DM native_DM;

class OptiCallback {
public:
  OptiCallback() {
  }
  OptiCallback(const OptiCallback& obj) {
    casadi_error("Callback objects cannot be copied");
  }
  virtual void call(int i) {
    userOut() << "This is a simple callback at iteration" << i << std::endl;
  }
  virtual ~OptiCallback() {}
};

class OptiSol;
class OptiDebug;
/** \brief A simplified interface for NLP modeling/solving

      This is the low-level base class.
      Direct usage of this class is not recommended unless for debugging.
      There are no guaranties API stability


      \date 2017
      \author Joris Gillis, Erik Lambrechts
*/
class CASADI_EXPORT OptiStack : public PrintableObject<OptiStack> {
  friend class InternalOptiCallback;
public:
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

  struct IndexAbstraction {
    IndexAbstraction() : start(0), stop(0) {}
    int start;
    int stop;
  };
  struct MetaCon : IndexAbstraction {
    MetaCon() :  n(1), flipped(false) {}
    MX original;  // original expression
    MX canon; // Canonical expression
    ConstraintType type;
    MX lb;
    MX ub;
    int n;
    bool flipped;
    MX dual_canon;
    MX dual;
  };
  struct MetaVar : IndexAbstraction {
    std::string attribute;
    int n;
    int m;
    VariableType type;
    int count;
    int i;
    Dict extra;
  };

  /** \brief Create Opti Context
  */
  OptiStack();

  /** \brief Create a decision variable (symbol)
  *
  * The order of creation matters.
  * The order will be reflected in the optimization problem.
  * It is not required for decision variables to actualy appear in the optimization problem.
  *
  * \param[in] n number of rows (default 1)
  * \param[in] m number of columnss (default 1)
  * \param[in] attribute: 'full' (default) or 'symmetric'
  */
  MX variable(int n=1, int m=1, const std::string& attribute="full");
  /** \brief Create a parameter (symbol); fixed during optimization
  *
  * The order of creation does not matter.
  * It is not required for parameter to actualy appear in the optimization problem.
  * Parameters that do appear, must be given a value before the problem can be solved.
  *
  * \param[in] n number of rows (default 1)
  * \param[in] m number of columnss (default 1)
  * \param[in] attribute: 'full' (default) or 'symmetric'
  */
  MX parameter(int n=1, int m=1, const std::string& attribute="full");

  /** \brief Set objective
  *
  * Objective must be a scalar. Default objective: 0
  * When method is called multiple times, the last call takes effect
  */
  void minimize(const MX& f);

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
  */
  void subject_to(const MX& g);
  void subject_to(const std::vector<MX>& g);
  /// @}

  /// Clear constraints
  void subject_to();

  /** \brief Set a solver
  *
  * \param[in] solver any of the nlpsol plugins can be used here
  *            In practice, not all nlpsol plugins may be supported yet
  * \param[in] options passed on to nlpsol
  *            No stability can be guaranteed about this part of the API
  */
  void solver(const std::string& solver, const Dict& options=Dict());

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
  */
  void set_value(const MX& x, const DM& v);
  void set_value(const std::vector<MX>& assignments);
  /// @}

  /// Crunch the numbers; solve the problem
  OptiSol solve();

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
  native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const { return x; }
  native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const {
    return DM::nan(x.sparsity());
  }
  /// @}

  /// Copy
  OptiStack copy() const { return *this; }

  /** \brief Get statistics
  *
  * nlpsol stats are passed as-is.
  * No stability can be guaranteed about this part of the API
  */
  Dict stats() const;

  /** \brief Get return status of solver
  *          passed as-is from nlpsol
  * No stability can be guaranteed about this part of the API
  */
  std::string return_status() const;

  /// Get the underlying CasADi solver of the Opti stack
  Function casadi_solver() const;

  /// get assignment expressions for initial values
  std::vector<MX> initial() const;

  /// get assignment expressions for latest values
  std::vector<MX> value_variables() const;
  std::vector<MX> value_parameters() const;

  void callback_class(OptiCallback* callback);
  void callback_class();

  /// return true if expression is only dependant on Opti parameters, not variables
  bool is_parametric(const MX& expr) const;

  /// @{
  /** \brief Get symbols present in expression
  *
  *  Returned vector is ordered according to the order of
  *  variable()/parameter() calls used to create the variables
  */
  std::vector<MX> symvar() const;
  std::vector<MX> symvar(const MX& expr) const;
  std::vector<MX> symvar(const MX& expr, VariableType type) const;
  /// @}

  /// Interpret an expression (for internal use only)
  OptiStack::MetaCon canon_expr(const MX& expr) const;

  /// Get meta-data of symbol (for internal use only)
  MetaVar get_meta(const MX& m) const;

  /// Get meta-data of symbol (for internal use only)
  MetaCon get_meta_con(const MX& m) const;

  /// Set meta-data of an expression
  void set_meta(const MX& m, const MetaVar& meta);

  /// Set meta-data of an expression
  void set_meta_con(const MX& m, const MetaCon& meta);

  /** \brief get the dual variable
  *
  * m must be a constraint expression.
  * The returned value is still a symbolic expression.
  * Use `value` on it to obtain the numerical value.
  */
  MX dual(const MX& m) const;

  void assert_active_symbol(const MX& m) const;

  std::vector<MX> active_symvar(OptiStack::VariableType type) const;
  std::vector<DM> active_values(OptiStack::VariableType type) const;


  void solve_prepare();
  DMDict solve_actual(const DMDict& args);

  DMDict arg() const { return arg_; }
  void res(const DMDict& res);
  std::vector<MX> constraints() const { return g_; }
  MX objective() const { return f_; }

  /// Number of (scalarised) decision variables
  int nx() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("x").size1();
  }

  /// Number of (scalarised) parameters
  int np() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("p").size1();
  }

  /// Number of (scalarised) constraints
  int ng() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("g").size1();
  }

  /// Get all (scalarised) decision variables as a symbolic column vector
  MX x() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("x");
  }

  /// Get all (scalarised) parameters as a symbolic column vector
  MX p() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("p");
  }

  /// Get all (scalarised) constraint expressions as a column vector
  MX g() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("g");
  }

  /// Get objective expression
  MX f() {
    if (problem_dirty()) internal_bake();
    return nlp_.at("f");
  }

  /** \brief Get all (scalarised) dual variables as a symbolic column vector
  *
  * Useful for obtaining the Lagrange Hessian:
  * \verbatim
  * sol.value(hessian(opti.f+opti.lam_g'*opti.g,opti.x)) % MATLAB
  * sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]) # Python
  * \endverbatim
  */
  MX lam_g() {
    if (problem_dirty()) internal_bake();
    return lam_;
  }
  void assert_empty() const;

  ///  Print representation
  void repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

  /// Print description
  void print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const;

  /// Fix the structure of the optimization problem
  void internal_bake();
private:

  static std::map<VariableType, std::string> VariableType2String_;
  std::string variable_type_to_string(VariableType vt) const;

  bool parse_opti_name(const std::string&name, VariableType& vt) const;
  void register_dual(MetaCon& meta);

  /// Set value of symbol
  void set_value_internal(const MX& x, const DM& v, std::vector<DM>& store);

  /** \brief decompose a chain of inequalities
  *
  * a<=b -> [a,b]
  * a<=b<=c [a,b,c]
  *
  * When flipped is set, [a,b,c] corresponds to c>=b>=a
  */
  static std::vector<MX> ineq_unchain(const MX& a, bool& SWIG_OUTPUT(flipped));

  /// Get meta-dat by const-ref
  const MetaVar& meta(const MX& m) const;
  /// Get meta-dat by ref
  MetaVar& meta(const MX& m);

  /// Get meta-dat by const-ref
  const MetaCon& meta_con(const MX& m) const;
  /// Get meta-dat by ref
  MetaCon& meta_con(const MX& m);

  /// Sort symbols according to Opti order
  std::vector<MX> sort(const std::vector<MX>& v) const;

  /// Throw an error
  void assert_has(const MX& m) const;

  bool has(const MX& m) const;

  /// Throw an error
  void assert_has_con(const MX& m) const;

  bool has_con(const MX& m) const;

  // Data members

  /// Map symbols to metadata
  std::map<MXNode*, MetaVar> meta_;
  /// map constraints to metadata
  std::map<MXNode*, MetaCon> meta_con_;

  /// Store references to all symbols
  std::vector<MX> symbols_;

  /// Symbol counter
  int count_;

  int count_var_;
  int count_par_;
  int count_dual_;

  /// Storing latest values for all parameters (including inactive)
  std::vector<DM> values_;
  /// Storing initial values for all variables (including inactive)
  std::vector<DM> initial_;
  /// Storing latest values for all variables (including inactive)
  std::vector<DM> latest_;

  /// Storing initial values for all duals (including inactive)
  std::vector<DM> initial_duals_;
  /// Storing latest values for all duals (including inactive)
  std::vector<DM> latest_duals_;

  /// Is symbol present in problem?
  std::vector<bool> symbol_active_;

  /// Solver
  Function solver_;

  /// Result of solver
  DMDict res_;
  DMDict arg_;
  MXDict nlp_;
  MX lam_;

  /// Bounds helper function: p -> lbg, ubg
  Function bounds_;

  /// Constraints verbatim as passed in with 'subject_to'
  std::vector<MX> g_;

  /// Objective verbatim as passed in with 'minimize'
  MX f_;

  std::vector<OptiCallback*> callbacks_;
  InternalOptiCallback* internal_callback_;
  Function callback_;


  std::string solver_name_;
  Dict solver_options_;
  bool solver_has_callback_;

  void assert_only_opti_symbols(const MX& e) const;
  void assert_only_opti_nondual(const MX& e) const;


  static int instance_count_;
  int instance_number_;

  std::string name_prefix() const;
public:



  bool problem_dirty_;
  void mark_problem_dirty(bool flag=true) { problem_dirty_=flag; mark_solver_dirty(); }
  bool problem_dirty() const { return problem_dirty_; }

  bool solver_dirty_;
  void mark_solver_dirty(bool flag=true) { solver_dirty_=flag; mark_solved(false); }
  bool solver_dirty() const { return solver_dirty_; }

  bool solved_;
  void mark_solved(bool flag=true) { solved_ = flag;}
  bool solved() const { return solved_; }

  void assert_solved() const;
  void assert_baked() const;
};

/** \brief A simplified interface for NLP modeling/solving

      This class offers a view with model description facilities
      The API is guaranteed to be stable.

      \date 2017
      \author Joris Gillis, Erik Lambrechts
*/
class CASADI_EXPORT Opti : private OptiStack {
  public:
    using OptiStack::getDescription;
    using OptiStack::getRepresentation;
    using OptiStack::repr;
    using OptiStack::print;
    using OptiStack::variable;
    using OptiStack::parameter;
    using OptiStack::minimize;
    using OptiStack::subject_to;
    using OptiStack::solver;
    using OptiStack::solve;
    OptiStack debug() { return *this; }
    Opti copy() { return *this; }
    using OptiStack::set_initial;
    using OptiStack::set_value;
    using OptiStack::initial;
    using OptiStack::constraints;
    using OptiStack::objective;
    using OptiStack::callback_class;
    using OptiStack::dual;

    using OptiStack::nx;
    using OptiStack::ng;
    using OptiStack::np;
    using OptiStack::x;
    using OptiStack::p;
    using OptiStack::f;
    using OptiStack::g;
    using OptiStack::lam_g;

    #ifndef SWIGMATLAB
    /** \brief Construct a double inequality
    *
    * Constructs:  lb(p) <= g(x,p) <= ub(p)
    *
    * Python prohibits such syntax directly
    */
    static MX bounded(const MX& lb, const MX& expr, const MX& ub) { return (lb<=expr)<= ub; }
    #endif

};

/** \brief A simplified interface for NLP modeling/solving

      This class offers a view with solution retrieval facilities
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
      \author Joris Gillis, Erik Lambrechts
*/
class CASADI_EXPORT OptiSol : private OptiStack {
  public:
    OptiSol(const OptiStack& opti) : OptiStack(opti) {}
    using OptiStack::getDescription;
    using OptiStack::getRepresentation;
    using OptiStack::repr;
    using OptiStack::print;
    using OptiStack::value;
    using OptiStack::value_variables;
    using OptiStack::value_parameters;
    using OptiStack::stats;
    OptiStack debug() { return *this; }
    Opti opti() { return *reinterpret_cast<Opti*>(this); }
};

//check if we can make splines work
// value_parameters -> parameters
// inline doxygen

} // namespace casadi

#endif // CASADI_OPTI_HPP
