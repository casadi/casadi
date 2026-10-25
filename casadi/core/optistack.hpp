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

    Opti collects decision variables, parameters, an objective and constraints written in
    natural syntax, builds an nlpsol (or qpsol) instance from them and gives numerical
    access to the solution. The bookkeeping of flattening variables into one vector and
    constraints into bounds is hidden.

    Workflow:
      - Declare symbols with variable() and parameter(); they are plain MX symbols.
      - State the problem with minimize() and subject_to().
      - Pick a numerical back-end with solver().
      - Give parameters a value with set_value(), optionally a guess with set_initial().
      - Call solve(); it returns an OptiSol snapshot whose value() evaluates expressions.

    Opti is stateful: every call above mutates the instance, and the problem is rebuilt
    lazily at the next solve(). An OptiSol is a copy taken at the time of solve(), so later
    changes to the Opti do not affect it. Only symbols created by this instance may appear
    in its expressions; raw MX.sym symbols or symbols of another Opti raise an error.

    The constructor argument selects the problem class: "nlp" (default, solved with nlpsol)
    or "conic" (quadratic programs and SOCPs, solved with qpsol).

    An iteration callback can be registered with opti.callback(fh); fh is called with the
    iteration number after every solver iteration, and opti.debug.value() inside it returns
    the current iterate. Calling opti.callback() without argument removes it.

    \lang{python}
    \doctest
    >>> opti = Opti()
    >>> x = opti.variable()
    >>> y = opti.variable()
    >>> opti.minimize((y-x**2)**2)
    >>> opti.subject_to(x**2+y**2==1)
    >>> opti.subject_to(x+y>=1)
    >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
    >>> sol = opti.solve()
    >>> sol.value(x), sol.value(y)
    (0.7862, 0.6180)
    >>> sol.value(opti.f)
    0
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> opti = Opti();
    >> x = opti.variable(); y = opti.variable();
    >> opti.minimize((y-x^2)^2);
    >> opti.subject_to(x^2+y^2==1);
    >> opti.subject_to(x+y>=1);
    >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
    >> sol = opti.solve();
    >> sol.value(x), sol.value(y)
    0.7862
    0.6180
    >> sol.value(opti.f)
    0
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> opti = ca.Opti();
    julia> x = ca.variable(opti);
    julia> y = ca.variable(opti);
    julia> ca.minimize(opti, (y - x^2)^2)
    julia> ca.subject_to(opti, x^2 + y^2 == 1)
    julia> ca.subject_to(opti, x + y >= 1)
    julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                      Dict("print_level" => 0, "sb" => "yes"))
    julia> sol = ca.solve(opti);
    julia> ca.value(sol, x), ca.value(sol, y)
    (0.7862, 0.6180)
    julia> ca.value(sol, ca.f(opti))
    0
    \enddoctest
    \endlang

    Parametric problems are re-solved by changing parameter values:

    \lang{python}
    \doctest
    >>> opti = Opti()
    >>> x = opti.variable(2)
    >>> p = opti.parameter()
    >>> opti.minimize(sumsqr(x-p))
    >>> opti.subject_to(1 <= (sum1(x) <= 2))
    >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
    >>> opti.set_value(p, 3)
    >>> sol = opti.solve()
    >>> sol.value(sum1(x))
    2.0
    >>> opti.set_value(p, 0)
    >>> sol = opti.solve()
    >>> sol.value(sum1(x))
    1.0
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> opti = Opti();
    >> x = opti.variable(2); p = opti.parameter();
    >> opti.minimize(sumsqr(x-p));
    >> opti.subject_to(1 <= sum1(x) <= 2);
    >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
    >> opti.set_value(p, 3); sol = opti.solve();
    >> sol.value(sum1(x))
    2
    >> opti.set_value(p, 0); sol = opti.solve();
    >> sol.value(sum1(x))
    1
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> opti = ca.Opti();
    julia> x = ca.variable(opti, 2);
    julia> p = ca.parameter(opti);
    julia> ca.minimize(opti, ca.sumsqr(x - p))
    julia> ca.subject_to(opti, 1 <= (sum(x) <= 2))
    julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                      Dict("print_level" => 0, "sb" => "yes"))
    julia> ca.set_value(opti, p, 3)
    julia> sol = ca.solve(opti);
    julia> ca.value(sol, sum(x))
    2.0
    julia> ca.set_value(opti, p, 0)
    julia> sol = ca.solve(opti);
    julia> ca.value(sol, sum(x))
    1.0
    \enddoctest
    \endlang

    An iteration callback:

    \lang{python}
    \doctest
    >>> opti = Opti()
    >>> x = opti.variable()
    >>> opti.minimize((x-1)**2)
    >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
    >>> history = []
    >>> opti.callback(lambda i: history.append((i, opti.debug.value(x))))
    >>> sol = opti.solve()
    >>> history
    [(0, 0.0), (1, 1.0)]
    >>> opti.callback()
    >>> sol = opti.solve()
    >>> len(history)
    2
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^2);
    >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
    >> opti.callback(@(i) fprintf('%d %.3f\n', i, opti.debug.value(x)));
    >> sol = opti.solve();
    0 0.000
    1 1.000
    >> opti.callback();
    >> sol = opti.solve();
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> opti = ca.Opti();
    julia> x = ca.variable(opti);
    julia> ca.minimize(opti, (x - 1)^2)
    julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                      Dict("print_level" => 0, "sb" => "yes"))
    julia> history = Tuple{Int,Float64}[];
    julia> ca.callback(opti, i -> push!(history, (i, ca.value(ca.debug(opti), x))))
    julia> sol = ca.solve(opti);
    julia> history
    2-element Vector{Tuple{Int64, Float64}}:
     (0, 0.0)
     (1, 1.0)
    julia> ca.callback(opti)
    julia> sol = ca.solve(opti);
    julia> length(history)
    2
    \enddoctest
    \endlang

    \sa OptiSol, OptiAdvanced, nlpsol, qpsol

    \date 2017
    \author Joris Gillis, Erik Lambrechts, Joel Andersson

    \identifier{16} */
class CASADI_EXPORT Opti
  : public SWIG_IF_ELSE(PrintableCommon, Printable<Opti>),
  public SharedObject {
  friend class InternalOptiCallback;
public:

  /** \brief Create an empty Opti context

      \param problem_type "nlp" (default): general nonlinear program, solved with an nlpsol
             plugin. "conic": quadratic program or SOCP, solved with a qpsol plugin; the
             objective must then be quadratic and the constraints linear or second-order cone.
             Any other string raises an error.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> opti.nx, opti.np, opti.ng
      (0, 0, 0)
      >>> opti = Opti("conic")
      >>> x = opti.variable(2)
      >>> opti.minimize(sumsqr(x-1))
      >>> opti.subject_to(x[0]+x[1] <= 1)
      >>> opti.solver("qrqp", {"print_time": False, "print_iter": False,
      ...                      "print_header": False, "print_info": False})
      >>> opti.solve().value(x)
      array([0.5, 0.5])
      >>> try:
      ...   Opti("lp")
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Specified problem type 'lp'unknown. Choose 'nlp' (default) or 'conic'.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); fprintf('%d %d %d\n', opti.nx, opti.np, opti.ng)
      0 0 0
      >> opti = Opti('conic'); x = opti.variable(2);
      >> opti.minimize(sumsqr(x-1)); opti.subject_to(x(1)+x(2) <= 1);
      >> opts = struct('print_time', false, 'print_iter', false, 'print_header', false);
      >> opts.print_info = false; opti.solver('qrqp', opts);
      >> sol = opti.solve(); disp(sol.value(x)')
      0.5000 0.5000
      >> try, Opti('lp'); catch e, disp(contains(e.message, 'Specified problem type')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> ca.nx(opti), ca.np(opti), ca.ng(opti)
      (0, 0, 0)
      julia> opti = ca.Opti("conic");
      julia> x = ca.variable(opti, 2);
      julia> ca.minimize(opti, ca.sumsqr(x - 1))
      julia> ca.subject_to(opti, x[1] + x[2] <= 1)
      julia> ca.solver(opti, "qrqp", Dict("print_time" => false, "print_iter" => false,
                                           "print_header" => false, "print_info" => false))
      julia> ca.value(ca.solve(opti), x)
      2-element Vector{Float64}:
       0.5
       0.5
      julia> try
                 ca.Opti("lp")
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Specified problem type 'lp'unknown. Choose 'nlp' (default) or 'conic'.
      \enddoctest
      \endlang

      \identifier{17} */
  Opti(const std::string& problem_type="nlp");

  /** \brief Create a decision variable (symbol)

      Returns an MX symbol of the requested shape, registered with this Opti instance.
      Variables are stacked into the solver's x vector in order of creation, so the order
      of variable() calls determines the layout of x and of lam_x. A variable that does not
      appear in the objective or the constraints is simply left out of the problem, but
      asking for its value after solve() raises an error.

      With attribute "symmetric" (square shapes only) the returned expression is a
      symmetric matrix built from n*(n+1)/2 underlying scalars.

      The default initial guess is 0; change it with set_initial().

      \param n Number of rows (default 1)
      \param m Number of columns (default 1)
      \param attribute "full" (default) or "symmetric"; anything else raises an error
      \return A dense MX symbol of shape n-by-m

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> X = opti.variable(3, 2)
      >>> x.shape, X.shape, X.is_symbolic()
      ((1, 1), (3, 2), True)
      >>> S = opti.variable(2, 2, "symmetric")
      >>> S.shape, S.is_symbolic()
      ((2, 2), False)
      >>> opti.minimize(sumsqr(X) + x**2 + sumsqr(S))
      >>> opti.nx
      10
      >>> try:
      ...   opti.variable(2, 3, "symmetric")
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      You specified attribute 'symmetric', while matrix is not even square, but 2-by-3.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); X = opti.variable(3, 2);
      >> size(x), size(X), X.is_symbolic()
      1 1
      3 2
      1
      >> S = opti.variable(2, 2, 'symmetric'); disp(size(S)), disp(S.is_symbolic())
      2 2
      0
      >> opti.minimize(sumsqr(X) + x^2 + sumsqr(S)); disp(opti.nx)
      10
      >> try, opti.variable(2, 3, 'symmetric'); catch e, disp(contains(e.message, 'not even')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> X = ca.variable(opti, 3, 2);
      julia> size(x), size(X), ca.is_symbolic(X)
      ((1, 1), (3, 2), true)
      julia> S = ca.variable(opti, 2, 2, "symmetric");
      julia> size(S), ca.is_symbolic(S)
      ((2, 2), false)
      julia> ca.minimize(opti, ca.sumsqr(X) + x^2 + ca.sumsqr(S))
      julia> ca.nx(opti)
      10
      julia> try
                 ca.variable(opti, 2, 3, "symmetric")
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      You specified attribute 'symmetric', while matrix is not even square, but 2-by-3.
      \enddoctest
      \endlang

      \identifier{18} */
  MX variable(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");
  /// Create a decision variable with a given sparsity pattern; only "full" is supported
  MX variable(const Sparsity& sp, const std::string& attribute="full");
  /// Register an existing MX symbol as a decision variable of this instance
  MX variable(const MX& symbol, const std::string& attribute="full");

  /** \brief Create a parameter (symbol); fixed during optimization

      Returns an MX symbol that is held constant by the solver. Each parameter that appears
      in the problem must be given a numerical value with set_value() before solve(),
      otherwise solve() raises an error naming the parameter. The order of creation does
      not matter for the result.

      \param n Number of rows (default 1)
      \param m Number of columns (default 1)
      \param attribute "full" (default) or "symmetric"
      \return A dense MX symbol of shape n-by-m

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter(2)
      >>> opti.minimize(sumsqr(x-p))
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print(str(e).splitlines()[-2])
      You have forgotten to assign a value to a parameter ('set_value'), or have set it to NaN/Inf:
      >>> opti.set_value(p, [1, 3])
      >>> opti.solve().value(x)
      2.0
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter(2);
      >> opti.minimize(sumsqr(x-p));
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> try, opti.solve(); catch e, disp(contains(e.message, 'forgotten to assign a value')); end
      1
      >> opti.set_value(p, [1; 3]); sol = opti.solve(); disp(sol.value(x))
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti, 2);
      julia> ca.minimize(opti, ca.sumsqr(x - p))
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
                 println(split(sprint(showerror, e), "\n")[end-1])
             end
      You have forgotten to assign a value to a parameter ('set_value'), or have set it to NaN/Inf:
      julia> ca.set_value(opti, p, [1, 3])
      julia> ca.value(ca.solve(opti), x)
      2.0
      \enddoctest
      \endlang

      \identifier{19} */
  MX parameter(casadi_int n=1, casadi_int m=1, const std::string& attribute="full");
  /// Create a parameter with a given sparsity pattern
  MX parameter(const Sparsity& sp, const std::string& attribute="full");
  /// Register an existing MX symbol as a parameter of this instance
  MX parameter(const MX& symbol, const std::string& attribute="full");

  /** \brief Set the objective

      The objective must be a scalar expression in the variables and parameters of this
      instance; a non-scalar raises an error. The default objective is 0. Each call replaces
      the previous objective. Dual variables (opti.dual, opti.lam_g) may not appear in it.

      \param f Scalar MX expression to minimize
      \param linear_scale The solver sees f divided by this factor (default 1)

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.minimize(sumsqr(x))
      >>> opti.minimize(x[0]**2, 2)
      >>> print(opti.f, opti.f_linear_scale)
      sq(opti..._x_1[0]) 2.0
      >>> try:
      ...   opti.minimize(x)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Objective must be scalar, got 2x1.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2);
      >> opti.minimize(sumsqr(x)); opti.minimize(x(1)^2, 2);
      >> opti.f, opti.f_linear_scale
      sq(opti..._x_1[0])
      2
      >> try, opti.minimize(x); catch e, disp(contains(e.message, 'Objective must be scalar')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.minimize(opti, ca.sumsqr(x))
      julia> opti.minimize(x[1]^2, 2)
      julia> print(ca.f(opti), " ", ca.f_linear_scale(opti))
      sq(opti..._x_1[0]) 2.0
      julia> try
                 ca.minimize(opti, x)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Objective must be scalar, got 2x1.
      \enddoctest
      \endlang

      \identifier{1a} */
  void minimize(const MX& f, double linear_scale=1);

  /// @{
  /** \brief Add constraints

      \p g is a comparison expression built with ==, <=, >= (< and > behave as <= and >=).
      Vector expressions give one constraint per element; both sides may be vectors of the
      same shape or a scalar that is broadcast. Each call appends to the constraint set.

      Recognized forms and how they reach the solver:
        - g(x,p) <= b(p), b(p) <= g(x,p), g(x,p) == b(p): one side free of decision variables
          becomes the bound (lbg/ubg) and the other side the constraint expression g.
        - lb(p) <= (g(x,p) <= ub(p)): a chained inequality whose outer sides are free of
          decision variables becomes a double-bounded constraint lb <= g <= ub (one row).
          Both bounds may be vectors. In Python the inner comparison must be parenthesized,
          because the plain chain 0 <= x <= 1 asks for the truth value of an MX and raises
          an error; bounded() builds the same expression.
        - g1(x,p) <= g2(x,p), g1(x,p) == g2(x,p): both sides depend on variables, giving
          g1-g2 <= 0 or g1-g2 == 0.
      A double-bounded constraint on a pure decision variable, such as -1 <= (x <= 1), is
      still passed as a constraint row; set the nlpsol option "detect_simple_bounds" in
      solver() to have it moved to the variable bounds lbx/ubx.

      Errors: a constant comparison, an empty or fully sparse expression, a symbol not
      declared by this instance, and dual variables all raise an error immediately. A square
      matrix comparison A >= B is interpreted as a semidefinite constraint, which "nlp" does
      not support: vectorize the matrix with vec(A) for an elementwise constraint. An
      expression without comparison operator is accepted but rejected at solve().

      \param g Constraint expression, or a list of them (each added in turn)
      \param linear_scale Scalar or same-shape matrix; the solver sees g divided by it
      \param options Dict with optional entries "meta" (dict stored in the user_dict of the
             constraint) and "stacktrace"; any other key raises an error

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(3)
      >>> y = opti.variable()
      >>> p = opti.parameter()
      >>> opti.subject_to(x >= 0)
      >>> opti.subject_to(x[0] + x[1] == p)
      >>> opti.subject_to(-1 <= (y <= 1))
      >>> opti.subject_to([y >= x[2], y*x[0] >= 1])
      >>> opti.ng
      7
      >>> print(opti.lbg, opti.ubg)
      vertcat(zeros(3x1), opti..._p_1, -1, -inf, 1) vertcat(inf(3x1), opti..._p_1, 1, 0, inf)
      >>> try:
      ...   opti.subject_to(0 <= y <= 1)
      ... except Exception as e:
      ...   print("Can only determine truth value of a numeric MX" in str(e))
      True
      >>> try:
      ...   opti.subject_to(MX(2) <= 3)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      You passed a constant to `subject_to`. You need a symbol to form a constraint.
      >>> try:
      ...   opti.subject_to(MX.sym("w") <= 3)
      ... except Exception as e:
      ...   print("MX symbol 'w' of shape 1x1, declared outside of Opti" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(3); y = opti.variable(); p = opti.parameter();
      >> opti.subject_to(x >= 0);
      >> opti.subject_to(x(1) + x(2) == p);
      >> opti.subject_to(-1 <= y <= 1);
      >> opti.subject_to({y >= x(3), y*x(1) >= 1});
      >> opti.ng
      7
      >> opti.lbg, opti.ubg
      vertcat(zeros(3x1), opti..._p_1, -1, -inf, 1)
      vertcat(inf(3x1), opti..._p_1, 1, 0, inf)
      >> try, opti.subject_to(MX(2) <= 3); catch e, disp(contains(e.message, 'a constant')); end
      1
      >> w = MX.sym('w');
      >> try, opti.subject_to(w <= 3); catch e, disp(contains(e.message, 'outside of Opti')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 3);
      julia> y = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.subject_to(opti, x >= 0)
      julia> ca.subject_to(opti, x[1] + x[2] == p)
      julia> ca.subject_to(opti, -1 <= (y <= 1))
      julia> ca.subject_to(opti, [y >= x[3], y*x[1] >= 1])
      julia> ca.ng(opti)
      7
      julia> print(ca.lbg(opti), " ", ca.ubg(opti))
      vertcat(zeros(3x1), opti..._p_1, -1, -inf, 1) vertcat(inf(3x1), opti..._p_1, 1, 0, inf)
      julia> try
                 ca.subject_to(opti, 0 <= y <= 1)
             catch e
                 println(typeof(e))
             end
      TypeError
      julia> try
                 ca.subject_to(opti, ca.MX(2) <= 3)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      You passed a constant to `subject_to`. You need a symbol to form a constraint.
      julia> try
                 ca.subject_to(opti, ca.MX.sym("w") <= 3)
             catch e
                 msg = "MX symbol 'w' of shape 1x1, declared outside of Opti"
                 println(occursin(msg, sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      Elementwise constraints on matrices need vec(); options attach metadata:

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> A = opti.variable(2, 2)
      >>> opti.subject_to(A >= 0)
      >>> opti.minimize(sumsqr(A))
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print("Psd constraints not implemented yet" in str(e))
      True
      >>> opti.subject_to()
      >>> con = vec(A) >= 0
      >>> opti.subject_to(con, {"meta": {"label": "nonneg"}})
      >>> opti.ng, opti.user_dict(con)["label"]
      (4, 'nonneg')
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); A = opti.variable(2, 2);
      >> opti.subject_to(A >= 0); opti.minimize(sumsqr(A));
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> try, opti.solve(); catch e, disp(contains(e.message, 'Psd constraints')); end
      1
      >> opti.subject_to(); con = vec(A) >= 0;
      >> opti.subject_to(con, struct('meta', struct('label', 'nonneg')));
      >> d = opti.user_dict(con); fprintf('%d %s\n', opti.ng, d.label)
      4 nonneg
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> A = ca.variable(opti, 2, 2);
      julia> ca.subject_to(opti, A >= 0)
      julia> ca.minimize(opti, ca.sumsqr(A))
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
                 println(occursin("Psd constraints not implemented yet", sprint(showerror, e)))
             end
      true
      julia> ca.subject_to(opti)
      julia> con = vec(A) >= 0;
      julia> ca.subject_to(opti, con, Dict("meta" => Dict("label" => "nonneg")))
      julia> ca.ng(opti), ca.user_dict(opti, con)["label"]
      (4, "nonneg")
      \enddoctest
      \endlang

      \sa bounded, subject_to(), dual, lbg, ubg, g

      \identifier{1b} */
  void subject_to(const MX& g, const Dict& options=Dict());
  void subject_to(const std::vector<MX>& g, const Dict& options=Dict());
  void subject_to(const MX& g, const DM& linear_scale, const Dict& options=Dict());
  void subject_to(const std::vector<MX>& g, const DM& linear_scale, const Dict& options=Dict());
  /// @}

  /** \brief Clear all constraints

      Removes every constraint added so far, including their dual variables, so the problem
      can be restated (typical when the constraint set changes between solves). The objective,
      variables, parameters, values and initial guesses are kept.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> con = x >= 1
      >>> opti.subject_to(con)
      >>> opti.minimize(x**2)
      >>> opti.ng
      1
      >>> opti.subject_to()
      >>> opti.ng
      0
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> try:
      ...   sol.value(opti.dual(con))
      ... except Exception as e:
      ...   print(str(e).splitlines()[-2])
      This expression has a dual for a constraint that is not given to Opti:
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); con = x >= 1;
      >> opti.subject_to(con); opti.minimize(x^2); disp(opti.ng)
      1
      >> opti.subject_to(); disp(opti.ng)
      0
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> try, sol.value(opti.dual(con)); catch e, disp(contains(e.message, 'not given')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> con = x >= 1;
      julia> ca.subject_to(opti, con)
      julia> ca.minimize(opti, x^2)
      julia> ca.ng(opti)
      1
      julia> ca.subject_to(opti)
      julia> ca.ng(opti)
      0
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> try
                 ca.value(sol, ca.dual(opti, con))
             catch e
                 println(split(sprint(showerror, e), "\n")[end-1])
             end
      This expression has a dual for a constraint that is not given to Opti:
      \enddoctest
      \endlang

      \identifier{2ms} */
  void subject_to();

  /** \brief Select the numerical solver and its options

      Must be called before solve(). The solver is (re)constructed lazily at the next solve()
      when this method or the problem structure changed.

      \param solver Name of an nlpsol plugin ("ipopt", "sqpmethod", "bonmin", ...) or,
             for a "conic" Opti, of a qpsol plugin ("qrqp", "qpoases", "osqp", ...)
      \param plugin_options Dict passed as-is to nlpsol/qpsol: CasADi-level options such as
             "print_time", "expand", "detect_simple_bounds", "error_on_fail",
             "iteration_callback", and plugin options that are registered as CasADi options
             (e.g. "qpsol" and "print_iteration" of sqpmethod, "print_iter" of qrqp)
      \param solver_options Dict of the back-end's native options; it is nested under the
             key \p solver, so {"max_iter": 50} here equals {"ipopt": {"max_iter": 50}} in
             \p plugin_options. Opti also sets the nlpsol options "equality" and "discrete"
             unless given explicitly.

      No stability is guaranteed for option names; they belong to the plugins.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-2)**2)
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print("You must call 'solver' on the Opti stack" in str(e))
      True
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> opti.solve().value(x)
      2.0
      >>> opti.solver("sqpmethod", {"print_time": False, "print_header": False,
      ...   "print_iteration": False, "print_status": False, "qpsol": "qrqp",
      ...   "qpsol_options": {"print_iter": False, "print_header": False, "print_info": False}})
      >>> opti.solve().value(x)
      2.0
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-2)^2);
      >> try, opti.solve(); catch e, disp(contains(e.message, 'You must call ''solver''')); end
      1
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); disp(sol.value(x))
      2
      >> qp = struct('print_iter', false, 'print_header', false, 'print_info', false);
      >> opts = struct('print_time', false, 'print_header', false, 'print_iteration', false);
      >> opts.print_status = false; opts.qpsol = 'qrqp'; opts.qpsol_options = qp;
      >> opti.solver('sqpmethod', opts); sol = opti.solve(); disp(sol.value(x))
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 2)^2)
      julia> try
                 ca.solve(opti)
             catch e
                 println(occursin("You must call 'solver' on the Opti stack", sprint(showerror, e)))
             end
      true
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> ca.value(ca.solve(opti), x)
      2.0
      julia> qp_opts = Dict("print_iter" => false, "print_header" => false, "print_info" => false);
      julia> ca.solver(opti, "sqpmethod", Dict("print_time" => false, "print_header" => false,
                 "print_iteration" => false, "print_status" => false, "qpsol" => "qrqp",
                 "qpsol_options" => qp_opts))
      julia> ca.value(ca.solve(opti), x)
      2.0
      \enddoctest
      \endlang

      \sa nlpsol, qpsol, solve

      \identifier{1c} */
  void solver(const std::string& solver,
              const Dict& plugin_options=Dict(),
              const Dict& solver_options=Dict());

  /// @{
  /** \brief Set the initial guess of decision variables or dual variables

      The guess is used as x0 (and lam_g0 for opti.lam_g or opti.dual(con)) at the next
      solve(); the default is 0. \p x may be a variable, a slice of one, or a linear mapping
      of variables such as 10*x[0]; a nonlinear expression raises an error. Parameters are
      refused with a hint to use set_value(). The list form takes assignment expressions
      var==value, as returned by value_variables() or initial(), which makes warm starting a
      one-liner. Calling this method invalidates value() until the next solve().

      \param x Variable (or linear mapping of variables) to initialize
      \param v Numerical value; scalars are broadcast to the shape of \p x

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> p = opti.parameter()
      >>> opti.set_initial(x, 1)
      >>> opti.set_initial(10*x[0], 2)
      >>> print(opti.initial())
      [MX(([0.2, 1]==opti..._x_1))]
      >>> try:
      ...   opti.set_initial(p, 1)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      You cannot set an initial value for a parameter. Did you mean 'set_value'?
      >>> try:
      ...   opti.set_initial(x[0]*x[1], 1)
      ... except Exception as e:
      ...   print("Use symbols or simple mappings of symbols" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2); p = opti.parameter();
      >> opti.set_initial(x, 1); opti.set_initial(10*x(1), 2);
      >> ini = opti.initial(); disp(ini{1})
      ([0.2, 1]==opti..._x_1)
      >> try, opti.set_initial(p, 1); catch e, disp(contains(e.message, 'mean ''set_value''')); end
      1
      >> try, opti.set_initial(x(1)*x(2), 1); catch e, disp(contains(e.message, 'simple')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> p = ca.parameter(opti);
      julia> ca.set_initial(opti, x, 1)
      julia> ca.set_initial(opti, 10*x[1], 2)
      julia> print(ca.initial(opti))
      MX[([0.2, 1]==opti..._x_1)]
      julia> try
                 ca.set_initial(opti, p, 1)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      You cannot set an initial value for a parameter. Did you mean 'set_value'?
      julia> try
                 ca.set_initial(opti, x[1]*x[2], 1)
             catch e
                 msg = "Use symbols or simple mappings of symbols"
                 println(occursin(msg, sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      Warm start from a previous solution:

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> y = opti.variable()
      >>> opti.minimize((y-x**2)**2)
      >>> opti.subject_to(x**2+y**2==1)
      >>> opti.subject_to(x+y>=1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> opti.set_initial(sol.value_variables())
      >>> opti.set_initial(opti.lam_g, sol.value(opti.lam_g))
      >>> opti.solve().stats()["iter_count"] < sol.stats()["iter_count"]
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); y = opti.variable();
      >> opti.minimize((y-x^2)^2); opti.subject_to(x^2+y^2==1); opti.subject_to(x+y>=1);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> opti.set_initial(sol.value_variables());
      >> opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
      >> sol2 = opti.solve(); s1 = sol.stats(); s2 = sol2.stats();
      >> s2.iter_count < s1.iter_count
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> y = ca.variable(opti);
      julia> ca.minimize(opti, (y - x^2)^2)
      julia> ca.subject_to(opti, x^2 + y^2 == 1)
      julia> ca.subject_to(opti, x + y >= 1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> opti.set_initial(ca.value_variables(sol))
      julia> ca.set_initial(opti, ca.lam_g(opti), ca.value(sol, ca.lam_g(opti)))
      julia> ca.stats(ca.solve(opti))["iter_count"] < ca.stats(sol)["iter_count"]
      true
      \enddoctest
      \endlang

      \sa set_value, initial, value_variables

      \identifier{2mt} */
  void set_initial(const MX& x, const DM& v);
  void set_initial(const std::vector<MX>& assignments);
  /// @}

  /// @{
  /** \brief Set the value of a parameter

      Every parameter that appears in the problem must have a value before solve(); the
      value stays until changed, so re-solving with new data only needs set_value() and
      solve(). \p x may be a parameter, a slice, or a linear mapping of parameters. Decision
      variables are refused with a hint to use set_initial(). The list form accepts
      assignment expressions par==value, as returned by value_parameters(). Calling this
      method invalidates value() of variables until the next solve().

      \param x Parameter (or linear mapping of parameters)
      \param v Numerical value; scalars are broadcast

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter(2)
      >>> opti.set_value(p, [1, 2])
      >>> opti.set_value(p[1], 5)
      >>> opti.value(p)
      array([1., 5.])
      >>> try:
      ...   opti.set_value(x, 1)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      You cannot set a value for a variable. Did you mean 'set_initial'?
      >>> opti.set_value([p==7])
      >>> opti.value(p)
      array([7., 7.])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter(2);
      >> opti.set_value(p, [1; 2]); opti.set_value(p(2), 5); disp(opti.value(p)')
      1 5
      >> try, opti.set_value(x, 1); catch e, disp(contains(e.message, 'mean ''set_initial''')); end
      1
      >> opti.set_value({p==7}); disp(opti.value(p)')
      7 7
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti, 2);
      julia> ca.set_value(opti, p, [1, 2])
      julia> ca.set_value(opti, p[2], 5)
      julia> ca.value(opti, p)
      2-element Vector{Float64}:
       1.0
       5.0
      julia> try
                 ca.set_value(opti, x, 1)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      You cannot set a value for a variable. Did you mean 'set_initial'?
      julia> ca.set_value(opti, [p == 7])
      julia> ca.value(opti, p)
      2-element Vector{Float64}:
       7.0
       7.0
      \enddoctest
      \endlang

      \sa set_initial, value_parameters

      \identifier{1d} */
  void set_value(const MX& x, const DM& v);
  void set_value(const std::vector<MX>& assignments);
  /// @}

  /** \brief Set the domain of a decision variable

      Marks all scalars of \p x as integer or real (default). The information is passed to
      the solver as the nlpsol/qpsol option "discrete"; a solver without integer support
      raises an error at solve().

      \param x A decision variable (pure symbol, not an expression)
      \param domain "real" or "integer"; anything else raises an error

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.set_domain(x, "integer")
      >>> opti.minimize((x-1.4)**2)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Discrete variables require a solver with integer support
      >>> try:
      ...   opti.set_domain(x, "binary")
      ... except Exception as e:
      ...   print("Unknown domain 'binary'. Known values are 'real', 'integer'." in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.set_domain(x, 'integer');
      >> opti.minimize((x-1.4)^2);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> try, opti.solve(); catch e, disp(contains(e.message, 'solver with integer support')); end
      1
      >> try, opti.set_domain(x, 'binary'); catch e, disp(contains(e.message, 'Unknown')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.set_domain(opti, x, "integer")
      julia> ca.minimize(opti, (x - 1.4)^2)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Discrete variables require a solver with integer support
      julia> try
                 ca.set_domain(opti, x, "binary")
             catch e
                 msg = "Unknown domain 'binary'. Known values are 'real', 'integer'."
                 println(occursin(msg, sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{27t} */
  void set_domain(const MX& x, const std::string& domain);

  /** \brief Set the linear scaling of a decision variable

      The solver works with (x-offset)/scale instead of x; all user-facing values (value(),
      initial guesses, duals) stay in the original units. Useful when variables differ by
      orders of magnitude.

      \param x Decision variable or linear mapping of variables
      \param scale Scalar or matrix of the shape of \p x (default 1)
      \param offset Scalar or matrix of the shape of \p x (default 0)

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> T = opti.variable()
      >>> opti.set_linear_scale(T, 20, 273.15)
      >>> opti.minimize((T-300)**2)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> sol.value(T), opti.x_linear_scale, opti.x_linear_scale_offset
      (300.0, DM(20), DM(273.15))
      >>> opti.debug.arg()["x0"]
      DM(-13.6575)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); T = opti.variable();
      >> opti.set_linear_scale(T, 20, 273.15); opti.minimize((T-300)^2);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> sol.value(T), opti.x_linear_scale, opti.x_linear_scale_offset
      300
      20
      273.15
      >> a = opti.debug.arg(); disp(a.x0)
      -13.6575
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> T = ca.variable(opti);
      julia> ca.set_linear_scale(opti, T, 20, 273.15)
      julia> ca.minimize(opti, (T - 300)^2)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, T)
      300.0
      julia> ca.x_linear_scale(opti), ca.x_linear_scale_offset(opti)
      (20, 273.15)
      \enddoctest
      \endlang

      \sa x_linear_scale, x_linear_scale_offset, scale_helper

      \identifier{2bs} */
  void set_linear_scale(const MX& x, const DM& scale, const DM& offset=0);

  /** \brief Solve the problem

      Builds (or reuses) the solver, evaluates the bounds for the current parameter values,
      runs the solver from the current initial guess and returns an OptiSol snapshot. Raises
      an error when the solver does not report success (including iteration or time limits,
      see solve_limited()); the latest iterate then remains available through
      opti.debug.value() and opti.value(). Calling solve() again without changing anything
      repeats the same solve from the same initial guess.

      \return OptiSol copy of the problem and its solution

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-1)**6)
      >>> opti.solver("ipopt", {"print_time": False},
      ...             {"print_level": 0, "max_iter": 2, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print("Solver failed" in str(e), "'Maximum_Iterations_Exceeded'" in str(e))
      True True
      >>> opti.debug.value(x), opti.return_status()
      (0.36, 'Maximum_Iterations_Exceeded')
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> sol.value(x), sol.stats()["success"]
      (1.0, True)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^6);
      >> ipopt = struct('print_level', 0, 'max_iter', 2, 'sb', 'yes');
      >> opti.solver('ipopt', struct('print_time', false), ipopt);
      >> try, opti.solve(); catch e, disp(contains(e.message, 'Maximum_Iterations_Exceeded')); end
      1
      >> opti.debug.value(x), opti.return_status()
      0.36
      Maximum_Iterations_Exceeded
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); s = sol.stats(); disp(sol.value(x)), disp(s.success)
      1.0
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 1)^6)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "max_iter" => 2, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
                 msg = sprint(showerror, e)
                 println(occursin("Solver failed", msg), " ",
                         occursin("'Maximum_Iterations_Exceeded'", msg))
             end
      true true
      julia> ca.value(ca.debug(opti), x), ca.return_status(opti)
      (0.36, "Maximum_Iterations_Exceeded")
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, x), ca.stats(sol)["success"]
      (1.0, true)
      \enddoctest
      \endlang

      \sa solve_limited, OptiSol, stats, return_status

      \identifier{2mu} */
  OptiSol solve();

  /** \brief Solve the problem, accepting a limit as outcome

      Same as solve(), but does not raise when the solver stopped because of an iteration or
      time limit (unified_return_status "SOLVER_RET_LIMITED"); other failures such as
      infeasibility still raise.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-1)**6)
      >>> opti.solver("ipopt", {"print_time": False},
      ...             {"print_level": 0, "max_iter": 2, "sb": "yes"})
      >>> sol = opti.solve_limited()
      >>> sol.stats()["unified_return_status"], sol.value(x)
      ('SOLVER_RET_LIMITED', 0.36)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^6);
      >> ipopt = struct('print_level', 0, 'max_iter', 2, 'sb', 'yes');
      >> opti.solver('ipopt', struct('print_time', false), ipopt);
      >> sol = opti.solve_limited(); s = sol.stats();
      >> s.unified_return_status, sol.value(x)
      SOLVER_RET_LIMITED
      0.36
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 1)^6)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "max_iter" => 2, "sb" => "yes"))
      julia> sol = ca.solve_limited(opti);
      julia> ca.stats(sol)["unified_return_status"], ca.value(sol, x)
      ("SOLVER_RET_LIMITED", 0.36)
      \enddoctest
      \endlang

      \identifier{1e} */
  OptiSol solve_limited();

  /// @{
  /** \brief Evaluate an expression at the current values

      Variables take the latest solver result, parameters their set_value(), dual variables
      the latest multipliers. Before the first solve(), or after set_initial()/set_value()
      changed the state, an expression involving variables raises an error unless every
      variable in it is overruled through \p values; parameter-only expressions always work.
      After a failed solve() the latest iterate is returned (the same holds for
      opti.debug.value()). A variable that does not appear in the problem raises an error.

      \param x Expression in the symbols of this instance (an MX). A DM is returned as is.
             SX expressions are not supported: a NaN matrix of the same sparsity is returned.
      \param values Optional list of assignments var==value (e.g. opti.initial()) that
             overrule the stored values for this evaluation only
      \return The numerical value

      \lang{python}
      A scalar becomes a float, a dense matrix a numpy array and a sparse matrix a
      scipy sparse matrix.

      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.minimize((x-p)**2)
      >>> opti.set_value(p, 3)
      >>> opti.value(2*p)
      6.0
      >>> try:
      ...   opti.value(x)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      This action is forbidden since you have not solved the Opti stack yet (with calling 'solve').
      >>> opti.value(x+p, [x==1]), opti.value(x**2, opti.initial())
      (4.0, 0.0)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> opti.value(x), sol.value(x), sol.value(x+p, [p==10])
      (3.0, 3.0, 13.0)
      >>> type(sol.value(vertcat(x, p))).__name__, sol.value(DM([1, 2]))
      ('ndarray', array([1., 2.]))
      \enddoctest
      \endlang

      \lang{matlab}
      The value is a double matrix, sparse when the expression is sparse.

      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.minimize((x-p)^2); opti.set_value(p, 3);
      >> opti.value(2*p)
      6
      >> try, opti.value(x); catch e, disp(contains(e.message, 'not solved the Opti stack')); end
      1
      >> opti.value(x+p, {x==1}), opti.value(x^2, opti.initial())
      4
      0
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> opti.value(x), sol.value(x), sol.value(x+p, {p==10})
      3
      3
      13
      >> v = sol.value(vertcat(x, p)); disp(class(v)), disp(size(v)), disp(sol.value(DM([1; 2]))')
      double
      2 1
      1 2
      \enddoctest
      \endlang

      \lang{julia}
      A scalar becomes a Float64, any other value a dense Matrix{Float64}.

      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.minimize(opti, (x - p)^2)
      julia> ca.set_value(opti, p, 3)
      julia> ca.value(opti, 2*p)
      6.0
      julia> try
                 ca.value(opti, x)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      This action is forbidden since you have not solved the Opti stack yet (with calling 'solve').
      julia> opti.value(x + p, [x == 1]), opti.value(x^2, ca.initial(opti))
      (4.0, 0.0)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(opti, x), ca.value(sol, x)
      (3.0, 3.0)
      julia> sol.value(x + p, [p == 10])
      13.0
      julia> ca.value(sol, [x; p]), ca.value(sol, ca.DM([1, 2]))
      ([3.0, 3.0], [1.0, 2.0])
      \enddoctest
      \endlang

      \sa initial, value_variables, OptiSol::value

      \identifier{2mv} */
  native_DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
  native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const;
  native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const;
  /// @}

  /** \brief Get the statistics of the last solve

      The dict returned by the nlpsol/qpsol plugin, passed as-is: always "success",
      "return_status" and "unified_return_status", plus plugin-specific entries such as
      "iter_count". Raises an error before the first solve(). No stability is guaranteed for
      the plugin-specific entries.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-1)**2)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> s = sol.stats()
      >>> s["success"], s["return_status"], s["unified_return_status"], s["iter_count"]
      (True, 'Solve_Succeeded', 'SOLVER_RET_SUCCESS', 1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^2);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); s = sol.stats();
      >> s.success, s.return_status, s.unified_return_status, s.iter_count
      1
      Solve_Succeeded
      SOLVER_RET_SUCCESS
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 1)^2)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> s = ca.stats(sol);
      julia> s["success"], s["return_status"], s["unified_return_status"], s["iter_count"]
      (true, "Solve_Succeeded", "SOLVER_RET_SUCCESS", 1)
      \enddoctest
      \endlang

      \identifier{1f} */
  Dict stats() const;

  /** \brief Get the return status of the last solve

      The plugin's "return_status" string (e.g. "Solve_Succeeded" or
      "Infeasible_Problem_Detected" for Ipopt), or "unknown" when no solve has happened.
      No stability is guaranteed for the strings.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize(x**2)
      >>> opti.return_status()
      'unknown'
      >>> opti.subject_to(x >= 1)
      >>> opti.subject_to(x <= -1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   print(opti.return_status())
      Infeasible_Problem_Detected
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize(x^2);
      >> opti.return_status()
      unknown
      >> opti.subject_to(x >= 1); opti.subject_to(x <= -1);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> try, opti.solve(); catch e, disp(opti.return_status()); end
      Infeasible_Problem_Detected
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, x^2)
      julia> ca.return_status(opti)
      "unknown"
      julia> ca.subject_to(opti, x >= 1)
      julia> ca.subject_to(opti, x <= -1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
                 println(ca.return_status(opti))
             end
      Infeasible_Problem_Detected
      \enddoctest
      \endlang

      \identifier{1g} */
  std::string return_status() const;

  /** \brief Get assignment expressions for the initial guesses

      One expression value==var per variable and dual variable, holding the current initial
      guess. Pass the list to value() to evaluate an expression at the initial guess, or to
      set_initial() of another Opti instance with the same symbols.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.set_initial(x, 2)
      >>> opti.minimize(x**2)
      >>> print(opti.initial())
      [MX((2==opti..._x_1))]
      >>> opti.value(x**2+1, opti.initial())
      5.0
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.set_initial(x, 2); opti.minimize(x^2);
      >> ini = opti.initial(); disp(ini{1})
      (2==opti..._x_1)
      >> opti.value(x^2+1, opti.initial())
      5
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.set_initial(opti, x, 2)
      julia> ca.minimize(opti, x^2)
      julia> print(ca.initial(opti))
      MX[(2==opti..._x_1)]
      julia> opti.value(x^2 + 1, ca.initial(opti))
      5.0
      \enddoctest
      \endlang

      \identifier{266} */
  std::vector<MX> initial() const;

  /** \brief Get assignment expressions for the latest values

      value_variables() gives value==var for every variable with its latest solver result;
      value_parameters() gives value==par for every parameter with its set value. Both feed
      set_initial()/set_value() to transfer a solution to the next solve or to another Opti.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.minimize((x-p)**2)
      >>> opti.set_value(p, 2)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> print(opti.value_variables(), opti.value_parameters())
      [MX((2==opti..._x_1))] [MX((2==opti..._p_1))]
      >>> opti.set_initial(sol.value_variables())
      >>> opti.value(x, opti.initial())
      2.0
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.minimize((x-p)^2); opti.set_value(p, 2);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> vv = opti.value_variables(); vp = opti.value_parameters(); disp(vv{1}), disp(vp{1})
      (2==opti..._x_1)
      (2==opti..._p_1)
      >> opti.set_initial(sol.value_variables()); disp(opti.value(x, opti.initial()))
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.minimize(opti, (x - p)^2)
      julia> ca.set_value(opti, p, 2)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> print(ca.value_variables(opti), " ", ca.value_parameters(opti))
      MX[(2==opti..._x_1)] MX[(2==opti..._p_1)]
      julia> opti.set_initial(ca.value_variables(sol))
      julia> opti.value(x, ca.initial(opti))
      2.0
      \enddoctest
      \endlang

      \identifier{267} */
  std::vector<MX> value_variables() const;
  std::vector<MX> value_parameters() const;

  /** \brief Scale a helper function constructed via opti.x, opti.g, ...

      Wraps a user-written derivative function so that it operates in the scaled space the
      solver sees (set_linear_scale(), minimize()/subject_to() scales). \p h must be named
      "nlp_jac_g" (inputs x, p; outputs g, jac_g_x) or "nlp_hess_l" (inputs x, p, lam_f,
      lam_g; output the Lagrangian Hessian); other names raise an error. The result can be
      handed to a solver that accepts custom derivative functions.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.set_linear_scale(x, 2, 1)
      >>> opti.minimize((x-p)**2)
      >>> opti.subject_to(x >= 0)
      >>> h = Function("nlp_jac_g", [opti.x, opti.p], [opti.g, jacobian(opti.g, opti.x)],
      ...              ["x", "p"], ["g", "jac_g_x"])
      >>> hs = opti.scale_helper(h)
      >>> hs(x=3, p=0)
      {'g': DM(7), 'jac_g_x': DM(2)}
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.set_linear_scale(x, 2, 1); opti.minimize((x-p)^2); opti.subject_to(x >= 0);
      >> J = jacobian(opti.g, opti.x);
      >> h = Function('nlp_jac_g', {opti.x, opti.p}, {opti.g, J}, {'x', 'p'}, {'g', 'jac_g_x'});
      >> hs = opti.scale_helper(h); r = hs('x', 3, 'p', 0); disp(r.g), disp(r.jac_g_x)
      7
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.set_linear_scale(opti, x, 2, 1)
      julia> ca.minimize(opti, (x - p)^2)
      julia> ca.subject_to(opti, x >= 0)
      julia> g, X = ca.g(opti), ca.x(opti);
      julia> h = ca.Function("nlp_jac_g", [X, ca.p(opti)], [g, ca.jacobian(g, X)],
                             ["x", "p"], ["g", "jac_g_x"]);
      julia> hs = ca.scale_helper(opti, h);
      julia> hs(x=3, p=0)
      Dict{String, DM} with 2 entries:
        "g" => 7
        "jac_g_x" => 2
      \enddoctest
      \endlang

      \identifier{2ci} */
  Function scale_helper(const Function& h) const;

  /** \brief Get the dual variable (Lagrange multiplier) of a constraint

      \p m must be the exact expression object that was passed to subject_to(). The result
      is symbolic; evaluate it with value(). For inequality constraints the dual is
      nonnegative (the sign of lam_g is removed); for equality constraints it equals the
      corresponding entries of lam_g. Keep a reference to the constraint expression when you
      need its dual later.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-3)**2)
      >>> con = x <= 1
      >>> opti.subject_to(con)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> sol.value(opti.dual(con)), sol.value(opti.lam_g)
      (4.0, 4.0)
      >>> opti.subject_to()
      >>> con = x >= 4
      >>> opti.subject_to(con)
      >>> sol = opti.solve()
      >>> sol.value(opti.dual(con)), sol.value(opti.lam_g)
      (2.0, -2.0)
      >>> try:
      ...   opti.dual(x >= 4)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Constraint not present in Opti stack.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-3)^2);
      >> con = x <= 1; opti.subject_to(con);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve();
      >> sol.value(opti.dual(con)), sol.value(opti.lam_g)
      4
      4
      >> opti.subject_to(); con = x >= 4; opti.subject_to(con); sol = opti.solve();
      >> sol.value(opti.dual(con)), sol.value(opti.lam_g)
      2
      -2
      >> try, opti.dual(x >= 4); catch e, disp(contains(e.message, 'Constraint not present')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 3)^2)
      julia> con = x <= 1;
      julia> ca.subject_to(opti, con)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, ca.dual(opti, con))
      4.0
      julia> ca.value(sol, ca.lam_g(opti))
      4.0
      julia> ca.subject_to(opti)
      julia> con = x >= 4;
      julia> ca.subject_to(opti, con)
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, ca.dual(opti, con))
      2.0
      julia> ca.value(sol, ca.lam_g(opti))
      -2.0
      julia> try
                 ca.dual(opti, x >= 4)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Constraint not present in Opti stack.
      \enddoctest
      \endlang

      \sa lam_g

      \identifier{1h} */
  MX dual(const MX& m) const;

  /** \brief Number of (scalarised) decision variables

      Counts the nonzeros of all variables that appear in the objective or constraints.
      In Python this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(3)
      >>> y = opti.variable()
      >>> opti.minimize(sumsqr(x))
      >>> opti.nx
      3
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(3); y = opti.variable();
      >> opti.minimize(sumsqr(x)); disp(opti.nx)
      3
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 3);
      julia> y = ca.variable(opti);
      julia> ca.minimize(opti, ca.sumsqr(x))
      julia> ca.nx(opti)
      3
      \enddoctest
      \endlang

      \identifier{268} */
  casadi_int nx() const;

  /** \brief Number of (scalarised) parameters

      Counts the nonzeros of all parameters that appear in the problem. In Python this is a
      read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter(2, 2)
      >>> opti.minimize(x*sum1(vec(p)))
      >>> opti.np
      4
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter(2, 2);
      >> opti.minimize(x*sum1(vec(p))); disp(opti.np)
      4
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti, 2, 2);
      julia> ca.minimize(opti, x*sum(vec(p)))
      julia> ca.np(opti)
      4
      \enddoctest
      \endlang

      \identifier{269} */
  casadi_int np() const;

  /** \brief Number of (scalarised) constraints

      Total number of constraint rows; a double-bounded constraint counts once. In Python
      this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.subject_to(0 <= (x <= 1))
      >>> opti.subject_to(x[0] == x[1])
      >>> opti.ng
      3
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2);
      >> opti.subject_to(0 <= x <= 1); opti.subject_to(x(1) == x(2)); disp(opti.ng)
      3
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.subject_to(opti, 0 <= (x <= 1))
      julia> ca.subject_to(opti, x[1] == x[2])
      julia> ca.ng(opti)
      3
      \enddoctest
      \endlang

      \identifier{26a} */
  casadi_int ng() const;

  /** \brief Get all (scalarised) decision variables as a symbolic column vector

      The vector the solver optimizes over, in order of variable() creation, restricted to
      the variables that appear in the problem. Use it with jacobian(), hessian() or
      to_function(). In Python this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> y = opti.variable()
      >>> opti.minimize(sumsqr(x) + y**2)
      >>> opti.x.shape
      (3, 1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> sol.value(jacobian(opti.f, opti.x))
      array([0., 0., 0.])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2); y = opti.variable();
      >> opti.minimize(sumsqr(x) + y^2); disp(size(opti.x))
      3 1
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); disp(sol.value(jacobian(opti.f, opti.x)))
      0 0 0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> y = ca.variable(opti);
      julia> ca.minimize(opti, ca.sumsqr(x) + y^2)
      julia> size(ca.x(opti))
      (3, 1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, ca.jacobian(ca.f(opti), ca.x(opti)))
      3-element Vector{Float64}:
       0.0
       0.0
       0.0
      \enddoctest
      \endlang

      \identifier{26b} */
  MX x() const;

  /** \brief Get all (scalarised) parameters as a symbolic column vector

      Restricted to the parameters that appear in the problem. In Python this is a read-only
      property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter(2)
      >>> opti.minimize(x*p[0])
      >>> opti.p.shape
      (2, 1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter(2);
      >> opti.minimize(x*p(1)); disp(size(opti.p))
      2 1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti, 2);
      julia> ca.minimize(opti, x*p[1])
      julia> size(ca.p(opti))
      (2, 1)
      \enddoctest
      \endlang

      \identifier{26c} */
  MX p() const;

  /** \brief Get all (scalarised) constraint expressions as a column vector

      The expression g passed to the solver, with lbg <= g <= ubg. In Python this is a
      read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.subject_to(x >= 1)
      >>> opti.subject_to(x**2 == 4)
      >>> opti.minimize(x)
      >>> print(opti.g)
      vertcat(opti..._x_1, sq(opti..._x_1))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable();
      >> opti.subject_to(x >= 1); opti.subject_to(x^2 == 4); opti.minimize(x);
      >> opti.g
      vertcat(opti..._x_1, sq(opti..._x_1))
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.subject_to(opti, x >= 1)
      julia> ca.subject_to(opti, x^2 == 4)
      julia> ca.minimize(opti, x)
      julia> ca.g(opti)
      MX(vertcat(opti33_x_1, sq(opti33_x_1)))
      \enddoctest
      \endlang

      \identifier{26d} */
  MX g() const;

  /** \brief Get the objective expression

      In Python this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.f
      MX(0)
      >>> opti.minimize(x**2)
      >>> print(opti.f)
      sq(opti..._x_1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); disp(opti.f)
      0
      >> opti.minimize(x^2); disp(opti.f)
      sq(opti..._x_1)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.f(opti)
      MX(0)
      julia> ca.minimize(opti, x^2)
      julia> ca.f(opti)
      MX(sq(opti34_x_1))
      \enddoctest
      \endlang

      \identifier{26e} */
  MX f() const;

  /** \brief Get all (scalarised) bounds on constraints as a column vector

      lbg() and ubg() are the lower and upper bounds matching g(); they are expressions in
      the parameters (constants otherwise), with -inf/inf for one-sided constraints. In
      Python these are read-only properties.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.subject_to(x >= 1)
      >>> opti.subject_to(x <= p)
      >>> opti.subject_to(x**2 == 4)
      >>> opti.minimize(x)
      >>> print(opti.lbg, opti.ubg)
      vertcat(1, -inf, 4) vertcat(inf, opti..._p_1, 4)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.subject_to(x >= 1); opti.subject_to(x <= p); opti.subject_to(x^2 == 4);
      >> opti.minimize(x); disp(opti.lbg), disp(opti.ubg)
      vertcat(1, -inf, 4)
      vertcat(inf, opti..._p_1, 4)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.subject_to(opti, x >= 1)
      julia> ca.subject_to(opti, x <= p)
      julia> ca.subject_to(opti, x^2 == 4)
      julia> ca.minimize(opti, x)
      julia> print(ca.lbg(opti), " ", ca.ubg(opti))
      vertcat(1, -inf, 4) vertcat(inf, opti..._p_1, 4)
      \enddoctest
      \endlang

      \identifier{26f} */
  MX lbg() const;
  MX ubg() const;

  /// Per-entry scale of x (see set_linear_scale); a Python property
  DM x_linear_scale() const;
  /// Per-entry offset of x (see set_linear_scale); a Python property
  DM x_linear_scale_offset() const;
  /// Per-row scale of g (see subject_to); a Python property
  DM g_linear_scale() const;
  /// Scale of the objective (see minimize); a Python property
  double f_linear_scale() const;

  /** \brief Get all (scalarised) dual variables as a symbolic column vector

      One multiplier per row of g(), as the solver reports them: with Lagrangian
      f + lam_g'*g, an entry is positive when the upper bound is active and negative when
      the lower bound is active, and -lam_g is the sensitivity of the optimal objective to
      the active bound. Useful for the Lagrangian Hessian, warm starting (set_initial) and
      as input or output of to_function(). In Python this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.minimize(x**2)
      >>> opti.subject_to(x >= p)
      >>> opti.set_value(p, 1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> sol.value(opti.lam_g)
      -2.0
      >>> sol.value(hessian(opti.f + dot(opti.lam_g, opti.g), opti.x)[0])
      2.0
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.minimize(x^2); opti.subject_to(x >= p); opti.set_value(p, 1);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); disp(sol.value(opti.lam_g))
      -2
      >> H = hessian(opti.f + dot(opti.lam_g, opti.g), opti.x); disp(sol.value(H))
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.minimize(opti, x^2)
      julia> ca.subject_to(opti, x >= p)
      julia> ca.set_value(opti, p, 1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> ca.value(sol, ca.lam_g(opti))
      -2.0
      julia> L = ca.f(opti) + ca.dot(ca.lam_g(opti), ca.g(opti));
      julia> ca.value(sol, ca.hessian(L, ca.x(opti))[1])
      2.0
      \enddoctest
      \endlang

      \sa dual

      \identifier{1i} */
  MX lam_g() const;

  /// @{
  /** \brief Create a CasADi Function from the Opti solver

      Packs bounds evaluation, the solver call and the evaluation of \p res into one
      Function, which is much cheaper to call repeatedly than set_value/solve/value. Inputs
      corresponding to parameters supply their values; inputs corresponding to variables or
      dual variables supply the initial guess; anything not listed keeps the value stored in
      the Opti at the time of the call. Outputs are evaluated at the solution.

      A solver failure does not raise by default: the last iterate is returned and the
      Function's stats() report success False. Pass "error_on_fail": True in the
      plugin_options of solver() to raise instead. Function options in \p opts are applied
      to the result; the entry "helper_options" is applied to the internal evaluation
      function.

      \param name Name of the resulting Function
      \param args List of parameters, decision variables and dual variables (or vertcat of
             them); each must be purely symbolic and independent, else an error is raised
      \param res List of expressions evaluated at the solution
      \param name_in Optional input names (same length as \p args)
      \param name_out Optional output names (same length as \p res)
      \param opts Function options

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> y = opti.variable()
      >>> p = opti.parameter()
      >>> opti.minimize((x-p)**2 + (y-x)**2)
      >>> opti.subject_to(x >= 0)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> F = opti.to_function("F", [p, x], [vertcat(x, y), opti.lam_g],
      ...                      ["p", "x0"], ["xy", "lam"])
      >>> F
      Function(F:(p,x0)->(xy[2],lam) MXFunction)
      >>> r = F(p=2, x0=1)
      >>> r["xy"], r["lam"]
      (DM([2, 2]), DM(0.0000))
      >>> F.stats()["success"]
      True
      >>> F = opti.to_function("F", {"p": p, "x": x}, ["p"], ["x"])
      >>> F(3)
      DM(3)
      >>> try:
      ...   opti.to_function("F", [p+1], [x])
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Argument 0 is not purely symbolic.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); y = opti.variable(); p = opti.parameter();
      >> opti.minimize((x-p)^2 + (y-x)^2); opti.subject_to(x >= 0);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> F = opti.to_function('F', {p, x}, {vertcat(x, y), opti.lam_g}, {'p', 'x0'}, {'xy', 'lam'});
      >> F
      F:(p,x0)->(xy[2],lam) MXFunction
      >> r = F('p', 2, 'x0', 1); disp(r.xy), disp(r.lam)
      [2, 2]
      0.0000
      >> s = F.stats(); disp(s.success)
      1
      >> F = opti.to_function('F', struct('p', p, 'x', x), {'p'}, {'x'}); disp(F(3))
      3
      >> try, opti.to_function('F', {p+1}, {x}); catch e, disp(contains(e.message, 'symbolic')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> y = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.minimize(opti, (x - p)^2 + (y - x)^2)
      julia> ca.subject_to(opti, x >= 0)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> F = ca.to_function(opti, "F", [p, x], [[x; y], ca.lam_g(opti)],
                                ["p", "x0"], ["xy", "lam"])
      Function(F:(p,x0)->(xy[2],lam) MXFunction)
      julia> r = F(p=2, x0=1);
      julia> r["xy"], r["lam"]
      ([2, 2], 0.0000)
      julia> ca.stats(F)["success"]
      true
      julia> F = ca.to_function(opti, "F", Dict("p" => p, "x" => x), ["p"], ["x"]);
      julia> F(3)
      DM(3)
      julia> try
                 ca.to_function(opti, "F", [p + 1], [x])
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Argument 0 is not purely symbolic.
      \enddoctest
      \endlang

      \sa nlpsol

      \identifier{1j} */
  Function to_function(const std::string& name,
      const std::vector<MX>& args, const std::vector<MX>& res,
      const Dict& opts = Dict());

  Function to_function(const std::string& name,
      const std::vector<MX>& args, const std::vector<MX>& res,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts = Dict());

  /// Same, with inputs and outputs picked by name from a dict of expressions
  Function to_function(const std::string& name,
      const std::map<std::string, MX>& dict,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts = Dict());
  /// @}

  #ifndef SWIGMATLAB
  /** \brief Construct a double inequality

      Returns the expression (lb <= expr) <= ub, which subject_to() recognizes as a
      double-bounded constraint lb <= expr <= ub when \p lb and \p ub are free of decision
      variables. Identical to writing lb <= (expr <= ub) by hand; provided because Python's
      chained comparison lb <= expr <= ub raises an error for MX. A static method; the
      instance is not needed.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.subject_to(opti.bounded(-1, x, 1))
      >>> opti.subject_to(Opti.bounded(0, x[0] + x[1], 1))
      >>> opti.ng, opti.lbg, opti.ubg
      (3, MX(vertcat(all_-1(2x1), 0)), MX(ones(3x1)))
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.subject_to(opti, ca.bounded(-1, x, 1))
      julia> ca.subject_to(opti, ca.Opti.bounded(0, x[1] + x[2], 1))
      julia> ca.ng(opti), ca.lbg(opti), ca.ubg(opti)
      (3, vertcat(all_-1(2x1), 0), ones(3x1))
      \enddoctest
      \endlang

      \identifier{1k} */
  static MX bounded(const MX& lb, const MX& expr, const MX& ub) { return (lb<=expr)<= ub; }
  #endif

  /** \brief Get a copy with advanced (debugging) functionality

      Returns an OptiAdvanced deep copy of the current state: same symbols, constraints,
      values and latest solver result, with extra methods such as value() at a non-converged
      iterate, show_infeasibilities(), x_describe()/g_describe() and casadi_solver().
      Changes to the copy do not affect the original. No API stability is guaranteed for the
      extra methods. In Python this is a read-only property: opti.debug.value(x).

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize((x-1)**6)
      >>> opti.solver("ipopt", {"print_time": False},
      ...             {"print_level": 0, "max_iter": 2, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   pass
      >>> opti.debug.value(x)
      0.36
      >>> print(opti.debug.x_describe(0))
      Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^6);
      >> ipopt = struct('print_level', 0, 'max_iter', 2, 'sb', 'yes');
      >> opti.solver('ipopt', struct('print_time', false), ipopt);
      >> try, opti.solve(); catch e, end
      >> opti.debug.value(x)
      0.36
      >> opti.debug.x_describe(1)
      Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, (x - 1)^6)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "max_iter" => 2, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
             end
      julia> ca.value(ca.debug(opti), x)
      0.36
      julia> print(ca.x_describe(ca.debug(opti), 0))
      Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \sa OptiAdvanced, advanced, copy

      \identifier{1l} */
  OptiAdvanced debug() const;

  /** \brief Get a copy with advanced functionality

      Same as debug(); the name signals use of the unstable OptiAdvanced API rather than
      debugging. In Python this is a read-only property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize(x**2)
      >>> print(opti.advanced.symvar())
      [MX(opti..._x_1)]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize(x^2);
      >> sv = opti.advanced.symvar(); disp(numel(sv)), disp(sv{1})
      1
      opti..._x_1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, x^2)
      julia> print(ca.symvar(ca.advanced(opti)))
      MX[opti..._x_1]
      \enddoctest
      \endlang

      \identifier{1m} */
  OptiAdvanced advanced() const;

  /** \brief Get a deep copy of the Opti instance

      The copy shares the symbols (so expressions in them stay valid for both) but has its own
      constraints, values, initial guesses and solver state. Updating the copy does not update
      the original.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize(x**2)
      >>> opti2 = opti.copy()
      >>> opti2.subject_to(x >= 2)
      >>> opti.ng, opti2.ng
      (0, 1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize(x^2);
      >> opti2 = opti.copy(); opti2.subject_to(x >= 2); fprintf('%d %d\n', opti.ng, opti2.ng)
      0 1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, x^2)
      julia> opti2 = opti.copy();
      julia> ca.subject_to(opti2, x >= 2)
      julia> ca.ng(opti), ca.ng(opti2)
      (0, 1)
      \enddoctest
      \endlang

      \identifier{1n} */
  Opti copy() const;

  /** \brief Add user data to a symbol or constraint

      Merges \p meta into the dict attached to a variable, parameter, dual variable (or every
      symbol of an expression) or to a constraint expression passed to subject_to(). Opti
      itself stores a "stacktrace" entry there in Python and MATLAB, which is what
      describe()/x_describe() print as "called from".

      \param m Symbol, expression of symbols, or constraint expression
      \param meta Dict entries to add (existing keys are overwritten)

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.update_user_dict(x, {"name": "speed", "unit": "m/s"})
      >>> d = opti.user_dict(x)
      >>> d["name"], d["unit"], "stacktrace" in d
      ('speed', 'm/s', True)
      >>> con = x >= 0
      >>> opti.subject_to(con)
      >>> opti.update_user_dict(con, {"name": "nonneg"})
      >>> opti.user_dict(con)["name"]
      'nonneg'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable();
      >> opti.update_user_dict(x, struct('name', 'speed', 'unit', 'm/s'));
      >> d = opti.user_dict(x); fprintf('%s %s %d\n', d.name, d.unit, isfield(d, 'stacktrace'))
      speed m/s 1
      >> con = x >= 0; opti.subject_to(con); opti.update_user_dict(con, struct('name', 'nonneg'));
      >> d = opti.user_dict(con); disp(d.name)
      nonneg
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.update_user_dict(opti, x, Dict("name" => "speed", "unit" => "m/s"))
      julia> d = ca.user_dict(opti, x);
      julia> d["name"], d["unit"], haskey(d, "stacktrace")
      ("speed", "m/s", true)
      julia> con = x >= 0;
      julia> ca.subject_to(opti, con)
      julia> ca.update_user_dict(opti, con, Dict("name" => "nonneg"))
      julia> ca.user_dict(opti, con)["name"]
      "nonneg"
      \enddoctest
      \endlang

      \identifier{1o} */
  void update_user_dict(const MX& m, const Dict& meta);
  void update_user_dict(const std::vector<MX>& m, const Dict& meta);
  /** \brief Get the user data of a symbol or constraint

      \param m A single symbol or a constraint expression; other expressions raise an error

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> sorted(opti.user_dict(x).keys())
      ['stacktrace']
      >>> try:
      ...   opti.user_dict(x+1)
      ... except Exception as e:
      ...   print(str(e).splitlines()[-1])
      Symbol expected, got expression.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); disp(fieldnames(opti.user_dict(x)))
      {'stacktrace'}
      >> try, opti.user_dict(x+1); catch e, disp(contains(e.message, 'Symbol expected')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> sort(collect(keys(ca.user_dict(opti, x))))
      1-element Vector{String}:
       "stacktrace"
      julia> try
                 ca.user_dict(opti, x + 1)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Symbol expected, got expression.
      \enddoctest
      \endlang

      \identifier{2mw} */
  Dict user_dict(const MX& m) const;

  /// Readable name of the class
  std::string type_name() const { return "Opti"; }

  /** \brief Print a summary of the instance

      Shows the instance number, the number of variables, parameters and constraints that
      appear in the problem, whether the solver is up to date and the last return status.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.minimize(sumsqr(x))
      >>> opti.subject_to(x[0] >= 1)
      >>> print(opti)
      Opti {
        instance #...
        #variables: 1 (nx = 2)
        #parameters: 0 (np = 0)
        #constraints: 1 (ng = 1)
        CasADi solver needs updating.
      }
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2);
      >> opti.minimize(sumsqr(x)); opti.subject_to(x(1) >= 1); disp(opti)
      Opti {
        instance #...
        #variables: 1 (nx = 2)
        #parameters: 0 (np = 0)
        #constraints: 1 (ng = 1)
        CasADi solver needs updating.
      }
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.minimize(opti, ca.sumsqr(x))
      julia> ca.subject_to(opti, x[1] >= 1)
      julia> opti
      Opti(Opti {
        instance #45
        #variables: 1 (nx = 2)
        #parameters: 0 (np = 0)
        #constraints: 1 (ng = 1)
        CasADi solver needs updating.
      })
      \enddoctest
      \endlang

      \identifier{2mx} */
  void disp(std::ostream& stream, bool more=false) const;

  /// Get string representation
  std::string get_str(bool more=false) const;

  ///@{
  /** \brief Helper methods for callback()

      Do not use directly; use opti.callback(fh) in Python/MATLAB, which wraps a callable in
      an OptiCallback and registers it here. The no-argument form removes the callback.

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

  /// How subject_to() classified a constraint (MetaCon::type)
  enum ConstraintType {
    OPTI_GENERIC_EQUALITY,  // g1(x,p) == g2(x,p)
    OPTI_GENERIC_INEQUALITY,  // g1(x,p) <= g2(x,p)
    OPTI_EQUALITY, // g(x,p) == bound(p)
    OPTI_INEQUALITY,  // g(x,p) <= bound(p)
    OPTI_DOUBLE_INEQUALITY,  // lb(p) <= g(x,p) <= ub(p)
    OPTI_PSD, // A(x,p) >= b(p)
    OPTI_UNKNOWN};
  /// Kind of an Opti symbol (MetaVar::type)
  enum VariableType {
    OPTI_VAR, // variable
    OPTI_PAR,  // parameter
    OPTI_DUAL_G // dual
  };
  /// Domain of a decision variable (see Opti::set_domain)
  enum DomainType {
    OPTI_DOMAIN_REAL,
    OPTI_DOMAIN_INTEGER
  };


  /// Range [start, stop) of a symbol or constraint in the flattened x or g vector
  struct IndexAbstraction {
    IndexAbstraction() : start(0), stop(0) {}
    casadi_int start;
    casadi_int stop;
  };
  /// Metadata of a constraint, as returned by OptiAdvanced::canon_expr/get_meta_con
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
  /// Metadata of a symbol, as returned by OptiAdvanced::get_meta
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


/** \brief Base class of iteration callbacks for Opti

    Subclass and override call(i), then register with Opti::callback_class. In Python and
    MATLAB, Opti::callback(fh) does this for a plain callable.

    \identifier{2my} */
class OptiCallback {
public:
  OptiCallback() {
  }
  OptiCallback(const OptiCallback& obj) {
    casadi_error("Callback objects cannot be copied");
  }
  /// Called after every solver iteration with the iteration number
  virtual void call(casadi_int i) {
    uout() << "This is a simple callback at iteration" << i << std::endl;
  }
  virtual ~OptiCallback() {}
};

/** \brief Opti with debugging and low-level access

    Obtained through Opti::debug or Opti::advanced; it is a deep copy of the Opti instance
    with the methods below added. No API stability is guaranteed for these methods. Most of
    them bake the problem (flatten variables and constraints) on the fly, so they also work
    on a copy that has not been solved.

    \identifier{2mz} */
class CASADI_EXPORT OptiAdvanced : public Opti {
  friend class InternalOptiCallback;
public:

  OptiAdvanced(const Opti& x);

  /** \brief Destructor

      \identifier{1t} */
  ~OptiAdvanced() {}


  /** \brief Get the underlying nlpsol/qpsol Function

      The solver Function constructed at the last solve() (null before that). Its inputs
      follow the nlpsol convention (x0, p, lbx, ubx, lbg, ubg, lam_x0, lam_g0), and its
      get_function() gives access to the derivative functions. In Python this is a read-only
      property.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.minimize(sumsqr(x))
      >>> opti.subject_to(x[0] >= 1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> S = opti.debug.casadi_solver
      >>> S.name(), S.class_name(), S.name_in()
      ('solver', 'IpoptInterface', ['x0', 'p', 'lbx', 'ubx', 'lbg', 'ubg', 'lam_x0', 'lam_g0'])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2);
      >> opti.minimize(sumsqr(x)); opti.subject_to(x(1) >= 1);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); S = opti.debug.casadi_solver();
      >> S.name(), S.class_name(), strjoin(cellstr(S.name_in()), ' ')
      solver
      IpoptInterface
      x0 p lbx ubx lbg ubg lam_x0 lam_g0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.minimize(opti, ca.sumsqr(x))
      julia> ca.subject_to(opti, x[1] >= 1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = ca.solve(opti);
      julia> S = opti.debug.casadi_solver;
      julia> ca.name(S), ca.class_name(S)
      ("solver", "IpoptInterface")
      julia> print(ca.name_in(S))
      ["x0", "p", "lbx", "ubx", "lbg", "ubg", "lam_x0", "lam_g0"]
      \enddoctest
      \endlang

      \identifier{2n0} */
  Function casadi_solver() const;

  /** \brief Check whether an expression depends on parameters only

      True when no decision variable of this instance appears in \p expr (constants count
      as parametric).

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.debug.is_parametric(2*p), opti.debug.is_parametric(x*p)
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> fprintf('%d %d\n', opti.debug.is_parametric(2*p), opti.debug.is_parametric(x*p))
      1 0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> d = ca.debug(opti);
      julia> ca.is_parametric(d, 2*p), ca.is_parametric(d, x*p)
      (true, false)
      \enddoctest
      \endlang

      \identifier{2n1} */
  bool is_parametric(const MX& expr) const;

  /// @{
  /** \brief Get the symbols of this instance

      Without argument: every symbol created so far (variables, parameters and the dual
      variables of constraints), whether or not it appears in the problem. With \p expr: the
      symbols of this instance that \p expr depends on, optionally filtered by \p type.
      The result is ordered by creation.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> y = opti.variable()
      >>> print(opti.debug.symvar())
      [MX(opti..._x_1), MX(opti..._p_1), MX(opti..._x_2)]
      >>> print(opti.debug.symvar(y*p + 1), opti.debug.symvar(y*p, OPTI_PAR))
      [MX(opti..._p_1), MX(opti..._x_2)] [MX(opti..._p_1)]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter(); y = opti.variable();
      >> sv = opti.debug.symvar(); disp(numel(sv)), disp(sv{2})
      3
      opti..._p_1
      >> sv = opti.debug.symvar(y*p + 1); disp(sv{1}), disp(sv{2})
      opti..._p_1
      opti..._x_2
      >> sv = opti.debug.symvar(y*p, OPTI_PAR); disp(numel(sv)), disp(sv{1})
      1
      opti..._p_1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> y = ca.variable(opti);
      julia> d = ca.debug(opti);
      julia> print(ca.symvar(d))
      MX[opti..._x_1, opti..._p_1, opti..._x_2]
      julia> print(ca.symvar(d, y*p + 1), " ", ca.symvar(d, y*p, ca.OPTI_PAR))
      MX[opti..._p_1, opti..._x_2] MX[opti..._p_1]
      \enddoctest
      \endlang

      \identifier{1u} */
  std::vector<MX> symvar() const;
  std::vector<MX> symvar(const MX& expr) const;
  std::vector<MX> symvar(const MX& expr, VariableType type) const;
  /// @}

  /** \brief Classify a constraint expression (for internal use only)

      Returns the MetaCon that subject_to() would store: the constraint type, canonical
      expression and bounds.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> c = opti.debug.canon_expr(0 <= (x <= 1))
      >>> print(c.type == OPTI_DOUBLE_INEQUALITY, c.lb, c.ub, c.canon)
      True 0 1 opti..._x_1
      >>> opti.debug.canon_expr(x + 1).type == OPTI_UNKNOWN
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable();
      >> c = opti.debug.canon_expr(0 <= x <= 1);
      >> c.type == OPTI_DOUBLE_INEQUALITY, c.lb, c.ub, c.canon
      1
      0
      1
      opti..._x_1
      >> c = opti.debug.canon_expr(x + 1); disp(c.type == OPTI_UNKNOWN)
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = opti.variable();
      julia> c = ca.canon_expr(opti.debug, 0 <= (x <= 1));
      julia> println(c.type == ca.OPTI_DOUBLE_INEQUALITY, " ", c.lb, " ", c.ub, " ", c.canon)
      true 0 1 opti..._x_1
      julia> ca.canon_expr(opti.debug, x + 1).type == ca.OPTI_UNKNOWN
      true
      \enddoctest
      \endlang

      \identifier{2n2} */
  MetaCon canon_expr(const MX& expr) const;

  /// Get meta-data of symbol (for internal use only)
  MetaVar get_meta(const MX& m) const;

  /// Get meta-data of symbol (for internal use only)
  MetaCon get_meta_con(const MX& m) const;

  /// Set meta-data of an expression
  void set_meta(const MX& m, const MetaVar& meta);

  /// Set meta-data of an expression
  void set_meta_con(const MX& m, const MetaCon& meta);

  /// Raise an error unless \p m is a symbol of this instance that appears in the problem
  void assert_active_symbol(const MX& m) const;

  /** \brief Get the symbols of a given type that appear in the problem

      Empty before the problem has been baked (solve, or any method that bakes).

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> y = opti.variable()
      >>> opti.minimize(x**2)
      >>> d = opti.debug
      >>> d.active_symvar(OPTI_VAR)
      []
      >>> d.bake()
      >>> print(d.active_symvar(OPTI_VAR), d.active_values(OPTI_VAR))
      [MX(opti..._x_1)] [DM(0)]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); y = opti.variable(); opti.minimize(x^2);
      >> d = opti.debug; disp(isempty(d.active_symvar(OPTI_VAR)))
      1
      >> d.bake(); a = d.active_symvar(OPTI_VAR); v = d.active_values(OPTI_VAR);
      >> a{1}, v{1}
      opti..._x_1
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> y = ca.variable(opti);
      julia> ca.minimize(opti, x^2)
      julia> d = ca.debug(opti);
      julia> ca.active_symvar(d, ca.OPTI_VAR)
      MX[]
      julia> ca.bake(d)
      julia> print(ca.active_symvar(d, ca.OPTI_VAR), " ", ca.active_values(d, ca.OPTI_VAR))
      MX[opti..._x_1] DM[0]
      \enddoctest
      \endlang

      \identifier{2n3} */
  std::vector<MX> active_symvar(VariableType type) const;
  /// Initial guesses (or parameter values) of the symbols returned by active_symvar
  std::vector<DM> active_values(VariableType type) const;

  /** \brief Get the variable that owns entry i of the flattened x vector

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> y = opti.variable()
      >>> opti.minimize(sumsqr(x) + y**2)
      >>> print(opti.debug.x_lookup(1), opti.debug.x_lookup(2))
      opti..._x_1 opti..._x_2
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2); y = opti.variable();
      >> opti.minimize(sumsqr(x) + y^2);
      >> opti.debug.x_lookup(2), opti.debug.x_lookup(3)
      opti..._x_1
      opti..._x_2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> y = ca.variable(opti);
      julia> ca.minimize(opti, ca.sumsqr(x) + y^2)
      julia> d = ca.debug(opti);
      julia> print(ca.x_lookup(d, 1), " ", ca.x_lookup(d, 2))
      opti..._x_1 opti..._x_2
      \enddoctest
      \endlang

      \identifier{2n4} */
  MX x_lookup(casadi_index i) const;
  /** \brief Get the constraint expression that owns row i of g

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(2)
      >>> opti.subject_to(x >= 0)
      >>> opti.subject_to(x[0] + x[1] == 1)
      >>> print(opti.debug.g_lookup(2))
      (1==(opti..._x_1[0]+opti..._x_1[1]))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(2);
      >> opti.subject_to(x >= 0); opti.subject_to(x(1) + x(2) == 1);
      >> opti.debug.g_lookup(3)
      (1==(opti..._x_1[0]+opti..._x_1[1]))
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 2);
      julia> ca.subject_to(opti, x >= 0)
      julia> ca.subject_to(opti, x[1] + x[2] == 1)
      julia> ca.g_lookup(ca.debug(opti), 2)
      MX((1==(opti52_x_1[0]+opti52_x_1[1])))
      \enddoctest
      \endlang

      \identifier{2n5} */
  MX g_lookup(casadi_index i) const;

  /// Row of the reduced g for row i of g when "detect_simple_bounds" is active (-1 if moved)
  casadi_index g_index_reduce_g(casadi_index i) const;
  /// Index into x when row i of g was reduced to a simple bound (-1 otherwise)
  casadi_index g_index_reduce_x(casadi_index i) const;
  /// Row of g for row i of the reduced g
  casadi_index g_index_unreduce_g(casadi_index i) const;

  /** \brief Describe entry i of the flattened x vector

      Human-readable description of the variable owning x[i] (name, shape, and where it was
      created in Python/MATLAB), plus the position inside the variable. Use it to trace
      NaN/Inf reports of a solver back to the source.

      \param opts "max_stacktrace_depth" (default 1)

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(3)
      >>> opti.minimize(sumsqr(x))
      >>> print(opti.debug.x_describe(2))
      Opti decision variable 'opti..._x_1' of shape 3x1, called from ...
      At nonzero 2.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(3); opti.minimize(sumsqr(x));
      >> opti.debug.x_describe(3)
      Opti decision variable 'opti..._x_1' of shape 3x1, called from ...
      At nonzero 3.
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 3);
      julia> ca.minimize(opti, ca.sumsqr(x))
      julia> print(ca.x_describe(ca.debug(opti), 2))
      Opti decision variable 'opti..._x_1' of shape 3x1, called from ...
      At nonzero 2.
      \enddoctest
      \endlang

      \identifier{2n6} */
  std::string x_describe(casadi_index i, const Dict& opts=Dict()) const;
  /** \brief Describe row i of the constraint vector g

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable(3)
      >>> opti.subject_to(x >= 0)
      >>> print(opti.debug.g_describe(1))
      Opti constraint of shape 3x1, called from ...
      At nonzero 1, part 0.
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(3); opti.subject_to(x >= 0);
      >> opti.debug.g_describe(2)
      Opti constraint of shape 3x1, called from ...
      At nonzero 2, part 0.
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti, 3);
      julia> ca.subject_to(opti, x >= 0)
      julia> print(ca.g_describe(ca.debug(opti), 1))
      Opti constraint of shape 3x1, called from ...
      At nonzero 1, part 0.
      \enddoctest
      \endlang

      \identifier{2n7} */
  std::string g_describe(casadi_index i, const Dict& opts=Dict()) const;
  /** \brief Describe an expression in terms of this instance

      For a symbol: its kind, name, shape and origin, or a note that it belongs to another
      Opti instance or was declared outside of Opti. For a constraint: its shape and origin.
      For another expression: the symbols it depends on.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.debug.describe(MX.sym("w"))
      "MX symbol 'w' of shape 1x1, declared outside of Opti."
      >>> print(opti.debug.describe(x*p))
      General expression, dependent on 2 symbols:
        Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
        Opti parameter 'opti..._p_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.debug.describe(MX.sym('w'))
      MX symbol 'w' of shape 1x1, declared outside of Opti.
      >> opti.debug.describe(x*p)
      General expression, dependent on 2 symbols:
        Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
        Opti parameter 'opti..._p_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> p = ca.parameter(opti);
      julia> ca.describe(ca.debug(opti), ca.MX.sym("w"))
      "MX symbol 'w' of shape 1x1, declared outside of Opti."
      julia> print(ca.describe(ca.debug(opti), x*p))
      General expression, dependent on 2 symbols:
        Opti decision variable 'opti..._x_1' of shape 1x1, called from ...
        Opti parameter 'opti..._p_1' of shape 1x1, called from ...
      \enddoctest
      \endlang

      \identifier{2n8} */
  std::string describe(const MX& x, casadi_index indent=0, const Dict& opts=Dict()) const;

  /** \brief Print the constraints violated at the latest iterate

      For every row of g whose violation max(g-ubg, lbg-g) is at least \p tol, prints the
      bounds, the value, the violation and the description of the constraint (also in scaled
      units when scaling is active). Requires a solve; typically called after an infeasible
      solve through opti.debug.

      \param tol Violation threshold (default 0, which also lists rounding-level violations)
      \param opts Passed to g_describe

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.minimize(x**2)
      >>> opti.subject_to(x >= 1)
      >>> opti.subject_to(x <= -1)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> try:
      ...   opti.solve()
      ... except Exception as e:
      ...   pass
      >>> opti.debug.show_infeasibilities(1e-6)
      Violated constraints (tol 1e-06), in order of declaration:
      ------- i = 0/2 ------
      1 <= ... <= inf (viol ...)
      Opti constraint of shape 1x1, called from ...
      ------- i = 1/2 ------
      -inf <= ... <= -1 (viol ...)
      Opti constraint of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); opti.minimize(x^2);
      >> opti.subject_to(x >= 1); opti.subject_to(x <= -1);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> try, opti.solve(); catch e, end
      >> opti.debug.show_infeasibilities(1e-6)
      Violated constraints (tol 1e-06), in order of declaration:
      ------- i = 1/2 ------
      1 <= ... <= inf (viol ...)
      Opti constraint of shape 1x1, called from ...
      ------- i = 2/2 ------
      -inf <= ... <= -1 (viol ...)
      Opti constraint of shape 1x1, called from ...
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.minimize(opti, x^2)
      julia> ca.subject_to(opti, x >= 1)
      julia> ca.subject_to(opti, x <= -1)
      julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                        Dict("print_level" => 0, "sb" => "yes"))
      julia> try
                 ca.solve(opti)
             catch e
             end
      julia> ca.show_infeasibilities(ca.debug(opti), 1e-6)
      Violated constraints (tol 1e-06), in order of declaration:
      ------- i = 0/2 ------
      1 <= ... <= inf (viol ...)
      Opti constraint of shape 1x1, called from ...
      ------- i = 1/2 ------
      -inf <= ... <= -1 (viol ...)
      Opti constraint of shape 1x1, called from ...
      \enddoctest
      \endlang

      \identifier{2n9} */
  void show_infeasibilities(double tol=0, const Dict& opts=Dict()) const;

  /** \brief Prepare the solver arguments (first half of solve)

      Checks the constraints and parameter values and fills arg() with x0, p, lam_g0, lbg
      and ubg. solve_actual() then runs the solver on such a dict, and res() stores its
      result so that value() works. Together they allow custom solver invocations.

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> p = opti.parameter()
      >>> opti.minimize((x-p)**2)
      >>> opti.subject_to(x <= 1)
      >>> opti.set_value(p, 3)
      >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
      >>> sol = opti.solve()
      >>> d = opti.debug
      >>> d.solve_prepare()
      >>> d.arg()
      {'lam_g0': DM(0), 'lbg': DM(-inf), 'p': DM(3), 'ubg': DM(1), 'x0': DM(0)}
      >>> r = d.solve_actual(d.arg())
      >>> d.res(r)
      >>> d.value(x), sorted(d.res().keys())
      (1.0, ['f', 'g', 'lam_g', 'lam_p', 'lam_x', 'x'])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable(); p = opti.parameter();
      >> opti.minimize((x-p)^2); opti.subject_to(x <= 1); opti.set_value(p, 3);
      >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
      >> sol = opti.solve(); d = opti.debug; d.solve_prepare();
      >> a = d.arg(); disp(fieldnames(a)'), disp(a.lbg), disp(a.ubg), disp(a.x0)
      {'lam_g0'} {'lbg'} {'p'} {'ubg'} {'x0'}
      -inf
      1
      0
      >> r = d.solve_actual(d.arg()); d.res(r);
      >> d.value(x), sort(fieldnames(d.res()))'
      1
      {'f'} {'g'} {'lam_g'} {'lam_p'} {'lam_x'} {'x'}
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = opti.variable();
      julia> p = opti.parameter();
      julia> opti.minimize((x - p)^2)
      julia> opti.subject_to(x <= 1)
      julia> opti.set_value(p, 3)
      julia> opti.solver("ipopt", Dict("print_time" => false), Dict("print_level" => 0, "sb" => "yes"))
      julia> sol = opti.solve();
      julia> d = opti.debug;
      julia> ca.solve_prepare(d)
      julia> sort(collect(ca.arg(d)), by=first)
      5-element Vector{Pair{String, DM}}:
       "lam_g0" => 0
          "lbg" => -inf
            "p" => 3
          "ubg" => 1
           "x0" => 0
      julia> r = ca.solve_actual(d, ca.arg(d));
      julia> ca.res(d, r)
      julia> ca.value(d, x), sort(collect(keys(ca.res(d))))
      (1.0, ["f", "g", "lam_g", "lam_p", "lam_x", "x"])
      \enddoctest
      \endlang

      \identifier{2na} */
  void solve_prepare();
  /// Run the solver on an argument dict (see solve_prepare)
  DMDict solve_actual(const DMDict& args);

  /// Solver arguments prepared by solve_prepare (x0, p, lam_g0, lbg, ubg)
  DMDict arg() const;
  /// Store a solver result dict as the latest values (marks the instance solved)
  void res(const DMDict& res);
  /// Result dict of the latest solve (x, f, g, lam_x, lam_g, lam_p)
  DMDict res() const;
  /** \brief Get the constraint expressions as passed to subject_to

      \lang{python}
      \doctest
      >>> opti = Opti()
      >>> x = opti.variable()
      >>> opti.subject_to(x >= 0)
      >>> opti.subject_to(1 <= (x**2 <= 4))
      >>> print(opti.debug.constraints())
      [MX((0<=opti..._x_1)), MX((1<=(sq(opti..._x_1)<=4)))]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opti = Opti(); x = opti.variable();
      >> opti.subject_to(x >= 0); opti.subject_to(1 <= x^2 <= 4);
      >> c = opti.debug.constraints(); disp(c{1}), disp(c{2})
      (0<=opti..._x_1)
      ((1<=sq(opti..._x_1))<=4)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opti = ca.Opti();
      julia> x = ca.variable(opti);
      julia> ca.subject_to(opti, x >= 0)
      julia> ca.subject_to(opti, 1 <= (x^2 <= 4))
      julia> print(ca.constraints(ca.debug(opti)))
      MX[(0<=opti..._x_1), (1<=(sq(opti..._x_1)<=4))]
      \enddoctest
      \endlang

      \identifier{2nb} */
  std::vector<MX> constraints() const;
  /// Objective expression; same as f()
  MX objective() const;

  /// Copy of the instance with the problem baked (see bake)
  OptiAdvanced baked_copy() const;

  /// Raise an error if an objective or constraints are present
  void assert_empty() const;


  /** \brief Fix the structure of the optimization problem

      Flattens the active variables, parameters and constraints into the x, p, g, lbg and ubg
      vectors and builds the bounds function. Called automatically by solve() and the
      accessors; raises an error when neither an objective nor a constraint is set.

      \identifier{2nc} */
  void bake();

  bool problem_dirty_;
  /// Flag that the problem structure changed since the last bake
  void mark_problem_dirty(bool flag=true);
  /// True when the problem must be baked again before solving
  bool problem_dirty() const;

  bool solver_dirty_;
  /// Flag that the solver must be reconstructed (solver() was called)
  void mark_solver_dirty(bool flag=true);
  /// True when the solver must be reconstructed
  bool solver_dirty() const;

  bool solved_;
  /// Set or clear the solved state
  void mark_solved(bool flag=true);
  /// True when a solver result is available for value()
  bool solved() const;

  /// Raise an error unless solved()
  void assert_solved() const;
  /// Raise an error when the problem is dirty
  void assert_baked() const;

  /// Sequence number of this Opti instance in the process; prefixes the symbol names
  casadi_int instance_number() const;

protected:
  OptiAdvanced() {}
};

/** \brief Snapshot of an Opti problem and its solution

    Returned by Opti::solve and Opti::solve_limited. It is a deep copy made at that moment:
    later set_value/set_initial/subject_to calls on the Opti do not change what value()
    returns, so several OptiSol objects can be kept and compared. The API is guaranteed to
    be stable.

    \lang{python}
    \doctest
    >>> opti = Opti()
    >>> x = opti.variable()
    >>> p = opti.parameter()
    >>> opti.minimize((x-p)**2)
    >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
    >>> opti.set_value(p, 1)
    >>> sol1 = opti.solve()
    >>> opti.set_value(p, 2)
    >>> sol2 = opti.solve()
    >>> sol1.value(x), sol2.value(x)
    (1.0, 2.0)
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> opti = Opti(); x = opti.variable(); p = opti.parameter();
    >> opti.minimize((x-p)^2);
    >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
    >> opti.set_value(p, 1); sol1 = opti.solve();
    >> opti.set_value(p, 2); sol2 = opti.solve();
    >> sol1.value(x), sol2.value(x)
    1
    2
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> opti = ca.Opti();
    julia> x = ca.variable(opti);
    julia> p = ca.parameter(opti);
    julia> ca.minimize(opti, (x - p)^2)
    julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                      Dict("print_level" => 0, "sb" => "yes"))
    julia> ca.set_value(opti, p, 1)
    julia> sol1 = ca.solve(opti);
    julia> ca.set_value(opti, p, 2)
    julia> sol2 = ca.solve(opti);
    julia> ca.value(sol1, x), ca.value(sol2, x)
    (1.0, 2.0)
    \enddoctest
    \endlang

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
    /** \brief Evaluate an expression at the solution

        Same as Opti::value on the snapshot: variables take their optimal values, parameters
        the values used for the solve, dual variables the multipliers. \p values overrules
        individual symbols for this evaluation only, e.g. opti.initial() to evaluate at the
        initial guess.

        \param x Expression (MX) in the symbols of the Opti; a DM is returned as is
        \param values Optional list of assignments var==value

        \lang{python}
        A scalar becomes a float, a dense matrix a numpy array and a sparse matrix a
        scipy sparse matrix.

        \doctest
        >>> opti = Opti()
        >>> x = opti.variable()
        >>> y = opti.variable()
        >>> opti.minimize((y-x**2)**2)
        >>> opti.subject_to(x**2+y**2==1)
        >>> opti.subject_to(x+y>=1)
        >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
        >>> sol = opti.solve()
        >>> sol.value(x), sol.value(vertcat(x, y)).shape
        (0.7862, (2,))
        >>> sol.value(opti.f, [y==2])
        1.9098
        >>> sol.value(x**2+y, opti.initial())
        0.0
        >>> sol.value(jacobian(opti.g, opti.x)).shape
        (2, 2)
        \enddoctest
        \endlang

        \lang{matlab}
        The value is a double matrix, sparse when the expression is sparse.

        \doctest
        >> opti = Opti(); x = opti.variable(); y = opti.variable();
        >> opti.minimize((y-x^2)^2); opti.subject_to(x^2+y^2==1); opti.subject_to(x+y>=1);
        >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
        >> sol = opti.solve();
        >> sol.value(x), size(sol.value(vertcat(x, y)))
        0.7862
        2 1
        >> sol.value(opti.f, {y==2})
        1.9098
        >> sol.value(x^2+y, opti.initial())
        0
        >> size(sol.value(jacobian(opti.g, opti.x)))
        2 2
        \enddoctest
        \endlang

        \lang{julia}
        A scalar becomes a Float64, any other value a dense Matrix{Float64}.

        \doctest
        julia> opti = ca.Opti();
        julia> x = ca.variable(opti);
        julia> y = ca.variable(opti);
        julia> ca.minimize(opti, (y - x^2)^2)
        julia> ca.subject_to(opti, x^2 + y^2 == 1)
        julia> ca.subject_to(opti, x + y >= 1)
        julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                          Dict("print_level" => 0, "sb" => "yes"))
        julia> sol = ca.solve(opti);
        julia> ca.value(sol, x), size(ca.value(sol, [x; y]))
        (0.7862, (2,))
        julia> sol.value(ca.f(opti), [y == 2])
        1.9098
        julia> sol.value(x^2 + y, ca.initial(opti))
        0.0
        julia> size(ca.value(sol, ca.jacobian(ca.g(opti), ca.x(opti))))
        (2, 2)
        \enddoctest
        \endlang

        \identifier{2nd} */
    native_DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const;
    /// @}

    /** \brief Get assignment expressions for the optimal solution

        value_variables() gives value==var for every variable, value_parameters() value==par
        for every parameter. Pass them to Opti::set_initial / Opti::set_value to warm start
        the next solve.

        \lang{python}
        \doctest
        >>> opti = Opti()
        >>> x = opti.variable()
        >>> p = opti.parameter()
        >>> opti.minimize((x-p)**2)
        >>> opti.set_value(p, 2)
        >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
        >>> sol = opti.solve()
        >>> print(sol.value_variables(), sol.value_parameters())
        [MX((2==opti..._x_1))] [MX((2==opti..._p_1))]
        >>> opti.set_initial(sol.value_variables())
        >>> opti.solve().stats()["iter_count"]
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> opti = Opti(); x = opti.variable(); p = opti.parameter();
        >> opti.minimize((x-p)^2); opti.set_value(p, 2);
        >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
        >> sol = opti.solve();
        >> vv = sol.value_variables(); vp = sol.value_parameters(); disp(vv{1}), disp(vp{1})
        (2==opti..._x_1)
        (2==opti..._p_1)
        >> opti.set_initial(sol.value_variables()); sol2 = opti.solve();
        >> s = sol2.stats(); disp(s.iter_count)
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> opti = ca.Opti();
        julia> x = ca.variable(opti);
        julia> p = ca.parameter(opti);
        julia> ca.minimize(opti, (x - p)^2)
        julia> ca.set_value(opti, p, 2)
        julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                          Dict("print_level" => 0, "sb" => "yes"))
        julia> sol = ca.solve(opti);
        julia> print(ca.value_variables(sol), " ", ca.value_parameters(sol))
        MX[(2==opti..._x_1)] MX[(2==opti..._p_1)]
        julia> opti.set_initial(ca.value_variables(sol))
        julia> ca.stats(ca.solve(opti))["iter_count"]
        0
        \enddoctest
        \endlang

        \identifier{2ne} */
    std::vector<MX> value_variables() const;
    std::vector<MX> value_parameters() const;

    /** \brief Get the solver statistics of this solve

        The dict returned by the nlpsol/qpsol plugin, passed as-is ("success",
        "return_status", "iter_count", ...). No stability is guaranteed for the entries.

        \lang{python}
        \doctest
        >>> opti = Opti()
        >>> x = opti.variable()
        >>> opti.minimize((x-1)**2)
        >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
        >>> sol = opti.solve()
        >>> sol.stats()["success"], sol.stats()["return_status"]
        (True, 'Solve_Succeeded')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^2);
        >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
        >> sol = opti.solve(); s = sol.stats(); disp(s.success), disp(s.return_status)
        1
        Solve_Succeeded
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> opti = ca.Opti();
        julia> x = ca.variable(opti);
        julia> ca.minimize(opti, (x - 1)^2)
        julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                          Dict("print_level" => 0, "sb" => "yes"))
        julia> sol = ca.solve(opti);
        julia> ca.stats(sol)["success"], ca.stats(sol)["return_status"]
        (true, "Solve_Succeeded")
        \enddoctest
        \endlang

        \identifier{1w} */
    Dict stats() const;

    /** \brief Get the Opti copy this solution belongs to

        An Opti holding the state at the time of the solve; opti.debug on it gives the
        OptiAdvanced accessors. In Python this is a read-only property.

        \lang{python}
        \doctest
        >>> opti = Opti()
        >>> x = opti.variable()
        >>> opti.minimize((x-1)**2)
        >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
        >>> sol = opti.solve()
        >>> sol.opti.return_status(), sol.opti.debug.value(x)
        ('Solve_Succeeded', 1.0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> opti = Opti(); x = opti.variable(); opti.minimize((x-1)^2);
        >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
        >> sol = opti.solve(); o = sol.opti;
        >> o.return_status(), o.debug.value(x)
        Solve_Succeeded
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> opti = ca.Opti();
        julia> x = ca.variable(opti);
        julia> ca.minimize(opti, (x - 1)^2)
        julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                          Dict("print_level" => 0, "sb" => "yes"))
        julia> sol = ca.solve(opti);
        julia> ca.return_status(ca.opti(sol)), ca.value(ca.debug(ca.opti(sol)), x)
        ("Solve_Succeeded", 1.0)
        \enddoctest
        \endlang

        \identifier{2nf} */
    Opti opti() const { return optistack_; } // NOLINT(cppcoreguidelines-slicing)

  protected:
    OptiSol(const Opti& opti);
    OptiAdvanced optistack_;
};


} // namespace casadi

#endif // CASADI_OPTI_HPP
