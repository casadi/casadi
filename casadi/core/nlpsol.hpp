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



#ifndef CASADI_NLPSOL_HPP
#define CASADI_NLPSOL_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_nlpsol Title
      \par

      Create an NLP solver

      Returns a Function that solves the parametric nonlinear program (NLP)
      \verbatim

      min          F(x, p)
      x

      subject to
      LBX <=   x    <= UBX
      LBG <= G(x, p) <= UBG
      p  == P

      nx: number of decision variables
      ng: number of constraints
      np: number of parameters

      \endverbatim

      The problem is given as a dictionary with fields "x" (decision variables, a dense vector),
      "f" (scalar objective), and optionally "p" (parameters, dense vector) and "g" (constraint
      expression, dense vector). All entries must be of one symbolic type (SX or MX). A missing "f"
      is taken as 0 and a missing "g" as an empty vector; any other field name is an error, as is a
      symbol in "f" or "g" that is not part of "x" or "p". The problem can also be given as a
      Function with inputs "x", "p" and outputs "f", "g" (for instance one created with external),
      as the name of a shared library or a ".c" file (compiled on the fly with clang) containing
      such a function named "nlp", as an Importer, or as an NlpBuilder (an .nl file).

      The returned Function has inputs x0, p, lbx, ubx, lbg, ubg, lam_x0, lam_g0 (see nlpsol_in)
      and outputs x, f, g, lam_x, lam_g, lam_p (see nlpsol_out). Call it with keyword arguments
      (returns a dictionary) or positionally with all 8 inputs (returns a tuple). Omitted inputs
      take their defaults: -inf for lbx and lbg, +inf for ubx and ubg, 0 otherwise (see
      nlpsol_default_in). A scalar passed for a vector input is broadcast to every entry. An input
      with lb > ub, lb = +inf or ub = -inf raises an "Ill-posed problem" error.

      Multiplier convention: lam_x and lam_g are the Lagrange multipliers of the bound constraints
      on x and g; a multiplier is negative when the lower bound is active and positive when the
      upper bound is active, such that grad_x f + jac_x(g)' lam_g + lam_x = 0 at the solution.
      lam_p is the multiplier of p == P, i.e. minus the gradient of the optimal objective with
      respect to p. lam_x0 and lam_g0 are only used by solvers that support warm starts.

      When the solver does not converge, the result still holds the last iterate and
      stats() reports "success" False and the plugin-specific "return_status"; with option
      "error_on_fail" True an exception is raised instead (default False for nlpsol).

      Options common to all plugins (the plugin options follow in the tables below):
      - "expand": convert an MX problem to SX before creating the derivative functions
        (default false). Usually faster evaluation; fails on nodes without SX support.
      - "print_time": print a timing table after each call (default true for nlpsol).
      - "error_on_fail": raise on an unsuccessful solve (default false).
      - "calc_lam_x", "calc_lam_p", "calc_f", "calc_g": recompute lam_x (from the stationarity
        condition), lam_p, f and g in the base class from the solution returned by the plugin
        (defaults false, true, false, false). "calc_multipliers" sets calc_lam_x and calc_lam_p.
      - "bound_consistency": after the solve, clip x and g to their bounds and set each multiplier
        to the sign of its active bound, or 0 when both bounds are infinite (default false).
      - "discrete": list of nx booleans marking integer variables; requires a plugin with integer
        support, otherwise an error at construction.
      - "equality": list of ng booleans hinting which constraints have equal bounds.
      - "warn_initial_bounds": warn when x0 violates lbx or ubx (default false).
      - "iteration_callback": a Function (typically a Callback) with the nlpsol output scheme as
        inputs, called every "iteration_callback_step" iterations with the current iterate.
      - "detect_simple_bounds": move entries of g that are bounds on a single variable into the
        simple bounds lbx/ubx (default false).
      - "verbose", "monitor", "show_eval_warnings", "common_options", "specific_options": control
        diagnostics and the options of the generated derivative functions.

      Plugin options are passed in a nested dictionary under the plugin name (e.g.
      {"ipopt": {"print_level": 0}}) or dotted ("ipopt.print_level"). The call to the returned
      Function can be embedded in an MX expression and differentiated: the sensitivities of the
      primal-dual solution with respect to any input follow from the KKT conditions.

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); y = SX.sym("y"); p = SX.sym("p")
      >>> nlp = {"x": vertcat(x, y), "p": p, "f": (1 - x)**2 + p*(y - x**2)**2, "g": x + y}
      >>> opts = {"ipopt.print_level": 0, "ipopt.sb": "yes", "print_time": False}
      >>> S = nlpsol("S", "ipopt", nlp, opts)
      >>> print(S)
      S:(x0[2],p,lbx[2],ubx[2],lbg,ubg,lam_x0[2],lam_g0)->(x[2],f,g,lam_x[2],lam_g,lam_p) ...
      >>> sol = S(x0=[2.5, 3.0], p=100, lbg=-1, ubg=1)
      >>> sol["f"], sol["x"]
      (DM(0.1456), DM([0.6188, 0.3812]))
      >>> S.stats()["success"], S.stats()["return_status"]
      (True, 'Solve_Succeeded')
      >>> x_opt, f_opt, g_opt, lam_x, lam_g, lam_p = S([2.5, 3.0], 100, -inf, inf, -1, 1, 0, 0)
      >>> f_opt, lam_g
      (DM(0.1456), DM(0.3407))
      >>> S(p=100)["x"]
      DM([1, 1])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); y = SX.sym('y'); p = SX.sym('p');
      >> nlp = struct('x', vertcat(x, y), 'p', p, 'f', (1 - x)^2 + p*(y - x^2)^2, 'g', x + y);
      >> opts = struct('print_time', false, 'ipopt', struct('print_level', 0, 'sb', 'yes'));
      >> S = nlpsol('S', 'ipopt', nlp, opts); disp(S)
      S:(x0[2],p,lbx[2],ubx[2],lbg,ubg,lam_x0[2],lam_g0)->(x[2],f,g,lam_x[2],lam_g,lam_p) ...
      >> sol = S('x0', [2.5; 3.0], 'p', 100, 'lbg', -1, 'ubg', 1);
      >> sol.f, sol.x
      0.1456
      [0.6188, 0.3812]
      >> st = S.stats(); disp(st.success), disp(st.return_status)
      1
      Solve_Succeeded
      >> [x_opt, f_opt, g_opt, lam_x, lam_g, lam_p] = S([2.5; 3.0], 100, -inf, inf, -1, 1, 0, 0);
      >> f_opt, lam_g
      0.1456
      0.3407
      >> sol = S('p', 100); disp(sol.x)
      [1, 1]
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); y = ca.SX.sym("y"); p = ca.SX.sym("p");
      julia> nlp = Dict("x" => [x; y], "p" => p, "f" => (1 - x)^2 + p*(y - x^2)^2, "g" => x + y);
      julia> opts = Dict("ipopt.print_level" => 0, "ipopt.sb" => "yes", "print_time" => false);
      julia> S = ca.nlpsol("S", "ipopt", nlp, opts);
      julia> print(ca.name_in(S))
      ["x0", "p", "lbx", "ubx", "lbg", "ubg", "lam_x0", "lam_g0"]
      julia> print(ca.name_out(S))
      ["x", "f", "g", "lam_x", "lam_g", "lam_p"]
      julia> sol = S(x0=[2.5, 3.0], p=100, lbg=-1, ubg=1);
      julia> sol["f"]
      DM(0.1456)
      julia> sol["x"]
      DM([0.6188, 0.3812])
      julia> ca.stats(S)["success"], ca.stats(S)["return_status"]
      (true, "Solve_Succeeded")
      julia> x_opt, f_opt, g_opt, lam_x, lam_g, lam_p = S([2.5, 3.0], 100, -Inf, Inf, -1, 1, 0, 0);
      julia> f_opt, lam_g
      (0.1456, 0.3407)
      julia> S(p=100)["x"]
      DM([1, 1])
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); p = SX.sym("p")
      >>> opts = {"ipopt.print_level": 0, "ipopt.sb": "yes", "print_time": False}
      >>> S = nlpsol("S", "ipopt", {"x": x, "p": p, "f": (x - p)**2 + 3*p, "g": 2*x}, opts)
      >>> sol = S(p=1, lbg=4)
      >>> sol["x"], sol["lam_g"]
      (DM(2), DM(-1))
      >>> sol["lam_p"]
      DM(-1)
      >>> S(p=1, ubx=0)["lam_x"]
      DM(2)
      >>> pm = MX.sym("p")
      >>> dxdp = Function("dxdp", [pm], [jacobian(S(p=pm)["x"], pm)])
      >>> dxdp(5)
      DM(1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); p = SX.sym('p');
      >> opts = struct('print_time', false, 'ipopt', struct('print_level', 0, 'sb', 'yes'));
      >> S = nlpsol('S', 'ipopt', struct('x', x, 'p', p, 'f', (x - p)^2 + 3*p, 'g', 2*x), opts);
      >> sol = S('p', 1, 'lbg', 4); sol.x, sol.lam_g
      2
      -1
      >> sol.lam_p
      -1
      >> sol = S('p', 1, 'ubx', 0); sol.lam_x
      2
      >> pm = MX.sym('p'); sol = S('p', pm);
      >> dxdp = Function('dxdp', {pm}, {jacobian(sol.x, pm)}); disp(dxdp(5))
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
      julia> opts = Dict("ipopt.print_level" => 0, "ipopt.sb" => "yes", "print_time" => false);
      julia> nlp = Dict("x" => x, "p" => p, "f" => (x - p)^2 + 3*p, "g" => 2*x);
      julia> S = ca.nlpsol("S", "ipopt", nlp, opts);
      julia> sol = S(p=1, lbg=4);
      julia> sol["x"], sol["lam_g"]
      (2, -1)
      julia> sol["lam_p"]
      DM(-1)
      julia> S(p=1, ubx=0)["lam_x"]
      DM(2)
      julia> pm = ca.MX.sym("p");
      julia> dxdp = ca.Function("dxdp", [pm], [ca.jacobian(S(p=pm)["x"], pm)]);
      julia> dxdp(5)
      DM(1)
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> x = SX.sym("x")
      >>> opts = {"ipopt.print_level": 0, "ipopt.sb": "yes", "print_time": False}
      >>> S = nlpsol("S", "ipopt", {"x": x, "f": x, "g": x}, opts)
      >>> sol = S(lbx=0, ubx=1, lbg=2, ubg=3)
      >>> S.stats()["success"], S.stats()["return_status"]
      (False, 'Infeasible_Problem_Detected')
      >>> S = nlpsol("S", "ipopt", {"x": x, "f": x, "g": x}, dict(opts, error_on_fail=True))
      >>> try:
      ...     S(lbx=0, ubx=1, lbg=2, ubg=3)
      ... except RuntimeError as e:
      ...     print("nlpsol process failed" in str(e))
      True
      >>> try:
      ...     S(lbx=1, ubx=0)
      ... except RuntimeError as e:
      ...     print(str(e).splitlines()[-1])
      Ill-posed problem detected: LBX[0] <= UBX[0] was violated. Got LBX[0]=1 and UBX[0] = 0.
      >>> try:
      ...     nlpsol("S", "ipopt", {"x": x, "f": x, "h": x}, opts)
      ... except RuntimeError as e:
      ...     print("No such field: h" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x');
      >> opts = struct('print_time', false, 'ipopt', struct('print_level', 0, 'sb', 'yes'));
      >> S = nlpsol('S', 'ipopt', struct('x', x, 'f', x, 'g', x), opts);
      >> sol = S('lbx', 0, 'ubx', 1, 'lbg', 2, 'ubg', 3);
      >> st = S.stats(); disp(st.success), disp(st.return_status)
      0
      Infeasible_Problem_Detected
      >> opts.error_on_fail = true;
      >> S = nlpsol('S', 'ipopt', struct('x', x, 'f', x, 'g', x), opts);
      >> try, S('lbx', 0, 'ubx', 1, 'lbg', 2, 'ubg', 3); catch e, msg = e.message; end
      >> contains(msg, 'nlpsol process failed')
      1
      >> try, S('lbx', 1, 'ubx', 0); catch e, msg = e.message; end
      >> contains(msg, 'Ill-posed problem detected: LBX[0] <= UBX[0] was violated')
      1
      >> nlp = struct('x', x, 'f', x, 'h', x);
      >> try, nlpsol('S', 'ipopt', nlp, opts); catch e, msg = e.message; end
      >> contains(msg, 'No such field: h')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x");
      julia> opts = Dict("ipopt.print_level" => 0, "ipopt.sb" => "yes", "print_time" => false);
      julia> S = ca.nlpsol("S", "ipopt", Dict("x" => x, "f" => x, "g" => x), opts);
      julia> sol = S(lbx=0, ubx=1, lbg=2, ubg=3);
      julia> ca.stats(S)["success"], ca.stats(S)["return_status"]
      (false, "Infeasible_Problem_Detected")
      julia> S = ca.nlpsol("S", "ipopt", Dict("x" => x, "f" => x, "g" => x),
                           merge(opts, Dict("error_on_fail" => true)));
      julia> try
                 S(lbx=0, ubx=1, lbg=2, ubg=3)
             catch e
                 println(occursin("nlpsol process failed", sprint(showerror, e)))
             end
      true
      julia> try
                 S(lbx=1, ubx=0)
             catch e
                 println(last(split(sprint(showerror, e), "\n")))
             end
      Ill-posed problem detected: LBX[0] <= UBX[0] was violated. Got LBX[0]=1 and UBX[0] = 0.
      julia> try
                 ca.nlpsol("S", "ipopt", Dict("x" => x, "f" => x, "h" => x), opts)
             catch e
                 println(occursin("No such field: h", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \generalsection{Nlpsol}
      \pluginssection{Nlpsol}

      \author Joel Andersson
      \date 2011-2015

      \identifier{21q} */

  /** \defgroup nlpsol Title
  * @copydoc main_nlpsol
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_nlpsol
  * \endif
  */
  ///@{
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const SXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const MXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const std::string& fname, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Importer& compiler, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const NlpBuilder& nl, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Function& nlp, const Dict& opts=Dict());
  ///@}

  /** \brief Get the input scheme of NLP solvers

      Returns the names of the 8 inputs of a Function created by nlpsol, in order:
      x0, p, lbx, ubx, lbg, ubg, lam_x0, lam_g0.

  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif

      \lang{python}
      \doctest
      >>> nlpsol_in()
      ['x0', 'p', 'lbx', 'ubx', 'lbg', 'ubg', 'lam_x0', 'lam_g0']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_in()
      x0
      p
      lbx
      ubx
      lbg
      ubg
      lam_x0
      lam_g0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_in()
      8-element Vector{String}:
       "x0"
       "p"
       "lbx"
       "ubx"
       "lbg"
       "ubg"
       "lam_x0"
       "lam_g0"
      \enddoctest
      \endlang

      \sa nlpsol_out, nlpsol_n_in, nlpsol_default_in

      \identifier{1sy} */
  CASADI_EXPORT std::vector<std::string> nlpsol_in();

  /** \brief Get the output scheme of NLP solvers

      Returns the names of the 6 outputs of a Function created by nlpsol, in order:
      x, f, g, lam_x, lam_g, lam_p.

  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif

      \lang{python}
      \doctest
      >>> nlpsol_out()
      ['x', 'f', 'g', 'lam_x', 'lam_g', 'lam_p']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_out()
      x
      f
      g
      lam_x
      lam_g
      lam_p
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_out()
      6-element Vector{String}:
       "x"
       "f"
       "g"
       "lam_x"
       "lam_g"
       "lam_p"
      \enddoctest
      \endlang

      \sa nlpsol_in, nlpsol_n_out

      \identifier{1sz} */
  CASADI_EXPORT std::vector<std::string> nlpsol_out();

  /** \brief Get the name of an NLP solver input by index

      \param ind Index in 0..nlpsol_n_in()-1; any other value gives an empty string

  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif

      \lang{python}
      \doctest
      >>> nlpsol_in(0), nlpsol_in(7), repr(nlpsol_in(8))
      ('x0', 'lam_g0', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_in(0), nlpsol_in(7), isempty(nlpsol_in(8))
      x0
      lam_g0
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_in(0), ca.nlpsol_in(7), ca.nlpsol_in(8)
      ("x0", "lam_g0", "")
      \enddoctest
      \endlang

      \identifier{1t0} */
  CASADI_EXPORT std::string nlpsol_in(casadi_int ind);

  /** \brief Get the name of an NLP solver output by index

      \param ind Index in 0..nlpsol_n_out()-1; any other value gives an empty string

  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif

      \lang{python}
      \doctest
      >>> nlpsol_out(0), nlpsol_out(5), repr(nlpsol_out(6))
      ('x', 'lam_p', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_out(0), nlpsol_out(5), isempty(nlpsol_out(6))
      x
      lam_p
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_out(0), ca.nlpsol_out(5), ca.nlpsol_out(6)
      ("x", "lam_p", "")
      \enddoctest
      \endlang

      \identifier{1t1} */
  CASADI_EXPORT std::string nlpsol_out(casadi_int ind);

  /** \brief Number of NLP solver inputs

      \lang{python}
      \doctest
      >>> nlpsol_n_in(), len(nlpsol_in())
      (8, 8)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> fprintf('%d %d\n', nlpsol_n_in(), numel(cellstr(nlpsol_in())))
      8 8
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_n_in(), length(ca.nlpsol_in())
      (8, 8)
      \enddoctest
      \endlang

      \identifier{1t2} */
  CASADI_EXPORT casadi_int nlpsol_n_in();

  /** \brief Number of NLP solver outputs

      \lang{python}
      \doctest
      >>> nlpsol_n_out(), len(nlpsol_out())
      (6, 6)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> fprintf('%d %d\n', nlpsol_n_out(), numel(cellstr(nlpsol_out())))
      6 6
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_n_out(), length(ca.nlpsol_out())
      (6, 6)
      \enddoctest
      \endlang

      \identifier{1t3} */
  CASADI_EXPORT casadi_int nlpsol_n_out();

  ///@{
  /** \brief Default value of an NLP solver input

      The value used for an input that is not passed when calling an NLP solver Function:
      -inf for lbx and lbg, +inf for ubx and ubg, 0 for x0, p, lam_x0 and lam_g0.

      \param ind Input index as in nlpsol_in; without argument, the list of all 8 defaults

      \lang{python}
      \doctest
      >>> nlpsol_default_in()
      [0.0, 0.0, -inf, inf, -inf, inf, 0.0, 0.0]
      >>> nlpsol_default_in(3)
      inf
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_default_in()
      0 0 -Inf Inf -Inf Inf 0 0
      >> nlpsol_default_in(3)
      Inf
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_default_in()
      8-element Vector{Float64}:
         0.0
         0.0
       -Inf
        Inf
       -Inf
        Inf
         0.0
         0.0
      julia> ca.nlpsol_default_in(3)
      Inf
      \enddoctest
      \endlang

      \identifier{1t4} */
  CASADI_EXPORT double nlpsol_default_in(casadi_int ind);
  CASADI_EXPORT std::vector<double> nlpsol_default_in();
  ///@}

  /** \brief Get all options of an NLP solver plugin

      Returns the names of the plugin-specific options only; the options shared by all plugins
      (e.g. "print_time", "calc_lam_x") are not included. Raises an error for an unknown plugin.

      \lang{python}
      \doctest
      >>> opts = nlpsol_options("sqpmethod")
      >>> "max_iter" in opts, "qpsol" in opts, "print_time" in opts
      (True, True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opts = cellstr(nlpsol_options('sqpmethod'));
      >> ismember({'max_iter', 'qpsol', 'print_time'}, opts)
      1 1 0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opts = ca.nlpsol_options("sqpmethod");
      julia> "max_iter" in opts, "qpsol" in opts, "print_time" in opts
      (true, true, false)
      \enddoctest
      \endlang

      \sa nlpsol_option_type, nlpsol_option_info

      \identifier{1t5} */
  CASADI_EXPORT std::vector<std::string> nlpsol_options(const std::string& name);

  /** \brief Get the type of an option of an NLP solver plugin

      \param name Plugin name
      \param op Option name; an unknown option raises an error
      \return Type name such as "OT_INT", "OT_BOOL", "OT_DICT"

      \lang{python}
      \doctest
      >>> nlpsol_option_type("sqpmethod", "max_iter"), nlpsol_option_type("ipopt", "ipopt")
      ('OT_INT', 'OT_DICT')
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_option_type('sqpmethod', 'max_iter')
      OT_INT
      >> nlpsol_option_type('ipopt', 'ipopt')
      OT_DICT
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_option_type("sqpmethod", "max_iter"), ca.nlpsol_option_type("ipopt", "ipopt")
      ("OT_INT", "OT_DICT")
      \enddoctest
      \endlang

      \identifier{1t6} */
  CASADI_EXPORT std::string nlpsol_option_type(const std::string& name, const std::string& op);

  /** \brief Get the documentation string of an option of an NLP solver plugin

      \param name Plugin name
      \param op Option name; an unknown option raises an error

      \lang{python}
      \doctest
      >>> nlpsol_option_info("sqpmethod", "max_iter")
      'Maximum number of SQP iterations'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> nlpsol_option_info('sqpmethod', 'max_iter')
      Maximum number of SQP iterations
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.nlpsol_option_info("sqpmethod", "max_iter")
      "Maximum number of SQP iterations"
      \enddoctest
      \endlang

      \identifier{1t7} */
  CASADI_EXPORT std::string nlpsol_option_info(const std::string& name, const std::string& op);

  /** \brief Check if an NLP solver plugin is available

      Tries to load the plugin \p name if it is not loaded yet; false when it is not found.

      \lang{python}
      \doctest
      >>> has_nlpsol("ipopt"), has_nlpsol("no_such_plugin")
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> has_nlpsol('ipopt'), has_nlpsol('no_such_plugin')
      1
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.has_nlpsol("ipopt"), ca.has_nlpsol("no_such_plugin")
      (true, false)
      \enddoctest
      \endlang

      \identifier{2pb} */
  CASADI_EXPORT bool has_nlpsol(const std::string& name);

  /** \brief Explicitly load an NLP solver plugin

      Normally plugins are loaded on first use; this raises an error when the plugin cannot be
      found.

      \lang{python}
      \doctest
      >>> try:
      ...     load_nlpsol("no_such_plugin")
      ... except RuntimeError as e:
      ...     print("no_such_plugin" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> try, load_nlpsol('no_such_plugin'); catch e, msg = e.message; end
      >> contains(msg, 'no_such_plugin')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> try
                 ca.load_nlpsol("no_such_plugin")
             catch e
                 println(occursin("no_such_plugin", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2pc} */
  CASADI_EXPORT void load_nlpsol(const std::string& name);

  /** \brief Get the documentation string of an NLP solver plugin

      The text includes the plugin-specific option table.

      \lang{python}
      \doctest
      >>> "max_iter" in doc_nlpsol("sqpmethod")
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> contains(doc_nlpsol('sqpmethod'), 'max_iter')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> occursin("max_iter", ca.doc_nlpsol("sqpmethod"))
      true
      \enddoctest
      \endlang

      \identifier{2pd} */
  CASADI_EXPORT std::string doc_nlpsol(const std::string& name);

  /** @} */

#ifndef SWIG
/// Input arguments of an NLP function
enum NLPInput {
  /// Decision variable
  NL_X,
  /// Fixed parameter
  NL_P,
  /// Number of NLP inputs
  NL_NUM_IN
};

/// Shortname for onput arguments of an NLP function
const std::vector<std::string> NL_INPUTS = {"x", "p"};

/// Output arguments of an NLP function
enum NLPOutput {
  /// Objective function
  NL_F,
  /// Constraint function
  NL_G,
  /// Number of NLP outputs
  NL_NUM_OUT
};

/// Shortname for output arguments of an NLP function
const std::vector<std::string> NL_OUTPUTS = {"f", "g"};

/// Input arguments of an NLP Solver
enum NlpsolInput {
  /// Decision variables, initial guess (nx x 1), default 0
  NLPSOL_X0,
  /// Value of fixed parameters (np x 1), default 0
  NLPSOL_P,
  /// Decision variables lower bound (nx x 1), default -inf
  NLPSOL_LBX,
  /// Decision variables upper bound (nx x 1), default +inf
  NLPSOL_UBX,
  /// Constraints lower bound (ng x 1), default -inf
  NLPSOL_LBG,
  /// Constraints upper bound (ng x 1), default +inf
  NLPSOL_UBG,
  /// Lagrange multipliers for bounds on X, initial guess (nx x 1), default 0
  NLPSOL_LAM_X0,
  /// Lagrange multipliers for bounds on G, initial guess (ng x 1), default 0
  NLPSOL_LAM_G0,
  NLPSOL_NUM_IN
};

/// Output arguments of an NLP Solver
enum NlpsolOutput {
  /// Decision variables at the optimal solution (nx x 1)
  NLPSOL_X,
  /// Cost function value at the optimal solution (1 x 1)
  NLPSOL_F,
  /// Constraints function at the optimal solution (ng x 1)
  NLPSOL_G,
  /// Lagrange multipliers for bounds on X at the solution (nx x 1); negative for an active
  /// lower bound, positive for an active upper bound
  NLPSOL_LAM_X,
  /// Lagrange multipliers for bounds on G at the solution (ng x 1); same sign convention
  NLPSOL_LAM_G,
  /// Lagrange multipliers for p == P at the solution (np x 1): minus the gradient of the
  /// optimal cost with respect to p
  NLPSOL_LAM_P,
  NLPSOL_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_NLPSOL_HPP
