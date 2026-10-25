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


#ifndef CASADI_CONIC_HPP
#define CASADI_CONIC_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_conic Title
      \par

      Create a QP solver

      Returns a Function that solves the convex quadratic program (QP)

      \verbatim
      min          1/2 x' H x + g' x
      x

      subject to
      LBA <= A x <= UBA
      LBX <= x   <= UBX

      resize(Q x, np, np) + P >= 0 (psd)

      with :
      H sparse (n x n) positive definite
      g dense  (n x 1)
      A sparse (nc x n)
      Q sparse symmetric (np^2 x n)
      P sparse symmetric (np x nq)

      n: number of decision variables (x)
      nc: number of constraints (A)
      nq: shape of psd constraint matrix

      \endverbatim

      Whether an indefinite H is detected depends on the plugin: some raise an error, others
      return a stationary point. Second-order cone constraints can be added as psd constraints
      through the helper function soc: || x ||_2 <= y is equivalent to soc(x, y) being psd, with
      soc(x, y) = [y*I x; x' y] (Schur complement). Only some plugins support psd constraints.

      conic is the low-level form: the problem is given by the sparsity patterns of "h" and "a"
      (and optionally "q", "p"), and the numerical values are passed at every call. At least one
      of "h" and "a" is needed; a missing "h" means a linear objective (LP) and a missing "a" means
      no linear constraints. "h" must have a square, symmetric sparsity pattern. Inputs are h, g,
      a, lba, uba, lbx, ubx, x0, lam_x0, lam_a0, q, p (see conic_in); outputs are x, cost, lam_a,
      lam_x (see conic_out). Omitted inputs default to -inf (lba, lbx), +inf (uba, ubx) or 0.

      qpsol is the symbolic form: a dictionary with fields "x" (dense vector), "f" (scalar),
      optionally "p" (parameters) and "g" (dense vector), as for nlpsol. The matrices H, g, A and
      the constant terms are obtained by symbolic differentiation and evaluated at every call at
      x0 and p; "f" must be quadratic and "g" linear in "x", this is NOT checked: other
      expressions are silently replaced by their second order (f) or first order (g) expansion
      around x0. The returned Function has the nlpsol input and output scheme (x0, p, lbx, ubx,
      lbg, ubg, lam_x0, lam_g0 -> x, f, g, lam_x, lam_g, lam_p); lam_p is always NaN.

      Multiplier convention as for nlpsol: a multiplier is negative when the lower bound is
      active and positive when the upper bound is active, such that H x + g + A' lam_a + lam_x = 0.

      To silence a solver, set "print_time" False plus the plugin's own options, e.g.
      {"print_iter": False, "print_header": False, "print_info": False} for qrqp or
      {"printLevel": "none"} for qpoases (its copyright notice is printed regardless).

      Options common to all plugins ("error_on_fail" defaults to true for QP solvers, so a failed
      solve raises; with false the plugin result is returned and stats() reports "success"):
      - "discrete": list of n booleans marking integer variables; needs a plugin with integer
        support (e.g. cplex, gurobi, highs), otherwise an error at construction.
      - "equality": list of nc booleans hinting which linear constraints have equal bounds.
      - "print_problem": print the numerical problem data at every call (default false).
      - "print_time": print a timing table after each call (default false).
      - "expand" (qpsol only): convert an MX problem to SX (default false).
      - "sens_min_lam", "sens_linsol", "sens_linsol_options": used when differentiating the
        solver Function (a call embedded in an MX expression): a bound is treated as active when
        its multiplier exceeds "sens_min_lam" in magnitude (default 0) and the KKT system is
        solved with linear solver "sens_linsol" (default "qr").
      The plugin options follow in the tables below.

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); y = SX.sym("y")
      >>> qp = {"x": vertcat(x, y), "f": x**2 + y**2 + x*y - 3*x, "g": x + y}
      >>> opts = {"print_iter": False, "print_header": False, "print_info": False}
      >>> S = qpsol("S", "qrqp", qp, opts)
      >>> S.name_in()
      ['x0', 'p', 'lbx', 'ubx', 'lbg', 'ubg', 'lam_x0', 'lam_g0']
      >>> S.name_out()
      ['x', 'f', 'g', 'lam_x', 'lam_g', 'lam_p']
      >>> sol = S(lbg=1, ubg=1, lbx=0)
      >>> sol["x"], sol["f"]
      (DM([1, 0]), DM(-2))
      >>> sol["lam_g"], sol["lam_x"]
      (DM(1), DM([0, -2]))
      >>> sol["g"], sol["lam_p"].shape
      (DM(1), (0, 0))
      >>> S.stats()["success"], S.stats()["return_status"]
      (True, 'success')
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); y = SX.sym('y');
      >> qp = struct('x', vertcat(x, y), 'f', x^2 + y^2 + x*y - 3*x, 'g', x + y);
      >> opts = struct('print_iter', false, 'print_header', false, 'print_info', false);
      >> S = qpsol('S', 'qrqp', qp, opts);
      >> S.name_in()
      x0
      p
      lbx
      ubx
      lbg
      ubg
      lam_x0
      lam_g0
      >> S.name_out()
      x
      f
      g
      lam_x
      lam_g
      lam_p
      >> sol = S('lbg', 1, 'ubg', 1, 'lbx', 0);
      >> sol.x, sol.f
      [1, 0]
      -2
      >> sol.lam_g, sol.lam_x
      1
      [0, -2]
      >> sol.g, size(sol.lam_p)
      1
      0 0
      >> st = S.stats(); disp(st.success), disp(st.return_status)
      1
      success
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
      julia> qp = Dict("x" => [x; y], "f" => x^2 + y^2 + x*y - 3*x, "g" => x + y);
      julia> opts = Dict("print_iter" => false, "print_header" => false, "print_info" => false);
      julia> S = ca.qpsol("S", "qrqp", qp, opts);
      julia> ca.name_in(S)
      8-element Vector{String}:
       "x0"
       "p"
       "lbx"
       "ubx"
       "lbg"
       "ubg"
       "lam_x0"
       "lam_g0"
      julia> ca.name_out(S)
      6-element Vector{String}:
       "x"
       "f"
       "g"
       "lam_x"
       "lam_g"
       "lam_p"
      julia> sol = S(lbg=1, ubg=1, lbx=0);
      julia> sol["x"], sol["f"]
      ([1, 0], -2)
      julia> sol["lam_g"], sol["lam_x"]
      (1, [0, -2])
      julia> sol["g"], size(sol["lam_p"])
      (1, (0, 0))
      julia> ca.stats(S)["success"], ca.stats(S)["return_status"]
      (true, "success")
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> H = DM([[2, 1], [1, 2]]); g = DM([-3, 0]); A = DM([[1, 1]])
      >>> opts = {"print_iter": False, "print_header": False, "print_info": False}
      >>> C = conic("C", "qrqp", {"h": H.sparsity(), "a": A.sparsity()}, opts)
      >>> print(C)
      C:(h[2x2],g[2],a[1x2],lba,uba,lbx[2],ubx[2],x0[2],...)->(x[2],cost,lam_a,lam_x[2]) Qrqp
      >>> sol = C(h=H, g=g, a=A, lba=1, uba=1, lbx=0)
      >>> sol["x"], sol["cost"]
      (DM([1, 0]), DM(-2))
      >>> sol["lam_a"], sol["lam_x"]
      (DM(1), DM([0, -2]))
      >>> LP = conic("LP", "qrqp", {"a": A.sparsity()}, opts)
      >>> LP(g=g, a=A, lba=1, uba=1, lbx=0, ubx=1)["x"]
      DM([1, 0])
      >>> try:
      ...     conic("C", "qrqp", {"h": Sparsity.upper(2)}, opts)
      ... except RuntimeError as e:
      ...     print("square & symmetric" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> H = DM([2 1; 1 2]); g = DM([-3; 0]); A = DM([1 1]);
      >> opts = struct('print_iter', false, 'print_header', false, 'print_info', false);
      >> C = conic('C', 'qrqp', struct('h', H.sparsity(), 'a', A.sparsity()), opts); disp(C)
      C:(h[2x2],g[2],a[1x2],lba,uba,lbx[2],ubx[2],x0[2],...)->(x[2],cost,lam_a,lam_x[2]) Qrqp
      >> sol = C('h', H, 'g', g, 'a', A, 'lba', 1, 'uba', 1, 'lbx', 0);
      >> sol.x, sol.cost
      [1, 0]
      -2
      >> sol.lam_a, sol.lam_x
      1
      [0, -2]
      >> LP = conic('LP', 'qrqp', struct('a', A.sparsity()), opts);
      >> sol = LP('g', g, 'a', A, 'lba', 1, 'uba', 1, 'lbx', 0, 'ubx', 1); disp(sol.x)
      [1, 0]
      >> qp = struct('h', Sparsity.upper(2));
      >> try, conic('C', 'qrqp', qp, opts); catch e, msg = e.message; end
      >> contains(msg, 'square & symmetric')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> H = ca.DM([2 1; 1 2]); g = ca.DM([-3, 0]); A = ca.DM([1 1]);
      julia> opts = Dict("print_iter" => false, "print_header" => false, "print_info" => false);
      julia> C = ca.conic("C", "qrqp", Dict("h" => ca.sparsity(H), "a" => ca.sparsity(A)), opts);
      julia> ca.class_name(C), ca.n_in(C), ca.name_out(C)
      ("Qrqp", 12, ["x", "cost", "lam_a", "lam_x"])
      julia> sol = C(h=H, g=g, a=A, lba=1, uba=1, lbx=0);
      julia> sol["x"], sol["cost"]
      ([1, 0], -2)
      julia> sol["lam_a"], sol["lam_x"]
      (1, [0, -2])
      julia> LP = ca.conic("LP", "qrqp", Dict("a" => ca.sparsity(A)), opts);
      julia> LP(g=g, a=A, lba=1, uba=1, lbx=0, ubx=1)["x"]
      DM([1, 0])
      julia> try
                 ca.conic("C", "qrqp", Dict("h" => ca.Sparsity.upper(2)), opts)
             catch e
                 println(occursin("square & symmetric", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> x = SX.sym("x")
      >>> opts = {"print_iter": False, "print_header": False, "print_info": False}
      >>> S = qpsol("S", "qrqp", {"x": x, "f": x**2, "g": x}, dict(opts, error_on_fail=False))
      >>> sol = S(lbx=0, ubx=1, lbg=2, ubg=3)
      >>> S.stats()["success"]
      False
      >>> S = qpsol("S", "qrqp", {"x": x, "f": x**2 - 4*x, "g": 2*x}, opts)
      >>> S(ubx=1)["lam_x"], S(lbx=3)["lam_x"], S(ubg=2)["lam_g"]
      (DM(2), DM(-2), DM(1))
      >>> try:
      ...     qpsol("S", "qrqp", {"x": x, "f": x**2}, dict(opts, discrete=[True]))
      ... except RuntimeError as e:
      ...     print("integer support" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); qp = struct('x', x, 'f', x^2, 'g', x);
      >> opts = struct('print_iter', false, 'print_header', false, 'print_info', false);
      >> S = qpsol('S', 'qrqp', qp, setfield(opts, 'error_on_fail', false));
      >> sol = S('lbx', 0, 'ubx', 1, 'lbg', 2, 'ubg', 3); st = S.stats(); disp(st.success)
      0
      >> S = qpsol('S', 'qrqp', struct('x', x, 'f', x^2 - 4*x, 'g', 2*x), opts);
      >> r1 = S('ubx', 1); r2 = S('lbx', 3); r3 = S('ubg', 2);
      >> r1.lam_x, r2.lam_x, r3.lam_g
      2
      -2
      1
      >> opts.discrete = true;
      >> try, qpsol('S', 'qrqp', qp, opts); catch e, msg = e.message; end
      >> contains(msg, 'integer support')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x");
      julia> opts = Dict("print_iter" => false, "print_header" => false, "print_info" => false);
      julia> qp = Dict("x" => x, "f" => x^2, "g" => x);
      julia> S = ca.qpsol("S", "qrqp", qp, merge(opts, Dict("error_on_fail" => false)));
      julia> sol = S(lbx=0, ubx=1, lbg=2, ubg=3);
      julia> ca.stats(S)["success"]
      false
      julia> S = ca.qpsol("S", "qrqp", Dict("x" => x, "f" => x^2 - 4*x, "g" => 2*x), opts);
      julia> S(ubx=1)["lam_x"], S(lbx=3)["lam_x"], S(ubg=2)["lam_g"]
      (2, -2, 1)
      julia> try
                 ca.qpsol("S", "qrqp", qp, merge(opts, Dict("discrete" => [true])))
             catch e
                 println(occursin("integer support", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \generalsection{Conic}
      \pluginssection{Conic}

      \author Joel Andersson
      \date 2011-2015

      \identifier{21n} */

  /** \defgroup conic Title
  * @copydoc main_conic
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_conic
  * \endif
  */
  ///@{
  CASADI_EXPORT Function conic(const std::string& name, const std::string& solver,
                               const SpDict& qp, const Dict& opts=Dict());
  CASADI_EXPORT Function qpsol(const std::string& name, const std::string& solver,
                               const SXDict& qp, const Dict& opts=Dict());
  CASADI_EXPORT Function qpsol(const std::string& name, const std::string& solver,
                               const MXDict& qp, const Dict& opts=Dict());
  ///@}

  /** \brief Get the input scheme of QP solvers

      Returns the names of the 12 inputs of a Function created by conic, in order:
      h, g, a, lba, uba, lbx, ubx, x0, lam_x0, lam_a0, q, p. A Function created by qpsol
      uses the nlpsol scheme instead.

      \lang{python}
      \doctest
      >>> conic_in()
      ['h', 'g', 'a', 'lba', 'uba', 'lbx', 'ubx', 'x0', 'lam_x0', 'lam_a0', 'q', 'p']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_in()
      h
      g
      a
      lba
      uba
      lbx
      ubx
      x0
      lam_x0
      lam_a0
      q
      p
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_in()
      12-element Vector{String}:
       "h"
       "g"
       "a"
       "lba"
       "uba"
       "lbx"
       "ubx"
       "x0"
       "lam_x0"
       "lam_a0"
       "q"
       "p"
      \enddoctest
      \endlang

      \sa conic_out, conic_n_in, nlpsol_in

      \identifier{1ee} */
  CASADI_EXPORT std::vector<std::string> conic_in();

  /** \brief Get the output scheme of QP solvers

      Returns the names of the 4 outputs of a Function created by conic, in order:
      x, cost, lam_a, lam_x.

      \lang{python}
      \doctest
      >>> conic_out()
      ['x', 'cost', 'lam_a', 'lam_x']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_out()
      x
      cost
      lam_a
      lam_x
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_out()
      4-element Vector{String}:
       "x"
       "cost"
       "lam_a"
       "lam_x"
      \enddoctest
      \endlang

      \sa conic_in, conic_n_out

      \identifier{1ef} */
  CASADI_EXPORT std::vector<std::string> conic_out();

  /** \brief Get the name of a QP solver input by index

      \param ind Index in 0..conic_n_in()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> conic_in(0), conic_in(11), repr(conic_in(12))
      ('h', 'p', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_in(0), conic_in(11), isempty(conic_in(12))
      h
      p
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_in(0), ca.conic_in(11), ca.conic_in(12)
      ("h", "p", "")
      \enddoctest
      \endlang

      \identifier{1eg} */
  CASADI_EXPORT std::string conic_in(casadi_int ind);

  /** \brief Get the name of a QP solver output by index

      \param ind Index in 0..conic_n_out()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> conic_out(1), repr(conic_out(4))
      ('cost', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_out(1), isempty(conic_out(4))
      cost
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_out(1), ca.conic_out(4)
      ("cost", "")
      \enddoctest
      \endlang

      \identifier{1eh} */
  CASADI_EXPORT std::string conic_out(casadi_int ind);

  /** \brief Get the number of QP solver inputs

      \lang{python}
      \doctest
      >>> conic_n_in()
      12
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_n_in()
      12
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_n_in()
      12
      \enddoctest
      \endlang

      \identifier{1ei} */
  CASADI_EXPORT casadi_int conic_n_in();

  /** \brief Get the number of QP solver outputs

      \lang{python}
      \doctest
      >>> conic_n_out()
      4
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_n_out()
      4
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_n_out()
      4
      \enddoctest
      \endlang

      \identifier{1ej} */
  CASADI_EXPORT casadi_int conic_n_out();

  /** \brief Get all options of a QP solver plugin

      Returns the names of the plugin-specific options only; the options shared by all plugins
      (e.g. "print_time", "discrete") are not included. Raises an error for an unknown plugin.

      \lang{python}
      \doctest
      >>> "max_iter" in conic_options("qrqp"), "discrete" in conic_options("qrqp")
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> opts = cellstr(conic_options('qrqp'));
      >> ismember({'max_iter', 'discrete'}, opts)
      1 0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> "max_iter" in ca.conic_options("qrqp"), "discrete" in ca.conic_options("qrqp")
      (true, false)
      \enddoctest
      \endlang

      \sa conic_option_type, conic_option_info

      \identifier{1ek} */
  CASADI_EXPORT std::vector<std::string> conic_options(const std::string& name);

  /** \brief Get the type of an option of a QP solver plugin

      \param name Plugin name
      \param op Option name; an unknown option raises an error
      \return Type name such as "OT_INT", "OT_BOOL", "OT_STRING"

      \lang{python}
      \doctest
      >>> conic_option_type("qrqp", "max_iter"), conic_option_type("qrqp", "print_iter")
      ('OT_INT', 'OT_BOOL')
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_option_type('qrqp', 'max_iter'), conic_option_type('qrqp', 'print_iter')
      OT_INT
      OT_BOOL
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_option_type("qrqp", "max_iter"), ca.conic_option_type("qrqp", "print_iter")
      ("OT_INT", "OT_BOOL")
      \enddoctest
      \endlang

      \identifier{1el} */
  CASADI_EXPORT std::string conic_option_type(const std::string& name, const std::string& op);

  /** \brief Get the documentation string of an option of a QP solver plugin

      \param name Plugin name
      \param op Option name; an unknown option raises an error

      \lang{python}
      \doctest
      >>> conic_option_info("qrqp", "max_iter")
      'Maximum number of iterations [1000].'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> conic_option_info('qrqp', 'max_iter')
      Maximum number of iterations [1000].
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.conic_option_info("qrqp", "max_iter")
      "Maximum number of iterations [1000]."
      \enddoctest
      \endlang

      \identifier{1em} */
  CASADI_EXPORT std::string conic_option_info(const std::string& name, const std::string& op);

  /** \brief Check if a QP solver plugin is available

      Tries to load the plugin \p name if it is not loaded yet; false when it is not found.

      \lang{python}
      \doctest
      >>> has_conic("qrqp"), has_conic("no_such_plugin")
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> has_conic('qrqp'), has_conic('no_such_plugin')
      1
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.has_conic("qrqp"), ca.has_conic("no_such_plugin")
      (true, false)
      \enddoctest
      \endlang

      \identifier{2od} */
  CASADI_EXPORT bool has_conic(const std::string& name);

  /** \brief Explicitly load a QP solver plugin

      Normally plugins are loaded on first use; this raises an error when the plugin cannot be
      found.

      \lang{python}
      \doctest
      >>> try:
      ...     load_conic("no_such_plugin")
      ... except RuntimeError as e:
      ...     print("no_such_plugin" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> try, load_conic('no_such_plugin'); catch e, msg = e.message; end
      >> contains(msg, 'no_such_plugin')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> try
                 ca.load_conic("no_such_plugin")
             catch e
                 println(occursin("no_such_plugin", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2oe} */
  CASADI_EXPORT void load_conic(const std::string& name);

  /** \brief Get the documentation string of a QP solver plugin

      The text includes the plugin-specific option table.

      \lang{python}
      \doctest
      >>> "active-set" in doc_conic("qrqp")
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> contains(doc_conic('qrqp'), 'active-set')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> occursin("active-set", ca.doc_conic("qrqp"))
      true
      \enddoctest
      \endlang

      \identifier{2of} */
  CASADI_EXPORT std::string doc_conic(const std::string& name);

  /** \brief Generate native code in the interfaced language for debugging

      Writes a standalone problem description in the solver's own language (e.g. a script for
      the solver's native API) for a Function created by conic. Only some plugins implement
      this; the others raise an error.

      \param f Function created by conic (not qpsol)
      \param filename Output file name

      \lang{python}
      \doctest
      >>> import os, tempfile
      >>> C = conic("C", "qrqp", {"h": Sparsity.dense(1)}, {"print_header": False})
      >>> d = tempfile.mkdtemp(); fn = os.path.join(d, "C_debug.txt")
      >>> try:
      ...     conic_debug(C, fn)
      ... except RuntimeError as e:
      ...     print("generateNativeCode not defined" in str(e))
      True
      >>> os.remove(fn); os.rmdir(d)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> C = conic('C', 'qrqp', struct('h', Sparsity.dense(1)), struct('print_header', false));
      >> fn = [tempname '.txt'];
      >> try, conic_debug(C, fn); catch e, msg = e.message; end
      >> contains(msg, 'generateNativeCode not defined')
      1
      >> delete(fn)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> opts = Dict("print_header" => false);
      julia> C = ca.conic("C", "qrqp", Dict("h" => ca.Sparsity.dense(1)), opts);
      julia> d = mktempdir(); fn = joinpath(d, "C_debug.txt");
      julia> try
                 ca.conic_debug(C, fn)
             catch e
                 println(occursin("generateNativeCode not defined", sprint(showerror, e)))
             end
      true
      julia> rm(d, recursive=true)
      \enddoctest
      \endlang

      \identifier{2og} */
  CASADI_EXPORT void conic_debug(const Function& f, const std::string &filename);

  /** \brief Generate native code in the interfaced language for debugging, to a stream

      \identifier{2oh} */
  CASADI_EXPORT void conic_debug(const Function& f, std::ostream &file);

  /** @} */

#ifndef SWIG
/// Input arguments of a QP problem
enum ConicInput {
  /// The square matrix H: sparse, (n x n). Only the lower triangular part is actually used.
  /// The matrix is assumed to be symmetrical.
  CONIC_H,
  /// The vector g: dense,  (n x 1)
  CONIC_G,
  /// The matrix A: sparse, (nc x n) - product with x must be dense.
  CONIC_A,
  /// dense, (nc x 1), default -inf
  CONIC_LBA,
  /// dense, (nc x 1), default +inf
  CONIC_UBA,
  /// dense, (n x 1), default -inf
  CONIC_LBX,
  /// dense, (n x 1), default +inf
  CONIC_UBX,
  /// dense, (n x 1), initial guess, default 0
  CONIC_X0,
  /// dense, (n x 1), initial guess for the multipliers of the simple bounds, default 0
  CONIC_LAM_X0,
  /// dense, (nc x 1), initial guess for the multipliers of the linear constraints, default 0
  CONIC_LAM_A0,
  /// The matrix Q: sparse symmetric, (np^2 x n)
  CONIC_Q,
  /// The matrix P: sparse symmetric, (np x np)
  CONIC_P,
  CONIC_NUM_IN};

/// Output arguments of an QP Solver
enum ConicOutput {
  /// The primal solution
  CONIC_X,
  /// The optimal cost
  CONIC_COST,
  /// The dual solution corresponding to linear bounds; negative for an active lower bound,
  /// positive for an active upper bound
  CONIC_LAM_A,
  /// The dual solution corresponding to simple bounds; same sign convention
  CONIC_LAM_X,
  CONIC_NUM_OUT};
#endif // SWIG

} // namespace casadi

#endif // CASADI_CONIC_HPP
