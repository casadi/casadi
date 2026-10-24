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


#ifndef CASADI_INTEGRATOR_HPP
#define CASADI_INTEGRATOR_HPP

#include "function.hpp"
#include "linsol.hpp"
#include "rootfinder.hpp"
#include "casadi_enum.hpp"

namespace casadi {

  /** \defgroup main_integrator Title
      \par

      Create an ODE/DAE integrator
      Solves an initial value problem (IVP) coupled to a terminal value problem
      with differential equation given as an implicit ODE coupled to an algebraic
      equation and a set of quadratures:

      \verbatim
      Initial conditions at t=t0
      x(t0)  = x0
      q(t0)  = 0

      Forward integration from t=t0 to t=tf
      der(x) = fx(x, z, p, u, t)                  Forward ODE
      0 = fz(x, z, p, u, t)                  Forward algebraic equations
      der(q) = fq(x, z, p, u, t)                  Forward quadratures

      Terminal conditions at t=tf
      rx(tf)  = rx0
      rq(tf)  = 0

      Backward integration from t=tf to t=t0
      der(rx) = gx(rx, rz, rp, x, z, p, t)        Backward ODE
      0 = gz(rx, rz, rp, x, z, p, t)        Backward algebraic equations
      der(rq) = gq(rx, rz, rp, x, z, p, t)        Backward quadratures

      where we assume that both the forward and backwards integrations are index-1
      (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
      gx, gz and gq have a linear dependency on rx, rz and rp.

      \endverbatim

      The DAE is given as a dictionary with the fields "x" (differential state), "ode"
      (its right-hand side, same size as "x") and optionally "t" (time, scalar), "z" (algebraic
      variables) with "alg" (residual, same size as "z"), "p" (parameters), "u" (controls, piecewise
      constant on the output grid), "quad" (quadrature right-hand side) and "zero" (event
      indicators, experimental). Every entry must be a dense vector; row vectors are flattened
      to columns. Any other field name is an error, as is a symbol in an output that is not an
      input. The DAE can also be given as a Function with the full input scheme t, x, z, p, u
      and output scheme ode, alg, quad, zero (see dyn_in and dyn_out), empty where unused.
      The pure ODE integrators (cvodes, rk, erk, rosenbrock) reject a problem with "z".

      Time grid: integrator(name, solver, dae, opts) integrates from t0 = 0 to tf = 1;
      integrator(name, solver, dae, t0, tf, opts) to a single output time tf; and
      integrator(name, solver, dae, t0, tout, opts) returns the solution at every time in the
      list \p tout. The grid must be increasing but need not contain t0; including it gives the
      initial state as the first output column. t0 > tf integrates backwards in time. The legacy
      options "t0", "tf", "grid" and "output_t0" are deprecated and emit a warning.

      The returned Function has inputs x0, z0, p, u, adj_xf, adj_zf, adj_qf (see integrator_in)
      and outputs xf, zf, qf, adj_x0, adj_z0, adj_p, adj_u (see integrator_out). With N output
      times, x0, z0 and p are column vectors, u has one column per output time (the control on
      the interval ending at that time; a single column is broadcast), and xf, zf, qf have one
      column per output time. z0 is the initial guess for the algebraic variables; it must be
      consistent for the Sundials solvers to converge. Omitted inputs are zero. The adj_* inputs
      and outputs are empty unless option "nadj" is positive: then adj_xf, adj_zf and adj_qf
      take one adjoint seed per output time and adj_x0, adj_p and adj_u return the corresponding
      adjoint sensitivities, i.e. adj_x0 is the gradient of sum(adj_xf .* xf) + sum(adj_qf .* qf)
      with respect to x0. Option "nfwd" adds, in the same way, forward sensitivity columns to
      every input and output. Without these options derivatives are obtained by embedding the
      call in an MX expression and differentiating (sensitivity equations are created on demand).

      Options common to all plugins (the plugin options follow in the tables below):
      - "nfwd", "nadj": number of forward and adjoint sensitivity directions built into the
        Function (default 0).
      - "expand": convert an MX DAE to SX before creating derivative functions (default false).
      - "print_stats": print integrator statistics after each call; stats() returns them.
      - "transition", "max_event_iter", "max_events", "event_tol": event handling for a DAE
        with "zero" indicators (experimental, see event_in and event_out).
      - "augmented_options": options passed to the integrator created for sensitivities.
      - "print_time", "verbose", "monitor", "common_options", "specific_options": diagnostics and
        options of the generated derivative functions, as for nlpsol.

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); p = SX.sym("p")
      >>> dae = {"x": x, "p": p, "ode": -p*x, "quad": x}
      >>> F = integrator("F", "cvodes", dae)
      >>> print(F)
      F:(x0,z0[0],p,u[0],adj_xf[],adj_zf[],adj_qf[])->(xf,zf[0],qf,adj_x0[],...,adj_u[]) ...
      >>> r = F(x0=1, p=1)
      >>> r["xf"], r["qf"]
      (DM(0.3679), DM(0.6321))
      >>> F = integrator("F", "rk", dae, 0, [1, 2, 3])
      >>> F(x0=1, p=1)["xf"]
      DM([[0.3679, 0.1353, 0.0498]])
      >>> xf, zf, qf, adj_x0, adj_z0, adj_p, adj_u = F(1, 0, 1, 0, 0, 0, 0)
      >>> xf.shape, zf.shape, adj_x0.shape
      ((1, 3), (0, 3), (0, 0))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); p = SX.sym('p');
      >> dae = struct('x', x, 'p', p, 'ode', -p*x, 'quad', x);
      >> F = integrator('F', 'cvodes', dae); disp(F)
      F:(x0,z0[0],p,u[0],adj_xf[],adj_zf[],adj_qf[])->(xf,zf[0],qf,adj_x0[],...,adj_u[]) ...
      >> r = F('x0', 1, 'p', 1); r.xf, r.qf
      0.3679
      0.6321
      >> F = integrator('F', 'rk', dae, 0, [1, 2, 3]);
      >> r = F('x0', 1, 'p', 1); r.xf
      [[0.3679, 0.1353, 0.0498]]
      >> [xf, zf, qf, adj_x0, adj_z0, adj_p, adj_u] = F(1, 0, 1, 0, 0, 0, 0);
      >> size(xf), size(zf), size(adj_x0)
      1 3
      0 3
      0 0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
      julia> dae = Dict("x" => x, "p" => p, "ode" => -p*x, "quad" => x);
      julia> F = ca.integrator("F", "cvodes", dae);
      julia> ca.n_in(F), ca.n_out(F), ca.size_in(F, 4)
      (7, 7, (0, 0))
      julia> r = F(x0=1, p=1);
      julia> r["xf"], r["qf"]
      (0.3679, 0.6321)
      julia> F = ca.integrator("F", "rk", dae, 0, [1, 2, 3]);
      julia> F(x0=1, p=1)["xf"]
      DM([[0.3679, 0.1353, 0.0498]])
      julia> xf, zf, qf, adj_x0, adj_z0, adj_p, adj_u = F(1, 0, 1, 0, 0, 0, 0);
      julia> size(xf), size(zf), size(adj_x0)
      ((1, 3), (0, 3), (0, 0))
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); z = SX.sym("z"); u = SX.sym("u")
      >>> F = integrator("F", "idas", {"x": x, "z": z, "ode": -z, "alg": z**2 - x})
      >>> r = F(x0=1, z0=1)
      >>> r["xf"], r["zf"]
      (DM(0.25), DM(0.5))
      >>> F(x0=1, z0=-1)["zf"]
      DM(-1.5)
      >>> F = integrator("F", "rk", {"x": x, "u": u, "ode": u}, 0, [1, 2, 3])
      >>> F(x0=0, u=[1, 2, 3])["xf"], F(x0=0, u=1)["xf"]
      (DM([[1, 3, 6]]), DM([[1, 2, 3]]))
      >>> try:
      ...     integrator("F", "cvodes", {"x": x, "z": z, "ode": -z, "alg": z - x})
      ... except RuntimeError as e:
      ...     print("does not support algebraic variables" in str(e))
      True
      >>> try:
      ...     integrator("F", "cvodes", {"x": x, "ode": vertcat(x, x)})
      ... except RuntimeError as e:
      ...     print("Dimension mismatch for 'ode'" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); z = SX.sym('z'); u = SX.sym('u');
      >> F = integrator('F', 'idas', struct('x', x, 'z', z, 'ode', -z, 'alg', z^2 - x));
      >> r = F('x0', 1, 'z0', 1); r.xf, r.zf
      0.25
      0.5
      >> r = F('x0', 1, 'z0', -1); r.zf
      -1.5
      >> F = integrator('F', 'rk', struct('x', x, 'u', u, 'ode', u), 0, [1, 2, 3]);
      >> r = F('x0', 0, 'u', [1 2 3]); disp(r.xf)
      [[1, 3, 6]]
      >> r = F('x0', 0, 'u', 1); disp(r.xf)
      [[1, 2, 3]]
      >> dae = struct('x', x, 'z', z, 'ode', -z, 'alg', z - x);
      >> try, integrator('F', 'cvodes', dae); catch e, msg = e.message; end
      >> contains(msg, 'does not support algebraic variables')
      1
      >> dae = struct('x', x, 'ode', [x; x]);
      >> try, integrator('F', 'cvodes', dae); catch e, msg = e.message; end
      >> contains(msg, 'Dimension mismatch for ''ode''')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); z = ca.SX.sym("z"); u = ca.SX.sym("u");
      julia> dae = Dict("x" => x, "z" => z, "ode" => -z, "alg" => z^2 - x);
      julia> F = ca.integrator("F", "idas", dae);
      julia> r = F(x0=1, z0=1);
      julia> r["xf"], r["zf"]
      (0.25, 0.5)
      julia> F(x0=1, z0=-1)["zf"]
      DM(-1.5)
      julia> F = ca.integrator("F", "rk", Dict("x" => x, "u" => u, "ode" => u), 0, [1, 2, 3]);
      julia> F(x0=0, u=[1 2 3])["xf"]
      DM([[1, 3, 6]])
      julia> F(x0=0, u=1)["xf"]
      DM([[1, 2, 3]])
      julia> try
                 ca.integrator("F", "cvodes", Dict("x" => x, "z" => z, "ode" => -z, "alg" => z - x))
             catch e
                 println(occursin("does not support algebraic variables", sprint(showerror, e)))
             end
      true
      julia> try
                 ca.integrator("F", "cvodes", Dict("x" => x, "ode" => [x; x]))
             catch e
                 println(occursin("Dimension mismatch for 'ode'", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); p = SX.sym("p")
      >>> dae = {"x": x, "p": p, "ode": -p*x, "quad": x}
      >>> F = integrator("F", "cvodes", dae, {"nadj": 1})
      >>> print(F)
      F:(x0,z0[0],p,u[0],adj_xf,adj_zf[0],adj_qf)->(xf,zf[0],qf,adj_x0,...,adj_p,adj_u[0]) ...
      >>> r = F(x0=1, p=1, adj_xf=1)
      >>> r["adj_x0"], r["adj_p"]
      (DM(0.3679), DM(-0.3679))
      >>> F = integrator("F", "cvodes", dae)
      >>> pm = MX.sym("p")
      >>> J = Function("J", [pm], [jacobian(F(x0=1, p=pm)["xf"], pm)])
      >>> J(1)
      DM(-0.3679)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); p = SX.sym('p');
      >> dae = struct('x', x, 'p', p, 'ode', -p*x, 'quad', x);
      >> F = integrator('F', 'cvodes', dae, struct('nadj', 1)); disp(F)
      F:(x0,z0[0],p,u[0],adj_xf,adj_zf[0],adj_qf)->(xf,zf[0],qf,adj_x0,...,adj_p,adj_u[0]) ...
      >> r = F('x0', 1, 'p', 1, 'adj_xf', 1); r.adj_x0, r.adj_p
      0.3679
      -0.3679
      >> F = integrator('F', 'cvodes', dae); pm = MX.sym('p'); r = F('x0', 1, 'p', pm);
      >> J = Function('J', {pm}, {jacobian(r.xf, pm)}); J(1)
      -0.3679
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
      julia> dae = Dict("x" => x, "p" => p, "ode" => -p*x, "quad" => x);
      julia> F = ca.integrator("F", "cvodes", dae, Dict("nadj" => 1));
      julia> ca.size_in(F, 4), ca.size_out(F, 5)
      ((1, 1), (1, 1))
      julia> r = F(x0=1, p=1, adj_xf=1);
      julia> r["adj_x0"], r["adj_p"]
      (0.3679, -0.3679)
      julia> F = ca.integrator("F", "cvodes", dae);
      julia> pm = ca.MX.sym("p");
      julia> J = ca.Function("J", [pm], [ca.jacobian(F(x0=1, p=pm)["xf"], pm)]);
      julia> J(1)
      DM(-0.3679)
      \enddoctest
      \endlang

      \generalsection{Integrator}
      \pluginssection{Integrator}

      \author Joel Andersson
      \date 2011-2015

      \identifier{21k} */
  /** \defgroup integrator Title
  * @copydoc main_integrator
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_integrator
  * \endif
  */
  ///@{
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, const std::vector<double>& tout, const Dict& opts=Dict());
#ifndef SWIGMATLAB
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, double tf, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, double tf, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, double tf, const Dict& opts=Dict());
#endif // SWIGMATLAB
  ///@}

  /** \brief Check if an integrator plugin is available

      Tries to load the plugin \p name if it is not loaded yet; false when it is not found.

      \lang{python}
      \doctest
      >>> has_integrator("rk"), has_integrator("no_such_plugin")
      (True, False)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> has_integrator('rk'), has_integrator('no_such_plugin')
      1
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.has_integrator("rk"), ca.has_integrator("no_such_plugin")
      (true, false)
      \enddoctest
      \endlang

      \identifier{2oi} */
  CASADI_EXPORT bool has_integrator(const std::string& name);

  /** \brief Explicitly load an integrator plugin

      Normally plugins are loaded on first use; this raises an error when the plugin cannot be
      found.

      \lang{python}
      \doctest
      >>> try:
      ...     load_integrator("no_such_plugin")
      ... except RuntimeError as e:
      ...     print("no_such_plugin" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> try, load_integrator('no_such_plugin'); catch e, msg = e.message; end
      >> contains(msg, 'no_such_plugin')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> try
                 ca.load_integrator("no_such_plugin")
             catch e
                 println(occursin("no_such_plugin", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2oj} */
  CASADI_EXPORT void load_integrator(const std::string& name);

  /** \brief Get the documentation string of an integrator plugin

      The text includes the plugin-specific option table.

      \lang{python}
      \doctest
      >>> "Runge-Kutta" in doc_integrator("rk")
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> contains(doc_integrator('rk'), 'Runge-Kutta')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> occursin("Runge-Kutta", ca.doc_integrator("rk"))
      true
      \enddoctest
      \endlang

      \identifier{2ok} */
  CASADI_EXPORT std::string doc_integrator(const std::string& name);

  /** \brief Get the input scheme of integrators

      Returns the names of the 7 inputs of a Function created by integrator, in order:
      x0, z0, p, u, adj_xf, adj_zf, adj_qf.

      \lang{python}
      \doctest
      >>> integrator_in()
      ['x0', 'z0', 'p', 'u', 'adj_xf', 'adj_zf', 'adj_qf']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_in()
      x0
      z0
      p
      u
      adj_xf
      adj_zf
      adj_qf
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_in()
      7-element Vector{String}:
       "x0"
       "z0"
       "p"
       "u"
       "adj_xf"
       "adj_zf"
       "adj_qf"
      \enddoctest
      \endlang

      \sa integrator_out, integrator_n_in

      \identifier{7b} */
  CASADI_EXPORT std::vector<std::string> integrator_in();

  /** \brief Get the output scheme of integrators

      Returns the names of the 7 outputs of a Function created by integrator, in order:
      xf, zf, qf, adj_x0, adj_z0, adj_p, adj_u.

      \lang{python}
      \doctest
      >>> integrator_out()
      ['xf', 'zf', 'qf', 'adj_x0', 'adj_z0', 'adj_p', 'adj_u']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_out()
      xf
      zf
      qf
      adj_x0
      adj_z0
      adj_p
      adj_u
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_out()
      7-element Vector{String}:
       "xf"
       "zf"
       "qf"
       "adj_x0"
       "adj_z0"
       "adj_p"
       "adj_u"
      \enddoctest
      \endlang

      \sa integrator_in, integrator_n_out

      \identifier{7c} */
  CASADI_EXPORT std::vector<std::string> integrator_out();

  /** \brief Get the name of an integrator input by index

      \param ind Index in 0..integrator_n_in()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> integrator_in(0), integrator_in(6), repr(integrator_in(7))
      ('x0', 'adj_qf', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_in(0), integrator_in(6), isempty(integrator_in(7))
      x0
      adj_qf
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_in(0), ca.integrator_in(6), ca.integrator_in(7)
      ("x0", "adj_qf", "")
      \enddoctest
      \endlang

      \identifier{7d} */
  CASADI_EXPORT std::string integrator_in(casadi_int ind);

  /** \brief Get the name of an integrator output by index

      \param ind Index in 0..integrator_n_out()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> integrator_out(0), integrator_out(6), repr(integrator_out(7))
      ('xf', 'adj_u', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_out(0), integrator_out(6), isempty(integrator_out(7))
      xf
      adj_u
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_out(0), ca.integrator_out(6), ca.integrator_out(7)
      ("xf", "adj_u", "")
      \enddoctest
      \endlang

      \identifier{7e} */
  CASADI_EXPORT std::string integrator_out(casadi_int ind);

  /** \brief Get the number of integrator inputs

      \lang{python}
      \doctest
      >>> integrator_n_in()
      7
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_n_in()
      7
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_n_in()
      7
      \enddoctest
      \endlang

      \identifier{7f} */
  CASADI_EXPORT casadi_int integrator_n_in();

  /** \brief Get the number of integrator outputs

      \lang{python}
      \doctest
      >>> integrator_n_out()
      7
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> integrator_n_out()
      7
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.integrator_n_out()
      7
      \enddoctest
      \endlang

      \identifier{7g} */
  CASADI_EXPORT casadi_int integrator_n_out();

  /** \brief Get the input scheme of a DAE function

      Returns the names of the 5 inputs of a DAE Function as expected by integrator, in order:
      t, x, z, p, u. These are also the dictionary keys accepted by integrator.

      \lang{python}
      \doctest
      >>> dyn_in()
      ['t', 'x', 'z', 'p', 'u']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_in()
      t
      x
      z
      p
      u
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_in()
      5-element Vector{String}:
       "t"
       "x"
       "z"
       "p"
       "u"
      \enddoctest
      \endlang

      \sa dyn_out, integrator

      \identifier{25p} */
  CASADI_EXPORT std::vector<std::string> dyn_in();

  /** \brief Get the output scheme of a DAE function

      Returns the names of the 4 outputs of a DAE Function as expected by integrator, in order:
      ode, alg, quad, zero.

      \lang{python}
      \doctest
      >>> dyn_out()
      ['ode', 'alg', 'quad', 'zero']
      >>> x = SX.sym("x"); t = SX.sym("t"); e = SX(0, 1)
      >>> dae = Function("dae", [t, x, e, e, e], [-x, e, x, e], dyn_in(), dyn_out())
      >>> integrator("F", "rk", dae)(x0=1)["xf"]
      DM(0.3679)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_out()
      ode
      alg
      quad
      zero
      >> x = SX.sym('x'); t = SX.sym('t'); e = SX(0, 1);
      >> dae = Function('dae', {t, x, e, e, e}, {-x, e, x, e}, dyn_in(), dyn_out());
      >> F = integrator('F', 'rk', dae); r = F('x0', 1); r.xf
      0.3679
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_out()
      4-element Vector{String}:
       "ode"
       "alg"
       "quad"
       "zero"
      julia> x = ca.SX.sym("x"); t = ca.SX.sym("t"); e = ca.SX(0, 1);
      julia> dae = ca.Function("dae", [t, x, e, e, e], [-x, e, x, e], ca.dyn_in(), ca.dyn_out());
      julia> ca.integrator("F", "rk", dae)(x0=1)["xf"]
      DM(0.3679)
      \enddoctest
      \endlang

      \sa dyn_in, integrator

      \identifier{25q} */
  CASADI_EXPORT std::vector<std::string> dyn_out();

  /** \brief Get the name of a DAE function input by index

      \param ind Index in 0..dyn_n_in()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> dyn_in(1), repr(dyn_in(5))
      ('x', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_in(1), isempty(dyn_in(5))
      x
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_in(1), ca.dyn_in(5)
      ("x", "")
      \enddoctest
      \endlang

      \identifier{25r} */
  CASADI_EXPORT std::string dyn_in(casadi_int ind);

  /** \brief Get the name of a DAE function output by index

      \param ind Index in 0..dyn_n_out()-1; any other value gives an empty string

      \lang{python}
      \doctest
      >>> dyn_out(0), repr(dyn_out(4))
      ('ode', "''")
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_out(0), isempty(dyn_out(4))
      ode
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_out(0), ca.dyn_out(4)
      ("ode", "")
      \enddoctest
      \endlang

      \identifier{25s} */
  CASADI_EXPORT std::string dyn_out(casadi_int ind);

  /** \brief Get the number of inputs of a DAE function

      \lang{python}
      \doctest
      >>> dyn_n_in()
      5
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_n_in()
      5
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_n_in()
      5
      \enddoctest
      \endlang

      \identifier{25t} */
  CASADI_EXPORT casadi_int dyn_n_in();

  /** \brief Get the number of outputs of a DAE function

      \lang{python}
      \doctest
      >>> dyn_n_out()
      4
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dyn_n_out()
      4
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.dyn_n_out()
      4
      \enddoctest
      \endlang

      \identifier{25u} */
  CASADI_EXPORT casadi_int dyn_n_out();

  /** \brief Get the input scheme of an event transition function

      Returns the names of the inputs of the Function passed as option "transition" to
      integrator: index (of the triggered "zero" indicator), t, x, z, p, u.

      \lang{python}
      \doctest
      >>> event_in()
      ['index', 't', 'x', 'z', 'p', 'u']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> event_in()
      index
      t
      x
      z
      p
      u
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.event_in()
      6-element Vector{String}:
       "index"
       "t"
       "x"
       "z"
       "p"
       "u"
      \enddoctest
      \endlang

      \sa event_out

      \identifier{2b4} */
  CASADI_EXPORT std::vector<std::string> event_in();

  /** \brief Get the output scheme of an event transition function

      Returns the names of the outputs of a transition Function: post_x, post_z, the values of
      x and z right after the event.

      \lang{python}
      \doctest
      >>> event_out()
      ['post_x', 'post_z']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> event_out()
      post_x
      post_z
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.event_out()
      2-element Vector{String}:
       "post_x"
       "post_z"
      \enddoctest
      \endlang

      \sa event_in

      \identifier{2b5} */
  CASADI_EXPORT std::vector<std::string> event_out();
  /** @} */

#ifndef SWIG
/// Inputs of the symbolic representation of the DAE
enum DynIn {
  DYN_T,
  DYN_X,
  DYN_Z,
  DYN_P,
  DYN_U,
  DYN_NUM_IN};

/// Outputs of the symbolic representation of the DAE
enum DynOut {
  DYN_ODE,
  DYN_ALG,
  DYN_QUAD,
  DYN_ZERO,
  DYN_NUM_OUT};

/// Inputs of an event transition function
enum EventIn {
  EVENT_INDEX,
  EVENT_T,
  EVENT_X,
  EVENT_Z,
  EVENT_P,
  EVENT_U,
  EVENT_NUM_IN};

/// Outputs of an event transition function
enum EventOut {
  EVENT_POST_X,
  EVENT_POST_Z,
  EVENT_NUM_OUT};

/// Input arguments of an integrator
enum IntegratorInput {
  /// Differential state at the initial time
  INTEGRATOR_X0,
  /// Initial guess for the algebraic variable at the initial time
  INTEGRATOR_Z0,
  /// Parameters
  INTEGRATOR_P,
  /// Piecewise constant control, a new control interval starts at each output time
  INTEGRATOR_U,
  /// Adjoint seeds corresponding to the states at the output times
  INTEGRATOR_ADJ_XF,
  /// Adjoint seeds corresponding to the algebraic variables at the output times
  INTEGRATOR_ADJ_ZF,
  /// Adjoint seeds corresponding to the quadratures at the output times
  INTEGRATOR_ADJ_QF,
  /// Number of input arguments of an integrator
  INTEGRATOR_NUM_IN
};

/// Output arguments of an integrator
enum IntegratorOutput {
  /// Differential state at all output times
  INTEGRATOR_XF,
  /// Algebraic variable at all output times
  INTEGRATOR_ZF,
  /// Quadrature state at all output times
  INTEGRATOR_QF,
  /// Adjoint sensitivities corresponding to the initial state
  INTEGRATOR_ADJ_X0,
  /// Adjoint sensitivities corresponding to the algebraic variable guess
  INTEGRATOR_ADJ_Z0,
  /// Adjoint sensitivities corresponding to the parameter vector
  INTEGRATOR_ADJ_P,
  /// Adjoint sensitivities corresponding to the control vector
  INTEGRATOR_ADJ_U,
  /// Number of output arguments of an integrator
  INTEGRATOR_NUM_OUT
};

///@{
/// Number of entries in enums
template<> struct enum_traits<DynIn> {
  static const size_t n_enum = DYN_NUM_IN;
};
template<> struct enum_traits<DynOut> {
  static const size_t n_enum = DYN_NUM_OUT;
};
template<> struct enum_traits<EventIn> {
  static const size_t n_enum = EVENT_NUM_IN;
};
template<> struct enum_traits<EventOut> {
  static const size_t n_enum = EVENT_NUM_OUT;
};
///@}

///@{
/// Convert to string
CASADI_EXPORT std::string to_string(DynIn v);
CASADI_EXPORT std::string to_string(DynOut v);
CASADI_EXPORT std::string to_string(EventIn v);
CASADI_EXPORT std::string to_string(EventOut v);
///@}

#endif // SWIG

} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
