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

    A DaeBuilder holds a list of uniquely named model variables with FMI-style meta data
    (causality, variability, type, start value, bounds, unit, ...) together with the model
    equations. The model is either built symbolically with `add`, `eq`, `when` and friends,
    or imported from an FMU / `modelDescription.xml` (constructor with a path,
    `load_fmi_description`). Numerical evaluation goes through the function factory
    `create`, which returns differentiable Function objects.

    Every variable belongs to at most one CasADi category, derived from its causality,
    variability and defining equation:

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

    The category names (`"x"`, `"ode"`, ...) are the strings accepted by `all`, `inputs`,
    `outputs`, `create` and the other category-based members. Internal helper variables
    (named `der(x)`, `__alg__0`, `__assign__..`, `__when__0`, ...) are created as equations
    are added; they show up in `all()`.

    \lang{python}
    \doctest
    >>> dae = DaeBuilder("rocket")
    >>> a = dae.add("a", "parameter", "tunable")
    >>> u = dae.add("u", "input")
    >>> h = dae.add("h", dict(start=0))
    >>> v = dae.add("v", dict(start=0, unit="m/s"))
    >>> m = dae.add("m", dict(start=1))
    >>> dae.eq(dae.der(h), v)
    >>> dae.eq(dae.der(v), (u - a*v**2)/m - 9.81)
    >>> dae.eq(dae.der(m), -u**2)
    >>> dae.x(), dae.u(), dae.p()
    (['h', 'v', 'm'], ['u'], ['a'])
    >>> f = dae.create("f", ["x", "u", "p"], ["ode"])
    >>> f
    Function(f:(x[3],u,p)->(ode[3]) MXFunction)
    >>> f([0, 0, 1], 1, 0.5)
    DM([0, -8.81, -1])
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> dae = DaeBuilder('rocket');
    >> a = dae.add('a', 'parameter', 'tunable'); u = dae.add('u', 'input');
    >> h = dae.add('h', struct('start', 0)); v = dae.add('v', struct('start', 0, 'unit', 'm/s'));
    >> m = dae.add('m', struct('start', 1));
    >> dae.eq(dae.der(h), v); dae.eq(dae.der(v), (u - a*v^2)/m - 9.81); dae.eq(dae.der(m), -u^2);
    >> dae.x()
    h
    v
    m
    >> dae.u(), dae.p()
    u
    a
    >> f = dae.create('f', {'x', 'u', 'p'}, {'ode'}); disp(f)
    f:(x[3],u,p)->(ode[3]) MXFunction
    >> f([0; 0; 1], 1, 0.5)
    [0, -8.81, -1]
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> dae = ca.DaeBuilder("rocket");
    julia> a = ca.add(dae, "a", "parameter", "tunable");
    julia> u = ca.add(dae, "u", "input");
    julia> h = ca.add(dae, "h", Dict("start" => 0));
    julia> v = ca.add(dae, "v", Dict("start" => 0, "unit" => "m/s"));
    julia> m = ca.add(dae, "m", Dict("start" => 1));
    julia> ca.eq(dae, ca.der(dae, h), v)
    julia> ca.eq(dae, ca.der(dae, v), (u - a*v^2)/m - 9.81)
    julia> ca.eq(dae, ca.der(dae, m), -u^2)
    julia> ca.x(dae), ca.u(dae), ca.p(dae)
    (["h", "v", "m"], ["u"], ["a"])
    julia> f = ca.create(dae, "f", ["x", "u", "p"], ["ode"])
    Function(f:(x[3],u,p)->(ode[3]) MXFunction)
    julia> f([0, 0, 1], 1, 0.5)
    DM([0, -8.81, -1])
    \enddoctest
    \endlang

    \date 2012-2021
    \author Joel Andersson

    \identifier{5c} */
class CASADI_EXPORT DaeBuilder
  : public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<DaeBuilder>) {
 public:

  /// Readable name of the class
  std::string type_name() const {return "DaeBuilder";}

  /// Default constructor: a null instance, not usable until assigned
  DaeBuilder();

  /** \brief Construct a DaeBuilder instance

      Without \p path, the instance starts empty and the model is built symbolically with
      `add` and `eq`. With \p path, the model description is imported immediately:

      - a directory containing `modelDescription.xml` (an unzipped FMU),
      - a `.fmu` archive (requires a build with `WITH_GHC_FILESYSTEM=ON`; it is unzipped to
        a temporary directory that is removed once nothing references it any more),
      - a `.bmo` Base Modelica file (requires the `lacemodelica` ModelicaParser plugin).

      An imported standard FMU is not symbolic (`symbolic()` is false): equations are
      evaluated by calling the FMU binary, so `create` is available but `oracle`, `eq`,
      `eliminate` etc. are not.

      \param name Instance name; used as default function name prefix by `create()`,
        `transition()` and `export_fmu`
      \param path Model to import, see above; empty for a symbolic model
      \param opts Options (any other key is an error):
        - `"debug"` (bool, false): print debug information during import and enable FMU
          logging
        - `"fmutol"` (double, 0): tolerance passed to the FMU
        - `"ignore_time"` (bool, false): do not expose the hidden time variable of an
          imported model that declares no independent variable
        - `"detect_quad"` (bool, false): classify differential states that do not appear in
          any right-hand side as quadrature states `q`
        - `"resource_serialize_mode"` (string): how an imported FMU is serialized with
          functions created from it, `"link"` or `"embed"`
        - `"enable_ls_dae"` (bool, true): honour the DAE layered standard of an FMU
        - `"enable_ls_serialization"` (bool, true): honour the serialization layered
          standard of an FMU

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("rocket")
      >>> dae.name(), dae.symbolic(), dae.all()
      ('rocket', True, [])
      >>> try:
      ...   DaeBuilder("m", "/nonexistent/path")
      ... except RuntimeError:
      ...   print("import failed")
      import failed
      >>> try:
      ...   DaeBuilder("m", "", dict(foo=True))
      ... except Exception as e:
      ...   print("No such option" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('rocket');
      >> dae.name(), isempty(dae.all())
      rocket
      1
      >> try, DaeBuilder('m', '/nonexistent/path'); failed = false; catch e, failed = true; end
      >> failed
      1
      >> try, DaeBuilder('m', '', struct('foo', true)); catch e, msg = e.message; end
      >> contains(msg, 'No such option')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("rocket");
      julia> ca.name(dae), ca.symbolic(dae), ca.all(dae)
      ("rocket", true, String[])
      julia> try
                 ca.DaeBuilder("m", "/nonexistent/path")
             catch e
                 println("import failed")
             end
      import failed
      julia> try
                 ca.DaeBuilder("m", "", Dict("foo" => true))
             catch e
                 println(occursin("No such option", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa load_fmi_description, symbolic

      \identifier{2ko} */
  explicit DaeBuilder(const std::string& name, const std::string& path = "",
    const Dict& opts = Dict());

  /** \brief Name of instance

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("rocket")
      >>> dae.name()
      'rocket'
      >>> dae.create()
      Function(rocket_dae:(t[0],x[0],z[0],p[0],u[0])->(ode[0],alg[0],quad[0],zero[0]) MXFunction)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('rocket');
      >> dae.name()
      rocket
      >> dae.create()
      rocket_dae:(t[0],x[0],z[0],p[0],u[0])->(ode[0],alg[0],quad[0],zero[0]) MXFunction
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("rocket");
      julia> ca.name(dae)
      "rocket"
      julia> ca.create(dae)
      Function(rocket_dae:(t[0],x[0],z[0],p[0],u[0])->(ode[0],alg[0],quad[0],zero[0]) MXFunction)
      \enddoctest
      \endlang

      \identifier{5d} */
  const std::string& name() const;

  /** @name Variables and equations */
  ///@{
  /** \brief Expression for independent variable (usually time)

      Only available when an independent variable has been declared explicitly with
      `add(name, "independent")` or imported from the model description; otherwise an error
      is raised. Without an explicit declaration a hidden variable named `"time"` exists in
      `all()` but is not categorized as `t`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> t = dae.add("t", "independent")
      >>> x = dae.add("x")
      >>> dae.time(), dae.has_t(), dae.t_new()
      (MX(t), True, ['t'])
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> dae.all(), dae.has_t(), dae.t_new()
      (['time', 'x'], False, [])
      >>> try:
      ...   dae.time()
      ... except Exception as e:
      ...   print("No explicit time variable" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); t = dae.add('t', 'independent'); x = dae.add('x');
      >> dae.time(), dae.has_t(), dae.t_new()
      t
      1
      t
      >> dae = DaeBuilder('m'); x = dae.add('x');
      >> dae.all(), dae.has_t(), isempty(dae.t_new())
      time
      x
      0
      1
      >> try, dae.time(); catch e, disp(contains(e.message, 'No explicit time variable')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> t = ca.add(dae, "t", "independent");
      julia> x = ca.add(dae, "x");
      julia> dae.time(), ca.has_t(dae), ca.t_new(dae)
      (t, true, ["t"])
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> ca.all(dae), ca.has_t(dae), ca.t_new(dae)
      (["time", "x"], false, String[])
      julia> try
                 dae.time()
             catch e
                 println(occursin("No explicit time variable", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa has_t, t_new

      \identifier{2by} */
  const MX& time() const;

  /** \brief Independent variable (usually time)

      Name of the independent variable as a list: empty when none has been declared,
      otherwise a single name. Equivalent to `all("t")`.

      \sa time, has_t

      \identifier{2bz} */
  std::vector<std::string> t_new() const {return all("t");}

  /** \brief Differential states

      Names of the variables of category `x`, in category order. Equivalent to `all("x")`.
      A continuous variable becomes a differential state once its time derivative is given
      a defining equation with `eq`; it is an algebraic variable `z` before that.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> dae.x(), dae.z()
      ([], ['x', 'y'])
      >>> dae.eq(dae.der(x), -y)
      >>> dae.x(), dae.z(), dae.nx()
      (['x'], ['y'], 1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> isempty(dae.x()), dae.z()
      1
      x
      y
      >> dae.eq(dae.der(x), -y);
      >> dae.x(), dae.z(), dae.nx()
      x
      y
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> ca.x(dae), ca.z(dae)
      (String[], ["x", "y"])
      julia> ca.eq(dae, ca.der(dae, x), -y)
      julia> ca.x(dae), ca.z(dae), ca.nx(dae)
      (["x"], ["y"], 1)
      \enddoctest
      \endlang

      \sa z, q, nx, all, inputs

      \identifier{5f} */
  std::vector<std::string> x() const {return all("x");}

  /** \brief Outputs

      Names of the variables with causality `"output"`, in category order. Equivalent to
      `all("y")`. Outputs are also members of another category (`x`, `w`, `z`, ...);
      `outputs("y")` gives their expressions.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", "output")
      >>> y = dae.add("y", "output")
      >>> w = dae.add("w")
      >>> dae.eq(dae.der(x), -w)
      >>> dae.eq(w, x**2)
      >>> dae.eq(y, 2*x)
      >>> dae.y(), dae.x(), dae.w()
      (['x', 'y'], ['x'], ['der(x)', 'w', 'y'])
      >>> dae.outputs("y")
      [MX(x), MX(y)]
      >>> dae.set_causality("y", "local")
      >>> dae.y(), dae.ny()
      (['x'], 1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', 'output'); y = dae.add('y', 'output');
      >> w = dae.add('w');
      >> dae.eq(dae.der(x), -w); dae.eq(w, x^2); dae.eq(y, 2*x);
      >> dae.y()
      x
      y
      >> dae.x(), dae.w()
      x
      der(x)
      w
      y
      >> cellfun(@disp, dae.outputs('y'))
      x
      y
      >> dae.set_causality('y', 'local');
      >> dae.y(), dae.ny()
      x
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", "output");
      julia> y = ca.add(dae, "y", "output");
      julia> w = ca.add(dae, "w");
      julia> ca.eq(dae, ca.der(dae, x), -w)
      julia> ca.eq(dae, w, x^2)
      julia> ca.eq(dae, y, 2*x)
      julia> ca.y(dae), ca.x(dae), ca.w(dae)
      (["x", "y"], ["x"], ["der(x)", "w", "y"])
      julia> ca.outputs(dae, "y")
      2-element Vector{MX}:
       MX(x)
       MX(y)
      julia> ca.set_causality(dae, "y", "local")
      julia> ca.y(dae), ca.ny(dae)
      (["x"], 1)
      \enddoctest
      \endlang

      \sa set_causality, set_all, ny, outputs

      \identifier{2fu} */
  std::vector<std::string> y() const {return all("y");}

  /** \brief Algebraic variables

      Names of the variables of category `z`, in category order. Equivalent to `all("z")`.
      Every continuous `local` or `output` variable starts in this category and leaves it
      when it gets a defining equation (`w`) or an equation for its derivative (`x`).

      \sa x, w, nz, eq

      \identifier{2fv} */
  std::vector<std::string> z() const {return all("z");}

  /** \brief Quadrature states

      Names of the variables of category `q`: differential states that do not appear in any
      right-hand side. Equivalent to `all("q")`. A state is only placed here by
      `set_category(name, "q")` or automatically with the constructor option
      `"detect_quad"`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> q = dae.add("q")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.eq(dae.der(q), x**2)
      >>> dae.x(), dae.q()
      (['x', 'q'], [])
      >>> dae.set_category("q", "q")
      >>> dae.x(), dae.q(), dae.outputs("quad")
      (['x'], ['q'], [MX(sq(x))])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); q = dae.add('q');
      >> dae.eq(dae.der(x), -x); dae.eq(dae.der(q), x^2);
      >> dae.x(), isempty(dae.q())
      x
      q
      1
      >> dae.set_category('q', 'q');
      >> disp(dae.x()), disp(dae.q()), cellfun(@disp, dae.outputs('quad'))
      x
      q
      sq(x)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> q = ca.add(dae, "q");
      julia> ca.eq(dae, ca.der(dae, x), -x)
      julia> ca.eq(dae, ca.der(dae, q), x^2)
      julia> ca.x(dae), ca.q(dae)
      (["x", "q"], String[])
      julia> ca.set_category(dae, "q", "q")
      julia> ca.x(dae), ca.q(dae), ca.outputs(dae, "quad")
      (["x"], ["q"], MX[sq(x)])
      \enddoctest
      \endlang

      \sa set_category, x, nq

      \identifier{2fw} */
  std::vector<std::string> q() const {return all("q");}

  /** \brief Input expressions for a specific category

      The symbolic expressions of the variables of an input category, in category order.

      \param cat One of `"t"`, `"c"`, `"p"`, `"d"`, `"w"`, `"u"`, `"x"`, `"z"`, `"q"`;
        an output category such as `"ode"` is an error
      \return List of (column vector) expressions, one per variable; empty when the
        category is empty

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> z = dae.add("z")
      >>> u = dae.add("u", "input")
      >>> dae.eq(dae.der(x), -x + z*u)
      >>> dae.eq(0, x - 2*z)
      >>> dae.inputs("x"), dae.inputs("z"), dae.inputs("u"), dae.inputs("p")
      ([MX(x)], [MX(z)], [MX(u)], [])
      >>> try:
      ...   dae.inputs("ode")
      ... except Exception as e:
      ...   print("not an input category" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); z = dae.add('z'); u = dae.add('u', 'input');
      >> dae.eq(dae.der(x), -x + z*u); dae.eq(0, x - 2*z);
      >> cellfun(@disp, dae.inputs('x')), cellfun(@disp, dae.inputs('z'))
      x
      z
      >> cellfun(@disp, dae.inputs('u')), disp(numel(dae.inputs('p')))
      u
      0
      >> try, dae.inputs('ode'); catch e, disp(contains(e.message, 'not an input category')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> z = ca.add(dae, "z");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.eq(dae, ca.der(dae, x), -x + z*u)
      julia> ca.eq(dae, 0, x - 2*z)
      julia> ca.inputs(dae, "x"), ca.inputs(dae, "z"), ca.inputs(dae, "u"), ca.inputs(dae, "p")
      (MX[x], MX[z], MX[u], MX[])
      julia> try
                 ca.inputs(dae, "ode")
             catch e
                 println(occursin("not an input category", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa outputs, var, all

      \identifier{2fx} */
  std::vector<MX> inputs(const std::string& cat) const;

  /** \brief Output expressions for a specific category

      The defining expressions of an equation category, in the order of the corresponding
      input category: `"ode"` (right-hand sides of the states `x`), `"alg"` (residuals),
      `"quad"` (right-hand sides of `q`), `"zero"` (zero-crossing functions of the when
      equations), `"ddef"` and `"wdef"` (binding expressions of `d` and `w`), `"y"`
      (expressions of the outputs).

      \param cat One of the output categories above; an input category is an error
      \return List of expressions, empty when the category is empty

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> z = dae.add("z")
      >>> u = dae.add("u", "input")
      >>> dae.eq(dae.der(x), -x + z*u)
      >>> dae.eq(0, x - 2*z)
      >>> dae.outputs("ode"), dae.outputs("alg"), dae.outputs("y")
      ([MX(((z*u)-x))], [MX((-(x-(2.*z))))], [])
      >>> try:
      ...   dae.outputs("x")
      ... except Exception as e:
      ...   print("No input category" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); z = dae.add('z'); u = dae.add('u', 'input');
      >> dae.eq(dae.der(x), -x + z*u); dae.eq(0, x - 2*z);
      >> cellfun(@disp, dae.outputs('ode')), cellfun(@disp, dae.outputs('alg'))
      ((z*u)-x)
      (-(x-(2.*z)))
      >> numel(dae.outputs('y'))
      0
      >> try, dae.outputs('x'); catch e, disp(contains(e.message, 'No input category')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> z = ca.add(dae, "z");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.eq(dae, ca.der(dae, x), -x + z*u)
      julia> ca.eq(dae, 0, x - 2*z)
      julia> ca.outputs(dae, "ode"), ca.outputs(dae, "alg"), ca.outputs(dae, "y")
      (MX[((z*u)-x)], MX[(-(x-(2.*z)))], MX[])
      julia> try
                 ca.outputs(dae, "x")
             catch e
                 println(occursin("No input category", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa inputs, wdef, ddef, create

      \identifier{2fy} */
  std::vector<MX> outputs(const std::string& cat) const;

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Replaced with outputs("ode")
  std::vector<MX> ode() const {return outputs("ode");}

  /// [DEPRECATED] Replaced with outputs("alg")
  std::vector<MX> alg() const {return outputs("alg");}

  /// [DEPRECATED] Replaced with outputs("quad")
  std::vector<MX> quad() const {return outputs("quad");}

  /// [DEPRECATED] Replaced with outputs("zero")
  std::vector<MX> zero() const {return outputs("zero");}

  /// [DEPRECATED] Replaced with outputs("y")
  std::vector<MX> ydef() const {return outputs("y");}

  /// [DEPRECATED] Replaced with set_all("y", name)
  void set_y(const std::vector<std::string>& name) {set_all("y", name);}
#endif // WITH_DEPRECATED_FEATURES

  /** \brief Free controls

      Names of the variables of category `u` (causality `"input"`), in category order.
      Equivalent to `all("u")`.

      \sa p, nu, set_category

      \identifier{5n} */
  std::vector<std::string> u() const {return all("u");}

  /** \brief Parameters

      Names of the variables of category `p` (causality `"parameter"`, variability
      `"tunable"`, floating-point type), in category order. Equivalent to `all("p")`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> c = dae.add("c", "parameter", "fixed")
      >>> u = dae.add("u", "input")
      >>> dae.p(), dae.c(), dae.u()
      (['p'], ['c'], ['u'])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); p = dae.add('p', 'parameter', 'tunable');
      >> c = dae.add('c', 'parameter', 'fixed'); u = dae.add('u', 'input');
      >> dae.p(), dae.c(), dae.u()
      p
      c
      u
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> c = ca.add(dae, "c", "parameter", "fixed");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.p(dae), ca.c(dae), ca.u(dae)
      (["p"], ["c"], ["u"])
      \enddoctest
      \endlang

      \sa c, u, np, set_variability

      \identifier{5o} */
  std::vector<std::string> p() const {return all("p");}

  /** \brief Named constants

      Names of the variables of category `c` (fixed parameters and constants), in category
      order. Equivalent to `all("c")`. Constants are never function inputs: `create`
      substitutes their current value, see `set`.

      \sa p, cdef, set, nc

      \identifier{5p} */
  std::vector<std::string> c() const {return all("c");}

  /** \brief Definitions of named constants

      The binding expressions of the constants `c`, one per entry of `c()`. Only constants
      that have a binding equation (imported from a symbolic model description) are
      supported; a constant without one raises an error.

      \sa c, ddef, wdef

      \identifier{5q} */
  std::vector<MX> cdef() const;

  /** \brief Dependent parameters

      Names of the variables of category `d` (causality `"calculatedParameter"`, or
      `"local"` with `"fixed"`/`"tunable"` variability), in category order. Equivalent to
      `all("d")`. They are defined by binding equations and are calculated by
      `dependent_fun`.

      \sa ddef, dependent_fun, nd

      \identifier{5r} */
  std::vector<std::string> d() const {return all("d");}

  /** \brief Definitions of dependent parameters

      The binding expressions of the dependent parameters, one per entry of `d()`.
      Interdependencies are allowed but must be non-cyclic (see `sort`). Equivalent to
      `outputs("ddef")`. A dependent parameter without a binding equation raises an error.

      \sa d, wdef, sort, eliminate

      \identifier{5s} */
  std::vector<MX> ddef() const;

  /** \brief Dependent variables

      Names of the variables of category `w`, in category order. Equivalent to `all("w")`.
      A variable becomes dependent when it is the left-hand side of an explicit equation
      `eq(v, expr)`; the derivative variables `der(x)` of the states are dependent variables
      too.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w = dae.add("w")
      >>> dae.eq(w, sin(x))
      >>> dae.eq(dae.der(x), -w)
      >>> dae.w()
      ['w', 'der(x)']
      >>> dae.wdef()
      [MX(sin(x)), MX((-w))]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w = dae.add('w');
      >> dae.eq(w, sin(x)); dae.eq(dae.der(x), -w);
      >> dae.w()
      w
      der(x)
      >> cellfun(@disp, dae.wdef())
      sin(x)
      (-w)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w = ca.add(dae, "w");
      julia> ca.eq(dae, w, sin(x))
      julia> ca.eq(dae, ca.der(dae, x), -w)
      julia> ca.w(dae)
      2-element Vector{String}:
       "w"
       "der(x)"
      julia> ca.wdef(dae)
      2-element Vector{MX}:
       MX(sin(x))
       MX((-w))
      \enddoctest
      \endlang

      \sa wdef, dependent_fun, eliminate, nw

      \identifier{5t} */
  std::vector<std::string> w() const {return all("w");}

  /** \brief Dependent variables and corresponding definitions

      The binding expressions of the dependent variables, one per entry of `w()`.
      Interdependencies are allowed but must be non-cyclic (see `sort`). Equivalent to
      `outputs("wdef")`.

      \sa w, ddef, sort, eliminate

      \identifier{5u} */
  std::vector<MX> wdef() const;

  /** \brief Initial conditions, left-hand-side

      The variable expressions that have an initial equation, in the order the equations
      were set with `set_init`; `init_rhs()` holds the matching right-hand sides.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.set_init("x", 2*p)
      >>> dae.init_lhs(), dae.init_rhs()
      ([MX(x)], [MX((2.*p))])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); p = dae.add('p', 'parameter', 'tunable');
      >> dae.eq(dae.der(x), -x); dae.set_init('x', 2*p);
      >> cellfun(@disp, dae.init_lhs()), cellfun(@disp, dae.init_rhs())
      x
      (2.*p)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.set_init(dae, "x", 2*p)
      julia> ca.init_lhs(dae), ca.init_rhs(dae)
      (MX[x], MX[(2.*p)])
      \enddoctest
      \endlang

      \sa init_rhs, set_init

      \identifier{2b1} */
  std::vector<MX> init_lhs() const;

  /** \brief Initial conditions, right-hand-side

      The right-hand sides of the initial equations, matching `init_lhs()`.

      \sa init_lhs, set_init

      \identifier{2b2} */
  std::vector<MX> init_rhs() const;

#ifdef WITH_DEPRECATED_FEATURES
  /// [DEPRECATED] Renamed "y"
  std::vector<std::string> outputs() const {return y();}

  /// [DEPRECATED] Renamed "der"
  std::vector<std::string> derivatives() const {return der();}
#endif // WITH_DEPRECATED_FEATURES

  /** \brief Model structure: initial unknowns

      Names of the variables listed under `InitialUnknowns` in the model structure of an
      imported model description. Empty for a symbolically constructed model.

      \sa load_fmi_description, der

      \identifier{63} */
  std::vector<std::string> initial_unknowns() const;

  /** @name Variables and equations */
  ///@{

  /** \brief Is there a time variable?

      True when an independent variable has been declared explicitly (or imported), i.e.
      when `t_new()` is non-empty.

      \sa time, t_new

      \identifier{64} */
  bool has_t() const;

  /** \brief Differential states

      Number of variables of category `x`, i.e. the length of `x()`. Vector-valued
      variables count once.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> pos = dae.add("pos", dict(dimension=[3]))
      >>> dae.eq(dae.der(pos), -pos)
      >>> dae.nx(), dae.create()
      (1, Function(m_dae:(t[0],x[3],z[0],p[0],u[0])->(ode[3],alg[0],quad[0],zero[0]) MXFunction))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); pos = dae.add('pos', struct('dimension', 3));
      >> dae.eq(dae.der(pos), -pos);
      >> dae.nx(), dae.create()
      1
      m_dae:(t[0],x[3],z[0],p[0],u[0])->(ode[3],alg[0],quad[0],zero[0]) MXFunction
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> pos = ca.add(dae, "pos", Dict("dimension" => [3]));
      julia> ca.eq(dae, ca.der(dae, pos), -pos)
      julia> ca.nx(dae), ca.create(dae)
      (1, m_dae:(t[0],x[3],z[0],p[0],u[0])->(ode[3],alg[0],quad[0],zero[0]) MXFunction)
      \enddoctest
      \endlang

      \sa x

      \identifier{65} */
  casadi_int nx() const;

  /** \brief Algebraic variables

      Number of variables of category `z`, i.e. the length of `z()`.

      \sa z

      \identifier{66} */
  casadi_int nz() const;

  /** \brief Quadrature states

      Number of variables of category `q`, i.e. the length of `q()`.

      \sa q

      \identifier{67} */
  casadi_int nq() const;

  /** \brief Zero-crossing functions

      Number of event indicators, one per `when` equation; the length of
      `outputs("zero")`.

      \sa when, transition

      \identifier{2cb} */
  casadi_int nzero() const;

  /** \brief Output variables

      Number of output variables, i.e. the length of `y()`.

      \sa y

      \identifier{68} */
  casadi_int ny() const;

  /** \brief Free controls

      Number of variables of category `u`, i.e. the length of `u()`.

      \sa u

      \identifier{69} */
  casadi_int nu() const;

  /** \brief Parameters

      Number of variables of category `p`, i.e. the length of `p()`.

      \sa p

      \identifier{6a} */
  casadi_int np() const;

  /** \brief Named constants

      Number of variables of category `c`, i.e. the length of `c()`.

      \sa c

      \identifier{6b} */
  casadi_int nc() const;

  /** \brief Dependent parameters

      Number of variables of category `d`, i.e. the length of `d()`.

      \sa d

      \identifier{6c} */
  casadi_int nd() const;

  /** \brief Dependent variables

      Number of variables of category `w`, i.e. the length of `w()`.

      \sa w

      \identifier{6d} */
  casadi_int nw() const;
  ///@}

  /** @name Symbolic modeling
   *  Formulate a dynamic system model
   */
  ///@{

  /** \brief Add a new model variable

      Creates a variable with FMI 3.0 causality and variability and returns its symbolic
      expression, a dense column vector with `numel` elements (the product of the
      dimensions). The category follows from causality and variability and may change later
      when equations are added:

      - `"parameter"`: `"tunable"` gives `p` (`c` for non-floating-point types),
        `"fixed"` gives `c`; other variabilities are an error
      - `"calculatedParameter"`: `d`; requires `"fixed"` or `"tunable"`
      - `"input"`: `u`; requires `"continuous"` or `"discrete"`
      - `"output"` and `"local"`: `"constant"` gives `c`, `"discrete"` gives a state `x`
        with zero derivative, `"continuous"` gives `z` until an equation defines the
        variable (`w`) or its derivative (`x`); `"local"` also accepts `"fixed"` and
        `"tunable"` (`d`)
      - `"independent"`: `t`; requires `"continuous"`, at most one per model

      Variables with causality `"output"` are additionally listed in `y()`. Unless the first
      variable is the independent one, a hidden variable `"time"` is created first.

      \param name Unique variable name; reusing a name is an error
      \param causality `"parameter"`, `"calculatedParameter"`, `"input"`, `"output"`,
        `"local"` or `"independent"`
      \param variability `"constant"`, `"fixed"`, `"tunable"`, `"discrete"` or
        `"continuous"`
      \param opts Attributes (any other key is an error):
        - `"dimension"` (list of int, `[1]`): shape; all entries must be positive
        - `"type"` (string, `"Float64"`): FMI 3 type, e.g. `"Float32"`, `"Int32"`,
          `"Boolean"`, `"String"`
        - `"start"` (double or list): start value, scalar is broadcast; also sets the
          current value, see `get`
        - `"min"`, `"max"`, `"nominal"` (double): bounds and scaling
        - `"initial"` (string): `"exact"`, `"approx"`, `"calculated"` or `"na"`;
          default per FMI from causality and variability
        - `"description"`, `"unit"`, `"display_unit"` (string): meta data
      \return Symbolic expression of the new variable, named like the variable

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> x, x.shape, dae.category("x"), dae.causality("x"), dae.variability("x")
      (MX(x), (1, 1), 'z', 'local', 'continuous')
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.category("p"), dae.variability("p")
      ('p', 'tunable')
      >>> b = dae.add("b", "input", dict(type="Boolean", start=False))
      >>> dae.category("b"), dae.variability("b"), dae.type("b")
      ('u', 'discrete', 'Boolean')
      >>> r = dae.add("r", dict(dimension=[3], start=[1, 2, 3], min=0, unit="m"))
      >>> r.shape, dae.start("r"), dae.min("r"), dae.unit("r")
      ((3, 1), [1.0, 2.0, 3.0], 0.0, 'm')
      >>> try:
      ...   dae.add("x")
      ... except Exception as e:
      ...   print("already exists" in str(e))
      True
      >>> try:
      ...   dae.add("q", "parameter", "continuous")
      ... except Exception as e:
      ...   print("requires 'fixed' or 'tunable'" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x');
      >> x, size(x), dae.category('x'), dae.causality('x')
      x
      1 1
      z
      local
      >> dae.variability('x')
      continuous
      >> p = dae.add('p', 'parameter', 'tunable');
      >> dae.category('p'), dae.variability('p')
      p
      tunable
      >> b = dae.add('b', 'input', struct('type', 'Boolean', 'start', false));
      >> dae.category('b'), dae.variability('b'), dae.type('b')
      u
      discrete
      Boolean
      >> r = dae.add('r', struct('dimension', 3, 'start', [1 2 3], 'min', 0, 'unit', 'm'));
      >> size(r), dae.start('r'), dae.unit('r')
      3 1
      1 2 3
      m
      >> try, dae.add('x'); catch e, disp(contains(e.message, 'already exists')); end
      1
      >> try, dae.add('q', 'parameter', 'continuous'); catch e, msg = e.message; end
      >> contains(msg, 'requires ''fixed'' or ''tunable''')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> x, size(x), ca.category(dae, "x"), ca.causality(dae, "x"), ca.variability(dae, "x")
      (x, (1, 1), "z", "local", "continuous")
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.category(dae, "p"), ca.variability(dae, "p")
      ("p", "tunable")
      julia> b = ca.add(dae, "b", "input", Dict("type" => "Boolean", "start" => false));
      julia> ca.category(dae, "b"), ca.variability(dae, "b"), ca.type(dae, "b")
      ("u", "discrete", "Boolean")
      julia> opts = Dict("dimension" => [3], "start" => [1, 2, 3], "min" => 0, "unit" => "m");
      julia> r = ca.add(dae, "r", opts);
      julia> size(r), ca.start(dae, "r"), dae.min("r"), ca.unit(dae, "r")
      ((3, 1), [1.0, 2.0, 3.0], 0.0, "m")
      julia> try
                 ca.add(dae, "x")
             catch e
                 println(occursin("already exists", sprint(showerror, e)))
             end
      true
      julia> try
                 ca.add(dae, "q", "parameter", "continuous")
             catch e
                 println(occursin("requires 'fixed' or 'tunable'", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa eq, der, var, category, set_category

      \identifier{2kp} */
  MX add(const std::string& name,
    const std::string& causality,
    const std::string& variability,
    const Dict& opts=Dict());

  /** \brief Add a new model variable, default variability

      The variability follows the FMI 3.0 default: `"fixed"` for parameters and
      calculated parameters, `"continuous"` for floating-point variables and `"discrete"`
      for any other type (see the `"type"` option).

      \identifier{2kq} */
  MX add(const std::string& name,
    const std::string& causality,
    const Dict& opts=Dict());

  /** \brief Add a new model variable, default variability and causality

      The causality defaults to `"local"`; the variability follows from the type as above,
      so a plain `add(name)` gives a continuous local variable of category `z`.

      \identifier{2kr} */
  MX add(const std::string& name,
    const Dict& opts=Dict());

  /** \brief Add a new model variable, symbolic expression already available

      Registers an existing symbolic primitive instead of creating a new one. \p expr must
      be a pure symbol (`is_symbolic()`) whose name equals \p name; anything else is an
      error. The dimension attribute is not taken from the expression.

      \identifier{2ks} */
  void add(const std::string& name,
    const std::string& causality,
    const std::string& variability,
    const MX& expr,
    const Dict& opts=Dict());

  /** \brief Add a simple equation

      Adds the equation `lhs == rhs` and reclassifies variables accordingly:

      - `lhs` is the derivative of a state (`der(x)`): the equation becomes the ODE of `x`,
        which moves from `z` to `x`; `der(x)` becomes a dependent variable `w`
      - `lhs` is a variable of category `z`: the equation is its binding equation and the
        variable moves to `w`; it must not have a binding equation already
      - `lhs` is a concatenation of such variables: split into one equation per component
      - `lhs` is anything else, typically `0`: an implicit (algebraic) equation
        `0 == lhs - rhs` is added to `outputs("alg")`; a second equation for a variable
        that already has a binding equation is treated this way too

      A derivative `der(x)` that appears in `rhs` before it was defined is replaced by a
      new algebraic variable `der_x` with the trivial ODE `der(x) == der_x`, keeping the
      DAE semi-explicit. Scalar sides are broadcast to the length of the other side;
      sparse sides are densified. Variables appearing in `rhs` are marked as used, so a
      quadrature state `q` referenced there becomes a regular state `x`.

      \param lhs Column vector expression; a parameter, constant or state as left-hand
        side is an error
      \param rhs Column vector expression in model variables only; dimensions must match
        after broadcasting
      \param opts Must be empty

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> z = dae.add("z")
      >>> w = dae.add("w")
      >>> dae.x(), dae.z(), dae.w()
      ([], ['x', 'y', 'z', 'w'], [])
      >>> dae.eq(dae.der(x), -x + z)
      >>> dae.x(), dae.z(), dae.w()
      (['x'], ['y', 'z', 'w'], ['der(x)'])
      >>> dae.eq(0, x - 2*z)
      >>> dae.z(), dae.outputs("alg")
      (['y', 'z', 'w'], [MX((-(x-(2.*z))))])
      >>> dae.eq(w, x*z)
      >>> dae.w(), dae.wdef()
      (['der(x)', 'w'], [MX((z-x)), MX((x*z))])
      >>> dae.eq(vertcat(dae.der(y), w), vertcat(w, 1))
      >>> dae.x(), dae.outputs("ode"), dae.outputs("alg")
      (['x', 'y'], [MX((z-x)), MX(w)], [MX((-(x-(2.*z)))), MX((-(w-1)))])
      >>> try:
      ...   dae.eq(x, 1)
      ... except Exception as e:
      ...   print("Cannot handle left-hand-side" in str(e))
      True
      >>> try:
      ...   dae.eq(vertcat(dae.der(x), dae.der(y)), vertcat(1, 2, 3))
      ... except Exception as e:
      ...   print("Mismatched dimensions" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> z = dae.add('z'); w = dae.add('w');
      >> isempty(dae.x()), dae.z(), isempty(dae.w())
      1
      x
      y
      z
      w
      1
      >> dae.eq(dae.der(x), -x + z);
      >> dae.x(), dae.z(), dae.w()
      x
      y
      z
      w
      der(x)
      >> dae.eq(0, x - 2*z);
      >> disp(dae.z()), cellfun(@disp, dae.outputs('alg'))
      y
      z
      w
      (-(x-(2.*z)))
      >> dae.eq(w, x*z);
      >> disp(dae.w()), cellfun(@disp, dae.wdef())
      der(x)
      w
      (z-x)
      (x*z)
      >> dae.eq([dae.der(y); w], [w; 1]);
      >> disp(dae.x()), cellfun(@disp, dae.outputs('ode')), cellfun(@disp, dae.outputs('alg'))
      x
      y
      (z-x)
      w
      (-(x-(2.*z)))
      (-(w-1))
      >> try, dae.eq(x, 1); catch e, disp(contains(e.message, 'Cannot handle left-hand-side')); end
      1
      >> try, dae.eq([dae.der(x); dae.der(y)], [1; 2; 3]); catch e, msg = e.message; end
      >> contains(msg, 'Mismatched dimensions')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> z = ca.add(dae, "z");
      julia> w = ca.add(dae, "w");
      julia> ca.x(dae), ca.z(dae), ca.w(dae)
      (String[], ["x", "y", "z", "w"], String[])
      julia> ca.eq(dae, ca.der(dae, x), -x + z)
      julia> ca.x(dae), ca.z(dae), ca.w(dae)
      (["x"], ["y", "z", "w"], ["der(x)"])
      julia> ca.eq(dae, 0, x - 2*z)
      julia> ca.z(dae), ca.outputs(dae, "alg")
      (["y", "z", "w"], MX[(-(x-(2.*z)))])
      julia> ca.eq(dae, w, x*z)
      julia> ca.w(dae), ca.wdef(dae)
      (["der(x)", "w"], MX[(z-x), (x*z)])
      julia> ca.eq(dae, [ca.der(dae, y); w], [w; 1])
      julia> ca.x(dae), ca.outputs(dae, "ode"), ca.outputs(dae, "alg")
      (["x", "y"], MX[(z-x), w], MX[(-(x-(2.*z))), (-(w-1))])
      julia> try
                 ca.eq(dae, x, 1)
             catch e
                 println(occursin("Cannot handle left-hand-side", sprint(showerror, e)))
             end
      true
      julia> try
                 ca.eq(dae, [ca.der(dae, x); ca.der(dae, y)], [1, 2, 3])
             catch e
                 println(occursin("Mismatched dimensions", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa der, when, set_init, add

      \identifier{2kt} */
  void eq(const MX& lhs, const MX& rhs, const Dict& opts=Dict());

  /** \brief Add when equations

      Adds an event: whenever the zero-crossing function of \p cond becomes negative, the
      assignments in \p eqs are executed. The condition must be a strict inequality
      `a < b` (or `a > b`); it is stored as the zero-crossing function `a - b` in
      `outputs("zero")`. The event transition can be obtained with `transition`.

      \param cond Strict inequality between model expressions; `<=` or any other
        expression is an error
      \param eqs Names returned by `reinit` or `assign`; may be empty
      \param opts Must be empty

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("ball")
      >>> h = dae.add("h", dict(start=5))
      >>> v = dae.add("v", dict(start=0))
      >>> dae.eq(dae.der(h), v)
      >>> dae.eq(dae.der(v), -9.81)
      >>> dae.when(h < 0, [dae.reinit("v", -0.8*v)])
      >>> dae.nzero(), dae.outputs("zero")
      (1, [MX(h)])
      >>> T = dae.transition()
      >>> T
      Function(ball_transition:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      >>> T(0, 0, [-0.1, -3], [], [], [])
      (DM([-0.1, 2.4]), DM(0x1))
      >>> try:
      ...   dae.when(h <= 0, [])
      ... except Exception as e:
      ...   print("Only strict inequality" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('ball');
      >> h = dae.add('h', struct('start', 5)); v = dae.add('v', struct('start', 0));
      >> dae.eq(dae.der(h), v); dae.eq(dae.der(v), -9.81);
      >> dae.when(h < 0, {dae.reinit('v', -0.8*v)});
      >> disp(dae.nzero()), cellfun(@disp, dae.outputs('zero'))
      1
      h
      >> T = dae.transition(); disp(T)
      ball_transition:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction
      >> [post_x, post_z] = T(0, 0, [-0.1; -3], [], [], []); disp(post_x), disp(post_z)
      [-0.1, 2.4]
      0x1
      >> try, dae.when(h <= 0, {}); catch e, msg = e.message; end
      >> contains(msg, 'Only strict inequality')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("ball");
      julia> h = ca.add(dae, "h", Dict("start" => 5));
      julia> v = ca.add(dae, "v", Dict("start" => 0));
      julia> ca.eq(dae, ca.der(dae, h), v)
      julia> ca.eq(dae, ca.der(dae, v), -9.81)
      julia> ca.when(dae, h < 0, [ca.reinit(dae, "v", -0.8*v)])
      julia> ca.nzero(dae), ca.outputs(dae, "zero")
      (1, MX[h])
      julia> T = ca.transition(dae)
      Function(ball_transition:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      julia> T(0, 0, [-0.1, -3], [], [], [])
      ([-0.1, 2.4], 0x1)
      julia> try
                 ca.when(dae, h <= 0, String[])
             catch e
                 println(occursin("Only strict inequality", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa reinit, assign, transition, nzero

      \identifier{2ku} */
  void when(const MX& cond, const std::vector<std::string>& eqs, const Dict& opts=Dict());

  /** \brief Assignment inside a when-equation or if-else equation

      Creates a helper variable holding the expression \p val assigned to the variable
      \p name when the event fires, and returns the helper's name for use in `when`.
      Intended for discrete-valued variables; `reinit` is the counterpart for continuous
      states.

      \param name Model variable that receives the value
      \param val Expression in model variables
      \return Name of the created assignment variable, e.g. `__assign__b__`

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(start=0))
      >>> b = dae.add("b", dict(type="Boolean", start=False))
      >>> dae.eq(dae.der(x), 1)
      >>> dae.assign("b", True)
      '__assign__b__'
      >>> dae.when(x > 1, [dae.assign("b", True)])
      >>> dae.disp(True)
      nt = 0, nx = 2, nz = 0, nq = 0, ny = 0, np = 0, nc = 0, nd = 0, nw = 1, nu = 0
      Model variables:
        x = [x, b]
        w = [der(x)]
      List of dependent variables (w):
        der(x) := 1
      Time derivatives of differential states (x):
        x: der(x)
        b: 0
      When equations:
        when (1-x) < 0 :
          b := 1
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('start', 0));
      >> b = dae.add('b', struct('type', 'Boolean', 'start', false));
      >> dae.eq(dae.der(x), 1);
      >> dae.assign('b', true)
      __assign__b__
      >> dae.when(x > 1, {dae.assign('b', true)});
      >> dae.disp(true)
      nt = 0, nx = 2, nz = 0, nq = 0, ny = 0, np = 0, nc = 0, nd = 0, nw = 1, nu = 0
      Model variables:
        x = [x, b]
        w = [der(x)]
      List of dependent variables (w):
        der(x) := 1
      Time derivatives of differential states (x):
        x: der(x)
        b: 0
      When equations:
        when (1-x) < 0 :
          b := 1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("start" => 0));
      julia> b = ca.add(dae, "b", Dict("type" => "Boolean", "start" => false));
      julia> ca.eq(dae, ca.der(dae, x), 1)
      julia> ca.assign(dae, "b", true)
      "__assign__b__"
      julia> ca.when(dae, x > 1, [ca.assign(dae, "b", true)])
      julia> ca.disp(dae, true)
      nt = 0, nx = 2, nz = 0, nq = 0, ny = 0, np = 0, nc = 0, nd = 0, nw = 1, nu = 0
      Model variables:
        x = [x, b]
        w = [der(x)]
      List of dependent variables (w):
        der(x) := 1
      Time derivatives of differential states (x):
        x: der(x)
        b: 0
      When equations:
        when (1-x) < 0 :
          b := 1
      \enddoctest
      \endlang

      \sa when, reinit

      \identifier{2kv} */
  std::string assign(const std::string& name, const MX& val);

  /** \brief Reinitialize a state inside when-equations

      Creates a helper variable holding the new value \p val of the state \p name after the
      event, and returns the helper's name for use in `when`. Use `pre(v)` for the value
      before the event (currently an identity mapping).

      \param name Model variable that is reinitialized
      \param val Expression in model variables
      \return Name of the created variable, e.g. `__reinit__v__0`

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("ball")
      >>> h = dae.add("h", dict(start=5))
      >>> v = dae.add("v", dict(start=0))
      >>> dae.eq(dae.der(h), v)
      >>> dae.eq(dae.der(v), -9.81)
      >>> r = dae.reinit("v", -0.8*v)
      >>> r, dae.category(r)
      ('__reinit__v__0', 'calculated')
      >>> dae.when(h < 0, [r])
      >>> dae.nzero()
      1
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('ball');
      >> h = dae.add('h', struct('start', 5)); v = dae.add('v', struct('start', 0));
      >> dae.eq(dae.der(h), v); dae.eq(dae.der(v), -9.81);
      >> r = dae.reinit('v', -0.8*v); disp(r), disp(dae.category(r))
      __reinit__v__0
      calculated
      >> dae.when(h < 0, {r}); disp(dae.nzero())
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("ball");
      julia> h = ca.add(dae, "h", Dict("start" => 5));
      julia> v = ca.add(dae, "v", Dict("start" => 0));
      julia> ca.eq(dae, ca.der(dae, h), v)
      julia> ca.eq(dae, ca.der(dae, v), -9.81)
      julia> r = ca.reinit(dae, "v", -0.8*v);
      julia> r, ca.category(dae, r)
      ("__reinit__v__0", "calculated")
      julia> ca.when(dae, h < 0, [r])
      julia> ca.nzero(dae)
      1
      \enddoctest
      \endlang

      \sa when, assign, pre

      \identifier{2kw} */
  std::string reinit(const std::string& name, const MX& val);

  /** \brief Specify the initial equation for a variable

      Sets `name == init_rhs` as initial equation, replacing any previous one for the same
      variable. An empty \p init_rhs removes the initial equation.

      \param name Model variable
      \param init_rhs Expression in model variables, or an empty expression to remove

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.init_lhs(), dae.init_rhs()
      ([], [])
      >>> dae.set_init("x", 2*p)
      >>> dae.init_lhs(), dae.init_rhs()
      ([MX(x)], [MX((2.*p))])
      >>> dae.set_init("x", MX())
      >>> dae.init_lhs(), dae.init_rhs()
      ([], [])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); p = dae.add('p', 'parameter', 'tunable');
      >> dae.eq(dae.der(x), -x);
      >> numel(dae.init_lhs()), numel(dae.init_rhs())
      0
      0
      >> dae.set_init('x', 2*p);
      >> cellfun(@disp, dae.init_lhs()), cellfun(@disp, dae.init_rhs())
      x
      (2.*p)
      >> dae.set_init('x', MX());
      >> numel(dae.init_lhs()), numel(dae.init_rhs())
      0
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.init_lhs(dae), ca.init_rhs(dae)
      (MX[], MX[])
      julia> ca.set_init(dae, "x", 2*p)
      julia> ca.init_lhs(dae), ca.init_rhs(dae)
      (MX[x], MX[(2.*p)])
      julia> ca.set_init(dae, "x", ca.MX())
      julia> ca.init_lhs(dae), ca.init_rhs(dae)
      (MX[], MX[])
      \enddoctest
      \endlang

      \sa init_lhs, init_rhs

      \identifier{2kx} */
  void set_init(const std::string& name, const MX& init_rhs);

  /** \brief Check if dimensions match

      Raises an error if there is more than one independent variable or if it is not
      scalar. Called by `disp` with `more` set and by `add_lc`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> dae.sanity_check()
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m');
      >> dae.sanity_check()
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> ca.sanity_check(dae)
      \enddoctest
      \endlang

      \identifier{2ky} */
  void sanity_check() const;
  ///@}

  /** \brief Reorder variables in a category

      Replaces the ordering of a category by \p v, which must be a permutation of the
      current members: a missing or extra name is an error. The ordering determines the
      order of the entries in the vectors used by `create` (e.g. of `x` and `ode`).

      \param cat Category name, e.g. `"x"`
      \param v All current members of the category, in the new order

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> dae.eq(dae.der(x), -y)
      >>> dae.eq(dae.der(y), x)
      >>> dae.x(), dae.outputs("ode")
      (['x', 'y'], [MX((-y)), MX(x)])
      >>> dae.reorder("x", ["y", "x"])
      >>> dae.x(), dae.outputs("ode")
      (['y', 'x'], [MX(x), MX((-y))])
      >>> try:
      ...   dae.reorder("x", ["y"])
      ... except Exception as e:
      ...   print("1 elements provided for 2 components" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> dae.eq(dae.der(x), -y); dae.eq(dae.der(y), x);
      >> disp(dae.x()), cellfun(@disp, dae.outputs('ode'))
      x
      y
      (-y)
      x
      >> dae.reorder('x', {'y', 'x'});
      >> disp(dae.x()), cellfun(@disp, dae.outputs('ode'))
      y
      x
      x
      (-y)
      >> try, dae.reorder('x', {'y'}); catch e, msg = e.message; end
      >> contains(msg, '1 elements provided for 2 components')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> ca.eq(dae, ca.der(dae, x), -y)
      julia> ca.eq(dae, ca.der(dae, y), x)
      julia> ca.x(dae), ca.outputs(dae, "ode")
      (["x", "y"], MX[(-y), x])
      julia> ca.reorder(dae, "x", ["y", "x"])
      julia> ca.x(dae), ca.outputs(dae, "ode")
      (["y", "x"], MX[x, (-y)])
      julia> try
                 ca.reorder(dae, "x", ["y"])
             catch e
                 println(occursin("1 elements provided for 2 components", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa set_all, all

      \identifier{2kz} */
  void reorder(const std::string& cat, const std::vector<std::string>& v);

  /** \brief Set all variables within a a category

      Makes \p name the complete, ordered content of category \p v, moving variables in and
      out as needed:

      - `"y"`: the listed variables get causality `"output"`, all other outputs become
        `"local"`
      - `"u"`, `"p"`, `"c"`: the listed variables are moved into the category with
        `set_category`; members not listed become constants `c` (for `"u"` and `"p"`)
      - `"x"`, `"q"`: the listed states are moved in; members not listed lose their
        category (they must not appear in any right-hand side)

      Other categories and duplicate names are an error. Names not permitted in the
      category raise the error of `set_category`.

      \param v Category name
      \param name Variables that should form the category, in order

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> u1 = dae.add("u1", "input")
      >>> p1 = dae.add("p1", "parameter", "tunable")
      >>> dae.eq(dae.der(x), -x + u1*p1)
      >>> dae.set_all("u", ["p1", "u1"])
      >>> dae.u(), dae.p(), dae.c()
      (['p1', 'u1'], [], [])
      >>> dae.set_all("p", ["u1"])
      >>> dae.u(), dae.p(), dae.c()
      (['p1'], ['u1'], [])
      >>> dae.set_all("y", ["x"])
      >>> dae.y(), dae.causality("x")
      (['x'], 'output')
      >>> try:
      ...   dae.set_all("x", ["x", "x"])
      ... except Exception as e:
      ...   print("Duplicate names" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); u1 = dae.add('u1', 'input');
      >> p1 = dae.add('p1', 'parameter', 'tunable');
      >> dae.eq(dae.der(x), -x + u1*p1);
      >> dae.set_all('u', {'p1', 'u1'});
      >> dae.u(), isempty(dae.p()), isempty(dae.c())
      p1
      u1
      1
      1
      >> dae.set_all('p', {'u1'});
      >> dae.u(), dae.p()
      p1
      u1
      >> dae.set_all('y', {'x'});
      >> dae.y(), dae.causality('x')
      x
      output
      >> try, dae.set_all('x', {'x', 'x'}); catch e, msg = e.message; end
      >> contains(msg, 'Duplicate names')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> u1 = ca.add(dae, "u1", "input");
      julia> p1 = ca.add(dae, "p1", "parameter", "tunable");
      julia> ca.set_all(dae, "u", ["p1", "u1"])
      julia> ca.u(dae), ca.p(dae), ca.c(dae)
      (["p1", "u1"], String[], String[])
      julia> ca.set_all(dae, "p", ["u1"])
      julia> ca.u(dae), ca.p(dae), ca.c(dae)
      (["p1"], ["u1"], String[])
      julia> ca.set_all(dae, "y", ["x"])
      julia> ca.y(dae), ca.causality(dae, "x")
      (["x"], "output")
      julia> try
                 ca.set_all(dae, "x", ["x", "x"])
             catch e
                 println(occursin("Duplicate names", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa reorder, set_category, set_causality

      \identifier{2l0} */
  void set_all(const std::string& v, const std::vector<std::string>& name);

  /** @name Manipulation
   *  Reformulate the dynamic optimization problem.
   */
  ///@{


  /** \brief Eliminate all dependent parameters

      Substitutes the variables of an acyclic category out of all model expressions:
      `"c"`, `"d"` or `"w"` (sorted first with `sort`). The variables keep existing but
      become `calculated` and the category is left empty. For `"q"` the quadrature states
      are turned into regular states `x` instead. Any other category is an error.

      \param cat `"c"`, `"d"`, `"w"` or `"q"`

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w1 = dae.add("w1")
      >>> w2 = dae.add("w2")
      >>> dae.eq(w2, w1**2)
      >>> dae.eq(w1, 2*x)
      >>> dae.eq(dae.der(x), -w2)
      >>> dae.eliminate("w")
      >>> dae.w(), dae.outputs("ode"), dae.category("w1")
      ([], [MX((-sq((2.*x))))], 'calculated')
      >>> try:
      ...   dae.eliminate("x")
      ... except Exception as e:
      ...   print("Elimination not supported" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w1 = dae.add('w1'); w2 = dae.add('w2');
      >> dae.eq(w2, w1^2); dae.eq(w1, 2*x); dae.eq(dae.der(x), -w2);
      >> dae.eliminate('w');
      >> disp(isempty(dae.w())), cellfun(@disp, dae.outputs('ode')), disp(dae.category('w1'))
      1
      (-sq((2.*x)))
      calculated
      >> try, dae.eliminate('x'); catch e, msg = e.message; end
      >> contains(msg, 'Elimination not supported')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w1 = ca.add(dae, "w1");
      julia> w2 = ca.add(dae, "w2");
      julia> ca.eq(dae, w2, w1^2)
      julia> ca.eq(dae, w1, 2*x)
      julia> ca.eq(dae, ca.der(dae, x), -w2)
      julia> ca.eliminate(dae, "w")
      julia> ca.w(dae), ca.outputs(dae, "ode"), ca.category(dae, "w1")
      (String[], MX[(-sq((2.*x)))], "calculated")
      julia> try
                 ca.eliminate(dae, "x")
             catch e
                 println(occursin("Elimination not supported", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa sort, dependent_fun, lift

      \identifier{2l1} */
  void eliminate(const std::string& cat);

  /** \brief Sort dependent parameters

      Reorders an acyclic category (`"c"`, `"d"` or `"w"`) so that every variable comes
      after the ones it depends on, using a block triangular (Dulmage-Mendelsohn)
      decomposition of the dependency pattern. Vector-valued variables are treated as a
      unit. Other categories are an error.

      \param cat `"c"`, `"d"` or `"w"`

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w1 = dae.add("w1")
      >>> w2 = dae.add("w2")
      >>> dae.eq(w2, w1**2)
      >>> dae.eq(w1, 2*x)
      >>> dae.eq(dae.der(x), -w2)
      >>> dae.w(), dae.wdef()
      (['w2', 'w1', 'der(x)'], [MX(sq(w1)), MX((2.*x)), MX((-w2))])
      >>> dae.sort("w")
      >>> dae.w(), dae.wdef()
      (['w1', 'w2', 'der(x)'], [MX((2.*x)), MX(sq(w1)), MX((-w2))])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w1 = dae.add('w1'); w2 = dae.add('w2');
      >> dae.eq(w2, w1^2); dae.eq(w1, 2*x); dae.eq(dae.der(x), -w2);
      >> disp(dae.w()), cellfun(@disp, dae.wdef())
      w2
      w1
      der(x)
      sq(w1)
      (2.*x)
      (-w2)
      >> dae.sort('w');
      >> disp(dae.w()), cellfun(@disp, dae.wdef())
      w1
      w2
      der(x)
      (2.*x)
      sq(w1)
      (-w2)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w1 = ca.add(dae, "w1");
      julia> w2 = ca.add(dae, "w2");
      julia> ca.eq(dae, w2, w1^2)
      julia> ca.eq(dae, w1, 2*x)
      julia> ca.eq(dae, ca.der(dae, x), -w2)
      julia> ca.w(dae), ca.wdef(dae)
      (["w2", "w1", "der(x)"], MX[sq(w1), (2.*x), (-w2)])
      julia> ca.sort(dae, "w")
      julia> ca.w(dae), ca.wdef(dae)
      (["w1", "w2", "der(x)"], MX[(2.*x), sq(w1), (-w2)])
      \enddoctest
      \endlang

      \sa eliminate, reorder

      \identifier{2l2} */
  void sort(const std::string& cat);

  /** \brief Lift problem formulation by extracting shared subexpressions

      Scans the right-hand sides of the ODEs, quadratures and algebraic equations and
      replaces shared subexpressions and/or embedded function calls by new dependent
      variables named `w_<k>`, appended to `w` with their definitions. The lifted
      definitions become part of `wdef`; `create` with `"w"` as input then exposes them,
      and `create` with `lifted_calls` accounts for the calls in `jac_wdef_w` and Hessians.
      A warning is issued when `w` is non-empty already, which is the case as soon as a
      state has an ODE.

      \param lift_shared Lift subexpressions that are used more than once
      \param lift_calls Lift calls to embedded Function objects, including their arguments

      \lang{python}
      \doctest
      >>> a = MX.sym("a")
      >>> F = Function("F", [a], [sin(a) + a**2])
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> dae.eq(0, x - F(y))
      >>> dae.eq(0, y - F(x)**2 - F(y))
      >>> dae.lift()
      >>> dae.w()
      ['w_0', 'w_1', 'w_2', 'w_3', 'w_4']
      >>> dae.wdef()
      [MX(x), MX(y), MX(F(w_1){0}), MX(F(w_0){0}), MX(F(w_1){0})]
      >>> dae.outputs("alg")
      [MX((-(w_0-w_2))), MX((-((w_1-sq(w_3))-w_4)))]
      >>> dae.create("f", ["z", "w"], ["alg", "wdef"])
      Function(f:(z[2],w[5])->(alg[2],wdef[5]) MXFunction)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> a = MX.sym('a'); F = Function('F', {a}, {sin(a) + a^2});
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> dae.eq(0, x - F(y)); dae.eq(0, y - F(x)^2 - F(y));
      >> dae.lift();
      >> dae.w()
      w_0
      w_1
      w_2
      w_3
      w_4
      >> cellfun(@disp, dae.wdef())
      x
      y
      F(w_1){0}
      F(w_0){0}
      F(w_1){0}
      >> cellfun(@disp, dae.outputs('alg'))
      (-(w_0-w_2))
      (-((w_1-sq(w_3))-w_4))
      >> dae.create('f', {'z', 'w'}, {'alg', 'wdef'})
      f:(z[2],w[5])->(alg[2],wdef[5]) MXFunction
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> a = ca.MX.sym("a");
      julia> F = ca.Function("F", [a], [sin(a) + a^2]);
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> ca.eq(dae, 0, x - F(y))
      julia> ca.eq(dae, 0, y - F(x)^2 - F(y))
      julia> ca.lift(dae)
      julia> ca.w(dae)
      5-element Vector{String}:
       "w_0"
       "w_1"
       "w_2"
       "w_3"
       "w_4"
      julia> ca.wdef(dae)
      5-element Vector{MX}:
       MX(x)
       MX(y)
       MX(F(w_1){0})
       MX(F(w_0){0})
       MX(F(w_1){0})
      julia> ca.outputs(dae, "alg")
      2-element Vector{MX}:
       MX((-(w_0-w_2)))
       MX((-((w_1-sq(w_3))-w_4)))
      julia> ca.create(dae, "f", ["z", "w"], ["alg", "wdef"])
      Function(f:(z[2],w[5])->(alg[2],wdef[5]) MXFunction)
      \enddoctest
      \endlang

      \sa eliminate, w, wdef

      \identifier{2l3} */
  void lift(bool lift_shared = true, bool lift_calls = true);

  /** \brief Prune unused controls

      Intended to drop the parameters `p` and/or controls `u` that do not appear in any
      model equation from their categories (the variables themselves remain).

      \param prune_p Consider the parameters `p`
      \param prune_u Consider the controls `u`

      \sa p, u

      \identifier{2l4} */
  void prune(bool prune_p = true, bool prune_u = true);

  /** \brief Identify iteration variables and residual equations using naming convention

      Variables named `res__<iv>` are interpreted as residual equations with iteration
      variable `<iv>` (optionally followed by `__<hold>` markers). The iteration variables
      are removed from `z` and `u` and appended to `z`, so that `z` and `alg` match up.
      Hold markers are not implemented and only produce a warning.

      \sa z, sort

      \identifier{2l5} */
  void tear();
  ///@}

  /** @name Functions
   *  Add or load auxiliary functions
   */
  ///@{

  /** \brief Add a function from loaded expressions

      Creates a Function mapping the model variables \p arg to the binding expressions of
      the dependent variables \p res and stores it in the instance. Stored functions are
      listed by `disp` and retrieved with `fun`.

      \param name Function name; must not exist yet
      \param arg Names of model variables used as inputs
      \param res Names of dependent variables `w`, whose definitions become the outputs;
        any other name is an error
      \param opts Function options
      \return The created function

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w = dae.add("w")
      >>> dae.eq(w, sin(x))
      >>> dae.eq(dae.der(x), -w)
      >>> f = dae.add_fun("wfun", ["x"], ["w"])
      >>> f
      Function(wfun:(x)->(w) MXFunction)
      >>> dae.has_fun("wfun"), dae.has_fun("nope")
      (True, False)
      >>> try:
      ...   dae.add_fun("bad", ["x"], ["x"])
      ... except Exception as e:
      ...   print("Cannot find dependent" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w = dae.add('w');
      >> dae.eq(w, sin(x)); dae.eq(dae.der(x), -w);
      >> f = dae.add_fun('wfun', {'x'}, {'w'}); disp(f)
      wfun:(x)->(w) MXFunction
      >> dae.has_fun('wfun'), dae.has_fun('nope')
      1
      0
      >> try, dae.add_fun('bad', {'x'}, {'x'}); catch e, msg = e.message; end
      >> contains(msg, 'Cannot find dependent')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w = ca.add(dae, "w");
      julia> ca.eq(dae, w, sin(x))
      julia> f = ca.add_fun(dae, "wfun", ["x"], ["w"])
      Function(wfun:(x)->(w) MXFunction)
      julia> ca.has_fun(dae, "wfun"), ca.has_fun(dae, "nope")
      (true, false)
      julia> try
                 ca.add_fun(dae, "bad", ["x"], ["x"])
             catch e
                 println(occursin("Cannot find dependent", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa fun, has_fun, gather_fun

      \identifier{2l6} */
  Function add_fun(const std::string& name,
                   const std::vector<std::string>& arg,
                   const std::vector<std::string>& res, const Dict& opts=Dict());

  /** \brief Add an already existing function

      Stores \p f under its own name; a stored function with the same name is an error.

      \lang{python}
      \doctest
      >>> a = MX.sym("a")
      >>> F = Function("F", [a], [a**2])
      >>> dae = DaeBuilder("m")
      >>> dae.add_fun(F)
      Function(F:(i0)->(o0) MXFunction)
      >>> [g.name() for g in dae.fun()]
      ['F']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> a = MX.sym('a'); F = Function('F', {a}, {a^2});
      >> dae = DaeBuilder('m');
      >> dae.add_fun(F)
      F:(i0)->(o0) MXFunction
      >> cellfun(@(g) disp(g.name()), dae.fun())
      F
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> a = ca.MX.sym("a");
      julia> F = ca.Function("F", [a], [a^2]);
      julia> dae = ca.DaeBuilder("m");
      julia> ca.add_fun(dae, F)
      Function(F:(i0)->(o0) MXFunction)
      julia> [ca.name(g) for g in ca.fun(dae)]
      1-element Vector{String}:
       "F"
      \enddoctest
      \endlang

      \identifier{2l7} */
  Function add_fun(const Function& f);

  /** \brief Add an external function

      Loads the function \p name from the compiled sources or shared library behind
      \p compiler (see `external`) and stores it.

      \identifier{2l8} */
  Function add_fun(const std::string& name, const Importer& compiler,
                   const Dict& opts=Dict());

  /** \brief Does a particular function already exist?

      True when a function with this name has been stored with `add_fun` or `gather_fun`.

      \sa fun, add_fun

      \identifier{2l9} */
  bool has_fun(const std::string& name) const;

  /** \brief Get function by name

      Returns the stored function with this name; an unknown name is an error.

      \lang{python}
      \doctest
      >>> a = MX.sym("a")
      >>> dae = DaeBuilder("m")
      >>> f = dae.add_fun(Function("F", [a], [a**2]))
      >>> dae.fun("F")
      Function(F:(i0)->(o0) MXFunction)
      >>> try:
      ...   dae.fun("nope")
      ... except Exception as e:
      ...   print("No such function" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> a = MX.sym('a'); dae = DaeBuilder('m');
      >> f = dae.add_fun(Function('F', {a}, {a^2}));
      >> dae.fun('F')
      F:(i0)->(o0) MXFunction
      >> try, dae.fun('nope'); catch e, disp(contains(e.message, 'No such function')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> a = ca.MX.sym("a");
      julia> dae = ca.DaeBuilder("m");
      julia> f = ca.add_fun(dae, ca.Function("F", [a], [a^2]));
      julia> ca.fun(dae, "F")
      Function(F:(i0)->(o0) MXFunction)
      julia> try
                 ca.fun(dae, "nope")
             catch e
                 println(occursin("No such function", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa has_fun, add_fun

      \identifier{2la} */
  Function fun(const std::string& name) const;

  /// Get all functions
  std::vector<Function> fun() const;

  /** \brief Collect embedded functions from the expression graph

      Finds the Function objects called inside the model equations and stores them, so
      that they are listed by `fun()`. Functions whose name is already stored are skipped
      with a warning.

      \param max_depth How deep to look into nested function calls; `-1` for unlimited

      \lang{python}
      \doctest
      >>> a = MX.sym("a")
      >>> F = Function("F", [a], [a**2])
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> dae.eq(dae.der(x), F(x))
      >>> [g.name() for g in dae.fun()]
      []
      >>> dae.gather_fun()
      >>> [g.name() for g in dae.fun()]
      ['F']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> a = MX.sym('a'); F = Function('F', {a}, {a^2});
      >> dae = DaeBuilder('m'); x = dae.add('x');
      >> dae.eq(dae.der(x), F(x));
      >> numel(dae.fun())
      0
      >> dae.gather_fun();
      >> cellfun(@(g) disp(g.name()), dae.fun())
      F
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> a = ca.MX.sym("a");
      julia> F = ca.Function("F", [a], [a^2]);
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> ca.eq(dae, ca.der(dae, x), F(x))
      julia> [ca.name(g) for g in ca.fun(dae)]
      String[]
      julia> ca.gather_fun(dae)
      julia> [ca.name(g) for g in ca.fun(dae)]
      1-element Vector{String}:
       "F"
      \enddoctest
      \endlang

      \sa fun, add_fun

      \identifier{2lb} */
  void gather_fun(casadi_int max_depth = -1);
///@}

  /** @name Import and export
   */
  ///@{
  /// Import existing problem from FMI/XML, alias of load_fmi_description
  void parse_fmi(const std::string& filename) {load_fmi_description(filename); }

  /** \brief Does the FMU provide support for analytic derivatives

      Reads the `providesDirectionalDerivatives` flag of an imported FMU. Calling it on a
      symbolic instance is an error.

      \sa symbolic, load_fmi_description

      \identifier{2lc} */
  bool provides_directional_derivatives() const;

  /// Does the FMU provide support for analytic derivatives (FMI 2 naming)
  bool provides_directional_derivative() const {return provides_directional_derivatives();}

  /** \brief Import problem description from FMI or XML

      Reads a `modelDescription.xml` (FMI 1.0, 2.0 or 3.0) into an empty instance: model
      variables with their attributes, the default experiment (`start_time` etc.), the
      model structure (`der`, `initial_unknowns`, dependencies used by `jac_sparsity`) and,
      for symbolic model descriptions, the equations. A standard FMU results in a
      non-symbolic instance whose functions call the FMU binary found in the FMU
      directory.
      The constructor with a path calls this for `<path>/modelDescription.xml`. An instance
      that already has variables, or a missing file, is an error.

      \param filename Path to the XML file

      \sa DaeBuilder, symbolic, provides_directional_derivatives

      \identifier{2ld} */
  void load_fmi_description(const std::string& filename);

  /** \brief Export instance into an FMU

      Generates the files of an FMI 3.0 model exchange FMU in the current directory: C code
      for the model equations and their derivatives (`<name>.c`, `<name>.h`), the FMI
      wrapper `<name>_wrap.c`, `buildDescription.xml`, `modelDescription.xml` and, with
      serialization, the oracle as `oracle.casadi` plus a layered-standard manifest.
      Quadrature states are not supported. The result maps each local file name to its
      path inside the FMU archive, ready for `compile_fmu` and `pack_fmu`.

      \param opts Options (any other key is an error):
        - `"no_warning"` (bool, false): suppress the "experimental" warning
        - `"with_serialization"` (bool, true): include the serialized oracle
      \return Map local file name to archive path

      \sa compile_fmu, pack_fmu, oracle

      \identifier{2le} */
  Dict export_fmu(const Dict& opts=Dict());

  /** \brief Compile the sources produced by export_fmu

      Compiles all `.c` entries of \p files (through one amalgamation translation unit)
      into a shared library `<name>.<ext>` in the current directory using an Importer
      plugin, and returns the map extended with the binary at its FMI 3 archive location
      `binaries/<platform>/`. The FMI 3 headers must be on the include path.

      \param files the {local_file -> archive path} map returned by export_fmu
      \param opts compile options: `"compiler"` (Importer plugin, default `"shell"`),
        `"compiler_options"` (options for that plugin), `"include_dirs"` (list of
        directories)
      \return the file map augmented with the amalgamation source and compiled binary

      \sa export_fmu, pack_fmu

      \identifier{2iq} */
  Dict compile_fmu(const Dict& files, const Dict& opts=Dict());

  /** \brief Pack files from export_fmu / compile_fmu into a single .fmu archive

      Zips the entries of \p files at their archive paths and removes the local files
      afterwards. Requires a build with `WITH_LIBZIP=ON`.

      \param files the {local_file -> archive path} map to pack
      \param opts packaging options: path (default '<name>.fmu')
      \return Path of the created archive

      \sa export_fmu, compile_fmu

      \identifier{2ir} */
  std::string pack_fmu(const Dict& files, const Dict& opts=Dict());

  /** \brief Add a named linear combination of output expressions

      Declares the name of a linear combination of output categories, so that `create` can
      produce its Hessian blocks `hess_<name>_<in1>_<in2>` with the multipliers
      `lam_<out>` as extra inputs.

      \param name Alphanumeric name
      \param f_out Output categories in the combination, e.g. `["ode", "alg"]`; must be
        non-empty without duplicates

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> z = dae.add("z")
      >>> dae.eq(dae.der(x), -x*z)
      >>> dae.eq(0, z - x**2)
      >>> dae.add_lc("gamma", ["ode", "alg"])
      >>> H = dae.create("H", ["x", "z", "lam_ode", "lam_alg"],
      ...                ["hess_gamma_x_x", "hess_gamma_x_z"])
      >>> H
      Function(H:(x,z,lam_ode,lam_alg)->(hess_gamma_x_x,hess_gamma_x_z) MXFunction)
      >>> H(1, 2, 1, 1)
      (DM(2), DM(-1))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); z = dae.add('z');
      >> dae.eq(dae.der(x), -x*z); dae.eq(0, z - x^2);
      >> dae.add_lc('gamma', {'ode', 'alg'});
      >> out = {'hess_gamma_x_x', 'hess_gamma_x_z'};
      >> H = dae.create('H', {'x', 'z', 'lam_ode', 'lam_alg'}, out); disp(H)
      H:(x,z,lam_ode,lam_alg)->(hess_gamma_x_x,hess_gamma_x_z) MXFunction
      >> [Hxx, Hxz] = H(1, 2, 1, 1); disp(Hxx), disp(Hxz)
      2
      -1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> z = ca.add(dae, "z");
      julia> ca.eq(dae, ca.der(dae, x), -x*z)
      julia> ca.eq(dae, 0, z - x^2)
      julia> ca.add_lc(dae, "gamma", ["ode", "alg"])
      julia> H = ca.create(dae, "H", ["x", "z", "lam_ode", "lam_alg"],
                           ["hess_gamma_x_x", "hess_gamma_x_z"])
      Function(H:(x,z,lam_ode,lam_alg)->(hess_gamma_x_x,hess_gamma_x_z) MXFunction)
      julia> H(1, 2, 1, 1)
      (2, -1)
      \enddoctest
      \endlang

      \sa create

      \identifier{2lf} */
  void add_lc(const std::string& name, const std::vector<std::string>& f_out);

  /** \brief Construct a function object, legacy syntax

      Same as the overload with options, but with the flags `sx` (expand to an SX function)
      and `lifted_calls` (account for function calls lifted into `w`, see `lift`) instead;
      neither is available for an imported FMU.

      \identifier{2lg} */
  Function create(const std::string& fname,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out, bool sx, bool lifted_calls = false) const;

  /** \brief  Construct a function object, names provided

      The function factory: builds a Function whose inputs are variable categories and
      whose outputs are equation categories or derivatives of them. Inputs are any of
      `t`, `p`, `d`, `u`, `x`, `z`, `q`, `w`; constants `c` are never inputs, their current
      values (see `set`) are substituted. Outputs are `ode`, `alg`, `quad`, `zero`, `ddef`,
      `wdef`, `y`. Dependent variables `w` and parameters `d` not given as inputs are
      substituted out. Derivative outputs use the Function factory naming with
      underscores: `jac_ode_x` (Jacobian block) or `hess_<lc>_x_x` for a linear
      combination declared with `add_lc`, with the multipliers `lam_ode`, `lam_alg`, ...
      as inputs. The created function is a snapshot: later changes to the instance do not
      affect it. For an imported FMU the function evaluates the FMU binary.

      \param name    Name assigned to the resulting function object
      \param name_in   Names of all the inputs
      \param name_out  Names of all the outputs
      \param opts    Optional settings: Function options for a symbolic model; for an FMU
        also `"aux"` (names of variables reported in `stats()` after evaluation),
        `"enable_ad"`, `"parallelization"` and the other FmuFunction options

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> z = dae.add("z")
      >>> u = dae.add("u", "input")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> c = dae.add("c", "parameter", "fixed", dict(start=3))
      >>> w = dae.add("w")
      >>> y = dae.add("y", "output")
      >>> dae.eq(dae.der(x), -c*x + z*u)
      >>> dae.eq(0, x - p*z**2)
      >>> dae.eq(w, x + z)
      >>> dae.eq(y, w**2)
      >>> f = dae.create("f", ["x", "z", "u", "p"], ["ode", "alg", "y"])
      >>> f
      Function(f:(x,z,u,p)->(ode,alg,y) MXFunction)
      >>> f(1, 2, 3, 4)
      (DM(3), DM(15), DM(9))
      >>> J = dae.create("J", ["x", "z", "u", "p"], ["jac_ode_x", "jac_alg_z"])
      >>> J(1, 2, 3, 4)
      (DM(-3), DM(16))
      >>> dae.create("g", ["x", "z", "w"], ["ode", "wdef"])
      Function(g:(x,z,w[3])->(ode,wdef[3]) MXFunction)
      >>> dae.create("s", ["x", "z", "u", "p"], ["ode"], True)
      Function(s:(x,z,u,p)->(ode) SXFunction)
      >>> try:
      ...   dae.create("f", ["x", "z", "c"], ["ode"])
      ... except Exception as e:
      ...   print('Cannot process "c" as input' in str(e))
      True
      >>> try:
      ...   dae.create("f", ["x", "z"], ["w"])
      ... except Exception as e:
      ...   print('Cannot process "w" as output' in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); z = dae.add('z'); u = dae.add('u', 'input');
      >> p = dae.add('p', 'parameter', 'tunable');
      >> c = dae.add('c', 'parameter', 'fixed', struct('start', 3));
      >> w = dae.add('w'); y = dae.add('y', 'output');
      >> dae.eq(dae.der(x), -c*x + z*u); dae.eq(0, x - p*z^2);
      >> dae.eq(w, x + z); dae.eq(y, w^2);
      >> f = dae.create('f', {'x', 'z', 'u', 'p'}, {'ode', 'alg', 'y'}); disp(f)
      f:(x,z,u,p)->(ode,alg,y) MXFunction
      >> [ode, alg, yval] = f(1, 2, 3, 4); disp(ode), disp(alg), disp(yval)
      3
      15
      9
      >> J = dae.create('J', {'x', 'z', 'u', 'p'}, {'jac_ode_x', 'jac_alg_z'});
      >> [J_ode_x, J_alg_z] = J(1, 2, 3, 4); disp(J_ode_x), disp(J_alg_z)
      -3
      16
      >> dae.create('g', {'x', 'z', 'w'}, {'ode', 'wdef'})
      g:(x,z,w[3])->(ode,wdef[3]) MXFunction
      >> dae.create('s', {'x', 'z', 'u', 'p'}, {'ode'}, true)
      s:(x,z,u,p)->(ode) SXFunction
      >> try, dae.create('f', {'x', 'z', 'c'}, {'ode'}); catch e, msg = e.message; end
      >> contains(msg, 'Cannot process "c" as input')
      1
      >> try, dae.create('f', {'x', 'z'}, {'w'}); catch e, msg = e.message; end
      >> contains(msg, 'Cannot process "w" as output')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> z = ca.add(dae, "z");
      julia> u = ca.add(dae, "u", "input");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> c = ca.add(dae, "c", "parameter", "fixed", Dict("start" => 3));
      julia> w = ca.add(dae, "w");
      julia> y = ca.add(dae, "y", "output");
      julia> ca.eq(dae, ca.der(dae, x), -c*x + z*u)
      julia> ca.eq(dae, 0, x - p*z^2)
      julia> ca.eq(dae, w, x + z)
      julia> ca.eq(dae, y, w^2)
      julia> f = ca.create(dae, "f", ["x", "z", "u", "p"], ["ode", "alg", "y"])
      Function(f:(x,z,u,p)->(ode,alg,y) MXFunction)
      julia> f(1, 2, 3, 4)
      (3, 15, 9)
      julia> J = ca.create(dae, "J", ["x", "z", "u", "p"], ["jac_ode_x", "jac_alg_z"]);
      julia> J(1, 2, 3, 4)
      (-3, 16)
      julia> ca.create(dae, "g", ["x", "z", "w"], ["ode", "wdef"])
      Function(g:(x,z,w[3])->(ode,wdef[3]) MXFunction)
      julia> ca.create(dae, "s", ["x", "z", "u", "p"], ["ode"], true)
      Function(s:(x,z,u,p)->(ode) SXFunction)
      julia> try
                 ca.create(dae, "f", ["x", "z", "c"], ["ode"])
             catch e
                 println(occursin("Cannot process \"c\" as input", sprint(showerror, e)))
             end
      true
      julia> try
                 ca.create(dae, "f", ["x", "z"], ["w"])
             catch e
                 println(occursin("Cannot process \"w\" as output", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa oracle, add_lc, dependent_fun, transition

      \identifier{6e} */
  Function create(const std::string& name,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const Dict& opts=Dict()) const;
  ///@}

  /** \brief Create a function with standard integrator DAE signature

      Shorthand for `create(fname, ["t", "x", "z", "p", "u"], ["ode", "alg", "quad",
      "zero"], opts)`: a function that can be passed directly to `integrator`. Empty
      categories give empty inputs or outputs.

      \param name    Name assigned to the resulting function object
      \param opts    Optional settings

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> u = dae.add("u", "input")
      >>> dae.eq(dae.der(x), -x + u)
      >>> f = dae.create("mydae")
      >>> f
      Function(mydae:(t[0],x,z[0],p[0],u)->(ode,alg[0],quad[0],zero[0]) MXFunction)
      >>> f(0, 1, [], [], 3)
      (DM(2), DM(0x1), DM(0x1), DM(0x1))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); u = dae.add('u', 'input');
      >> dae.eq(dae.der(x), -x + u);
      >> f = dae.create('mydae'); disp(f)
      mydae:(t[0],x,z[0],p[0],u)->(ode,alg[0],quad[0],zero[0]) MXFunction
      >> [ode, alg, quad, zero] = f(0, 1, [], [], 3); disp(ode), disp(alg), disp(quad), disp(zero)
      2
      0x1
      0x1
      0x1
      >> res = f('t', 0, 'x', 1, 'u', 3); disp(res.ode)
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.eq(dae, ca.der(dae, x), -x + u)
      julia> f = ca.create(dae, "mydae")
      Function(mydae:(t[0],x,z[0],p[0],u)->(ode,alg[0],quad[0],zero[0]) MXFunction)
      julia> f(0, 1, [], [], 3)
      (2, 0x1, 0x1, 0x1)
      julia> f(t=0, x=1, u=3)["ode"]
      DM(2)
      \enddoctest
      \endlang

      \sa transition

      \identifier{2c0} */
  Function create(const std::string& fname, const Dict& opts=Dict()) const;

  /** \brief Create a function with standard integrator DAE signature, default naming

      The function is named `<name>_dae`.

      \identifier{2c1} */
  Function create() const {return create(name() + "_dae");}

  /** \brief Construct a function for evaluating dependent parameters

      Builds a Function that calculates the dependent parameters `d` and/or dependent
      variables `w` from their binding equations, substituting the interdependencies out,
      with the given categories as inputs.

      \param fname Function name
      \param s_in Input categories; a category that is also an output is an error, and
        every variable the definitions depend on must be covered
      \param s_out `"d"`, `"w"` or both; anything else is an error
      \return Function with one (concatenated) output per entry of \p s_out

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w1 = dae.add("w1")
      >>> w2 = dae.add("w2")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.eq(w1, 2*x)
      >>> dae.eq(w2, w1**2 + p)
      >>> dae.eq(dae.der(x), -w2)
      >>> f = dae.dependent_fun("f", ["x", "p"], ["w"])
      >>> f
      Function(f:(x,p)->(w[3]) MXFunction)
      >>> f(1, 3)
      DM([2, 7, -7])
      >>> try:
      ...   dae.dependent_fun("f", ["x", "w"], ["w"])
      ... except Exception as e:
      ...   print("cannot be both input and output" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w1 = dae.add('w1'); w2 = dae.add('w2');
      >> p = dae.add('p', 'parameter', 'tunable');
      >> dae.eq(w1, 2*x); dae.eq(w2, w1^2 + p); dae.eq(dae.der(x), -w2);
      >> f = dae.dependent_fun('f', {'x', 'p'}, {'w'}); disp(f)
      f:(x,p)->(w[3]) MXFunction
      >> f(1, 3)
      [2, 7, -7]
      >> try, dae.dependent_fun('f', {'x', 'w'}, {'w'}); catch e, msg = e.message; end
      >> contains(msg, 'cannot be both input and output')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w1 = ca.add(dae, "w1");
      julia> w2 = ca.add(dae, "w2");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.eq(dae, w1, 2*x)
      julia> ca.eq(dae, w2, w1^2 + p)
      julia> ca.eq(dae, ca.der(dae, x), -w2)
      julia> f = ca.dependent_fun(dae, "f", ["x", "p"], ["w"])
      Function(f:(x,p)->(w[3]) MXFunction)
      julia> f(1, 3)
      DM([2, 7, -7])
      julia> try
                 ca.dependent_fun(dae, "f", ["x", "w"], ["w"])
             catch e
                 println(occursin("cannot be both input and output", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa w, d, create, eliminate

      \identifier{2lh} */
  Function dependent_fun(const std::string& fname,
      const std::vector<std::string>& s_in,
      const std::vector<std::string>& s_out) const;

  /** \brief Construct a function describing transition at a specific events

      Builds a Function with inputs `t`, `x`, `z`, `p`, `u` and outputs `post_x`,
      `post_z`: the states after the event with index \p index (0-based, in the order the
      `when` equations were added). Dependent variables are substituted out.

      \param fname Function name
      \param index Event index; out of range is an error

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("ball")
      >>> h = dae.add("h", dict(start=5))
      >>> v = dae.add("v", dict(start=0))
      >>> dae.eq(dae.der(h), v)
      >>> dae.eq(dae.der(v), -9.81)
      >>> dae.when(h < 0, [dae.reinit("v", -0.8*v)])
      >>> T = dae.transition("tr", 0)
      >>> T
      Function(tr:(t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      >>> T(0, [-0.1, -3], [], [], [])
      (DM([-0.1, 2.4]), DM(0x1))
      >>> try:
      ...   dae.transition("tr", 1)
      ... except Exception as e:
      ...   print("Illegal event index" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('ball');
      >> h = dae.add('h', struct('start', 5)); v = dae.add('v', struct('start', 0));
      >> dae.eq(dae.der(h), v); dae.eq(dae.der(v), -9.81);
      >> dae.when(h < 0, {dae.reinit('v', -0.8*v)});
      >> T = dae.transition('tr', 0); disp(T)
      tr:(t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction
      >> [post_x, post_z] = T(0, [-0.1; -3], [], [], []); disp(post_x), disp(post_z)
      [-0.1, 2.4]
      0x1
      >> try, dae.transition('tr', 1); catch e, msg = e.message; end
      >> contains(msg, 'Illegal event index')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("ball");
      julia> h = ca.add(dae, "h", Dict("start" => 5));
      julia> v = ca.add(dae, "v", Dict("start" => 0));
      julia> ca.eq(dae, ca.der(dae, h), v)
      julia> ca.eq(dae, ca.der(dae, v), -9.81)
      julia> ca.when(dae, h < 0, [ca.reinit(dae, "v", -0.8*v)])
      julia> T = ca.transition(dae, "tr", 0)
      Function(tr:(t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      julia> T(0, [-0.1, -3], [], [], [])
      ([-0.1, 2.4], 0x1)
      julia> try
                 ca.transition(dae, "tr", 1)
             catch e
                 println(occursin("Illegal event index", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa when, create

      \identifier{2li} */
  Function transition(const std::string& fname, casadi_int index) const;

  /** \brief Construct a function describing transition at any events

      Builds the event function expected by the `"transition"` option of `integrator`,
      with an event `index` as first input: for a single event the function of
      `transition(fname, 0)` with a dummy index input, for several events a `switch` over
      the per-event functions. Without when equations a null Function is returned.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("ball")
      >>> h = dae.add("h", dict(start=5))
      >>> v = dae.add("v", dict(start=0))
      >>> dae.eq(dae.der(h), v)
      >>> dae.eq(dae.der(v), -9.81)
      >>> dae.transition().is_null()
      True
      >>> dae.when(h < 0, [dae.reinit("v", -0.8*v)])
      >>> dae.transition("tr")
      Function(tr:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('ball');
      >> h = dae.add('h', struct('start', 5)); v = dae.add('v', struct('start', 0));
      >> dae.eq(dae.der(h), v); dae.eq(dae.der(v), -9.81);
      >> dae.transition().is_null()
      1
      >> dae.when(h < 0, {dae.reinit('v', -0.8*v)});
      >> dae.transition('tr')
      tr:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("ball");
      julia> h = ca.add(dae, "h", Dict("start" => 5));
      julia> v = ca.add(dae, "v", Dict("start" => 0));
      julia> ca.eq(dae, ca.der(dae, h), v)
      julia> ca.eq(dae, ca.der(dae, v), -9.81)
      julia> ca.is_null(ca.transition(dae))
      true
      julia> ca.when(dae, h < 0, [ca.reinit(dae, "v", -0.8*v)])
      julia> ca.transition(dae, "tr")
      Function(tr:(index[],t[0],x[2],z[0],p[0],u[0])->(post_x[2],post_z[0]) MXFunction)
      \enddoctest
      \endlang

      \sa when, create

      \identifier{2lj} */
  Function transition(const std::string& fname) const;

  /** \brief Construct an event transition function, default naming

      The function is named `<name>_transition`.

      \identifier{2lk} */
  Function transition() const {return transition(name() + "_transition");}

  ///@{
  /** \brief Get variable expression by name

      The symbolic expression of a model variable: a dense column vector with `numel`
      entries, named after the variable. An unknown name is an error. Also available via
      the call operator.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(dimension=[2]))
      >>> dae.var("x"), dae.var("x").shape
      (MX(x), (2, 1))
      >>> dae.has("x"), dae.has("y")
      (True, False)
      >>> try:
      ...   dae.var("y")
      ... except Exception as e:
      ...   print("No such variable" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('dimension', 2));
      >> dae.var('x'), size(dae.var('x'))
      x
      2 1
      >> dae.has('x'), dae.has('y')
      1
      0
      >> try, dae.var('y'); catch e, disp(contains(e.message, 'No such variable')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("dimension" => [2]));
      julia> ca.var(dae, "x"), size(ca.var(dae, "x"))
      (x, (2, 1))
      julia> ca.has(dae, "x"), ca.has(dae, "y")
      (true, false)
      julia> try
                 ca.var(dae, "y")
             catch e
                 println(occursin("No such variable", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa has, all, inputs

      \identifier{2ll} */
  MX var(const std::string& name) const;
  MX operator()(const std::string& name) const {return var(name);}
  ///@}

  /** \brief Model structure: All time derivatives

      Names of all derivative variables, in creation order: the `der(x)` variables created
      by `der` for a symbolic model, or the variables marked as derivatives in an imported
      model structure.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> dae.der()
      []
      >>> dae.eq(dae.der(x), -y)
      >>> dae.eq(dae.der(y), x)
      >>> dae.der()
      ['der(x)', 'der(y)']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> isempty(dae.der())
      1
      >> dae.eq(dae.der(x), -y); dae.eq(dae.der(y), x);
      >> dae.der()
      der(x)
      der(y)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> ca.der(dae)
      String[]
      julia> ca.eq(dae, ca.der(dae, x), -y)
      julia> ca.eq(dae, ca.der(dae, y), x)
      julia> ca.der(dae)
      2-element Vector{String}:
       "der(x)"
       "der(y)"
      \enddoctest
      \endlang

      \identifier{2fz} */
  std::vector<std::string> der() const;

  /** \brief Get the time derivative of model variables

      The names of the derivative variables of the given variables; a variable without a
      derivative variable is an error.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.eq(dae.der(x), -y)
      >>> dae.eq(dae.der(y), x)
      >>> dae.der("x"), dae.der(["x", "y"])
      ('der(x)', ['der(x)', 'der(y)'])
      >>> try:
      ...   dae.der("p")
      ... except Exception as e:
      ...   print("No derivative expression" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> p = dae.add('p', 'parameter', 'tunable');
      >> dae.eq(dae.der(x), -y); dae.eq(dae.der(y), x);
      >> dae.der('x'), dae.der({'x', 'y'})
      der(x)
      der(x)
      der(y)
      >> try, dae.der('p'); catch e, disp(contains(e.message, 'No derivative expression')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.eq(dae, ca.der(dae, x), -y)
      julia> ca.eq(dae, ca.der(dae, y), x)
      julia> ca.der(dae, "x"), ca.der(dae, ["x", "y"])
      ("der(x)", ["der(x)", "der(y)"])
      julia> try
                 ca.der(dae, "p")
             catch e
                 println(occursin("No derivative expression", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2lm} */
  std::vector<std::string> der(const std::vector<std::string>& name) const;

  ///@{
  /** \brief Differentiate an expression with respect to time

      For a model variable `x` the derivative variable `der(x)` is returned, created on
      first use (as an algebraic variable until `eq` gives it a defining equation). For a
      concatenation of variables the derivatives are concatenated; for a general expression
      the chain rule is applied through the derivatives of the variables it depends on. The
      derivative of the independent variable is 1 and that of parameters, constants,
      inputs and discrete variables is 0. The expression must be a column vector in model
      variables.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> p = dae.add("p", "parameter", "tunable")
      >>> dae.der(x)
      MX(der(x))
      >>> dae.der(vertcat(x, y))
      MX(vertcat(der(x), der(y)))
      >>> dae.der(x*y + p)
      MX(((y*der(x))+(x*der(y))))
      >>> dae.der(p)
      MX(0)
      >>> dae.z()
      ['x', 'y', 'der(x)', 'der(y)']
      >>> try:
      ...   dae.der(vertcat(x, y).T)
      ... except Exception as e:
      ...   print("must be a vector" in str(e))
      True
      >>> try:
      ...   dae.der(MX.sym("foo"))
      ... except Exception as e:
      ...   print("No such variable" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> p = dae.add('p', 'parameter', 'tunable');
      >> dae.der(x)
      der(x)
      >> dae.der([x; y])
      vertcat(der(x), der(y))
      >> dae.der(x*y + p)
      ((y*der(x))+(x*der(y)))
      >> dae.der(p)
      0
      >> dae.z()
      x
      y
      der(x)
      der(y)
      >> try, dae.der([x; y]'); catch e, disp(contains(e.message, 'must be a vector')); end
      1
      >> try, dae.der(MX.sym('foo')); catch e, disp(contains(e.message, 'No such variable')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> ca.der(dae, x)
      MX(der(x))
      julia> ca.der(dae, [x; y])
      MX(vertcat(der(x), der(y)))
      julia> ca.der(dae, x*y + p)
      MX(((y*der(x))+(x*der(y))))
      julia> ca.der(dae, p)
      MX(0)
      julia> ca.z(dae)
      4-element Vector{String}:
       "x"
       "y"
       "der(x)"
       "der(y)"
      julia> try
                 ca.der(dae, [x; y]')
             catch e
                 println(occursin("must be a vector", sprint(showerror, e)))
             end
      true
      julia> try
                 ca.der(dae, ca.MX.sym("foo"))
             catch e
                 println(occursin("No such variable", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa eq, der, var

      \identifier{2ln} */
  MX der(const MX& v) const;
  MX der(const MX& v);
  ///@}

  /** \brief Get the pre-variables of model variables

      Not implemented: returns the names unchanged, with a warning the first time.

      \sa reinit

      \identifier{2lo} */
  std::vector<std::string> pre(const std::vector<std::string>& name) const;

  /** \brief Get the pre-expression given variable expression

      Not implemented: returns the expression unchanged, with a warning the first time.
      Meant for the value of a variable just before an event, as in
      `reinit("v", -0.8*pre(v))`.

      \sa reinit, when

      \identifier{2lp} */
  MX pre(const MX& v) const;

  /** \brief Does a variable have a binding equation?

      A binding equation is the right-hand side of an explicit equation `eq(v, expr)`
      (or an imported one); dependent variables, dependent parameters and defined
      derivatives have one.

      \sa beq, eq

      \identifier{2lq} */
  bool has_beq(const std::string& name) const;

  /** \brief Get the binding equation for a variable

      The right-hand side of the explicit equation defining the variable; a variable
      without binding equation is an error.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w = dae.add("w")
      >>> dae.eq(w, sin(x))
      >>> dae.eq(dae.der(x), -w)
      >>> dae.beq("w"), dae.beq("der(x)")
      (MX(sin(x)), MX((-w)))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w = dae.add('w');
      >> dae.eq(w, sin(x)); dae.eq(dae.der(x), -w);
      >> dae.beq('w'), dae.beq('der(x)')
      sin(x)
      (-w)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w = ca.add(dae, "w");
      julia> ca.eq(dae, w, sin(x))
      julia> ca.beq(dae, "w")
      MX(sin(x))
      \enddoctest
      \endlang

      \sa has_beq, wdef

      \identifier{2lr} */
  MX beq(const std::string& name) const;

  ///@{
  /** \brief Get/set value reference

      The FMI value reference of a variable: imported from the model description, or by
      default the index of the variable in the instance.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y")
      >>> dae.value_reference("x"), dae.value_reference("y")
      (1, 2)
      >>> dae.set_value_reference("x", 100)
      >>> dae.value_reference("x")
      100
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y');
      >> dae.value_reference('x'), dae.value_reference('y')
      1
      2
      >> dae.set_value_reference('x', 100);
      >> dae.value_reference('x')
      100
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y");
      julia> ca.value_reference(dae, "x"), ca.value_reference(dae, "y")
      (1, 2)
      julia> ca.set_value_reference(dae, "x", 100)
      julia> ca.value_reference(dae, "x")
      100
      \enddoctest
      \endlang

      \identifier{2ls} */
  casadi_int value_reference(const std::string& name) const;
  void set_value_reference(const std::string& name, casadi_int val);
  ///@}

  ///@{
  /** \brief Get/set description

      Free-text description of a variable, also settable with the `"description"` option
      of `add`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(description="position"))
      >>> dae.description("x")
      'position'
      >>> dae.set_description("x", "height")
      >>> dae.description("x")
      'height'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('description', 'position'));
      >> dae.description('x')
      position
      >> dae.set_description('x', 'height');
      >> dae.description('x')
      height
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("description" => "position"));
      julia> ca.description(dae, "x")
      "position"
      julia> ca.set_description(dae, "x", "height")
      julia> ca.description(dae, "x")
      "height"
      \enddoctest
      \endlang

      \identifier{2lt} */
  std::string description(const std::string& name) const;
  void set_description(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /** \brief Get/set the type

      The variable type in FMI 3 naming (`"Float64"`, `"Float32"`, `"Int32"`, `"Boolean"`,
      `"String"`, ...) or, with \p fmi_version 2, in FMI 2 naming (`"real"`, `"integer"`,
      `"boolean"`, `"string"`, `"enum"`). Setting accepts FMI 3 names; other versions or
      names are an error.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> n = dae.add("n", dict(type="Int32"))
      >>> dae.type("x"), dae.type("x", 2), dae.type("n"), dae.type("n", 2)
      ('Float64', 'real', 'Int32', 'integer')
      >>> dae.set_type("x", "Float32")
      >>> dae.type("x")
      'Float32'
      >>> try:
      ...   dae.set_type("x", "bogus")
      ... except Exception as e:
      ...   print("No such enum" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); n = dae.add('n', struct('type', 'Int32'));
      >> dae.type('x'), dae.type('x', 2), dae.type('n'), dae.type('n', 2)
      Float64
      real
      Int32
      integer
      >> dae.set_type('x', 'Float32');
      >> dae.type('x')
      Float32
      >> try, dae.set_type('x', 'bogus'); catch e, disp(contains(e.message, 'No such enum')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> n = ca.add(dae, "n", Dict("type" => "Int32"));
      julia> ca.type(dae, "x"), ca.type(dae, "x", 2), ca.type(dae, "n"), ca.type(dae, "n", 2)
      ("Float64", "real", "Int32", "integer")
      julia> ca.set_type(dae, "x", "Float32")
      julia> ca.type(dae, "x")
      "Float32"
      julia> try
                 ca.set_type(dae, "x", "bogus")
             catch e
                 println(occursin("No such enum", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{2lu} */
  std::string type(const std::string& name, casadi_int fmi_version = 3) const;
  void set_type(const std::string& name, const std::string& val);
  ///@}

  /** \brief Get the causality

      One of `"parameter"`, `"calculatedParameter"`, `"input"`, `"output"`, `"local"`,
      `"independent"`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> t = dae.add("t", "independent")
      >>> u = dae.add("u", "input")
      >>> x = dae.add("x")
      >>> dae.causality("t"), dae.causality("u"), dae.causality("x")
      ('independent', 'input', 'local')
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); t = dae.add('t', 'independent'); u = dae.add('u', 'input');
      >> x = dae.add('x');
      >> dae.causality('t'), dae.causality('u'), dae.causality('x')
      independent
      input
      local
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> t = ca.add(dae, "t", "independent");
      julia> u = ca.add(dae, "u", "input");
      julia> x = ca.add(dae, "x");
      julia> ca.causality(dae, "t"), ca.causality(dae, "u"), ca.causality(dae, "x")
      ("independent", "input", "local")
      \enddoctest
      \endlang

      \sa set_causality, variability, category

      \identifier{2lv} */
  std::string causality(const std::string& name) const;

  /** \brief Which categories are possible for a variable?

      The categories `set_category` accepts for the variable given its current category
      and usage; an empty string stands for "no category". Controls, parameters and
      constants can be exchanged freely; states not used in any right-hand side can
      become quadrature states or lose their category.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> u = dae.add("u", "input")
      >>> x = dae.add("x")
      >>> q = dae.add("q")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.eq(dae.der(q), x**2)
      >>> dae.categories("u"), dae.categories("x"), dae.categories("q")
      (['p', 'u', 'c'], ['x'], ['x', 'q', ''])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); u = dae.add('u', 'input'); x = dae.add('x'); q = dae.add('q');
      >> dae.eq(dae.der(x), -x); dae.eq(dae.der(q), x^2);
      >> dae.categories('u'), dae.categories('x')
      p
      u
      c
      x
      >> strjoin(cellstr(dae.categories('q')), ',')
      x,q,
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> u = ca.add(dae, "u", "input");
      julia> x = ca.add(dae, "x");
      julia> q = ca.add(dae, "q");
      julia> ca.eq(dae, ca.der(dae, x), -x)
      julia> ca.eq(dae, ca.der(dae, q), x^2)
      julia> ca.categories(dae, "u"), ca.categories(dae, "x"), ca.categories(dae, "q")
      (["p", "u", "c"], ["x"], ["x", "q", ""])
      \enddoctest
      \endlang

      \sa category, set_category

      \identifier{2lw} */
  std::vector<std::string> categories(const std::string& name) const;

  /** \brief Set the causality, if permitted

      The following changes are permitted:

      - For controls `u` (variability `"continuous"`, causality `"input"`), free
        parameters `p` (variability `"tunable"`, causality `"parameter"`) and fixed
        parameters `c` (variability `"fixed"`, causality `"parameter"`), causality can only
        be changed indirectly, by updating the variability
      - Add or remove an output `y` by setting the causality to `"output"` or `"local"`,
        respectively

      No other changes are permitted.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> y = dae.add("y", "output")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.eq(y, 2*x)
      >>> dae.set_causality("y", "local")
      >>> dae.set_causality("x", "output")
      >>> dae.y(), dae.causality("y")
      (['x'], 'local')
      >>> try:
      ...   dae.set_causality("x", "input")
      ... except Exception as e:
      ...   print("Cannot change causality" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); y = dae.add('y', 'output');
      >> dae.eq(dae.der(x), -x); dae.eq(y, 2*x);
      >> dae.set_causality('y', 'local'); dae.set_causality('x', 'output');
      >> dae.y(), dae.causality('y')
      x
      local
      >> try, dae.set_causality('x', 'input'); catch e, msg = e.message; end
      >> contains(msg, 'Cannot change causality')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> y = ca.add(dae, "y", "output");
      julia> ca.eq(dae, y, 2*x)
      julia> ca.set_causality(dae, "y", "local")
      julia> ca.set_causality(dae, "x", "output")
      julia> ca.y(dae), ca.causality(dae, "y")
      (["x"], "local")
      julia> try
                 ca.set_causality(dae, "x", "input")
             catch e
                 println(occursin("Cannot change causality", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa causality, set_variability, set_category, set_all

      \identifier{2c2} */
  void set_causality(const std::string& name, const std::string& val);

  /** \brief Get the variability

      One of `"constant"`, `"fixed"`, `"tunable"`, `"discrete"`, `"continuous"`.

      \sa set_variability, causality

      \identifier{2lx} */
  std::string variability(const std::string& name) const;

  /** \brief Set the variability, if permitted

      For controls `u` (variability `"continuous"`, causality `"input"`), free parameters
      `p` (variability `"tunable"`, causality `"parameter"`) and fixed parameters `c`
      (variability `"fixed"`, causality `"parameter"`), update variability in order to
      change the category. Causality is updated accordingly.

      Other changes are not permitted.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> u = dae.add("u", "input")
      >>> x = dae.add("x")
      >>> dae.eq(dae.der(x), -u*x)
      >>> dae.set_variability("u", "fixed")
      >>> dae.u(), dae.c(), dae.causality("u")
      ([], ['u'], 'parameter')
      >>> dae.set_variability("u", "tunable")
      >>> dae.p(), dae.causality("u")
      (['u'], 'parameter')
      >>> dae.set_variability("u", "continuous")
      >>> dae.u(), dae.causality("u")
      (['u'], 'input')
      >>> try:
      ...   dae.set_variability("x", "discrete")
      ... except Exception as e:
      ...   print("Cannot change variability" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); u = dae.add('u', 'input'); x = dae.add('x');
      >> dae.eq(dae.der(x), -u*x);
      >> dae.set_variability('u', 'fixed');
      >> isempty(dae.u()), dae.c(), dae.causality('u')
      1
      u
      parameter
      >> dae.set_variability('u', 'tunable');
      >> dae.p(), dae.causality('u')
      u
      parameter
      >> dae.set_variability('u', 'continuous');
      >> dae.u(), dae.causality('u')
      u
      input
      >> try, dae.set_variability('x', 'discrete'); catch e, msg = e.message; end
      >> contains(msg, 'Cannot change variability')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> u = ca.add(dae, "u", "input");
      julia> x = ca.add(dae, "x");
      julia> ca.set_variability(dae, "u", "fixed")
      julia> ca.u(dae), ca.c(dae), ca.causality(dae, "u")
      (String[], ["u"], "parameter")
      julia> ca.set_variability(dae, "u", "tunable")
      julia> ca.p(dae), ca.causality(dae, "u")
      (["u"], "parameter")
      julia> ca.set_variability(dae, "u", "continuous")
      julia> ca.u(dae), ca.causality(dae, "u")
      (["u"], "input")
      julia> try
                 ca.set_variability(dae, "x", "discrete")
             catch e
                 println(occursin("Cannot change variability", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa variability, set_causality, set_category

      \identifier{2c3} */
  void set_variability(const std::string& name, const std::string& val);

  /** \brief Get the variable category

      One of `"t"`, `"c"`, `"p"`, `"d"`, `"w"`, `"u"`, `"x"`, `"z"`, `"q"`,
      `"calculated"` (internal helper variables), or an empty string for no category.

      \sa categories, set_category, causality, variability

      \identifier{2ly} */
  std::string category(const std::string& name) const;

  /** \brief Set the variable category, if permitted

      The following changes are permitted:

      - Controls `u` can be changed to/from tunable parameters `p` or fixed parameters `c`
      - Differential states that do not appear in the right-hand-sides can be changed
        between regular states `x`, quadrature states `q` and no category (empty string)
      - A state can be made the independent variable `t`, replacing the current one

      Other changes are not permitted. Causality and variability is updated accordingly.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> u = dae.add("u", "input")
      >>> x = dae.add("x")
      >>> q = dae.add("q")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.eq(dae.der(q), x**2)
      >>> dae.set_category("q", "q")
      >>> dae.x(), dae.q()
      (['x'], ['q'])
      >>> dae.set_category("q", "")
      >>> dae.x(), dae.q(), dae.category("q")
      (['x'], [], '')
      >>> dae.set_category("u", "p")
      >>> dae.u(), dae.p(), dae.causality("u"), dae.variability("u")
      ([], ['u'], 'parameter', 'tunable')
      >>> try:
      ...   dae.set_category("x", "q")
      ... except Exception as e:
      ...   print("Cannot change category of x from 'x' to 'q'" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); u = dae.add('u', 'input'); x = dae.add('x'); q = dae.add('q');
      >> dae.eq(dae.der(x), -x); dae.eq(dae.der(q), x^2);
      >> dae.set_category('q', 'q');
      >> dae.x(), dae.q()
      x
      q
      >> dae.set_category('q', '');
      >> dae.x(), isempty(dae.q()), isempty(dae.category('q'))
      x
      1
      1
      >> dae.set_category('u', 'p');
      >> isempty(dae.u()), dae.p()
      1
      u
      >> dae.causality('u'), dae.variability('u')
      parameter
      tunable
      >> try, dae.set_category('x', 'q'); catch e, msg = e.message; end
      >> contains(msg, 'Cannot change category of x from ''x'' to ''q''')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> u = ca.add(dae, "u", "input");
      julia> x = ca.add(dae, "x");
      julia> q = ca.add(dae, "q");
      julia> ca.eq(dae, ca.der(dae, x), -x)
      julia> ca.eq(dae, ca.der(dae, q), x^2)
      julia> ca.set_category(dae, "q", "q")
      julia> ca.x(dae), ca.q(dae)
      (["x"], ["q"])
      julia> ca.set_category(dae, "q", "")
      julia> ca.x(dae), ca.q(dae), ca.category(dae, "q")
      (["x"], String[], "")
      julia> ca.set_category(dae, "u", "p")
      julia> ca.u(dae), ca.p(dae), ca.causality(dae, "u"), ca.variability(dae, "u")
      (String[], ["u"], "parameter", "tunable")
      julia> try
                 ca.set_category(dae, "x", "q")
             catch e
                 println(occursin("Cannot change category of x from 'x' to 'q'", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa category, categories, set_all

      \identifier{2c4} */
  void set_category(const std::string& name, const std::string& val);

  ///@{
  /** \brief Get/set the initial property

      The FMI `initial` attribute: `"exact"`, `"approx"`, `"calculated"` or `"na"`. The
      default follows from causality and variability per the FMI specification; it can be
      given with the `"initial"` option of `add`. Other values are an error.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(initial="approx"))
      >>> p = dae.add("p", "parameter", "tunable")
      >>> u = dae.add("u", "input")
      >>> dae.initial("x"), dae.initial("p"), dae.initial("u")
      ('approx', 'exact', 'na')
      >>> dae.set_initial("x", "exact")
      >>> dae.initial("x")
      'exact'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('initial', 'approx'));
      >> p = dae.add('p', 'parameter', 'tunable'); u = dae.add('u', 'input');
      >> dae.initial('x'), dae.initial('p'), dae.initial('u')
      approx
      exact
      na
      >> dae.set_initial('x', 'exact');
      >> dae.initial('x')
      exact
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("initial" => "approx"));
      julia> p = ca.add(dae, "p", "parameter", "tunable");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.initial(dae, "x"), ca.initial(dae, "p"), ca.initial(dae, "u")
      ("approx", "exact", "na")
      julia> dae.set_initial("x", "exact")
      julia> ca.initial(dae, "x")
      "exact"
      \enddoctest
      \endlang

      \identifier{2lz} */
  std::string initial(const std::string& name) const;
  void set_initial(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /** \brief Get/set the unit

      Unit string of a variable, e.g. `"m/s"`; empty by default. Also settable with the
      `"unit"` option of `add`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> v = dae.add("v", dict(unit="m/s"))
      >>> dae.unit("v")
      'm/s'
      >>> dae.set_unit("v", "km/h")
      >>> dae.unit("v")
      'km/h'
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); v = dae.add('v', struct('unit', 'm/s'));
      >> dae.unit('v')
      m/s
      >> dae.set_unit('v', 'km/h');
      >> dae.unit('v')
      km/h
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> v = ca.add(dae, "v", Dict("unit" => "m/s"));
      julia> ca.unit(dae, "v")
      "m/s"
      julia> ca.set_unit(dae, "v", "km/h")
      julia> ca.unit(dae, "v")
      "km/h"
      \enddoctest
      \endlang

      \identifier{2m0} */
  std::string unit(const std::string& name) const;
  void set_unit(const std::string& name, const std::string& val);
  ///@}

  ///@{
  /** \brief Get/set the display unit

      Display unit string of a variable; empty by default. Also settable with the
      `"display_unit"` option of `add`.

      \identifier{2m1} */
  std::string display_unit(const std::string& name) const;
  void set_display_unit(const std::string& name, const std::string& val);
  ///@}

  /** \brief Get the number of elements of a variable

      Product of the dimensions, i.e. the length of the variable's column vector
      expression.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> m = dae.add("m", dict(dimension=[2, 3]))
      >>> dae.numel("x"), dae.numel("m"), dae.dimension("m"), dae.var("m").shape
      (1, 6, [2, 3], (6, 1))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); m = dae.add('m', struct('dimension', [2 3]));
      >> dae.numel('x'), dae.numel('m')
      1
      6
      >> dae.dimension('m'), size(dae.var('m'))
      2 3
      6 1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> m = ca.add(dae, "m", Dict("dimension" => [2, 3]));
      julia> ca.numel(dae, "x"), ca.numel(dae, "m"), ca.dimension(dae, "m"), size(ca.var(dae, "m"))
      (1, 6, [2, 3], (6, 1))
      \enddoctest
      \endlang

      \sa dimension, var

      \identifier{2m2} */
  casadi_int numel(const std::string& name) const;

  /** \brief Get the dimensions of a variable

      The `"dimension"` list given to `add` (default `[1]`). The expression is always
      flattened to a column vector.

      \sa numel

      \identifier{2m3} */
  std::vector<casadi_int> dimension(const std::string& name) const;

  /** \brief Get the start time

      Start time of the default experiment: imported from the model description or set
      with `set_start_time`; nan when unknown.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> dae.start_time(), dae.stop_time(), dae.tolerance(), dae.step_size()
      (nan, nan, nan, nan)
      >>> dae.set_start_time(0)
      >>> dae.set_stop_time(10)
      >>> dae.set_tolerance(1e-6)
      >>> dae.set_step_size(0.01)
      >>> dae.start_time(), dae.stop_time(), dae.tolerance(), dae.step_size()
      (0.0, 10.0, 1e-06, 0.01)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m');
      >> fprintf('%g %g %g %g\n', dae.start_time, dae.stop_time, dae.tolerance, dae.step_size)
      NaN NaN NaN NaN
      >> dae.set_start_time(0); dae.set_stop_time(10);
      >> dae.set_tolerance(1e-6); dae.set_step_size(0.01);
      >> fprintf('%g %g %g %g\n', dae.start_time, dae.stop_time, dae.tolerance, dae.step_size)
      0 10 1e-06 0.01
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> ca.start_time(dae), ca.stop_time(dae), ca.tolerance(dae), ca.step_size(dae)
      (NaN, NaN, NaN, NaN)
      julia> ca.set_start_time(dae, 0)
      julia> ca.set_stop_time(dae, 10)
      julia> ca.set_tolerance(dae, 1e-6)
      julia> ca.set_step_size(dae, 0.01)
      julia> ca.start_time(dae), ca.stop_time(dae), ca.tolerance(dae), ca.step_size(dae)
      (0.0, 10.0, 1.0e-6, 0.01)
      \enddoctest
      \endlang

      \sa set_start_time, stop_time

      \identifier{2m4} */
  double start_time() const;

  /// Set the start time of the default experiment
  void set_start_time(double val);

  /// Get the stop time of the default experiment, nan when unknown
  double stop_time() const;

  /// Set the stop time of the default experiment
  void set_stop_time(double val);

  /// Get the tolerance of the default experiment, nan when unknown
  double tolerance() const;

  /// Set the tolerance of the default experiment
  void set_tolerance(double val);

  /// Get the step size of the default experiment, nan when unknown
  double step_size() const;

  /// Set the step size of the default experiment
  void set_step_size(double val);

  // The following routines are not needed in MATLAB and would cause ambiguity
  // Note that a multirow strings can be interpreted as a vector of strings
#if !(defined(SWIG) && defined(SWIGMATLAB))
  /** \brief Get the time derivative of model variables, single variable

      The name of the derivative variable of a single variable, e.g. `der(x)` for `x`.

      \identifier{2m5} */
  std::string der(const std::string& name) const;

  /// Get the pre-variables of model variables, single variable (not implemented)
  std::string pre(const std::string& name) const;

  /** \brief Get an attribute, single variable

      Numerical attribute of a scalar variable, selected by name: `"min"`, `"max"`,
      `"nominal"`, `"start"` or `"value"` (the current value, see `get`). Any other
      attribute name, or a non-scalar variable, is an error.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(start=1, min=-2, nominal=10))
      >>> v = dae.add("v", dict(dimension=[2], start=[1, 2]))
      >>> dae.attribute("start", "x"), dae.attribute("min", "x")
      (1.0, -2.0)
      >>> dae.attribute("start", ["x", "v"]), dae.attribute("nominal", ["x"])
      ([1.0, 1.0, 2.0], [10.0])
      >>> dae.set_attribute("min", "x", -1)
      >>> dae.set_attribute("start", ["x", "v"], [9, 8, 7])
      >>> dae.min("x"), dae.start("v")
      (-1.0, [8.0, 7.0])
      >>> try:
      ...   dae.attribute("unit", "x")
      ... except Exception as e:
      ...   print("No such enum" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('start', 1, 'min', -2, 'nominal', 10));
      >> v = dae.add('v', struct('dimension', 2, 'start', [1 2]));
      >> dae.attribute('start', 'x'), dae.attribute('min', 'x')
      1
      -2
      >> dae.attribute('start', {'x', 'v'}), dae.attribute('nominal', {'x'})
      1 1 2
      10
      >> dae.set_attribute('min', 'x', -1); dae.set_attribute('start', {'x', 'v'}, [9 8 7]);
      >> dae.min('x'), dae.start('v')
      -1
      8 7
      >> try, dae.attribute('unit', 'x'); catch e, disp(contains(e.message, 'No such enum')); end
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("start" => 1, "min" => -2, "nominal" => 10));
      julia> v = ca.add(dae, "v", Dict("dimension" => [2], "start" => [1, 2]));
      julia> ca.attribute(dae, "start", "x"), ca.attribute(dae, "min", "x")
      (1.0, -2.0)
      julia> ca.attribute(dae, "start", ["x", "v"]), ca.attribute(dae, "nominal", ["x"])
      ([1.0, 1.0, 2.0], [10.0])
      julia> ca.set_attribute(dae, "min", "x", -1)
      julia> ca.set_attribute(dae, "start", ["x", "v"], [9, 8, 7])
      julia> dae.min("x"), ca.start(dae, "v")
      (-1.0, [8.0, 7.0])
      julia> try
                 ca.attribute(dae, "unit", "x")
             catch e
                 println(occursin("No such enum", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa min, max, nominal, start, get

      \identifier{2m6} */
  double attribute(const std::string& a, const std::string& name) const;

  /// Set an attribute, single variable (vector-valued variables get all elements set)
  void set_attribute(const std::string& a, const std::string& name, double val);

  /** \brief Get the lower bound, single variable

      The `"min"` attribute; `-inf` by default. Only scalar variables are supported.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(min=-2, max=2, nominal=10))
      >>> dae.min("x"), dae.max("x"), dae.nominal("x")
      (-2.0, 2.0, 10.0)
      >>> dae.set_min("x", 0)
      >>> dae.set_max(["x"], [5])
      >>> dae.set_nominal("x", 2)
      >>> dae.min(["x"]), dae.max("x"), dae.nominal(["x"])
      ([0.0], 5.0, [2.0])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('min', -2, 'max', 2, 'nominal', 10));
      >> dae.min('x'), dae.max('x'), dae.nominal('x')
      -2
      2
      10
      >> dae.set_min('x', 0); dae.set_max({'x'}, 5); dae.set_nominal('x', 2);
      >> dae.min({'x'}), dae.max('x'), dae.nominal({'x'})
      0
      5
      2
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("min" => -2, "max" => 2, "nominal" => 10));
      julia> dae.min("x"), dae.max("x"), ca.nominal(dae, "x")
      (-2.0, 2.0, 10.0)
      julia> ca.set_min(dae, "x", 0)
      julia> ca.set_max(dae, ["x"], [5])
      julia> ca.set_nominal(dae, "x", 2)
      julia> dae.min(["x"]), dae.max("x"), ca.nominal(dae, ["x"])
      ([0.0], 5.0, [2.0])
      \enddoctest
      \endlang

      \sa set_min, max, nominal

      \identifier{2m7} */
  double min(const std::string& name) const;

  /// Set the lower bound, single variable
  void set_min(const std::string& name, double val);

  /// Get the upper bound, single variable (`inf` by default)
  double max(const std::string& name) const;

  /// Set the upper bound, single variable
  void set_max(const std::string& name, double val);

  /// Get the nominal value, single variable (1 by default)
  double nominal(const std::string& name) const;

  /// Set the nominal value, single variable
  void set_nominal(const std::string& name, double val);

  /** \brief Get the start attribute, single variable

      The start values of a variable as a list with `numel` entries, also for scalars.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x", dict(start=1))
      >>> v = dae.add("v", dict(dimension=[2]))
      >>> dae.start("x"), dae.start("v"), dae.start(["x", "v"])
      ([1.0], [0.0, 0.0], [1.0, 0.0, 0.0])
      >>> dae.set_start("v", [3, 4])
      >>> dae.set_start(["x", "v"], [0, 1])
      >>> dae.start(["x", "v"])
      [0.0, 1.0, 1.0]
      >>> dae.set_start(["x", "v"], [1, 2, 3])
      >>> dae.start(["x", "v"])
      [1.0, 2.0, 3.0]
      >>> try:
      ...   dae.set_start("v", [1, 2, 3])
      ... except Exception as e:
      ...   print("expected number of elements (2) or number of variables (1)" in str(e))
      True
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x', struct('start', 1));
      >> v = dae.add('v', struct('dimension', 2));
      >> dae.start('x'), dae.start('v'), dae.start({'x', 'v'})
      1
      0 0
      1 0 0
      >> dae.set_start('v', [3 4]); dae.set_start({'x', 'v'}, [0 1]);
      >> dae.start({'x', 'v'})
      0 1 1
      >> dae.set_start({'x', 'v'}, [1 2 3]);
      >> dae.start({'x', 'v'})
      1 2 3
      >> try, dae.set_start('v', [1 2 3]); catch e, msg = e.message; end
      >> contains(msg, 'expected number of elements (2) or number of variables (1)')
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x", Dict("start" => 1));
      julia> v = ca.add(dae, "v", Dict("dimension" => [2]));
      julia> ca.start(dae, "x"), ca.start(dae, "v"), ca.start(dae, ["x", "v"])
      ([1.0], [0.0, 0.0], [1.0, 0.0, 0.0])
      julia> ca.set_start(dae, "v", [3, 4])
      julia> ca.set_start(dae, ["x", "v"], [0, 1])
      julia> ca.start(dae, ["x", "v"])
      3-element Vector{Float64}:
       0.0
       1.0
       1.0
      julia> ca.set_start(dae, ["x", "v"], [1, 2, 3])
      julia> ca.start(dae, ["x", "v"])
      3-element Vector{Float64}:
       1.0
       2.0
       3.0
      julia> try
                 ca.set_start(dae, "v", [1, 2, 3])
             catch e
                 msg = "expected number of elements (2) or number of variables (1)"
                 println(occursin(msg, sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \sa set_start, get

      \identifier{2m8} */
  std::vector<double> start(const std::string& name) const;

  /// Set the start attribute, single variable (all elements for a vector-valued variable)
  void set_start(const std::string& name, double val);

  /// Set the start attribute, vector argument (length `numel`, or 1 to broadcast)
  void set_start(const std::string& name, const std::vector<double>& val);

  /** \brief Clear all set values

      Resets the current value of every variable to nan (empty for string variables); the
      start attributes are kept.

      \sa set, get

      \identifier{2m9} */
  void reset();

  /** \brief Set the current value, single value

      The current value is what `get` returns and what `create` substitutes for the
      constants `c`; set it before creating functions. A variable's current value is
      initialized from its `"start"` option and is nan when unset. A vector-valued variable
      gets all elements set.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> c = dae.add("c", "parameter", "fixed")
      >>> x = dae.add("x", dict(start=1))
      >>> v = dae.add("v", dict(dimension=[2]))
      >>> s = dae.add("s", dict(type="String"))
      >>> dae.get("c"), dae.get("x"), dae.get("v"), repr(dae.get("s"))
      (nan, 1.0, [nan, nan], "''")
      >>> dae.set("c", 2.5)
      >>> dae.set("v", 1)
      >>> dae.set("s", "hello")
      >>> dae.get(["c", "v", "s"])
      [2.5, [1.0, 1.0], 'hello']
      >>> dae.set(["c", "x"], [3, 4])
      >>> dae.set(["v"], [5, 6])
      >>> dae.get(["c", "x", "v"])
      [3.0, 4.0, [5.0, 6.0]]
      >>> dae.reset()
      >>> dae.get(["c", "x"])
      [nan, nan]
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); c = dae.add('c', 'parameter', 'fixed');
      >> x = dae.add('x', struct('start', 1)); v = dae.add('v', struct('dimension', 2));
      >> s = dae.add('s', struct('type', 'String'));
      >> val = dae.get({'c', 'x', 'v', 's'});
      >> val{1}, val{2}, val{3}, isempty(val{4})
      NaN
      1
      NaN NaN
      1
      >> dae.set('c', 2.5); dae.set('v', 1); dae.set('s', 'hello');
      >> val = dae.get({'c', 'v', 's'}); disp(val{1}), disp(val{2}), disp(val{3})
      2.5000
      1 1
      hello
      >> dae.set({'c', 'x'}, [3 4]); dae.set({'v'}, [5 6]);
      >> val = dae.get({'c', 'x', 'v'}); disp(val{1}), disp(val{2}), disp(val{3})
      3
      4
      5 6
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> c = ca.add(dae, "c", "parameter", "fixed");
      julia> x = ca.add(dae, "x", Dict("start" => 1));
      julia> v = ca.add(dae, "v", Dict("dimension" => [2]));
      julia> s = ca.add(dae, "s", Dict("type" => "String"));
      julia> ca.get(dae, "c"), ca.get(dae, "x"), ca.get(dae, "v"), ca.get(dae, "s")
      (NaN, 1.0, [NaN, NaN], "")
      julia> ca.set(dae, "c", 2.5)
      julia> ca.set(dae, "v", 1)
      julia> ca.set(dae, "s", "hello")
      julia> ca.get(dae, ["c", "v", "s"])
      3-element Vector{Any}:
       2.5
        [1.0, 1.0]
        "hello"
      julia> ca.set(dae, ["c", "x"], [3, 4])
      julia> ca.set(dae, ["v"], [5, 6])
      julia> ca.get(dae, ["c", "x", "v"])
      3-element Vector{Any}:
       3.0
       4.0
        [5.0, 6.0]
      julia> ca.reset(dae)
      julia> ca.get(dae, ["c", "x"])
      2-element Vector{Any}:
       NaN
       NaN
      \enddoctest
      \endlang

      The value of a constant is baked into the functions created afterwards:

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> c = dae.add("c", "parameter", "fixed")
      >>> x = dae.add("x")
      >>> dae.eq(dae.der(x), -c*x)
      >>> dae.set("c", 3)
      >>> dae.create("f", ["x"], ["ode"])(2)
      DM(-6)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); c = dae.add('c', 'parameter', 'fixed'); x = dae.add('x');
      >> dae.eq(dae.der(x), -c*x);
      >> dae.set('c', 3);
      >> f = dae.create('f', {'x'}, {'ode'}); disp(f(2))
      -6
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> c = ca.add(dae, "c", "parameter", "fixed");
      julia> x = ca.add(dae, "x");
      julia> ca.eq(dae, ca.der(dae, x), -c*x)
      julia> ca.set(dae, "c", 3)
      julia> ca.create(dae, "f", ["x"], ["ode"])(2)
      DM(-6)
      \enddoctest
      \endlang

      \sa get, reset, start, c

      \identifier{2ma} */
  void set(const std::string& name, double val);

  /// Set the current value, single value (string), for variables of type `"String"`
  void set(const std::string& name, const std::string& val);

  /** \brief Evaluate the values for a set of variables at the initial time, single value

      For a symbolic model the current value (see `set`): a number for a scalar, a list
      for a vector-valued variable, a string for a string variable; nan when unset. For an
      imported FMU the value is obtained by evaluating the FMU.

      \sa set, start

      \identifier{2mb} */
  GenericType get(const std::string& name) const;

  /** \brief Symbolic DAE instance?

      True for a symbolically constructed model (or a symbolic model description), false
      for an imported standard FMU whose equations live in a binary.

      \lang{python}
      \doctest
      >>> DaeBuilder("m").symbolic()
      True
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.symbolic(ca.DaeBuilder("m"))
      true
      \enddoctest
      \endlang

      \sa load_fmi_description, oracle

      \identifier{2mc} */
  bool symbolic() const;

#endif  // !SWIGMATLAB

  /// Get an attribute for several variables, concatenated (`"start"`/`"value"` per element)
  std::vector<double> attribute(const std::string& a, const std::vector<std::string>& name) const;

  /// Set an attribute, one value per variable or one value per element
  void set_attribute(const std::string& a, const std::vector<std::string>& name,
    const std::vector<double>& val);

  /// Get the lower bound for several (scalar) variables
  std::vector<double> min(const std::vector<std::string>& name) const;

  /// Set the lower bound for several variables
  void set_min(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the upper bound for several (scalar) variables
  std::vector<double> max(const std::vector<std::string>& name) const;

  /// Set the upper bound for several variables
  void set_max(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the nominal value for several (scalar) variables
  std::vector<double> nominal(const std::vector<std::string>& name) const;

  /// Set the nominal value for several variables
  void set_nominal(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Get the start attribute for several variables, concatenated
  std::vector<double> start(const std::vector<std::string>& name) const;

  /** \brief Set the start attribute for several variables

      \p val holds either one value per variable (broadcast over the elements) or one value
      per element of all variables concatenated; any other length is an error.

      \identifier{2md} */
  void set_start(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Set the current value for several variables (one per variable or one per element)
  void set(const std::vector<std::string>& name, const std::vector<double>& val);

  /// Set the current value (string) for several string variables
  void set(const std::vector<std::string>& name, const std::vector<std::string>& val);

  /// Evaluate the values for a set of variables at the initial time, one entry per name
  std::vector<GenericType> get(const std::vector<std::string>& name) const;

  /** \brief Check if a particular variable exists

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> dae.has("x"), dae.has("der(x)"), dae.has("time")
      (True, False, True)
      >>> dae.der(x), dae.has("der(x)")
      (MX(der(x)), True)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x');
      >> dae.has('x'), dae.has('der(x)'), dae.has('time')
      1
      0
      1
      >> dae.der(x), dae.has('der(x)')
      der(x)
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> ca.has(dae, "x"), ca.has(dae, "der(x)"), ca.has(dae, "time")
      (true, false, true)
      julia> ca.der(dae, x), ca.has(dae, "der(x)")
      (der(x), true)
      \enddoctest
      \endlang

      \sa all, var

      \identifier{2me} */
  bool has(const std::string& name) const;

  /** \brief Get a list of all variables

      Names of all variables in creation order, including the hidden `"time"` variable
      and the internal helper variables created by `eq`, `when`, `assign` and `reinit`.

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> dae.eq(dae.der(x), -x)
      >>> dae.all()
      ['time', 'x', 'der(x)', '__assign__der(x)__']
      >>> dae.all("x"), dae.all("w"), dae.all("calculated")
      (['x'], ['der(x)'], ['__assign__der(x)__'])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); dae.eq(dae.der(x), -x);
      >> dae.all()
      time
      x
      der(x)
      __assign__der(x)__
      >> dae.all('x'), dae.all('w'), dae.all('calculated')
      x
      der(x)
      __assign__der(x)__
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> ca.eq(dae, ca.der(dae, x), -x)
      julia> ca.all(dae)
      4-element Vector{String}:
       "time"
       "x"
       "der(x)"
       "__assign__der(x)__"
      julia> ca.all(dae, "x"), ca.all(dae, "w"), ca.all(dae, "calculated")
      (["x"], ["der(x)"], ["__assign__der(x)__"])
      \enddoctest
      \endlang

      \sa has, x, w

      \identifier{2mf} */
  std::vector<std::string> all() const;

  /** \brief Get a list of all variables of a particular category

      Names of the members of a category in category order; empty for an empty category,
      an unknown category name is an error. `x()`, `z()`, `u()` etc. are shorthands.

      \identifier{2mg} */
  std::vector<std::string> all(const std::string& cat) const;

  /** \brief Get the (cached) oracle, SX or MX

      The Function behind `create` for a symbolic model: inputs are all input categories
      except `c` (whose current values are substituted), outputs are all output categories,
      each as one concatenated vector (empty vectors for empty categories). The oracle is
      cached and rebuilt when variables or equations change. Not available for an imported
      FMU.

      \param sx Return the SX expansion of the oracle instead of the MX graph
      \param elim_w Substitute the dependent variables `w` out and drop them as input
      \param lifted_calls Replace embedded function calls in `wdef` by zeros, for use with
        `create` with `lifted_calls`; incompatible with \p elim_w

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> w = dae.add("w")
      >>> dae.eq(w, sin(x))
      >>> dae.eq(dae.der(x), -w)
      >>> f = dae.oracle()
      >>> f.name(), f.name_in()
      ('mx_oracle', ['t', 'p', 'd', 'w', 'u', 'x', 'z', 'q'])
      >>> f.name_out()
      ['ode', 'alg', 'quad', 'zero', 'ddef', 'wdef', 'y']
      >>> f.size1_in("w"), f.size1_in("p")
      (2, 0)
      >>> dae.oracle(True).class_name(), dae.oracle(False, True).name_in()
      ('SXFunction', ['t', 'p', 'd', 'u', 'x', 'z', 'q'])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); w = dae.add('w');
      >> dae.eq(w, sin(x)); dae.eq(dae.der(x), -w);
      >> f = dae.oracle();
      >> f.name(), strjoin(cellstr(f.name_in()), ' ')
      mx_oracle
      t p d w u x z q
      >> strjoin(cellstr(f.name_out()), ' ')
      ode alg quad zero ddef wdef y
      >> f.size1_in('w'), f.size1_in('p')
      2
      0
      >> g = dae.oracle(true); disp(g.class_name())
      SXFunction
      >> g = dae.oracle(false, true); disp(strjoin(cellstr(g.name_in()), ' '))
      t p d u x z q
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> w = ca.add(dae, "w");
      julia> ca.eq(dae, w, sin(x))
      julia> ca.eq(dae, ca.der(dae, x), -w)
      julia> f = ca.oracle(dae);
      julia> ca.name(f), ca.name_in(f)
      ("mx_oracle", ["t", "p", "d", "w", "u", "x", "z", "q"])
      julia> ca.name_out(f)
      7-element Vector{String}:
       "ode"
       "alg"
       "quad"
       "zero"
       "ddef"
       "wdef"
       "y"
      julia> ca.size1_in(f, "w"), ca.size1_in(f, "p")
      (2, 0)
      julia> ca.class_name(ca.oracle(dae, true)), ca.name_in(ca.oracle(dae, false, true))
      ("SXFunction", ["t", "p", "d", "u", "x", "z", "q"])
      \enddoctest
      \endlang

      \sa create, symbolic

      \identifier{2mh} */
  Function oracle(bool sx = false, bool elim_w = false, bool lifted_calls = false) const;

  /** \brief Get Jacobian sparsity

      Structural dependency pattern of the variables \p onames with respect to \p inames,
      one row per output name and one column per input name, taken from the dependency
      information of the model structure. This information comes from an imported model
      description or is computed by `export_fmu`; for a symbolically constructed model it
      is otherwise absent and the pattern has no nonzeros.

      \param onames Names of (derivative or output) variables
      \param inames Names of variables to differentiate with respect to
      \return Sparsity pattern of size `len(onames)` by `len(inames)`

      \lang{python}
      \doctest
      >>> dae = DaeBuilder("m")
      >>> x = dae.add("x")
      >>> u = dae.add("u", "input")
      >>> dae.eq(dae.der(x), -x*u)
      >>> dae.jac_sparsity(["der(x)"], ["x", "u"])
      Sparsity(1x2,0nz)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> dae = DaeBuilder('m'); x = dae.add('x'); u = dae.add('u', 'input');
      >> dae.eq(dae.der(x), -x*u);
      >> dae.jac_sparsity({'der(x)'}, {'x', 'u'})
      1x2,0nz
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> dae = ca.DaeBuilder("m");
      julia> x = ca.add(dae, "x");
      julia> u = ca.add(dae, "u", "input");
      julia> ca.eq(dae, ca.der(dae, x), -x*u)
      julia> ca.jac_sparsity(dae, ["der(x)"], ["x", "u"])
      Sparsity(1x2,0nz)
      \enddoctest
      \endlang

      \sa load_fmi_description, export_fmu

      \identifier{6g} */
  Sparsity jac_sparsity(const std::vector<std::string>& onames,
    const std::vector<std::string>& inames) const;

#ifndef SWIG
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
