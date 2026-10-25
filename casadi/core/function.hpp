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

#ifndef CASADI_FUNCTION_HPP
#define CASADI_FUNCTION_HPP

#include "sx_elem.hpp"
#include "mx.hpp"
#include "printable.hpp"
#include <exception>
#include <stack>

namespace casadi {

#ifndef SWIG
  /** Forward declaration of internal class */
  class FunctionInternal;
  class SerializingStream;
  class DeserializingStream;
#endif // SWIG

  /** \brief Function object

      A Function instance is a general multiple-input, multiple-output function
      where each input and output can be a sparse matrix.\n

      For an introduction to this class, see the CasADi user guide.\n

      Function is a reference counted and immutable class; copying a class instance
      is very cheap and its behavior (with some exceptions) is not affected by
      calling its member functions.\n

      Each input and output has a name (name_in, name_out) and a fixed sparsity pattern
      (sparsity_in, sparsity_out). The same object can be evaluated numerically (DM
      arguments) or symbolically (SX or MX arguments), differentiated (jacobian, forward,
      reverse, factory), transformed (map, mapaccum, fold, expand), code-generated
      (generate), serialized (serialize, save) and inspected (print_dimensions, stats, info).

      Instances are created from expressions (the SX and MX constructors), from a C body
      (jit), by combining other functions (conditional, if_else, map, mapaccum) or by the
      solver constructors nlpsol, qpsol, rootfinder, integrator, external. Printing shows
      the signature and the class name, e.g. f:(x,y[2])->(r[2],s) SXFunction.

      The exceptions to immutability are options changed with change_option and the cache
      of derivative functions.

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); y = SX.sym("y", 2)
      >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
      >>> f
      Function(f:(x,y[2])->(r[2],s) SXFunction)
      >>> r, s = f(2, [3, 4])
      >>> r, s
      (DM([6, 8]), DM(0.909297))
      >>> f(x=2, y=[3, 4])
      {'r': DM([6, 8]), 's': DM(0.909297)}
      >>> z = MX.sym("z")
      >>> f(z, [3, 4])
      (MX(f(z, [3, 4]){0}), MX(f(z, [3, 4]){1}))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); y = SX.sym('y', 2);
      >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
      >> f
      f:(x,y[2])->(r[2],s) SXFunction
      >> [r, s] = f(2, [3; 4]); disp(r), disp(s)
      [6, 8]
      0.909297
      >> res = f('x', 2, 'y', [3; 4]); disp(res.r), disp(res.s)
      [6, 8]
      0.909297
      >> z = MX.sym('z'); [r, s] = f(z, [3; 4]); disp(r), disp(s)
      f(z, [3, 4]){0}
      f(z, [3, 4]){1}
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
      julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
      Function(f:(x,y[2])->(r[2],s) SXFunction)
      julia> r, s = f(2, [3, 4])
      ([6, 8], 0.909297)
      julia> f(x=2, y=[3, 4])
      Dict{String, DM} with 2 entries:
        "r" => [6, 8]
        "s" => 0.909297
      julia> z = ca.MX.sym("z");
      julia> f(z, [3, 4])
      (f(z, [3, 4]){0}, f(z, [3, 4]){1})
      \enddoctest
      \endlang

      \author Joel Andersson
      \date 2010-2017

      \identifier{1uw} */
  class CASADI_EXPORT Function :
    public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Function>) {
  public:
    /** \brief Get the type name, always "Function"

        \lang{python}
        \doctest
        >>> Function.type_name()
        'Function'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Function.type_name()
        Function
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.Function.type_name()
        "Function"
        \enddoctest
        \endlang

        \identifier{1ux} */
    static std::string type_name() {return "Function";}

    /** \brief Default constructor, creates a null object

        A null Function has name "null", is_null() is true, and any other member
        function raises an error. Use it only as a placeholder to be assigned later.

        \lang{python}
        \doctest
        >>> f = Function()
        >>> f.is_null(), f.name()
        (True, 'null')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function(); disp(f.is_null()), disp(f.name())
        1
        null
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function();
        julia> ca.is_null(f), ca.name(f)
        (true, "null")
        \enddoctest
        \endlang

        \identifier{1uy} */
    Function();

#ifndef SWIG
    /** \brief Construct from a file; not implemented, raises an error

        Use Function::load to read a file written by Function::save.

        \sa load, save

        \identifier{1uz} */
    Function(const std::string& fname);
#endif // SWIG

    ///@{
    /** \brief Construct a function from SX input and output expressions

        The result is an SXFunction: a flat list of scalar instructions that is cheap to
        evaluate but has no notion of matrix structure (see expand for the difference with
        an MXFunction).

        Every symbolic primitive that appears in \p ex_out must be present in \p ex_in,
        otherwise construction fails with "Initialization failed since variables [...] are
        free" (the symbol was forgotten in the input list). Pass option "allow_free": true
        to construct anyway; such a function can only be called symbolically.
        Each input must be a purely symbolic matrix (symbols or a concatenation of symbols);
        an output may be any expression, including constants and empty matrices. The
        sparsity patterns of the expressions become sparsity_in and sparsity_out.

        \param name Display name, used in error messages and generated code. Must start
          with a letter and contain only letters, digits and non-consecutive underscores
          (see check_name); otherwise an error is raised.
        \param ex_in List of input expressions, each a symbolic SX
        \param ex_out List of output expressions
        \param name_in Names of the inputs; defaults to i0, i1, ...; duplicates raise an
          error unless option "allow_duplicate_io_names" is set
        \param name_out Names of the outputs; defaults to o0, o1, ...
        \param opts Options, see the table below. The common ones are "cse" (eliminate
          common subexpressions), "live_variables", "is_diff_in"/"is_diff_out",
          "default_in", "allow_free", "jit"

        The dict overload takes a single map from name to expression; \p name_in and
        \p name_out select which keys are inputs and outputs, in that order. A key that
        is in neither list raises "Unknown dictionary entry". A missing key gives an
        empty (0-by-0) input or output.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y, sin(x)])
        >>> f
        Function(f:(i0,i1[2])->(o0[2],o1) SXFunction)
        >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
        >>> f.name_in(), f.name_out()
        (['x', 'y'], ['r', 's'])
        >>> g = Function("g", {"x": x, "y": y, "r": x*y}, ["x", "y"], ["r"])
        >>> g
        Function(g:(x,y[2])->(r[2]) SXFunction)
        >>> p = SX.sym("p")
        >>> try:
        ...   Function("h", [x], [x*p])
        ... except RuntimeError as e:
        ...   print("Initialization failed" in str(e), "[p] are free" in str(e))
        True True
        >>> h = Function("h", [x], [x*p], {"allow_free": True})
        >>> h.get_free()
        ['p']
        >>> try:
        ...   Function("my func", [x], [x])
        ... except RuntimeError as e:
        ...   print("Function name is not valid" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y, sin(x)}); disp(f)
        f:(i0,i1[2])->(o0[2],o1) SXFunction
        >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
        >> f.name_in(), f.name_out()
        x
        y
        r
        s
        >> g = Function('g', struct('x', x, 'y', y, 'r', x*y), {'x', 'y'}, {'r'}); disp(g)
        g:(x,y[2])->(r[2]) SXFunction
        >> p = SX.sym('p');
        >> try, Function('h', {x}, {x*p}); catch e, disp(contains(e.message, '[p] are free')); end
        1
        >> h = Function('h', {x}, {x*p}, struct('allow_free', true)); disp(h.get_free())
        p
        >> try, Function('my func', {x}, {x}); catch e, disp(contains(e.message, 'not valid')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)])
        Function(f:(i0,i1[2])->(o0[2],o1) SXFunction)
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
        julia> ca.name_in(f), ca.name_out(f)
        (["x", "y"], ["r", "s"])
        julia> g = ca.Function("g", Dict("x" => x, "y" => y, "r" => x*y), ["x", "y"], ["r"])
        Function(g:(x,y[2])->(r[2]) SXFunction)
        julia> p = ca.SX.sym("p");
        julia> try
                   ca.Function("h", [x], [x*p])
               catch e
                   msg = sprint(showerror, e)
                   println(occursin("Initialization failed", msg), " ",
                           occursin("[p] are free", msg))
               end
        true true
        julia> h = ca.Function("h", [x], [x*p], Dict("allow_free" => true));
        julia> ca.get_free(h)
        1-element Vector{String}:
         "p"
        julia> try
                   ca.Function("my func", [x], [x])
               catch e
                   println(occursin("Function name is not valid", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa expand, check_name, has_free

        \identifier{1v0} */
    Function(const std::string& name,
             const std::vector<SX>& ex_in,
             const std::vector<SX>& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<SX>& ex_in,
             const std::vector<SX>& ex_out,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, SX>& dict,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief Construct a function from MX input and output expressions

        The result is an MXFunction: a graph of matrix-valued nodes that may embed calls
        to other Function objects (solvers, integrators, SX functions, ...). Compared with
        an SXFunction it uses less memory and supports every operation, but is slower per
        node to evaluate; call expand to convert when all nodes are expandable.

        The rules on free variables, naming and the dict overload are the same as for the
        SX constructor. Additional MX-only options: "post_expand" (expand right after
        construction) and "post_expand_options".

        \param name Display name, must satisfy check_name
        \param ex_in List of symbolic MX inputs
        \param ex_out List of MX output expressions
        \param name_in Input names; defaults to i0, i1, ...
        \param name_out Output names; defaults to o0, o1, ...
        \param opts Options, see the table below

        \lang{python}
        \doctest
        >>> x = MX.sym("x", 2); A = MX.sym("A", 2, 2)
        >>> f = Function("f", [A, x], [A @ x, norm_2(x)], ["A", "x"], ["Ax", "nrm"])
        >>> f
        Function(f:(A[2x2],x[2])->(Ax[2],nrm) MXFunction)
        >>> f.class_name(), f.n_instructions()
        ('MXFunction', 7)
        >>> f.expand()
        Function(f:(A[2x2],x[2])->(Ax[2],nrm) SXFunction)
        >>> g = Function("g", [x], [x], {"post_expand": True})
        >>> g.class_name()
        'SXFunction'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x', 2); A = MX.sym('A', 2, 2);
        >> f = Function('f', {A, x}, {A*x, norm_2(x)}, {'A', 'x'}, {'Ax', 'nrm'}); disp(f)
        f:(A[2x2],x[2])->(Ax[2],nrm) MXFunction
        >> f.class_name(), f.n_instructions()
        MXFunction
        7
        >> f.expand()
        f:(A[2x2],x[2])->(Ax[2],nrm) SXFunction
        >> g = Function('g', {x}, {x}, struct('post_expand', true)); disp(g.class_name())
        SXFunction
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x", 2); A = ca.MX.sym("A", 2, 2);
        julia> f = ca.Function("f", [A, x], [A * x, ca.norm_2(x)], ["A", "x"], ["Ax", "nrm"])
        Function(f:(A[2x2],x[2])->(Ax[2],nrm) MXFunction)
        julia> ca.class_name(f), ca.n_instructions(f)
        ("MXFunction", 7)
        julia> ca.expand(f)
        Function(f:(A[2x2],x[2])->(Ax[2],nrm) SXFunction)
        julia> g = ca.Function("g", [x], [x], Dict("post_expand" => true));
        julia> ca.class_name(g)
        "SXFunction"
        \enddoctest
        \endlang

        \sa expand, has_free

        \identifier{1v1} */
    Function(const std::string& name,
             const std::vector<MX>& ex_in,
             const std::vector<MX>& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name,
             const std::vector<MX>& ex_in,
             const std::vector<MX>& ex_out,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const std::map<std::string, MX>& dict,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    ///@}

    ///@{
    /** \brief To resolve ambiguity on some compilers

        \identifier{1v2} */
#ifndef SWIG
    Function(const std::string& name, SXIList ex_in,
      const SXVector& ex_out, const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& ex_in,
      SXIList ex_out, const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in,
      SXIList ex_out, const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in, const SXVector& ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const SXVector& ex_in, SXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, SXIList ex_in, SXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, const MXVector& ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& ex_in, MXIList ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, MXIList ex_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, const MXVector& ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, const MXVector& ex_in, MXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
    Function(const std::string& name, MXIList ex_in, MXIList ex_out,
             const StringVector& name_in, const StringVector& name_out,
             const Dict& opts=Dict());
#endif // SWIG
    ///@}

    ///@{
    /** \brief Create a just-in-time compiled function from a C language string

        \p body is the body of a C function that is compiled with the "compiler" plugin
        (option "compiler", default "shell" i.e. the system C compiler on most builds, or
        "clang") and loaded as an external function. This needs a working C compiler at
        run time; the source and library files are removed afterwards unless
        "jit_cleanup" is false.

        Inside the body, the nonzeros of input k are available through a pointer named
        name_in[k] (type const casadi_real*) and the nonzeros of output k must be written
        through a pointer named name_out[k] (type casadi_real*), in column-major order of
        the sparsity pattern. casadi_real is double unless the code generator is told
        otherwise. By default (option "buffered": true) inputs and outputs are copied
        through a work vector, so the pointers are never null; with "buffered": false a
        null input means all-zero and a null output may be skipped. Write "return 1;" to
        signal an evaluation failure.

        The result has no derivatives unless option "enable_fd": true (finite
        differences) is passed or a C body for the Jacobian is supplied with option "jac"
        (and "hess" for its Jacobian).

        \param name Function name
        \param body C code of the function body, without the enclosing braces
        \param name_in Input names, also the C variable names
        \param name_out Output names, also the C variable names
        \param sparsity_in Sparsity of each input; all scalar when omitted
        \param sparsity_out Sparsity of each output; all scalar when omitted
        \param opts Options; "compiler", "jit_options", "buffered", "jac", "hess" and the
          common Function options

        Example (not run here because it invokes a compiler):
        f = Function.jit("f", "r[0] = x[0]*x[0] + y[0];", ["x", "y"], ["r"]) gives a
        function f:(x,y)->(r) that evaluates like any other Function.

        \sa external

        \identifier{1v3} */
    static Function jit(const std::string& name, const std::string& body,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const Dict& opts=Dict());
    static Function jit(const std::string& name, const std::string& body,
             const std::vector<std::string>& name_in,
             const std::vector<std::string>& name_out,
             const std::vector<Sparsity>& sparsity_in,
             const std::vector<Sparsity>& sparsity_out,
             const Dict& opts=Dict());
    ///@}

    /** \brief  Destructor

        \identifier{1v4} */
    ~Function();

    /** \brief Convert an MXFunction into an equivalent SXFunction

        Evaluates the function symbolically with SX inputs and builds an SXFunction from
        the result: every matrix node is unrolled into scalar instructions. Evaluation
        of the result is typically faster per operation, at the cost of memory
        proportional to the number of scalar operations. The name, input/output names
        and the options of the original function (other than the SX-specific ones given in
        \p opts) are kept.

        Fails when the graph contains a node that has no scalar expansion, with
        "eval_sx not defined for ...": for example calls to external functions or
        linear solves (solve, inv) with a numerical plugin; use the "symbolicqr" plugin for
        those. Fails when the function has free variables. An SXFunction is returned
        unchanged (as a copy).

        \param name Name of the new function; defaults to the current name
        \param opts Options for the SXFunction constructor, e.g. "live_variables", "cse"
        \return SXFunction

        \lang{python}
        \doctest
        >>> x = MX.sym("x", 2)
        >>> f = Function("f", [x], [x.T @ x, 2*x], ["x"], ["q", "y"])
        >>> f, f.n_instructions()
        (Function(f:(x[2])->(q,y[2]) MXFunction), 7)
        >>> g = f.expand()
        >>> g, g.n_instructions()
        (Function(f:(x[2])->(q,y[2]) SXFunction), 10)
        >>> f([1, 2]), g([1, 2])
        ((DM(5), DM([2, 4])), (DM(5), DM([2, 4])))
        >>> f.expand("f_sx", {"live_variables": False})
        Function(f_sx:(x[2])->(q,y[2]) SXFunction)
        >>> A = MX.sym("A", 2, 2)
        >>> h = Function("h", [A, x], [solve(A, x)])
        >>> try:
        ...   h.expand()
        ... except RuntimeError as e:
        ...   print("eval_sx not defined" in str(e))
        True
        >>> Function("h", [A, x], [solve(A, x, "symbolicqr")]).expand().class_name()
        'SXFunction'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x', 2);
        >> f = Function('f', {x}, {x'*x, 2*x}, {'x'}, {'q', 'y'});
        >> f, f.n_instructions()
        f:(x[2])->(q,y[2]) MXFunction
        7
        >> g = f.expand(); disp(g), disp(g.n_instructions())
        f:(x[2])->(q,y[2]) SXFunction
        10
        >> [q, y] = g([1; 2]); disp(q), disp(y)
        5
        [2, 4]
        >> f.expand('f_sx', struct('live_variables', false))
        f_sx:(x[2])->(q,y[2]) SXFunction
        >> A = MX.sym('A', 2, 2); h = Function('h', {A, x}, {solve(A, x)});
        >> try, h.expand(); catch e, disp(contains(e.message, 'eval_sx not defined')); end
        1
        >> Function('h', {A, x}, {solve(A, x, 'symbolicqr')}).expand().class_name()
        SXFunction
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x", 2);
        julia> f = ca.Function("f", [x], [x' * x, 2x], ["x"], ["q", "y"]);
        julia> println(f); ca.n_instructions(f)
        f:(x[2])->(q,y[2]) MXFunction
        7
        julia> g = ca.expand(f);
        julia> println(g); ca.n_instructions(g)
        f:(x[2])->(q,y[2]) SXFunction
        10
        julia> f([1, 2]), g([1, 2])
        ((5, [2, 4]), (5, [2, 4]))
        julia> ca.expand(f, "f_sx", Dict("live_variables" => false))
        Function(f_sx:(x[2])->(q,y[2]) SXFunction)
        julia> A = ca.MX.sym("A", 2, 2);
        julia> h = ca.Function("h", [A, x], [A \ x]);
        julia> try
                   ca.expand(h)
               catch e
                   println(occursin("eval_sx not defined", sprint(showerror, e)))
               end
        true
        julia> ca.class_name(ca.expand(ca.Function("h", [A, x], [ca.solve(A, x, "symbolicqr")])))
        "SXFunction"
        \enddoctest
        \endlang

        \sa Function::Function, transform

        \identifier{1v5} */
    ///@{
    Function expand() const;
    Function expand(const std::string& name,
                    const Dict& opts=Dict()) const;
    ///@}

    /** \brief Rebuild the function after graph simplification or other passes

        Only defined for SXFunction and MXFunction. The result has the same signature,
        names and options; only the expression graph differs.

        The dict form runs a single "simplify" pass whose tasks are selected with boolean
        options, all defaulting to on except combine_terms:
        - empty_inputs: an input whose symbols do not occur in any output is replaced by
          an all-zero sparse matrix of the same shape (its nnz becomes 0)
        - cse: common subexpression elimination
        - ref_count: drop unreferenced nodes, repeated until a fixed point
        - const_folding: evaluate constant subexpressions
        - combine_terms: merge like terms of sums

        Unknown keys raise an error.

        The passes form runs an explicit list of passes in order. Each pass is a list whose
        first entry is the verb:
        - ["simplify", task, ...] with tasks from the list above; an integer before a task
          is its run count (N>0: N times, 0: until the graph stops shrinking)
        - ["expand"]: expand to an SXFunction
        - ["external", library, operation] or [..., opts]: a transform from a plugin
          library

        Here \p opts accepts only "verbose"; passing the boolean task options raises an
        error. With "verbose": true the node and instruction counts after each pass are
        printed.

        \param fname Name of the result; defaults to the current name
        \param passes Explicit list of passes
        \param opts Dict form: task booleans and "verbose"; passes form: "verbose" only

        \lang{python}
        \doctest
        >>> z = MX.sym("z"); w = MX.sym("w")
        >>> f = Function("f", [z, w], [sin(z)*cos(z) + sin(z)*cos(z)], ["z", "w"], ["r"])
        >>> f, f.n_instructions()
        (Function(f:(z,w)->(r) MXFunction), 9)
        >>> g = f.transform()
        >>> g, g.n_instructions()
        (Function(f:(z,w[1x1,0nz])->(r) MXFunction), 6)
        >>> f.transform([["simplify", "cse"]]).n_instructions()
        6
        >>> f.transform("fe", [["expand"]])
        Function(fe:(z,w)->(r) SXFunction)
        >>> try:
        ...   f.transform([["nope"]])
        ... except RuntimeError as e:
        ...   print("unknown pass verb" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> z = MX.sym('z'); w = MX.sym('w');
        >> f = Function('f', {z, w}, {sin(z)*cos(z) + sin(z)*cos(z)}, {'z', 'w'}, {'r'});
        >> f, f.n_instructions()
        f:(z,w)->(r) MXFunction
        9
        >> g = f.transform(); disp(g), disp(g.n_instructions())
        f:(z,w[1x1,0nz])->(r) MXFunction
        6
        >> f.transform({{'simplify', 'cse'}}).n_instructions()
        6
        >> f.transform('fe', {{'expand'}})
        fe:(z,w)->(r) SXFunction
        >> try, f.transform({{'nope'}}); catch e, disp(contains(e.message, 'unknown pass')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> z = ca.MX.sym("z"); w = ca.MX.sym("w");
        julia> f = ca.Function("f", [z, w], [sin(z)*cos(z) + sin(z)*cos(z)], ["z", "w"], ["r"]);
        julia> println(f); ca.n_instructions(f)
        f:(z,w)->(r) MXFunction
        9
        julia> g = ca.transform(f);
        julia> println(g); ca.n_instructions(g)
        f:(z,w[1x1,0nz])->(r) MXFunction
        6
        julia> ca.n_instructions(ca.transform(f, [["simplify", "cse"]]))
        6
        julia> ca.transform(f, "fe", [["expand"]])
        Function(fe:(z,w)->(r) SXFunction)
        julia> try
                   ca.transform(f, [["nope"]])
               catch e
                   println(occursin("unknown pass verb", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa expand

        \identifier{2ig} */
    ///@{
    Function transform(const Dict& opts = Dict()) const;
    Function transform(const std::string& fname, const Dict& opts = Dict()) const;
    Function transform(const std::vector<std::vector<GenericType> >& passes,
                       const Dict& opts = Dict()) const;
    Function transform(const std::string& fname,
                       const std::vector<std::vector<GenericType> >& passes,
                       const Dict& opts = Dict()) const;
    ///@}

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{1v6} */
    static Function create(FunctionInternal* node);

    /** \brief  Create from node and initialize

        \identifier{1v7} */
    static Function create(FunctionInternal* node, const Dict& opts);
#endif // SWIG
    /// \endcond

    /** \brief Get the number of function inputs

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> f = Function("f", [x, y], [x+y, x-y, x*y])
        >>> f.n_in(), f.n_out()
        (2, 3)
        >>> Function("g", [], [DM(1)]).n_in()
        0
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> f = Function('f', {x, y}, {x+y, x-y, x*y});
        >> fprintf('%d %d\n', f.n_in(), f.n_out())
        2 3
        >> Function('g', {}, {DM(1)}).n_in()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> f = ca.Function("f", [x, y], [x + y, x - y, x*y]);
        julia> ca.n_in(f), ca.n_out(f)
        (2, 3)
        julia> ca.n_in(ca.Function("g", [], [ca.DM(1)]))
        0
        \enddoctest
        \endlang

        \sa n_out, name_in, size_in

        \identifier{1v8} */
    casadi_int n_in() const;

    /** \brief Get the number of function outputs

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> Function("f", [x], [x, 2*x]).n_out()
        2
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> Function('f', {x}, {x, 2*x}).n_out()
        2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.n_out(ca.Function("f", [x], [x, 2x]))
        2
        \enddoctest
        \endlang

        \sa n_in, name_out, size_out

        \identifier{1v9} */
    casadi_int n_out() const;

    ///@{
    /** \brief Get the shape of an input: rows, columns or both

        size1_in is the number of rows, size2_in the number of columns and size_in the
        pair (rows, columns). The input is selected by index or by name; an unknown name
        raises an error listing the available names, an index out of range raises an
        error.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("A", 2, 3), SX.sym("b", 2)], [SX(0)], ["A", "b"], ["c"])
        >>> f.size_in(0), f.size1_in("A"), f.size2_in("A")
        ((2, 3), 2, 3)
        >>> f.size_in("b")
        (2, 1)
        >>> try:
        ...   f.size_in("nope")
        ... except RuntimeError as e:
        ...   print("Available names are: [A, b]" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('A', 2, 3), SX.sym('b', 2)}, {SX(0)}, {'A', 'b'}, {'c'});
        >> disp(f.size_in(0)), fprintf('%d %d\n', f.size1_in('A'), f.size2_in('A'))
        2 3
        2 3
        >> f.size_in('b')
        2 1
        >> try, f.size_in('nope'); catch e, disp(contains(e.message, 'Available names')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> A = ca.SX.sym("A", 2, 3); b = ca.SX.sym("b", 2);
        julia> f = ca.Function("f", [A, b], [ca.SX(0)], ["A", "b"], ["c"]);
        julia> ca.size_in(f, 0), ca.size1_in(f, "A"), ca.size2_in(f, "A")
        ((2, 3), 2, 3)
        julia> ca.size_in(f, "b")
        (2, 1)
        julia> try
                   ca.size_in(f, "nope")
               catch e
                   println(occursin("Available names are: [A, b]", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa size_out, sparsity_in, nnz_in, numel_in

        \identifier{1va} */
    casadi_int size1_in(casadi_int ind) const;
    casadi_int size1_in(const std::string& iname) const {return size1_in(index_in(iname));}
    casadi_int size2_in(casadi_int ind) const;
    casadi_int size2_in(const std::string& iname) const {return size2_in(index_in(iname));}
    std::pair<casadi_int, casadi_int> size_in(casadi_int ind) const;
    std::pair<casadi_int, casadi_int> size_in(const std::string& iname) const {
      return size_in(index_in(iname));
    }
    ///@}

    ///@{
    /** \brief Get the shape of an output: rows, columns or both

        Same conventions as size_in, size1_in and size2_in.

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 3)
        >>> f = Function("f", [x], [x.T, SX(2, 2)])
        >>> f.size_out(0), f.size_out("o1"), f.size1_out(0), f.size2_out(1)
        ((1, 3), (2, 2), 1, 2)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 3);
        >> f = Function('f', {x}, {x', SX(2, 2)});
        >> f.size_out(0), f.size_out('o1')
        1 3
        2 2
        >> fprintf('%d %d\n', f.size1_out(0), f.size2_out(1))
        1 2
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 3);
        julia> f = ca.Function("f", [x], [x', ca.SX(2, 2)]);
        julia> ca.size_out(f, 0), ca.size_out(f, "o1"), ca.size1_out(f, 0), ca.size2_out(f, 1)
        ((1, 3), (2, 2), 1, 2)
        \enddoctest
        \endlang

        \sa size_in, sparsity_out

        \identifier{1vb} */
    casadi_int size1_out(casadi_int ind) const;
    casadi_int size1_out(const std::string& oname) const {return size1_out(index_out(oname));}
    casadi_int size2_out(casadi_int ind) const;
    casadi_int size2_out(const std::string& oname) const {return size2_out(index_out(oname));}
    std::pair<casadi_int, casadi_int> size_out(casadi_int ind) const;
    std::pair<casadi_int, casadi_int> size_out(const std::string& oname) const {
      return size_out(index_out(oname));
    }
    ///@}

    ///@{
    /** \brief Get the number of structural nonzeros of one input, or of all inputs

        Without argument, the sum over all inputs. The nonzero count is what numerical
        evaluation reads for the input; it is smaller than numel_in for sparse inputs.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x", 2), SX.sym("D", Sparsity.diag(3))], [SX(1)])
        >>> f.nnz_in(0), f.nnz_in(1), f.nnz_in()
        (2, 3, 5)
        >>> f.numel_in(0), f.numel_in(1), f.numel_in()
        (2, 9, 11)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x', 2), SX.sym('D', Sparsity.diag(3))}, {SX(1)});
        >> fprintf('%d %d %d\n', f.nnz_in(0), f.nnz_in(1), f.nnz_in())
        2 3 5
        >> fprintf('%d %d %d\n', f.numel_in(0), f.numel_in(1), f.numel_in())
        2 9 11
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> D = ca.SX.sym("D", ca.Sparsity.diag(3));
        julia> f = ca.Function("f", [ca.SX.sym("x", 2), D], [ca.SX(1)]);
        julia> ca.nnz_in(f, 0), ca.nnz_in(f, 1), ca.nnz_in(f)
        (2, 3, 5)
        julia> ca.numel_in(f, 0), ca.numel_in(f, 1), ca.numel_in(f)
        (2, 9, 11)
        \enddoctest
        \endlang

        \sa nnz_out, numel_in, sparsity_in

        \identifier{1vc} */
    casadi_int nnz_in() const;
    casadi_int nnz_in(casadi_int ind) const;
    casadi_int nnz_in(const std::string& iname) const {return nnz_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief Get the number of structural nonzeros of one output, or of all outputs

        Without argument, the sum over all outputs.

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 2)
        >>> f = Function("f", [x], [x, diag(x), SX(1, 1)])
        >>> f.nnz_out(0), f.nnz_out(1), f.nnz_out(2), f.nnz_out()
        (2, 2, 0, 4)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 2);
        >> f = Function('f', {x}, {x, diag(x), SX(1, 1)});
        >> fprintf('%d %d %d %d\n', f.nnz_out(0), f.nnz_out(1), f.nnz_out(2), f.nnz_out())
        2 2 0 4
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 2);
        julia> f = ca.Function("f", [x], [x, ca.diag(x), ca.SX(1, 1)]);
        julia> ca.nnz_out(f, 0), ca.nnz_out(f, 1), ca.nnz_out(f, 2), ca.nnz_out(f)
        (2, 2, 0, 4)
        \enddoctest
        \endlang

        \sa nnz_in, numel_out, sparsity_out

        \identifier{1vd} */
    casadi_int nnz_out() const;
    casadi_int nnz_out(casadi_int ind) const;
    casadi_int nnz_out(const std::string& oname) const {return nnz_out(index_out(oname));}
    ///@}

    /** \brief Output signal activity induced by a given input activity
     *
     * \p arg is a flat activity mask over all input nonzeros (size nnz_in(),
     * true = active = possibly nonzero); the result is the corresponding mask over
     * all output nonzeros (size nnz_out()). An output that comes out inactive (false)
     * is guaranteed zero for any value of the active inputs. Unlike sparsity
     * propagation this exploits annihilation: an inactive operand of a multiply kills
     * the product.

        \identifier{2ih} */
    std::vector<bool> activity(const std::vector<bool>& arg) const;

    ///@{
    /** \brief Get the number of elements (rows times columns) of one input, or of all

        Counts structural zeros too, unlike nnz_in.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("D", Sparsity.diag(3))], [SX(1)])
        >>> f.numel_in(0), f.nnz_in(0), f.numel_in("i0")
        (9, 3, 9)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('D', Sparsity.diag(3))}, {SX(1)});
        >> fprintf('%d %d %d\n', f.numel_in(0), f.nnz_in(0), f.numel_in('i0'))
        9 3 9
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("D", ca.Sparsity.diag(3))], [ca.SX(1)]);
        julia> ca.numel_in(f, 0), ca.nnz_in(f, 0), ca.numel_in(f, "i0")
        (9, 3, 9)
        \enddoctest
        \endlang

        \sa nnz_in, numel_out

        \identifier{1ve} */
    casadi_int numel_in() const;
    casadi_int numel_in(casadi_int ind) const;
    casadi_int numel_in(const std::string& iname) const {return numel_in(index_in(iname));}
    ///@}

    ///@{
    /** \brief Get the number of elements (rows times columns) of one output, or of all

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 2)
        >>> f = Function("f", [x], [SX(2, 2), x])
        >>> f.numel_out(0), f.numel_out(1), f.numel_out()
        (4, 2, 6)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 2);
        >> f = Function('f', {x}, {SX(2, 2), x});
        >> fprintf('%d %d %d\n', f.numel_out(0), f.numel_out(1), f.numel_out())
        4 2 6
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 2);
        julia> f = ca.Function("f", [x], [ca.SX(2, 2), x]);
        julia> ca.numel_out(f, 0), ca.numel_out(f, 1), ca.numel_out(f)
        (4, 2, 6)
        \enddoctest
        \endlang

        \sa nnz_out, numel_in

        \identifier{1vf} */
    casadi_int numel_out() const;
    casadi_int numel_out(casadi_int ind) const;
    casadi_int numel_out(const std::string& oname) const {return numel_out(index_out(oname));}
    ///@}

    /** \brief Get the names of all inputs, in order

        Names come from the constructor; when none were given they are i0, i1, ...
        Solver objects have fixed schemes, e.g. x0, p, lbx, ubx, lbg, ubg, lam_x0, lam_g0
        for nlpsol.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> Function("f", [x, y], [x+y]).name_in()
        ['i0', 'i1']
        >>> f = Function("f", [x, y], [x+y], ["x", "y"], ["r"])
        >>> f.name_in(), f.name_in(1), f.name_out()
        (['x', 'y'], 'y', ['r'])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> Function('f', {x, y}, {x+y}).name_in()
        i0
        i1
        >> f = Function('f', {x, y}, {x+y}, {'x', 'y'}, {'r'});
        >> f.name_in(), f.name_in(1), f.name_out()
        x
        y
        y
        r
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> ca.name_in(ca.Function("f", [x, y], [x + y]))
        2-element Vector{String}:
         "i0"
         "i1"
        julia> f = ca.Function("f", [x, y], [x + y], ["x", "y"], ["r"]);
        julia> ca.name_in(f), ca.name_in(f, 1), ca.name_out(f)
        (["x", "y"], "y", ["r"])
        \enddoctest
        \endlang

        \sa name_out, index_in, has_in

        \identifier{1vg} */
    const std::vector<std::string>& name_in() const;

    /** \brief Get the names of all outputs, in order

        Defaults are o0, o1, ... when the constructor was given no names.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> Function("f", [x], [x, 2*x]).name_out()
        ['o0', 'o1']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> Function('f', {x}, {x, 2*x}).name_out()
        o0
        o1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.name_out(ca.Function("f", [x], [x, 2x]))
        2-element Vector{String}:
         "o0"
         "o1"
        \enddoctest
        \endlang

        \sa name_in, index_out, has_out

        \identifier{1vh} */
    const std::vector<std::string>& name_out() const;

    /** \brief Get the name of input \p ind; out of range raises an error

        \identifier{1vi} */
    const std::string& name_in(casadi_int ind) const;

    /** \brief Get the name of output \p ind; out of range raises an error

        \identifier{1vj} */
    const std::string& name_out(casadi_int ind) const;

    /** \brief Get the index of the input with a given name

        Raises an error listing the available names when no input is called \p name.
        Use has_in to test without raising.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x"), SX.sym("y")], [SX(0)], ["x", "y"], ["r"])
        >>> f.index_in("y"), f.index_out("r")
        (1, 0)
        >>> try:
        ...   f.index_in("z")
        ... except RuntimeError as e:
        ...   print('could not find entry "z". Available names are: [x, y]' in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x'), SX.sym('y')}, {SX(0)}, {'x', 'y'}, {'r'});
        >> fprintf('%d %d\n', f.index_in('y'), f.index_out('r'))
        1 0
        >> try, f.index_in('z'); catch e, disp(contains(e.message, 'find entry "z"')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> f = ca.Function("f", [x, y], [ca.SX(0)], ["x", "y"], ["r"]);
        julia> ca.index_in(f, "y"), ca.index_out(f, "r")
        (1, 0)
        julia> try
                   ca.index_in(f, "z")
               catch e
                   msg = sprint(showerror, e)
                   println(occursin("could not find entry \"z\". Available names are: [x, y]", msg))
               end
        true
        \enddoctest
        \endlang

        \sa index_out, has_in, name_in

        \identifier{1vk} */
    casadi_int index_in(const std::string &name) const;

    /** \brief Get the index of the output with a given name; unknown names raise an error

        \sa index_in, has_out, name_out

        \identifier{1vl} */
    casadi_int index_out(const std::string &name) const;

    /** \brief Check whether an input with the given name exists

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)], ["x"], ["r"])
        >>> f.has_in("x"), f.has_in("r"), f.has_out("r"), f.has_out("x")
        (True, False, True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)}, {'x'}, {'r'});
        >> fprintf('%d %d %d %d\n', f.has_in('x'), f.has_in('r'), f.has_out('r'), f.has_out('x'))
        1 0 1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)], ["x"], ["r"]);
        julia> ca.has_in(f, "x"), ca.has_in(f, "r"), ca.has_out(f, "r"), ca.has_out(f, "x")
        (true, false, true, false)
        \enddoctest
        \endlang

        \sa has_out, index_in

        \identifier{2c9} */
    bool has_in(const std::string &name) const;
    /** \brief Check whether an output with the given name exists

        \sa has_in, index_out

        \identifier{2ca} */
    bool has_out(const std::string &name) const;

    /** \brief Get the value used for input \p ind when a call by name omits it

        When calling with a dict (or keyword arguments), every input that is not given is
        filled with this scalar, broadcast to the input shape. For functions built from
        expressions it is 0 unless option "default_in" was set at construction. Solvers
        define their own: nlpsol uses -inf for lbx and lbg, inf for ubx and ubg, 0 for the
        rest.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x + y], ["x", "y"], ["r"])
        >>> f.default_in(0), f(x=3)
        (0.0, {'r': DM([3, 3])})
        >>> g = Function("g", [x, y], [x + y], ["x", "y"], ["r"], {"default_in": [1, 2]})
        >>> g.default_in(1), g(x=3), g(y=[0, 1])
        (2.0, {'r': DM([5, 5])}, {'r': DM([1, 2])})
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x + y}, {'x', 'y'}, {'r'});
        >> disp(f.default_in(0)), res = f('x', 3); disp(res.r)
        0
        [3, 3]
        >> opts = struct('default_in', [1, 2]);
        >> g = Function('g', {x, y}, {x + y}, {'x', 'y'}, {'r'}, opts);
        >> disp(g.default_in(1)), res = g('x', 3); disp(res.r), res = g('y', [0; 1]); disp(res.r)
        2
        [5, 5]
        [1, 2]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x + y], ["x", "y"], ["r"]);
        julia> ca.default_in(f, 0), f(x=3)["r"]
        (0.0, [3, 3])
        julia> opts = Dict("default_in" => [1, 2]);
        julia> g = ca.Function("g", [x, y], [x + y], ["x", "y"], ["r"], opts);
        julia> ca.default_in(g, 1), g(x=3)["r"], g(y=[0, 1])["r"]
        (2.0, [5, 5], [1, 2])
        \enddoctest
        \endlang

        \sa max_in, min_in, nominal_in, convert_in

        \identifier{1vm} */
    double default_in(casadi_int ind) const;

    /** \brief Get the largest admissible value of input \p ind (inf unless a solver

        restricts it)

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.max_in(0), f.min_in(0)
        (inf, -inf)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> fprintf('%g %g\n', f.max_in(0), f.min_in(0))
        Inf -Inf
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.max_in(f, 0), ca.min_in(f, 0)
        (Inf, -Inf)
        \enddoctest
        \endlang

        \sa min_in, default_in

        \identifier{1vn} */
    double max_in(casadi_int ind) const;

    /** \brief Get the smallest admissible value of input \p ind (-inf unless a solver

        restricts it)

        \sa max_in, default_in

        \identifier{1vo} */
    double min_in(casadi_int ind) const;

    /** \brief Get the nominal value of each nonzero of input \p ind

        A list of length nnz_in(ind), all ones unless the class defines scaling.

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 2)
        >>> f = Function("f", [x], [vertcat(x, 1)])
        >>> f.nominal_in(0), f.nominal_out(0)
        ([1.0, 1.0], [1.0, 1.0, 1.0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 2);
        >> f = Function('f', {x}, {vertcat(x, 1)});
        >> f.nominal_in(0), f.nominal_out(0)
        1 1
        1 1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 2);
        julia> f = ca.Function("f", [x], [[x; 1]]);
        julia> ca.nominal_in(f, 0), ca.nominal_out(f, 0)
        ([1.0, 1.0], [1.0, 1.0, 1.0])
        \enddoctest
        \endlang

        \sa nominal_out, default_in

        \identifier{1vp} */
    std::vector<double> nominal_in(casadi_int ind) const;

    /** \brief Get the nominal value of each nonzero of output \p ind; all ones by default

        \sa nominal_in

        \identifier{1vq} */
    std::vector<double> nominal_out(casadi_int ind) const;

    /** \brief Get the sparsity pattern of an input, by index or name

        Numerical arguments are projected onto this pattern before evaluation: entries
        outside it are dropped silently. Symbolic inputs keep the pattern of the expression
        passed to the constructor.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("D", Sparsity.diag(2))], [SX(0)], ["D"], ["r"])
        >>> f.sparsity_in(0)
        Sparsity(2x2,2nz)
        >>> f.sparsity_in("D").is_diag()
        True
        >>> f.sparsity_out(0)
        Sparsity(1x1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('D', Sparsity.diag(2))}, {SX(0)}, {'D'}, {'r'});
        >> f.sparsity_in(0)
        2x2,2nz
        >> f.sparsity_in('D').is_diag()
        1
        >> f.sparsity_out(0)
        1x1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> D = ca.SX.sym("D", ca.Sparsity.diag(2));
        julia> f = ca.Function("f", [D], [ca.SX(0)], ["D"], ["r"]);
        julia> ca.sparsity_in(f, 0)
        Sparsity(2x2,2nz)
        julia> ca.is_diag(ca.sparsity_in(f, "D"))
        true
        julia> ca.sparsity_out(f, 0)
        Sparsity(1x1)
        \enddoctest
        \endlang

        \sa sparsity_out, size_in, nnz_in

        \identifier{1vr} */
    /// @{
    const Sparsity& sparsity_in(casadi_int ind) const;
    const Sparsity& sparsity_in(const std::string& iname) const;
    /// @}

    /** \brief Get the sparsity pattern of an output, by index or name

        Results of a numerical call are DM objects with exactly this pattern.

        \sa sparsity_in, size_out, nnz_out

        \identifier{1vs} */
    /// @{
    const Sparsity& sparsity_out(casadi_int ind) const;
    const Sparsity& sparsity_out(const std::string& iname) const;
    /// @}

    /** \brief Check which inputs or outputs are marked differentiable

        All true by default; set with options "is_diff_in" and "is_diff_out". A
        non-differentiable input or output is skipped by jacobian, forward, reverse and
        factory: the corresponding derivative blocks are all-zero sparse. The signature
        printout marks them with a '#'. Without argument, the list for all inputs or
        outputs is returned.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); n = SX.sym("n")
        >>> f = Function("f", [x, n], [x*n], ["x", "n"], ["r"], {"is_diff_in": [True, False]})
        >>> f
        Function(f:(x,#n)->(r) SXFunction)
        >>> f.is_diff_in(), f.is_diff_in(1), f.is_diff_out()
        ([True, False], False, [True])
        >>> f.jacobian()
        Function(jac_f:(x,#n,out_r[1x1,0nz])->(jac_r_x,#jac_r_n[1x1,0nz]) SXFunction)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); n = SX.sym('n');
        >> opts = struct('is_diff_in', [true, false]);
        >> f = Function('f', {x, n}, {x*n}, {'x', 'n'}, {'r'}, opts); disp(f)
        f:(x,#n)->(r) SXFunction
        >> f.is_diff_in(), f.is_diff_in(1), f.is_diff_out()
        1 0
        0
        1
        >> f.jacobian()
        jac_f:(x,#n,out_r[1x1,0nz])->(jac_r_x,#jac_r_n[1x1,0nz]) SXFunction
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); n = ca.SX.sym("n");
        julia> opts = Dict("is_diff_in" => [true, false]);
        julia> f = ca.Function("f", [x, n], [x*n], ["x", "n"], ["r"], opts);
        julia> f
        Function(f:(x,#n)->(r) SXFunction)
        julia> ca.is_diff_in(f), ca.is_diff_in(f, 1), ca.is_diff_out(f)
        (Bool[1, 0], false, Bool[1])
        julia> ca.jacobian(f)
        Function(jac_f:(x,#n,out_r[1x1,0nz])->(jac_r_x,#jac_r_n[1x1,0nz]) SXFunction)
        \enddoctest
        \endlang

        \sa jacobian

        \identifier{1vt} */
    /// @{
    bool is_diff_in(casadi_int ind) const;
    bool is_diff_out(casadi_int ind) const;
    std::vector<bool> is_diff_in() const;
    std::vector<bool> is_diff_out() const;
    /// @}

    // A linear combination of inputs
    typedef std::map<std::string, std::vector<std::string> > AuxOut;

    /** \brief Build a function computing derivatives by naming inputs and outputs

        Only defined for SXFunction and MXFunction (other classes raise an error).
        The new function takes the inputs listed in \p s_in and returns the outputs
        listed in \p s_out, in that order. Each entry is either an existing name or a
        derivative request built from existing names with ':' separators:
        - inputs: "fwd:x" a forward seed for input x, "adj:r" an adjoint seed for output
          r, "lam:r" a multiplier for output r (used by \p aux)
        - outputs: "fwd:r" forward sensitivity, "adj:x" adjoint sensitivity, "jac:r:x"
          the Jacobian of r with respect to x, "grad:r:x" its gradient (r must be scalar),
          "hess:r:x:y" the Hessian of scalar r with respect to x and y
        - an output may be prefixed with "triu:", "tril:" or "sym:" to return only the
          upper or lower triangle of a Hessian, or to symmetrize it

        The output names of the result have ':' replaced by '_' (jac:r:x becomes jac_r_x).
        Non-differentiable inputs and outputs (is_diff_in/is_diff_out) give all-zero
        blocks. An input that is needed by an output but not listed in \p s_in is
        treated as zero. An unknown name raises an error listing the available names.

        \param name Name of the new function
        \param s_in Input requests, see above
        \param s_out Output requests, see above
        \param aux Extra scalar outputs defined as linear combinations: {"gamma": ["f",
          "g"]} defines gamma = lam_f'*f + lam_g'*g, so that "grad:gamma:x" and
          "hess:gamma:x:x" are the Lagrangian gradient and Hessian with inputs "lam:f"
          and "lam:g"
        \param opts "helper_options" and "final_options" (dicts for the intermediate and
          the returned function) and any Function option for the result

        \lang{python}
        \doctest
        >>> x = SX.sym("x", 2); p = SX.sym("p")
        >>> f = Function("f", [x, p], [p*dot(x, x), x[0]*x[1]], ["x", "p"], ["f", "g"])
        >>> g = f.factory("g", ["x", "p"], ["jac:f:x", "hess:f:x:x", "jac:g:x"])
        >>> g
        Function(g:(x[2],p)->(jac_f_x[1x2],hess_f_x_x[2x2,2nz],jac_g_x[1x2]) SXFunction)
        >>> g([1, 2], 3)
        (DM([[6, 12]]), DM(
        [[6, 00],
         [00, 6]]), DM([[2, 1]]))
        >>> L = f.factory("L", ["x", "p", "lam:f", "lam:g"], ["grad:gamma:x"],
        ...               {"gamma": ["f", "g"]})
        >>> L, L([1, 2], 3, 1, 2)
        (Function(L:(x[2],p,lam_f,lam_g)->(grad_gamma_x[2]) SXFunction), DM([10, 14]))
        >>> h = f.factory("h", ["x", "p", "fwd:x", "fwd:p"], ["fwd:f"])
        >>> h([1, 2], 3, [1, 0], 0)
        DM(6)
        >>> try:
        ...   f.factory("bad", ["x"], ["jac:q:x"])
        ... except RuntimeError as e:
        ...   print('Cannot process "q" as output. Available: f,g.' in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x', 2); p = SX.sym('p');
        >> f = Function('f', {x, p}, {p*dot(x, x), x(1)*x(2)}, {'x', 'p'}, {'f', 'g'});
        >> g = f.factory('g', {'x', 'p'}, {'jac:f:x', 'hess:f:x:x', 'jac:g:x'}); disp(g)
        g:(x[2],p)->(jac_f_x[1x2],hess_f_x_x[2x2,2nz],jac_g_x[1x2]) SXFunction
        >> [jac_f_x, hess_f_x_x, jac_g_x] = g([1; 2], 3);
        >> jac_f_x, hess_f_x_x, jac_g_x
        [[6, 12]]
        [[6, 00],
         [00, 6]]
        [[2, 1]]
        >> aux = struct('gamma', {{'f', 'g'}});
        >> L = f.factory('L', {'x', 'p', 'lam:f', 'lam:g'}, {'grad:gamma:x'}, aux);
        >> L, L([1; 2], 3, 1, 2)
        L:(x[2],p,lam_f,lam_g)->(grad_gamma_x[2]) SXFunction
        [10, 14]
        >> h = f.factory('h', {'x', 'p', 'fwd:x', 'fwd:p'}, {'fwd:f'});
        >> h([1; 2], 3, [1; 0], 0)
        6
        >> try, f.factory('bad', {'x'}, {'jac:q:x'}); catch e, msg = e.message; end
        >> contains(msg, 'Cannot process "q" as output. Available: f,g.')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 2); p = ca.SX.sym("p");
        julia> f = ca.Function("f", [x, p], [p * ca.dot(x, x), x[1] * x[2]],
                               ["x", "p"], ["f", "g"]);
        julia> g = ca.factory(f, "g", ["x", "p"], ["jac:f:x", "hess:f:x:x", "jac:g:x"]);
        julia> g
        Function(g:(x[2],p)->(jac_f_x[1x2],hess_f_x_x[2x2,2nz],jac_g_x[1x2]) SXFunction)
        julia> jac_f_x, hess_f_x_x, jac_g_x = g([1, 2], 3);
        julia> jac_f_x
        DM([[6, 12]])
        julia> hess_f_x_x
        DM(
        [[6, 00],
         [00, 6]])
        julia> jac_g_x
        DM([[2, 1]])
        julia> L = ca.factory(f, "L", ["x", "p", "lam:f", "lam:g"], ["grad:gamma:x"],
                             Dict("gamma" => ["f", "g"]));
        julia> println(L); L([1, 2], 3, 1, 2)
        L:(x[2],p,lam_f,lam_g)->(grad_gamma_x[2]) SXFunction
        DM([10, 14])
        julia> h = ca.factory(f, "h", ["x", "p", "fwd:x", "fwd:p"], ["fwd:f"]);
        julia> h([1, 2], 3, [1, 0], 0)
        DM(6)
        julia> try
                   ca.factory(f, "bad", ["x"], ["jac:q:x"])
               catch e
                   msg = sprint(showerror, e)
                   println(occursin("Cannot process \"q\" as output. Available: f,g.", msg))
               end
        true
        \enddoctest
        \endlang

        \sa jacobian, forward, reverse, which_depends

        \identifier{2mi} */
    Function factory(const std::string& name,
                     const std::vector<std::string>& s_in,
                     const std::vector<std::string>& s_out,
                     const AuxOut& aux=AuxOut(),
                     const Dict& opts=Dict()) const;

    /** \brief Get the user-supplied problem function behind a solver object

        For nlpsol, qpsol, rootfinder and integrator objects this is the Function with the
        problem expressions (e.g. nlp:(x,p)->(f,g)); for other classes an error is raised.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> opts = {"qpsol": "qrqp", "print_header": False, "print_iteration": False,
        ...         "print_time": False, "qpsol_options": {"print_header": False}}
        >>> S = nlpsol("S", "sqpmethod", {"x": x, "f": x**2}, opts)
        >>> S.oracle()
        Function(nlp:(x,p[])->(f,g[0]) SXFunction)
        >>> try:
        ...   S.oracle().oracle()
        ... except RuntimeError as e:
        ...   print("'oracle' not defined for SXFunction" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> opts = struct('qpsol', 'qrqp', 'print_header', false, 'print_iteration', false);
        >> opts.print_time = false; opts.qpsol_options = struct('print_header', false);
        >> S = nlpsol('S', 'sqpmethod', struct('x', x, 'f', x^2), opts);
        >> S.oracle()
        nlp:(x,p[])->(f,g[0]) SXFunction
        >> try, S.oracle().oracle(); catch e, disp(contains(e.message, 'not defined')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> opts = Dict("qpsol" => "qrqp", "print_header" => false, "print_iteration" => false,
                           "print_time" => false, "qpsol_options" => Dict("print_header" => false));
        julia> S = ca.nlpsol("S", "sqpmethod", Dict("x" => x, "f" => x^2), opts);
        julia> ca.oracle(S)
        Function(nlp:(x,p[])->(f,g[0]) SXFunction)
        julia> try
                   ca.oracle(ca.oracle(S))
               catch e
                   println(occursin("'oracle' not defined for SXFunction", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa get_function, find_functions

        \identifier{1vu} */
    Function oracle() const;

    /** \brief Wrap in an MXFunction whose graph is a single call to this function

        The wrapper has the same signature and names, carries over the derivative options
        (ad_weight, ad_weight_sp, max_num_dir, is_diff_in/out, derivative_of) and is
        cached by name. Classes without their own derivative or slice implementation use
        it internally; it is also a way to give any function a plain MX graph.

        \param name Name of the wrapper; defaults to "wrap_" + name()

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [sin(x)], ["x"], ["y"])
        >>> w = f.wrap()
        >>> w, w.n_instructions(), w.get_function()
        (Function(wrap_f:(x)->(y) MXFunction), 3, ['f'])
        >>> f.wrap("g").name()
        'g'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {sin(x)}, {'x'}, {'y'});
        >> w = f.wrap(); disp(w), disp(w.n_instructions()), disp(w.get_function())
        wrap_f:(x)->(y) MXFunction
        3
        f
        >> f.wrap('g').name()
        g
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [sin(x)], ["x"], ["y"]);
        julia> w = ca.wrap(f);
        julia> println(w); ca.n_instructions(w), ca.get_function(w)
        wrap_f:(x)->(y) MXFunction
        (3, ["f"])
        julia> ca.name(ca.wrap(f, "g"))
        "g"
        \enddoctest
        \endlang

        \sa wrap_as_needed, expand

        \identifier{1vv} */
    /// @{
    Function wrap() const;
    Function wrap(const std::string& name) const;
    /// @}

    /** \brief Wrap in an MXFunction constructed with the given options

        Like wrap, but the wrapper is built with \p opts (merged with the derivative
        options of this function) and is not cached. The named overload returns this
        function itself when \p opts is empty and \p name equals its current name; the
        unnamed overload always wraps, under the name "wrap_" + name().

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [sin(x)])
        >>> f.wrap_as_needed("f", {}).class_name(), f.wrap_as_needed({}).class_name()
        ('SXFunction', 'MXFunction')
        >>> f.wrap_as_needed("g", {"ad_weight": 1}).name()
        'g'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {sin(x)});
        >> f.wrap_as_needed('f', struct()).class_name()
        SXFunction
        >> f.wrap_as_needed(struct()).class_name()
        MXFunction
        >> f.wrap_as_needed('g', struct('ad_weight', 1)).name()
        g
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [sin(x)]);
        julia> ca.class_name(ca.wrap_as_needed(f, "f", Dict()))
        "SXFunction"
        julia> ca.class_name(ca.wrap_as_needed(f, Dict()))
        "MXFunction"
        julia> ca.name(ca.wrap_as_needed(f, "g", Dict("ad_weight" => 1)))
        "g"
        \enddoctest
        \endlang

        \sa wrap

        \identifier{1vw} */
    /// @{
    Function wrap_as_needed(const Dict& opts) const;
    Function wrap_as_needed(const std::string& name, const Dict& opts) const;
    /// @}

    /** \brief Find which entries of an input the outputs depend on, linearly or not

        Returns one boolean per nonzero of input \p s_in: true when at least one of the
        outputs \p s_out depends on that entry (\p order 1) or depends on it nonlinearly
        (\p order 2). With \p tr true the relation is flipped: one boolean per nonzero of
        the vertically concatenated outputs, true when it depends on any entry of \p s_in.
        The analysis is structural (it follows the graph, not the values).

        \param s_in Name of one input
        \param s_out Names of the outputs to inspect
        \param order 1 (any dependence) or 2 (nonlinear dependence); other values raise
        \param tr Transpose the relation, see above

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y[0], y[1]**2 + x, 3], ["x", "y"], ["a", "b", "c"])
        >>> f.which_depends("y", ["a", "b"])
        [True, True]
        >>> f.which_depends("y", ["a", "b"], 2)
        [False, True]
        >>> f.which_depends("y", ["a", "b", "c"], 1, True)
        [True, True, False]
        >>> f.which_depends("x", ["a"], 2)
        [False]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y(1), y(2)^2 + x, 3}, {'x', 'y'}, {'a', 'b', 'c'});
        >> f.which_depends('y', {'a', 'b'})
        1 1
        >> f.which_depends('y', {'a', 'b'}, 2)
        0 1
        >> f.which_depends('y', {'a', 'b', 'c'}, 1, true)
        1 1 0
        >> f.which_depends('x', {'a'}, 2)
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y[1], y[2]^2 + x, 3], ["x", "y"], ["a", "b", "c"]);
        julia> ca.which_depends(f, "y", ["a", "b"])
        2-element Vector{Bool}:
         1
         1
        julia> ca.which_depends(f, "y", ["a", "b"], 2)
        2-element Vector{Bool}:
         0
         1
        julia> ca.which_depends(f, "y", ["a", "b", "c"], 1, true)
        3-element Vector{Bool}:
         1
         1
         0
        julia> ca.which_depends(f, "x", ["a"], 2)
        1-element Vector{Bool}:
         0
        \enddoctest
        \endlang

        \sa factory, jac_sparsity

        \identifier{1vx} */
    std::vector<bool> which_depends(const std::string& s_in,
                                    const std::vector<std::string>& s_out,
                                    casadi_int order=1, bool tr=false) const;

    /** \brief Print the name and shape of every input and output

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x", 2), SX.sym("p")], [SX(1, 3)], ["x", "p"], ["r"])
        >>> f.print_dimensions()
        Number of inputs: 2
        Input 0 ("x"): 2x1
        Input 1 ("p"): 1x1
        Number of outputs: 1
        Output 0 ("r"): 1x3
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x', 2), SX.sym('p')}, {SX(1, 3)}, {'x', 'p'}, {'r'});
        >> f.print_dimensions()
        Number of inputs: 2
        Input 0 ("x"): 2x1
        Input 1 ("p"): 1x1
        Number of outputs: 1
        Output 0 ("r"): 1x3
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x", 2); p = ca.SX.sym("p");
        julia> f = ca.Function("f", [x, p], [ca.SX(1, 3)], ["x", "p"], ["r"]);
        julia> ca.print_dimensions(f)
         Number of inputs: 2
          Input 0 ("x"): 2x1
          Input 1 ("p"): 1x1
         Number of outputs: 1
          Output 0 ("r"): 1x3
        \enddoctest
        \endlang

        \sa size_in, sparsity_in

        \identifier{1vy} */
    void print_dimensions(std::ostream &stream=casadi::uout()) const;

    /** \brief Print the name, type and description of every option of this class

        The list depends on the class (SXFunction, MXFunction, a solver plugin, ...).
        Options are set at construction; see change_option for the few that can be
        changed afterwards.

        \sa print_option, has_option

        \identifier{1vz} */
    void print_options(std::ostream &stream=casadi::uout()) const;

    /** \brief Print the type and description of one option; unknown names print a

        notice

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.print_option("verbose")
        > "verbose"          [OT_BOOL]      "Verbose evaluation -- for debugging"
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> f.print_option('verbose')
        > "verbose"          [OT_BOOL]      "Verbose evaluation -- for debugging"
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.print_option(f, "verbose")
        > "verbose"          [OT_BOOL]      "Verbose evaluation -- for debugging"
        \enddoctest
        \endlang

        \sa print_options, has_option

        \identifier{1w0} */
    void print_option(const std::string &name, std::ostream &stream = casadi::uout()) const;

    /** \brief Check whether this class accepts an option with the given name

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.has_option("live_variables"), f.has_option("ipopt")
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> fprintf('%d %d\n', f.has_option('live_variables'), f.has_option('ipopt'))
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.has_option(f, "live_variables"), ca.has_option(f, "ipopt")
        (true, false)
        \enddoctest
        \endlang

        \sa print_options, change_option

        \identifier{1w1} */
    bool has_option(const std::string &option_name) const;

    /** \brief Change a debugging option on an existing object

        Only options that do not change the numerical result can be changed:
        "verbose", "regularity_check", "print_in", "print_out", "dump", "dump_in",
        "dump_out", "dump_trace", "dump_dir", "dump_format", "print_time", "record_time"
        and "print_instructions". Any other existing option raises "cannot be changed";
        an unknown name raises "does not exist".

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.change_option("verbose", False)
        >>> try:
        ...   f.change_option("cse", True)
        ... except RuntimeError as e:
        ...   print("Option 'cse' cannot be changed" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> f.change_option('verbose', false)
        >> try, f.change_option('cse', true); catch e, msg = e.message; end
        >> contains(msg, 'Option ''cse'' cannot be changed')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.change_option(f, "verbose", false)
        julia> try
                   ca.change_option(f, "cse", true)
               catch e
                   println(occursin("Option 'cse' cannot be changed", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa has_option

        \identifier{1w2} */
    void change_option(const std::string& option_name, const GenericType& option_value);

    /** \brief Reset to zero the counter that numbers the files written by dump_in,

        dump_out and dump_trace

        \identifier{2dy} */
    void reset_dump_count();

    /** \brief Check whether the derivative functions read the nondifferentiated outputs

        The functions returned by jacobian, forward and reverse always accept the
        original outputs as extra inputs (out_r, ...). For most classes those inputs are
        ignored and this returns false; it is true for classes that reuse the outputs,
        for example to avoid recomputing a solver solution.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.uses_output()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> f.uses_output()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.uses_output(f)
        false
        \enddoctest
        \endlang

        \sa jacobian, forward, reverse

        \identifier{1w3} */
    bool uses_output() const;

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Jacobian of output \p oind w.r.t. input \p iind, followed by

        all outputs; use factory with "jac:out:in" instead

        \identifier{1w4} */
    Function jacobian_old(casadi_int iind, casadi_int oind) const;

    /** \brief [DEPRECATED] Hessian and gradient of scalar output \p oind w.r.t. input

        \p iind, followed by all outputs; use factory with "hess:out:in:in" instead

        \identifier{1w5} */
    Function hessian_old(casadi_int iind, casadi_int oind) const;

    ///@{
    /// [DEPRECATED] Sparsity of the Jacobian of output \p oind w.r.t. input \p iind;
    /// use jac_sparsity(oind, iind) (note the reversed argument order)
    const Sparsity sparsity_jac(casadi_int iind, casadi_int oind,
                                bool compact=false, bool symmetric=false) const;
    const Sparsity sparsity_jac(const std::string &iind, casadi_int oind=0,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(index_in(iind), oind, compact, symmetric);
    }
    const Sparsity sparsity_jac(casadi_int iind, const std::string &oind,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(iind, index_out(oind), compact, symmetric);
    }
    const Sparsity sparsity_jac(const std::string &iind, const std::string &oind,
                                bool compact=false, bool symmetric=false) const {
      return sparsity_jac(index_in(iind), index_out(oind), compact, symmetric);
    }
    ///@}
#endif // WITH_DEPRECATED_FEATURES

  /** \brief Get a function computing all Jacobian blocks at once

      The result is named jac_<name> and takes all inputs followed by all outputs
      (named out_<output>); the output values are only read when uses_output() is true
      and may be passed as 0 otherwise. It returns one block per (output, input) pair,
      ordered by output first: for f:(x,y)->(r,s) the outputs are jac_r_x, jac_r_y,
      jac_s_x, jac_s_y. Blocks involving a non-differentiable input or output are
      all-zero sparse. Each block is a size_out(o)-times-size_in(i) matrix where the
      output and input are vectorized column-major; blocks that are structurally zero
      have an all-zero sparsity pattern.

      The result is cached: calling jacobian again returns the same object as long as it
      is alive. The derivative itself is built by symbolic AD for SX and MX functions, by
      the solver's own sensitivity method for solver objects, or raises "Derivatives
      cannot be calculated" for classes without any derivative support. For a single
      block, factory("...", ..., ["jac:r:x"]) avoids computing the others.

      \lang{python}
      \doctest
      >>> x = SX.sym("x"); y = SX.sym("y", 2)
      >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
      >>> J = f.jacobian()
      >>> print(J)
      jac_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz])->(jac_r_x[2],jac_r_y[2x2,2nz],...) SXFunction
      >>> J(1, [2, 3], 0, 0)
      (DM([2, 3]), DM(
      [[1, 00],
       [00, 1]]), DM(0.540302), DM([[00, 00]]))
      >>> f.jacobian().name_out()
      ['jac_r_x', 'jac_r_y', 'jac_s_x', 'jac_s_y']
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> x = SX.sym('x'); y = SX.sym('y', 2);
      >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
      >> J = f.jacobian(); disp(J)
      jac_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz])->(jac_r_x[2],jac_r_y[2x2,2nz],...) SXFunction
      >> [jac_r_x, jac_r_y, jac_s_x, jac_s_y] = J(1, [2; 3], 0, 0);
      >> jac_r_x, jac_r_y, jac_s_x, jac_s_y
      [2, 3]
      [[1, 00],
       [00, 1]]
      0.540302
      [[00, 00]]
      >> J.name_out()
      jac_r_x
      jac_r_y
      jac_s_x
      jac_s_y
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
      julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
      julia> J = ca.jacobian(f);
      julia> ca.name(J), ca.n_in(J), ca.n_out(J)
      ("jac_f", 4, 4)
      julia> jac_r_x, jac_r_y, jac_s_x, jac_s_y = J(1, [2, 3], 0, 0);
      julia> jac_r_x, jac_s_x, jac_s_y
      ([2, 3], 0.540302, [[00, 00]])
      julia> jac_r_y
      DM(
      [[1, 00],
       [00, 1]])
      julia> ca.name_out(J)
      4-element Vector{String}:
       "jac_r_x"
       "jac_r_y"
       "jac_s_x"
       "jac_s_y"
      \enddoctest
      \endlang

      \sa factory, forward, reverse, jac_sparsity, uses_output

      \identifier{1w6} */
    Function jacobian() const;

    ///@{
    /** \brief Evaluate with a list or a dict of arguments; the result has the same form

        This is the general evaluation entry point; the call operator (f(...) in Python,
        C++ and MATLAB) is a shorthand for it. Argument lists must have exactly n_in()
        entries. A dict maps input names to values; missing inputs take default_in,
        unknown names raise an error. The result is a list of n_out() values, or a dict
        keyed by output name.

        The argument type selects the evaluation mode:
        - DM (or anything convertible to DM, such as numbers and numeric arrays): numerical
          evaluation, returns DM
        - SX: symbolic evaluation producing an SX expression; an MX function is
          expanded on the fly (fails if it contains non-expandable nodes)
        - MX: symbolic evaluation producing a call node f(args){k} in the MX graph,
          unless the function is inlined

        Shapes: each argument must have the input shape, or be a scalar (broadcast to all
        entries), a transposed vector, an empty matrix (treated as all zeros), or
        N-by-(P*M) for an N-by-M input (P evaluations, outputs stacked horizontally;
        other arguments are repeated). Anything else raises an error that lists these
        rules. Sparsity is projected: entries outside sparsity_in are ignored. When all
        inputs are scalar and the arguments are matrices of one common shape, the function
        is applied elementwise. Numerical evaluation of a function with free variables
        raises "Cannot evaluate ... since variables [...] are free".

        \param arg List or dict of arguments, all of the same type
        \param res Output: list or dict of results
        \param always_inline MX only: expand the call into the caller's graph instead of
          creating a call node (SX functions are then inlined as scalar operations);
          combined with never_inline raises an error
        \param never_inline SX only: keep the call as an opaque call node inside the
          SX expression instead of inlining it

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
        >>> f.call([2, [3, 4]])
        [DM([6, 8]), DM(0.909297)]
        >>> f.call({"x": 2})
        {'r': DM([0, 0]), 's': DM(0.909297)}
        >>> f(2, [3, 4])
        (DM([6, 8]), DM(0.909297))
        >>> f(x=2, y=[3, 4])
        {'r': DM([6, 8]), 's': DM(0.909297)}
        >>> f(2, DM([[3, 5], [4, 6]]))
        (DM(
        [[6, 10],
         [8, 12]]), DM([[0.909297, 0.909297]]))
        >>> f(DM(0, 0), 3)
        (DM([0, 0]), DM(0))
        >>> z = MX.sym("z"); w = MX.sym("w", 2)
        >>> f.call([z, w])
        [MX(f(z, w){0}), MX(f(z, w){1})]
        >>> f.call([z, w], True)
        [MX(vertcat((z*w[0]), (z*w[1]))), MX(sin(z))]
        >>> a = SX.sym("a")
        >>> f.call([a, [1, 2]])
        [SX([a, (2.*a)]), SX(sin(a))]
        >>> try:
        ...   f(1)
        ... except RuntimeError as e:
        ...   print("Incorrect number of inputs: Expected 2, got 1" in str(e))
        True
        >>> try:
        ...   f(1, [1, 2, 3])
        ... except RuntimeError as e:
        ...   print("Input 1 (y) has mismatching shape. Got 3-by-1" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
        >> res = f.call({2, [3; 4]}); disp(res{1}), disp(res{2})
        [6, 8]
        0.909297
        >> res = f.call(struct('x', 2)); disp(res.r), disp(res.s)
        [0, 0]
        0.909297
        >> [r, s] = f(2, [3; 4]); disp(r), disp(s)
        [6, 8]
        0.909297
        >> res = f('x', 2, 'y', [3; 4]); disp(res.r), disp(res.s)
        [6, 8]
        0.909297
        >> [r, s] = f(2, DM([3, 5; 4, 6])); disp(r), disp(s)
        [[6, 10],
         [8, 12]]
        [[0.909297, 0.909297]]
        >> [r, s] = f(DM(0, 0), 3); disp(r), disp(s)
        [0, 0]
        0
        >> z = MX.sym('z'); w = MX.sym('w', 2);
        >> res = f.call({z, w}); disp(res{1}), disp(res{2})
        f(z, w){0}
        f(z, w){1}
        >> res = f.call({z, w}, true); disp(res{1}), disp(res{2})
        vertcat((z*w[0]), (z*w[1]))
        sin(z)
        >> a = SX.sym('a'); res = f.call({a, [1; 2]}); disp(res{1}), disp(res{2})
        [a, (2.*a)]
        sin(a)
        >> try, f(1); catch e, disp(contains(e.message, 'Expected 2, got 1')); end
        1
        >> try, f(1, [1; 2; 3]); catch e, disp(contains(e.message, 'Got 3-by-1')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
        julia> ca.call(f, [2, [3, 4]])
        2-element Vector{DM}:
         DM([6, 8])
         DM(0.909297)
        julia> ca.call(f, Dict("x" => 2))
        Dict{String, DM} with 2 entries:
          "r" => [0, 0]
          "s" => 0.909297
        julia> f(2, [3, 4])
        ([6, 8], 0.909297)
        julia> f(x=2, y=[3, 4])
        Dict{String, DM} with 2 entries:
          "r" => [6, 8]
          "s" => 0.909297
        julia> r, s = f(2, [3 5; 4 6]);
        julia> r
        DM(
        [[6, 10],
         [8, 12]])
        julia> s
        DM([[0.909297, 0.909297]])
        julia> f(ca.DM(0, 0), 3)
        ([0, 0], 0)
        julia> z = ca.MX.sym("z"); w = ca.MX.sym("w", 2);
        julia> ca.call(f, [z, w])
        2-element Vector{MX}:
         MX(f(z, w){0})
         MX(f(z, w){1})
        julia> ca.call(f, [z, w], true)
        2-element Vector{MX}:
         MX(vertcat((z*w[0]), (z*w[1])))
         MX(sin(z))
        julia> a = ca.SX.sym("a");
        julia> ca.call(f, [a, [1, 2]])
        2-element Vector{SX}:
         SX([a, (2.*a)])
         SX(sin(a))
        julia> try
                   f(1)
               catch e
                   msg = sprint(showerror, e)
                   println(occursin("Incorrect number of inputs: Expected 2, got 1", msg))
               end
        true
        julia> try
                   f(1, [1, 2, 3])
               catch e
                   msg = sprint(showerror, e)
                   println(occursin("Input 1 (y) has mismatching shape. Got 3-by-1", msg))
               end
        true
        \enddoctest
        \endlang

        \sa default_in, convert_in, map, mapsum

        \identifier{1w7} */
    void call(const std::vector<DM> &arg, std::vector<DM>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const std::vector<SX> &arg, std::vector<SX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const std::vector<MX> &arg, std::vector<MX>& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const DMDict& arg, DMDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const SXDict& arg, SXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    void call(const MXDict& arg, MXDict& SWIG_OUTPUT(res),
              bool always_inline=false, bool never_inline=false) const;
    ///@}

#ifndef SWIG
    /// Check if same as another function
    bool operator==(const Function& f) const;

    ///@{
    /// Functor shorthand for evaluation
    std::vector<DM> operator()(const std::vector<DM>& arg) const;
    std::vector<SX> operator()(const std::vector<SX>& arg) const;
    std::vector<MX> operator()(const std::vector<MX>& arg) const;
    const DMDict operator()(const DMDict& arg) const;
    const SXDict operator()(const SXDict& arg) const;
    const MXDict operator()(const MXDict& arg) const;
    ///@}

    ///@{
    /** \brief Evaluate with temporary memory allocation

        \identifier{1w8} */
    void operator()(std::vector<const double*> arg, std::vector<double*> res) const;
    void operator()(std::vector<const bvec_t*> arg, std::vector<bvec_t*> res) const;
    void operator()(std::vector<const SXElem*> arg, std::vector<SXElem*> res) const;
    template<typename D> void call_gen(std::vector<const D*> arg, std::vector<D*> res) const;
    ///@}

    ///@{
    /** \brief Supported arguments for numerical evaluation and converters

        \identifier{1w9} */
    typedef const std::vector<std::vector<double>>& VecArg;
    std::vector<const double*> buf_in(VecArg arg) const;
    typedef std::vector<std::vector<double>>& VecRes;
    std::vector<double*> buf_out(VecRes res) const;
    typedef std::vector<std::vector<double>*> VPrRes;
    std::vector<double*> buf_out(VPrRes res) const;

    typedef const std::map<std::string, std::vector<double>>& MapArg;
    std::vector<const double*> buf_in(MapArg arg) const;
    typedef std::map<std::string, std::vector<double>>& MapRes;
    std::vector<double*> buf_out(MapRes res) const;
    typedef std::map<std::string, std::vector<double>*> MPrRes;
    std::vector<double*> buf_out(MPrRes res) const;
    ///@}

    ///@{
    /** \brief Numerical evaluation

        \identifier{1wa} */
    void operator()(VecArg arg, VecRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MapRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, VPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(VecArg arg, MPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }

    void operator()(MapArg arg, VecRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MapRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, VPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    void operator()(MapArg arg, MPrRes res) const { (*this)(buf_in(arg), buf_out(res)); }
    ///@}

    ///@{
    /// Functor shorthand for evaluation, single argument (only C++)
    std::vector<DM> operator()(const DM& arg0) const {
      return operator()(std::vector<DM>{arg0});
    }
    std::vector<SX> operator()(const SX& arg0) const {
      return operator()(std::vector<SX>{arg0});
    }
    std::vector<MX> operator()(const MX& arg0) const {
      return operator()(std::vector<MX>{arg0});
    }
    ///@}

    /** \brief Evaluate memory-less, numerically

        \identifier{1wb} */
    int operator()(const double** arg, double** res,
        casadi_int* iw, double* w, int mem) const;

    /** \brief Evaluate numerically with checkout/release

        \identifier{1wc} */
    int operator()(const double** arg, double** res,
        casadi_int* iw, double* w) const;

    /** \brief Evaluate memory-less SXElem

        Same syntax as the double version, allowing use in templated code

        \identifier{1wd} */
    int operator()(const SXElem** arg, SXElem** res,
        casadi_int* iw, SXElem* w, int mem=0) const;

    /** \brief  Propagate sparsity forward

        \identifier{1we} */
    int operator()(const bvec_t** arg, bvec_t** res,
        casadi_int* iw, bvec_t* w, int mem=0) const;

    /** \brief Propagate signal activity forward (bit set = active (possibly nonzero))

        \identifier{2ii} */
    int eval_activity(const bvec_t** arg, bvec_t** res,
        casadi_int* iw, bvec_t* w, int mem=0) const;

    /** \brief  Propagate sparsity backward

        \identifier{1wf} */
    int rev(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, int mem=0) const;

    /** \brief Propagate sparsity backward with temporary memory allocation

        \identifier{1wg} */
    int rev(std::vector<bvec_t*> arg, std::vector<bvec_t*> res) const;

#endif // SWIG

    /** \brief Evaluate on horizontally stacked MX arguments and sum the results

        Each argument is either the input itself or n copies stacked horizontally
        (N-by-(n*M) for an N-by-M input); the largest such n sets the number of
        evaluations and non-stacked arguments are shared across evaluations. Returns the
        sum over the n evaluations of each output, as MX expressions calling an
        internal map-and-reduce function.

        \param x List of n_in() MX arguments
        \param parallelization "serial" only; other values raise an error

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> f = Function("f", [x, y], [x*y, x - y], ["x", "y"], ["p", "d"])
        >>> p, d = f.mapsum([MX.sym("a"), MX.sym("b", 1, 3)])
        >>> p
        MX(mapsum(a, b){0})
        >>> f.mapsum([2, DM([[1, 2, 3]])])
        [MX(mapsum(2, [[1, 2, 3]]){0}), MX(mapsum(2, [[1, 2, 3]]){1})]
        >>> evalf(f.mapsum([2, DM([[1, 2, 3]])])[0])
        DM(12)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> f = Function('f', {x, y}, {x*y, x - y}, {'x', 'y'}, {'p', 'd'});
        >> res = f.mapsum({MX.sym('a'), MX.sym('b', 1, 3)}); disp(res{1})
        mapsum(a, b){0}
        >> res = f.mapsum({2, DM([1, 2, 3])}); disp(res{1}), disp(res{2})
        mapsum(2, [[1, 2, 3]]){0}
        mapsum(2, [[1, 2, 3]]){1}
        >> evalf(res{1})
        12
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> f = ca.Function("f", [x, y], [x*y, x - y], ["x", "y"], ["p", "d"]);
        julia> p, d = ca.mapsum(f, [ca.MX.sym("a"), ca.MX.sym("b", 1, 3)]);
        julia> p
        MX(mapsum(a, b){0})
        julia> ca.mapsum(f, [2, ca.DM([1 2 3])])
        2-element Vector{MX}:
         MX(mapsum(2, [[1, 2, 3]]){0})
         MX(mapsum(2, [[1, 2, 3]]){1})
        julia> ca.evalf(ca.mapsum(f, [2, ca.DM([1 2 3])])[1])
        DM(12)
        \enddoctest
        \endlang

        \sa map, call

        \identifier{1wh} */
    std::vector<MX> mapsum(const std::vector<MX > &x,
                           const std::string& parallelization="serial") const;

    ///@{
    /** \brief Create a function that applies this function N times, feeding outputs back

        as inputs

        The first \p n_accum inputs are accumulators: in every step they take the value
        of the first \p n_accum outputs of the previous step. For a step
        f:(x, u)->(x_next, y) the result is
        F:(x0, U)->(X, Y) with U = horzcat(u0, ..., u_(N-1)), X = horzcat(x1, ..., x_N),
        Y = horzcat(y0, ..., y_(N-1)) and x_(k+1), y_k = f(x_k, u_k). Non-accumulated
        inputs and all outputs are stacked horizontally over the N steps; the accumulator
        inputs keep their original shape. An accumulator output may be wider than its
        input (e.g. 1-by-2 feeding a 1-by-1 input): only the last columns are fed back.

        Compared with writing the loop in the host language the construction is much
        faster and the graph is a tower of nested unrolled loops of at most \p base
        iterations (option "base", default 10; -1 unrolls completely). This keeps
        memory use proportional to base rather than N at a small evaluation cost.

        The overloads differ in how accumulators are selected: a count \p n_accum
        (default 1), or explicit index/name lists \p accum_in and \p accum_out of equal
        length that may pick arbitrary inputs and outputs. The unnamed overload names the
        result "mapaccum_" + name(). fold is mapaccum returning only the final
        accumulator state.

        \param name Name of the result
        \param N Number of steps, must be positive
        \param n_accum Number of leading inputs/outputs that are accumulated; at most
          min(n_in, n_out)
        \param opts "base" plus options for the resulting MXFunction

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); u = SX.sym("u")
        >>> step = Function("step", [x, u], [x + u, x*u], ["x", "u"], ["x_next", "y"])
        >>> F = step.mapaccum(3)
        >>> F
        Function(mapaccum_step:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction)
        >>> F(1, DM([[1, 2, 3]]))
        (DM([[2, 4, 7]]), DM([[1, 4, 12]]))
        >>> step.mapaccum("F", 3, ["x"], ["x_next"])
        Function(F:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction)
        >>> step.fold(3)
        Function(fold_step:(x,u[1x3])->(x_next,y[1x3]) MXFunction)
        >>> step.fold(3)(1, DM([[1, 2, 3]]))
        (DM(7), DM([[1, 4, 12]]))
        >>> try:
        ...   step.mapaccum(0)
        ... except RuntimeError as e:
        ...   print("mapaccum: N must be positive" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); u = SX.sym('u');
        >> step = Function('step', {x, u}, {x + u, x*u}, {'x', 'u'}, {'x_next', 'y'});
        >> F = step.mapaccum(3); disp(F)
        mapaccum_step:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction
        >> [X, Y] = F(1, DM([1, 2, 3])); disp(X), disp(Y)
        [[2, 4, 7]]
        [[1, 4, 12]]
        >> step.mapaccum('F', 3, {'x'}, {'x_next'})
        F:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction
        >> G = step.fold(3); disp(G)
        fold_step:(x,u[1x3])->(x_next,y[1x3]) MXFunction
        >> [X, Y] = G(1, DM([1, 2, 3])); disp(X), disp(Y)
        7
        [[1, 4, 12]]
        >> try, step.mapaccum(0); catch e, disp(contains(e.message, 'N must be positive')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); u = ca.SX.sym("u");
        julia> step = ca.Function("step", [x, u], [x + u, x*u], ["x", "u"], ["x_next", "y"]);
        julia> F = ca.mapaccum(step, 3);
        julia> F
        Function(mapaccum_step:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction)
        julia> F(1, [1 2 3])
        ([[2, 4, 7]], [[1, 4, 12]])
        julia> ca.mapaccum(step, "F", 3, ["x"], ["x_next"])
        Function(F:(x,u[1x3])->(x_next[1x3],y[1x3]) MXFunction)
        julia> ca.fold(step, 3)
        Function(fold_step:(x,u[1x3])->(x_next,y[1x3]) MXFunction)
        julia> ca.fold(step, 3)(1, [1 2 3])
        (7, [[1, 4, 12]])
        julia> try
                   ca.mapaccum(step, 0)
               catch e
                   println(occursin("mapaccum: N must be positive", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa map, fold

        \identifier{1wi} */
    Function mapaccum(const std::string& name, casadi_int N, const Dict& opts = Dict()) const;
    Function mapaccum(const std::string& name, casadi_int N, casadi_int n_accum,
                      const Dict& opts = Dict()) const;
    Function mapaccum(const std::string& name, casadi_int n,
                      const std::vector<casadi_int>& accum_in,
                      const std::vector<casadi_int>& accum_out,
                      const Dict& opts=Dict()) const;
    Function mapaccum(const std::string& name, casadi_int n,
                      const std::vector<std::string>& accum_in,
                      const std::vector<std::string>& accum_out,
                      const Dict& opts=Dict()) const;
    Function mapaccum(casadi_int N, const Dict& opts = Dict()) const;
    Function fold(casadi_int N, const Dict& opts = Dict()) const;
    ///@}

    /** \brief Create a function that evaluates this function n times on stacked arguments

        For f:(a, p)->(s) the result is F:(A, P)->(S) with A = horzcat(a0, ..., a_(n-1)),
        the same for P, and S = horzcat(f(a0, p0), ..., f(a_(n-1), p_(n-1))). Every input
        and output is widened to n times its number of columns; when calling, a
        non-widened argument is repeated for all evaluations. The graph size and
        construction time are independent of n. For n == 1 the function itself is
        returned; n < 1 raises "Degenerate map operation".

        \param n Number of evaluations
        \param parallelization How the evaluations are run:
          "serial" (a Map node evaluated in a loop; cached, so repeated calls return the
          same object), "openmp" (OpenMP parallel loop, requires an OpenMP build),
          "thread" (std::thread per evaluation), "unroll" (n explicit calls in an
          MXFunction), "inline" (like unroll with the calls inlined). Other values raise
          an error.
        \param max_num_threads With "thread" or "openmp": split the n evaluations into
          at most this many parallel chunks, each a serial map

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); p = SX.sym("p")
        >>> f = Function("f", [x, p], [p*x**2], ["x", "p"], ["y"])
        >>> F = f.map(3)
        >>> F
        Function(map3_f:(x[1x3],p[1x3])->(y[1x3]) Map)
        >>> F(DM([[1, 2, 3]]), 2)
        DM([[2, 8, 18]])
        >>> f.map(4, "thread", 2)
        Function(threadmap2_map2_f:(x[1x4],p[1x4])->(y[1x4]) ThreadMap)
        >>> f.map(3, "unroll")
        Function(f_3:(x[1x3],p[1x3])->(y[1x3]) MXFunction)
        >>> f.map(1) is f, f.map(1).name()
        (False, 'f')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); p = SX.sym('p');
        >> f = Function('f', {x, p}, {p*x^2}, {'x', 'p'}, {'y'});
        >> F = f.map(3); disp(F)
        map3_f:(x[1x3],p[1x3])->(y[1x3]) Map
        >> F(DM([1, 2, 3]), 2)
        [[2, 8, 18]]
        >> f.map(4, 'thread', 2)
        threadmap2_map2_f:(x[1x4],p[1x4])->(y[1x4]) ThreadMap
        >> f.map(3, 'unroll')
        f_3:(x[1x3],p[1x3])->(y[1x3]) MXFunction
        >> f.map(1).name()
        f
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
        julia> f = ca.Function("f", [x, p], [p*x^2], ["x", "p"], ["y"]);
        julia> F = f.map(3);
        julia> F
        Function(map3_f:(x[1x3],p[1x3])->(y[1x3]) Map)
        julia> F([1 2 3], 2)
        DM([[2, 8, 18]])
        julia> f.map(4, "thread", 2)
        Function(threadmap2_map2_f:(x[1x4],p[1x4])->(y[1x4]) ThreadMap)
        julia> f.map(3, "unroll")
        Function(f_3:(x[1x3],p[1x3])->(y[1x3]) MXFunction)
        julia> f.map(1) === f, ca.name(f.map(1))
        (false, "f")
        \enddoctest
        \endlang

        \sa mapaccum, mapsum, call

        \identifier{1wj} */
    Function map(casadi_int n, const std::string& parallelization="serial") const;
    Function map(casadi_int n, const std::string& parallelization,
      casadi_int max_num_threads) const;

    ///@{
    /** \brief Map where some inputs are shared and some outputs are summed

        Inputs listed in \p reduce_in are not widened: the same value is used in all n
        evaluations. Outputs listed in \p reduce_out are summed over the evaluations
        instead of stacked. The index and name overloads build an MXFunction around a
        plain map (\p parallelization as in map; \p opts is currently ignored). The
        boolean-mask overload creates a dedicated serial MapSum node wrapped in an
        MXFunction built with \p opts; masks must have length n_in and n_out (an empty
        \p reduce_out means no reduction).

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); p = SX.sym("p")
        >>> f = Function("f", [x, p], [p*x, x], ["x", "p"], ["y", "z"])
        >>> F = f.map("F", "serial", 3, ["p"], ["y"])
        >>> F
        Function(F:(x[1x3],p)->(y,z[1x3]) MXFunction)
        >>> F(DM([[1, 2, 3]]), 2)
        (DM(12), DM([[1, 2, 3]]))
        >>> f.map(3, [False, True], [True, False])(DM([[1, 2, 3]]), 2)
        (DM(12), DM([[1, 2, 3]]))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); p = SX.sym('p');
        >> f = Function('f', {x, p}, {p*x, x}, {'x', 'p'}, {'y', 'z'});
        >> F = f.map('F', 'serial', 3, {'p'}, {'y'}); disp(F)
        F:(x[1x3],p)->(y,z[1x3]) MXFunction
        >> [y, z] = F(DM([1, 2, 3]), 2); disp(y), disp(z)
        12
        [[1, 2, 3]]
        >> G = f.map(3, [false, true], [true, false]);
        >> [y, z] = G(DM([1, 2, 3]), 2); disp(y), disp(z)
        12
        [[1, 2, 3]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
        julia> f = ca.Function("f", [x, p], [p*x, x], ["x", "p"], ["y", "z"]);
        julia> F = f.map("F", "serial", 3, ["p"], ["y"]);
        julia> F
        Function(F:(x[1x3],p)->(y,z[1x3]) MXFunction)
        julia> F([1 2 3], 2)
        (12, [[1, 2, 3]])
        julia> f.map(3, [false, true], [true, false])([1 2 3], 2)
        (12, [[1, 2, 3]])
        \enddoctest
        \endlang

        \sa map, mapsum

        \identifier{1wk} */
    Function map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<casadi_int>& reduce_in,
      const std::vector<casadi_int>& reduce_out,
      const Dict& opts=Dict()) const;
    Function map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<std::string>& reduce_in,
      const std::vector<std::string>& reduce_out,
      const Dict& opts=Dict()) const;
    Function map(casadi_int n,
      const std::vector<bool>& reduce_in,
      const std::vector<bool>& reduce_out=std::vector<bool>(),
      const Dict& opts=Dict()) const;
    ///@}

    /** \brief Create a function with a subset or reordering of the inputs and outputs

        The expressions are unchanged; only the lists of inputs and outputs are
        rearranged. Indices may be repeated or omitted. An omitted input that an output
        depends on makes the result have free variables (construction then fails unless
        "allow_free" is in \p opts). Classes other than SXFunction and MXFunction are
        wrapped first.

        \param name Name of the result
        \param order_in Indices of the original inputs to keep, in the new order
        \param order_out Indices of the original outputs to keep, in the new order
        \param opts Options for the new function

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> f = Function("f", [x, y], [x + y, x - y, 2*x], ["x", "y"], ["a", "b", "c"])
        >>> g = f.slice("g", [1, 0], [2, 0])
        >>> g, g(1, 2)
        (Function(g:(y,x)->(c,a) SXFunction), (DM(4), DM(3)))
        >>> f.slice("h", [0], [2])
        Function(h:(x)->(c) SXFunction)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> f = Function('f', {x, y}, {x + y, x - y, 2*x}, {'x', 'y'}, {'a', 'b', 'c'});
        >> g = f.slice('g', [1, 0], [2, 0]); disp(g)
        g:(y,x)->(c,a) SXFunction
        >> [c, a] = g(1, 2); disp(c), disp(a)
        4
        3
        >> f.slice('h', 0, 2)
        h:(x)->(c) SXFunction
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> f = ca.Function("f", [x, y], [x + y, x - y, 2x], ["x", "y"], ["a", "b", "c"]);
        julia> g = ca.slice(f, "g", [1, 0], [2, 0]);
        julia> println(g); g(1, 2)
        g:(y,x)->(c,a) SXFunction
        (4, 3)
        julia> ca.slice(f, "h", [0], [2])
        Function(h:(x)->(c) SXFunction)
        \enddoctest
        \endlang

        \sa factory

        \identifier{1wl} */
    Function slice(const std::string& name, const std::vector<casadi_int>& order_in,
                   const std::vector<casadi_int>& order_out, const Dict& opts=Dict()) const;

    /** \brief Create a function that selects one of several functions by an index input

        The result takes an index as first input followed by the inputs of the branch
        functions. At evaluation the index is truncated to an integer k; branch f[k] is
        evaluated when 0 <= k < len(f), otherwise \p f_def. All branches must have the
        same number of inputs and outputs; the signature (names, shapes) is taken from
        \p f_def and arguments are projected when a branch has a different sparsity.
        Derivatives are the derivatives of the selected branch. Symbolic MX evaluation
        gives a call node; unlike if_else on expressions, only the selected branch is
        evaluated numerically.

        \param name Name of the result
        \param f Branch functions, selected by index 0, 1, ...
        \param f_def Default branch for indices out of range
        \param opts Function options

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f0 = Function("f0", [x], [x]); f1 = Function("f1", [x], [x**2])
        >>> fdef = Function("fdef", [x], [-x])
        >>> sw = Function.conditional("sw", [f0, f1], fdef)
        >>> sw
        Function(sw:(i0,i1)->(o0) Switch)
        >>> sw(0, 3), sw(1, 3), sw(1.9, 3), sw(2, 3), sw(-1, 3)
        (DM(3), DM(9), DM(9), DM(-3), DM(-3))
        >>> sw(MX.sym("k"), MX.sym("v"))
        MX(sw(k, v){0})
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f0 = Function('f0', {x}, {x}); f1 = Function('f1', {x}, {x^2});
        >> fdef = Function('fdef', {x}, {-x});
        >> sw = Function.conditional('sw', {f0, f1}, fdef); disp(sw)
        sw:(i0,i1)->(o0) Switch
        >> [sw(0, 3), sw(1, 3), sw(1.9, 3), sw(2, 3), sw(-1, 3)]
        [[3, 9, 9, -3, -3]]
        >> sw(MX.sym('k'), MX.sym('v'))
        sw(k, v){0}
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f0 = ca.Function("f0", [x], [x]); f1 = ca.Function("f1", [x], [x^2]);
        julia> fdef = ca.Function("fdef", [x], [-x]);
        julia> sw = ca.Function.conditional("sw", [f0, f1], fdef);
        julia> sw
        Function(sw:(i0,i1)->(o0) Switch)
        julia> sw(0, 3), sw(1, 3), sw(1.9, 3), sw(2, 3), sw(-1, 3)
        (3, 9, 9, -3, -3)
        julia> sw(ca.MX.sym("k"), ca.MX.sym("v"))
        MX(sw(k, v){0})
        \enddoctest
        \endlang

        \sa if_else

        \identifier{1wm} */
    static Function conditional(const std::string& name, const std::vector<Function>& f,
                                const Function& f_def, const Dict& opts=Dict());

    /** \brief Create a function that calls \p f when its first input is nonzero and

        returns zeros otherwise

        Shorthand for if_else with a zero-returning function of the same signature as
        the false branch; evaluating with a zero condition skips \p f entirely.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [x**2])
        >>> c = Function.conditional("c", f)
        >>> c, c(0, 3), c(1, 3)
        (Function(c:(i0,i1)->(o0) Switch), DM(0), DM(9))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {x^2});
        >> c = Function.conditional('c', f); disp(c)
        c:(i0,i1)->(o0) Switch
        >> [c(0, 3), c(1, 3)]
        [[0, 9]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [x^2]);
        julia> c = ca.Function.conditional("c", f);
        julia> println(c); c(0, 3), c(1, 3)
        c:(i0,i1)->(o0) Switch
        (0, 9)
        \enddoctest
        \endlang

        \sa if_else

        \identifier{1wn} */
    static Function conditional(const std::string& name,
      const Function& f, const Dict& opts=Dict());

    /** \brief Create a function evaluating a B-spline with known coefficients

        The result maps a point of dimension len(knots) (one column vector input) to the
        m-dimensional spline value. Coefficients are stored flat in \p coeffs with the
        value dimension m fastest, then the first spline dimension, etc.; the number of
        coefficients per dimension is len(knots[d]) - degree[d] - 1. A mismatch between
        the lengths of \p knots and \p degree raises an error. For data-driven splines
        or general interpolation see interpolant.

        \param name Name of the result
        \param knots Knot vector per dimension, non-decreasing, with repeated end knots
        \param coeffs Flat coefficient list
        \param degree Polynomial degree per dimension
        \param m Dimension of the spline value
        \param opts "lookup_mode" (per dimension: "linear", "exact", "binary") and
          Function options

        \lang{python}
        \doctest
        >>> f = Function.bspline("f", [[0, 0, 0, 1, 2, 2, 2]], [0, 1, 2, 3], [2])
        >>> f
        Function(f:(i0)->(o0) MXFunction)
        >>> f(0), f(0.5), f(2)
        (DM(0), DM(0.875), DM(3))
        >>> g = Function.bspline("g", [[0, 0, 1, 1], [0, 0, 1, 1]], [0, 1, 2, 3], [1, 1])
        >>> g([0.5, 0.5])
        DM(1.5)
        >>> try:
        ...   Function.bspline("f", [[0, 0, 0, 1, 2, 2, 2]], [0, 1, 2, 3], [2, 2])
        ... except RuntimeError as e:
        ...   print("must match knot list length" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> knots = {[0, 0, 0, 1, 2, 2, 2]};
        >> f = Function.bspline('f', knots, [0, 1, 2, 3], 2); disp(f)
        f:(i0)->(o0) MXFunction
        >> [f(0), f(0.5), f(2)]
        [[0, 0.875, 3]]
        >> g = Function.bspline('g', {[0, 0, 1, 1], [0, 0, 1, 1]}, [0, 1, 2, 3], [1, 1]);
        >> g([0.5; 0.5])
        1.5
        >> try, Function.bspline('f', knots, [0, 1, 2, 3], [2, 2]); catch e, msg = e.message; end
        >> contains(msg, 'must match knot list length')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function.bspline("f", [[0, 0, 0, 1, 2, 2, 2]], [0, 1, 2, 3], [2]);
        julia> f
        Function(f:(i0)->(o0) MXFunction)
        julia> f(0), f(0.5), f(2)
        (0, 0.875, 3)
        julia> g = ca.Function.bspline("g", [[0, 0, 1, 1], [0, 0, 1, 1]], [0, 1, 2, 3], [1, 1]);
        julia> g([0.5, 0.5])
        DM(1.5)
        julia> try
                   ca.Function.bspline("f", [[0, 0, 0, 1, 2, 2, 2]], [0, 1, 2, 3], [2, 2])
               catch e
                   println(occursin("must match knot list length", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa interpolant, MX::bspline

        \identifier{1wo} */
    static Function bspline(const std::string &name,
      const std::vector< std::vector<double> >& knots, const std::vector<double>& coeffs,
      const std::vector<casadi_int>& degree, casadi_int m=1, const Dict& opts=Dict());

    /** \brief Create a function that evaluates \p f_true or \p f_false depending on

        a condition input

        The result takes the condition as first input followed by the inputs of the
        branches; a nonzero condition (truncated to integer) selects \p f_true, zero
        selects \p f_false. Only the selected branch is evaluated numerically. The
        signature is taken from \p f_false (see conditional).

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function.if_else("f", Function("t", [x], [x**2]), Function("e", [x], [x]))
        >>> f
        Function(f:(i0,i1)->(o0) Switch)
        >>> f(1, 3), f(0, 3), f(5, 3)
        (DM(9), DM(3), DM(9))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f_true = Function('t', {x}, {x^2}); f_false = Function('e', {x}, {x});
        >> f = Function.if_else('f', f_true, f_false); disp(f)
        f:(i0,i1)->(o0) Switch
        >> [f(1, 3), f(0, 3), f(5, 3)]
        [[9, 3, 9]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f_true = ca.Function("t", [x], [x^2]); f_false = ca.Function("e", [x], [x]);
        julia> f = ca.Function.if_else("f", f_true, f_false);
        julia> f
        Function(f:(i0,i1)->(o0) Switch)
        julia> f(1, 3), f(0, 3), f(5, 3)
        (9, 3, 9)
        \enddoctest
        \endlang

        \sa conditional

        \identifier{1wp} */
    static Function if_else(const std::string& name, const Function& f_true,
                            const Function& f_false, const Dict& opts=Dict());

    /** \brief Get a function computing \p nfwd forward directional derivatives

        The result is named fwd<nfwd>_<name> and has n_in + n_out + n_in inputs: the
        original inputs, the original outputs (out_<name>, read only when uses_output()
        is true) and the forward seeds fwd_<input>, each with nfwd times the columns of
        the input. Its n_out outputs fwd_<output> are the forward sensitivities, nfwd
        stacked horizontally. Uses the class' own forward mode when available ("enable_forward"),
        finite differences when "enable_fd" is set, or otherwise the full Jacobian; raises
        "Derivatives cannot be calculated" when none applies. Results are cached while
        alive. \p nfwd must be positive.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
        >>> F = f.forward(1)
        >>> print(F)
        fwd1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],fwd_x,fwd_y[2])->(...) SXFunction
        >>> F(1, [2, 3], 0, 0, 1, [0, 0])
        (DM([2, 3]), DM(0.540302))
        >>> print(f.forward(2))
        fwd2_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],fwd_x[1x2],fwd_y[2x2])->(...) SXFunction
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
        >> F = f.forward(1); disp(F)
        fwd1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],fwd_x,fwd_y[2])->(...) SXFunction
        >> [fwd_r, fwd_s] = F(1, [2; 3], 0, 0, 1, [0; 0]); disp(fwd_r), disp(fwd_s)
        [2, 3]
        0.540302
        >> f.forward(2)
        fwd2_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],fwd_x[1x2],fwd_y[2x2])->(...) SXFunction
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
        julia> F = ca.forward(f, 1);
        julia> F
        Function(fwd1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],fwd_x,fwd_y[2])->(fwd_r[2],fwd_s) SXFunction)
        julia> F(1, [2, 3], 0, 0, 1, [0, 0])
        ([2, 3], 0.540302)
        julia> F2 = ca.forward(f, 2);
        julia> ca.name(F2), ca.size_in(F2, 4), ca.size_in(F2, 5)
        ("fwd2_f", (1, 2), (2, 2))
        \enddoctest
        \endlang

        \sa reverse, jacobian, factory

        \identifier{1wq} */
    Function forward(casadi_int nfwd) const;

    /** \brief Get a function computing \p nadj reverse (adjoint) directional derivatives

        The result is named adj<nadj>_<name> and has n_in + n_out + n_out inputs: the
        original inputs, the original outputs (out_<name>) and the adjoint seeds
        adj_<output>, each with nadj times the columns of the output. Its n_in outputs
        adj_<input> are the adjoint sensitivities, nadj stacked horizontally. Falls back
        to the full Jacobian when the class has no reverse mode; raises when it has no
        derivatives at all. Results are cached while alive. \p nadj must be positive.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
        >>> R = f.reverse(1)
        >>> R
        Function(adj1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],adj_r[2],adj_s)->(adj_x,adj_y[2])
        SXFunction)
        >>> R(1, [2, 3], 0, 0, [1, 1], 0)
        (DM(5), DM([1, 1]))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
        >> R = f.reverse(1); disp(R)
        adj1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],adj_r[2],adj_s)->(adj_x,adj_y[2]) SXFunction
        >> [adj_x, adj_y] = R(1, [2; 3], 0, 0, [1; 1], 0); disp(adj_x), disp(adj_y)
        5
        [1, 1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
        julia> R = f.reverse(1);
        julia> R
        Function(adj1_f:(x,y[2],out_r[2x1,0nz],out_s[1x1,0nz],adj_r[2],adj_s)->(adj_x,adj_y[2]) SXFunction)
        julia> R(1, [2, 3], 0, 0, [1, 1], 0)
        (5, [1, 1])
        \enddoctest
        \endlang

        \sa forward, jacobian, factory

        \identifier{1wr} */
    Function reverse(casadi_int nadj) const;

    /** \brief Get the sparsity patterns of all Jacobian blocks

        One pattern per (output, input) pair, ordered by output first, as in jacobian.
        Patterns are computed by sparsity propagation the first time and cached.

        \param compact When true, the pattern is over the nonzeros of output and input
          (nnz_out-by-nnz_in) rather than over all elements (numel_out-by-numel_in)

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); D = SX.sym("D", Sparsity.diag(2))
        >>> f = Function("f", [x, D], [x*D, sin(x)])
        >>> f.jac_sparsity()
        [Sparsity(4x1,2nz), Sparsity(4x4,2nz), Sparsity(1x1), Sparsity(1x4,0nz)]
        >>> f.jac_sparsity(True)
        [Sparsity(2x1), Sparsity(2x2,2nz), Sparsity(1x1), Sparsity(1x2,0nz)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); D = SX.sym('D', Sparsity.diag(2));
        >> f = Function('f', {x, D}, {x*D, sin(x)});
        >> sp = f.jac_sparsity(); disp(numel(sp)), cellfun(@disp, sp)
        4
        4x1,2nz
        4x4,2nz
        1x1
        1x4,0nz
        >> cellfun(@disp, f.jac_sparsity(true))
        2x1
        2x2,2nz
        1x1
        1x2,0nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); D = ca.SX.sym("D", ca.Sparsity.diag(2));
        julia> f = ca.Function("f", [x, D], [x*D, sin(x)]);
        julia> ca.jac_sparsity(f)
        4-element Vector{Sparsity}:
         Sparsity(4x1,2nz)
         Sparsity(4x4,2nz)
         Sparsity(1x1)
         Sparsity(1x4,0nz)
        julia> ca.jac_sparsity(f, true)
        4-element Vector{Sparsity}:
         Sparsity(2x1)
         Sparsity(2x2,2nz)
         Sparsity(1x1)
         Sparsity(1x2,0nz)
        \enddoctest
        \endlang

        \sa jacobian, which_depends

        \identifier{1ws} */
    const std::vector<Sparsity>& jac_sparsity(bool compact = false) const;

    /** \brief Get the sparsity pattern of the Jacobian of output \p oind w.r.t. input

        \p iind

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); D = SX.sym("D", Sparsity.diag(2))
        >>> f = Function("f", [x, D], [x*D, sin(x)])
        >>> f.jac_sparsity(0, 1), f.jac_sparsity(0, 1, True), f.jac_sparsity(1, 1)
        (Sparsity(4x4,2nz), Sparsity(2x2,2nz), Sparsity(1x4,0nz))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); D = SX.sym('D', Sparsity.diag(2));
        >> f = Function('f', {x, D}, {x*D, sin(x)});
        >> f.jac_sparsity(0, 1), f.jac_sparsity(0, 1, true), f.jac_sparsity(1, 1)
        4x4,2nz
        2x2,2nz
        1x4,0nz
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); D = ca.SX.sym("D", ca.Sparsity.diag(2));
        julia> f = ca.Function("f", [x, D], [x*D, sin(x)]);
        julia> ca.jac_sparsity(f, 0, 1), ca.jac_sparsity(f, 0, 1, true), ca.jac_sparsity(f, 1, 1)
        (4x4,2nz, 2x2,2nz, 1x4,0nz)
        \enddoctest
        \endlang

        \sa jacobian

        \identifier{1wt} */
    Sparsity jac_sparsity(casadi_int oind, casadi_int iind, bool compact = false) const;

    /** \brief Generate self-contained C code for the function

        Writes \p fname (a ".c" suffix is added when missing; the base name must satisfy
        check_name) in the current directory and returns the file name. The generated
        file exports, for a function named f, the evaluation routine
        int f(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w,
        int mem) plus f_n_in, f_n_out, f_name_in, f_name_out, f_sparsity_in,
        f_sparsity_out, f_work, f_checkout, f_release, f_incref, f_decref and
        f_default_in. The arrays hold the nonzeros of each input and output in
        column-major order; work vector sizes come from f_work. Functions embedded in the
        graph (solvers, maps, ...) are generated too when they support code generation;
        otherwise an error is raised. The same file can hold several functions via the
        CodeGenerator class.

        \param fname File name, e.g. "f.c" or "f"
        \param opts Code generator options: "with_header" (also write a .h file),
          "main" (add a main() reading inputs from a file written by generate_in),
          "casadi_real" (scalar type, default "double"), "casadi_int" (integer type),
          "with_mem" (simplified memory API), "mex", "cpp", "indent", "prefix",
          "verbose", "with_export", "with_import", "include_math", "infinity", "nan",
          "real_min", "unroll_args", "codegen_scalars", "avoid_stack", "thread_safe",
          "force_canonical"

        \lang{python}
        \doctest
        >>> import os
        >>> x = SX.sym("x"); y = SX.sym("y", 2)
        >>> f = Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"])
        >>> fname = f.generate("f.c")
        >>> fname, os.path.isfile("f.c")
        ('f.c', True)
        >>> src = open("f.c").read()
        >>> print("int f(const casadi_real** arg, casadi_real** res" in src, "f_n_in" in src)
        True True
        >>> os.remove("f.c")
        >>> f.generate("g", {"with_header": True}), sorted(os.listdir(".")).count("g.h")
        ('g.c', 1)
        >>> os.remove("g.c"); os.remove("g.h")
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y', 2);
        >> f = Function('f', {x, y}, {x*y, sin(x)}, {'x', 'y'}, {'r', 's'});
        >> fname = f.generate('f.c'); disp(fname), disp(isfile('f.c'))
        f.c
        1
        >> src = fileread('f.c');
        >> disp(contains(src, 'int f(const casadi_real** arg, casadi_real** res'))
        1
        >> contains(src, 'f_n_in')
        1
        >> delete('f.c');
        >> f.generate('g', struct('with_header', true)), isfile('g.h')
        g.c
        1
        >> delete('g.c'); delete('g.h');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y", 2);
        julia> f = ca.Function("f", [x, y], [x*y, sin(x)], ["x", "y"], ["r", "s"]);
        julia> fname = ca.generate(f, "f.c")
        "f.c"
        julia> isfile("f.c")
        true
        julia> src = read("f.c", String);
        julia> occursin("int f(const casadi_real** arg, casadi_real** res", src)
        true
        julia> occursin("f_n_in", src)
        true
        julia> rm("f.c")
        julia> ca.generate(f, "g", Dict("with_header" => true)), isfile("g.h")
        ("g.c", true)
        julia> rm("g.c"); rm("g.h")
        \enddoctest
        \endlang

        \sa CodeGenerator, generate_in, generate_out, jit, external

        \identifier{1wu} */
    std::string generate(const std::string& fname, const Dict& opts=Dict()) const;

    /** \brief Generate C code into a file named after the function (name() + ".c")

        \identifier{1wv} */
    std::string generate(const Dict& opts=Dict()) const;

    /** \brief Generate C code for the problem functions of a solver object

        Only defined for solver classes built on an oracle (nlpsol, rootfinder,
        integrator, ...): writes code for the helper functions (e.g. nlp_f, nlp_jac_g)
        into \p fname and returns the file name. Other classes raise an error. The
        library compiled from this file can be loaded with external and passed back to
        the solver constructor in place of the symbolic problem.

        \param fname File name
        \param opts Code generator options, as for generate

        \sa generate, external

        \identifier{1ww} */
    std::string generate_dependencies(const std::string& fname, const Dict& opts=Dict()) const;

    /** \brief Write or read a text file with the nonzeros of all inputs

        With arguments, writes the nonzeros of \p arg (one value per line, in input
        order) to \p fname; this is the input format read by code generated with option
        "main". Without arguments, reads such a file back into a list of DM with the
        input sparsities. The arguments follow the shape rules of call (scalars are
        broadcast); a file with the wrong number of values raises a dimension error.

        \lang{python}
        \doctest
        >>> import os
        >>> f = Function("f", [SX.sym("x"), SX.sym("y", 2)], [SX(0)])
        >>> f.generate_in("f_in.txt", [1, [2, 3]])
        >>> open("f_in.txt").read().split()
        ['1.0000000000000000e+00', '2.0000000000000000e+00', '3.0000000000000000e+00']
        >>> f.generate_in("f_in.txt")
        [DM(1), DM([2, 3])]
        >>> os.remove("f_in.txt")
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x'), SX.sym('y', 2)}, {SX(0)});
        >> f.generate_in('f_in.txt', {1, [2; 3]});
        >> fileread('f_in.txt')
        1.0000000000000000e+00
        2.0000000000000000e+00
        3.0000000000000000e+00
        >> arg = f.generate_in('f_in.txt'); disp(arg{1}), disp(arg{2})
        1
        [2, 3]
        >> delete('f_in.txt');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x"), ca.SX.sym("y", 2)], [ca.SX(0)]);
        julia> ca.generate_in(f, "f_in.txt", [1, [2, 3]])
        julia> split(read("f_in.txt", String))
        3-element Vector{SubString{String}}:
         "1.0000000000000000e+00"
         "2.0000000000000000e+00"
         "3.0000000000000000e+00"
        julia> ca.generate_in(f, "f_in.txt")
        2-element Vector{DM}:
         DM(1)
         DM([2, 3])
        julia> rm("f_in.txt")
        \enddoctest
        \endlang

        \sa generate_out, nz_from_in, nz_to_in, generate

        \identifier{1wx} */
    /// @{
    void generate_in(const std::string& fname, const std::vector<DM>& arg);
    std::vector<DM> generate_in(const std::string& fname);
    /// @}

    /** \brief Write or read a text file with the nonzeros of all outputs

        Same format as generate_in, for the outputs: the file written by generated code
        with option "main" can be read back and compared with the result of call.

        \lang{python}
        \doctest
        >>> import os
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [x, 2*x])
        >>> f.generate_out("f_out.txt", f.call([1.5]))
        >>> f.generate_out("f_out.txt")
        [DM(1.5), DM(3)]
        >>> os.remove("f_out.txt")
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {x, 2*x});
        >> f.generate_out('f_out.txt', f.call({1.5}));
        >> res = f.generate_out('f_out.txt'); disp(res{1}), disp(res{2})
        1.5
        3
        >> delete('f_out.txt');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [x, 2x]);
        julia> ca.generate_out(f, "f_out.txt", ca.call(f, [1.5]))
        julia> ca.generate_out(f, "f_out.txt")
        2-element Vector{DM}:
         DM(1.5)
         DM(3)
        julia> rm("f_out.txt")
        \enddoctest
        \endlang

        \sa generate_in, nz_from_out, nz_to_out

        \identifier{1wy} */
    /// @{
    void generate_out(const std::string& fname, const std::vector<DM>& arg);
    std::vector<DM> generate_out(const std::string& fname);
    /// @}

    /** \brief Export the algorithm as source code in another language

        Only "matlab" is supported, and only for SXFunction and MXFunction without free
        variables; other classes or languages raise an error. The code is a plain MATLAB
        function with the same name that evaluates the algorithm on numeric arguments,
        without any CasADi dependency. The file overload writes to \p fname, the other
        returns the code as a string.

        \param lang "matlab"
        \param options "indent_level" (integer)

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [sin(x)])
        >>> f.export_code("matlab").splitlines()[0]
        'function [varargout] = f(varargin)'
        >>> try:
        ...   f.export_code("python")
        ... except RuntimeError as e:
        ...   print("Only matlab language supported" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {sin(x)});
        >> code = strsplit(f.export_code('matlab'), newline); disp(code{1})
        function [varargout] = f(varargin)
        >> try, f.export_code('python'); catch e, disp(contains(e.message, 'Only matlab')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [sin(x)]);
        julia> first(split(ca.export_code(f, "matlab"), "\n"))
        "function [varargout] = f(varargin)"
        julia> try
                   ca.export_code(f, "python")
               catch e
                   println(occursin("Only matlab language supported", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa generate

        \identifier{1wz} */
    ///@{
    void export_code(const std::string& lang,
      const std::string &fname, const Dict& options=Dict()) const;

#ifndef SWIG
    /** \brief Serialize

        \identifier{1x0} */
    void serialize(std::ostream &stream, const Dict& opts=Dict()) const;

    /** \brief Serialize an object

        \identifier{1x1} */
    void serialize(SerializingStream &s) const;
#endif

    /** \brief Serialize to a string that deserialize turns back into an equal function

        Any Function can be serialized, including solvers and functions with embedded
        calls; a just-in-time compiled function is stored according to option
        "jit_serialize". The string is binary-safe text; it is also what the Python
        pickle protocol and MATLAB save use.

        \param opts "debug" (bool) adds consistency markers to the stream

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [sin(x)], ["x"], ["y"])
        >>> s = f.serialize()
        >>> g = Function.deserialize(s)
        >>> g, g(1)
        (Function(f:(x)->(y) SXFunction), DM(0.841471))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {sin(x)}, {'x'}, {'y'});
        >> s = f.serialize(); g = Function.deserialize(s); disp(g), disp(g(1))
        f:(x)->(y) SXFunction
        0.841471
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [sin(x)], ["x"], ["y"]);
        julia> s = ca.serialize(f);
        julia> g = ca.Function.deserialize(s);
        julia> println(g); g(1)
        f:(x)->(y) SXFunction
        DM(0.841471)
        \enddoctest
        \endlang

        \sa deserialize, save, load

        \identifier{1x2} */
    std::string serialize(const Dict& opts=Dict()) const;

    /** \brief Export a Function graph as .html, .dot, or .casadi_viz.
     *
     * HTML and .casadi_viz embed the serialized Function, decoded by the browser viewer;
     * any serializable Function works, and helpers such as nlp_hess_l can be explored.
     * include_functions (default true) shows called and helper Functions in the viewer.
     * DOT is standalone Graphviz of an SX/MX root graph only. Options show_matrix_contents and
     * show_matrix_sizes (both default true) control matrix detail and dimensions.
     * Load dump_trace JSONL files in the viewer to replay numerical evaluation.
     * Options: view (function, expression; default function),
     * direction (LR, RL, TB, BT; expression view only; default TB), viz_js (local Viz.js
     * global bundle to embed; by default uses the viewer package renderer).
     * viewer_url (string) overrides the HTML viewer ESM module URL; default is
     * https://unpkg.com/@casadi/casadi-viz@MAJOR.MINOR/dist/index.js for this CasADi version.

        \identifier{2k3} */
    void export_graph(const std::string& fname, const Dict& opts=Dict()) const;

    /** \brief Return the .casadi_viz JSON bundle without writing a file.

     * Uses the same graph options as the filename overload.

        \identifier{2k4} */
    std::string export_graph(const Dict& opts=Dict()) const;

    /** \brief Save the function to a file that load reads back

        Writes a serialization of the function to \p fname (conventionally with the
        extension .casadi). The file is portable across platforms and languages
        (Python, MATLAB, C++); it is compatible across CasADi versions as far as the
        serialization format allows.

        \param fname File name
        \param opts Serializer options, e.g. "debug"

        \lang{python}
        \doctest
        >>> import os
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [x**2])
        >>> f.save("f.casadi")
        >>> g = Function.load("f.casadi")
        >>> g, g(3)
        (Function(f:(i0)->(o0) SXFunction), DM(9))
        >>> os.remove("f.casadi")
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {x^2});
        >> f.save('f.casadi'); g = Function.load('f.casadi'); disp(g), disp(g(3))
        f:(i0)->(o0) SXFunction
        9
        >> delete('f.casadi');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [x^2]);
        julia> ca.save(f, "f.casadi")
        julia> g = ca.Function.load("f.casadi");
        julia> println(g); g(3)
        f:(i0)->(o0) SXFunction
        DM(9)
        julia> rm("f.casadi")
        \enddoctest
        \endlang

        \sa load, serialize, FileSerializer

        \identifier{240} */
    void save(const std::string &fname, const Dict& opts=Dict()) const;

    std::string export_code(const std::string& lang, const Dict& options=Dict()) const;
#ifndef SWIG
    void export_code(const std::string& lang,
      std::ostream &stream, const Dict& options=Dict()) const;
#endif // SWIG
    ///@}
#ifndef SWIG
    /// \cond INTERNAL
    /// Get a const pointer to the node
    FunctionInternal* get() const;

    /// Get a pointer and typecast
    template<typename T>
    T* get() const {
      T* ret = dynamic_cast<T*>(get());
      casadi_assert_dev(ret!=nullptr);
      return ret;
    }

    /** \brief  Const access functions of the node

        \identifier{1x3} */
    FunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);
    /// \endcond
#endif // SWIG

    /** \brief Get the statistics recorded by the last numerical evaluation

        The content depends on the class: solvers report e.g. "success",
        "return_status", "iter_count"; with option "record_time" every function reports
        "t_wall_total", "t_proc_total" and "n_call_total". Raises an error when the
        function has not been evaluated numerically yet.

        \param mem Memory object (from checkout) whose statistics to return; 0 is the
          memory used by ordinary calls

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [x**2], {"record_time": True})
        >>> try:
        ...   f.stats()
        ... except RuntimeError as e:
        ...   print("first evaluate it numerically" in str(e))
        True
        >>> r = f(2)
        >>> sorted(f.stats().keys())
        ['n_call_total', 't_proc_total', 't_wall_total']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {x^2}, struct('record_time', true));
        >> try, f.stats(); catch e, disp(contains(e.message, 'first evaluate it numerically')); end
        1
        >> r = f(2); disp(strjoin(sort(fieldnames(f.stats()))', ' '))
        n_call_total t_proc_total t_wall_total
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [x^2], Dict("record_time" => true));
        julia> try
                   ca.stats(f)
               catch e
                   println(occursin("first evaluate it numerically", sprint(showerror, e)))
               end
        true
        julia> r = f(2);
        julia> sort(collect(keys(ca.stats(f))))
        3-element Vector{String}:
         "n_call_total"
         "t_proc_total"
         "t_wall_total"
        \enddoctest
        \endlang

        \sa checkout

        \identifier{2mj} */
    Dict stats(int mem=0) const;

    ///@{
    /** \brief Create fresh symbols with the names and sparsity patterns of the inputs

        Convenient for calling the function symbolically or for building a wrapper:
        f(f.mx_in()) is an MX expression of the outputs. Each call creates new symbols;
        they are not the symbols the function was constructed from. All inputs (no
        argument), one by index, or one by name.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("a"), SX.sym("B", Sparsity.diag(2))], [SX(0)],
        ...              ["a", "B"], ["r"])
        >>> f.sx_in()
        [SX(a), SX(
        [[B_0, 00],
         [00, B_1]])]
        >>> f.mx_in(1), f.mx_in("a").shape
        (MX(B), (1, 1))
        >>> f.sx_out(), f.mx_out(0)
        ([SX(r)], MX(r))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> a = SX.sym('a'); B = SX.sym('B', Sparsity.diag(2));
        >> f = Function('f', {a, B}, {SX(0)}, {'a', 'B'}, {'r'});
        >> in = f.sx_in(); disp(in{1}), disp(in{2})
        a
        [[B_0, 00],
         [00, B_1]]
        >> f.mx_in(1), size(f.mx_in('a'))
        B
        1 1
        >> out = f.sx_out(); disp(out{1}), disp(f.mx_out(0))
        r
        r
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> a = ca.SX.sym("a"); B = ca.SX.sym("B", ca.Sparsity.diag(2));
        julia> f = ca.Function("f", [a, B], [ca.SX(0)], ["a", "B"], ["r"]);
        julia> a_in, B_in = ca.sx_in(f);
        julia> a_in
        SX(a)
        julia> B_in
        SX(
        [[B_0, 00],
         [00, B_1]])
        julia> ca.mx_in(f, 1), size(ca.mx_in(f, "a"))
        (B, (1, 1))
        julia> ca.sx_out(f), ca.mx_out(f, 0)
        (SX[r], r)
        \enddoctest
        \endlang

        \sa sx_out, mx_out, sparsity_in

        \identifier{1x4} */
    const SX sx_in(casadi_int iind) const;
    const SX sx_in(const std::string& iname) const {
      return sx_in(index_in(iname));
    }
    const std::vector<SX> sx_in() const;
    const MX mx_in(casadi_int ind) const;
    const MX mx_in(const std::string & iname) const {
      return mx_in(index_in(iname));
    }
    const std::vector<MX> mx_in() const;
#ifndef SWIG
    template<typename T>
    const T sym_in(casadi_int iind) const;
    template<typename T>
    const T sym_in(const std::string& iname) const {
      return sym_in<T>(index_in(iname));
    }
    template<typename T>
    const std::vector<T> sym_in() const;
#endif // SWIG
    ///@}

    ///@{
    /** \brief Create fresh symbols with the names and sparsity patterns of the outputs

        Useful as adjoint seeds or placeholders; see sx_in.

        \identifier{1x5} */
    const SX sx_out(casadi_int oind) const;
    const SX sx_out(const std::string& oname) const {
      return sx_out(index_out(oname));
    }
    const std::vector<SX> sx_out() const;
    const MX mx_out(casadi_int ind) const;
    const MX mx_out(const std::string& oname) const {
      return mx_out(index_out(oname));
    }
    const std::vector<MX> mx_out() const;
    ///@}

    /** \brief Convert between a list of input/output matrices and one flat list of

        nonzeros

        nz_from_in flattens a list of n_in() DM arguments (shape rules as in call) into
        the concatenated nonzeros, in input order and column-major within each input;
        nz_to_in does the reverse, the list must have exactly nnz_in() entries. nz_from_out
        and nz_to_out do the same for the outputs with nnz_out(). This is the memory
        layout that generated C code and FunctionBuffer use.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x"), SX.sym("D", Sparsity.diag(2))], [SX(0)])
        >>> f.nz_from_in([1, DM([[2, 0], [0, 3]])])
        [1.0, 2.0, 3.0]
        >>> f.nz_to_in([1, 2, 3])
        [DM(1), DM(
        [[2, 00],
         [00, 3]])]
        >>> try:
        ...   f.nz_to_in([1, 2])
        ... except RuntimeError as e:
        ...   print("Dimension mismatch. Expecting 3, got 2" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x'), SX.sym('D', Sparsity.diag(2))}, {SX(0)});
        >> f.nz_from_in({1, DM([2, 0; 0, 3])})
        1 2 3
        >> arg = f.nz_to_in([1, 2, 3]); disp(arg{1}), disp(arg{2})
        1
        [[2, 00],
         [00, 3]]
        >> try, f.nz_to_in([1, 2]); catch e, disp(contains(e.message, 'Expecting 3, got 2')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> D = ca.SX.sym("D", ca.Sparsity.diag(2));
        julia> f = ca.Function("f", [ca.SX.sym("x"), D], [ca.SX(0)]);
        julia> ca.nz_from_in(f, [1, ca.DM([2 0; 0 3])])
        3-element Vector{Float64}:
         1.0
         2.0
         3.0
        julia> x_in, D_in = ca.nz_to_in(f, [1, 2, 3]);
        julia> x_in
        DM(1)
        julia> D_in
        DM(
        [[2, 00],
         [00, 3]])
        julia> try
                   ca.nz_to_in(f, [1, 2])
               catch e
                   println(occursin("Dimension mismatch. Expecting 3, got 2", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa convert_in, generate_in, FunctionBuffer

        \identifier{1x6} */
    /// @{
    std::vector<double> nz_from_in(const std::vector<DM>& arg) const;
    std::vector<double> nz_from_out(const std::vector<DM>& arg) const;
    std::vector<DM> nz_to_in(const std::vector<double>& arg) const;
    std::vector<DM> nz_to_out(const std::vector<double>& arg) const;
    ///@}

    /** \brief Convert between a list and a dict of inputs or outputs

        convert_in turns a list of n_in() values into a dict keyed by input name, and a
        dict into a list where absent inputs get default_in. convert_out does the same
        for outputs; absent outputs become nan. An unknown key or a list of the wrong
        length raises an error. Values are passed through unchanged (no shape or sparsity
        check), for any of DM, SX and MX. Useful to prepare arguments for call or
        generate_in.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x"), SX.sym("y")], [SX(0), SX(1)], ["x", "y"], ["a", "b"])
        >>> f.convert_in([1, 2])
        {'x': DM(1), 'y': DM(2)}
        >>> f.convert_in({"y": 2})
        [DM(0), DM(2)]
        >>> f.convert_out({"b": 5})
        [DM(nan), DM(5)]
        >>> try:
        ...   f.convert_in({"z": 1})
        ... except RuntimeError as e:
        ...   print('could not find entry "z"' in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x'), SX.sym('y')}, {SX(0), SX(1)}, {'x', 'y'}, {'a', 'b'});
        >> s = f.convert_in({1, 2}); disp(s.x), disp(s.y)
        1
        2
        >> c = f.convert_in(struct('y', 2)); disp(c{1}), disp(c{2})
        0
        2
        >> c = f.convert_out(struct('b', 5)); disp(c{1}), disp(c{2})
        nan
        5
        >> try, f.convert_in(struct('z', 1)); catch e, msg = e.message; end
        >> contains(msg, 'could not find entry "z"')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x"), ca.SX.sym("y")], [ca.SX(0), ca.SX(1)],
                               ["x", "y"], ["a", "b"]);
        julia> ca.convert_in(f, Dict("y" => 2))
        2-element Vector{DM}:
         DM(0)
         DM(2)
        julia> ca.convert_out(f, Dict("b" => 5))
        2-element Vector{DM}:
         DM(nan)
         DM(5)
        julia> try
                   ca.convert_in(f, Dict("z" => 1))
               catch e
                   println(occursin("could not find entry \"z\"", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa call, default_in, nz_from_in

        \identifier{1x7} */
    /// @{
    DMDict convert_in(const std::vector<DM>& arg) const;
    std::vector<DM> convert_in(const DMDict& arg) const;
    DMDict convert_out(const std::vector<DM>& arg) const;
    std::vector<DM> convert_out(const DMDict& arg) const;
    SXDict convert_in(const std::vector<SX>& arg) const;
    std::vector<SX> convert_in(const SXDict& arg) const;
    SXDict convert_out(const std::vector<SX>& arg) const;
    std::vector<SX> convert_out(const SXDict& arg) const;
    MXDict convert_in(const std::vector<MX>& arg) const;
    std::vector<MX> convert_in(const MXDict& arg) const;
    MXDict convert_out(const std::vector<MX>& arg) const;
    std::vector<MX> convert_out(const MXDict& arg) const;
    /// @}

    /** \brief Check whether the outputs contain symbols that are not inputs

        Only possible for SXFunction and MXFunction constructed with "allow_free".
        Such a function can be called symbolically but not numerically, expanded,
        code-generated or differentiated into a numerically usable form.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); p = SX.sym("p")
        >>> f = Function("f", [x], [x*p], {"allow_free": True})
        >>> f.has_free(), f.get_free(), f.free_sx()
        (True, ['p'], [SX(p)])
        >>> try:
        ...   f(1)
        ... except RuntimeError as e:
        ...   print("since variables [p] are free" in str(e))
        True
        >>> f(SX.sym("q"))
        SX((q*p))
        >>> Function("g", [x], [x]).has_free()
        False
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); p = SX.sym('p');
        >> f = Function('f', {x}, {x*p}, struct('allow_free', true));
        >> disp(f.has_free()), disp(f.get_free()), fr = f.free_sx(); disp(fr{1})
        1
        p
        p
        >> try, f(1); catch e, disp(contains(e.message, 'variables [p] are free')); end
        1
        >> f(SX.sym('q'))
        (q*p)
        >> Function('g', {x}, {x}).has_free()
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); p = ca.SX.sym("p");
        julia> f = ca.Function("f", [x], [x*p], Dict("allow_free" => true));
        julia> ca.has_free(f), ca.get_free(f), ca.free_sx(f)
        (true, ["p"], SX[p])
        julia> try
                   f(1)
               catch e
                   println(occursin("since variables [p] are free", sprint(showerror, e)))
               end
        true
        julia> f(ca.SX.sym("q"))
        SX((q*p))
        julia> ca.has_free(ca.Function("g", [x], [x]))
        false
        \enddoctest
        \endlang

        \sa get_free, free_sx, free_mx

        \identifier{1x8} */
    bool has_free() const;

    /** \brief Get the names of the free variables; empty when there are none

        \sa has_free, free_sx, free_mx

        \identifier{1x9} */
    std::vector<std::string> get_free() const;

    /** \brief Get the free variables of an SXFunction as SX symbols; other classes raise

        \sa has_free, free_mx

        \identifier{1xa} */
    std::vector<SX> free_sx() const;

    /** \brief Get the free variables of an MXFunction as MX symbols; other classes raise

        \lang{python}
        \doctest
        >>> x = MX.sym("x"); p = MX.sym("p")
        >>> f = Function("f", [x], [x*p], {"allow_free": True})
        >>> f.free_mx()
        [MX(p)]
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x'); p = MX.sym('p');
        >> f = Function('f', {x}, {x*p}, struct('allow_free', true));
        >> fr = f.free_mx(); disp(fr{1})
        p
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x"); p = ca.MX.sym("p");
        julia> f = ca.Function("f", [x], [x*p], Dict("allow_free" => true));
        julia> ca.free_mx(f)
        1-element Vector{MX}:
         MX(p)
        \enddoctest
        \endlang

        \sa has_free, free_sx

        \identifier{1xb} */
    std::vector<MX> free_mx() const;

    /** \brief Extract the functions needed for the Lifted Newton method (MXFunction only)

        Returns two functions: one defining the lifted variables and one for their
        initialization. Raises for other classes.

        \identifier{1xc} */
    void generate_lifted(Function& SWIG_OUTPUT(vdef_fcn),
                         Function& SWIG_OUTPUT(vinit_fcn)) const;

    /** \brief Get the number of nodes in the expression graph (SXFunction and MXFunction)

        Counts the unique nodes of the output expressions, inputs included; a cheap
        measure of graph size for comparing transformations. Other classes raise.

        \lang{python}
        \doctest
        >>> x = MX.sym("x")
        >>> Function("f", [x], [sin(x) + sin(x)]).n_nodes()
        4
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x');
        >> Function('f', {x}, {sin(x) + sin(x)}).n_nodes()
        4
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x");
        julia> ca.n_nodes(ca.Function("f", [x], [sin(x) + sin(x)]))
        4
        \enddoctest
        \endlang

        \sa n_instructions

        \identifier{1xd} */
    casadi_int n_nodes() const;

    /** \brief Get the number of instructions of the algorithm (SXFunction and MXFunction)

        The algorithm is the sequence of operations evaluated for a call, including the
        input and output instructions; instruction_id, instruction_input,
        instruction_output, instruction_constant and instruction_MX inspect them by
        index. Other classes raise.

        \lang{python}
        \doctest
        >>> x = SX.sym("x"); y = SX.sym("y")
        >>> f = Function("f", [x, y], [x*y + 1])
        >>> f.n_instructions()
        6
        >>> for k in range(f.n_instructions()):
        ...   print(f.instruction_id(k), f.instruction_input(k), f.instruction_output(k))
        45 [0, 0] [0]
        45 [1, 0] [1]
        3 [0, 1] [0]
        44 [] [1]
        1 [0, 1] [0]
        46 [0] [0, 0]
        >>> f.instruction_id(3) == OP_CONST, f.instruction_constant(3)
        (True, 1.0)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x'); y = SX.sym('y');
        >> f = Function('f', {x, y}, {x*y + 1});
        >> f.n_instructions()
        6
        >> arrayfun(@(k) f.instruction_id(k), 0:5)
        45 45 3 44 1 46
        >> f.instruction_input(2), f.instruction_output(2)
        0 1
        0
        >> f.instruction_id(3) == OP_CONST, f.instruction_constant(3)
        1
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x"); y = ca.SX.sym("y");
        julia> f = ca.Function("f", [x, y], [x*y + 1]);
        julia> ca.n_instructions(f)
        6
        julia> for k in 0:ca.n_instructions(f)-1
                   println(ca.instruction_id(f, k), " ", ca.instruction_input(f, k), " ",
                           ca.instruction_output(f, k))
               end
        45 [0, 0] [0]
        45 [1, 0] [1]
        3 [0, 1] [0]
        44 Int64[] [1]
        1 [0, 1] [0]
        46 [0] [0, 0]
        julia> ca.instruction_id(f, 3) == ca.OP_CONST, ca.instruction_constant(f, 3)
        (true, 1.0)
        \enddoctest
        \endlang

        \sa instruction_id, instructions_sx, n_nodes

        \identifier{1xe} */
    casadi_int n_instructions() const;

    /** \brief Get the operation code of instruction \p k (SXFunction and MXFunction)

        The codes are the OP_* constants (OP_INPUT, OP_OUTPUT, OP_CONST, OP_ADD,
        OP_MUL, OP_SIN, ...; OP_CALL for a call node in an MXFunction).

        \sa n_instructions, instruction_input, instruction_output

        \identifier{1xf} */
    casadi_int instruction_id(casadi_int k) const;

    /** \brief Get the work-vector locations read by instruction \p k

        For an input instruction (OP_INPUT) the pair is [input index, nonzero index].
        (SXFunction and MXFunction)

        \sa instruction_output, n_instructions

        \identifier{1xg} */
    std::vector<casadi_int> instruction_input(casadi_int k) const;

    /** \brief Get the value of a constant instruction (OP_CONST) of an SXFunction

        Only meaningful for constant instructions; MXFunction raises.

        \sa instruction_id, n_instructions

        \identifier{1xh} */
    double instruction_constant(casadi_int k) const;

    /** \brief Get the work-vector locations written by instruction \p k

        For an output instruction (OP_OUTPUT) the pair is [output index, nonzero index].
        (SXFunction and MXFunction)

        \sa instruction_input, n_instructions

        \identifier{1xi} */
    std::vector<casadi_int> instruction_output(casadi_int k) const;

    /** \brief Get the MX expression computed by instruction \p k of an MXFunction

        SXFunction raises; use instructions_sx instead.

        \lang{python}
        \doctest
        >>> x = MX.sym("x")
        >>> f = Function("f", [x], [2*sin(x)])
        >>> [str(f.instruction_MX(k)) for k in range(f.n_instructions())]
        ['input[0][0]', 'sin(x)', '(2.*sin(x))', 'output[0][0]']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x');
        >> f = Function('f', {x}, {2*sin(x)});
        >> for k = 0:f.n_instructions()-1, disp(f.instruction_MX(k)); end
        input[0][0]
        sin(x)
        (2.*sin(x))
        output[0][0]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x");
        julia> f = ca.Function("f", [x], [2sin(x)]);
        julia> [string(ca.instruction_MX(f, k)) for k in 0:ca.n_instructions(f)-1]
        4-element Vector{String}:
         "input[0][0]"
         "sin(x)"
         "(2.*sin(x))"
         "output[0][0]"
        \enddoctest
        \endlang

        \sa n_instructions, instructions_sx

        \identifier{1xj} */
    MX instruction_MX(casadi_int k) const;

    /** \brief Get the SX expression of every instruction of an SXFunction, as one vector

        Entry k is the expression computed by instruction k; input and output
        instructions have no expression and are nan. MXFunction raises.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [2*sin(x)])
        >>> f.instructions_sx()
        SX(@1=nan, @2=sin(x), [@1, @2, (2.*@2), @1])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {2*sin(x)});
        >> f.instructions_sx()
        @1=nan, @2=sin(x), [@1, @2, (2.*@2), @1]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [2sin(x)]);
        julia> ca.instructions_sx(f)
        SX(@1=nan, @2=sin(x), [@1, @2, (2.*@2), @1])
        \enddoctest
        \endlang

        \sa n_instructions, instruction_MX

        \identifier{1xk} */
    SX instructions_sx() const;

    ///@{
    /** \brief Check whether the class supports forward (has_spfwd) or reverse (has_sprev)

        sparsity propagation

        Both are true for SXFunction and MXFunction; classes without it fall back to
        dense Jacobian sparsity patterns.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.has_spfwd(), f.has_sprev()
        (True, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> fprintf('%d %d\n', f.has_spfwd(), f.has_sprev())
        1 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.has_spfwd(f), ca.has_sprev(f)
        (true, true)
        \enddoctest
        \endlang

        \sa jac_sparsity

        \identifier{1xl} */
    bool has_spfwd() const;
    bool has_sprev() const;
    ///@}

    /** \brief Get the length of the input pointer array needed to evaluate (at least n_in)

        sz_arg, sz_res, sz_iw and sz_w are the sizes of the arg, res, iw and w arrays of
        the low-level evaluation interface used by generated code and FunctionBuffer.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x"), SX.sym("y")], [SX(0)])
        >>> f.sz_arg(), f.sz_res(), f.sz_iw(), f.sz_w()
        (2, 1, 0, 1)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x'), SX.sym('y')}, {SX(0)});
        >> fprintf('%d %d %d %d\n', f.sz_arg(), f.sz_res(), f.sz_iw(), f.sz_w())
        2 1 0 1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x"), ca.SX.sym("y")], [ca.SX(0)]);
        julia> ca.sz_arg(f), ca.sz_res(f), ca.sz_iw(f), ca.sz_w(f)
        (2, 1, 0, 1)
        \enddoctest
        \endlang

        \sa sz_res, sz_iw, sz_w, FunctionBuffer

        \identifier{1xm} */
    size_t sz_arg() const;

    /** \brief Get the length of the output pointer array needed to evaluate (at least

        n_out)

        \identifier{1xn} */
    size_t sz_res() const;

    /** \brief Get the length of the integer work vector needed to evaluate

        \identifier{1xo} */
    size_t sz_iw() const;

    /** \brief Get the length of the real work vector needed to evaluate

        \identifier{1xp} */
    size_t sz_w() const;

#ifndef SWIG
    /** \brief Get number of temporary variables needed

        \identifier{1xq} */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Set the (persistent) work vectors

        \identifier{1xr} */
    void set_work(const double**& arg, double**& res,
      casadi_int*& iw, double*& w, int mem=0) const;

    /** \brief Set the (temporary) work vectors

        \identifier{1xs} */
    void set_temp(const double** arg, double** res,
        casadi_int* iw, double* w, int mem=0) const;

    /** \brief Set the (persistent and temporary) work vectors

        \identifier{1xt} */
    void setup(const double** arg, double** res, casadi_int* iw, double* w, int mem=0) const;

    /** \brief Call using a map

        \identifier{1xu} */
    template<typename M>
    void call_gen(const std::map<std::string, M>& arg, std::map<std::string, M>& res,
               bool always_inline, bool never_inline) const;

    /** \brief List merge opportunitities

        \identifier{2b6} */
    void merge(const std::vector<MX>& arg,
        std::vector<MX>& subs_from, std::vector<MX>& subs_to) const;
#endif // SWIG
    /// \endcond

    /** \brief Get the name given at construction ("null" for a null object)

        \lang{python}
        \doctest
        >>> Function("f", [SX.sym("x")], [SX(0)]).name()
        'f'
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Function('f', {SX.sym('x')}, {SX(0)}).name()
        f
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.name(ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]))
        "f"
        \enddoctest
        \endlang

        \sa check_name

        \identifier{1xv} */
    const std::string& name() const;

    /** \brief Check whether the function is of a given class

        \p type is a class name as printed after the signature: "SXFunction",
        "MXFunction", "Map", "Switch", "External", "Nlpsol", "Conic", "Integrator",
        "Rootfinder" and the plugin names (e.g. "Ipopt"). With \p recursive true (default)
        base classes match as well ("Nlpsol" for an Ipopt solver, "Map" for a ThreadMap,
        "FunctionInternal" for everything).

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> f.is_a("SXFunction"), f.is_a("MXFunction"), f.class_name()
        (True, False, 'SXFunction')
        >>> f.map(2, "thread").is_a("Map"), f.map(2, "thread").is_a("Map", False)
        (True, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> fprintf('%d %d %s\n', f.is_a('SXFunction'), f.is_a('MXFunction'), f.class_name())
        1 0 SXFunction
        >> F = f.map(2, 'thread'); fprintf('%d %d\n', F.is_a('Map'), F.is_a('Map', false))
        1 0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> ca.is_a(f, "SXFunction"), ca.is_a(f, "MXFunction"), ca.class_name(f)
        (true, false, "SXFunction")
        julia> F = f.map(2, "thread");
        julia> ca.is_a(F, "Map"), ca.is_a(F, "Map", false)
        (true, false)
        \enddoctest
        \endlang

        \sa class_name

        \identifier{1xw} */
    bool is_a(const std::string& type, bool recursive=true) const;

    /** \brief Check whether a string is a valid function name

        A valid name starts with a letter and continues with letters, digits and
        underscores, with no two consecutive underscores; it may not be one of the
        keywords "null", "jac", "hess". Constructors and generate raise for invalid
        names; names must be valid C identifiers because they appear in generated code.

        \lang{python}
        \doctest
        >>> Function.check_name("f_1"), Function.check_name("1f")
        (True, False)
        >>> Function.check_name("a__b"), Function.check_name("jac"), Function.check_name("")
        (False, False, False)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fprintf('%d %d\n', Function.check_name('f_1'), Function.check_name('1f'))
        1 0
        >> fprintf('%d %d\n', Function.check_name('a__b'), Function.check_name('jac'))
        0 0
        >> Function.check_name('')
        0
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> F = ca.Function;
        julia> F.check_name("f_1"), F.check_name("1f")
        (true, false)
        julia> F.check_name("a__b"), F.check_name("jac"), F.check_name("")
        (false, false, false)
        \enddoctest
        \endlang

        \sa fix_name

        \identifier{1xx} */
    static bool check_name(const std::string& name);

    /** \brief Turn any string into a valid function name as defined by check_name

        Runs of non-alphanumeric characters become a single underscore, a leading
        non-letter gets an "a" prepended, the empty string becomes "unnamed" and a
        keyword gets "1" appended. A valid name is returned unchanged.

        \lang{python}
        \doctest
        >>> Function.fix_name("my func!"), Function.fix_name("1f"), Function.fix_name("")
        ('my_func_', 'a1f', 'unnamed')
        >>> Function.fix_name("jac"), Function.fix_name("ok")
        ('jac1', 'ok')
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> Function.fix_name('my func!'), Function.fix_name('1f')
        my_func_
        a1f
        >> Function.fix_name(''), Function.fix_name('jac')
        unnamed
        jac1
        >> Function.fix_name('ok')
        ok
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> F = ca.Function;
        julia> F.fix_name("my func!"), F.fix_name("1f"), F.fix_name("")
        ("my_func_", "a1f", "unnamed")
        julia> F.fix_name("jac"), F.fix_name("ok")
        ("jac1", "ok")
        \enddoctest
        \endlang

        \sa check_name

        \identifier{1xy} */
    static std::string fix_name(const std::string& name);

    /** \brief Build a function from a serialization stream

        \identifier{1xz} */
    static Function deserialize(std::istream& stream);

    /** \brief Build a function from a string produced by serialize

        Raises "DeserializingStream sanity check failed" for data that is not a
        serialized CasADi object.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function.deserialize(Function("f", [x], [2*x]).serialize())
        >>> f(3)
        DM(6)
        >>> try:
        ...   Function.deserialize("not a function")
        ... except RuntimeError as e:
        ...   print("sanity check failed" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function.deserialize(Function('f', {x}, {2*x}).serialize());
        >> f(3)
        6
        >> try, Function.deserialize('not a function'); catch e, msg = e.message; end
        >> contains(msg, 'sanity check failed')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function.deserialize(ca.serialize(ca.Function("f", [x], [2x])));
        julia> f(3)
        DM(6)
        julia> try
                   ca.Function.deserialize("not a function")
               catch e
                   println(occursin("sanity check failed", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa serialize, load

        \identifier{1y0} */
    static Function deserialize(const std::string& s);

    /** \brief Read a function from a file written by save

        The file must contain exactly one Function (as written by save); a file
        written by FileSerializer with other objects raises an error pointing to
        FileDeserializer, a missing file raises "Could not open".

        \sa save, deserialize, FileDeserializer

        \identifier{1y1} */
    static Function load(const std::string& filename);

    /** \brief Build a function from a DeserializingStream (C++ serialization API)

        \identifier{1y2} */
    static Function deserialize(DeserializingStream& s);

    /** \brief Raise an error unless input \p i is \p nrow-by-\p ncol

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x", 2)], [SX(1, 3)])
        >>> f.assert_size_in(0, 2, 1)
        >>> try:
        ...   f.assert_size_out(0, 2, 1)
        ... except RuntimeError as e:
        ...   print("Expected 2-by-1 but got 1-by-3" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x', 2)}, {SX(1, 3)});
        >> f.assert_size_in(0, 2, 1)
        >> try, f.assert_size_out(0, 2, 1); catch e, msg = e.message; end
        >> contains(msg, 'Expected 2-by-1 but got 1-by-3')
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x", 2)], [ca.SX(1, 3)]);
        julia> ca.assert_size_in(f, 0, 2, 1)
        julia> try
                   ca.assert_size_out(f, 0, 2, 1)
               catch e
                   println(occursin("Expected 2-by-1 but got 1-by-3", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa assert_size_out, size_in

        \identifier{2mk} */
    void assert_size_in(casadi_int i, casadi_int nrow, casadi_int ncol) const;

    /// Raise an error unless output \p i is \p nrow-by-\p ncol
    void assert_size_out(casadi_int i, casadi_int nrow, casadi_int ncol) const;

    /** \brief Raise an error unless output \p i has the sparsity \p sp repeated \p n

        times horizontally

        With \p allow_all_zero_sparse, an output without nonzeros of the right shape
        passes. Used to check derivative functions.

        \identifier{2ml} */
    void assert_sparsity_out(casadi_int i, const Sparsity& sp,
        casadi_int n = 1, bool allow_all_zero_sparse = true) const;

    /** \brief Reserve a memory object for evaluation and return its index

        Each memory object holds the state of one evaluation context (solver workspace,
        statistics); ordinary calls use index 0. Reserving separate objects allows
        concurrent evaluation from several threads with the low-level interface and
        separate stats. Release it with release when done.

        \lang{python}
        \doctest
        >>> f = Function("f", [SX.sym("x")], [SX(0)])
        >>> m = f.checkout()
        >>> m
        1
        >>> f.release(m)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> f = Function('f', {SX.sym('x')}, {SX(0)});
        >> m = f.checkout(); disp(m)
        1
        >> f.release(m)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> f = ca.Function("f", [ca.SX.sym("x")], [ca.SX(0)]);
        julia> m = ca.checkout(f)
        1
        julia> ca.release(f, m)
        \enddoctest
        \endlang

        \sa release, stats, FunctionBuffer

        \identifier{2mm} */
    casadi_int checkout() const;

    /// Give back a memory object obtained with checkout
    void release(int mem) const;

#ifndef SWIG
    /// Get memory object
    void* memory(int ind) const;

    static std::vector<SX> order(const std::vector<SX>& expr);
    static std::vector<MX> order(const std::vector<MX>& expr);
#endif // SWIG

    /** \brief Get the derivative and map functions currently cached for this function

        jacobian, forward, reverse, wrap and serial map store their result in a cache
        keyed by name (jac_f, fwd1_f, adj1_f, wrap_f, map3_f, ...). The cache holds weak
        references: an entry stays only while some other reference keeps the function
        alive. Option "cache" prepopulates it, e.g. with a hand-written Jacobian.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> f = Function("f", [x], [sin(x)])
        >>> f.cache()
        {}
        >>> J = f.jacobian()
        >>> list(f.cache().keys())
        ['jac_f']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> f = Function('f', {x}, {sin(x)});
        >> isempty(fieldnames(f.cache()))
        1
        >> J = f.jacobian(); disp(strjoin(fieldnames(f.cache()), ' '))
        jac_f
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> f = ca.Function("f", [x], [sin(x)]);
        julia> ca.cache(f)
        Dict{String, Any}()
        julia> J = ca.jacobian(f);
        julia> collect(keys(ca.cache(f)))
        1-element Vector{String}:
         "jac_f"
        \enddoctest
        \endlang

        \sa jacobian, forward, reverse

        \identifier{26i} */
    Dict cache() const;

    /** \brief Get the names of the functions this function directly owns

        For an SXFunction or MXFunction: the functions called in its graph (an MXFunction
        called with MX arguments is inlined and does not appear; SX functions, solvers
        and calls made with never_inline do); for a map: the mapped function; for a
        solver: its helper functions (e.g. nlp_f, nlp_jac_g for nlpsol, which can also be
        code-generated with generate_dependencies). Empty for leaf classes. Note that
        has_function only answers true for solver objects; for SX and MX functions it is
        false even for names in this list.

        \lang{python}
        \doctest
        >>> x = MX.sym("x")
        >>> f = Function("f", [x], [sin(x)]).expand()
        >>> g = Function("g", [x], [f(x) + 1])
        >>> g.get_function(), g.has_function("f")
        (['f'], True)
        >>> g.get_function("f")
        Function(f:(i0)->(o0) SXFunction)
        >>> try:
        ...   g.get_function("h")
        ... except RuntimeError as e:
        ...   print("No such function 'h'" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x');
        >> f = Function('f', {x}, {sin(x)}).expand();
        >> g = Function('g', {x}, {f(x) + 1});
        >> g.get_function(), g.has_function('f')
        f
        1
        >> g.get_function('f')
        f:(i0)->(o0) SXFunction
        >> try, g.get_function('h'); catch e, disp(contains(e.message, 'No such function')); end
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x");
        julia> f = ca.expand(ca.Function("f", [x], [sin(x)]));
        julia> g = ca.Function("g", [x], [f(x) + 1]);
        julia> ca.get_function(g), ca.has_function(g, "f")
        (["f"], true)
        julia> ca.get_function(g, "f")
        Function(f:(i0)->(o0) SXFunction)
        julia> try
                   ca.get_function(g, "h")
               catch e
                   println(occursin("No such function 'h'", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \sa has_function, find_functions, oracle

        \identifier{1y3} */
    std::vector<std::string> get_function() const;

    /** \brief Get a directly owned function by name; unknown names raise an error

        \sa has_function, find_function

        \identifier{1y4} */
    Function get_function(const std::string &name) const;

    /** \brief Check if a particular dependency exists

        True if \p fname is among the names returned by get_function().

        \sa get_function

        \identifier{1y5} */
    bool has_function(const std::string& fname) const;

    /** \brief Get all functions embedded in this function, recursively

        Walks the owned functions (get_function) and their owned functions in turn,
        returning each distinct Function once, in depth-first discovery order.

        \param max_depth Maximum nesting depth to descend; 0 returns only the directly
          owned functions, a negative number means no limit
        \return Unique list of embedded functions

        \lang{python}
        \doctest
        >>> x = MX.sym("x")
        >>> inner = Function("inner", [x], [sin(x)]).expand()
        >>> middle = Function("middle", [x], [inner(x) + 1])
        >>> outer = Function("outer", [x], [middle.call([x], False, True)[0]*2])
        >>> [g.name() for g in outer.find_functions()]
        ['middle', 'inner']
        >>> [g.name() for g in outer.find_functions(0)]
        ['middle']
        >>> outer.find_function("inner")
        Function(inner:(i0)->(o0) SXFunction)
        >>> Function("o2", [x], [middle(x)*2]).get_function()
        ['inner']
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = MX.sym('x');
        >> inner = Function('inner', {x}, {sin(x)}).expand();
        >> middle = Function('middle', {x}, {inner(x) + 1});
        >> res = middle.call({x}, false, true); outer = Function('outer', {x}, {res{1}*2});
        >> fs = outer.find_functions(); cellfun(@(g) disp(g.name()), fs)
        middle
        inner
        >> fs = outer.find_functions(0); disp(fs{1}.name())
        middle
        >> outer.find_function('inner')
        inner:(i0)->(o0) SXFunction
        >> Function('o2', {x}, {middle(x)*2}).get_function()
        inner
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.MX.sym("x");
        julia> inner = ca.expand(ca.Function("inner", [x], [sin(x)]));
        julia> middle = ca.Function("middle", [x], [inner(x) + 1]);
        julia> outer = ca.Function("outer", [x], [ca.call(middle, [x], false, true)[1] * 2]);
        julia> [ca.name(g) for g in ca.find_functions(outer)]
        2-element Vector{String}:
         "middle"
         "inner"
        julia> [ca.name(g) for g in ca.find_functions(outer, 0)]
        1-element Vector{String}:
         "middle"
        julia> ca.find_function(outer, "inner")
        Function(inner:(i0)->(o0) SXFunction)
        julia> ca.get_function(ca.Function("o2", [x], [middle(x) * 2]))
        1-element Vector{String}:
         "inner"
        \enddoctest
        \endlang

        \sa find_function, get_function

        \identifier{1y6} */
    std::vector<Function> find_functions(casadi_int max_depth = -1) const;

    /** \brief Get an embedded function by name, searching recursively

        \param name Name of the function; raises "'name' not found" when absent
        \param max_depth Maximum nesting depth, negative for no limit

        \sa find_functions

        \identifier{1y7} */
    Function find_function(const std::string &name, casadi_int max_depth=-1) const;

    /** \brief Get class-specific information as a dict

        Empty for most classes; e.g. a Switch (conditional, if_else) reports its branch
        functions.

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> Function("f", [x], [x]).info()
        {}
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> isempty(fieldnames(Function('f', {x}, {x}).info()))
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.info(ca.Function("f", [x], [x]))
        Dict{String, Any}()
        \enddoctest
        \endlang

        \identifier{2mn} */
    Dict info() const;

#ifndef SWIG
    protected:
    ///@{
    /** \brief Called by constructors

        \identifier{1y8} */
    void construct(const std::string& name,
                   const std::vector<SX>& ex_in, const std::vector<SX>& ex_out,
                   const std::vector<std::string>& name_in,
                   const std::vector<std::string>& name_out,
                   const Dict& opts);
    void construct(const std::string& name,
                   const std::vector<MX>& ex_in, const std::vector<MX>& ex_out,
                   const std::vector<std::string>& name_in,
                   const std::vector<std::string>& name_out,
                   const Dict& opts);
    template<typename M>
      void construct(const std::string& name, const std::map<std::string, M>& dict,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const Dict& opts);
    ///@}

    /// Helper function for parsing .casadi files
    static bool proceed_to(std::istream& file, const std::string& str);

    /// Helper function for mapaccum
    Function mapaccum(const std::string& name, const std::vector<Function>& chain,
                      casadi_int n_accum=1, const Dict& opts = Dict()) const;

#ifdef WITH_EXTRA_CHECKS
    public:
    // How many times have we passed through
    // operator()(const double** arg, double** res, casadi_int* iw, double* w, int mem)?
    static thread_local casadi_int call_depth_;
#endif


#endif // SWIG



  };


/** \brief Evaluate a Function repeatedly with minimal overhead

    A FunctionBuffer holds the work vectors and a memory object of a Function and
    evaluates it directly on caller-owned contiguous double buffers, skipping argument
    conversion and allocation. Bind each input with set_arg and each output with set_res
    (nonzeros in column-major order, see nz_from_in), then evaluate; the results appear
    in the bound output buffers and ret() gives the return code (0 on success). The
    buffers must stay alive and unchanged in size while the FunctionBuffer is used.

    \lang{python}
    Function.buffer() returns a pair (FunctionBuffer, callable): the callable runs one
    evaluation; bind numpy arrays with memoryview.

    \doctest
    >>> x = SX.sym("x"); y = SX.sym("y", 2)
    >>> f = Function("f", [x, y], [x*y], ["x", "y"], ["r"])
    >>> fb, evaluate = f.buffer()
    >>> a = np.array([2.0]); b = np.array([3.0, 4.0]); r = np.zeros(2)
    >>> fb.set_arg(0, memoryview(a)); fb.set_arg(1, memoryview(b))
    >>> fb.set_res(0, memoryview(r))
    >>> evaluate()
    >>> r, fb.ret()
    (array([6., 8.]), 0)
    >>> a[0] = 5; evaluate(); print(r)
    [15. 20.]
    >>> try:
    ...   fb.set_arg(1, memoryview(a))
    ... except RuntimeError as e:
    ...   print("Buffer is not large enough" in str(e))
    True
    \enddoctest
    \endlang

    \sa Function::nz_from_in, Function::checkout, Function::sz_w

    \identifier{1y9} */
class CASADI_EXPORT FunctionBuffer {
  Function f_;
  std::vector<double> w_;
  std::vector<casadi_int> iw_;
  std::vector<const double*> arg_;
  std::vector<double*> res_;
  FunctionInternal* f_node_;
  casadi_int mem_;
  void *mem_internal_;
  int ret_;
public:
  /** \brief Create a buffer for \p f, allocating its work vectors and a memory object

      \identifier{1ya} */
  FunctionBuffer(const Function& f);
#ifndef SWIG
  ~FunctionBuffer();
  FunctionBuffer(const FunctionBuffer& f);
  FunctionBuffer& operator=(const FunctionBuffer& f);
#endif // SWIG

  /** \brief Bind a contiguous buffer holding the nonzeros of input \p i

      \p size is the buffer size in bytes; it must be at least nnz_in(i)*sizeof(double)
      or an error is raised. The values are read in column-major (Fortran) order of the
      input sparsity.

      \lang{python}
      Pass a memoryview of a numpy array: set_arg(0, memoryview(a)).
      \endlang

      \identifier{1yb} */
  void set_arg(casadi_int i, const double* a, casadi_int size);

  /** \brief Bind a contiguous buffer that receives the nonzeros of output \p i

      \p size is the buffer size in bytes, at least nnz_out(i)*sizeof(double). Values are
      written in column-major order of the output sparsity.

      \lang{python}
      Pass a memoryview of a numpy array: set_res(0, memoryview(r)).
      \endlang

      \identifier{1yc} */
  void set_res(casadi_int i, double* a, casadi_int size);
  /// Get the return code of the last evaluation (0 on success)
  int ret();
  void _eval();
  void* _self() { return this; }
  /// Get the statistics of the memory object used by this buffer, see Function::stats
  Dict stats() const;
};

void CASADI_EXPORT _function_buffer_eval(void* raw);

template<>
const SX CASADI_EXPORT Function::sym_in(casadi_int iind) const;
template<>
const MX CASADI_EXPORT Function::sym_in(casadi_int iind) const;
template<>
const std::vector<SX> CASADI_EXPORT Function::sym_in() const;
template<>
const std::vector<MX> CASADI_EXPORT Function::sym_in() const;

} // namespace casadi

#include "casadi_interrupt.hpp"
#include "runtime/shared.hpp"

#endif // CASADI_FUNCTION_HPP
