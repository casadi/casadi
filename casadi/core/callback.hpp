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

#ifndef CASADI_CALLBACK_HPP
#define CASADI_CALLBACK_HPP

#include "function.hpp"

namespace casadi {
  /** Forward declaration of internal class */
  class CallbackInternal;

  /** \brief User-defined Function implemented by subclassing

   A Callback is a Function whose numerical evaluation (and optionally its derivatives) is
   written by the user in Python, MATLAB or C++ by overriding virtual methods. It can be
   called like any Function, embedded in MX graphs, mapped, and used inside NLPs and
   integrators.

   A subclass must call the base constructor and then construct(name, opts) from its own
   constructor; construct() queries the virtual methods below (get_n_in, get_n_out,
   get_sparsity_in/out, get_name_in/out) and calls init() and finalize(), so all data
   those methods need must be set before construct() runs. Only eval() is mandatory. Defaults:
   one scalar input "i0" and one scalar output "o0".

   Object lifetime: the Function side holds a raw pointer to the user object and does not
   own it. The user must keep the Callback object alive (in Python: keep a reference) for
   as long as any Function, expression or solver that uses it exists; otherwise evaluation
   fails with "Callback object has been deleted". The same holds for objects returned from
   get_jacobian/get_forward/get_reverse. Callback objects cannot be copied.

   Derivatives: by default a Callback is not differentiable, so jacobian(), forward/reverse
   sweeps and NLP solvers that need gradients fail with "Derivatives cannot be calculated".
   Either override has_jacobian/get_jacobian, has_forward/get_forward or
   has_reverse/get_reverse, or pass the option "enable_fd": True to construct() to use
   finite differences ("fd_method": "central" (default), "forward", "backward" or
   "smoothing"; "fd_options" such as {"h": 1e-6} tune the step). A user-provided Jacobian
   is itself differentiated with finite differences only when it is a Callback with
   "enable_fd"; a Function built from SX/MX gives exact second derivatives.

   \lang{python}
   \doctest
   >>> class MyCallback(Callback):
   ...   def __init__(self, name, d, opts={}):
   ...     Callback.__init__(self)
   ...     self.d = d
   ...     self.construct(name, opts)
   ...   def get_n_in(self): return 1
   ...   def get_n_out(self): return 1
   ...   def eval(self, arg):
   ...     return [sin(self.d*arg[0])]
   >>> f = MyCallback("f", 0.5)
   >>> f
   Callback(f:(i0)->(o0) CallbackInternal)
   >>> f(2)
   DM(0.8415)
   >>> x = MX.sym("x")
   >>> f(x)
   MX(f(x){0})
   >>> try:
   ...   jacobian(f(x), x)
   ... except Exception as e:
   ...   print("Derivatives cannot be calculated for f" in str(e))
   True
   >>> f = MyCallback("f", 0.5, {"enable_fd": True})
   >>> J = Function("J", [x], [jacobian(f(x), x)])
   >>> J(2)
   DM(0.2702)
   \enddoctest
   \endlang

   \lang{matlab}
   \doctest
   >> c = {'classdef MyCallback < casadi.Callback', 'properties', 'd', 'end', 'methods'};
   >> c = [c, {'function self = MyCallback(name, d, varargin)'}];
   >> c = [c, {'self@casadi.Callback(); self.d = d; construct(self, name, varargin{:});'}];
   >> c = [c, {'end'}];
   >> c = [c, {'function r = eval(self, arg)', 'r = {sin(self.d*arg{1})};'}];
   >> c = [c, {'end', 'end', 'end'}];
   >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'MyCallback.m'), 'w');
   >> fprintf(fid, '%s\n', c{:}); fclose(fid);
   >> f = MyCallback('f', 0.5); disp(f)
   f:(i0)->(o0) CallbackInternal
   >> f(2)
   0.8415
   >> x = MX.sym('x'); disp(f(x))
   f(x){0}
   >> try, jacobian(f(x), x); catch e, disp(contains(e.message, 'Derivatives cannot be')); end
   1
   >> f = MyCallback('f', 0.5, struct('enable_fd', true));
   >> J = Function('J', {x}, {jacobian(f(x), x)}); J(2)
   0.2702
   >> rmpath(d); rmdir(d, 's');
   \enddoctest
   \endlang

   \lang{julia}
   \doctest
   julia> struct MyCallback <: ca.Callback
              d::Float64
          end
   julia> ca.get_n_in(::MyCallback) = 1
   julia> ca.get_n_out(::MyCallback) = 1
   julia> ca.eval_(cb::MyCallback, arg) = [sin(cb.d * arg[1])]
   julia> f = ca.callback(MyCallback(0.5), "f");
   julia> ca.Function(f)
   Function(f:(i0)->(o0) CallbackInternal)
   julia> f(2)
   DM(0.8415)
   julia> x = ca.MX.sym("x");
   julia> f(x)
   MX(f(x){0})
   julia> try
              ca.jacobian(f(x), x)
          catch e
              println(occursin("Derivatives cannot be calculated for f", sprint(showerror, e)))
          end
   true
   julia> f = ca.callback(MyCallback(0.5), "f", Dict("enable_fd" => true));
   julia> J = ca.Function("J", [x], [ca.jacobian(f(x), x)]);
   julia> J(2)
   DM(0.2702)
   \enddoctest
   \endlang

   The callback object must outlive every use of the Function:

   \lang{python}
   \doctest
   >>> class Sq(Callback):
   ...   def __init__(self, name, opts={}):
   ...     Callback.__init__(self)
   ...     self.construct(name, opts)
   ...   def eval(self, arg):
   ...     return [arg[0]**2]
   >>> x = MX.sym("x")
   >>> f = Sq("f")
   >>> F = Function("F", [x], [f(x)])
   >>> F(3)
   DM(9)
   >>> del f
   >>> try:
   ...   F(3)
   ... except Exception as e:
   ...   print(str(e).splitlines()[-1])
   Callback object has been deleted
   \enddoctest
   \endlang

   \lang{matlab}
   \doctest
   >> c = {'classdef SqCb < casadi.Callback', 'methods'};
   >> c = [c, {'function self = SqCb(name)', 'self@casadi.Callback(); construct(self, name);'}];
   >> c = [c, {'end'}];
   >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end', 'end', 'end'}];
   >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqCb.m'), 'w');
   >> fprintf(fid, '%s\n', c{:}); fclose(fid);
   >> x = MX.sym('x'); f = SqCb('f'); F = Function('F', {x}, {f(x)}); disp(F(3))
   9
   >> delete(f)
   >> try, F(3); catch e, disp(contains(e.message, 'Callback object has been deleted')); end
   ...
   1
   >> rmpath(d); rmdir(d, 's');
   \enddoctest
   \endlang

   \lang{julia}
   \doctest
   julia> struct Sq <: ca.Callback end
   julia> ca.eval_(::Sq, arg) = [arg[1]^2]
   julia> x = ca.MX.sym("x");
   julia> f = ca.callback(Sq(), "f");
   julia> F = ca.Function("F", [x], [f(x)]);
   julia> F(3)
   DM(9)
   julia> f = nothing; GC.gc()
   julia> try
              F(3)
          catch e
              println(last(split(sprint(showerror, e), "\n")))
          end
   Callback object has been deleted
   \enddoctest
   \endlang

   Use in an NLP with finite differences:

   \lang{python}
   \doctest
   >>> class Sq(Callback):
   ...   def __init__(self, name, opts={}):
   ...     Callback.__init__(self)
   ...     self.construct(name, opts)
   ...   def eval(self, arg):
   ...     return [arg[0]**2]
   >>> f = Sq("f", {"enable_fd": True})
   >>> opti = Opti()
   >>> v = opti.variable()
   >>> opti.minimize(f(v-1))
   >>> opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
   >>> opti.solve().value(v)
   1.0
   \enddoctest
   \endlang

   \lang{matlab}
   \doctest
   >> c = {'classdef SqFd < casadi.Callback', 'methods'};
   >> c = [c, {'function self = SqFd(name, opts)'}];
   >> c = [c, {'self@casadi.Callback(); construct(self, name, opts);', 'end'}];
   >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end', 'end', 'end'}];
   >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqFd.m'), 'w');
   >> fprintf(fid, '%s\n', c{:}); fclose(fid);
   >> f = SqFd('f', struct('enable_fd', true));
   >> opti = Opti(); v = opti.variable(); opti.minimize(f(v-1));
   >> opti.solver('ipopt', struct('print_time', false), struct('print_level', 0, 'sb', 'yes'));
   >> sol = opti.solve(); disp(sol.value(v))
   1.0000
   >> rmpath(d); rmdir(d, 's');
   \enddoctest
   \endlang

   \lang{julia}
   \doctest
   julia> struct Sq <: ca.Callback end
   julia> ca.eval_(::Sq, arg) = [arg[1]^2]
   julia> f = ca.callback(Sq(), "f", Dict("enable_fd" => true));
   julia> opti = ca.Opti();
   julia> v = ca.variable(opti);
   julia> ca.minimize(opti, f(v - 1))
   julia> ca.solver(opti, "ipopt", Dict("print_time" => false),
                     Dict("print_level" => 0, "sb" => "yes"))
   julia> ca.value(ca.solve(opti), v)
   1.0
   \enddoctest
   \endlang

   \sa Function, external

   \author Joris Gillis, Joel Andersson
   \date 2015

      \identifier{o0} */
  class CASADI_EXPORT Callback : public Function {
  public:
    /** \brief Get type name

        \identifier{o1} */
    static std::string type_name() {return "Callback";}

    /** \brief Default constructor

        Creates a null Function; construct() must be called before use.

        \identifier{o2} */
    Callback();

    /** \brief Copy constructor (throws an error)

        Callback objects cannot be copied; share the Function instead.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        >>> f = Sq("f")
        >>> try:
        ...   Callback(f)
        ... except Exception as e:
        ...   print("Callback objects cannot be copied" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqCopy < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqCopy(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqCopy.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqCopy('f');
        >> try, Callback(f); catch e, disp(contains(e.message, 'cannot be copied')); end
        1
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> f = ca.callback(Sq(), "f");
        julia> try
                   ca.CallbackFunction(f)
               catch e
                   println(e isa MethodError)
               end
        true
        \enddoctest
        \endlang

        \identifier{o3} */
    Callback(const Callback& obj);

    /** \brief  Destructor

        Detaches the user object; a Function that still references it fails afterwards.

        \identifier{o4} */
    virtual ~Callback();

    /** \brief Construct the internal Function object

        The class constructor only creates a null pointer; this call creates the Function.
        It must be called once from the subclass constructor, after all data needed by
        get_n_in/get_n_out/get_sparsity_in/get_sparsity_out/get_name_in/get_name_out is set,
        because those methods are evaluated here. init() and finalize() are called at the
        end. Calling it twice raises an error.

        \param name Name of the Function
        \param opts Function options, e.g. "enable_fd", "fd_method", "fd_options",
               "enable_forward", "enable_reverse", "enable_jacobian"

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        >>> f = Sq("f", {"enable_fd": True, "fd_method": "forward"})
        >>> f.name(), f.n_in(), f.n_out()
        ('f', 1, 1)
        >>> try:
        ...   f.construct("g")
        ... except Exception as e:
        ...   print("Internal class already created" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqOpts < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqOpts(name, opts)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name, opts);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqOpts.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqOpts('f', struct('enable_fd', true, 'fd_method', 'forward'));
        >> fprintf('%s %d %d\n', f.name(), f.n_in(), f.n_out())
        f 1 1
        >> try, f.construct('g'); catch e, disp(contains(e.message, 'already created')); end
        1
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> f = ca.callback(Sq(), "f", Dict("enable_fd" => true, "fd_method" => "forward"));
        julia> ca.name(f), ca.n_in(f), ca.n_out(f)
        ("f", 1, 1)
        julia> try
                   ca.construct(f, "g")
               catch e
                   println(occursin("Internal class already created", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{o5} */
    void construct(const std::string& name, const Dict& opts=Dict());

    /** \brief Initialize the object

        Called from construct() after the inputs and outputs are known and the options are
        processed, before finalize(). Override to allocate data; the default does nothing.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     print("before construct")
        ...     self.construct(name, opts)
        ...     print("after construct")
        ...   def init(self): print("init", self.n_in())
        ...   def finalize(self): print("finalize")
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        >>> f = Sq("f")
        before construct
        init 1
        finalize
        after construct
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqInit < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqInit(name)'}];
        >> c = [c, {'self@casadi.Callback(); disp(''before construct'');'}];
        >> c = [c, {'construct(self, name); disp(''after construct'');', 'end'}];
        >> c = [c, {'function init(self)', 'fprintf(''init %d\n'', self.n_in());', 'end'}];
        >> c = [c, {'function finalize(self)', 'disp(''finalize'');', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqInit.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqInit('f');
        before construct
        init 1
        finalize
        after construct
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.init(::Sq) = println("init")
        julia> ca.finalize(::Sq) = println("finalize")
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> begin
                   println("before construct")
                   f = ca.callback(Sq(), "f")
                   println("after construct")
               end
        before construct
        init
        finalize
        after construct
        \enddoctest
        \endlang

        \identifier{o6} */
    virtual void init() {}

    /** \brief Finalize the object

        Called from construct() after init(), as the last step before the Function can be
        evaluated. The default does nothing.

        \identifier{o7} */
    virtual void finalize() {}

    /** \brief Evaluate numerically

        The method a subclass must implement. Receives one DM per input (with the sparsity
        of get_sparsity_in) and must return a list with one entry per output; each entry is
        converted to a DM and must match the shape of get_sparsity_out, otherwise the call
        raises an error. Exceptions raised inside propagate to the caller. The default
        implementation raises an error.

        This signature is not thread-safe. For guaranteed thread-safety, use eval_buffer.

        \param arg Input values
        \return Output values

        \lang{python}
        \doctest
        >>> class Scale(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_n_in(self): return 2
        ...   def get_sparsity_in(self, i):
        ...     return Sparsity.dense(2, 1) if i == 0 else Sparsity.dense(1, 1)
        ...   def get_sparsity_out(self, i): return Sparsity.dense(2, 1)
        ...   def eval(self, arg):
        ...     return [arg[0]*arg[1]]
        >>> f = Scale("f")
        >>> f([1, 2], 3)
        DM([3, 6])
        >>> class Bad(Scale):
        ...   def eval(self, arg):
        ...     return [arg[1]]
        >>> try:
        ...   Bad("b")([1, 2], 3)
        ... except Exception as e:
        ...   print("Shape mismatch" in str(e))
        True
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef Scale < casadi.Callback', 'methods'};
        >> c = [c, {'function self = Scale(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function n = get_n_in(self)', 'n = 2;', 'end'}];
        >> c = [c, {'function sp = get_sparsity_in(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(1, 1);'}];
        >> c = [c, {'if i == 0, sp = casadi.Sparsity.dense(2, 1); end', 'end'}];
        >> c = [c, {'function sp = get_sparsity_out(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(2, 1);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}*arg{2}};'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'Scale.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = Scale('f'); disp(f([1; 2], 3))
        [3, 6]
        >> c = {'classdef BadScale < Scale', 'methods', 'function self = BadScale(name)'};
        >> c = [c, {'self@Scale(name);', 'end', 'function r = eval(self, arg)', 'r = {arg{2}};'}];
        >> c = [c, {'end', 'end', 'end'}]; fid = fopen(fullfile(d, 'BadScale.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid); g = BadScale('b');
        >> try, g([1; 2], 3); catch e, disp(contains(e.message, 'Shape mismatch')); end
        1
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> abstract type AbstractScale <: ca.Callback end
        julia> struct Scale <: AbstractScale end
        julia> ca.get_n_in(::AbstractScale) = 2
        julia> ca.get_sparsity_in(::AbstractScale, i) = ca.Sparsity.dense(i == 0 ? 2 : 1, 1)
        julia> ca.get_sparsity_out(::AbstractScale, i) = ca.Sparsity.dense(2, 1)
        julia> ca.eval_(::Scale, arg) = [arg[1] * arg[2]]
        julia> f = ca.callback(Scale(), "f");
        julia> f([1, 2], 3)
        DM([3, 6])
        julia> struct Bad <: AbstractScale end
        julia> ca.eval_(::Bad, arg) = [arg[2]]
        julia> try
                   ca.callback(Bad(), "b")([1, 2], 3)
               catch e
                   println(occursin("Shape mismatch", sprint(showerror, e)))
               end
        true
        \enddoctest
        \endlang

        \identifier{o8} */
    virtual std::vector<DM> eval(const std::vector<DM>& arg) const;

    /** \brief A copy-free low level interface

        Alternative to eval() that works on the raw nonzero buffers, avoiding DM allocation.
        Override has_eval_buffer() to return true to activate it. In Python the method
        receives two tuples of memoryview objects (inputs and outputs), each holding only
        the structural nonzeros of the corresponding input or output in column-major order;
        write the result into the output buffers and return 0.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_sparsity_in(self, i): return Sparsity.dense(3, 1)
        ...   def get_sparsity_out(self, i): return Sparsity.dense(3, 1)
        ...   def has_eval_buffer(self): return True
        ...   def eval_buffer(self, arg, res):
        ...     a = np.frombuffer(arg[0], dtype=np.float64)
        ...     r = np.frombuffer(res[0], dtype=np.float64)
        ...     r[:] = a**2
        ...     return 0
        >>> f = Sq("f")
        >>> f([1, 2, 3])
        DM([1, 4, 9])
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.get_sparsity_in(::Sq, i) = ca.Sparsity.dense(3, 1)
        julia> ca.get_sparsity_out(::Sq, i) = ca.Sparsity.dense(3, 1)
        julia> ca.has_eval_buffer(::Sq) = true
        julia> ca.eval_buffer(::Sq, arg, res) = (res[1] .= arg[1] .^ 2; 0)
        julia> f = ca.callback(Sq(), "f");
        julia> f([1, 2, 3])
        DM([1, 4, 9])
        \enddoctest
        \endlang

        \identifier{o9} */
    virtual int eval_buffer(const double **arg, const std::vector<casadi_int>& sizes_arg,
                              double **res, const std::vector<casadi_int>& sizes_res) const;
    /** \brief Does the Callback class support a copy-free low level interface ?

        Return true to have eval_buffer() called instead of eval(). Default: false.

        \identifier{265} */
    virtual bool has_eval_buffer() const;

    /** \brief Get the number of inputs

        Called during construct(). Default: 1.

        \identifier{oa} */
    virtual casadi_int get_n_in();

    /** \brief Get the number of outputs

        Called during construct(). Default: 1.

        \lang{python}
        \doctest
        >>> class Split(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_n_in(self): return 1
        ...   def get_n_out(self): return 2
        ...   def eval(self, arg):
        ...     return [arg[0]**2, -arg[0]]
        >>> f = Split("f")
        >>> f(3)
        (DM(9), DM(-3))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SplitCb < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SplitCb(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function n = get_n_out(self)', 'n = 2;', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2, -arg{1}};'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SplitCb.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SplitCb('f'); [a, b] = f(3); disp(a), disp(b)
        9
        -3
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Split <: ca.Callback end
        julia> ca.get_n_in(::Split) = 1
        julia> ca.get_n_out(::Split) = 2
        julia> ca.eval_(::Split, arg) = [arg[1]^2, -arg[1]]
        julia> f = ca.callback(Split(), "f");
        julia> f(3)
        (9, -3)
        \enddoctest
        \endlang

        \identifier{ob} */
    virtual casadi_int get_n_out();

    /** \brief Get the sparsity of an input

        Called during construct() for i in 0..n_in-1. Default: scalar (1-by-1 dense).
        Inputs are padded or projected to this pattern when the Function is called.

        \lang{python}
        \doctest
        >>> class Trace(Callback):
        ...   def __init__(self, name, n, opts={}):
        ...     Callback.__init__(self)
        ...     self.n = n
        ...     self.construct(name, opts)
        ...   def get_sparsity_in(self, i): return Sparsity.diag(self.n)
        ...   def eval(self, arg):
        ...     return [sum1(diag(arg[0]))]
        >>> f = Trace("f", 3)
        >>> f.sparsity_in(0), f(DM.ones(3, 3))
        (Sparsity(3x3,3nz), DM(3))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef TraceCb < casadi.Callback', 'properties', 'n', 'end', 'methods'};
        >> c = [c, {'function self = TraceCb(name, n)'}];
        >> c = [c, {'self@casadi.Callback(); self.n = n; construct(self, name);', 'end'}];
        >> c = [c, {'function sp = get_sparsity_in(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.diag(self.n);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {sum1(diag(arg{1}))};'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'TraceCb.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = TraceCb('f', 3); disp(f.sparsity_in(0)), disp(f(DM.ones(3, 3)))
        3x3,3nz
        3
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Trace <: ca.Callback
                   n::Int
               end
        julia> ca.get_sparsity_in(cb::Trace, i) = ca.Sparsity.diag(cb.n)
        julia> ca.eval_(::Trace, arg) = [ca.sum1(ca.diag(arg[1]))]
        julia> f = ca.callback(Trace(3), "f");
        julia> ca.sparsity_in(f, 0), f(ca.DM.ones(3, 3))
        (3x3,3nz, 3)
        \enddoctest
        \endlang

        \identifier{oc} */
    virtual Sparsity get_sparsity_in(casadi_int i);

    /** \brief Get the sparsity of an output

        Called during construct() for i in 0..n_out-1. Default: scalar (1-by-1 dense).
        The value returned by eval() for output i must have this shape.

        \identifier{od} */
    virtual Sparsity get_sparsity_out(casadi_int i);

    /** \brief Get the name of an input

        Called during construct(). Default: "i0", "i1", ... Names allow calling the Function
        with keyword arguments.

        \lang{python}
        \doctest
        >>> class Lin(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_n_in(self): return 2
        ...   def get_name_in(self, i): return ["x", "a"][i]
        ...   def get_name_out(self, i): return "y"
        ...   def eval(self, arg):
        ...     return [arg[1]*arg[0]]
        >>> f = Lin("f")
        >>> f
        Callback(f:(x,a)->(y) CallbackInternal)
        >>> f(x=2, a=3)
        {'y': DM(6)}
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef Lin < casadi.Callback', 'methods'};
        >> c = [c, {'function self = Lin(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function n = get_n_in(self)', 'n = 2;', 'end'}];
        >> c = [c, {'function s = get_name_in(self, i)'}];
        >> c = [c, {'names = {''x'', ''a''}; s = names{i+1};', 'end'}];
        >> c = [c, {'function s = get_name_out(self, i)', 's = ''y'';', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{2}*arg{1}};'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'Lin.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = Lin('f'); disp(f)
        f:(x,a)->(y) CallbackInternal
        >> r = f('x', 2, 'a', 3); disp(r.y)
        6
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Lin <: ca.Callback end
        julia> ca.get_n_in(::Lin) = 2
        julia> ca.get_name_in(::Lin, i) = ["x", "a"][i + 1]
        julia> ca.get_name_out(::Lin, i) = "y"
        julia> ca.eval_(::Lin, arg) = [arg[2] * arg[1]]
        julia> f = ca.callback(Lin(), "f");
        julia> ca.Function(f)
        Function(f:(x,a)->(y) CallbackInternal)
        julia> f(x=2, a=3)
        Dict{String, DM} with 1 entry:
          "y" => 6
        \enddoctest
        \endlang

        \identifier{oe} */
    virtual std::string get_name_in(casadi_int i);

    /** \brief Get the name of an output

        Called during construct(). Default: "o0", "o1", ...

        \identifier{of} */
    virtual std::string get_name_out(casadi_int i);

    /** \brief Do the derivative functions need nondifferentiated outputs?

        Derivative functions obtained from get_jacobian/get_forward/get_reverse receive the
        nominal outputs as extra inputs. When this returns false (default) those inputs are
        passed as empty (all-zero sparsity) placeholders; return true to receive the actual
        output values, e.g. to reuse them in the derivative computation.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        ...   def has_jacobian(self): return True
        ...   def get_jacobian(self, name, inames, onames, opts):
        ...     x = SX.sym("x")
        ...     y = SX.sym("y", self.sparsity_out(0) if self.uses_output() else Sparsity(1, 1))
        ...     return Function(name, [x, y], [2*x], inames, onames, opts)
        >>> class SqU(Sq):
        ...   def uses_output(self): return True
        >>> Sq("f").jacobian().sparsity_in(1), SqU("f").jacobian().sparsity_in(1)
        (Sparsity(1x1,0nz), Sparsity(1x1))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqJac < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqJac(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end'}];
        >> c = [c, {'function b = has_jacobian(self)', 'b = true;', 'end'}];
        >> c = [c, {'function J = get_jacobian(self, name, inames, onames, opts)'}];
        >> c = [c, {'x = casadi.SX.sym(''x''); y = casadi.SX.sym(''y'', casadi.Sparsity(1, 1));'}];
        >> c = [c, {'if self.uses_output(), y = casadi.SX.sym(''y'', self.sparsity_out(0)); end'}];
        >> c = [c, {'J = casadi.Function(name, {x, y}, {2*x}, inames, onames, opts);'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqJac.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> c = {'classdef SqJacU < SqJac', 'methods', 'function self = SqJacU(name)'};
        >> c = [c, {'self@SqJac(name);', 'end', 'function b = uses_output(self)', 'b = true;'}];
        >> c = [c, {'end', 'end', 'end'}]; fid = fopen(fullfile(d, 'SqJacU.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqJac('f'); J = f.jacobian(); disp(J.sparsity_in(1))
        1x1,0nz
        >> g = SqJacU('g'); J = g.jacobian(); disp(J.sparsity_in(1))
        1x1
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> ca.has_jacobian(::Sq) = true
        julia> function ca.get_jacobian(::Sq, name, inames, onames, opts)
                   x = ca.SX.sym("x"); y = ca.SX.sym("y", ca.Sparsity(1, 1))
                   ca.Function(name, [x, y], [2x], inames, onames, opts)
               end
        julia> ca.sparsity_in(ca.jacobian(ca.callback(Sq(), "f")), 1)
        Sparsity(1x1,0nz)
        \enddoctest
        \endlang

        \identifier{og} */
    virtual bool uses_output() const;

    /** \brief Customize calls to the function factory

        Called when a derived function is requested by name through Function::factory;
        the default delegates to the standard factory.

        \identifier{2de} */
    virtual Function get_factory(const std::string& name,
        const std::vector<std::string>& s_in,
        const std::vector<std::string>& s_out,
        const Function::AuxOut& aux,
        const Dict& opts) const;

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements

        Override has_jacobian() to return true and get_jacobian() to return a Function with
        inputs [all nominal inputs, all nominal outputs] and one output per (output, input)
        pair, each the dense or sparse Jacobian block of that output with respect to that
        input (row-major over the pairs). Use the \p inames, \p onames and \p opts passed
        in: they carry the expected names and the "derivative_of" option. The returned
        Function may be an SX/MX Function (then higher derivatives are exact) or another
        Callback (keep a reference to it).

        \lang{python}
        \doctest
        >>> class F43(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_sparsity_in(self, i): return Sparsity.dense(4, 1)
        ...   def get_sparsity_out(self, i): return Sparsity.dense(3, 1)
        ...   def eval(self, arg):
        ...     a, b, c, d = vertsplit(arg[0])
        ...     return [vertcat(sin(c)*d+d**2, 2*a+c, b**2+5*c)]
        ...   def has_jacobian(self): return True
        ...   def get_jacobian(self, name, inames, onames, opts):
        ...     print(name, inames, onames)
        ...     x = SX.sym("x", 4)
        ...     y = SX.sym("y", Sparsity(3, 1))
        ...     a, b, c, d = vertsplit(x)
        ...     J = blockcat([[0, 0, d*cos(c), sin(c)+2*d], [2, 0, 1, 0], [0, 2*b, 5, 0]])
        ...     return Function(name, [x, y], [J], inames, onames, opts)
        >>> f = F43("f")
        >>> x = MX.sym("x", 4)
        >>> J = Function("J", [x], [jacobian(f(x), x)])
        jac_f ['i0', 'out_o0'] ['jac_o0_i0']
        >>> J([1, 2, 0, 3])
        DM(
        [[0, 0, 3, 6],
         [2, 0, 1, 0],
         [0, 4, 5, 0]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef F43 < casadi.Callback', 'methods'};
        >> c = [c, {'function self = F43(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function sp = get_sparsity_in(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(4, 1);', 'end'}];
        >> c = [c, {'function sp = get_sparsity_out(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(3, 1);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)'}];
        >> c = [c, {'v = arg{1}; a = v(1); b = v(2); c = v(3); d = v(4);'}];
        >> c = [c, {'r = {[sin(c)*d+d^2; 2*a+c; b^2+5*c]};', 'end'}];
        >> c = [c, {'function b = has_jacobian(self)', 'b = true;', 'end'}];
        >> c = [c, {'function J = get_jacobian(self, name, inames, onames, opts)'}];
        >> c = [c, {'disp(name); disp(inames); disp(onames);'}];
        >> c = [c, {'x = casadi.SX.sym(''x'', 4); a = x(1); b = x(2); c = x(3); d = x(4);'}];
        >> c = [c, {'y = casadi.SX.sym(''y'', casadi.Sparsity(3, 1));'}];
        >> c = [c, {'Jm = [0, 0, d*cos(c), sin(c)+2*d; 2, 0, 1, 0; 0, 2*b, 5, 0];'}];
        >> c = [c, {'J = casadi.Function(name, {x, y}, {Jm}, inames, onames, opts);'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'F43.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = F43('f'); x = MX.sym('x', 4); J = Function('J', {x}, {jacobian(f(x), x)});
        jac_f
        i0
        out_o0
        jac_o0_i0
        >> J([1; 2; 0; 3])
        [[0, 0, 3, 6],
         [2, 0, 1, 0],
         [0, 4, 5, 0]]
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct F43 <: ca.Callback end
        julia> ca.get_sparsity_in(::F43, i) = ca.Sparsity.dense(4, 1)
        julia> ca.get_sparsity_out(::F43, i) = ca.Sparsity.dense(3, 1)
        julia> function ca.eval_(::F43, arg)
                   a, b, c, d = arg[1][1], arg[1][2], arg[1][3], arg[1][4]
                   [[sin(c)*d + d^2; 2a + c; b^2 + 5c]]
               end
        julia> ca.has_jacobian(::F43) = true
        julia> function ca.get_jacobian(::F43, name, inames, onames, opts)
                   println(name, " ", inames, " ", onames)
                   x = ca.SX.sym("x", 4); y = ca.SX.sym("y", ca.Sparsity(3, 1))
                   a, b, c, d = x[1], x[2], x[3], x[4]
                   J = [0 0 d*cos(c) sin(c)+2d; 2 0 1 0; 0 2b 5 0]
                   ca.Function(name, [x, y], [J], inames, onames, opts)
               end
        julia> f = ca.callback(F43(), "f");
        julia> x = ca.MX.sym("x", 4);
        julia> J = ca.Function("J", [x], [ca.jacobian(f(x), x)]);
        jac_f ["i0", "out_o0"] ["jac_o0_i0"]
        julia> J([1, 2, 0, 3])
        DM(
        [[0, 0, 3, 6],
         [2, 0, 1, 0],
         [0, 4, 5, 0]])
        \enddoctest
        \endlang

        \identifier{oh} */
    virtual bool has_jacobian() const;
    virtual Function get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives

        Override has_forward(nfwd) to return true for the numbers of directions you support
        (nfwd==1 is enough; CasADi then loops) and get_forward() to return a Function with
        inputs [nominal inputs, nominal outputs, one forward seed per input] and one forward
        sensitivity per output. For nfwd>1 the seeds and sensitivities are horizontally
        stacked (n-by-nfwd*m). forward(nfwd) on the Function caches the result.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        ...   def has_forward(self, nfwd): return nfwd == 1
        ...   def get_forward(self, nfwd, name, inames, onames, opts):
        ...     print(nfwd, name, inames, onames)
        ...     x, y, xd = SX.sym("x"), SX.sym("y", Sparsity(1, 1)), SX.sym("xd")
        ...     return Function(name, [x, y, xd], [2*x*xd], inames, onames, opts)
        >>> f = Sq("f")
        >>> x = MX.sym("x")
        >>> J = Function("J", [x], [jacobian(f(x), x)])
        1 fwd1_f ['i0', 'out_o0', 'fwd_i0'] ['fwd_o0']
        >>> J(3)
        DM(6)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqFwd < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqFwd(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end'}];
        >> c = [c, {'function b = has_forward(self, nfwd)', 'b = nfwd == 1;', 'end'}];
        >> c = [c, {'function F = get_forward(self, nfwd, name, inames, onames, opts)'}];
        >> c = [c, {'disp(nfwd); disp(name); disp(inames); disp(onames);'}];
        >> c = [c, {'x = casadi.SX.sym(''x''); y = casadi.SX.sym(''y'', casadi.Sparsity(1, 1));'}];
        >> c = [c, {'xd = casadi.SX.sym(''xd'');'}];
        >> c = [c, {'F = casadi.Function(name, {x, y, xd}, {2*x*xd}, inames, onames, opts);'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqFwd.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqFwd('f'); x = MX.sym('x'); J = Function('J', {x}, {jacobian(f(x), x)});
        1
        fwd1_f
        i0
        out_o0
        fwd_i0
        fwd_o0
        >> J(3)
        6
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> ca.has_forward(::Sq, nfwd) = nfwd == 1
        julia> function ca.get_forward(::Sq, nfwd, name, inames, onames, opts)
                   println(nfwd, " ", name, " ", inames, " ", onames)
                   x = ca.SX.sym("x"); y = ca.SX.sym("y", ca.Sparsity(1, 1)); xd = ca.SX.sym("xd")
                   ca.Function(name, [x, y, xd], [2x*xd], inames, onames, opts)
               end
        julia> f = ca.callback(Sq(), "f");
        julia> x = ca.MX.sym("x");
        julia> J = ca.Function("J", [x], [ca.jacobian(f(x), x)]);
        1 fwd1_f ["i0", "out_o0", "fwd_i0"] ["fwd_o0"]
        julia> J(3)
        DM(6)
        \enddoctest
        \endlang

        \identifier{oi} */
    virtual bool has_forward(casadi_int nfwd) const;
    virtual Function get_forward(casadi_int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives

        Override has_reverse(nadj) to return true for the numbers of directions you support
        and get_reverse() to return a Function with inputs [nominal inputs, nominal outputs,
        one adjoint seed per output] and one adjoint sensitivity per input. For nadj>1 the
        seeds and sensitivities are horizontally stacked. Reverse mode is what gradient-based
        NLP solvers use for the objective gradient. reverse(nadj) caches the result.

        \lang{python}
        \doctest
        >>> class Sq(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        ...   def has_reverse(self, nadj): return nadj == 1
        ...   def get_reverse(self, nadj, name, inames, onames, opts):
        ...     print(nadj, name, inames, onames)
        ...     x, y, yb = SX.sym("x"), SX.sym("y", Sparsity(1, 1)), SX.sym("yb")
        ...     return Function(name, [x, y, yb], [2*x*yb], inames, onames, opts)
        >>> f = Sq("f")
        >>> x = MX.sym("x")
        >>> G = Function("G", [x], [gradient(f(x), x)])
        1 adj1_f ['i0', 'out_o0', 'adj_o0'] ['adj_i0']
        >>> G(3)
        DM(6)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef SqAdj < casadi.Callback', 'methods'};
        >> c = [c, {'function self = SqAdj(name)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}^2};', 'end'}];
        >> c = [c, {'function b = has_reverse(self, nadj)', 'b = nadj == 1;', 'end'}];
        >> c = [c, {'function F = get_reverse(self, nadj, name, inames, onames, opts)'}];
        >> c = [c, {'disp(nadj); disp(name); disp(inames); disp(onames);'}];
        >> c = [c, {'x = casadi.SX.sym(''x''); y = casadi.SX.sym(''y'', casadi.Sparsity(1, 1));'}];
        >> c = [c, {'yb = casadi.SX.sym(''yb'');'}];
        >> c = [c, {'F = casadi.Function(name, {x, y, yb}, {2*x*yb}, inames, onames, opts);'}];
        >> c = [c, {'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'SqAdj.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = SqAdj('f'); x = MX.sym('x'); G = Function('G', {x}, {gradient(f(x), x)});
        1
        adj1_f
        i0
        out_o0
        adj_o0
        adj_i0
        >> G(3)
        6
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Sq <: ca.Callback end
        julia> ca.eval_(::Sq, arg) = [arg[1]^2]
        julia> ca.has_reverse(::Sq, nadj) = nadj == 1
        julia> function ca.get_reverse(::Sq, nadj, name, inames, onames, opts)
                   println(nadj, " ", name, " ", inames, " ", onames)
                   x = ca.SX.sym("x"); y = ca.SX.sym("y", ca.Sparsity(1, 1)); yb = ca.SX.sym("yb")
                   ca.Function(name, [x, y, yb], [2x*yb], inames, onames, opts)
               end
        julia> f = ca.callback(Sq(), "f");
        julia> x = ca.MX.sym("x");
        julia> G = ca.Function("G", [x], [ca.gradient(f(x), x)]);
        1 adj1_f ["i0", "out_o0", "adj_o0"] ["adj_i0"]
        julia> G(3)
        DM(6)
        \enddoctest
        \endlang

        \identifier{oj} */
    virtual bool has_reverse(casadi_int nadj) const;
    virtual Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return sparsity of Jacobian of all input elements

        with respect to all output elements. Override has_jac_sparsity(oind, iind) to
        return true and get_jac_sparsity() to return the pattern of the Jacobian block of
        output \p oind with respect to input \p iind. Without it the block is assumed dense.
        Finite differencing exploits the pattern (fewer evaluations) and the result of
        jacobian() carries it.

        \lang{python}
        \doctest
        >>> class Diag(Callback):
        ...   def __init__(self, name, opts={}):
        ...     Callback.__init__(self)
        ...     self.construct(name, opts)
        ...   def get_sparsity_in(self, i): return Sparsity.dense(3, 1)
        ...   def get_sparsity_out(self, i): return Sparsity.dense(3, 1)
        ...   def eval(self, arg):
        ...     return [arg[0]**2]
        ...   def has_jac_sparsity(self, oind, iind): return True
        ...   def get_jac_sparsity(self, oind, iind, symmetric): return Sparsity.diag(3)
        >>> f = Diag("f", {"enable_fd": True})
        >>> x = MX.sym("x", 3)
        >>> J = jacobian(f(x), x)
        >>> J.sparsity()
        Sparsity(3x3,3nz)
        >>> Function("J", [x], [J])([1, 2, 3])
        DM(
        [[2, 00, 00],
         [00, 4, 00],
         [00, 00, 6]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> c = {'classdef DiagCb < casadi.Callback', 'methods'};
        >> c = [c, {'function self = DiagCb(name, opts)'}];
        >> c = [c, {'self@casadi.Callback(); construct(self, name, opts);', 'end'}];
        >> c = [c, {'function sp = get_sparsity_in(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(3, 1);', 'end'}];
        >> c = [c, {'function sp = get_sparsity_out(self, i)'}];
        >> c = [c, {'sp = casadi.Sparsity.dense(3, 1);', 'end'}];
        >> c = [c, {'function r = eval(self, arg)', 'r = {arg{1}.^2};', 'end'}];
        >> c = [c, {'function b = has_jac_sparsity(self, oind, iind)', 'b = true;', 'end'}];
        >> c = [c, {'function sp = get_jac_sparsity(self, oind, iind, symmetric)'}];
        >> c = [c, {'sp = casadi.Sparsity.diag(3);', 'end', 'end', 'end'}];
        >> d = tempname; mkdir(d); addpath(d); fid = fopen(fullfile(d, 'DiagCb.m'), 'w');
        >> fprintf(fid, '%s\n', c{:}); fclose(fid);
        >> f = DiagCb('f', struct('enable_fd', true)); x = MX.sym('x', 3); J = jacobian(f(x), x);
        >> J.sparsity()
        3x3,3nz
        >> JF = Function('J', {x}, {J}); disp(JF([1; 2; 3]))
        [[2, 00, 00],
         [00, 4, 00],
         [00, 00, 6]]
        >> rmpath(d); rmdir(d, 's');
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> struct Diag <: ca.Callback end
        julia> ca.get_sparsity_in(::Diag, i) = ca.Sparsity.dense(3, 1)
        julia> ca.get_sparsity_out(::Diag, i) = ca.Sparsity.dense(3, 1)
        julia> ca.eval_(::Diag, arg) = [arg[1] .^ 2]
        julia> ca.has_jac_sparsity(::Diag, oind, iind) = true
        julia> ca.get_jac_sparsity(::Diag, oind, iind, symmetric) = ca.Sparsity.diag(3)
        julia> f = ca.callback(Diag(), "f", Dict("enable_fd" => true));
        julia> x = ca.MX.sym("x", 3);
        julia> J = ca.jacobian(f(x), x);
        julia> ca.sparsity(J)
        Sparsity(3x3,3nz)
        julia> ca.Function("J", [x], [J])([1, 2, 3])
        DM(
        [[2, 00, 00],
         [00, 4, 00],
         [00, 00, 6]])
        \enddoctest
        \endlang

        \identifier{ok} */
    virtual bool has_jac_sparsity(casadi_int oind, casadi_int iind) const { return false; }
    virtual Sparsity get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const {
      return Sparsity(); }
    ///@}
  };

} // namespace casadi

#endif // CASADI_CALLBACK_HPP
