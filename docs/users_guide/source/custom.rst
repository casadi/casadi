.. _sec-user-defined:

User-defined function objects
=============================

.. include:: defs.rst

There are situations when rewriting user-functions using |casadi| symbolics is not
possible or practical. To tackle this, |casadi| provides a number of ways to
embed a call to a "black box" function defined in the language |casadi| is being
used from (C++, MATLAB or Python) or in C.
That being said, the recommendation is always to try to avoid this when possible,
even if it means investing a lot of time reimplementing existing code.
Functions defined using |casadi| symbolics are almost always more
efficient, especially when derivative calculation is involved, since a lot more
structure can typically be exploited.

Depending on the circumstances, the user can implement custom |Function|
objects in a number of different ways, which will be elaborated on in the following sections:


* Subclassing :class:`FunctionInternal`: :numref:`sec-function_internal`
* Subclassing :class:`Callback`: :numref:`sec-callback`
* Importing a function with ``external``: :numref:`sec-external`
* Just-in-time compile a C language string: :numref:`sec-jit_function`
* Replace the function call with a lookup table: :numref:`sec-lookup`

.. _sec-function_internal:

Subclassing |FunctionInternal|
------------------------------

All function objects presented in :numref:`Chapter %s <sec-function>` are implemented
in |casadi| as C++ classes inheriting from the |FunctionInternal| abstract
base class. In principle, a user with familiarity with C++ programming, can
implement a class inheriting from |FunctionInternal|,
overloading the virtual methods of this class. The best reference for doing so
is the C++ API documentation, choosing "switch to internal" to expose the internal
API.

Since |FunctionInternal| is not considered part of the stable, public API,
we advice against this in general, unless the plan is to contribution to |casadi|'s
source.

.. _sec-callback:

Subclassing |Callback|
----------------------
The |Callback| class provides a public API to |FunctionInternal|
and inheriting from this class has the same effect as inheriting directly from
|FunctionInternal|. Thanks to *cross-language polymorphism*, it
is possible to implement the exposed methods of |Callback| from either
Python, MATLAB/Octave or C++.

The derived class consists of the following parts:

* A constructor or a static function replacing the constructor
* A number of *virtual* functions, all optional, that can be overloaded in order to get the desired behavior. This includes the number of inputs and outputs using ``get_n_in`` and ``get_n_out``, their names using ``get_name_in`` and ``get_name_out`` and their sparsity patterns ``get_sparsity_in`` and ``get_sparsity_out``.
* An optional ``init`` function called duing the object construction.
* A function for numerical evaluation.
* Optional functions for derivatives. You can choose to work with finite differences (``enable_fd``), supply a full Jacobian (``has_jacobian``, ``get_jacobian``), or choose to supply forward/reverse sensitivities (``has_forward``, ``get_forward``,  ``has_reverse``, ``get_reverse``).

For a complete list of functions, see the C++ API documentation for |Callback|.
Also see the ``callback.py`` example.

The usage from the different languages are described in the following.

Python
^^^^^^

.. exec-block:: python

    class MyCallback(Callback):
      def __init__(self, name, d, opts={}):
        Callback.__init__(self)
        self.d = d
        self.construct(name, opts)

      # Number of inputs and outputs
      def get_n_in(self): return 1
      def get_n_out(self): return 1

      # Initialize the object
      def init(self):
         print('initializing object')

      # Evaluate numerically
      def eval(self, arg):
        x = arg[0]
        f = sin(self.d*x)
        return [f]

The implementation should include a constructor, which should begin with a call to
the base class constructor using
``Callback.__init__(self)`` and end with a call to
initialize object construction using ``self.construct(name, opts)``.

This function can be used as any built-in |casadi| function with the important
caveat that when embedded in graphs, the ownership of the class will *not*
be shared between all references. So it is important that the user does not
allow the Python class to go out of scope while it is still needed in
calculations.

.. exec-block:: python

    class MyCallback(Callback): [hidden]
      def __init__(self, name, d, opts={}): [hidden]
        Callback.__init__(self) [hidden]
        self.d = d [hidden]
        self.construct(name, opts) [hidden]

      # Number of inputs and outputs [hidden]
      def get_n_in(self): return 1 [hidden]
      def get_n_out(self): return 1 [hidden]

      # Initialize the object [hidden]
      def init(self): [hidden]
         print('initializing object') [hidden]

      # Evaluate numerically [hidden]
      def eval(self, arg): [hidden]
        x = arg[0] [hidden]
        f = sin(self.d*x) [hidden]
        return [f] [hidden]

    
    # Use the function
    f = MyCallback('f', 0.5)
    print(f(2))

    x = MX.sym("x")
    print(f(x))


MATLAB
^^^^^^
In MATLAB, a custom function class can be defined as follows, in a file
``MyCallback.m``:

.. code-block:: octave

    classdef MyCallback < casadi.Callback
      properties
        d
      end
      methods
        function self = MyCallback(name, d)
          self@casadi.Callback();
          self.d = d;
          construct(self, name);
        end

        % Number of inputs and outputs
        function v=get_n_in(self)
          v=1;
        end
        function v=get_n_out(self)
          v=1;
        end

        % Initialize the object
        function init(self)
          disp('initializing object')
        end

        % Evaluate numerically
        function arg = eval(self, arg)
          x = arg{1};
          f = sin(self.d * x);
          arg = {f};
        end
      end
    end

This function can be used as any built-in |casadi| function, but as for Python,
the ownership of the class will *not* be shared between all references.
So the user must not allow a class instance to get deleted while it is still
in use, e.g. by making it ``persistent``.

..
    f = fopen('MyCallback.m','w'); [hidden]
    fprintf(f,'classdef MyCallback < casadi.Callback\n'); [hidden]
    fprintf(f,'  properties\n'); [hidden]
    fprintf(f,'    d\n'); [hidden]
    fprintf(f,'  end\n'); [hidden]
    fprintf(f,'  methods\n'); [hidden]
    fprintf(f,'    function self = MyCallback(name, d)\n'); [hidden]
    fprintf(f,'      self@casadi.Callback();\n'); [hidden]
    fprintf(f,'      self.d = d;\n'); [hidden]
    fprintf(f,'      construct(self, name);\n'); [hidden]
    fprintf(f,'    end [hidden]\n');
    fprintf(f,'\n');
    fprintf(f,'    %% Number of inputs and outputs\n'); [hidden]
    fprintf(f,'    function v=get_n_in(self)\n'); [hidden]
    fprintf(f,'      v=1; [hidden]\n');
    fprintf(f,'    end [hidden]\n');
    fprintf(f,'    function v=get_n_out(self)\n'); [hidden]
    fprintf(f,'      v=1; [hidden]\n');
    fprintf(f,'    end [hidden]\n');
    fprintf(f,'\n');
    fprintf(f,'    %% Initialize the object\n'); [hidden]
    fprintf(f,'    function init(self)\n'); [hidden]
    fprintf(f,'      disp(''initializing object'')\n'); [hidden]
    fprintf(f,'    end\n'); [hidden]
    fprintf(f,'\n');
    fprintf(f,'    %% Evaluate numerically\n'); [hidden]
    fprintf(f,'    function arg = eval(self, arg)\n'); [hidden]
    fprintf(f,'      x = arg{1};\n'); [hidden]
    fprintf(f,'      f = sin(self.d * x);\n'); [hidden]
    fprintf(f,'      arg = {f};\n'); [hidden]
    fprintf(f,'    end\n'); [hidden]
    fprintf(f,'  end\n'); [hidden]
    fprintf(f,'end\n'); [hidden]
    fclose(f); [hidden]

.. code-block:: octave

    % Use the function
    f = MyCallback('f', 0.5);
    res = f(2);
    disp(res)

    x = MX.sym('x');
    disp(f(x))

C++
^^^

In C++, the syntax is as follows:

.. code-block:: cpp

    #include "casadi/casadi.hpp"
    using namespace casadi;
    class MyCallback : public Callback {
      // Data members
      double d;
    public:
      // Constructor
      MyCallback(const std::string& name, double d,
                 const Dict& opts=Dict()) : d(d) {
        construct(name, opts);
      }

      // Destructor
      ~MyCallback() override {}

      // Number of inputs and outputs
      casadi_int get_n_in() override { return 1;}
      casadi_int get_n_out() override { return 1;}

      // Initialize the object
      void init override() {
        std::cout << "initializing object" << std::endl;
      }

      // Evaluate numerically
      std::vector<DM> eval(const std::vector<DM>& arg) const override {
        DM x = arg.at(0);
        DM f = sin(d*x);
        return {f};
      }
    };

A class created this way can be used as any other |Function| instance,
but with the important difference that the user is responsible to managing
the memory of this class.

.. code-block:: cpp

  int main() {
    MyCallback f("f", 0.5);
    std::vector<DM> arg = {2};
    std::vector<DM> res = f(arg);
    std::cout << res << std::endl;
    return 0;
  }


.. _sec-external:


Importing a function with ``external``
--------------------------------------

The basic usage of |casadi|'s ``external`` function was demonstrated in
:numref:`sec-using_codegen` in the context of using autogenerated code. The
same function can also be used for importing a user-defined function, as long as
it also uses the C API described in :numref:`sec-c_api`.

The following sections expands on this.

Default functions
^^^^^^^^^^^^^^^^^

It is usually *not* necessary to define all the functions defined in
:numref:`sec-c_api`. If ``fname_incref`` and ``fname_decref``
are absent, it is assumed that no memory management is needed. If no
names of inputs and outputs are provided, they will be given default names.
Sparsity patterns are in general assumed to be scalar by default, unless the
function corresponds to a derivative of another function (see below), in which
case they are assumed to be dense and of the correct dimension.

Furthermore, work vectors are assumed not to be needed if ``fname_work`` has
not been implemented.

Meta information as comments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you rely on |casadi|'s just-in-time compiler, you can provide meta information
as a comment in the C code instead of implementing the actual callback function.

The structure of such meta information should be as follows:

.. code-block:: c

    /*CASADIMETA
    :fname_N_IN 1
    :fname_N_OUT 2
    :fname_NAME_IN[0] x
    :fname_NAME_OUT[0] r
    :fname_NAME_OUT[1] s
    :fname_SPARSITY_IN[0] 2 1 0 2
    */

Derivatives
^^^^^^^^^^^

The external function can be made differentiable by providing functions for
calculating derivatives. During derivative calculations, |casadi| will look for
symbols in the same file/shared library that follows a certain
*naming convention*. For example, you can specify a Jacobian for all the
outputs with respect to all inputs for a function named ``fname`` by
implementing a function named ``jac_fname``. Similary, you can specify
a function for calculating one forward directional derivative by providing a
function named ``fwd1_fname``, where 1 can be replaced by 2, 4, 8, 16,
32 or 64 for calculating multiple forward directional derivatives at once.
For reverse mode directional derivatives, replace ``fwd`` with ``adj``.

This is an experimental feature.

.. _sec-jit_function:

Just-in-time compile a C language string
----------------------------------------

In the previous section we showed how to specify a C file with functions for numerical
evaluation and meta information. As was shown, this file can be just-in-time
compiled by CasADi's interface to Clang. There exists a shorthand for this approach,
where the user simply specifies the source code as a C language string.

.. side-by-side::
    .. exec-block:: python

        body =\
        'r[0] = x[0];'+\
        'while (r[0]<s[0]) {'+\
        ' r[0] *= r[0];'+\
        '}'

        f = Function.jit('f',body,\
              ['x','s'],['r'])
        print(f)
    &&

    .. exec-block:: octave

        body =[...
        'r[0] = x[0];',...
        'while (r[0]<s[0]) {',...
        ' r[0] *= r[0];',...
        '}'];

        f = Function.jit('f',body,...
              {'x','s'},{'r'});
        disp(f)


These four arguments of ``Function.jit`` are mandatory:
The function name, the C source as a string and the names of inputs and outputs.
In the C source, the input/output names correspond to arrays of type ``casadi_real_t``
containing the nonzero elements of the function inputs and outputs. By default,
all inputs and outputs are scalars (i.e. 1-by-1 and dense). To specify a different sparsity pattern, provide two additional function arguments containing vectors/lists
of the sparsity patterns:

.. side-by-side::

    .. exec-block:: python

        body =\ [hidden]
        'r[0] = x[0];'+\ [hidden]
        'while (r[0]<s[0]) {'+\ [hidden]
        ' r[0] *= r[0];'+\ [hidden]
        '}\n' [hidden]

        sp = Sparsity.scalar()
        f = Function.jit('f',body,\
             ['x','s'], ['r'],\
             [sp,sp], [sp])
        print(f)
    &&

    .. exec-block:: octave

        body =[... [hidden]
        'r[0] = x[0];',... [hidden]
        'while (r[0]<s[0]) {',... [hidden]
        ' r[0] *= r[0];',... [hidden]
        '}']; [hidden]

        sp = Sparsity.scalar();
        f = Function.jit('f',body,...
             {'x','s'}, {'r'},...
             {sp,sp}, {sp});
        disp(f)


Both variants accept an optional 5th (or 7th) argument in the form of an
options dictionary.


.. _sec-lookup:

Using lookup-tables
-------------------
Lookup-tables can be created using |casadi|'s ``interpolant`` function. Different interpolating schemes are implemented as *plugins*, similar to ``nlpsol`` or ``integrator`` objects. In addition to the identifier name and plugin, the ``interpolant`` function expects a set of grid points with the corresponding numerical values.

The result of an ``interpolant`` call is a |casadi| Function object that is differentiable, and can be embedded into |casadi| computational graphs by calling with MX arguments. Furthermore, C code generation is fully supported for such graphs.

Currently, two plugins exist for ``interpolant``: ``'linear'`` and ``'bspline'``. They are intended to behave simiarly to MATLAB/Octave's ``interpn`` with the method set to ``'linear'`` or ``'spline'`` -- corresponding to a multilinear interpolation and a (by default cubic) spline interpolation with
not-a-knot boundary conditions.

In the case of ``bspline``, coefficients will be sought at construction time that fit the provided data. Alternatively, you may also use the more low-level ``Function.bspline`` to supply the coefficients yourself. The default degree of the bspline is 3 in each dimension. You may deviate from this default by passing a ``degree`` option.

We will walk through the syntax of ``interpolant`` for the 1D and 2D versions, but the syntax in fact generalizes to an arbitrary number of dimensions.

1D lookup tables
^^^^^^^^^^^^^^^^
A 1D spline fit can be done in |casadi|/Python as follows, compared with the corresponding method in SciPy:

.. exec-block:: python

    import casadi as ca
    import numpy as np
    xgrid = np.linspace(1,6,6)
    V = [-1,-1,-2,-3,0,2]
    lut = ca.interpolant('LUT','bspline',[xgrid],V)
    print(lut(2.5))
    # Using SciPy
    import scipy.interpolate as ip
    interp = ip.InterpolatedUnivariateSpline(xgrid, V)
    print(interp(2.5))

In MATLAB/Octave, the corresponding code reads:

.. exec-block:: octave

    xgrid = 1:6;
    V = [-1 -1 -2 -3 0 2];
    lut = casadi.interpolant('LUT','bspline',{xgrid},V);
    lut(2.5)
    % Using MATLAB/Octave builtin
    interp1(xgrid,V,2.5,'spline')

Note in particular that the ``grid`` and ``values`` arguments to ``interpolant`` must be numerical in nature.

2D lookup tables
^^^^^^^^^^^^^^^^

In two dimensions, we get the following in Python, also compared to SciPy for reference:

.. exec-block:: python

    import casadi as ca [hidden]
    import numpy as np [hidden]
    import scipy.interpolate as ip [hidden]

    xgrid = np.linspace(-5,5,11)
    ygrid = np.linspace(-4,4,9)
    X,Y = np.meshgrid(xgrid,ygrid,indexing='ij')
    R = np.sqrt(5*X**2 + Y**2)+ 1
    data = np.sin(R)/R
    data_flat = data.ravel(order='F')
    lut = ca.interpolant('name','bspline',[xgrid,ygrid],data_flat)
    print(lut([0.5,1]))
    # Using Scipy

    interp = ip.RectBivariateSpline(xgrid, ygrid, data)
    print(interp.ev(0.5,1))


or, in MATLAB/Octave compared to the built-in functions:

.. exec-block:: octave

    xgrid = -5:1:5;
    ygrid = -4:1:4;
    [X,Y] = ndgrid(xgrid, ygrid);
    R = sqrt(5*X.^2 + Y.^2)+ 1;
    V = sin(R)./R;
    lut = interpolant('LUT','bspline',{xgrid, ygrid},V(:));
    lut([0.5 1])
    % Using Matlab builtin
    interpn(X,Y,V,0.5,1,'spline')

In particular note how the ``values`` argument had to be flatten to a one-dimensional array.

.. _sec-fd:

Derivative calculation using finite differences
-----------------------------------------------

CasADi 3.3 introduced support for finite difference calculation for all function
objects, in particular including external functions defined as outlined in :numref:`sec-callback`, :numref:`sec-external` or :numref:`sec-jit_function`
(for lookup tables, :numref:`sec-lookup`, analytical derivatives are available).

Finite difference derivative are disabled by default, with the exception of ``Function.jit``, and to enable it, you must set the option ``'enable_fd'``
to ``True``/``true``:

.. side-by-side::
    .. exec-block:: python

        from os import system [hidden]

        x = MX.sym('x') [hidden]
        f = Function('f',[x],[sin(x)]) [hidden]
        f.generate('gen.c') [hidden]
        system('gcc -fPIC -shared gen.c -o gen.so') [hidden]
        f = external('f', './gen.so',\
           dict(enable_fd=True))

        e = jacobian(f(x),x)
        D = Function('d',[x],[e])
        print(D(0))

    &&

    .. exec-block:: octave

        x = MX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]
        f.generate('gen.c'); [hidden]
        system('gcc -fPIC -shared gen.c -o gen.so'); [hidden]
        f = external('f', './gen.so',...
            struct('enable_fd',true));

        e = jacobian(f(x),x);
        D = Function('d',{x},{e});
        disp(D(0))

cf. :numref:`sec-codegen_syntax`.

The ``'enable_fd'`` options enables |casadi| to use finite differences, *if*
analytical derivatives are not available. To force |casadi| to use finite differences,
you can set ``'enable_forward'``, ``'enable_reverse'`` and ``'enable_jacobian'`` to ``False``/``false``, corresponding to the three types of analytical derivative information that |casadi| works with.

The default method is central differences with a step size determined by estimates of roundoff errors and truncation errors of the function. You can change the method by setting the option ``'fd_method'`` to ``'forward'`` (corresponding to first order forward differences), ``'backward'`` (corresponding to first order backward differences) and ``'smoothing'`` for a second-order accurate discontinuity avoiding scheme, suitable when derivatives need to be calculated at the edges of a domain. Additional algorithmic options for the finite differences are available by setting ``'fd_options'`` option.

