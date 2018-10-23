.. _sec-function:

Function objects
================

.. include:: defs.rst

|casadi| allows the user to create function objects, in C++ terminology often referred to as *functors*. This includes functions that are defined by a symbolic expression, ODE/DAE integrators, QP solvers, NLP solvers etc.

Function objects are typically created with the syntax:

.. code-block:: none

      f = functionname(name, arguments, ..., [options])

The name is mainly a display name that will show up in e.g. error messages or as comments in generated C code. This is followed by a set of arguments, which is class dependent. Finally, the user can pass an options structure for customizing the behavior of the class. The options structure is a dictionary type in Python, a struct in MATLAB or |casadi|'s ``Dict`` type in C++.

A ``Function`` can be constructed by passing a list of input expressions and a list of output expressions:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym('x',2)
        y = SX.sym('y')
        f = Function('f',[x,y],\
                   [x,sin(y)*x])
        print(f)
    &&

    .. exec-block:: octave

        x = SX.sym('x',2);
        y = SX.sym('y');
        f = Function('f',{x,y},...
                   {x,sin(y)*x});
        disp(f)

which defines a function
:math:`f : \mathbb{R}^{2} \times \mathbb{R} \rightarrow \mathbb{R}^{2} \times \mathbb{R}^{2}, \quad (x,y) \mapsto (x,\sin(y) x)`.
Note that all function objects in |casadi|, including the above, are multiple matrix-valued input, multiple, matrix-valued output.

|MX| expression graphs work the same way:

.. side-by-side::
    .. exec-block:: python

        x = MX.sym('x',2)
        y = MX.sym('y')
        f = Function('f',[x,y],\
                   [x,sin(y)*x])
        print(f)
    &&

    .. exec-block:: octave

        x = MX.sym('x',2);
        y = MX.sym('y');
        f = Function('f',{x,y},...
                   {x,sin(y)*x});
        disp(f)

When creating a ``Function`` from expressions like that, it is always advisory to *name* the inputs and outputs as follows:


.. side-by-side::
    .. exec-block:: python

        x = MX.sym('x',2) [hidden]
        y = MX.sym('y') [hidden]
        f = Function('f',[x,y],\
              [x,sin(y)*x],\
              ['x','y'],['r','q'])
        print(f)
    &&

    .. exec-block:: octave

        x = MX.sym('x',2); [hidden]
        y = MX.sym('y'); [hidden]
        f = Function('f',{x,y},...
              {x,sin(y)*x},...
              {'x','y'},{'r','q'});
        disp(f)

Naming inputs and outputs is preferred for a number of reasons:

* No need to remember the number or order of arguments
* Inputs or outputs that are absent can be left unset
* More readable and less error prone syntax.

For ``Function`` instances -- to be encountered later -- that are *not* created directly from expressions,
the inputs and outputs are named automatically.

Calling function objects
------------------------

|MX| expressions may contain calls to ``Function``-derived functions. Calling a function object is both done for the numerical evaluation and, by passing symbolic arguments, for embedding a *call* to the function object into an expression graph (cf. also :numref:`sec-integrator`).

To call a function object, you either pass the argument in the correct order:

.. side-by-side::
    .. exec-block:: python

        x = MX.sym('x',2) [hidden]
        y = MX.sym('y') [hidden]
        f = Function('f',[x,y],[x,sin(y)*x],['x','y'],['r','q']) [hidden]
        r0, q0 = f(1.1,3.3)
        print('r0:',r0)
        print('q0:',q0)
    &&

    .. exec-block:: octave

        x = MX.sym('x',2); [hidden]
        y = MX.sym('y'); [hidden]
        f = Function('f',{x,y},{x,sin(y)*x},{'x','y'},{'r','q'}); [hidden]
        [r0, q0] = f(1.1,3.3);
        disp(r0)
        disp(q0)

or the arguments and their names as follows, which will result in a dictionary (``dict`` in Python, ``struct`` in MATLAB and ``std::map<std::string, MatrixType>`` in C++):

.. side-by-side::
    .. exec-block:: python

        x = MX.sym('x',2) [hidden]
        y = MX.sym('y') [hidden]
        f = Function('f',[x,y],[x,sin(y)*x],['x','y'],['r','q']) [hidden]
        res = f(x=1.1, y=3.3)
        print('res:', res)
    &&

    .. exec-block:: octave

        x = MX.sym('x',2); [hidden]
        y = MX.sym('y'); [hidden]
        f = Function('f',{x,y},{x,sin(y)*x},{'x','y'},{'r','q'}); [hidden]
        res = f('x',1.1,'y',3.3);
        disp(res)

When calling a function object, the dimensions (but not necessarily the sparsity patterns) of the evaluation arguments have to match those of the function inputs, with two exceptions:

* A row vector can be passed instead of a column vector and vice versa.
* A scalar argument can always be passed, regardless of the input dimension. This has the meaning of setting all elements of the input matrix to that value.

When the number of inputs to a function object is large or changing, an alternative syntax to the above is to use the *call* function which takes a Python list / MATLAB cell array or, alternatively, a Python dict / MATLAB struct. The return value will have the same type:

.. side-by-side::
    .. exec-block:: python

        x = MX.sym('x',2) [hidden]
        y = MX.sym('y') [hidden]
        f = Function('f',[x,y],[x,sin(y)*x],['x','y'],['r','q']) [hidden]
        arg = [1.1,3.3]
        res = f.call(arg)
        print('res:', res)
        arg = {'x':1.1,'y':3.3}
        res = f.call(arg)
        print('res:', res)
    &&

    .. exec-block:: octave

        x = MX.sym('x',2); [hidden]
        y = MX.sym('y'); [hidden]
        f = Function('f',{x,y},{x,sin(y)*x},{'x','y'},{'r','q'}); [hidden]

        arg = {1.1,3.3};
        res = f.call(arg);
        disp(res)
        arg = struct('x',1.1,'y',3.3);
        res = f.call(arg);
        disp(res)

Converting |MX| to |SX|
-----------------------
A function object defined by an |MX| graph that only contains built-in operations (e.g. element-wise operations such as addition, square root, matrix multiplications and calls to |SX| functions, can be converted into a function defined purely by an |SX| graph using the syntax:

.. code-block:: python

     sx_function = mx_function.expand()


This might speed up the calculations significantly, but might also cause extra memory overhead.


.. _sec-rootfinder:

Nonlinear root-finding problems
-------------------------------

Consider the following system of equations:

.. math::

    \begin{aligned}
    &g_0(z, x_1, x_2, \ldots, x_n) &&= 0 \\
    &g_1(z, x_1, x_2, \ldots, x_n) &&= y_1 \\
    &g_2(z, x_1, x_2, \ldots, x_n) &&= y_2 \\
    &\qquad \vdots \qquad &&\qquad \\
    &g_m(z, x_1, x_2, \ldots, x_n) &&= y_m,
    \end{aligned}

where the first equation uniquely defines :math:`z` as a function of :math:`x_1`, \ldots, :math:`x_n` by the *implicit function theorem*
and the remaining equations define the auxiliary outputs :math:`y_1`, \ldots, :math:`y_m`.

Given a function :math:`g` for evaluating :math:`g_0`, \ldots, :math:`g_m`, we can use |casadi| to automatically formulate a function
:math:`G: \{z_{\text{guess}}, x_1, x_2, \ldots, x_n\} \rightarrow \{z, y_1, y_2, \ldots, y_m\}`.
This function includes a guess for :math:`z` to handle the case when the solution is non-unique.
The syntax for this, assuming :math:`n=m=1` for simplicity, is:


.. side-by-side::
    .. exec-block:: python

        nz = 1 [hidden]
        nx = 1 [hidden]

        z = SX.sym('x',nz)
        x = SX.sym('x',nx)
        g0 = sin(x+z)
        g1 = cos(x-z)
        g = Function('g',[z,x],[g0,g1])
        G = rootfinder('G','newton',g)
        print(G)
    &&

    .. exec-block:: octave

        nz = 1; [hidden]
        nx = 1; [hidden]

        z = SX.sym('x',nz);
        x = SX.sym('x',nx);
        g0 = sin(x+z);
        g1 = cos(x-z);
        g = Function('g',{z,x},{g0,g1});
        G = rootfinder('G','newton',g);
        disp(G)

where the ``rootfinder`` function expects a display name, the name of a solver plugin
(here a simple full-step Newton method) and the residual function.

Rootfinding objects in |casadi| are differential objects and derivatives can be calculated exactly to arbitrary order.

.. seealso:: `API of rootfinder <http://casadi.org/python-api/#rootfinding>`_


.. _sec-integrator:

Initial-value problems and sensitivity analysis
-----------------------------------------------

|casadi| can be used to solve initial-value problems in ODE or DAE. The problem formulation used
is a DAE of semi-explicit form with quadratures:

.. math::

    \begin{aligned}
     \dot{x} &= f_{\text{ode}}(t,x,z,p), \qquad x(0) = x_0 \\
          0  &= f_{\text{alg}}(t,x,z,p) \\
     \dot{q} &= f_{\text{quad}}(t,x,z,p), \qquad q(0) = 0
    \end{aligned}

For solvers of *ordinary* differential equations, the second equation and the algebraic variables :math:`z` must be absent.

An integrator in |casadi| is a function that takes the state at the initial time ``x0``, a set of parameters ``p``, and a guess for the algebraic variables (only for DAEs) ``z0`` and returns the state vector ``xf``, algebraic variables ``zf`` and the quadrature state ``qf``, all at the final time.

The freely available `SUNDIALS suite <https://computation.llnl.gov/casc/sundials/description/description.html>`_ (distributed along with |casadi|) contains the two popular integrators CVodes and IDAS for ODEs and DAEs respectively. These integrators have support for forward and adjoint sensitivity analysis and when used via |casadi|'s Sundials interface, |casadi| will automatically formulate the Jacobian information, which is needed by the backward differentiation formula (BDF) that CVodes and IDAS use. Also automatically formulated will be the forward and adjoint sensitivity equations.

Creating integrators
^^^^^^^^^^^^^^^^^^^^

Integrators are created using |casadi|'s ``integrator`` function. Different integrators schemes and interfaces are implemented as *plugins*, essentially shared libraries that are loaded at runtime.

Consider for example the DAE:

.. math::

  \begin{aligned}
   \dot{x} &= z+p, \\
        0  &= z \, \cos(z)-x
  \end{aligned}

An integrator, using the "idas" plugin, can be created using the syntax:

.. exec-block:: python

    x = SX.sym('x'); z = SX.sym('z'); p = SX.sym('p')
    dae = {'x':x, 'z':z, 'p':p, 'ode':z+p, 'alg':z*cos(z)-x}
    F = integrator('F', 'idas', dae)
    print(F)

.. exec-block:: octave

    x = SX.sym('x'); z = SX.sym('z'); p = SX.sym('p');
    dae = struct('x',x,'z',z,'p',p,'ode',z+p,'alg',z*cos(z)-x);
    F = integrator('F', 'idas', dae);
    disp(F)


Integrating this DAE from 0 to 1 with :math:`x(0)=0`, :math:`p=0.1` and using the guess :math:`z(0)=0`, can
be done by evaluating the created function object:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym('x'); z = SX.sym('z'); p = SX.sym('p') [hidden]
        dae = {'x':x, 'z':z, 'p':p, 'ode':z+p, 'alg':z*cos(z)-x} [hidden]
        F = integrator('F', 'idas', dae) [hidden]

        r = F(x0=0, z0=0, p=0.1)
        print(r['xf'])
    &&

    .. exec-block:: octave

        x = SX.sym('x'); z = SX.sym('z'); p = SX.sym('p'); [hidden]
        dae = struct('x',x,'z',z,'p',p,'ode',z+p,'alg',z*cos(z)-x); [hidden]
        F = integrator('F', 'idas', dae); [hidden]

        r = F('x0',0,'z0',0,'p',0.1);
        disp(r.xf)

The time horizon is assumed to be fixed [#f1]_ and can be changed from its default [0, 1] by setting the options "t0" and "tf".

Sensitivity analysis
^^^^^^^^^^^^^^^^^^^^

From a usage point of view, an integrator behaves just like the function objects created from expressions earlier in the chapter.
You can use member functions in the Function class to generate new function objects corresponding to directional derivatives (forward or reverse mode) or complete Jacobians. Then evaluate these function objects numerically to obtain sensitivity information. The documented example "sensitivity_analysis" (available in |casadi|'s example collection for Python, MATLAB and C++) demonstrate how |casadi| can be used to calculate first and second order derivative information (forward-over-forward, forward-over-adjoint, adjoint-over-adjoint) for a simple DAE.


.. _sec-nlpsol:

Nonlinear programming
---------------------

.. note:: This section assumes familiarity with much of the material that comes above.
             There is also a higher-level interface in :numref:`Chapter %s <sec-opti>`. That interface can be learned stand-alone.


The NLP solvers distributed with or interfaced to |casadi| solves parametric NLPs of the following form:

.. math::
    :label: nlp

    \begin{array}{cc}
    \begin{array}{c}
    \text{minimize:} \\
    x
    \end{array}
    &
    f(x,p)
    \\
    \begin{array}{c}
    \text{subject to:}
    \end{array}
    &
    \begin{array}{rcl}
      x_{\textrm{lb}} \le &  x   & \le x_{\textrm{ub}} \\
      g_{\textrm{lb}} \le &g(x,p)& \le g_{\textrm{ub}}
    \end{array}
    \end{array}


where :math:`x \in \mathbb{R}^{nx}` is the decision variable and :math:`p \in \mathbb{R}^{np}` is a known parameter vector.

An NLP solver in |casadi| is a function that takes the parameter value (``p``), the bounds (``lbx``, ``ubx``, ``lbg``, ``ubg``) and a guess for the primal-dual solution (``x0``, ``lam_x0``, ``lam_g0``) and returns the optimal solution. Unlike integrator objects, NLP solver functions are currently not differentiable functions in |casadi|.

There are several NLP solvers interfaced with |casadi|. The most popular one is IPOPT, an open-source primal-dual interior point method which is included in |casadi| installations. Others, that require the installation of third-party software, include SNOPT, WORHP and KNITRO. Whatever the NLP solver used, the interface will automatically generate the information that it needs to solve the NLP, which may be solver and option dependent. Typically an NLP solver will need a function that gives the Jacobian of the constraint function and a Hessian of the Lagrangian function (:math:`L(x,\lambda) = f(x) + \lambda^{\text{T}} \, g(x))` with respect to :math:`x`.

Creating NLP solvers
^^^^^^^^^^^^^^^^^^^^

NLP solvers are created using |casadi|'s ``nlpsol`` function. Different solvers and interfaces are implemented as *plugins*.
Consider the following form of the so-called Rosenbrock problem:

.. math::

    \begin{array}{cc}
    \begin{array}{c}
    \text{minimize:} \\
    x,y,z
    \end{array}
    &
    x^2 + 100 \, z^2  \\
    \begin{array}{c}
    \text{subject to:}
    \end{array}
    &  z+(1-x)^2-y = 0
    \end{array}


A solver for this problem, using the "ipopt" plugin, can be created using the syntax:

.. exec-block:: python

    x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z')
    nlp = {'x':vertcat(x,y,z), 'f':x**2+100*z**2, 'g':z+(1-x)**2-y}
    S = nlpsol('S', 'ipopt', nlp)
    print(S)

.. exec-block:: octave

    x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
    nlp = struct('x',[x;y;z], 'f',x^2+100*z^2, 'g',z+(1-x)^2-y);
    S = nlpsol('S', 'ipopt', nlp);
    disp(S)


Once the solver has been created, we can solve the NLP, using ``[2.5,3.0,0.75]`` as an initial guess, by evaluating the
function ``S``:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z') [hidden]
        nlp = {'x':vertcat(x,y,z), 'f':x**2+100*z**2, 'g':z+(1-x)**2-y} [hidden]
        S = nlpsol('S', 'ipopt', nlp) [hidden]

        r = S(x0=[2.5,3.0,0.75],\
              lbg=0, ubg=0)
        x_opt = r['x']
        print('x_opt: ', x_opt)


    &&

    .. exec-block:: octave

        x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z'); [hidden]
        nlp = struct('x',[x;y;z], 'f',x^2+100*z^2, 'g',z+(1-x)^2-y); [hidden]
        S = nlpsol('S', 'ipopt', nlp); [hidden]

        r = S('x0',[2.5,3.0,0.75],...
              'lbg',0,'ubg',0);
        x_opt = r.x;
        disp(x_opt)


.. _sec-qpsol:

Quadratic programming
---------------------

|casadi| provides interfaces to solve quadratic programs (QPs). Supported solvers are the open-source solvers qpOASES (distributed with |casadi|) and
OOQP as well as the commercial solvers CPLEX and GUROBI.

There are two different ways to solve QPs in |casadi|, using a high-level interface and a low-level interface. They are described in the following.

High-level interface
^^^^^^^^^^^^^^^^^^^^

The high-level interface for quadratic programming mirrors that of nonlinear programming, i.e. expects a problem of the form :eq:`nlp`,
with the restriction that objective function :math:`f(x,p)` must be a convex quadratic function in :math:`x` and the constraint function :math:`g(x,p)` must be linear in :math:`x`.
If the functions are not quadratic and linear, respectively, the solution is done at the current linearization point, given by the "initial guess" for :math:`x`.

If the objective function is not convex, the solver may or may not fail to find a solution or the solution may not be unique.

To illustrate the syntax, we consider the following convex QP:

.. math::
  :label: simple_qp

    \begin{array}{cc}
    \begin{array}{c}
    \text{minimize:} \\
    x,y
    \end{array}
    &
    x^2 + y^2  \\
    \begin{array}{c}
    \text{subject to:}
    \end{array}
    & x+y-10 \ge 0
    \end{array}


To solve this problem with the high-level interface, we simply replace ``nlpsol`` with ``qpsol`` and use a QP solver plugin such as the with |casadi| distributed qpOASES:

.. exec-block:: python

    x = SX.sym('x'); y = SX.sym('y')
    qp = {'x':vertcat(x,y), 'f':x**2+y**2, 'g':x+y-10}
    S = qpsol('S', 'qpoases', qp)
    print(S)

.. exec-block:: octave

    x = SX.sym('x'); y = SX.sym('y');
    qp = struct('x',[x;y], 'f',x^2+y^2, 'g',x+y-10);
    S = qpsol('S', 'qpoases', qp);
    disp(S)

The created solver object ``S`` will have the same input and output signature as the solver objects
created with ``nlpsol``. Since the solution is unique, it is less important to provide an initial guess:


.. side-by-side::
    .. exec-block:: python

        x = SX.sym('x'); y = SX.sym('y') [hidden]
        qp = {'x':vertcat(x,y), 'f':x**2+y**2, 'g':x+y-10} [hidden]
        S = qpsol('S', 'qpoases', qp) [hidden]

        r = S(lbg=0)
        x_opt = r['x']
        print('x_opt: ', x_opt)

    &&

    .. exec-block:: octave

        x = SX.sym('x'); y = SX.sym('y'); [hidden]
        qp = struct('x',[x;y], 'f',x^2+y^2, 'g',x+y-10); [hidden]
        S = qpsol('S', 'qpoases', qp); [hidden]

        r = S('lbg',0);
        x_opt = r.x;
        disp(x_opt)


Low-level interface
^^^^^^^^^^^^^^^^^^^

The low-level interface, on the other hand, solves QPs of the following form:

.. math::
  :label: qp

    \begin{array}{cc}
    \begin{array}{c}
    \text{minimize:} \\
    x
    \end{array}
    &
    \frac{1}{2} x^T \, H \, x + g^T \, x
    \\
    \begin{array}{c}
    \text{subject to:}
    \end{array}
    &
    \begin{array}{rcl}
      x_{\textrm{lb}} \le &  x   & \le x_{\textrm{ub}} \\
      a_{\textrm{lb}} \le & A \, x& \le a_{\textrm{ub}}
    \end{array}
    \end{array}


Encoding problem :eq:`simple_qp` in this form, omitting bounds that are infinite, is straightforward:


.. side-by-side::
    .. code-block:: python

        H = 2*DM.eye(2)
        A = DM.ones(1,2)
        g = DM.zeros(2)
        lba = 10.


    .. code-block:: octave

        H = 2*DM.eye(2);
        A = DM.ones(1,2);
        g = DM.zeros(2);
        lba = 10;

To create a solver instance, instead of passing symbolic expressions for the QP, we now pass the sparsity patterns of the matrices :math:`H` and :math:`A`.
Since we used |casadi|'s |DM|-type above, we can simply query the sparsity patterns:

.. side-by-side::
    .. exec-block:: python

        H = 2*DM.eye(2) [hidden]
        A = DM.ones(1,2) [hidden]
        g = DM.zeros(2) [hidden]
        lba = 10. [hidden]

        qp = {}
        qp['h'] = H.sparsity()
        qp['a'] = A.sparsity()
        S = conic('S','qpoases',qp)
        print(S)


    .. exec-block:: octave

        H = 2*DM.eye(2); [hidden]
        A = DM.ones(1,2); [hidden]
        g = DM.zeros(2); [hidden]
        lba = 10; [hidden]

        qp = struct;
        qp.h = H.sparsity();
        qp.a = A.sparsity();
        S = conic('S','qpoases',qp);
        disp(S)


The returned ``Function`` instance will have a *different* input/output signature compared to the high-level interface, one that includes the matrices :math:`H` and :math:`A`:

.. side-by-side::
    .. exec-block:: python

        H = 2*DM.eye(2) [hidden]
        A = DM.ones(1,2) [hidden]
        g = DM.zeros(2) [hidden]
        lba = 10. [hidden]

        qp = {} [hidden]
        qp['h'] = H.sparsity() [hidden]
        qp['a'] = A.sparsity() [hidden]
        S = conic('S','qpoases',qp) [hidden]

        r = S(h=H, g=g, \
              a=A, lba=lba)
        x_opt = r['x']
        print('x_opt: ', x_opt)

    &&

    .. exec-block:: octave

        H = 2*DM.eye(2); [hidden]
        A = DM.ones(1,2); [hidden]
        g = DM.zeros(2); [hidden]
        lba = 10; [hidden]

        qp = struct; [hidden]
        qp.h = H.sparsity(); [hidden]
        qp.a = A.sparsity(); [hidden]
        S = conic('S','qpoases',qp); [hidden]

        r = S('h', H, 'g', g,...
              'a', A, 'lba', lba);
        x_opt = r.x;
        disp(x_opt)



For-loop equivalents
--------------------

When modeling using expression graphs in |casadi|, it is a common pattern to use of for-loop constructs of the host language (C++/Python/Matlab).

The graph size will grow linearly with the loop size :math:`n`, and so will the construction time of the expression graph and the initialization time of functions using that expression.

We offer some special constructs that improve on this complexity.

Map
^^^

Suppose you are interested in computing a function :math:`f : \mathbb{R}^{n} \rightarrow \mathbb{R}^{m}` repeatedly on all columns of a matrix :math:`X \in \mathbb{R}^{n \times N}`, and aggregating all results in a result matrix :math:`Y \in \mathbb{R}^{m \times N}`:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym("x") [hidden]
        f = Function("f",[x],[sin(x)]) [hidden]

        N = 4
        X = MX.sym("X",1,N)

        print(f)

        ys = []
        for i in range(N):
          ys.append(f(X[:,i]))

        Y = hcat(ys)
        F = Function('F',[X],[Y])
        print(F)

    &&

    .. exec-block:: octave

        x = SX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]

        N = 4;
        X = MX.sym('X',1,N);

        disp(f)

        ys = {};
        for i=1:N
          ys{end+1} = f(X(:,i));
        end
        Y = [ys{:}];
        F = Function('F',{X},{Y});
        disp(F)

The aggregate function :math:`F : \mathbb{R}^{n \times N} \rightarrow \mathbb{R}^{m \times N}` can be obtained with the ``map`` construct:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym("x") [hidden]
        f = Function("f",[x],[sin(x)]) [hidden]
        N = 4 [hidden]


        F = f.map(N)
        print(F)

    &&

    .. exec-block:: octave

        x = SX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]
        N = 4; [hidden]

        F = f.map(N);
        disp(F)

|casadi| can be instructed to parallelize when :math:`F` gets evaluated. In the following example, we dedicate 2 threads for the ``map`` task.

.. side-by-side::
    .. exec-block:: python

        x = SX.sym("x") [hidden]
        f = Function("f",[x],[sin(x)]) [hidden]
        N = 4 [hidden]

        F = f.map(N,"thread",2)
        print(F)

    &&

    .. exec-block:: octave

        x = SX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]
        N = 4; [hidden]

        F = f.map(N,'thread',2);
        disp(F)

The ``map`` operation supports primitive functions :math:`f` with multiple inputs/outputs which can also be matrices. Aggregation will always happen horizontally.

The ``map`` operation exhibits constant graph size and initialization time.

Fold
^^^^

In case each for-loop iteration depends on the result from the previous iteration, the ``fold`` construct applies. In the following, the ``x`` variable acts as an accumulater that is initialized by :math:`x_0 \in \mathbb{R}^{n}`:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym("x") [hidden]
        f = Function("f",[x],[sin(x)]) [hidden]
        x0 = MX.sym("x0") [hidden]
        N = 4 [hidden]

        x = x0
        for i in range(N):
          x = f(x)

        F = Function('F',[x0],[x])
        print(F)

    &&

    .. exec-block:: octave

        x = SX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]
        x0 = MX.sym('x0'); [hidden]
        N = 4; [hidden]

        x = x0;
        for i=1:N
          x = f(x);
        end

        F = Function('F',{x0},{x});
        disp(F)

For a given function :math:`f : \mathbb{R}^{n} \rightarrow \mathbb{R}^{n}`, the result function :math:`F : \mathbb{R}^{n} \rightarrow \mathbb{R}^{n}` can be obtained with the ``fold`` construct:

.. side-by-side::
    .. exec-block:: python

        x = SX.sym("x") [hidden]
        f = Function("f",[x],[sin(x)]) [hidden]
        N = 4 [hidden]

        F = f.fold(N)
        print(F)
    &&

    .. exec-block:: octave

        x = SX.sym('x'); [hidden]
        f = Function('f',{x},{sin(x)}); [hidden]
        N = 4; [hidden]

        F = f.fold(N);
        disp(F)

In case intermediate accumulator values are desired as output (:math:`\mathbb{R}^{n} \rightarrow \mathbb{R}^{n \times N}`), use ``mapaccum`` instead of ``fold``.

The ``fold``/``mapaccum`` operation supports primitive functions :math:`f` with multiple inputs/outputs which can also be matrices.
The first input and output are used for accumulating, while the remainder inputs are read column-wise over the iterations.

The ``map``/``mapaccum`` operation exhibits a graph size and initialization time that scales logarithmically with :math:`n`.

.. rubric:: Footnotes

.. [#f1] for problems with free end time, you can always scale time by introducing an extra parameter and substitute :math:`t` for a dimensionless time variable that goes from 0 to 1


