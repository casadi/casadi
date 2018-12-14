.. _sec-opti:

Opti stack
==========

.. include:: defs.rst

The Opti stack is a collection of |casadi| helper classes that provides a close correspondence between mathematical NLP notation, e.g.

.. math::
  :label: simple_nlp

    \begin{array}{cc}
    \begin{array}{c}
    \text{minimize} \\
    x,y
    \end{array}
    &
    (y-x^2)^2  \\
    \begin{array}{c}
    \text{subject to}
    \end{array}
    & x^2+y^2=1 \\
    & x+y \ge 1,
    \end{array}


and computer code:

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti()

        x = opti.variable()
        y = opti.variable()

        opti.minimize(  (y-x**2)**2   )
        opti.subject_to( x**2+y**2==1 )
        opti.subject_to(       x+y>=1 )

        opti.solver('ipopt')
        opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]


        sol = opti.solve()

        print(sol.value(x))
        print(sol.value(y))
    &&

    .. exec-block:: octave

        opti = casadi.Opti();

        x = opti.variable();
        y = opti.variable();

        opti.minimize(  (y-x^2)^2   );
        opti.subject_to( x^2+y^2==1 );
        opti.subject_to(     x+y>=1 );

        opti.solver('ipopt');
        options = struct('print_time',false); [hidden]
        ipopt_options = struct('print_level',1); [hidden]

        opti.solver('ipopt',options,ipopt_options); [hidden]
        sol = opti.solve();

        sol.value(x)
        sol.value(y)

The main characteristics of the Opti stack are:


* Allows *natural* syntax for constraints.
* Indexing/bookkeeping of decision variables is hidden.
* Closer mapping of numerical data-type to the host language: no encounter with |DM|.

Problem specification
---------------------

**Variables**: Declare any amount of decision variables:

* ``x = opti.variable()``: scalar
* ``x = opti.variable(5)``: column vector
* ``x = opti.variable(5,3)``: matrix
* ``x = opti.variable(5,5,'symmetric')``: symmetric matrix


The order in which you declare the variables is respected by the solver.
Note that the variables are in fact plain MX symbols.
You may perform any CasADi MX operations on them, e.g. embedding integrator calls.

**Parameters**: Declare any amount of parameters. You must fix them to a specific numerical value before solving, and you may overwrite this value at any time.

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]

        p = opti.parameter()
        opti.set_value(p, 3)
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]

        p = opti.parameter();
        opti.set_value(p, 3)


**Objective**: Declare an objective using an expression that may involve all variables or parameters. Calling the command again discards the old objective.

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        x = opti.variable() [hidden]
        y = opti.variable() [hidden]
        p = opti.parameter() [hidden]

        opti.minimize(   sin(x*(y-p))   )
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        x = opti.variable(); [hidden]
        y = opti.variable(); [hidden]
        p = opti.parameter(); [hidden]

        opti.minimize(   sin(x*(y-p))   )

**Constraints**: Declare any amount of equality/inequality constraints:

* ``opti.subject_to( sqrt(x+y) >= 1)``: inequality
* ``opti.subject_to( sqrt(x+y) > 1)``: same as above
* ``opti.subject_to( 1<= sqrt(x+y) )``: same as above
* ``opti.subject_to( 5*x+y==1 )``: equality

You may also throw in several constraints at once:

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        x = opti.variable() [hidden]
        y = opti.variable() [hidden]

        opti.subject_to([x*y>=1,x==3])
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        x = opti.variable(); [hidden]
        y = opti.variable(); [hidden]

        opti.subject_to({x*y>=1,x==3});

You may declare double inequalities:

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        x = opti.variable() [hidden]

        opti.subject_to( opti.bounded(0,x,1) )
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        x = opti.variable(); [hidden]

        opti.subject_to( 0<=x<=1 );

When the bounds of the double inequalities are free of variables, the constraint will be passed on efficiently to solvers that support them (notably IPOPT).

You may make element-wise (in)equalities with vectors:

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        p = opti.parameter() [hidden]
        x = opti.variable(5,1)

        opti.subject_to( x*p<=3 )
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        p = opti.parameter(); [hidden]
        x = opti.variable(5,1);

        opti.subject_to( x*p<=3 )

Elementwise (in)equalities for matrices are not supported with a natural syntax,
since there is an ambiguity with semi-definiteness constraints.
The workaround is to vectorize first:

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        A = opti.variable(5,5)
        opti.subject_to( vec(A)<=3 )
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        A = opti.variable(5,5);
        opti.subject_to( A(:)<=3 );

Each ``subject_to`` command adds to the set of constraints in the problem specification.
Use ``subject_to()`` to empty this set and start over.

**Solver**:
You must always declare the ``solver`` (numerical back-end).
An optional dictionary of CasADi plugin options can be given as second argument.
An optional dictionary of ``solver`` options can be given as third argument.

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]

        x = opti.variable() [hidden]
        y = opti.variable() [hidden]

        opti.minimize(  (y-x**2)**2   ) [hidden]
        opti.subject_to( x**2+y**2==1 ) [hidden]
        opti.subject_to(       x+y>=1 ) [hidden]

        p_opts = {"expand":True}
        s_opts = {"max_iter": 100}
        opti.solver("ipopt",p_opts,
                            s_opts)
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]

        x = opti.variable(); [hidden]
        y = opti.variable(); [hidden]

        opti.minimize(  (y-x^2)^2   ); [hidden]
        opti.subject_to( x^2+y^2==1 ); [hidden]
        opti.subject_to(     x+y>=1 ); [hidden]

        p_opts = struct('expand',true);
        s_opts = struct('max_iter',100);
        opti.solver('ipopt',p_opts,
                            s_opts);

**Initial guess**:
You may provide initial guesses for decision variables (or simple mappings of decision variables). When no initial guess is provided, numerical zero is assumed.

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]
        x = opti.variable() [hidden]

        opti.set_initial(x, 2)
        opti.set_initial(10*x[0], 2)
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]
        x = opti.variable(); [hidden]

        opti.set_initial(x, 2);
        opti.set_initial(10*x(1), 2)

Problem solving and retrieving
------------------------------

Solving
^^^^^^^

After setting up the problem, you may call the solve method, which constructs a CasADi ``nlpsol`` and calls it.

.. side-by-side::
    .. exec-block:: python

        opti = casadi.Opti() [hidden]

        x = opti.variable() [hidden]
        y = opti.variable() [hidden]

        opti.minimize(  (y-x**2)**2   ) [hidden]
        opti.subject_to( x**2+y**2==1 ) [hidden]
        opti.subject_to(       x+y>=1 ) [hidden]

        opti.solver("ipopt") [hidden]

        sol = opti.solve()
    &&

    .. exec-block:: octave

        opti = casadi.Opti(); [hidden]

        x = opti.variable(); [hidden]
        y = opti.variable(); [hidden]

        opti.minimize(  (y-x^2)^2   ); [hidden]
        opti.subject_to( x^2+y^2==1 ); [hidden]
        opti.subject_to(     x+y>=1 ); [hidden]

        opti.solver('ipopt'); [hidden]

        sol = opti.solve();

The call will fail with an error if the solver fails to convergence. You may still inspect the non-converged solution (see :numref:`sec-opti-extra`).

You may call ``solve`` any number of times. You will always get an immutable copy of the problem specification and its solution.
Consecutively calling ``solve`` will not help the convergence of the problem.

To warm start a solver, you need to explicitly transfer the solution of one problem to the initial value of the next.


.. exec-block:: python

    opti = casadi.Opti() [hidden]

    x = opti.variable() [hidden]
    y = opti.variable() [hidden]

    opti.minimize(  (y-x**2)**2   ) [hidden]
    opti.subject_to( x**2+y**2==1 ) [hidden]
    opti.subject_to(       x+y>=1 ) [hidden]

    opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]


    sol1 = opti.solve()
    print(sol1.stats()["iter_count"])
    
    # Solving again makes no difference
    sol1 = opti.solve()
    print(sol1.stats()["iter_count"])

    # Passing initial makes a difference
    opti.set_initial(sol1.value_variables())
    sol2 = opti.solve()
    print(sol2.stats()["iter_count"])


.. exec-block:: octave

    opti = casadi.Opti(); [hidden]

    x = opti.variable(); [hidden]
    y = opti.variable(); [hidden]

    opti.minimize(  (y-x^2)^2   ); [hidden]
    opti.subject_to( x^2+y^2==1 ); [hidden]
    opti.subject_to(     x+y>=1 ); [hidden]

    options = struct('print_time',false); [hidden]
    ipopt_options = struct('print_level',1); [hidden]

    opti.solver('ipopt',options,ipopt_options); [hidden]

    sol1 = opti.solve();
    sol1.stats.iter_count

    % Solving again makes no difference
    sol1 = opti.solve();
    sol1.stats.iter_count

    % Passing initial makes a difference
    opti.set_initial(sol1.value_variables());
    sol2 = opti.solve();
    sol2.stats.iter_count

In order to initialize the dual variables, e.g. when solving a set of similar optimization problems, you can use the following syntax:


.. exec-block:: python

    opti = casadi.Opti() [hidden]

    x = opti.variable() [hidden]
    y = opti.variable() [hidden]

    opti.minimize(  (y-x**2)**2   ) [hidden]
    opti.subject_to( x**2+y**2==1 ) [hidden]
    opti.subject_to(       x+y>=1 ) [hidden]

    opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]

    sol = opti.solve()
    lam_g0 = sol.value(opti.lam_g)
    opti.set_initial(opti.lam_g, lam_g0)

.. exec-block:: octave

    opti = casadi.Opti(); [hidden]

    x = opti.variable(); [hidden]
    y = opti.variable(); [hidden]

    opti.minimize(  (y-x^2)^2   ); [hidden]
    opti.subject_to( x^2+y^2==1 ); [hidden]
    opti.subject_to(     x+y>=1 ); [hidden]

    options = struct('print_time',false); [hidden]
    ipopt_options = struct('print_level',1); [hidden]

    opti.solver('ipopt',options,ipopt_options); [hidden]

    sol = opti.solve();
    lam_g0 = sol.value(opti.lam_g);
    opti.set_initial(opti.lam_g, lam_g0);

Numerical value at the solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Afterwards, you may retrieve the numerical values of variables (or expressions of those variables) at the solution:


* ``sol.value(x)``: value of a decision variable
* ``sol.value(p)``: value of a parameter
* ``sol.value(sin(x+p))``: value of an expression
* ``sol.value(jacobian(opti.g,opti.x))``: value of constraint jacobian

Note that the return type of ``value`` is sparse when applicable.

Numerical value at other points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may pass a list of overruling assignment expressions to ``value``.
In the following code, we are asking for the value of the objective,
using all optimal values at the solution, except for ``y``, which we set equal to 2.
Note that such statement does not modify the actual optimal value of ``y`` in a permanent way.


.. exec-block:: python

    opti = casadi.Opti() [hidden]

    x = opti.variable() [hidden]
    y = opti.variable() [hidden]

    opti.subject_to( x**2+y**2==1 ) [hidden]
    opti.subject_to(       x+y>=1 ) [hidden]

    obj = (y-x**2)**2

    opti.minimize(obj)

    opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]
    sol = opti.solve() [hidden]

    print(sol.value(obj,[y==2]))

.. exec-block:: octave

    opti = casadi.Opti(); [hidden]

    x = opti.variable(); [hidden]
    y = opti.variable(); [hidden]

    opti.subject_to( x^2+y^2==1 ); [hidden]
    opti.subject_to(     x+y>=1 ); [hidden]

    obj = (y-x^2)^2;

    opti.minimize(obj);

    options = struct('print_time',false); [hidden]
    ipopt_options = struct('print_level',1); [hidden]

    opti.solver('ipopt',options,ipopt_options); [hidden]
    sol = opti.solve(); [hidden]

    sol.value(obj,{y==2})

A related usage pattern is to evaluate an expression at the initial guess:


.. exec-block:: python

    opti = casadi.Opti() [hidden]

    x = opti.variable() [hidden]
    y = opti.variable() [hidden]

    opti.minimize(  (y-x**2)**2   ) [hidden]
    opti.subject_to( x**2+y**2==1 ) [hidden]
    opti.subject_to(       x+y>=1 ) [hidden]

    opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]
    sol = opti.solve() [hidden]

    print(sol.value(x**2+y,opti.initial()))

.. exec-block:: octave

    opti = casadi.Opti(); [hidden]

    x = opti.variable(); [hidden]
    y = opti.variable(); [hidden]

    opti.minimize(  (y-x^2)^2   ); [hidden]
    opti.subject_to( x^2+y^2==1 ); [hidden]
    opti.subject_to(     x+y>=1 ); [hidden]

    options = struct('print_time',false); [hidden]
    ipopt_options = struct('print_level',1); [hidden]

    opti.solver('ipopt',options,ipopt_options); [hidden]
    sol = opti.solve(); [hidden]

    sol.value(x**2+y,opti.initial())

Dual variables
^^^^^^^^^^^^^^

In order to obtain dual variables (Lagrange multipliers) of constraints, make sure you save the constraint expression first:


.. exec-block:: python
    :hide-output:

    opti = casadi.Opti() [hidden]

    x = opti.variable() [hidden]
    y = opti.variable() [hidden]

    opti.minimize(  (y-x**2)**2   ) [hidden]

    opti.solver("ipopt",dict(print_time=False),dict(print_level=False)) [hidden]

    con = sin(x+y)>=1
    opti.subject_to(con)
    sol = opti.solve()

    print(sol.value(opti.dual(con)))

.. exec-block:: octave
    :hide-output:

    opti = casadi.Opti(); [hidden]

    x = opti.variable(); [hidden]
    y = opti.variable(); [hidden]

    opti.minimize(  (y-x^2)^2   ); [hidden]

    options = struct('print_time',false); [hidden]
    ipopt_options = struct('print_level',1); [hidden]

    opti.solver('ipopt',options,ipopt_options); [hidden]
    sol = opti.solve(); [hidden]

    con = sin(x+y)>=1;
    opti.subject_to(con);
    sol = opti.solve();

    sol.value(opti.dual(con))


.. _sec-opti-extra:

Extras
------


It may well happen that the solver does not find an optimal solution.
In such cases, you may still access the non-converged solution through debug mode:

.. code-block:: python

     opti.debug.value(x)

Related, you may inspect the value of an expression, at the initial guess that you supplied to the solver:

.. code-block:: python

      opti.debug.value(x,opti.initial())


In case the solver stops due to problem infeasibility, you may identify the problematic constraints with:

.. code-block:: python

      opti.debug.show_infeasibilities()

In case the solver reports NaN/Inf at a certain location, you may find out which constraint or variable is to blame by looking at its description:

.. code-block:: python

    opti.debug.x_describe(index)
    opti.debug.g_describe(index)


You may specify a callback function; it will be called at each iteration of the solver, with the current iteration number as argument.
To plot the progress of the solver, you may access the non-converged solution through debug mode:

.. code-block:: python

    opti.callback(lambda i: plot(opti.debug.value(x)))

.. code-block:: octave

    opti.callback(@(i) plot(opti.debug.value(x)))

The callback may be cleared from the Opti stack by calling the |Callback| function without arguments.


