.. include:: defs.rst

The |DaeBuilder| class
======================

.. _sec-daebuilder:

The |DaeBuilder| class in |casadi| is an auxiliary class intended to
facilitate the modeling as well as import and export of complex dynamical
systems. This class has lower level of abstraction than
a physical modeling language such as Modelica (cf. :numref:`sec-modelica`),
while still being higher level than working directly with |casadi| symbolic
expressions. In particular, it can be used to import or export physical models
defined using the Functional Mock-up Interface (FMI) format. It can also be used as a building block for constructing a domain specific modeling environments.

The syntax and scope of the class has changed significantly since its initial inception. In the following, we provide an overview of the syntax introduced in CasADi 3.7.

.. _sec-model_variables:

|DaeBuilder| model variables and their categorization
-----------------------------------------------------
Each |DaeBuilder| instance maintains a list of model variables. The variables must be uniquely named and hold import meta information about the variable such as its type (e.g. floating point, integer or sting-valued), minimum and maximum values, nominal values, dimension, description, unit etc. The format for model variables closely mimics that of the FMI standard, version 3.0, and to be able to use the |DaeBuilder| class efficently, you should probably be somewhat familar with the format. As in FMI, variables are classified according to their *causality* and *variability* (consult the FMI specification if you are unfamilar with these terms). They are also given a CasADi-specific categorization into the following categories:

:math:`t`
  The independent variable, usually time

:math:`c`
  Constants

:math:`p`
  Independent parameters

:math:`d`
  Dependent parameters, calculated from :math:`c` and :math:`p` and, acyclically, other :math:`d`

:math:`u`
  Input variables

:math:`x`
  Differential states, each defined by an explicit ordinary differential equation (ODE)

:math:`q`
  Quadrature states, i.e. a differential states that do not appear in
  the right-hand-sides allowing them to be calculated by so-called quadrature formulas

:math:`z`
  Algebraic variables, defined implicitly using residual equations

:math:`w`
  Dependent variables, calculated explicitly from :math:`c`, :math:`p`, :math:`d`, :math:`u`, :math:`x`, :math:`z` and, acyclically, other :math:`w`

:math:`y`
  Output variables, calculated explicitly from :math:`c`, :math:`p`, :math:`d`, :math:`u`, :math:`x`, :math:`q`, :math:`z` and :math:`w`

There are also additional categories of internal character, for certain derived expressions (residual variables, event indicators, derivative expressions etc.)

The category for a variable is determined from its causality and variability according to the following table, where the rows correspond to the causality and the columns to the variability (cf. Table 17 in FMI specification, 3.0.2):

+-------------+------------+-------------+-----------+--------------+-----------+--------------+
||            || parameter || calculated || input    || output      || local    || independent |
||            ||           || Parameter  ||          ||             ||          ||             |
+=============+============+=============+===========+==============+===========+==============+
| constant    | n/a        | n/a         | n/a       | :math:`c`    | :math:`c` | n/a          |
+-------------+------------+-------------+-----------+--------------+-----------+--------------+
| fixed       | :math:`c`  | :math:`d`   | n/a       | n/a          | :math:`d` | n/a          |
+-------------+------------+-------------+-----------+--------------+-----------+--------------+
| tunable     | :math:`p`  | :math:`d`   | n/a       | n/a          | :math:`d` | n/a          |
+-------------+------------+-------------+-----------+--------------+-----------+--------------+
| discrete    | n/a        | n/a         | :math:`u` | :math:`y`    | :math:`x` | n/a          |
+-------------+------------+-------------+-----------+--------------+-----------+--------------+
| continuous  | n/a        | n/a         | :math:`u` | :math:`u` or | :math:`x` | :math:`t`    |
|             |            |             |           | :math:`q`    |           |              |
+-------------+------------+-------------+-----------+--------------+-----------+--------------+

Not all combinations of causality and variability are permitted, as explained in the FMI specification. These combinations are marked "n/a" in the table. There may also be multiple combinations mapping to the same variable category. For example, discrete variables are treated as continuous variables with time derivative zero in CasADi. For variables with output causality and continuous variability, the category is either :math:`y` or :math:`q` depending on whether there exists a defining equation for the variable itself or its time derivative.

.. _sec-model_equations:

|DaeBuilder| model equations
----------------------------

At the time of this writing, |DaeBuilder| instances supported the following types of equations, leaving out :math:`c`, :math:`d` and :math:`w` for simplicity:

  * Ordinary differential equation, e.g.: :math:`\dot{x} = f_{\text{ode}}(t,x,z,u,p)`
  * Output equation, e.g.: :math:`y = f_{\text{ydef}}(t,x,z,u,p)`
  * Algebraic equation, e.g.: :math:`0 = f_{\text{alg}}(t,x,z,u,p)`
  * Quadrature equation, e.g.: :math:`\dot{q} = f_{\text{quad}}(t,x,z,u,p)`
  * When equation, e.g.: when :math:`f_\text{zero}(t,x,z,u,p) < 0` then :math:`x := f_\text{transition}(t,x,z,u,p)`

Future versions of |DaeBuilder| may support more types of equations, in particular related to proper handling of initial equations and event handling, guided by the Modelica standard.

The functions in the list above are associated with additional, internal, variable categories, in addition to the categories listed in :numref:`sec-model_variables`. For example ordinary differential equations and quadrature equations are associated with derivative variables and algebraic equations are associated with residual variables. When equations are associated with zero-crossing conditions and a set of assignment operations. We will return to when equations in :numref:`sec-hybrid` below.

Model equations may be given either as symbolic expressions, as described in :numref:`sec-daebuilder_symbolic`, from Modelica source code as described in :numref:`sec-modelica` or via standard Functional Mock-up Unit (FMU) shared libraries, as described in :numref:`sec-fmi`.

.. _sec-daebuilder_symbolic:

Constructing a |DaeBuilder| instance symbolically
-------------------------------------------------

We can create |DaeBuilder| instances using CasADi symbolic expressions. The current ambition is to support an acausal model approach taking elements from the Functional Mock-up Interface (FMI) standard as well as the proposed Base Modelica standard. Specifically, model variables (cf. :numref:`sec-model_variables`) are defined using the conventions defined by the FMI standard, version 3.0. The model equations (cf. :numref:`sec-model_equations`), on the other hand, where standard FMI relies on a black-box C API, we try to follow the proposed Base Modelica standard. We will discuss in :numref:`sec-modelica` how this support can be used to enable symbolic import of models formulated in actual Modelica.

To illustrate the usage, consider the following simple DAE corresponding to a controlled rocket subject to quadratic air friction term and gravity, which loses mass as it uses up fuel:

.. math::

    \begin{aligned}
     \dot{h} &= v,                    \qquad &h(0) = 0 \\
     \dot{v} &= (u - a \, v^2)/m - g, \qquad &v(0) = 0 \\
     \dot{m} &= -b \, u^2,            \qquad &m(0) = 1
    \end{aligned}

where the three states correspond to height, velocity and mass, respectively.
:math:`u` is the thrust of the rocket and :math:`(a,b)` are parameters.

To construct a DAE formulation for this problem, start with an empty |DaeBuilder| instance and add the model variables, including any meta information, and model equations step-by-step as follows.

.. side-by-side::
    .. exec-block:: python

        dae = DaeBuilder('rocket')
        # Add model variables
        a = dae.add('a', 'parameter', 'tunable')
        b = dae.add('b', 'parameter', 'tunable')
        u = dae.add('u', 'input')
        h = dae.add('h', dict(start=0))
        v = dae.add('v', dict(start=0))
        m = dae.add('m', dict(start=1))
        # Change some meta information
        dae.set_unit('h', 'm')
        dae.set_unit('v', 'm/s')
        dae.set_unit('m', 'kg')
        # Add model equations
        dae.eq(dae.der(h), v)
        dae.eq(dae.der(v), (u-a*v**2)/m-9.81)
        dae.eq(dae.der(m), -b*u**2)

    &&

    .. exec-block:: octave

        dae = DaeBuilder('rocket')
        % Add model variables
        a = dae.add('a', 'parameter', 'tunable');
        b = dae.add('b', 'parameter', 'tunable');
        u = dae.add('u', 'input');
        h = dae.add('h', struct('start', 0));
        v = dae.add('v', struct('start', 0));
        m = dae.add('m', struct('start', 1));
        % Change some meta information
        dae.set_unit('h','m');
        dae.set_unit('v','m/s');
        dae.set_unit('m','kg');
        % Add model equations
        dae.eq(dae.der(h), v);
        dae.eq(dae.der(v), (u-a*v^2)/m-9.81);
        dae.eq(dae.der(m), -b*u^2);

Other input and output expressions can be added in an analogous way. For a full
list of functions, see the C++ API documentation for |DaeBuilder|.

.. _sec-modelica:

Creating a |DaeBuilder| instance from Modelica, symbolically
------------------------------------------------------------

The original goal of the |DaeBuilder| class was to enable symbolic import of models formulated in Modelica. This was made possible by supporting the import models in the "FMUX" format, an XML format derived from FMI 1.0. In this format, the model variables (cf. :numref:`sec-model_variables`) are defined using standard FMI 1.0 XML syntax, while the model equations (c.f. :numref:`sec-model_equations`) are provided in the same XML instead of the C API defined by the standard. The support for FMUX import in CasADi never reached a mature state, in part because there was never a mature way to generate  FMUX files using available Modelica tools. While originally supported by the JModelica.org tool, this support was later dropped in favor of a tighter, SWIG-based CasADi interface. FMUX export has also been implemented in the OpenModelica compiler.

There is still legacy support for FMUX import in the |DaeBuilder| class, but since it is not actively maintained and very incomplete, it may require diving into the C++ source code to get it to work for nontrivial models.

We also note that JModelica.org's closed-source successor code, OCT, does retain CasADi interoperability. Details of how to use OCT to generate CasADi expressions is best described in OCT's user guide. As of this writing, this support did not use the DaeBuilder class in CasADi, although it is possible to create DaeBuilder instances from OCT-generated CasADi expressions, using the standard symbolic API described in :numref:`sec-daebuilder_symbolic`.

In future versions of CasADi, we hope to provide a more mature support for symbolic import of Modelica models based on the emerging Base Modelica standard. Base Modelica is intended to be a subset of the full Modelica language, which can be generated from any Modelica models, but removing many of the high-level features associated to object oriented physical modeling.

Unlike previous attempts at symbolic import, this support will also include support for hybrid (event-driven) modeling, while still retaining the differentiability needed by gradient-based optimization algorithms. Specifically, we try to ensure that the zero-crossing functions used to trigger events are differentiable and, secondly, that the state transitions are provided as differentiable functions. The Base Modelica effort currently does not include a Modelica parser, instead we are working with the open-source Rumoca translator, which can be used to convert Base Modelica source code into Pythons scripts that construct DaeBuilder instances symbolically.

.. _sec-fmi:

Constructing a |DaeBuilder| instance from an FMU
------------------------------------------------
The DaeBuilder class can be used to import standard FMUs, adhering to the FMI standard version 2.0 or 3.0. This is a dedicated, native, foreign function interface (FFI) of CasADi that communicates directly with the FMI C API. The evaluation takes place in CasADi function objects created using the function factory described in :numref:`sec-daebuilder_factory`. Importantly, the binary FMI interface enables the efficient calculation of first and second derivatives of the model equations using a hybrid symbolic-numeric approach, relying on analytical derivatives for the first order derivatives and efficient finite differences implementations for the second order derivatives. Sparsity and parallelization is also exploited, whenever possible. For an overview of the functionality, we refer to the fmu_demo.ipynb tutorial in the CasADi repository.

.. _sec-hybrid:

Hybrid modeling with |DaeBuilder|
---------------------------------

CasADi 3.7 introduces initial support for hybrid modeling and sensitivity analysis. In |DaeBuilder|, state events can be formulated with when equations, mirroring the syntax in (Base) Modelica. Whenever a when equation is added, a differential zero crossing condition is created automatically, using a form compatible with an extension of the CasADi integrator class to enable automatic sensitivity analysis of hybrid systems. While the event support is still at an early stage, it has been successfully used to model simple hybrid systems, such as the following bouncing ball example:

.. side-by-side::
    .. exec-block:: python

        dae = DaeBuilder('bouncing_ball')
        # Add model variables
        h = dae.add('h', dict(start=5))
        v = dae.add('v', dict(start=0))
        # Add model equations
        dae.eq(dae.der(h), v)
        dae.eq(dae.der(v), -9.81)
        dae.when(h < 0, \
          [dae.reinit('v', -0.8*dae.pre(v))])

    &&

    .. exec-block:: octave

        dae = DaeBuilder('bouncing_ball')
        % Add model variables
        h = dae.add('h', struct('start', 5));
        v = dae.add('v', struct('start', 0));
        % Add model equations
        dae.eq(dae.der(h), v);
        dae.eq(dae.der(v), -9.81);
        dae.when(h < 0, ...
          {dae.reinit('v', -0.8*dae.pre(v))});


.. _sec-daebuilder_factory:

Function factory
----------------
Once a |DaeBuilder| has been formulated and possibly reformulated to
a satisfactory form, we can generate |casadi| functions corresponding to the
input and output expressions outlined in :numref:`sec-model_variables` and :numref:`sec-model_equations`.
For example, to create a function for the ODE right-hand-side for the rocket
model in :numref:`sec-daebuilder_symbolic`, simply provide a display
name of the function being created, a list of input expressions
and a list of output expressions:

.. side-by-side::
    .. code-block:: python

        f = dae.create('f',\
             ['x','u','p'],['ode'])

    &&

    .. code-block:: octave

        f = dae.create('f',...
             {'x','u','p'},{'ode'});


Using a naming convention, we can also create Jacobians, e.g. for the 'ode'
output with respect to 'x':

.. side-by-side::
    .. code-block:: python

        f = dae.create('f',\
             ['x','u','p'],\
             ['jac_ode_x'])

    &&

    .. code-block:: octave

        f = dae.create('f',...
             {'x','u','p'},
             {'jac_ode_x'});

