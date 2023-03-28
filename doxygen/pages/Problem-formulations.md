# Problem formulations {#page_problem_formulations}

[TOC]

## General NLP formulation {#problem_formulations_general}

Most alpaqa solvers deal with problems in the following form:

@f[
\begin{equation}\tag{P}\label{eq:problem_main}
    \begin{aligned}
        & \underset{x}{\text{minimize}}
        & & f(x) &&&& f : \Rn \rightarrow \R \\
        & \text{subject to}
        & & \underline{x} \le x \le \overline{x} \\
        &&& \underline{z} \le g(x) \le \overline{z} &&&& g : \Rn \rightarrow \Rm
    \end{aligned}
\end{equation}
@f]

The solver needs to be able to evaluate the following required functions and
derivatives:
  - @ref alpaqa::TypeErasedProblem::eval_f "eval_f"                     @f$ : \Rn \to \R : x \mapsto f(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_f "eval_grad_f"           @f$ : \Rn \to \Rn : x \mapsto \nabla f(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_g "eval_g"                     @f$ : \Rn \to \Rm : x \mapsto g(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_g_prod "eval_grad_g_prod" @f$ : \Rn \times \Rm \to \Rn : (x, y) \mapsto \nabla g(x)\, y @f$

Usually, [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)
(AD) is used to evaluate the gradients and gradient-vector products. Many AD
software packages are available, see e.g. <https://autodiff.org/> for an overview.

Additionally, the solver needs to be able to project onto the rectangular sets
@f[
\begin{equation}
    \begin{aligned}
        C &\;\defeq\; \defset{x \in \Rn}{\underline{x} \le x \le \overline{x}}, \\
        D &\;\defeq\; \defset{z \in \Rm}{\underline{z} \le z \le \overline{z}}.
    \end{aligned}
\end{equation}
@f]

<!-- Given two boxes @f$ C @f$ and @f$ D @f$, the @ref alpaqa::BoxConstrProblem "BoxConstrProblem"
class provides default implementations for the necessary projections:
@ref alpaqa::TypeErasedProblem::eval_proj_diff_g "eval_proj_diff_g",
@ref alpaqa::TypeErasedProblem::eval_proj_multipliers "eval_proj_multipliers" and
@ref alpaqa::TypeErasedProblem::eval_prox_grad_step "eval_prox_grad_step". -->

## Problem API

The alpaqa solvers access the problem functions through the API outlined in
@ref alpaqa::TypeErasedProblem.  
Usually, problems are defined using C++ structs, providing the evaluations
described above as public member functions. These problem structs are
[structurally typed](https://en.wikipedia.org/wiki/Structural_type_system),
which means that they only need to provide member functions with the correct
names and signatures. Inheriting from a common base class is not required.

As an example, the following struct defines a problem that can be passed to the
alpaqa solvers. Detailed descriptions of each function can be found in the
@ref alpaqa::TypeErasedProblem documentation.

```cpp
struct RosenbrockProblem {
    USING_ALPAQA_CONFIG(alpaqa::DefaultConfig);
    // Problem dimensions
    length_t get_n() const; // number of unknowns
    length_t get_m() const; // number of general constraints

    // Cost
    real_t eval_f(crvec x) const;
    // Gradient of cost
    void eval_grad_f(crvec x, rvec grad_fx) const;

    // Constraints
    void eval_g(crvec x, rvec gx) const;
    // Gradient-vector product of constraints
    void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const;

    // Proximal gradient step
    real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad, rvec x̂, rvec p) const;
    // Projecting difference onto constraint set D
    void eval_proj_diff_g(crvec z, rvec p) const;
    // Projection of Lagrange multipliers
    void eval_proj_multipliers(rvec y, real_t max_y) const;
};
```

### Base classes for common use cases

Convenience classes with default implementations of some of these functions are
provided for common use cases:
  - @ref alpaqa::BoxConstrProblem "BoxConstrProblem" defines the 
    @ref alpaqa::TypeErasedProblem::eval_proj_diff_g "eval_proj_diff_g",
    @ref alpaqa::TypeErasedProblem::eval_proj_multipliers "eval_proj_multipliers" and
    @ref alpaqa::TypeErasedProblem::eval_prox_grad_step "eval_prox_grad_step" functions
    for the specific case where @f$ C @f$ and @f$ D @f$ are rectangular boxes,
    as in @f$ \eqref{eq:problem_main} @f$.
  - @ref alpaqa::UnconstrProblem "UnconstrProblem" defines those same functions,
    and additional empty implementations of constraint-related functions to
    allow a more concise formulation of unconstrained problems.

The user can simply inherit from these classes to inject the default
implementations into their problem definition, as demonstrated in the following
examples.
  - @ref C++/CustomCppProblem/main.cpp
  - @ref C++/SimpleUnconstrProblem/main.cpp

It is highly recommended to study the @ref C++/CustomCppProblem/main.cpp example
now to see how optimization problems can be formulated in practice, before we
continue with some more specialized use cases.

### Second-order derivatives

Some solvers can exploit information about the Hessian of the (augmented)
Lagrangian of the problem. To use these solvers, some of the following functions
are required, they should be added as member functions to your problem struct.
  - @ref alpaqa::TypeErasedProblem::eval_jac_g "eval_jac_g"
  - @ref alpaqa::TypeErasedProblem::get_jac_g_num_nonzeros "get_jac_g_num_nonzeros"
  - @ref alpaqa::TypeErasedProblem::eval_grad_gi "eval_grad_gi"
  - @ref alpaqa::TypeErasedProblem::eval_hess_L_prod "eval_hess_L_prod"
  - @ref alpaqa::TypeErasedProblem::eval_hess_L "eval_hess_L"
  - @ref alpaqa::TypeErasedProblem::get_hess_L_num_nonzeros "get_hess_L_num_nonzeros"
  - @ref alpaqa::TypeErasedProblem::eval_hess_ψ_prod "eval_hess_ψ_prod"
  - @ref alpaqa::TypeErasedProblem::eval_hess_ψ "eval_hess_ψ"
  - @ref alpaqa::TypeErasedProblem::get_hess_ψ_num_nonzeros "get_hess_ψ_num_nonzeros"

Sparse matrices are stored in [compressed column storage](https://www.eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseIntro)
(CCS) format. For symmetric Hessian matrices, only the upper triangle is stored.
Upon initialization, the number of nonzeros is queried by the solver using the
`get_xyz_num_nonzeros()` function, and storage is allocated for the arrays of
column/row indices and for the nonzero values. Then the `eval_xyz()` function is
called once with an empty `values` argument (`values.size() == 0`), indicating
that the column/row indices representing the sparsity should be initialized.
Subsequent calls to `eval_xyz()` then pass a non-empty `values` argument, in
addition to the initialized column/row indices, and the user should then
overwrite the nonzero values of the matrix.

If the matrix is dense, `get_xyz_num_nonzeros()` should return zero, the
column/row indices are not used, and the `values` argument to `eval_xyz()`
provides storage for a dense column-major matrix.

@note   Currently, symmetric dense matrices should store the full matrix, not
just the upper triangular part. This is different from the sparse case, and we
might want to change this in the future.

Some solvers do not require the full Hessian matrices, but use Hessian-vector
products only, for example when using Newton-CG. These products can often be
computed efficiently using automatic differentiation, at a computational cost
that's not much higher than a gradient evaluation.

The @ref alpaqa::TypeErasedProblem "TypeErasedProblem" class provides functions
to query which optional problem functions are available. For example,
@ref alpaqa::TypeErasedProblem::provides_eval_jac_g "provides_eval_jac_g"
returns true if the problem provides an implementation for
@ref alpaqa::TypeErasedProblem::eval_jac_g "eval_jac_g". Calling an optional
function that is not provided results in an @ref alpaqa::not_implemented_error
exception being thrown.

### Specialized combined evaluations

In practice, the solvers do not always evaluate the functions @f$ f(x) @f$ and
@f$ g(x) @f$ directly. Instead, they evaluate the Lagrangian and augmented
Lagrangian functions of the problem. In many applications, such as
single-shooting optimal control problems, some computations are common to the
evaluation of both @f$ f(x) @f$ and @f$ g(x) @f$, and significant speedups can
be achieved by providing implementations that evaluate both at the same time,
or even compute the (augmented) Lagrangian directly. Similarly, when using
automatic differentiation, evaluation of the gradient @f$ \nabla f(x) @f$
produces the function value @f$ f(x) @f$ as a byproduct, motivating the
simultaneous evaluation of these quantities as well.

The full list of these combined evaluations can be found in the @ref alpaqa::TypeErasedProblem "TypeErasedProblem"
documentation. They can be provided in the same fashion as `eval_f` above.
  - @ref alpaqa::TypeErasedProblem::eval_f_grad_f "eval_f_grad_f": @f$ f(x) @f$ and @f$ \nabla f(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_f_g "eval_f_g": @f$ f(x) @f$ and @f$ g(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_f_grad_g_prod "eval_grad_f_grad_g_prod": @f$ \nabla f(x) @f$ and @f$ \nabla g(x)\,y @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_L "eval_grad_L": gradient of Lagrangian: @f$ \nabla L(x, y) = \nabla f(x) + \nabla g(x)\,y @f$
  - @ref alpaqa::TypeErasedProblem::eval_ψ "eval_ψ": augmented Lagrangian: @f$ \psi(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_ψ_from_ŷ "eval_grad_ψ_from_ŷ": gradient of augmented Lagrangian: @f$ \nabla \psi(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_grad_ψ "eval_grad_ψ": gradient of augmented Lagrangian: @f$ \nabla \psi(x) @f$
  - @ref alpaqa::TypeErasedProblem::eval_ψ_grad_ψ "eval_ψ_grad_ψ": augmented Lagrangian and gradient: @f$ \psi(x) @f$ and @f$ \nabla \psi(x) @f$
