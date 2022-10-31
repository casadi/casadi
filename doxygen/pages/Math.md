# Math {#page_math}

[TOC]

## Augmented Lagrangian method

### Definitions

The nonlinear program we wish to solve is of the form:
@f[
\begin{equation}
    \begin{aligned}
        & \underset{x}{\text{minimize}}
        & & f(x) &&&& f : \Rn \rightarrow \R \\
        & \text{subject to}
        & & \underline{x} \le x \le \overline{x} \\
        &&& \underline{z} \le g(x) \le \overline{z} &&&& g : \Rn \rightarrow \R^m
    \end{aligned}
    \label{eq:problem-orig}
    \tag{P}
\end{equation}
@f]

Define the convex sets @f$C@f$ and @f$D@f$ as
@f[
\begin{equation}
    \begin{aligned}
        C &= \left\{ x \in \Rn \mid \underline x \le x \le \overline x \right\} \\
        D &= \left\{ z \in \R^m \mid \underline z \le z \le \overline z \right\}. \\
    \end{aligned}
    \label{eq:setCD}
\end{equation}
@f]
These rectangular boxes can be decomposed as Cartesian products of 1-dimensional closed intervals:
@f[
\begin{equation}
    \begin{aligned}
        C &= C_1 \times C_2 \times \dots \times C_n \quad&&\text{where } C_i = \left[ \underline x_i, \overline x_i \right] \\
        D &= D_1 \times D_2 \times \dots \times D_m \quad&&\text{where } D_i = \left[ \underline z_i, \overline z_i \right] \\
    \end{aligned}
    \label{eq:setCDcartprod}
\end{equation}
@f]
Using these definitions, problem @f$\eqref{eq:problem-orig}@f$ can equivalently be expressed as
@f[
\begin{equation}
    \begin{aligned}
        & \underset{x\in C}{\text{minimize}}
        & & f(x) \\
        & \text{subject to}
        & & g(x) \in D.
    \end{aligned}
    \label{eq:problem-in-setCD}
\end{equation}
@f]

After introduction of a slack variable @f$z@f$, problem @f$\eqref{eq:problem-in-setCD}@f$ can be stated as 

@f[
\begin{equation}
    \label{eq:problem-origCD-alm}
    \tag{P-ALM}
    \begin{aligned}
        & \underset{x\in C,\; z\in D}{\text{minimize}}
        & & f(x) \\
        & \text{subject to}
        & & g(x) - z = 0.
    \end{aligned}
\end{equation}
@f]

The Lagrangian function of problem @f$\eqref{eq:problem-origCD-alm}@f$ is given by:
@f[
\begin{equation}\label{eq:def-lagr}
    \begin{aligned}
        \Lagr : \Rn \times \R^m \times \R^m \rightarrow \R : (x, z, y) &\mapsto \Lagr(x, z, y) \\
        &\,\triangleq\, f(x) + \left\langle g(x) - z,\; y\right\rangle.
    \end{aligned}
\end{equation}
@f]
The vector @f$y \in \R^m@f$ is called the vector of Lagrange multipliers.

The augmented Lagrangian function with penalty factor @f$\Sigma@f$ of the problem @f$\eqref{eq:problem-origCD-alm}@f$ is defined as the sum of the Lagrangian function and a quadratic 
term that penalizes the constraint violation:
@f[
\begin{equation}\label{eq:def-auglagr}
    \begin{aligned}
        \Lagr_\Sigma : \Rn \times \R^m \times \R^m \rightarrow \R : (x, z, y) &\mapsto \Lagr_\Sigma(x, z, y) \\
        &\,\triangleq\, \Lagr(x, z, y) + \tfrac{1}{2} \left\|g(x) - z\right\|^2_\Sigma,
    \end{aligned}
\end{equation}
@f]
where @f$\Sigma@f$ is a symmetric positive definite @f$m\times m@f$ matrix that defines a norm on @f$\R^m@f$, @f$\|z\|^2_\Sigma \triangleq z^\top \Sigma z@f$.

### The augmented Lagrangian method algorithm

The augmented Lagrangian method for solving Problem @f$\eqref{eq:problem-origCD-alm}@f$
consists of the successive minimization of @f$\Lagr_\Sigma@f$ 
with respect to the decision variables @f$x@f$ and the slack variables @f$z@f$ (1),
after which the Lagrange multipliers @f$y@f$ are updated (2),
and the penalty factors @f$\Sigma_{ii}@f$ corresponding to constraints with high
violation are increased (3).

The augmented Lagrangian function is used as an exact penalty function 
for problem @f$\eqref{eq:problem-origCD-alm}@f$, it is equivalent to the shifted quadratic penalty method with shift @f$\Sigma^{-1}y@f$.

#### 1. Minimization of the augmented Lagrangian 

Using some algebraic manipulations, the augmented Lagrangian defined in @f$\eqref{eq:def-auglagr}@f$ can be expressed as
@f[
\begin{equation}
    \label{eq:auglagr2}
    \Lagr_\Sigma(x, z, y) = f(x) \;\;+\;\; \frac{1}{2} \Big\Vert g(x) - z + \Sigma^{-1} y \Big\Vert_\Sigma^2
    \;\;-\;\; \frac{1}{2}\Big\Vert y \Big\Vert_{\Sigma^{-1}.}^2
\end{equation}
@f]

At each iteration @f$\nu@f$ of the ALM algorithm, the following minimization problem is solved:
@f[
\begin{equation}
    \label{eq:alm-step-1-argmin}
    (x^\nu, z^\nu) \;\;=\;\; \underset{x\in C,\; z\in D}{\text{argmin}}\quad \Lagr_{\Sigma^{\nu-1}}(x,\, z;\, y^{\nu-1})
\end{equation}
@f]

#### 2. Update of the Lagrange multipliers

The update of the Lagrange multipliers corrects the shift @f$\Sigma^{-1} y@f$ in @f$\eqref{eq:auglagr2}@f$:
if the constraint violation @f$g(x^\nu) - z^\nu@f$ is positive, 
the shift is increased, in an attempt to drive the next iterate
towards a smaller constraint violation @f$g(x^\nu) - z^\nu@f$. The following update rule formalizes that idea:
@f[
\begin{equation}\label{eq:lagr-update-explanation}
    y^\nu \leftarrow y^{\nu-1} + \Sigma^{\nu-1} \left(g(x^\nu) - z^\nu\right)
\end{equation}
@f]
When the constraint violation becomes zero, the Lagrange multipliers are no longer updated.

As the penalty factors @f$\Sigma@f$ tend towards infinity, the shift @f$\Sigma^{-1} y@f$ has to vanish, because in that case,
the quadratic penalty method without shifts solves the problem exactly. 
For @f$\Sigma^{-1} y@f$ to vanish, the Lagrange multipliers must be bounded, which is achieved by the following 
projection:

Let @f$M > 0@f$ be some large but finite bound.
@f[
    \begin{gather}
        \underline y_i \triangleq \begin{cases}
            0 & \underline z_i = -\infty \\
            -M & \text{otherwise},
        \end{cases}
        \quad\quad\quad\quad\quad\quad 
        \overline y_i \triangleq \begin{cases}
            0 & \overline z_i = +\infty \\
            +M & \text{otherwise}
        \end{cases} \\[0.66em]
        Y \triangleq [\underline y_1, \overline y_1] \times \dots \times [\underline y_m, \overline y_m]
    \end{gather}
@f]
The result of @f$\eqref{eq:lagr-update-explanation}@f$ is therefore clamped as follows:
@f[
\begin{equation}\label{eq:lagr-update-explanation-clamp}
    y^\nu \leftarrow \Pi_Y\left(y^{\nu-1} + \Sigma^{\nu-1} \left(g(x^\nu) - z^\nu\right)\right)
\end{equation}
@f]

#### 3. Update of the penalty factors

When the penalty factor for the @f$i@f$-th constraint, @f$\Sigma_{ii}@f$ is increased, minimizing the violation of this 
constraint becomes more important in @f$\eqref{eq:alm-step-1-argmin}@f$. Therefore, if the constraint violation cannot be reduced by updating the 
shifts alone, the penalty factors are increased. 

Selecting when and by how much each penalty factor should be increased is more of a heuristic. The strategy 
used here is to compare the violation at the current iterate with the violation at the previous iterate, it is the same strategy as used in [QPALM](https://arxiv.org/abs/2010.02653). 
Denote the vector of constraint violations as @f$e^\nu \triangleq g(x^\nu) - z^\nu@f$.
Let @f$\theta \in (0, 1)@f$.  
If @f$|e^\nu_i| \le \theta |e^{\nu-1}_i|@f$, meaning that the constraint violation has decreased by at least a factor @f$\theta@f$ compared to the previous iteration, then the penalty factor is not updated.  
If the constraint violation did not decrease sufficiently, then the penalty factor @f$\Sigma_{ii}@f$ is increased by a factor
@f[
\begin{equation}\label{eq:multipliers-update-factor}
    \Delta \dfrac{|e^\nu_i|}{\|e^\nu\|_\infty},
\end{equation}
@f]
where @f$\Delta > 1@f$ is a tuning parameter. The violation of each individual constraint is scaled by
the maximum violation of all constraints,
such that the penalty factors of constraints with a large violation are increased 
more aggressively. If the factor in @f$\eqref{eq:multipliers-update-factor}@f$ is less than one, the penalty factor is not 
updated (otherwise it would result in a reduction of the penalty).

## PANOC

PANOC is an algorithm that solves optimization problems of the form:
@f[
\begin{equation}
    \begin{aligned}
        & \underset{x}{\text{minimize}}
        & & \psi(x) + h(x),
    \end{aligned}
    \label{eq:problem-panoc}
    \tag{P-PANOC}
\end{equation}
@f]
where @f$\psi : \Rn \rightarrow \R @f$ has Lipschitz gradient, and
@f$h : \Rn \rightarrow \overline \R @f$ allows efficient computation of the
proximal operator.

Recall the inner minimization problem @f$\eqref{eq:alm-step-1-argmin}@f$ 
in the first step of the ALM algorithm. It can be simplified to:
@f[
\begin{equation}
    \label{eq:problem-inner}
    \begin{aligned}
        \min_{x\in C,\, z\in D} \Lagr_\Sigma(x, z, y)
        &= -\tfrac{1}{2} \lVert y \rVert_{\Sigma^{-1}}^2
        + \min_{x\in C}\left\{ f(x)
        + \min_{z\in D} \left\{
        \tfrac{1}{2} \left\Vert z - \left(g(x) + \Sigma^{-1} y\right)\right\Vert_\Sigma^2
        \right\}
        \right\} \\
        &= -\tfrac{1}{2} \lVert y \rVert_{\Sigma^{-1}}^2
        + \min_{x\in C}\;\;
        \bigg\{
        \underbrace{
        f(x)
        + \tfrac{1}{2} \text{dist}_\Sigma^2 \left(
        g(x) + \Sigma^{-1}y, \ D
        \right)}_{\triangleq\, \psi_{\Sigma}(x;\,y)}
        \bigg\}
        \end{aligned}
\end{equation}
@f]
Within the PANOC algorithm,
the parameters @f$y@f$ and @f$\Sigma@f$ remain constant, and will be omitted from the function names to ease notation:
@f[
\begin{equation}
    \begin{aligned}
        \psi(x) &= f(x)
        + \tfrac{1}{2} \text{dist}_\Sigma^2 \left(
        g(x) + \Sigma^{-1}y, \ D
        \right)
    \end{aligned}
    \label{eq:psi-inner-panoc}
\end{equation}
@f]
The inner problem in @f$\eqref{eq:problem-inner}@f$ has the same minimizers as the following problem that will be solved using the PANOC algorithm:
@f[
\begin{equation}
    \begin{aligned}
        & \underset{x\in C}{\text{minimize}}
        & & \psi(x),
    \end{aligned}
    \label{eq:problem-inner-panoc}
\end{equation}
@f]
This problem is an instance of problem @f$\eqref{eq:problem-panoc}@f$ where the nonsmooth term @f$h@f$ is the indicator of the 
set @f$C@f$, @f$h(x) = \delta_C(x)@f$.

### Evaluation

The following is a list of symbols and formulas that are used in the 
implementation of the PANOC algorithm.

@f[
    \DeclareMathOperator{\prox}{\mathbf{prox}}
    \DeclareMathOperator*{\argmin}{\mathbf{argmin}}
    \DeclareMathOperator{\dist}{\mathbf{dist}}
    \DeclareMathOperator*{\minimize}{\mathbf{minimize}\;\;}

    \begin{aligned}
    y &\in \R^m &\text{Current Lagrange multipliers} \\
    \Sigma &\in \text{diag}(\R^m_{>0}) &\text{Current penalty factor} \\
    x^k &\in \Rn &\text{Current PANOC iterate} \\
    \gamma_k &\in \R_{>0} &\text{Current proximal gradient step size} \\[1em]
    \zeta^k &\triangleq g(x^k) + \Sigma^{-1}y &\text{Shifted constraint value}\\
    \hat{z}^k &\triangleq \Pi_D\left(g(x^k) + \Sigma^{-1}y\right) &\text{Closest feasible value for slack variable } z \\
    &= \Pi_D(\zeta^k) \\
    d^k &\triangleq \zeta^k - \Pi_D(\zeta^k) &\text{How far the shifted constraint value }\zeta \\
    &= \zeta^k - \hat{z}^k  &\text{is from the feasible set}\\
    e^k &\triangleq g(x^k) - \hat z^k &\text{Constraint violation} \\
    \hat{y}^k &\triangleq \Sigma\, d^k &\text{Candidate Lagrange multipliers,}\\
    &= \Sigma\, \left(g(x^k) + \Sigma^{-1}y - \Pi_D\left(g(x^k) + \Sigma^{-1}y\right)\right) &\text{see \eqref{eq:lagr-update-explanation}}\\
    &= y + \Sigma\,\left(g(x^k) - \hat z^k\right) \\
    &= y + \Sigma\, e^k \\[1em]
    \psi(x^k) &= \Lagr_\Sigma(x^k, \hat z^k, y) + \tfrac{1}{2} \lVert y \rVert_{\Sigma^{-1}}^2 &\text{PANOC objective function} \\    
    &= f(x^k) + \tfrac{1}{2} \dist_\Sigma^2\left(g(x^k) + \Sigma^{-1}y,\;D\right) \\
    &= f(x^k) + \tfrac{1}{2} \left\|\left(g(x^k) + \Sigma^{-1}y\right) -\Pi_D\left(g(x^k) + \Sigma^{-1}y\right)\right\|_\Sigma^2 \\
    &= f(x^k) + \tfrac{1}{2} \left\|\zeta^k - \hat{z}^k\right\|_\Sigma^2 \\
    &= f(x^k) + \tfrac{1}{2} \langle d^k, \hat{y}^k \rangle \\[1em]
    \nabla \psi(x^k) &= \nabla f(x^k)
        + \nabla g(x^k)\, \Sigma \left(g(x^k) + \Sigma^{-1}y - \Pi_D(g(x^k) + \Sigma^{-1}y)\right) &\text{Gradient of the objective} \\
    &= \nabla f(x^k) + \nabla g(x^k)\, \Sigma\, \big(\zeta^k - \hat{z}^k\big) \\
    &= \nabla f(x^k) + \nabla g(x^k)\, \hat{y}^k \\[1em]
    \hat{x}^k &\triangleq T_{\gamma^k}\left(x^k\right) &\text{Next proximal gradient iterate}\\
    &= \Pi_C\left(x^k - \gamma^k \nabla \psi(x^k)\right) \\
    p^k &\triangleq \hat{x}^k - x^k &\text{Proximal gradient step}\\
    r^k &\triangleq \tfrac{1}{\gamma^k} p^k &\text{Fixed-point residual (FPR)}\\[1em]
    \varphi_{\gamma^k}(x^k) &= \psi(x^k) + h(\hat{x}^k) + \tfrac{1}{2\gamma^k} \lVert \hat{x}^k - x^k \rVert^2 + \nabla\psi(x^k)^\top (\hat{x}^k - x^k) &\text{Forward-backward envelope (FBE)}\\
    &= \psi(x^k) + \tfrac{1}{2\gamma^k} \lVert p^k \rVert^2 + \nabla\psi(x^k)^\top p^k \\[1em]
    q^k &\triangleq H_k r^k &\text{Quasi-Newton step} \\
    x^{k+1} &= x^k + (1-\tau) p^k + \tau q^k &\text{Next PANOC iterate} \\
    \end{aligned}
@f]

Note that many of the intermediate values depend on the value of @f$x^k@f$, it
is sometimes easiest to define them as functions:
@f[
    \begin{aligned}
    \zeta(x) &\triangleq g(x) + \Sigma^{-1}y\\
    \hat{z}(x) &\triangleq \Pi_D\left(g(x) + \Sigma^{-1}y\right)  \\
    &= \Pi_D(\zeta(x)) \\
    d(x) &\triangleq \zeta(x) - \Pi_D(\zeta(x)) \\
    &= \zeta(x) - \hat{z}(x)  \\
    \hat{y}(x) &\triangleq \Sigma\, d(x) \\
    &= \Sigma\, \left(g(x) + \Sigma^{-1}y - \Pi_D\left(g(x) + \Sigma^{-1}y\right)\right)\\
    e(x) &\triangleq g(x) - \hat z(x)
    \end{aligned}
@f]

The result of the PANOC algorithm is the triple @f$(\hat x^k,\;\hat y(\hat x^k),\;\hat z(\hat x^k))@f$.

The following graph visualizes the dependencies between the different values
used in a PANOC iteration.

@image html expression-dep.gv.svg