# Inner solvers

This page provides a high level overview of the inner solvers supported by alpaqa, and how to use them.

### Forward-backward iterations

Both PANOC @cite stella2017panoc and PANTR solve problems of the form
@f[
\begin{equation}
    \begin{aligned}
        \underset{x}{\text{minimize}}\; \psi(x) + h(x) \quad \quad \psi : \Rn \rightarrow \R, \; \; h : \Rn \rightarrow \overline{\R}
    \end{aligned}
    \label{eq:problem-orig}
    \tag{P}
\end{equation}
@f]
where @f$\psi@f$ is (locally) Lipschitz-smooth and @f$h@f$ is proper, lowersemicontinuous and convex.
Note that the (possibly nonsmooth) term @f$h@f$ can be used to encode constraints or include a regularization term.

A well-known strategy for solving @f$\eqref{eq:problem-orig}@f$ consists of iteratively applying the forward-backward operator
@f[
\begin{equation}
    T_\gamma(x) = \prox_{\gamma h}(x - \gamma \nabla \psi(x)),
    \label{eq:fbs}
    \tag{FB}
\end{equation}
@f]
where @f$\prox_{\gamma h}(x) = \underset{u}{\text{argmin}} \{ h(u) + \frac{1}{2\gamma} \Vert u - x \Vert^2 \}@f$ denotes the proximal operator of @f$h@f$ with step size @f$\gamma@f$.
Remark that this scheme only requires evaluations of @f$\nabla \psi@f$ and @f$\prox_{\gamma h}@f$, which are assumed to be efficiently evaluated.
Using the same oracle, both PANOC and PANTR aim to accelerate the standard forward-backward iterations @f$\eqref{eq:fbs}@f$.

## PANOC

PANOC @cite stella2017panoc combines forward-backward iterations @f$\eqref{eq:fbs}@f$ with a quasi-Newton linesearch procedure to attain fast asymptotic convergence.

TO DO: API

## PANTR

PANTR is similar to PANOC, but replaces the quasi-Newton directions by solutions to trust-region subproblems, which can be seen as regularized Newton updates.

TO DO: API