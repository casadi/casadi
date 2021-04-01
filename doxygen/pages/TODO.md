# TODO

- ~~Implement cautious BFGS updates~~
- ~~In the 2D QP test, why does L-BFGS only have a correct inverse Hessian 
  estimate after 3 iterations instead of 2?~~
- Why is a margin needed in the Lipschitz check of PANOC?
- Implement memory pool to optimize vector allocations
- Handle NaN correctly in min/max functions
- ~~Check if L-BFGS should use r or p~~
- ~~Adapt SpecializedLBFGS to use new generic API~~
- LBFGSParam
- xₖ₋₁ is incorrect when previous SL-BFGS update was rejected :(
- Automatically use unconstrained L-BFGS inner solver if the problem has no
  box constraints on x
- MPC-specific tests (mass-spring benchmark or ask Joris Gillis (?) or 
  Ben Hermans)
- Implement better heuristic scaling/preconditioning
- Investigate how Ipopt and other established solvers handle this scaling
- ALM must detect if PANOC fails to reach tolerance, it makes little sense to 
  increase the penalty based on a bad result from PANOC
- Anderson acceleration 
- Alternative fast directions
- Find “optimal” parameters for ALM, first order should be less aggressive than
  e.g. QPALM.
- Look at Christian Kanzow's (?) solver and 
  [Lancelot solver](https://www.numerical.rl.ac.uk/lancelot/blurb.html)
- Combine function and gradient evaluations?

## Open questions

- ALM: when can we be certain that we'll converge to a feasible point?
- Investigate second order PANOC (specifically how to compute 6.1.c efficiently)
- How do multiple shooting formulations compare to single shooting?
- How can we balance ε and δ? Keep λ constant if ε is too large?
- Can we accept unit step size if ‖Rγ(xₖ₊₁)‖ < c ‖Rγ(xₖ)‖, c ∈ (0, 1)
- When PANOC suddenly converges after stalling, why?

## Questions for Panos

- Why no accelerated proximal gradient in PANOC?  
  yₖ = xₖ + βₖ(xₖ - xₖ₋₁)  
  xₖ₊₁ = proxₕ(yₖ - γₖ ∇ψ(yₖ))
- More efficient stopping criterion AA-PGA?
- Regularization of least squares Anderson acceleration?
- Underdetermined system Anderson acceleration?
- QR with column pivoting? (Much more expensive, faster alternative?)
- Incremental condition number estimation for dropping columns?
- Is checking that ψₖ == ψₖ₋₁ for N iterations a valid check instead of
  xₖ == xₖ₋₁?
- Do we really have to flush L-BFGS if γ changes? 
  Since Rγ(x) = γ⁻¹ (x - Tγ(x)) = γ⁻¹ (x - Π(x - γ∇ψ(x))), and if we're not 
  at the boundaries of the feasible set, Π(x - γ∇ψ(x)) = x - γ∇ψ(x), and 
  Rγ(x) = γ⁻¹ (-γ∇ψ(x)) = ∇ψ(x) which is independent of γ.
  This might be useful for Anderson as well.  
  _Update: I tried this for Anderson, and it nicely speeds up the algorithm and results in more test problems converging._
- Trust region? Lancelot?
- Does it make sense to include the Wolfe curvature condition in our line
  search? Maybe this results in fewer L-BFGS updates being rejected? This
  condition can be evaluated relatively inexpensively.