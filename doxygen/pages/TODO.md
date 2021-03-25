# TODO

 - Implement careful BFGS updates
 - In the 2D QP test, why does L-BFGS only have a correct inverse Hessian 
   estimate after 3 iterations instead of 2?
 - Why is a margin needed in the Lipschitz check of PANOC?
 - Implement memory pool to optimize vector allocations
 - Handle NaN correctly in min/max functions
 - Check if L-BFGS should use r or p
 - Adapt SpecializedLBFGS to use new generic API
 - LBFGSParam

 # Questions for Panos

- Why no accelerated proximal gradient in PANOC?  
  yₖ = xₖ + βₖ(xₖ - xₖ₋₁)  
  xₖ₊₁ = proxₕ(yₖ - γₖ ∇ψ(yₖ))
- More efficient stopping criterion AA-PGA?
- Regularization of least squares Anderson acceleration?
- Underdetermined system Anderson acceleration?
- QR with column pivoting? (Much more expensive, faster alternative?)
- Incremental condition number estimation for dropping columns?