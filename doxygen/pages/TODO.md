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