# Test results

PDF files contain the figures, XLSX files contain the Excel sheets with the full
results (one sheet per solver).

## Standard PANOC -vs- 2nd order PANOC (L-BFGS, finite differences)

[[PDF](panoc-lbfgs-23-4-ls-margin-vs-panoc-2nd-lbfgs-23-4-fd5.pdf)] 
[[XLSX](panoc-lbfgs-23-4-ls-margin-vs-panoc-2nd-lbfgs-23-4-fd5.xlsx)]

Hessian of the Augmented Lagrangian is approximated using L-BFGS on the gradient,
Hessian-vector product is computed using finite differences on the gradient.

Problems converged: 155 vs 160 out of 219  
Total time: 893s vs 880s  
Time of converged problems: 53s vs 74s

Conclusions: Second order PANOC is faster and requires fewer iterations.  
It is somewhat surprising that the average τ is so much smaller.

## 2nd order PANOC (L-BFGS, AD) -vs- 2nd order PANOC (L-BFGS, finite differences)

[[PDF](panoc-2nd-lbfgs-23-4-aug-hess-vs-panoc-2nd-lbfgs-23-4-fd5.pdf)] 
[[XLSX](panoc-2nd-lbfgs-23-4-aug-hess-vs-panoc-2nd-lbfgs-23-4-fd5.pdf)]

Hessian-vector product of the Lagrangian is computed using AD, then ∇g_i(x) is
computed using AD as well to compute the Hessian-vector product of the Augmented
Lagrangian. This basically corresponds to computing the entire Jacobian of g, 
which is expensive.

Problems converged: 159 vs 160 out of 219  
Total time: 1442s vs 880s  
Time of converged problems: 128s vs 74s

Conclusions: The number of iterations is roughly the same, but FD is faster, 
especially for larger problems. FD does require more objective gradient
evaluations, but this is compensated by Hessian-vector products and constraints
Jacobians, which are much more expensive.

## L-BFGS-B -vs- 2nd order PANOC (L-BFGS, finite differences)

[[PDF](lbfgsbpp-23-4-vs-panoc-2nd-lbfgs-23-4-fd5.pdf)] 
[[XLSX](lbfgsbpp-23-4-vs-panoc-2nd-lbfgs-23-4-fd5.pdf)]

Problems converged: 164 vs 160 out of 219  
Total time: 538s vs 880s  
Time of converged problems: 85s vs 74s

Conclusions: PANOC requires many more iterations than L-BFGS-B, but PANOC
iterations are cheaper. For small problems, PANOC is faster, for large problems,
L-BFGS-B is much faster than PANOC.  
PANOC also requires significantly fewer ALM iterations for some reason.

|     | Problems converged | Total time [s] | Time of converged problems [s] |
|:----|:------------------:|:----------:|:--------------------------:|
| Standard PANOC | 155 | 893 | 53 |
| 2nd order PANOC (L-BFGS, FD) | 160 | 880 | 74 |
| 2nd order PANOC (L-BFGS, AD) | 169 | 1442 | 128 |
| L-BFGS-B | 164 | 538 | 85 |