/// \cond INTERNAL

extern "C" {

  extern void sqic(
   const int *m, // Number of constraints + 1 (for the objective)
   const int* n, // Number of decision variables
   const int* nnzA, // Number of nonzeros in objective-augmented linear constraint matrix  A
   const int *indA, // colind of Compressed Column Storage A , length: nnzA
   const int *locA, // row of  Compressed Column Storage A, length n + 1
   const double *valA, // Values of A
   const double* bl, // Lower bounds to decision variables + objective
   const double* bu, // Upper bounds to decision variables + objective
   const int *hEtype, // ?
   const int *hs, // ?
   double *x,  // Decision variables + evaluated linear constraints ((initial + optimal), length n+m
   double *pi, // ?
   double *rc, // Multipliers (initial + optimal), length n+m
   const int* nnzH, // Number of nonzeros in full hessian H
   const int* indH, // colind of Compressed Column Storage H , length: nnzH
   const int* locH, // row of  Compressed Column Storage H, length n + 1
   double* valH
   );

  extern void sqicSolve(
   double* Obj // Output: hessian part of the resulting objective
  );

  extern void sqicSolveStabilized(
   double* Obj, // Output: hessian part of the resulting objective
   double *mu,
   int *lenpi,
   double* piE
  );

  extern void sqicDestroy();
}
/// \endcond
