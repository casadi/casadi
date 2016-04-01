// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 


#ifndef MINOTAURCHOL_H
#define MINOTAURCHOL_H

#include "LinearFunction.h"
#include "QuadraticFunction.h"

namespace Minotaur {


  // /**
  // The CholCalculator class is used to calculate Cholesky factorization 
  //  of given quadratic function
  // */
  class CholCalculator {
    public:

      /// default constructor
      CholCalculator();

      // /**
      // Construct using a quadratic function.
      // */
      CholCalculator(ConstQuadraticFunctionPtr qf);

      // /**
      // Destroy
      // */
      ~CholCalculator() {};

    private:

      // /**
      // The quadratic function for whose Hessian we wish to find the eigen values.
      // */
      ConstQuadraticFunctionPtr qf_;

      // /**
      // Dimension of the square matrix
      // */
      UInt n_;

      // /**
      // The square matrix is stored as a single array. The element A[i,j] can
      // be accessed at A_[i+j*n_]. And element A_[i] = A[i mod n_, i/n_].
      // */
      double *A_;

      ///
      double abstol_;

  };


}
#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
