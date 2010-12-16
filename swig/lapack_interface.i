%{
#include "interfaces/lapack/lapack_lu_dense.hpp"
#include "interfaces/lapack/lapack_qr_dense.hpp"
%}

namespace CasADi{
  
  /// Linear solver using dense LU factorization
  class LapackLUDense : public LinearSolver{
  
    public:
      
      /// Default constructor
      LapackLUDense();

      /// Create a linear solver given a sparsity pattern
      LapackLUDense(int nrow, int ncol, int nrhs=1);
  
  };

  /// Linear solver using dense QR factorization
  class LapackQRDense : public LinearSolver{
  
    public:
      
      /// Default constructor
      LapackQRDense();

      /// Create a linear solver given a sparsity pattern
      LapackQRDense(int nrow, int ncol, int nrhs=1);
  
  };
} // namespace CasADi

