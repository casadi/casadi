%{
#include "superlu_interface/superlu.hpp"
%}

namespace CasADi{
  
  /// Type of solve call
  enum Factorization{DOFACT, SAMEPATTERN, SAMEPATTERN_SAMEROWPERM, FACTORED};

  /// LU factorizer
  class SuperLU : public FX{
  
    public:

      /// Create a linear solver given a sparsity pattern
      SuperLU(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1);
  
      /// Solve
      void solve(Factorization fact);
  
  };
} // namespace CasADi

