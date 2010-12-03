%{
#include "superlu_interface/superlu.hpp"
%}

namespace CasADi{
  
  /// LU factorizer
  class SuperLU : public LinearSolver{
  
    public:
      
      /// Default constructor
      SuperLU();

      /// Create a linear solver given a sparsity pattern
      SuperLU(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1);
  
  };
} // namespace CasADi

