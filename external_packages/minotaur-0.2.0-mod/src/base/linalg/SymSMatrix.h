
// Symmetric Sparse Matrix

#include "SymMatrix.h"
#include "SVector.h"

namespace Minotaur {
  namespace LinearAlgebra {
    //   Row major
    template <class T> class SymCSRMatrix : public SymMatrix<T> {
    private:
      int nrow;
      int ncol;
      T * dat;

    public:
      SymCSRMatrix();
    };

    //   Col major
    template <class T> class SymCSCMatrix : public SymMatrix<T> {
    private:
      int nrow;
      int ncol;
      T * dat;

    public:
      SymCSCMatrix();
    };

    //  Outer product matrix

    template <class T> class SymSOMatrix : public SymMatrix<T> {
    private:
      SVector<T> v;
    };
  };
};

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
