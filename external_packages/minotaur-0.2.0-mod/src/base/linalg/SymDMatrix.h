
// Symmetric Dense Matrix

#include "SymMatrix.h"
#include "DVector.h"

namespace Minotaur {
  namespace LinearAlgebra {
    //   Row major
    template <class T> class SymDRMatrix : public SymMatrix<T> {
    private:
      int nrow;
      int ncol;
      T * dat;

    public:
      SymDRMatrix();
    };

    //   Col major
    template <class T> class SymDCMatrix : public SymMatrix<T> {
    private:
      int nrow;
      int ncol;
      T * dat;

    public:
      SymDCMatrix();
    };

    template <class T> class SymDOMatrix : public SymMatrix<T> {
    private:
      DVector<T> v;
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
