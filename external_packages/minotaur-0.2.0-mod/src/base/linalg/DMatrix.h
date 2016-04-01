
// Dense Matrix

#include "Types.h"
#include "Matrix.h"
#include "DVector.h"

namespace Minotaur {
  namespace LinearAlgebra {
    //   Row major
    template <class T> class DRMatrix : public Matrix<T> {
    private:
      UInt nrow;
      UInt ncol;
      T *dat;

    public:
      DRMatrix();
    };

    //   Col major
    template <class T> class DCMatrix : public Matrix<T> {
    private:
      int nrow;
      int ncol;
      T * dat;

    public:
      DCMatrix();
    };

    //  Outer product matrix
    template <class T> class DOMatrix : public Matrix<T> {
    private:
      DVector<T> v1;
      DVector<T> v2;
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
