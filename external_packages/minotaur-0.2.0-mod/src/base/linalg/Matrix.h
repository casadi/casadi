/*
 *     MINOTAUR -- It's only 1/2 bull
 *
 *     (C)opyright 2008-- The MINOTAUR Team
 */

#ifndef MINOTAURMATRIX_H
#define MINOTAURMATRIX_H

#include "Types.h"

//
// yell at ashu if you dont like this. Default Matrix class defined here
// assumes fully dense matrix. Derived classes are sparse, symmetric etc.
//
// The matrix should:
//     # provide methods to get number of rows, cols, nonzeros.
//     # get element (i,j).
//     # change the values of elements or report that they cannot be changed.
//     # tell what kind of matrix is it (sparse, dense, symmetric ... ).
//     # tell if it is ordered by rows or columns.
//     # add your wish here.
//

namespace Minotaur {
  namespace LinearAlgebra {
    template <class T> class Matrix {
      public:
        UInt getNumRows() const { return nrows_; }
        UInt getNumCols() const { return ncols_; }
        UInt getNumElements() const { return ncols_; } 

      private:
        UInt nrows_;
        UInt ncols_;
        UInt nelements_;
    };
    typedef boost::shared_ptr< Matrix<double> > DoubleMatrixPtr;
  }
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
