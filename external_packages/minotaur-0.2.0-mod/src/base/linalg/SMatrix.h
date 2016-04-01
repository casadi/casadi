// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 -- 2010 The MINOTAUR Team.
// 

// /**
// \file SMatrix.h
// \brief Declare sparse matrices
// \author Ashutosh Mahajan, Argonne National Laboratory
//
// This file contains declarations for methods of sparse matrices.
// */


#ifndef MINOTAURSMATRIX_H
#define MINOTAURSMATRIX_H

#include "Types.h"
#include "Matrix.h"
#include "SVector.h"

namespace Minotaur {
  namespace LinearAlgebra {

    typedef enum {
      ColumnOrdered,
      RowOrdered
    } MatrixOrderType;

    /// Class for Compressed Sparse matrix. We probably don't need different classes for
    /// column and row orders; instead a flag will suffice.
    template <class T> class CSMatrix : public Matrix<T> {
      public:
        /// Default constructor
        CSMatrix(MatrixOrderType order=ColumnOrdered);

        /// Append a new row
        void addRow(std::vector<const T > &row);

        /// Append a new column
        void addCol(std::vector<const T > &col);

        /// Get major-size
        int getNumMajors() { return maj_size_; };

        /// Get minor-size
        int getNumMinors() { return min_size_; };

        /// Get number of columns
        int getNumCols() 
        { return (order_t_==ColumnOrdered) ? maj_size_ : min_size_; };

        /// Get number of rows
        int getNumRows() 
        { return (order_t_==RowOrdered) ? maj_size_ : min_size_; };

        /// Get number of nonzeros
        int getNumNzs() { return nzs_; };

        /// Get the type of ordering
        MatrixOrderType getOrderType() { return order_t_; };

      private:
        /// Size of the major ordering
        int maj_size_;

        /// Size of the minor ordering
        int min_size_;

        /// Number of nonzeros
        int nzs_;

        /// Type of ordering
        MatrixOrderType order_t_;
    };



    //  Compressed Spare Row Matrix
    // template <class T> class CSRMatrix : public Matrix<T> {
    // private:
    //   int nrow;
    //   int ncol;
    //   T * dat;

    // public:
    //   CSRMatrix();
    // };

    // //  Compressed Sparse Col Matrix
    // template <class T> class CSCMatrix : public Matrix<T> {
    // private:
    //   int nrow;
    //   int ncol;
    //   T * dat;

    // public:
    //   CSCMatrix();
    // };

    // // Sparse Outer produce matrix
    // template <class T> class SOMatrix : public Matrix<T> {
    //   SVector<T> v1;
    //   SVector<T> v2;
    // };
  };
};

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
