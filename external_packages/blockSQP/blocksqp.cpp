/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp.cpp
 * \author Dennis Janka, Joel Andersson
 * \date 2012-2015, 2016
 *
 */


#include "blocksqp.hpp"

namespace blocksqp {


  /**
   * Compute scalar product aTb
   */
  double adotb(const Matrix &a, const Matrix &b) {
    double norm = 0.0;

    if (a.n != 1 || b.n != 1) {
        printf("a or b is not a vector!\n");
    } else if (a.m != b.m) {
        printf("a and b must have the same dimension!\n");
    } else {
      for (int k=0; k<a.m; k++)
        norm += a(k) * b(k);
    }
    return norm;
  }

  /**
   * Compute the matrix vector product for a column-compressed sparse matrix A with a vector b and store it in result
   */
  void Atimesb(double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result) {
    int nCol = b.m;
    int nRow = result.m;
    int i, k;

    for (i=0; i<nRow; i++)
      result(i) = 0.0;

    for (i=0; i<nCol; i++) {
        // k runs over all elements in one column
        for (k=AIndCol[i]; k<AIndCol[i+1]; k++)
          result(AIndRow[k]) += Anz[k] * b(i);
    }
  }

  /**
   * Compute the matrix vector product A*b and store it in result
   */
  void Atimesb(const Matrix &A, const Matrix &b, Matrix &result) {
    result.Initialize(0.0);
    for (int i=0; i<A.m; i++)
      for (int k=0; k<A.n; k++)
        result(i) += A(i, k) * b(k);
  }

  double l1VectorNorm(const Matrix &v) {
    double norm = 0.0;

    if (v.n != 1) {
        printf("v is not a vector!\n");
    } else {
        for (int k=0; k<v.m; k++)
          norm += fabs(v(k));
    }
    return norm;
  }

  double l2VectorNorm(const Matrix &v) {
    double norm = 0.0;
    if (v.n != 1) {
      printf("v is not a vector!\n");
    } else {
      for (int k=0; k<v.m; k++)
        norm += v(k)* v(k);
    }
    return sqrt(norm);
  }

  double lInfVectorNorm(const Matrix &v) {
    double norm = 0.0;
    if (v.n != 1) {
        printf("v is not a vector!\n");
    } else {
      for (int k=0; k<v.m; k++)
        if (fabs(v(k)) > norm)
          norm = fabs(v(k));
    }
    return norm;
  }

  /**
   * Calculate l_Infinity norm of constraint violations
   */
  double lInfConstraintNorm(const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl) {
    double norm = 0.0;
    int i;
    int nVar = xi.m;
    int nCon = constr.m;

    // Violation of simple bounds
    for (i=0; i<nVar; i++) {
      if (xi(i) - bu(i) > norm)
        norm = xi(i) - bu(i);
      else if (bl(i) - xi(i) > norm)
        norm = bl(i) - xi(i);
    }

    // Find out the largest constraint violation
    for (i=0; i<nCon; i++) {
      if (constr(i) - bu(nVar+i) > norm)
        norm = constr(i) - bu(nVar+i);
      if (bl(nVar+i) - constr(i) > norm)
        norm = bl(nVar+i) - constr(i);
    }

    return norm;
  }

  void Error(const char *F) {
    printf("Error: %s\n", F);
  }

  /* ----------------------------------------------------------------------- */

  int Matrix::malloc() {
    int len;

    if (tflag)
      Error("malloc cannot be called with Submatrix");

    if (ldim < m)
      ldim = m;

    len = ldim*n;

    if (len == 0)
      array = 0;
    else
      if ((array = new double[len]) == 0)
        Error("'new' failed");

    return 0;
  }


  int Matrix::free() {
    if (tflag) Error("free cannot be called with Submatrix");
    if (array != 0) delete[] array;
    return 0;
  }

  double &Matrix::operator()(int i, int j) {
    return array[i+j*ldim];
  }

  double &Matrix::operator()(int i, int j) const {
    return array[i+j*ldim];
  }

  double &Matrix::operator()(int i) {
    return array[i];
  }

  double &Matrix::operator()(int i) const {
    return array[i];
  }

  Matrix::Matrix(int M, int N, int LDIM) {
    m = M;
    n = N;
    ldim = LDIM;
    tflag = 0;
    malloc();
  }

  Matrix::Matrix(int M, int N, double *ARRAY, int LDIM) {
    m = M;
    n = N;
    array = ARRAY;
    ldim = LDIM;
    tflag = 0;
    if (ldim < m) ldim = m;
  }

  Matrix::Matrix(const Matrix &A) {
    m = A.m;
    n = A.n;
    ldim = A.ldim;
    tflag = 0;
    malloc();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n ; j++)
        (*this)(i, j) = A(i, j);
  }

  Matrix &Matrix::operator=(const Matrix &A) {
    if (this != &A) {
        if (!tflag) {
            free();

            m = A.m;
            n = A.n;
            ldim = A.ldim;

            malloc();

            for (int i = 0; i < m; i++)
              for (int j = 0; j < n ; j++)
                (*this)(i, j) = A(i, j);
          } else {
            if (m != A.m || n != A.n)
              Error("= operation not allowed");

            for (int i = 0; i < m; i++)
              for (int j = 0; j < n ; j++)
                (*this)(i, j) = A(i, j);
          }
      }
    return *this;
  }

  Matrix::~Matrix() {
    if (!tflag) free();
  }

  Matrix &Matrix::Dimension(int M, int N, int LDIM) {
    if (M != m || N != n || (LDIM != ldim && LDIM != -1)) {
      if (tflag) {
        Error("Cannot set new dimension for Submatrix");
      } else {
        free();
        m = M;
        n = N;
        ldim = LDIM;

        malloc();
      }
    }
    return *this;
  }

  Matrix &Matrix::Initialize(double (*f)(int, int)) {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        (*this)(i, j) = f(i, j);

    return *this;
  }

  Matrix &Matrix::Initialize(double val) {
    for (int i=0; i < m; i++)
      for (int j=0; j < n; j++)
        (*this)(i, j) = val;
    return *this;
  }

  Matrix &Matrix::Submatrix(const Matrix &A, int M, int N, int i0, int j0) {
    if (i0 + M > A.m || j0 + N > A.n) Error("Cannot create Submatrix");
    if (!tflag) free();
    tflag = 1;
    m = M;
    n = N;
    array = &A.array[i0+j0*A.ldim];
    ldim = A.ldim;
    return *this;
  }

  int SymMatrix::malloc() {
    int len = m*(m+1)/2.0;
    if (len == 0) {
      array = 0;
    } else {
      if ((array = new double[len]) == 0) Error("'new' failed");
    }
    return 0;
  }

  int SymMatrix::free() {
    if (array != 0) delete[] array;
    return 0;
  }

  double &SymMatrix::operator()(int i, int j) {
    int pos;

    if (i < j)//reference to upper triangular part
      pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
      pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
  }


  double &SymMatrix::operator()(int i, int j) const {
    int pos;

    if (i < j)//reference to upper triangular part
      pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
      pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
  }

  double &SymMatrix::operator()(int i) {
    return array[i];
  }

  double &SymMatrix::operator()(int i) const {
    return array[i];
  }

  SymMatrix::SymMatrix(int M) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;
    malloc();
  }

  SymMatrix::SymMatrix(int M, double *ARRAY) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;
    malloc();
    array = ARRAY;
  }


  SymMatrix::SymMatrix(int M, int N, int LDIM) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;
    malloc();
  }


  SymMatrix::SymMatrix(int M, int N, double *ARRAY, int LDIM) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;
    malloc();
    array = ARRAY;
  }


  SymMatrix::SymMatrix(const Matrix &A) {
    m = A.m;
    n = A.m;
    ldim = A.m;
    tflag = 0;
    malloc();
    for (int j=0; j<m; j++)
      for (int i=j; i<m; i++)
        (*this)(i, j) = A(i, j);
  }


  SymMatrix::SymMatrix(const SymMatrix &A) {
    m = A.m;
    n = A.n;
    ldim = A.ldim;
    tflag = 0;
    malloc();
    for (int j=0; j<m; j++)
      for (int i=j; i<m; i++)
        (*this)(i, j) = A(i, j);
  }


  SymMatrix::~SymMatrix() {
    if (!tflag)
      free();
  }

  SymMatrix &SymMatrix::Dimension(int M) {
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
  }


  SymMatrix &SymMatrix::Dimension(int M, int N, int LDIM) {
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
  }


  SymMatrix &SymMatrix::Initialize(double (*f)(int, int)) {
    for (int j=0; j<m; j++)
      for (int i=j; i<n ; i++)
        (*this)(i, j) = f(i, j);

    return *this;
  }


  SymMatrix &SymMatrix::Initialize(double val) {
    for (int j=0; j<m; j++)
      for (int i=j; i<n; i++)
        (*this)(i, j) = val;

    return *this;
  }

  SymMatrix &SymMatrix::Submatrix(const Matrix &A, int M, int N, int i0, int j0) {
    Error("SymMatrix doesn't support Submatrix");
    return *this;
  }

} // namespace blocksqp
