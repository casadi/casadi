/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp.hpp
 * \author Dennis Janka, Joel Andersson
 * \date 2012-2015
 */

#ifndef BLOCKSQP_HPP
#define BLOCKSQP_HPP

#include "math.h"
#include "stdio.h"

namespace blocksqp {
  /**
   * \brief Class for easy access of elements of a dense matrix.
   * \author Dennis Janka
   * \date 2012-2015
   */
  class Matrix {
  public:
    int m;                                                        ///< internal number of rows
    int n;                                                        ///< internal number of columns
    int ldim;                                                     ///< internal leading dimension not necesserily equal to m or n
    double *array;                                                ///< array of how the matrix is stored in the memory
    int tflag;                                                    ///< 1 if it is a Teilmatrix
  private:
    int malloc();                                           ///< memory allocation
    int free();                                             ///< memory free
  public:
    Matrix(int = 1, int = 1, int = -1);                         ///< constructor with standard arguments
    Matrix(int, int, double*, int = -1);
    Matrix(const Matrix& A);
    virtual ~Matrix();

    virtual double &operator()(int i, int j);                   ///< access element i,j of the matrix
    virtual double &operator()(int i, int j) const;
    virtual double &operator()(int i);                          ///< access element i of the matrix (columnwise)
    virtual double &operator()(int i) const;
    virtual Matrix &operator=(const Matrix &A);                 ///< assignment operator

    Matrix &Dimension(int, int = 1, int = -1);                  ///< set dimension (rows, columns, leading dimension)
    Matrix &Initialize(double (*)(int, int));                 ///< set matrix elements i,j to f(i,j)
    Matrix &Initialize(double);                                 ///< set all matrix elements to a constant

    /// Returns just a pointer to the full matrix
    Matrix& Submatrix(const Matrix&, int, int, int = 0, int = 0);
  };

  /**
   * \brief Class for easy access of elements of a dense symmetric matrix.
   * \author Dennis Janka
   * \date 2012-2015
   */
  class SymMatrix : public Matrix {
  protected:
    int malloc();
    int free();

  public:
    SymMatrix(int = 1);
    SymMatrix(int, double*);
    SymMatrix(int, int, int);
    SymMatrix(int, int, double*, int = -1);
    SymMatrix(const Matrix& A);
    SymMatrix(const SymMatrix& A);
    virtual ~SymMatrix();

    virtual double &operator()(int i, int j);
    virtual double &operator()(int i, int j) const;
    virtual double &operator()(int i);
    virtual double &operator()(int i) const;

    SymMatrix &Dimension(int = 1);
    SymMatrix &Dimension(int, int, int);
    SymMatrix &Initialize(double (*)(int, int));
    SymMatrix &Initialize(double);

    SymMatrix& Submatrix(const Matrix&, int, int, int = 0, int = 0);
  };

  //  Declaration of general purpose routines for matrix and vector computations
  double l1VectorNorm(const Matrix &v);
  double l2VectorNorm(const Matrix &v);
  double lInfVectorNorm(const Matrix &v);
  double lInfConstraintNorm(const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl);

  double adotb(const Matrix &a, const Matrix &b);
  void Atimesb(const Matrix &A, const Matrix &b, Matrix &result);
  void Atimesb(double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result);

} // namespace blocksqp

#endif // BLOCKSQP_HPP
