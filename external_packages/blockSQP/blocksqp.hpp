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
#include "string.h"
#include <qpOASES.hpp>
#include <set>

namespace blocksqp {

  typedef char PATHSTR[4096];

  /**
   * \brief Class for easy access of elements of a dense matrix.
   * \author Dennis Janka
   * \date 2012-2015
   */
  class Matrix {
    private:
    int malloc( void );                                           ///< memory allocation
    int free( void );                                             ///< memory free

  public:
    int m;                                                        ///< internal number of rows
    int n;                                                        ///< internal number of columns
    int ldim;                                                     ///< internal leading dimension not necesserily equal to m or n
    double *array;                                                ///< array of how the matrix is stored in the memory
    int tflag;                                                    ///< 1 if it is a Teilmatrix

    Matrix( int = 1, int = 1, int = -1 );                         ///< constructor with standard arguments
    Matrix( int, int, double*, int = -1 );
    Matrix( const Matrix& A );
    virtual ~Matrix( void );

    int M( void ) const;                                          ///< number of rows
    int N( void ) const;                                          ///< number of columns
    int LDIM( void ) const;                                       ///< leading dimensions
    double *ARRAY( void ) const;                                  ///< returns pointer to data array
    int TFLAG( void ) const;                                      ///< returns this->tflag (1 if it is a submatrix and does not own the memory and 0 otherwise)

    virtual double &operator()( int i, int j );                   ///< access element i,j of the matrix
    virtual double &operator()( int i, int j ) const;
    virtual double &operator()( int i );                          ///< access element i of the matrix (columnwise)
    virtual double &operator()( int i ) const;
    virtual Matrix &operator=( const Matrix &A );                 ///< assignment operator

    Matrix &Dimension( int, int = 1, int = -1 );                  ///< set dimension (rows, columns, leading dimension)
    Matrix &Initialize( double (*)( int, int ) );                 ///< set matrix elements i,j to f(i,j)
    Matrix &Initialize( double );                                 ///< set all matrix elements to a constant

    /// Returns just a pointer to the full matrix
    Matrix& Submatrix( const Matrix&, int, int, int = 0, int = 0 );
    /// Matrix that points on <tt>ARRAY</tt>
    Matrix& Arraymatrix( int M, int N, double* ARRAY, int LDIM = -1 );

    /** Flag == 0: bracket output
     * Flag == 1: Matlab output
     * else: plain output */
    const Matrix &Print( FILE* = stdout,   ///< file for output
                         int = 13,       ///< number of digits
                         int = 1         ///< Flag for format
                         ) const;
  };

  /**
   * \brief Class for easy access of elements of a dense symmetric matrix.
   * \author Dennis Janka
   * \date 2012-2015
   */
  class SymMatrix : public Matrix {
  protected:
    int malloc( void );
    int free( void );

  public:
    SymMatrix( int = 1 );
    SymMatrix( int, double* );
    SymMatrix( int, int, int );
    SymMatrix( int, int, double*, int = -1 );
    SymMatrix( const Matrix& A );
    SymMatrix( const SymMatrix& A );
    virtual ~SymMatrix( void );

    virtual double &operator()( int i, int j );
    virtual double &operator()( int i, int j ) const;
    virtual double &operator()( int i );
    virtual double &operator()( int i ) const;

    SymMatrix &Dimension( int = 1 );
    SymMatrix &Dimension( int, int, int );
    SymMatrix &Initialize( double (*)( int, int ) );
    SymMatrix &Initialize( double );

    SymMatrix& Submatrix( const Matrix&, int, int, int = 0, int = 0 );
    SymMatrix& Arraymatrix( int, double* );
    SymMatrix& Arraymatrix( int, int, double*, int = -1 );
  };

  Matrix Transpose( const Matrix& A); ///< Overwrites \f$ A \f$ with its transpose \f$ A^T \f$
  Matrix &Transpose( const Matrix &A, Matrix &T ); ///< Computes \f$ T = A^T \f$
  double delta( int, int );

  /**
   * \brief Base class for problem specification as required by SQPMethod.
   * \author Dennis Janka
   * \date 2012-2015
   */
  class Problemspec {
    /*
     * VARIABLES
     */
  public:
    int         nVar;               ///< number of variables
    int         nCon;               ///< number of constraints
    int         nnCon;              ///< number of nonlinear constraints

    double      objLo;              ///< lower bound for objective
    double      objUp;              ///< upper bound for objective
    Matrix      bl;                 ///< lower bounds of variables and constraints
    Matrix      bu;                 ///< upper bounds of variables and constraints

    int         nBlocks;            ///< number of separable blocks of Lagrangian
    int*        blockIdx;           ///< [blockwise] index in the variable vector where a block starts

    /*
     * METHODS
     */
  public:
    Problemspec( ){};
    virtual ~Problemspec( ){};

    /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (dense version).
    virtual void initialize( Matrix &xi,            ///< optimization variables
                             Matrix &lambda,        ///< Lagrange multipliers
                             Matrix &constrJac      ///< constraint Jacobian (dense)
                             ){};

    /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (sparse version).
    virtual void initialize( Matrix &xi,            ///< optimization variables
                             Matrix &lambda,        ///< Lagrange multipliers
                             double *&jacNz,        ///< nonzero elements of constraint Jacobian
                             int *&jacIndRow,       ///< row indices of nonzero elements
                             int *&jacIndCol        ///< starting indices of columns
                             ){};

    /// Evaluate objective, constraints, and derivatives (dense version).
    virtual void evaluate( const Matrix &xi,        ///< optimization variables
                           const Matrix &lambda,    ///< Lagrange multipliers
                           double *objval,          ///< objective function value
                           Matrix &constr,          ///< constraint function values
                           Matrix &gradObj,         ///< gradient of objective
                           Matrix &constrJac,       ///< constraint Jacobian (dense)
                           SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                           int dmode,               ///< derivative mode
                           int *info                ///< error flag
                           ){};

    /// Evaluate objective, constraints, and derivatives (sparse version).
    virtual void evaluate( const Matrix &xi,        ///< optimization variables
                           const Matrix &lambda,    ///< Lagrange multipliers
                           double *objval,          ///< objective function value
                           Matrix &constr,          ///< constraint function values
                           Matrix &gradObj,         ///< gradient of objective
                           double *&jacNz,          ///< nonzero elements of constraint Jacobian
                           int *&jacIndRow,         ///< row indices of nonzero elements
                           int *&jacIndCol,         ///< starting indices of columns
                           SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                           int dmode,               ///< derivative mode
                           int *info                ///< error flag
                           ){};

    /// Short cut if no derivatives are needed
    virtual void evaluate( const Matrix &xi,        ///< optimization variables
                           double *objval,          ///< objective function value
                           Matrix &constr,          ///< constraint function values
                           int *info                ///< error flag
                           );

    /*
     * Optional Methods
     */
    /// Problem specific heuristic to reduce constraint violation
    virtual void reduceConstrVio( Matrix &xi,       ///< optimization variables
                                  int *info         ///< error flag
                                  ){ *info = 1; };

    /// Print information about the current problem
    virtual void printInfo(){};
  };

  //  Declaration of general purpose routines for matrix and vector computations
  double l1VectorNorm( const Matrix &v );
  double l2VectorNorm( const Matrix &v );
  double lInfVectorNorm( const Matrix &v );

  double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl, const Matrix &weights );
  double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
  double l2ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
  double lInfConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );

  double adotb( const Matrix &a, const Matrix &b );
  void Atimesb( const Matrix &A, const Matrix &b, Matrix &result );
  void Atimesb( double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result );

  int calcEigenvalues( const Matrix &B, Matrix &ev );
  double estimateSmallestEigenvalue( const Matrix &B );
  int inverse( const Matrix &A, Matrix &Ainv );

} // namespace blocksqp

#endif // BLOCKSQP_HPP
