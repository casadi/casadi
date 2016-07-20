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

// LAPACK routines
extern "C" {
  void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
               double *w, double *work, int *lwork, int *info,
               int strlen_jobz, int strlen_uplo );

  void dspev_( char *jobz, char *uplo, int *n, double *ap, double *w, double *z, int *ldz,
               double *work, int *info, int strlen_jobz, int strlen_uplo );

  void dgetrf_( int *m, int *n, double *a, int *lda, int *ipiv, int *info );

  void dgetri_( int *n, double *a, int *lda,
                int *ipiv, double *work, int *lwork, int *info );
}

namespace blocksqp {


  /**
   * Compute the inverse of a matrix
   * using LU decomposition (DGETRF and DGETRI)
   */
  int inverse( const Matrix &A, Matrix &Ainv ) {
    int i, j;
    int n, ldim, lwork, info = 0;
    int *ipiv;
    double *work;

    for (i=0; i<A.N(); i++ )
      for (j=0; j<A.M(); j++ )
        Ainv( j,i ) = A( j,i );

    n = Ainv.N();
    ldim = Ainv.LDIM();
    ipiv = new int[n];
    lwork = n*n;
    work = new double[lwork];

    // Compute LU factorization
    dgetrf_( &n, &n, Ainv.ARRAY(), &ldim, ipiv, &info );
    if ( info != 0 )
      printf( "WARNING: DGETRF returned info=%i\n", info );
    // Compute inverse
    dgetri_( &n, Ainv.ARRAY(), &ldim, ipiv, work, &lwork, &info );
    if ( info != 0 )
      printf( "WARNING: DGETRI returned info=%i\n", info );

    return info;
  }

  /**
   * Compute eigenvalues of a symmetric matrix by DSPEV
   */
  int calcEigenvalues( const SymMatrix &B, Matrix &ev ) {
    int n;
    SymMatrix temp;
    double *work, *dummy = 0;
    int info, iDummy = 1;

    n = B.M();
    ev.Dimension( n ).Initialize( 0.0 );
    work = new double[3*n];

    // copy Matrix, will be overwritten
    temp = SymMatrix( B );

    // DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
    // real symmetric matrix A in packed storage.
    dspev_( const_cast<char*>("N"), const_cast<char*>("L"), &n, temp.ARRAY(), ev.ARRAY(), dummy, &iDummy,
            work, &info, strlen("N"), strlen("L") );

    delete[] work;

    return info;
  }

  /**
   * Estimate the smalles eigenvalue of a sqare matrix
   * with the help og Gershgorin's circle theorem
   */
  double estimateSmallestEigenvalue( const Matrix &B )
  {
    int i, j;
    double radius;
    int dim = B.M();
    double lambdaMin = 0.0;

    // For each row, sum up off-diagonal elements
    for (i=0; i<dim; i++ )
      {
        radius = 0.0;
        for (j=0; j<dim; j++ )
          if (j != i )
            radius += fabs( B( i,j ) );

        if (B( i,i ) - radius < lambdaMin )
          lambdaMin = B( i,i ) - radius;
      }

    return lambdaMin;
  }


  /**
   * Compute scalar product aTb
   */
  double adotb( const Matrix &a, const Matrix &b ) {
    double norm = 0.0;

    if (a.N() != 1 || b.N() != 1 )
      {
        printf("a or b is not a vector!\n");
      }
    else if (a.M() != b.M() )
      {
        printf("a and b must have the same dimension!\n");
      }
    else
      {
        for (int k=0; k<a.M(); k++ )
          norm += a(k) * b(k);
      }

    return norm;
  }

  /**
   * Compute the matrix vector product for a column-compressed sparse matrix A with a vector b and store it in result
   */
  void Atimesb( double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result ) {
    int nCol = b.M();
    int nRow = result.M();
    int i, k;

    for (i=0; i<nRow; i++ )
      result( i ) = 0.0;

    for (i=0; i<nCol; i++ )
      {
        // k runs over all elements in one column
        for (k=AIndCol[i]; k<AIndCol[i+1]; k++ )
          result( AIndRow[k] ) += Anz[k] * b( i );
      }

  }

  /**
   * Compute the matrix vector product A*b and store it in result
   */
  void Atimesb( const Matrix &A, const Matrix &b, Matrix &result ) {
    result.Initialize( 0.0 );
    for (int i=0; i<A.M(); i++ )
      for (int k=0; k<A.N(); k++ )
        result( i ) += A( i, k ) * b( k );
  }

  double l1VectorNorm( const Matrix &v ) {
    double norm = 0.0;

    if (v.N() != 1 )
      {
        printf("v is not a vector!\n");
      }
    else
      {
        for (int k=0; k<v.M(); k++ )
          norm += fabs(v( k ));
      }

    return norm;
  }

  double l2VectorNorm( const Matrix &v ) {
    double norm = 0.0;

    if (v.N() != 1 )
      {
        printf("v is not a vector!\n");
      }
    else
      {
        for (int k=0; k<v.M(); k++ )
          norm += v( k )* v( k );
      }

    return sqrt(norm);
  }

  double lInfVectorNorm( const Matrix &v ) {
    double norm = 0.0;

    if (v.N() != 1 )
      {
        printf("v is not a vector!\n");
      }
    else
      {
        for (int k=0; k<v.M(); k++ )
          if (fabs(v( k )) > norm )
            norm = fabs(v( k ));
      }

    return norm;
  }


  /**
   * Calculate weighted l1 norm of constraint violations
   */
  double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl, const Matrix &weights ) {
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    if (weights.M() < constr.M() + nVar )
      {
        printf("Weight vector too short!\n");
        return 0.0;
      }

    // Weighted violation of simple bounds
    for (i=0; i<nVar; i++ )
      {
        if (xi( i ) > bu( i ) )
          norm += fabs(weights( i )) * (xi( i ) - bu( i ));
        else if (xi( i ) < bl( i ) )
          norm += fabs(weights( i )) * (bl( i ) - xi( i ));
      }

    // Calculate weighted sum of constraint violations
    for (i=0; i<constr.M(); i++ )
      {
        if (constr( i ) > bu( nVar+i ) )
          norm += fabs(weights( nVar+i )) * (constr( i ) - bu( nVar+i ));
        else if (constr( i ) < bl( nVar+i ) )
          norm += fabs(weights( nVar+i )) * (bl( nVar+i ) - constr( i ));
      }

    return norm;
  }


  /**
   * Calculate l1 norm of constraint violations
   */
  double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl ) {
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for (i=0; i<nVar; i++ )
      {
        if (xi( i ) > bu( i ) )
          norm += xi( i ) - bu( i );
        else if (xi( i ) < bl( i ) )
          norm += bl( i ) - xi( i );
      }

    // Calculate sum of constraint violations
    for (i=0; i<constr.M(); i++ )
      {
        if (constr( i ) > bu( nVar+i ) )
          norm += constr( i ) - bu( nVar+i );
        else if (constr( i ) < bl( nVar+i ) )
          norm += bl( nVar+i ) - constr( i );
      }

    return norm;
  }


  /**
   * Calculate l2 norm of constraint violations
   */
  double l2ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl ) {
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for (i=0; i<nVar; i++ )
      if (xi( i ) > bu( i ) )
        norm += xi( i ) - bu( i );
    if (xi( i ) < bl( i ) )
      norm += bl( i ) - xi( i );

    // Calculate sum of constraint violations
    for (i=0; i<constr.M(); i++ )
      if (constr( i ) > bu( nVar+i ) )
        norm += pow(constr( i ) - bu( nVar+i ), 2);
      else if (constr( i ) < bl( nVar+i ) )
        norm += pow(bl( nVar+i ) - constr( i ), 2);

    return sqrt(norm);
  }


  /**
   * Calculate l_Infinity norm of constraint violations
   */
  double lInfConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl ) {
    double norm = 0.0;
    int i;
    int nVar = xi.M();
    int nCon = constr.M();

    // Violation of simple bounds
    for (i=0; i<nVar; i++ ) {
      if (xi( i ) - bu( i ) > norm )
        norm = xi( i ) - bu( i );
      else if (bl( i ) - xi( i ) > norm )
        norm = bl( i ) - xi( i );
    }

    // Find out the largest constraint violation
    for (i=0; i<nCon; i++ ) {
      if (constr( i ) - bu( nVar+i ) > norm )
        norm = constr( i ) - bu( nVar+i );
      if (bl( nVar+i ) - constr( i ) > norm )
        norm = bl( nVar+i ) - constr( i );
    }

    return norm;
  }

  void Error( const char *F ) {
    printf("Error: %s\n", F );
    //exit( 1 );
  }

  /* ----------------------------------------------------------------------- */

  int Matrix::malloc( void ) {
    int len;

    if ( tflag )
      Error("malloc cannot be called with Submatrix");

    if ( ldim < m )
      ldim = m;

    len = ldim*n;

    if ( len == 0 )
      array = 0;
    else
      if ( ( array = new double[len] ) == 0 )
        Error("'new' failed");

    return 0;
  }


  int Matrix::free( void ) {
    if ( tflag )
      Error("free cannot be called with Submatrix");

    if ( array != 0 )
      delete[] array;

    return 0;
  }


  double &Matrix::operator()( int i, int j ) {
    return array[i+j*ldim];
  }

  double &Matrix::operator()( int i, int j ) const {
    return array[i+j*ldim];
  }

  double &Matrix::operator()( int i ) {
    return array[i];
  }

  double &Matrix::operator()( int i ) const {
    return array[i];
  }

  Matrix::Matrix( int M, int N, int LDIM ) {
    m = M;
    n = N;
    ldim = LDIM;
    tflag = 0;

    malloc();
  }


  Matrix::Matrix( int M, int N, double *ARRAY, int LDIM ) {
    m = M;
    n = N;
    array = ARRAY;
    ldim = LDIM;
    tflag = 0;

    if ( ldim < m )
      ldim = m;
  }


  Matrix::Matrix( const Matrix &A ) {
    int i, j;

    m = A.m;
    n = A.n;
    ldim = A.ldim;
    tflag = 0;

    malloc();

    for ( i = 0; i < m; i++ )
      for ( j = 0; j < n ; j++ )
        (*this)(i,j) = A(i,j);
    //(*this)(i,j) = A.a(i,j);
  }

  Matrix &Matrix::operator=( const Matrix &A ) {
    int i, j;

    if ( this != &A )
      {
        if ( !tflag )
          {
            free();

            m = A.m;
            n = A.n;
            ldim = A.ldim;

            malloc();

            for ( i = 0; i < m; i++ )
              for ( j = 0; j < n ; j++ )
                (*this)(i,j) = A(i,j);
          }
        else
          {
            if ( m != A.m || n != A.n )
              Error("= operation not allowed");

            for ( i = 0; i < m; i++ )
              for ( j = 0; j < n ; j++ )
                (*this)(i,j) = A(i,j);
          }
      }

    return *this;
  }


  Matrix::~Matrix( void ) {
    if ( !tflag )
      free();
  }

  /* ----------------------------------------------------------------------- */

  int Matrix::M( void ) const {
    return m;
  }


  int Matrix::N( void ) const {
    return n;
  }


  int Matrix::LDIM( void ) const {
    return ldim;
  }


  double *Matrix::ARRAY( void ) const {
    return array;
  }


  int Matrix::TFLAG( void ) const {
    return tflag;
  }

  /* ----------------------------------------------------------------------- */

  Matrix &Matrix::Dimension( int M, int N, int LDIM ) {
    if ( M != m || N != n || ( LDIM != ldim && LDIM != -1 ) )
      {
        if ( tflag )
          Error("Cannot set new dimension for Submatrix");
        else
          {
            free();
            m = M;
            n = N;
            ldim = LDIM;

            malloc();
          }
      }

    return *this;
  }

  Matrix &Matrix::Initialize( double (*f)( int, int ) ) {
    int i, j;

    for ( i = 0; i < m; i++ )
      for ( j = 0; j < n; j++ )
        (*this)(i,j) = f(i,j);

    return *this;
  }


  Matrix &Matrix::Initialize( double val ) {
    int i, j;

    for ( i = 0; i < m; i++ )
      for ( j = 0; j < n; j++ )
        (*this)(i,j) = val;

    return *this;
  }


  /* ----------------------------------------------------------------------- */

  Matrix &Matrix::Submatrix( const Matrix &A, int M, int N, int i0, int j0 ) {
    if ( i0 + M > A.m || j0 + N > A.n )
      Error("Cannot create Submatrix");

    if ( !tflag )
      free();

    tflag = 1;

    m = M;
    n = N;
    array = &A.array[i0+j0*A.ldim];
    ldim = A.ldim;

    return *this;
  }


  Matrix &Matrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM ) {
    if ( !tflag )
      free();

    tflag = 1;

    m = M;
    n = N;
    array = ARRAY;
    ldim = LDIM;

    if ( ldim < m )
      ldim = m;

    return *this;
  }


  const Matrix &Matrix::Print( FILE *f, int DIGITS, int flag ) const {
    int i, j;
    double x;

    // Flag == 1: Matlab output
    // else: plain output

    if ( flag == 1 )
      fprintf( f, "[" );

    for ( i = 0; i < m; i++ )
      {
        for ( j = 0; j < n; j++ )
          {
            x = (*this)(i,j);
            //x = a(i,j);

            if ( flag == 1 )
              {
                fprintf( f, j == 0 ? " " : ", " );
                fprintf( f, "%.*le", DIGITS, x );
              }
            else
              {
                fprintf( f, j == 0 ? "" : "  " );
                fprintf( f, "% .*le", DIGITS, x );
              }
          }
        if ( flag == 1 )
          {
            if ( i < m-1 )
              fprintf( f, ";\n" );
          }
        else
          {
            if ( i < m-1 )
              fprintf( f, "\n" );
          }
      }

    if ( flag == 1 )
      fprintf( f, " ];\n" );
    else
      fprintf( f, "\n" );

    return *this;
  }


  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */



  int SymMatrix::malloc( void ) {
    int len;

    len = m*(m+1)/2.0;

    if ( len == 0 )
      array = 0;
    else
      if ( ( array = new double[len] ) == 0 )
        Error("'new' failed");

    return 0;
  }


  int SymMatrix::free( void ) {
    if (array != 0 )
      delete[] array;

    return 0;
  }


  double &SymMatrix::operator()( int i, int j ) {
    int pos;

    if (i < j )//reference to upper triangular part
      pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
      pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
  }


  double &SymMatrix::operator()( int i, int j ) const {
    int pos;

    if (i < j )//reference to upper triangular part
      pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
      pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
  }


  double &SymMatrix::operator()( int i ) {
    return array[i];
  }


  double &SymMatrix::operator()( int i ) const {
    return array[i];
  }

  SymMatrix::SymMatrix( int M ) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
  }

  SymMatrix::SymMatrix( int M, double *ARRAY ) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
  }


  SymMatrix::SymMatrix( int M, int N, int LDIM ) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
  }


  SymMatrix::SymMatrix( int M, int N, double *ARRAY, int LDIM ) {
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
  }


  SymMatrix::SymMatrix( const Matrix &A ) {
    int i, j;

    m = A.M();
    n = A.M();
    ldim = A.M();
    tflag = 0;

    malloc();

    for ( j=0; j<m; j++ )//columns
      for ( i=j; i<m; i++ )//rows
        (*this)(i,j) = A(i,j);
    //(*this)(i,j) = A.a(i,j);
  }


  SymMatrix::SymMatrix( const SymMatrix &A ) {
    int i, j;

    m = A.m;
    n = A.n;
    ldim = A.ldim;
    tflag = 0;

    malloc();

    for ( j=0; j<m; j++ )//columns
      for ( i=j; i<m; i++ )//rows
        (*this)(i,j) = A(i,j);
    //(*this)(i,j) = A.a(i,j);
  }


  SymMatrix::~SymMatrix( void ) {
    if (!tflag )
      free();
  }



  SymMatrix &SymMatrix::Dimension( int M ) {
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
  }


  SymMatrix &SymMatrix::Dimension( int M, int N, int LDIM ) {
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
  }


  SymMatrix &SymMatrix::Initialize( double (*f)( int, int ) ) {
    int i, j;

    for ( j=0; j<m; j++ )
      for ( i=j; i<n ; i++ )
        (*this)(i,j) = f(i,j);

    return *this;
  }


  SymMatrix &SymMatrix::Initialize( double val ) {
    int i, j;

    for ( j=0; j<m; j++ )
      for ( i=j; i<n ; i++ )
        (*this)(i,j) = val;

    return *this;
  }


  SymMatrix &SymMatrix::Submatrix( const Matrix &A, int M, int N, int i0, int j0) {
    Error("SymMatrix doesn't support Submatrix");
    return *this;
  }


  SymMatrix &SymMatrix::Arraymatrix( int M, double *ARRAY ) {
    if (!tflag )
      free();

    tflag = 1;
    m = M;
    n = M;
    ldim = M;
    array = ARRAY;

    return *this;
  }


  SymMatrix &SymMatrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM ) {
    if (!tflag )
      free();

    tflag = 1;
    m = M;
    n = M;
    ldim = M;
    array = ARRAY;

    return *this;
  }


  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */


  double delta( int i, int j ) {
    return (i == j) ? 1.0 : 0.0;
  }


  Matrix Transpose( const Matrix &A ) {
    int i, j;
    double *array;

    if ( ( array = new double[A.N()*A.M()] ) == 0 )
      Error("'new' failed");

    for ( i = 0; i < A.N(); i++ )
      for ( j = 0; j < A.M(); j++ )
        array[i+j*A.N()] = A(j,i);
    //array[i+j*A.N()] = A.a(j,i);

    return Matrix( A.N(), A.M(), array, A.N() );
  }


  Matrix &Transpose( const Matrix &A, Matrix &T ) {
    int i, j;

    T.Dimension( A.N(), A.M() );

    for ( i = 0; i < A.N(); i++ )
      for ( j = 0; j < A.M(); j++ )
        T(i,j) = A(j,i);
    //T(i,j) = A.a(j,i);

    return T;
  }

  void Problemspec::evaluate( const Matrix &xi, double *objval, Matrix &constr, int *info ) {
    Matrix lambdaDummy, gradObjDummy;
    SymMatrix *hessDummy;
    int dmode = 0;

    Matrix constrJacDummy;
    double *jacNzDummy;
    int *jacIndRowDummy, *jacIndColDummy;
    *info = 0;

    // Try sparse version first
    evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, jacNzDummy, jacIndRowDummy, jacIndColDummy, hessDummy, dmode, info );

    // If sparse version is not implemented, try dense version
    if (info )
      evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, constrJacDummy, hessDummy, dmode, info );
  }

} // namespace blocksqp
