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

  /**
   * Standard Constructor:
   * Default settings
   */
  SQPoptions::SQPoptions() {
    /* qpOASES: dense (0), sparse (1), or Schur (2)
     * Choice of qpOASES method:
     * 0: dense Hessian and Jacobian, dense factorization of reduced Hessian
     * 1: sparse Hessian and Jacobian, dense factorization of reduced Hessian
     * 2: sparse Hessian and Jacobian, Schur complement approach (recommended) */
    sparseQP = 2;

    // 0: no output, 1: normal output, 2: verbose output
    printLevel = 2;

    /* 0: no debug output, 1: print one line per iteration to file,
       2: extensive debug output to files (impairs performance) */
    debugLevel = 0;

    //eps = 2.2204e-16;
    eps = 1.0e-16;
    opttol = 1.0e-6;
    nlinfeastol = 1.0e-6;

    // 0: no globalization, 1: filter line search
    globalization = 1;

    // 0: no feasibility restoration phase 1: if line search fails, start feasibility restoration phase
    restoreFeas = 1;

    // 0: globalization is always active, 1: take a full step at first SQP iteration, no matter what
    skipFirstGlobalization = false;

    // 0: one update for large Hessian, 1: apply updates blockwise, 2: 2 blocks: 1 block updates, 1 block Hessian of obj.
    blockHess = 1;

    // after too many consecutive skipped updates, Hessian block is reset to (scaled) identity
    maxConsecSkippedUpdates = 100;

    // for which blocks should second derivatives be provided by the user:
    // 0: none, 1: for the last block, 2: for all blocks
    whichSecondDerv = 0;

    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: geometric mean of 1 and 2
    // 4: centered Oren-Luenberger sizing according to Tapia paper
    hessScaling = 2;
    fallbackScaling = 4;
    iniHessDiag = 1.0;

    // Activate damping strategy for BFGS (if deactivated, BFGS might yield indefinite updates!)
    hessDamp = 1;

    // Damping factor for Powell modification of BFGS updates ( between 0.0 and 1.0 )
    hessDampFac = 0.2;

    // 0: constant, 1: SR1, 2: BFGS (damped), 3: [not used] , 4: finiteDiff, 5: Gauss-Newton
    hessUpdate = 1;
    fallbackUpdate = 2;

    //
    convStrategy = 0;

    // How many ADDITIONAL (convexified) QPs may be solved per iteration?
    maxConvQP = 1;

    // 0: full memory updates 1: limited memory
    hessLimMem = 1;

    // memory size for L-BFGS/L-SR1 updates
    hessMemsize = 20;

    // maximum number of line search iterations
    maxLineSearch = 20;

    // if step has to be reduced in too many consecutive iterations, feasibility restoration phase is invoked
    maxConsecReducedSteps = 100;

    // maximum number of second-order correction steps
    maxSOCiter = 3;

    // maximum number of QP iterations per QP solve
    maxItQP = 5000;
    // maximum time (in seconds) for one QP solve
    maxTimeQP = 10000.0;

    // Oren-Luenberger scaling parameters
    colEps = 0.1;
    colTau1 = 0.5;
    colTau2 = 1.0e4;

    // Filter line search parameters
    gammaTheta = 1.0e-5;
    gammaF = 1.0e-5;
    kappaSOC = 0.99;
    kappaF = 0.999;
    thetaMax = 1.0e7;       // reject steps if constr viol. is larger than thetaMax
    thetaMin = 1.0e-5;      // if constr viol. is smaller than thetaMin require Armijo cond. for obj.
    delta = 1.0;
    sTheta = 1.1;
    sF = 2.3;
    eta = 1.0e-4;

    // Inertia correction for filter line search and indefinite Hessians
    kappaMinus = 0.333;
    kappaPlus = 8.0;
    kappaPlusMax = 100.0;
    deltaH0 = 1.0e-4;
  }


  /**
   * Some options cannot be set together, resolve here
   */
  void SQPoptions::optionsConsistency() {
    // If we compute second constraints derivatives switch to finite differences Hessian (convenience)
    if (whichSecondDerv == 2 )
      {
        hessUpdate = 4;
        blockHess = 1;
      }

    // If we don't use limited memory BFGS we need to store only one vector.
    if (!hessLimMem )
      hessMemsize = 1;

    if (sparseQP != 2 && hessUpdate == 1 )
      {
        printf( "SR1 update only works with qpOASES Schur complement version. Using BFGS updates instead.\n" );
        hessUpdate = 2;
        hessScaling = fallbackScaling;
      }
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
