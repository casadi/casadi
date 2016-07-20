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

  SQPiterate::SQPiterate( Problemspec* prob, SQPoptions* param, bool full ) {
    int maxblocksize = 1;

    // Set nBlocks structure according to if we use block updates or not
    if (param->blockHess == 0 || prob->nBlocks == 1 ) {
      nBlocks = 1;
      blockIdx = new int[2];
      blockIdx[0] = 0;
      blockIdx[1] = prob->nVar;
      maxblocksize = prob->nVar;
      param->whichSecondDerv = 0;
    } else if (param->blockHess == 2 && prob->nBlocks > 1 ) {
      // hybrid strategy: 1 block for constraints, 1 for objective
      nBlocks = 2;
      blockIdx = new int[3];
      blockIdx[0] = 0;
      blockIdx[1] = prob->blockIdx[prob->nBlocks-1];
      blockIdx[2] = prob->nVar;
    } else {
      nBlocks = prob->nBlocks;
      blockIdx = new int[nBlocks+1];
      for (int k=0; k<nBlocks+1; k++ ) {
        blockIdx[k] = prob->blockIdx[k];
        if (k > 0 )
          if (blockIdx[k] - blockIdx[k-1] > maxblocksize )
            maxblocksize = blockIdx[k] - blockIdx[k-1];
      }
    }

    if (param->hessLimMem && param->hessMemsize == 0 )
      param->hessMemsize = maxblocksize;

    allocMin( prob );

    if (!param->sparseQP ) {
      constrJac.Dimension( prob->nCon, prob->nVar ).Initialize( 0.0 );
      hessNz = new double[prob->nVar*prob->nVar];
    } else {
      hessNz = 0;
    }

    jacNz = 0;
    jacIndCol = 0;
    jacIndRow = 0;

    hessIndCol = 0;
    hessIndRow = 0;
    hessIndLo = 0;
    hess = 0;
    hess1 = 0;
    hess2 = 0;

    noUpdateCounter = 0;

    if (full ) {
      allocHess( param );
      allocAlg( prob, param );
    }
  }


  SQPiterate::SQPiterate( const SQPiterate &iter ) {
    int i;

    nBlocks = iter.nBlocks;
    blockIdx = new int[nBlocks+1];
    for (i=0; i<nBlocks+1; i++ )
      blockIdx[i] = iter.blockIdx[i];

    xi = iter.xi;
    lambda = iter.lambda;
    constr = iter.constr;
    gradObj = iter.gradObj;
    gradLagrange = iter.gradLagrange;

    constrJac = iter.constrJac;
    if (iter.jacNz != 0 ) {
      int nVar = xi.M();
      int nnz = iter.jacIndCol[nVar];

      jacNz = new double[nnz];
      for (i=0; i<nnz; i++ )
        jacNz[i] = iter.jacNz[i];

      jacIndRow = new int[nnz + (nVar+1) + nVar];
      for (i=0; i<nnz + (nVar+1) + nVar; i++ )
        jacIndRow[i] = iter.jacIndRow[i];
      jacIndCol = jacIndRow + nnz;
    } else {
      jacNz = 0;
      jacIndRow = 0;
      jacIndCol = 0;
    }

    noUpdateCounter = 0;
    hessNz = 0;
    hessIndCol = 0;
    hessIndRow = 0;
    hessIndLo = 0;
    hess = 0;
    hess1 = 0;
    hess2 = 0;
  }


  /**
   * Allocate memory for variables
   * required by all optimization
   * algorithms except for the Jacobian
   */
  void SQPiterate::allocMin( Problemspec *prob ) {
    // current iterate
    xi.Dimension( prob->nVar ).Initialize( 0.0 );

    // dual variables (for general constraints and variable bounds)
    lambda.Dimension( prob->nVar + prob->nCon ).Initialize( 0.0 );

    // constraint vector with lower and upper bounds
    // (Box constraints are not included in the constraint list)
    constr.Dimension( prob->nCon ).Initialize( 0.0 );

    // gradient of objective
    gradObj.Dimension( prob->nVar ).Initialize( 0.0 );

    // gradient of Lagrangian
    gradLagrange.Dimension( prob->nVar ).Initialize( 0.0 );
  }


  void SQPiterate::allocHess( SQPoptions *param ) {
    int iBlock, varDim;

    // Create one Matrix for one diagonal block in the Hessian
    hess1 = new SymMatrix[nBlocks];
    for (iBlock=0; iBlock<nBlocks; iBlock++ )
      {
        varDim = blockIdx[iBlock+1] - blockIdx[iBlock];
        hess1[iBlock].Dimension( varDim ).Initialize( 0.0 );
      }

    // For SR1 or finite differences, maintain two Hessians
    if (param->hessUpdate == 1 || param->hessUpdate == 4 ) {
      hess2 = new SymMatrix[nBlocks];
      for (iBlock=0; iBlock<nBlocks; iBlock++ ) {
        varDim = blockIdx[iBlock+1] - blockIdx[iBlock];
        hess2[iBlock].Dimension( varDim ).Initialize( 0.0 );
      }
    }

    // Set Hessian pointer
    hess = hess1;
  }

  /**
   * Convert diagonal block Hessian to double array.
   * Assumes that hessNz is already allocated.
   */
  void SQPiterate::convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_ ) {
    if (hessNz == 0 ) return;
    int count = 0;
    int blockCnt = 0;
    for (int i=0; i<prob->nVar; i++ )
      for (int j=0; j<prob->nVar; j++ )
        {
          if (i == blockIdx[blockCnt+1] )
            blockCnt++;
          if (j >= blockIdx[blockCnt] && j < blockIdx[blockCnt+1] )
            hessNz[count++] = hess[blockCnt]( i - blockIdx[blockCnt], j - blockIdx[blockCnt] );
          else
            hessNz[count++] = 0.0;
        }
  }

  /**
   * Convert array *hess to a single symmetric sparse matrix in
   * Harwell-Boeing format (as used by qpOASES)
   */
  void SQPiterate::convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_,
                                   double *&hessNz_, int *&hessIndRow_, int *&hessIndCol_, int *&hessIndLo_ ) {
    int iBlock, count, colCountTotal, rowOffset, i, j;
    int nnz, nCols, nRows;

    // 1) count nonzero elements
    nnz = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++ )
      for (i=0; i<hess_[iBlock].N(); i++ )
        for (j=i; j<hess_[iBlock].N(); j++ )
          if (fabs(hess_[iBlock]( i,j )) > eps ) {
            nnz++;
            if (i != j ) {
              // off-diagonal elements count twice
              nnz++;
            }
          }

    if (hessNz_ != 0 ) delete[] hessNz_;
    if (hessIndRow_ != 0 ) delete[] hessIndRow_;

    hessNz_ = new double[nnz];
    hessIndRow_ = new int[nnz + (prob->nVar+1) + prob->nVar];
    hessIndCol_ = hessIndRow_ + nnz;
    hessIndLo_ = hessIndCol_ + (prob->nVar+1);

    // 2) store matrix entries columnwise in hessNz
    count = 0; // runs over all nonzero elements
    colCountTotal = 0; // keep track of position in large matrix
    rowOffset = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++ ) {
      nCols = hess_[iBlock].N();
      nRows = hess_[iBlock].M();

      for (i=0; i<nCols; i++ ) {
        // column 'colCountTotal' starts at element 'count'
        hessIndCol_[colCountTotal] = count;

        for (j=0; j<nRows; j++ )
          if (fabs(hess_[iBlock]( i,j )) > eps )
            {
              hessNz_[count] = hess_[iBlock]( i, j );
              hessIndRow_[count] = j + rowOffset;
              count++;
            }
        colCountTotal++;
      }

      rowOffset += nRows;
    }
    hessIndCol_[colCountTotal] = count;

    // 3) Set reference to lower triangular matrix
    for (j=0; j<prob->nVar; j++ ) {
      for (i=hessIndCol_[j]; i<hessIndCol_[j+1] && hessIndRow_[i]<j; i++);
      hessIndLo_[j] = i;
    }

    if (count != nnz )
      printf( "Error in convertHessian: %i elements processed, should be %i elements!\n", count, nnz );
  }


  /**
   * Allocate memory for additional variables
   * needed by the algorithm
   */
  void SQPiterate::allocAlg( Problemspec *prob, SQPoptions *param ) {
    int iBlock;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // current step
    deltaMat.Dimension( nVar, param->hessMemsize, nVar ).Initialize( 0.0 );
    deltaXi.Submatrix( deltaMat, nVar, 1, 0, 0 );
    // trial step (temporary variable, for line search)
    trialXi.Dimension( nVar, 1, nVar ).Initialize( 0.0 );

    // bounds for step (QP subproblem)
    deltaBl.Dimension( nVar+nCon ).Initialize( 0.0 );
    deltaBu.Dimension( nVar+nCon ).Initialize( 0.0 );

    // product of constraint Jacobian with step (deltaXi)
    AdeltaXi.Dimension( nCon ).Initialize( 0.0 );

    // dual variables of QP (simple bounds and general constraints)
    lambdaQP.Dimension( nVar+nCon ).Initialize( 0.0 );

    // line search parameters
    deltaH.Dimension( nBlocks ).Initialize( 0.0 );

    // filter as a set of pairs
    filter = new std::set< std::pair<double,double> >;

    // difference of Lagrangian gradients
    gammaMat.Dimension( nVar, param->hessMemsize, nVar ).Initialize( 0.0 );
    gamma.Submatrix( gammaMat, nVar, 1, 0, 0 );

    // Scalars that are used in various Hessian update procedures
    noUpdateCounter = new int[nBlocks];
    for (iBlock=0; iBlock<nBlocks; iBlock++ )
      noUpdateCounter[iBlock] = -1;

    // For selective sizing: for each block save sTs, sTs_, sTy, sTy_
    deltaNorm.Dimension( nBlocks ).Initialize( 1.0 );
    deltaNormOld.Dimension( nBlocks ).Initialize( 1.0 );
    deltaGamma.Dimension( nBlocks ).Initialize( 0.0 );
    deltaGammaOld.Dimension( nBlocks ).Initialize( 0.0 );
  }


  void SQPiterate::initIterate( SQPoptions* param ) {
    alpha = 1.0;
    nSOCS = 0;
    reducedStepCount = 0;
    steptype = 0;

    obj = param->inf;
    tol = param->inf;
    cNorm = param->thetaMax;
    gradNorm = param->inf;
    lambdaStepNorm = 0.0;
  }

  SQPiterate::~SQPiterate( void ) {
    if (blockIdx != 0 )
      delete[] blockIdx;
    if (noUpdateCounter != 0 )
      delete[] noUpdateCounter;
    if (jacNz != 0 )
      delete[] jacNz;
    if (jacIndRow != 0 )
      delete[] jacIndRow;
    if (hessNz != 0 )
      delete[] hessNz;
    if (hessIndRow != 0 )
      delete[] hessIndRow;
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
    // 1: (some) colorful output
    printColor = 1;

    /* 0: no debug output, 1: print one line per iteration to file,
       2: extensive debug output to files (impairs performance) */
    debugLevel = 0;

    //eps = 2.2204e-16;
    eps = 1.0e-16;
    inf = 1.0e20;
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
