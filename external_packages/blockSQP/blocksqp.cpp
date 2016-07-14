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

  Blocksqp::Blocksqp( Problemspec *problem, SQPoptions *parameters, SQPstats *statistics ) {
    prob = problem;
    param = parameters;
    stats = statistics;

    // Check if there are options that are infeasible and set defaults accordingly
    param->optionsConsistency();

    vars = new SQPiterate( prob, param, 1 );

    if ( param->sparseQP < 2 ) {
      qp = new qpOASES::SQProblem( prob->nVar, prob->nCon );
      qpSave = new qpOASES::SQProblem( prob->nVar, prob->nCon );
    } else {
      qp = new qpOASES::SQProblemSchur( prob->nVar, prob->nCon, qpOASES::HST_UNKNOWN, 50 );
      qpSave = new qpOASES::SQProblemSchur( prob->nVar, prob->nCon, qpOASES::HST_UNKNOWN, 50 );
    }

    initCalled = false;
  }

  Blocksqp::~Blocksqp() {
    delete qp;
    delete qpSave;
    delete vars;
  }


  void Blocksqp::init() {
    // Print header and information about the algorithmic parameters
    printInfo( param->printLevel );

    // Open output files
    stats->initStats( param );
    vars->initIterate( param );

    // Initialize filter with pair ( maxConstrViolation, objLowerBound )
    initializeFilter();

    // Set initial values for all xi and set the Jacobian for linear constraints
    if (param->sparseQP )
      prob->initialize( vars->xi, vars->lambda, vars->jacNz, vars->jacIndRow, vars->jacIndCol );
    else
      prob->initialize( vars->xi, vars->lambda, vars->constrJac );

    initCalled = true;
  }


  int Blocksqp::run( int maxIt, int warmStart ) {
    int it, infoQP = 0, infoEval = 0;
    bool skipLineSearch = false;
    bool hasConverged = false;
    int whichDerv = param->whichSecondDerv;

    if (!initCalled ) {
      printf("init() must be called before run(). Aborting.\n");
      return -1;
    }

    if (warmStart == 0 || stats->itCount == 0 ) {
      // SQP iteration 0

      /// Set initial Hessian approximation
      calcInitialHessian();

      /// Evaluate all functions and gradients for xi_0
      if (param->sparseQP )
        prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                        vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->hess, 1+whichDerv, &infoEval );
      else
        prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                        vars->constrJac, vars->hess, 1+whichDerv, &infoEval );
      stats->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol();
      stats->printProgress( prob, vars, param, hasConverged );
      if (hasConverged )
        return 0;

      stats->itCount++;
    }

    /*
     * SQP Loop: during first iteration, stats->itCount = 1
     */
    for (it=0; it<maxIt; it++ ) {
      /// Solve QP subproblem with qpOASES or QPOPT
      updateStepBounds( 0 );
      infoQP = solveQP( vars->deltaXi, vars->lambdaQP );

      if (infoQP == 1 )
        {// 1.) Maximum number of iterations reached
          printf("***Warning! Maximum number of QP iterations exceeded.***\n");
          ;// just continue ...
        }
      else if (infoQP == 2 || infoQP > 3 )
        {// 2.) QP error (e.g., unbounded), solve again with pos.def. diagonal matrix (identity)
          printf("***QP error. Solve again with identity matrix.***\n");
          resetHessian();
          infoQP = solveQP( vars->deltaXi, vars->lambdaQP );
          if (infoQP ) {
            // If there is still an error, terminate.
            printf( "***QP error. Stop.***\n" );
            return -1;
          } else {
            vars->steptype = 1;
          }
        } else if (infoQP == 3 ) {
        // 3.) QP infeasible, try to restore feasibility
        bool qpError = true;
        skipLineSearch = true; // don't do line search with restoration step

        // Try to reduce constraint violation by heuristic
        if (vars->steptype < 2 ) {
          printf("***QP infeasible. Trying to reduce constraint violation...");
          qpError = feasibilityRestorationHeuristic();
          if (!qpError ) {
            vars->steptype = 2;
            printf("Success.***\n");
          } else {
            printf("Failed.***\n");
          }
        }

        // Invoke feasibility restoration phase
        //if (qpError && vars->steptype < 3 && param->restoreFeas )
        if (qpError && param->restoreFeas && vars->cNorm > 0.01 * param->nlinfeastol ) {
          printf("***Start feasibility restoration phase.***\n");
          vars->steptype = 3;
          qpError = feasibilityRestorationPhase();
        }

        // If everything failed, abort.
        if (qpError ) {
          printf( "***QP error. Stop.***\n" );
          return -1;
        }
      }

      /// Determine steplength alpha
      if (param->globalization == 0 || (param->skipFirstGlobalization && stats->itCount == 1) ) {
        // No globalization strategy, but reduce step if function cannot be evaluated
        if (fullstep()) {
          printf( "***Constraint or objective could not be evaluated at new point. Stop.***\n" );
          return -1;
        }
        vars->steptype = 0;
      } else if (param->globalization == 1 && !skipLineSearch ) {
        // Filter line search based on Waechter et al., 2006 (Ipopt paper)
        if (filterLineSearch() || vars->reducedStepCount > param->maxConsecReducedSteps ) {
          // Filter line search did not produce a step. Now there are a few things we can try ...
          bool lsError = true;

          // Heuristic 1: Check if the full step reduces the KKT error by at least kappaF, if so, accept the step.
          lsError = kktErrorReduction( );
          if (!lsError )
            vars->steptype = -1;

          // Heuristic 2: Try to reduce constraint violation by closing continuity gaps to produce an admissable iterate
          if (lsError && vars->cNorm > 0.01 * param->nlinfeastol && vars->steptype < 2 ) {
            // Don't do this twice in a row!

            printf("***Warning! Steplength too short. Trying to reduce constraint violation...");

            // Integration over whole time interval
            lsError = feasibilityRestorationHeuristic( );
            if (!lsError )
              {
                vars->steptype = 2;
                printf("Success.***\n");
              }
            else
              printf("Failed.***\n");
          }

          // Heuristic 3: Recompute step with a diagonal Hessian
          if (lsError && vars->steptype != 1 && vars->steptype != 2 ) {
            // After closing continuity gaps, we already take a step with initial Hessian. If this step is not accepted then this will cause an infinite loop!

            printf("***Warning! Steplength too short. Trying to find a new step with identity Hessian.***\n");
            vars->steptype = 1;

            resetHessian();
            continue;
          }

          // If this does not yield a successful step, start restoration phase
          if (lsError && vars->cNorm > 0.01 * param->nlinfeastol && param->restoreFeas ) {
            printf("***Warning! Steplength too short. Start feasibility restoration phase.***\n");
            vars->steptype = 3;

            // Solve NLP with minimum norm objective
            lsError = feasibilityRestorationPhase( );
          }

          // If everything failed, abort.
          if (lsError ) {
            printf( "***Line search error. Stop.***\n" );
            return -1;
          }
        } else {
          vars->steptype = 0;
        }
      }

      /// Calculate "old" Lagrange gradient: gamma = dL(xi_k, lambda_k+1)
      calcLagrangeGradient( vars->gamma, 0 );

      /// Evaluate functions and gradients at the new xi
      if (param->sparseQP ) {
        prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                        vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->hess, 1+whichDerv, &infoEval );
      } else {
        prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                        vars->constrJac, vars->hess, 1+whichDerv, &infoEval );
      }
      stats->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol();

      /// Print one line of output for the current iteration
      stats->printProgress( prob, vars, param, hasConverged );
      if (hasConverged && vars->steptype < 2 ) {
        stats->itCount++;
        if (param->debugLevel > 2 ) {
          //printf("Computing finite differences Hessian at the solution ... \n");
          //calcFiniteDiffHessian( );
          //stats->printHessian( prob->nBlocks, vars->hess );
          stats->dumpQPCpp( prob, vars, qp, param->sparseQP );
        }
        return 0; //Convergence achieved!
      }

      /// Calculate difference of old and new Lagrange gradient: gamma = -gamma + dL(xi_k+1, lambda_k+1)
      calcLagrangeGradient( vars->gamma, 1 );

      /// Revise Hessian approximation
      if (param->hessUpdate < 4 && !param->hessLimMem ) {
        calcHessianUpdate( param->hessUpdate, param->hessScaling );
      } else if (param->hessUpdate < 4 && param->hessLimMem ) {
        calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling );
      } else if (param->hessUpdate == 4 ) {
        calcFiniteDiffHessian( );
      }

      // If limited memory updates  are used, set pointers deltaXi and gamma to the next column in deltaMat and gammaMat
      updateDeltaGamma();

      stats->itCount++;
      skipLineSearch = false;
    }

    return 1;
  }


  void Blocksqp::finish() {
    if (initCalled ) {
      initCalled = false;
    } else {
      printf("init() must be called before finish().\n");
      return;
    }

    stats->finish( param );
  }


  /**
   * Compute gradient of Lagrangian or difference of Lagrangian gradients (sparse version)
   *
   * flag == 0: output dL(xi, lambda)
   * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
   * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
   */
  void Blocksqp::calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, double *jacNz, int *jacIndRow, int *jacIndCol,
                                       Matrix &gradLagrange, int flag ) {
    int iVar, iCon;

    // Objective gradient
    if (flag == 0 ) {
      for (iVar=0; iVar<prob->nVar; iVar++ ) {
        gradLagrange( iVar ) = gradObj( iVar );
      }
    } else if (flag == 1 ) {
      for (iVar=0; iVar<prob->nVar; iVar++ ) {
        gradLagrange( iVar ) = gradObj( iVar ) - gradLagrange( iVar );
      }
    } else {
      gradLagrange.Initialize( 0.0 );
    }

    // - lambdaT * constrJac
    for (iVar=0; iVar<prob->nVar; iVar++ )
      for (iCon=jacIndCol[iVar]; iCon<jacIndCol[iVar+1]; iCon++ )
        gradLagrange( iVar ) -= lambda( prob->nVar + jacIndRow[iCon] ) * jacNz[iCon];

    // - lambdaT * simpleBounds
    for (iVar=0; iVar<prob->nVar; iVar++ ) gradLagrange( iVar ) -= lambda( iVar );
  }


  /**
   * Compute gradient of Lagrangian or difference of Lagrangian gradients (dense version)
   *
   * flag == 0: output dL(xi, lambda)
   * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
   * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
   */
  void Blocksqp::calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, const Matrix &constrJac,
                                       Matrix &gradLagrange, int flag ) {
    int iVar, iCon;

    // Objective gradient
    if (flag == 0 ) {
      for (iVar=0; iVar<prob->nVar; iVar++ ) {
        gradLagrange( iVar ) = gradObj( iVar );
      }
    } else if (flag == 1 ) {
      for (iVar=0; iVar<prob->nVar; iVar++ ) {
        gradLagrange( iVar ) = gradObj( iVar ) - gradLagrange( iVar );
      }
    } else {
      gradLagrange.Initialize( 0.0 );
    }

    // - lambdaT * constrJac
    for (iVar=0; iVar<prob->nVar; iVar++ )
      for (iCon=0; iCon<prob->nCon; iCon++ )
        gradLagrange( iVar ) -= lambda( prob->nVar + iCon ) * constrJac( iCon, iVar );

    // - lambdaT * simpleBounds
    for (iVar=0; iVar<prob->nVar; iVar++ ) {
      gradLagrange( iVar ) -= lambda( iVar );
    }
  }


  /**
   * Wrapper if called with standard arguments
   */
  void Blocksqp::calcLagrangeGradient( Matrix &gradLagrange, int flag ) {
    if (param->sparseQP ) {
      calcLagrangeGradient( vars->lambda, vars->gradObj, vars->jacNz, vars->jacIndRow, vars->jacIndCol, gradLagrange, flag );
    } else {
      calcLagrangeGradient( vars->lambda, vars->gradObj, vars->constrJac, gradLagrange, flag );
    }
  }


  /**
   * Compute optimality conditions:
   * ||gradLagrange(xi,lambda)||_infty / (1 + ||lambda||_infty) <= TOL
   * and
   * ||constrViolation||_infty / (1 + ||xi||_infty) <= TOL
   */
  bool Blocksqp::calcOptTol() {
    // scaled norm of Lagrangian gradient
    calcLagrangeGradient( vars->gradLagrange, 0 );
    vars->gradNorm = lInfVectorNorm( vars->gradLagrange );
    vars->tol = vars->gradNorm /( 1.0 + lInfVectorNorm( vars->lambda ) );

    // norm of constraint violation
    vars->cNorm  = lInfConstraintNorm( vars->xi, vars->constr, prob->bu, prob->bl );
    vars->cNormS = vars->cNorm /( 1.0 + lInfVectorNorm( vars->xi ) );

    if (vars->tol <= param->opttol && vars->cNormS <= param->nlinfeastol )
      return true;
    else
      return false;
  }

  void Blocksqp::printInfo( int printLevel ) {
    char hessString1[100];
    char hessString2[100];
    char globString[100];
    char qpString[100];

    if (printLevel == 0 )
      return;

    /* QP Solver */
    if (param->sparseQP == 0 )
      strcpy( qpString, "dense, reduced Hessian factorization" );
    else if (param->sparseQP == 1 )
      strcpy( qpString, "sparse, reduced Hessian factorization" );
    else if (param->sparseQP == 2 )
      strcpy( qpString, "sparse, Schur complement approach" );

    /* Globalization */
    if (param->globalization == 0 )
      strcpy( globString, "none (full step)" );
    else if (param->globalization == 1 )
      strcpy( globString, "filter line search" );

    /* Hessian approximation */
    if (param->blockHess && (param->hessUpdate == 1 || param->hessUpdate == 2) )
      strcpy( hessString1, "block " );
    else
      strcpy( hessString1, "" );

    if (param->hessLimMem && (param->hessUpdate == 1 || param->hessUpdate == 2) )
      strcat( hessString1, "L-" );

    /* Fallback Hessian */
    if (param->hessUpdate == 1 || param->hessUpdate == 4 || (param->hessUpdate == 2 && !param->hessDamp) )
      {
        strcpy( hessString2, hessString1 );

        /* Fallback Hessian update type */
        if (param->fallbackUpdate == 0 )
          strcat( hessString2, "Id" );
        else if (param->fallbackUpdate == 1 )
          strcat( hessString2, "SR1" );
        else if (param->fallbackUpdate == 2 )
          strcat( hessString2, "BFGS" );
        else if (param->fallbackUpdate == 4 )
          strcat( hessString2, "Finite differences" );

        /* Fallback Hessian scaling */
        if (param->fallbackScaling == 1 )
          strcat( hessString2, ", SP" );
        else if (param->fallbackScaling == 2 )
          strcat( hessString2, ", OL" );
        else if (param->fallbackScaling == 3 )
          strcat( hessString2, ", mean" );
        else if (param->fallbackScaling == 4 )
          strcat( hessString2, ", selective sizing" );
      }
    else
      strcpy( hessString2, "-" );

    /* First Hessian update type */
    if (param->hessUpdate == 0 )
      strcat( hessString1, "Id" );
    else if (param->hessUpdate == 1 )
      strcat( hessString1, "SR1" );
    else if (param->hessUpdate == 2 )
      strcat( hessString1, "BFGS" );
    else if (param->hessUpdate == 4 )
      strcat( hessString1, "Finite differences" );

    /* First Hessian scaling */
    if (param->hessScaling == 1 )
      strcat( hessString1, ", SP" );
    else if (param->hessScaling == 2 )
      strcat( hessString1, ", OL" );
    else if (param->hessScaling == 3 )
      strcat( hessString1, ", mean" );
    else if (param->hessScaling == 4 )
      strcat( hessString1, ", selective sizing" );

    printf( "\n+---------------------------------------------------------------+\n");
    printf( "| Starting blockSQP with the following algorithmic settings:    |\n");
    printf( "+---------------------------------------------------------------+\n");
    printf( "| qpOASES flavor            | %-34s|\n", qpString );
    printf( "| Globalization             | %-34s|\n", globString );
    printf( "| 1st Hessian approximation | %-34s|\n", hessString1 );
    printf( "| 2nd Hessian approximation | %-34s|\n", hessString2 );
    printf( "+---------------------------------------------------------------+\n\n");
  }

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
    dspev_( "N", "L", &n, temp.ARRAY(), ev.ARRAY(), dummy, &iDummy,
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
    for (i=0; i<nVar; i++ )
      {
        if (xi( i ) - bu( i ) > norm )
          norm = xi( i ) - bu( i );
        else if (bl( i ) - xi( i ) > norm )
          norm = bl( i ) - xi( i );
      }

    // Find out the largest constraint violation
    for (i=0; i<nCon; i++ )
      {
        if (constr( i ) - bu( nVar+i ) > norm )
          norm = constr( i ) - bu( nVar+i );
        if (bl( nVar+i ) - constr( i ) > norm )
          norm = bl( nVar+i ) - constr( i );
      }

    return norm;
  }
  void Blocksqp::acceptStep( double alpha ) {
    acceptStep( vars->deltaXi, vars->lambdaQP, alpha, 0 );
  }

  void Blocksqp::acceptStep( const Matrix &deltaXi, const Matrix &lambdaQP, double alpha, int nSOCS ) {
    int k;
    double lStpNorm;

    // Current alpha
    vars->alpha = alpha;
    vars->nSOCS = nSOCS;

    // Set new xi by accepting the current trial step
    for (k=0; k<vars->xi.M(); k++ )
      {
        vars->xi( k ) = vars->trialXi( k );
        vars->deltaXi( k ) = alpha * deltaXi( k );
      }

    // Store the infinity norm of the multiplier step
    vars->lambdaStepNorm = 0.0;
    for (k=0; k<vars->lambda.M(); k++ )
      if ((lStpNorm = fabs( alpha*lambdaQP( k ) - alpha*vars->lambda( k ) )) > vars->lambdaStepNorm )
        vars->lambdaStepNorm = lStpNorm;

    // Set new multipliers
    for (k=0; k<vars->lambda.M(); k++ )
      vars->lambda( k ) = (1.0 - alpha)*vars->lambda( k ) + alpha*lambdaQP( k );

    // Count consecutive reduced steps
    if (vars->alpha < 1.0 )
      vars->reducedStepCount++;
    else
      vars->reducedStepCount = 0;
  }

  void Blocksqp::reduceStepsize( double *alpha ) {
    *alpha = (*alpha) * 0.5;
  }

  void Blocksqp::reduceSOCStepsize( double *alphaSOC ) {
    int i;
    int nVar = prob->nVar;

    // Update bounds on linearized constraints for the next SOC QP:
    // That is different from the update for the first SOC QP!
    for (i=0; i<prob->nCon; i++ )
      {
        if (prob->bl( nVar+i ) != param->inf )
          vars->deltaBl( nVar+i ) = (*alphaSOC)*vars->deltaBl( nVar+i ) - vars->constr( i );
        else
          vars->deltaBl( nVar+i ) = param->inf;

        if (prob->bu( nVar+i ) != param->inf )
          vars->deltaBu( nVar+i ) = (*alphaSOC)*vars->deltaBu( nVar+i ) - vars->constr( i );
        else
          vars->deltaBu( nVar+i ) = param->inf;
      }

    *alphaSOC = (*alphaSOC) * 0.5;
  }


  /**
   * Take a full Quasi-Newton step, except when integrator fails:
   * xi = xi + deltaXi
   * lambda = lambdaQP
   */
  int Blocksqp::fullstep() {
    double alpha;
    double objTrial, cNormTrial;
    int i, k, info;
    int nVar = prob->nVar;

    // Backtracking line search
    alpha = 1.0;
    for (k=0; k<10; k++ ) {
      // Compute new trial point
      for (i=0; i<nVar; i++ )
        vars->trialXi( i ) = vars->xi( i ) + alpha * vars->deltaXi( i );

      // Compute problem functions at trial point
      prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
      stats->nFunCalls++;
      cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
      // Reduce step if evaluation fails, if lower bound is violated or if objective or a constraint is NaN
      if (info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) ) {
        printf("info=%i, objTrial=%g\n", info, objTrial );
        // evaluation error, reduce stepsize
        reduceStepsize( &alpha );
        continue;
      } else {
        acceptStep( alpha );
        return 0;
      }
    }
    return 1;
  }


  /**
   *
   * Backtracking line search based on a filter
   * as described in Ipopt paper (Waechter 2006)
   *
   */
  int Blocksqp::filterLineSearch() {
    double alpha = 1.0;
    double cNorm, cNormTrial, objTrial, dfTdeltaXi;

    int i, k, info;
    int nVar = prob->nVar;

    // Compute ||constr(xi)|| at old point
    cNorm = lInfConstraintNorm( vars->xi, vars->constr, prob->bu, prob->bl );

    // Backtracking line search
    for (k=0; k<param->maxLineSearch; k++ ) {
      // Compute new trial point
      for (i=0; i<nVar; i++ )
        vars->trialXi( i ) = vars->xi( i ) + alpha * vars->deltaXi( i );

      // Compute grad(f)^T * deltaXi
      dfTdeltaXi = 0.0;
      for (i=0; i<nVar; i++ )
        dfTdeltaXi += vars->gradObj( i ) * vars->deltaXi( i );

      // Compute objective and at ||constr(trialXi)||_1 at trial point
      prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
      stats->nFunCalls++;
      cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
      // Reduce step if evaluation fails, if lower bound is violated or if objective is NaN
      if (info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
        {
          // evaluation error, reduce stepsize
          reduceStepsize( &alpha );
          continue;
        }

      // Check acceptability to the filter
      if (pairInFilter( cNormTrial, objTrial )) {
        // Trial point is in the prohibited region defined by the filter, try second order correction
        if (secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 0, k )) {
          break; // SOC yielded suitable alpha, stop
        } else {
          reduceStepsize( &alpha );
          continue;
        }
      }

      // Check sufficient decrease, case I:
      // If we are (almost) feasible and a "switching condition" is satisfied
      // require sufficient progress in the objective instead of bi-objective condition
      if (cNorm <= param->thetaMin ) {
        // Switching condition, part 1: grad(f)^T * deltaXi < 0 ?
        if (dfTdeltaXi < 0 )
          // Switching condition, part 2: alpha * ( - grad(f)^T * deltaXi )**sF > delta * cNorm**sTheta ?
          if (alpha * pow( (-dfTdeltaXi), param->sF ) > param->delta * pow( cNorm, param->sTheta ) ) {
            // Switching conditions hold: Require satisfaction of Armijo condition for objective
            if (objTrial > vars->obj + param->eta*alpha*dfTdeltaXi ) {
              // Armijo condition violated, try second order correction
              if (secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 1, k ) ) {
                break; // SOC yielded suitable alpha, stop
              } else {
                reduceStepsize( &alpha );
                continue;
              }
            } else {
              // found suitable alpha, stop
              acceptStep( alpha );
              break;
            }
          }
      }

      // Check sufficient decrease, case II:
      // Bi-objective (filter) condition
      if (cNormTrial < (1.0 - param->gammaTheta) * cNorm || objTrial < vars->obj - param->gammaF * cNorm ) {
        // found suitable alpha, stop
        acceptStep( alpha );
        break;
      } else {
        // Trial point is dominated by current point, try second order correction
        if (secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 0, k ) ) {
          break; // SOC yielded suitable alpha, stop
        } else {
          reduceStepsize( &alpha );
          continue;
        }
      }
    }

    // No step could be found by the line search
    if (k == param->maxLineSearch ) return 1;

    // Augment the filter if switching condition or Armijo condition does not hold
    if (dfTdeltaXi >= 0 ) {
      augmentFilter( cNormTrial, objTrial );
    } else if (alpha * pow( (-dfTdeltaXi), param->sF ) > param->delta * pow( cNorm, param->sTheta ) ) { // careful with neg. exponents!
      augmentFilter( cNormTrial, objTrial );
    } else if (objTrial <= vars->obj + param->eta*alpha*dfTdeltaXi ) {
      augmentFilter( cNormTrial, objTrial );
    }

    return 0;
  }


  /**
   *
   * Perform a second order correction step, i.e. solve the QP:
   *
   * min_d d^TBd + d^TgradObj
   * s.t.  bl <= A^Td + constr(xi+alpha*deltaXi) - A^TdeltaXi <= bu
   *
   */
  bool Blocksqp::secondOrderCorrection( double cNorm, double cNormTrial, double dfTdeltaXi, bool swCond, int it ) {
    // Perform SOC only on the first iteration of backtracking line search
    if (it > 0 ) return false;
    // If constraint violation of the trialstep is lower than the current one skip SOC
    if (cNormTrial < cNorm ) return false;

    int nSOCS = 0;
    double cNormTrialSOC, cNormOld, objTrialSOC;
    int i, k, info;
    int nVar = prob->nVar;
    Matrix deltaXiSOC, lambdaQPSOC;

    // vars->constr contains result at first trial point: c(xi+deltaXi)
    // vars->constrJac, vars->AdeltaXi and vars->gradObj are unchanged so far.

    // First SOC step
    deltaXiSOC.Dimension( vars->deltaXi.M() ).Initialize( 0.0 );
    lambdaQPSOC.Dimension( vars->lambdaQP.M() ).Initialize( 0.0 );

    // Second order correction loop
    cNormOld = cNorm;
    for (k=0; k<param->maxSOCiter; k++ ) {
      nSOCS++;

      // Update bounds for SOC QP
      updateStepBounds( 1 );

      // Solve SOC QP to obtain new, corrected deltaXi
      // (store in separate vector to avoid conflict with original deltaXi -> need it in linesearch!)
      info = solveQP( deltaXiSOC, lambdaQPSOC, false );
      if (info != 0 ) return false; // Could not solve QP, abort SOC

      // Set new SOC trial point
      for (i=0; i<nVar; i++ ) {
        vars->trialXi( i ) = vars->xi( i ) + deltaXiSOC( i );
      }

      // Compute objective and ||constr(trialXiSOC)||_1 at SOC trial point
      prob->evaluate( vars->trialXi, &objTrialSOC, vars->constr, &info );
      stats->nFunCalls++;
      cNormTrialSOC = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
      if (info != 0 || objTrialSOC < prob->objLo || objTrialSOC > prob->objUp || !(objTrialSOC == objTrialSOC) || !(cNormTrialSOC == cNormTrialSOC) )
        return false; // evaluation error, abort SOC

      // Check acceptability to the filter (in SOC)
      if (pairInFilter( cNormTrialSOC, objTrialSOC ) ) return false; // Trial point is in the prohibited region defined by the filter, abort SOC

      // Check sufficient decrease, case I (in SOC)
      // (Almost feasible and switching condition holds for line search alpha)
      if (cNorm <= param->thetaMin && swCond ) {
        if (objTrialSOC > vars->obj + param->eta*dfTdeltaXi ) {
          // Armijo condition does not hold for SOC step, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > param->kappaSOC * cNormOld ) {
            return false;
          } else {
            cNormOld = cNormTrialSOC;
          }
          continue;
        } else {
          // found suitable alpha during SOC, stop
          acceptStep( deltaXiSOC, lambdaQPSOC, 1.0, nSOCS );
          return true;
        }
      }

      // Check sufficient decrease, case II (in SOC)
      if (cNorm > param->thetaMin || !swCond ) {
        if (cNormTrialSOC < (1.0 - param->gammaTheta) * cNorm || objTrialSOC < vars->obj - param->gammaF * cNorm ) {
          // found suitable alpha during SOC, stop
          acceptStep( deltaXiSOC, lambdaQPSOC, 1.0, nSOCS );
          return true;
        } else {
          // Trial point is dominated by current point, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > param->kappaSOC * cNormOld ) {
            return false;
          } else {
            cNormOld = cNormTrialSOC;
          }
          continue;
        }
      }
    }

    return false;
  }

  /**
   * Minimize constraint violation by solving an NLP with minimum norm objective
   *
   * "The dreaded restoration phase" -- Nick Gould
   */
  int Blocksqp::feasibilityRestorationPhase() {
    // No Feasibility restoration phase
    if (param->restoreFeas == 0 ) return -1;

    stats->nRestPhaseCalls++;

    int ret, it, i, k, info;
    int maxRestIt = 100;
    int warmStart;
    double cNormTrial, objTrial, lStpNorm;
    RestorationProblem *restProb;
    Blocksqp *restMethod;
    SQPoptions *restOpts;
    SQPstats *restStats;

    // Create a min(constrVio) NLP, an options and a stats object
    restProb = new RestorationProblem( prob, vars->xi );
    restOpts = new SQPoptions();
    restStats = new SQPstats( stats->outpath );

    // Set options for the SQP method for this problem
    restOpts->globalization = 1;
    restOpts->whichSecondDerv = 0;
    restOpts->restoreFeas = 0;
    restOpts->hessUpdate = 2;
    restOpts->hessLimMem = 1;
    restOpts->hessScaling = 2;
    restOpts->opttol = param->opttol;
    restOpts->nlinfeastol = param->nlinfeastol;

    // Create and initialize the SQP method for the minimum norm NLP
    restMethod = new Blocksqp( restProb, restOpts, restStats );
    restMethod->init();

    // Iterate until a point acceptable to the filter is found
    warmStart = 0;
    for (it=0; it<maxRestIt; it++) {
      // One iteration for minimum norm NLP
      ret = restMethod->run( 1, warmStart );
      warmStart = 1;

      // If restMethod yields error, stop restoration phase
      if (ret == -1 ) break;

      // Get new xi from the restoration phase
      for (i=0; i<prob->nVar; i++ ) {
        vars->trialXi( i ) = restMethod->vars->xi( i );
      }

      // Compute objective at trial point
      prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
      stats->nFunCalls++;
      cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
      if (info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) ) continue;

      // Is this iterate acceptable for the filter?
      if (!pairInFilter( cNormTrial, objTrial )) {
        // success
        printf("Found a point acceptable for the filter.\n");
        ret = 0;
        break;
      }

      // If minimum norm NLP has converged, declare local infeasibility
      if (restMethod->vars->tol < param->opttol && restMethod->vars->cNormS < param->nlinfeastol ) {
        ret = 1;
        break;
      }
    }

    // Success or locally infeasible
    if (ret == 0 || ret == 1 ) {
      // Store the infinity norm of the multiplier step
      vars->lambdaStepNorm = 0.0;
      // Compute restoration step
      for (k=0; k<prob->nVar; k++ ) {
        vars->deltaXi( k ) = vars->xi( k );

        vars->xi( k ) = vars->trialXi( k );

        // Store lInf norm of dual step
        if ((lStpNorm = fabs( restMethod->vars->lambda( k ) - vars->lambda( k ) )) > vars->lambdaStepNorm )
          vars->lambdaStepNorm = lStpNorm;
        vars->lambda( k ) = restMethod->vars->lambda( k );
        vars->lambdaQP( k ) = restMethod->vars->lambdaQP( k );

        vars->deltaXi( k ) -= vars->xi( k );
      }
      for (k=prob->nVar; k<prob->nVar+prob->nCon; k++ ) {
        // skip the dual variables for the slack variables in the restoration problem
        if ((lStpNorm = fabs( restMethod->vars->lambda( 2*prob->nCon + k ) - vars->lambda( k ) )) > vars->lambdaStepNorm )
          vars->lambdaStepNorm = lStpNorm;
        vars->lambda( k ) = restMethod->vars->lambda( 2*prob->nCon + k );
        vars->lambdaQP( k ) = restMethod->vars->lambdaQP( 2*prob->nCon + k );
      }
      vars->alpha = 1.0;
      vars->nSOCS = 0;

      // reset reduced step counter
      vars->reducedStepCount = 0;

      // reset Hessian and limited memory information
      resetHessian();
    }

    if (ret == 1 ) {
      stats->printProgress( prob, vars, param, 0 );
      printf("The problem seems to be locally infeasible. Infeasibilities minimized.\n");
    }

    // Clean up
    delete restMethod;
    delete restOpts;
    delete restStats;
    delete restProb;

    return ret;
  }


  /**
   * Try to (partly) improve constraint violation by satisfying
   * the (pseudo) continuity constraints, i.e. do a single shooting
   * iteration with the current controls and measurement weights q and w
   */
  int Blocksqp::feasibilityRestorationHeuristic() {
    stats->nRestHeurCalls++;

    int info, k;
    double cNormTrial;

    info = 0;

    // Call problem specific heuristic to reduce constraint violation.
    // For shooting methods that means setting consistent values for shooting nodes by one forward integration.
    for (k=0; k<prob->nVar; k++ ) // input: last successful step
      vars->trialXi( k ) = vars->xi( k );
    prob->reduceConstrVio( vars->trialXi, &info );
    if (info ) {
      // If an error occured in restoration heuristics, abort
      return -1;
    }

    // Compute objective and constraints at the new (hopefully feasible) point
    prob->evaluate( vars->trialXi, &vars->obj, vars->constr, &info );
    stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
    if (info != 0 || vars->obj < prob->objLo || vars->obj > prob->objUp || !(vars->obj == vars->obj) || !(cNormTrial == cNormTrial) )
      return -1;

    // Is the new point acceptable for the filter?
    if (pairInFilter( cNormTrial, vars->obj)) {
      // point is in the taboo region, restoration heuristic not successful!
      return -1;
    }

    // If no error occured in the integration all shooting variables now
    // have the values obtained by a single shooting integration.
    // This is done instead of a Newton-like step in the current SQP iteration

    vars->alpha = 1.0;
    vars->nSOCS = 0;

    // reset reduced step counter
    vars->reducedStepCount = 0;

    // Reset lambda
    vars->lambda.Initialize( 0.0 );
    vars->lambdaQP.Initialize( 0.0 );

    // Compute the "step" taken by closing the continuity conditions
    /// \note deltaXi is reset by resetHessian(), so this doesn't matter
    for (k=0; k<prob->nVar; k++ ) {
      //vars->deltaXi( k ) = vars->trialXi( k ) - vars->xi( k );
      vars->xi( k ) = vars->trialXi( k );
    }

    // reduce Hessian and limited memory information
    resetHessian();

    return 0;
  }


  /**
   * If the line search fails, check if the full step reduces the KKT error by a factor kappaF.
   */
  int Blocksqp::kktErrorReduction( ) {
    int i, info = 0;
    double objTrial, cNormTrial, trialGradNorm, trialTol;
    Matrix trialConstr, trialGradLagrange;

    // Compute new trial point
    for (i=0; i<prob->nVar; i++ )
      vars->trialXi( i ) = vars->xi( i ) + vars->deltaXi( i );

    // Compute objective and ||constr(trialXi)|| at trial point
    trialConstr.Dimension( prob->nCon ).Initialize( 0.0 );
    prob->evaluate( vars->trialXi, &objTrial, trialConstr, &info );
    stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm( vars->trialXi, trialConstr, prob->bu, prob->bl );
    if (info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) ) {
      // evaluation error
      return 1;
    }

    // Compute KKT error of the new point

    // scaled norm of Lagrangian gradient
    trialGradLagrange.Dimension( prob->nVar ).Initialize( 0.0 );
    if (param->sparseQP ) {
      calcLagrangeGradient( vars->lambdaQP, vars->gradObj, vars->jacNz,
                            vars->jacIndRow, vars->jacIndCol, trialGradLagrange, 0 );
    } else {
      calcLagrangeGradient( vars->lambdaQP, vars->gradObj, vars->constrJac,
                            trialGradLagrange, 0 );
    }

    trialGradNorm = lInfVectorNorm( trialGradLagrange );
    trialTol = trialGradNorm /( 1.0 + lInfVectorNorm( vars->lambdaQP ) );

    if (fmax( cNormTrial, trialTol ) < param->kappaF * fmax( vars->cNorm, vars->tol ) ) {
      acceptStep( 1.0 );
      return 0;
    } else {
      return 1;
    }
  }

  /**
   * Check if current entry is accepted to the filter:
   * (cNorm, obj) in F_k
   */
  bool Blocksqp::pairInFilter( double cNorm, double obj ) {
    std::set< std::pair<double,double> >::iterator iter;
    std::set< std::pair<double,double> > *filter;
    filter = vars->filter;

    /*
     * A pair is in the filter if:
     * - it increases the objective and
     * - it also increases the constraint violation
     * The second expression in the if-clause states that we exclude
     * entries that are within the feasibility tolerance, e.g.
     * if an entry improves the constraint violation from 1e-16 to 1e-17,
     * but increases the objective considerably we also think of this entry
     * as dominated
     */

    for (iter=filter->begin(); iter!=filter->end(); iter++ )
      if ((cNorm >= (1.0 - param->gammaTheta) * iter->first ||
           (cNorm < 0.01 * param->nlinfeastol && iter->first < 0.01 * param->nlinfeastol ) ) &&
          obj >= iter->second - param->gammaF * iter->first ) {
        return 1;
      }

    return 0;
  }


  void Blocksqp::initializeFilter() {
    std::set< std::pair<double,double> >::iterator iter;
    std::pair<double,double> initPair ( param->thetaMax, prob->objLo );

    // Remove all elements
    iter=vars->filter->begin();
    while (iter != vars->filter->end()) {
      std::set< std::pair<double,double> >::iterator iterToRemove = iter;
      iter++;
      vars->filter->erase( iterToRemove );
    }

    // Initialize with pair ( maxConstrViolation, objLowerBound );
    vars->filter->insert( initPair );
  }


  /**
   * Augment the filter:
   * F_k+1 = F_k U { (c,f) | c > (1-gammaTheta)cNorm and f > obj-gammaF*c
   */
  void Blocksqp::augmentFilter( double cNorm, double obj ) {
    std::set< std::pair<double,double> >::iterator iter;
    std::pair<double,double> entry ( (1.0 - param->gammaTheta)*cNorm, obj - param->gammaF*cNorm );

    // Augment filter by current element
    vars->filter->insert( entry );

    // Remove dominated elements
    iter=vars->filter->begin();
    while (iter != vars->filter->end()) {
      //printf(" iter->first=%g, entry.first=%g, iter->second=%g, entry.second=%g\n",iter->first, entry.first, iter->second, entry.second);
      if (iter->first > entry.first && iter->second > entry.second ) {
        std::set< std::pair<double,double> >::iterator iterToRemove = iter;
        iter++;
        vars->filter->erase( iterToRemove );
      } else {
        iter++;
      }
    }
    // Print current filter
    //for (iter=vars->filter->begin(); iter!=vars->filter->end(); iter++ )
    //printf("(%g,%g), ", iter->first, iter->second);
    //printf("\n");
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
      hessNz = NULL;
    }

    jacNz = NULL;
    jacIndCol = NULL;
    jacIndRow = NULL;

    hessIndCol = NULL;
    hessIndRow = NULL;
    hessIndLo = NULL;
    hess = NULL;
    hess1 = NULL;
    hess2 = NULL;

    noUpdateCounter = NULL;

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
    if (iter.jacNz != NULL ) {
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
      jacNz = NULL;
      jacIndRow = NULL;
      jacIndCol = NULL;
    }

    noUpdateCounter = NULL;
    hessNz = NULL;
    hessIndCol = NULL;
    hessIndRow = NULL;
    hessIndLo = NULL;
    hess = NULL;
    hess1 = NULL;
    hess2 = NULL;
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
    if (hessNz == NULL ) return;
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

    if (hessNz_ != NULL ) delete[] hessNz_;
    if (hessIndRow_ != NULL ) delete[] hessIndRow_;

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
    if (blockIdx != NULL )
      delete[] blockIdx;
    if (noUpdateCounter != NULL )
      delete[] noUpdateCounter;
    if (jacNz != NULL )
      delete[] jacNz;
    if (jacIndRow != NULL )
      delete[] jacIndRow;
    if (hessNz != NULL )
      delete[] hessNz;
    if (hessIndRow != NULL )
      delete[] hessIndRow;
  }

  /**
   * Initial Hessian: Identity matrix
   */
  void Blocksqp::calcInitialHessian()
  {
    int iBlock;

    for (iBlock=0; iBlock<vars->nBlocks; iBlock++ )
      //if objective derv is computed exactly, don't set the last block!
      if (!(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks-1) )
        calcInitialHessian( iBlock );
  }


  /**
   * Initial Hessian for one block: Identity matrix
   */
  void Blocksqp::calcInitialHessian( int iBlock ) {
    vars->hess[iBlock].Initialize( 0.0 );

    // Each block is a diagonal matrix
    for (int i=0; i<vars->hess[iBlock].M(); i++ )
      vars->hess[iBlock]( i, i ) = param->iniHessDiag;

    // If we maintain 2 Hessians, also reset the second one
    if (vars->hess2 != NULL ) {
      vars->hess2[iBlock].Initialize( 0.0 );
      for (int i=0; i<vars->hess2[iBlock].M(); i++ )
        vars->hess2[iBlock]( i, i ) = param->iniHessDiag;
    }
  }


  void Blocksqp::resetHessian() {
    for (int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
      //if objective derv is computed exactly, don't set the last block!
      if (!(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks - 1) )
        resetHessian( iBlock );
  }


  void Blocksqp::resetHessian( int iBlock ) {
    Matrix smallDelta, smallGamma;
    int nVarLocal = vars->hess[iBlock].M();

    // smallGamma and smallDelta are either subvectors of gamma and delta
    // or submatrices of gammaMat, deltaMat, i.e. subvectors of gamma and delta from m prev. iterations (for L-BFGS)
    smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
    smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

    // Remove past information on Lagrangian gradient difference
    smallGamma.Initialize( 0.0 );

    // Remove past information on steps
    smallDelta.Initialize( 0.0 );

    // Remove information on old scalars (used for COL sizing)
    vars->deltaNorm( iBlock ) = 1.0;
    vars->deltaGamma( iBlock ) = 0.0;
    vars->deltaNormOld( iBlock ) = 1.0;
    vars->deltaGammaOld( iBlock ) = 0.0;

    vars->noUpdateCounter[iBlock] = -1;

    calcInitialHessian( iBlock );
  }

  /**
   * Approximate Hessian by finite differences
   */
  int Blocksqp::calcFiniteDiffHessian() {
    int iVar, jVar, k, iBlock, maxBlock, info, idx, idx1, idx2;
    double dummy, lowerVio, upperVio;
    Matrix pert;
    SQPiterate varsP = SQPiterate( *vars );

    const double myDelta = 1.0e-4;
    const double minDelta = 1.0e-6;

    pert.Dimension( prob->nVar );

    info = 0;

    // Find out the largest block
    maxBlock = 0;
    for (iBlock=0; iBlock<vars->nBlocks; iBlock++ )
      if (vars->blockIdx[iBlock+1] - vars->blockIdx[iBlock] > maxBlock )
        maxBlock = vars->blockIdx[iBlock+1] - vars->blockIdx[iBlock];

    // Compute original Lagrange gradient
    calcLagrangeGradient( vars->lambda, vars->gradObj, vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->gradLagrange, 0 );

    for (iVar = 0; iVar<maxBlock; iVar++ ) {
      pert.Initialize( 0.0 );

      // Perturb all blocks simultaneously
      for (iBlock=0; iBlock<vars->nBlocks; iBlock++ ) {
        idx = vars->blockIdx[iBlock] + iVar;
        // Skip blocks that have less than iVar variables
        if (idx < vars->blockIdx[iBlock+1] ) {
          pert( idx ) = myDelta * fabs( vars->xi( idx ) );
          pert( idx ) = fmax( pert( idx ), minDelta );

          // If perturbation violates upper bound, try to perturb with negative
          upperVio = vars->xi( idx ) + pert( idx ) - prob->bu( idx );
          if (upperVio > 0 ) {
            lowerVio = prob->bl( idx ) -  ( vars->xi( idx ) - pert( idx ) );
            // If perturbation violates also lower bound, take the largest perturbation possible
            if (lowerVio > 0 ) {
              if (lowerVio > upperVio )
                pert( idx ) = -lowerVio;
              else
                pert( idx ) = upperVio;
            } else {
              // If perturbation does not violate lower bound, take -computed perturbation
              pert( idx ) = -pert( idx );
            }
          }
        }
      }

      // Add perturbation
      for (k=0; k<prob->nVar; k++ )
        vars->xi( k ) += pert( k );

      // Compute perturbed Lagrange gradient
      if (param->sparseQP ) {
        prob->evaluate( vars->xi, vars->lambda, &dummy, varsP.constr, varsP.gradObj,
                        varsP.jacNz, varsP.jacIndRow, varsP.jacIndCol, vars->hess, 1, &info );
        calcLagrangeGradient( vars->lambda, varsP.gradObj, varsP.jacNz, varsP.jacIndRow,
                              varsP.jacIndCol, varsP.gradLagrange, 0 );
      } else {
        prob->evaluate( vars->xi, vars->lambda, &dummy, varsP.constr, varsP.gradObj, varsP.constrJac, vars->hess, 1, &info );
        calcLagrangeGradient( vars->lambda, varsP.gradObj, varsP.constrJac, varsP.gradLagrange, 0 );
      }

      // Compute finite difference approximations: one column in every block
      for (iBlock=0; iBlock<vars->nBlocks; iBlock++ ) {
        idx1 = vars->blockIdx[iBlock] + iVar;
        // Skip blocks that have less than iVar variables
        if (idx1 < vars->blockIdx[iBlock+1] ) {
          for (jVar=iVar; jVar<vars->blockIdx[iBlock+1]-vars->blockIdx[iBlock]; jVar++ )
            {// Take symmetrized matrices
              idx2 = vars->blockIdx[iBlock] + jVar;
              vars->hess[iBlock]( iVar, jVar ) =  ( varsP.gradLagrange( idx1 ) - vars->gradLagrange( idx2 ) );
              vars->hess[iBlock]( iVar, jVar ) += ( varsP.gradLagrange( idx2 ) - vars->gradLagrange( idx1 ) );
              vars->hess[iBlock]( iVar, jVar ) *= 0.5 / pert( idx1 );
            }
        }
      }

      // Subtract perturbation
      for (k=0; k<prob->nVar; k++ ) vars->xi( k ) -= pert( k );
    }

    return info;
  }


  void Blocksqp::sizeInitialHessian( const Matrix &gamma, const Matrix &delta, int iBlock, int option ) {
    int i, j;
    double scale;
    double myEps = 1.0e3 * param->eps;

    if (option == 1 ) {
      // Shanno-Phua
      scale = adotb( gamma, gamma ) / fmax( adotb( delta, gamma ), myEps );
    }
    else if (option == 2 ) {
      // Oren-Luenberger
      scale = adotb( delta, gamma ) / fmax( adotb( delta, delta ), myEps );
      scale = fmin( scale, 1.0 );
    } else if (option == 3 ) {
      // Geometric mean of 1 and 2
      scale = sqrt( adotb( gamma, gamma ) / fmax( adotb( delta, delta ), myEps ) );
    } else {
      // Invalid option, ignore
      return;
    }

    if (scale > 0.0 ) {
      scale = fmax( scale, myEps );
      for (i=0; i<vars->hess[iBlock].M(); i++ )
        for (j=i; j<vars->hess[iBlock].M(); j++ )
          vars->hess[iBlock]( i,j ) *= scale;
    } else {
      scale = 1.0;
    }

    // statistics: average sizing factor
    stats->averageSizingFactor += scale;
  }


  void Blocksqp::sizeHessianCOL( const Matrix &gamma, const Matrix &delta, int iBlock ) {
    int i, j;
    double theta, scale, myEps = 1.0e3 * param->eps;
    double deltaNorm, deltaNormOld, deltaGamma, deltaGammaOld, deltaBdelta;

    // Get sTs, sTs_, sTy, sTy_, sTBs
    deltaNorm = vars->deltaNorm(iBlock);
    deltaGamma = vars->deltaGamma(iBlock);
    deltaNormOld = vars->deltaNormOld(iBlock);
    deltaGammaOld = vars->deltaGammaOld(iBlock);
    deltaBdelta = 0.0;
    for (i=0; i<delta.M(); i++ )
      for (j=0; j<delta.M(); j++ )
        deltaBdelta += delta( i ) * vars->hess[iBlock]( i, j ) * delta( j );

    // Centered Oren-Luenberger factor
    if (vars->noUpdateCounter[iBlock] == -1 ) {
      // in the first iteration, this should equal the OL factor
      theta = 1.0;
    } else {
      theta = fmin( param->colTau1, param->colTau2 * deltaNorm );
    }
    if (deltaNorm > myEps && deltaNormOld > myEps ) {
      scale = (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaBdelta / deltaNorm;
      if (scale > param->eps )
        scale = ( (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaGamma / deltaNorm ) / scale;
    } else {
      scale = 1.0;
    }

    // Size only if factor is between zero and one
    if (scale < 1.0 && scale > 0.0 ) {
      scale = fmax( param->colEps, scale );
      //printf("Sizing value (COL) block %i = %g\n", iBlock, scale );
      for (i=0; i<vars->hess[iBlock].M(); i++ )
        for (j=i; j<vars->hess[iBlock].M(); j++ )
          vars->hess[iBlock]( i,j ) *= scale;

      // statistics: average sizing factor
      stats->averageSizingFactor += scale;
    } else {
      stats->averageSizingFactor += 1.0;
    }
  }

  /**
   * Apply BFGS or SR1 update blockwise and size blocks
   */
  void Blocksqp::calcHessianUpdate( int updateType, int hessScaling ) {
    int iBlock, nBlocks;
    int nVarLocal;
    Matrix smallGamma, smallDelta;
    bool firstIter;

    //if objective derv is computed exactly, don't set the last block!
    if (param->whichSecondDerv == 1 && param->blockHess )
      nBlocks = vars->nBlocks - 1;
    else
      nBlocks = vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->averageSizingFactor = 0.0;

    for (iBlock=0; iBlock<nBlocks; iBlock++ ) {
      nVarLocal = vars->hess[iBlock].M();

      // smallGamma and smallDelta are subvectors of gamma and delta, corresponding to partially separability
      smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
      smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

      // Is this the first iteration or the first after a Hessian reset?
      firstIter = ( vars->noUpdateCounter[iBlock] == -1 );

      // Update sTs, sTs_ and sTy, sTy_
      vars->deltaNormOld(iBlock) = vars->deltaNorm(iBlock);
      vars->deltaGammaOld(iBlock) = vars->deltaGamma(iBlock);
      vars->deltaNorm(iBlock) = adotb( smallDelta, smallDelta );
      vars->deltaGamma(iBlock) = adotb( smallDelta, smallGamma );

      // Sizing before the update
      if (hessScaling < 4 && firstIter )
        sizeInitialHessian( smallGamma, smallDelta, iBlock, hessScaling );
      else if (hessScaling == 4 )
        sizeHessianCOL( smallGamma, smallDelta, iBlock );

      // Compute the new update
      if (updateType == 1 ) {
        calcSR1( smallGamma, smallDelta, iBlock );

        // Prepare to compute fallback update as well
        vars->hess = vars->hess2;

        // Sizing the fallback update
        if (param->fallbackScaling < 4 && firstIter )
          sizeInitialHessian( smallGamma, smallDelta, iBlock, param->fallbackScaling );
        else if (param->fallbackScaling == 4 )
          sizeHessianCOL( smallGamma, smallDelta, iBlock );

        // Compute fallback update
        if (param->fallbackUpdate == 2 )
          calcBFGS( smallGamma, smallDelta, iBlock );

        // Reset pointer
        vars->hess = vars->hess1;
      } else if (updateType == 2 ) {
        calcBFGS( smallGamma, smallDelta, iBlock );
      }

      // If an update is skipped to often, reset Hessian block
      if (vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates ) {
        resetHessian( iBlock );
      }
    }

    // statistics: average sizing factor
    stats->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::calcHessianUpdateLimitedMemory( int updateType, int hessScaling ) {
    int iBlock, nBlocks, nVarLocal;
    Matrix smallGamma, smallDelta;
    Matrix gammai, deltai;
    int i, m, pos, posOldest, posNewest;
    int hessDamped, hessSkipped;
    double averageSizingFactor;

    //if objective derv is computed exactly, don't set the last block!
    if (param->whichSecondDerv == 1 && param->blockHess ) {
      nBlocks = vars->nBlocks - 1;
    } else {
      nBlocks = vars->nBlocks;
    }

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->hessSkipped = 0;
    stats->averageSizingFactor = 0.0;

    for (iBlock=0; iBlock<nBlocks; iBlock++ ) {
      nVarLocal = vars->hess[iBlock].M();

      // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
      // i.e. subvectors of gamma and delta from m prev. iterations
      smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
      smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

      // Memory structure
      if (stats->itCount > smallGamma.N() ) {
        m = smallGamma.N();
        posOldest = stats->itCount % m;
        posNewest = (stats->itCount-1) % m;
      } else {
        m = stats->itCount;
        posOldest = 0;
        posNewest = m-1;
      }

      // Set B_0 (pretend it's the first step)
      calcInitialHessian( iBlock );
      vars->deltaNorm( iBlock ) = 1.0;
      vars->deltaNormOld( iBlock ) = 1.0;
      vars->deltaGamma( iBlock ) = 0.0;
      vars->deltaGammaOld( iBlock ) = 0.0;
      vars->noUpdateCounter[iBlock] = -1;

      // Size the initial update, but with the most recent delta/gamma-pair
      gammai.Submatrix( smallGamma, nVarLocal, 1, 0, posNewest );
      deltai.Submatrix( smallDelta, nVarLocal, 1, 0, posNewest );
      sizeInitialHessian( gammai, deltai, iBlock, hessScaling );

      for (i=0; i<m; i++ ) {
        pos = (posOldest+i) % m;

        // Get new vector from list
        gammai.Submatrix( smallGamma, nVarLocal, 1, 0, pos );
        deltai.Submatrix( smallDelta, nVarLocal, 1, 0, pos );

        // Update sTs, sTs_ and sTy, sTy_
        vars->deltaNormOld(iBlock) = vars->deltaNorm(iBlock);
        vars->deltaGammaOld(iBlock) = vars->deltaGamma(iBlock);
        vars->deltaNorm(iBlock) = adotb( deltai, deltai );
        vars->deltaGamma(iBlock) = adotb( gammai, deltai );

        // Save statistics, we want to record them only for the most recent update
        averageSizingFactor = stats->averageSizingFactor;
        hessDamped = stats->hessDamped;
        hessSkipped = stats->hessSkipped;

        // Selective sizing before the update
        if (hessScaling == 4 )
          sizeHessianCOL( gammai, deltai, iBlock );

        // Compute the new update
        if (updateType == 1 )
          calcSR1( gammai, deltai, iBlock );
        else if (updateType == 2 )
          calcBFGS( gammai, deltai, iBlock );

        stats->nTotalUpdates++;

        // Count damping statistics only for the most recent update
        if (pos != posNewest ) {
          stats->hessDamped = hessDamped;
          stats->hessSkipped = hessSkipped;
          if (hessScaling == 4 )
            stats->averageSizingFactor = averageSizingFactor;
        }
      }

      // If an update is skipped to often, reset Hessian block
      if (vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates ) {
        resetHessian( iBlock );
      }
    }//blocks
    stats->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::calcBFGS( const Matrix &gamma, const Matrix &delta, int iBlock ) {
    int i, j, k, dim = gamma.M();
    Matrix Bdelta;
    SymMatrix *B;
    double h1 = 0.0;
    double h2 = 0.0;
    double thetaPowell = 0.0;
    int damped;

    /* Work with a local copy of gamma because damping may need to change gamma.
     * Note that vars->gamma needs to remain unchanged!
     * This may be important in a limited memory context:
     * When information is "forgotten", B_i-1 is different and the
     *  original gamma might lead to an undamped update with the new B_i-1! */
    Matrix gamma2 = gamma;

    B = &vars->hess[iBlock];

    // Bdelta = B*delta (if sizing is enabled, B is the sized B!)
    // h1 = delta^T * B * delta
    // h2 = delta^T * gamma
    Bdelta.Dimension( dim ).Initialize( 0.0 );
    for (i=0; i<dim; i++ )
      {
        for (k=0; k<dim; k++ )
          Bdelta( i ) += (*B)( i,k ) * delta( k );

        h1 += delta( i ) * Bdelta( i );
        //h2 += delta( i ) * gamma( i );
      }
    h2 = vars->deltaGamma( iBlock );

    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    damped = 0;
    if (param->hessDamp )
      if (h2 < param->hessDampFac * h1 / vars->alpha && fabs( h1 - h2 ) > 1.0e-12 ) {
        // At the first iteration h1 and h2 are equal due to COL scaling

        thetaPowell = (1.0 - param->hessDampFac)*h1 / ( h1 - h2 );

        // Redefine gamma and h2 = delta^T * gamma
        h2 = 0.0;
        for (i=0; i<dim; i++ ) {
          gamma2( i ) = thetaPowell*gamma2( i ) + (1.0 - thetaPowell)*Bdelta( i );
          h2 += delta( i ) * gamma2( i );
        }

        // Also redefine deltaGamma for computation of sizing factor in the next iteration
        vars->deltaGamma( iBlock ) = h2;

        damped = 1;
      }

    // For statistics: count number of damped blocks
    stats->hessDamped += damped;

    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e2 * param->eps;
    if (fabs( h1 ) < myEps || fabs( h2 ) < myEps ) {
      // don't perform update because of bad condition, might introduce negative eigenvalues
      vars->noUpdateCounter[iBlock]++;
      stats->hessDamped -= damped;
      stats->hessSkipped++;
      stats->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++ )
        for (j=i; j<dim; j++ )
          (*B)( i,j ) = (*B)( i,j ) - Bdelta( i ) * Bdelta( j ) / h1
            + gamma2( i ) * gamma2( j ) / h2;

      vars->noUpdateCounter[iBlock] = 0;
    }
  }


  void Blocksqp::calcSR1( const Matrix &gamma, const Matrix &delta, int iBlock ) {
    int i, j, k, dim = gamma.M();
    Matrix gmBdelta;
    SymMatrix *B;
    double myEps = 1.0e2 * param->eps;
    double r = 1.0e-8;
    double h = 0.0;

    B = &vars->hess[iBlock];

    // gmBdelta = gamma - B*delta
    // h = (gamma - B*delta)^T * delta
    gmBdelta.Dimension( dim );
    for (i=0; i<dim; i++ ) {
      gmBdelta( i ) = gamma( i );
      for (k=0; k<dim; k++ )
        gmBdelta( i ) -= ( (*B)( i,k ) * delta( k ) );

      h += ( gmBdelta( i ) * delta( i ) );
    }

    // B_k+1 = B_k + gmBdelta * gmBdelta^T / h
    if (fabs( h ) < r * l2VectorNorm( delta ) * l2VectorNorm( gmBdelta ) || fabs( h ) < myEps ) {
      // Skip update if denominator is too small
      vars->noUpdateCounter[iBlock]++;
      stats->hessSkipped++;
      stats->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++ )
        for (j=i; j<dim; j++ )
          (*B)( i,j ) = (*B)( i,j ) + gmBdelta( i ) * gmBdelta( j ) / h;
      vars->noUpdateCounter[iBlock] = 0;
    }
  }


  /**
   * Set deltaXi and gamma as a column in the matrix containing
   * the m most recent delta and gamma
   */
  void Blocksqp::updateDeltaGamma() {
    int nVar = vars->gammaMat.M();
    int m = vars->gammaMat.N();

    if (m == 1 )
      return;

    vars->deltaXi.Submatrix( vars->deltaMat, nVar, 1, 0, stats->itCount % m );
    vars->gamma.Submatrix( vars->gammaMat, nVar, 1, 0, stats->itCount % m );
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
      array = NULL;
    else
      if ( ( array = new double[len] ) == NULL )
        Error("'new' failed");

    return 0;
  }


  int Matrix::free( void ) {
    if ( tflag )
      Error("free cannot be called with Submatrix");

    if ( array != NULL )
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
      array = NULL;
    else
      if ( ( array = new double[len] ) == NULL )
        Error("'new' failed");

    return 0;
  }


  int SymMatrix::free( void ) {
    if (array != NULL )
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

    if ( ( array = new double[A.N()*A.M()] ) == NULL )
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

  void Blocksqp::computeNextHessian( int idx, int maxQP ) {
    // Compute fallback update only once
    if (idx == 1 )
      {
        // Switch storage
        vars->hess = vars->hess2;

        // If last block contains exact Hessian, we need to copy it
        if (param->whichSecondDerv == 1 )
          for (int i=0; i<vars->hess[prob->nBlocks-1].M(); i++ )
            for (int j=i; j<vars->hess[prob->nBlocks-1].N(); j++ )
              vars->hess2[prob->nBlocks-1]( i,j ) = vars->hess1[prob->nBlocks-1]( i,j );

        // Limited memory: compute fallback update only when needed
        if (param->hessLimMem )
          {
            stats->itCount--;
            int hessDampSave = param->hessDamp;
            param->hessDamp = 1;
            calcHessianUpdateLimitedMemory( param->fallbackUpdate, param->fallbackScaling );
            param->hessDamp = hessDampSave;
            stats->itCount++;
          }
        /* Full memory: both updates must be computed in every iteration
         * so switching storage is enough */
      }

    // 'Nontrivial' convex combinations
    if (maxQP > 2 )
      {
        /* Convexification parameter: mu_l = l / (maxQP-1).
         * Compute it only in the first iteration, afterwards update
         * by recursion: mu_l/mu_(l-1) */
        double idxF = (double) idx;
        double mu = (idx==1) ? 1.0 / (maxQP-1) : idxF / (idxF - 1.0);
        double mu1 = 1.0 - mu;
        for (int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
          for (int i=0; i<vars->hess[iBlock].M(); i++ )
            for (int j=i; j<vars->hess[iBlock].N(); j++ )
              {
                vars->hess2[iBlock]( i,j ) *= mu;
                vars->hess2[iBlock]( i,j ) += mu1 * vars->hess1[iBlock]( i,j );
              }
      }
  }


  /**
   * Inner loop of SQP algorithm:
   * Solve a sequence of QPs until pos. def. assumption (G3*) is satisfied.
   */
  int Blocksqp::solveQP( Matrix &deltaXi, Matrix &lambdaQP, bool matricesChanged ) {
    Matrix jacT;
    int maxQP, l;
    if (param->globalization == 1 &&
        param->hessUpdate == 1 &&
        matricesChanged &&
        stats->itCount > 1 )
      {
        maxQP = param->maxConvQP + 1;
      }
    else
      maxQP = 1;

    /*
     * Prepare for qpOASES
     */

    // Setup QProblem data
    qpOASES::Matrix *A;
    qpOASES::SymmetricMatrix *H;
    if (matricesChanged )
      {
        if (param->sparseQP )
          {
            A = new qpOASES::SparseMatrix( prob->nCon, prob->nVar,
                                           vars->jacIndRow, vars->jacIndCol, vars->jacNz );
          }
        else
          {
            // transpose Jacobian (qpOASES needs row major arrays)
            Transpose( vars->constrJac, jacT );
            A = new qpOASES::DenseMatrix( prob->nCon, prob->nVar, prob->nVar, jacT.ARRAY() );
          }
      }
    double *g = vars->gradObj.ARRAY();
    double *lb = vars->deltaBl.ARRAY();
    double *lu = vars->deltaBu.ARRAY();
    double *lbA = vars->deltaBl.ARRAY() + prob->nVar;
    double *luA = vars->deltaBu.ARRAY() + prob->nVar;

    // qpOASES options
    qpOASES::Options opts;
    if (matricesChanged && maxQP > 1 )
      opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    else
      opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    qp->setOptions( opts );

    if (maxQP > 1 )
      {
        // Store last successful QP in temporary storage
        (*qpSave) = *qp;
        /** \todo Storing the active set would be enough but then the QP object
         *        must be properly reset after unsuccessful (SR1-)attempt.
         *        Moreover, passing a guessed active set doesn't yield
         *        exactly the same result as hotstarting a QP. This has
         *        something to do with how qpOASES handles user-given
         *        active sets (->check qpOASES source code). */
      }

    // Other variables for qpOASES
    double cpuTime = matricesChanged ? param->maxTimeQP : 0.1*param->maxTimeQP;
    int maxIt = matricesChanged ? param->maxItQP : 0.1*param->maxItQP;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue ret;

    /*
     * QP solving loop for convex combinations (sequential)
     */
    for (l=0; l<maxQP; l++ )
      {
        /*
         * Compute a new Hessian
         */
        if (l > 0 )
          {// If the solution of the first QP was rejected, consider second Hessian
            stats->qpResolve++;
            *qp = *qpSave;

            computeNextHessian( l, maxQP );
          }

        if (l == maxQP-1 )
          {// Enable inertia correction for supposedly convex QPs, just in case
            opts.enableInertiaCorrection = qpOASES::BT_TRUE;
            qp->setOptions( opts );
          }

        /*
         * Prepare the current Hessian for qpOASES
         */
        if (matricesChanged )
          {
            if (param->sparseQP )
              {
                // Convert block-Hessian to sparse format
                vars->convertHessian( prob, param->eps, vars->hess, vars->hessNz,
                                      vars->hessIndRow, vars->hessIndCol, vars->hessIndLo );
                H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                               vars->hessIndRow, vars->hessIndCol, vars->hessNz );
                dynamic_cast<qpOASES::SymSparseMat*>(H)->createDiagInfo();
              }
            else
              {
                // Convert block-Hessian to double array
                vars->convertHessian( prob, param->eps, vars->hess );
                H = new qpOASES::SymDenseMat( prob->nVar, prob->nVar, prob->nVar, vars->hessNz );
              }
          }

        /*
         * Call qpOASES
         */
        if (param->debugLevel > 2 ) stats->dumpQPCpp( prob, vars, qp, param->sparseQP );
        if (matricesChanged )
          {
            maxIt = param->maxItQP;
            cpuTime = param->maxTimeQP;
            if (qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                qp->getStatus() == qpOASES::QPS_SOLVED )
              {
                ret = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
              }
            else
              {
                ret = qp->init( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
              }
          }
        else if (!matricesChanged ) // Second order correction: H and A do not change
          {
            maxIt = 0.1*param->maxItQP;
            cpuTime = 0.1*param->maxTimeQP;
            ret = qp->hotstart( g, lb, lu, lbA, luA, maxIt, &cpuTime );
          }

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if (l < maxQP-1 && matricesChanged )
          {
            if (ret == qpOASES::SUCCESSFUL_RETURN )
              {
                if (param->sparseQP == 2 )
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints( dynamic_cast<qpOASES::SQProblemSchur*>(qp) );
                else
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints( qp );
              }

            if (ret == qpOASES::SUCCESSFUL_RETURN )
              {// QP was solved successfully and curvature is positive after removing bounds
                stats->qpIterations = maxIt + 1;
                break; // Success!
              }
            else
              {// QP solution is rejected, save statistics
                if (ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
                  stats->qpIterations2++;
                else
                  stats->qpIterations2 += maxIt + 1;
                stats->rejectedSR1++;
              }
          }
        else // Convex QP was solved, no need to check assumption (G3*)
          stats->qpIterations += maxIt + 1;

      } // End of QP solving loop

    /*
     * Post-processing
     */

    // Get solution from qpOASES
    qp->getPrimalSolution( deltaXi.ARRAY() );
    qp->getDualSolution( lambdaQP.ARRAY() );
    vars->qpObj = qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    if (param->sparseQP )
      Atimesb( vars->jacNz, vars->jacIndRow, vars->jacIndCol, deltaXi, vars->AdeltaXi );
    else
      Atimesb( vars->constrJac, deltaXi, vars->AdeltaXi );

    // Print qpOASES error code, if any
    if (ret != qpOASES::SUCCESSFUL_RETURN && matricesChanged )
      printf( "qpOASES error message: \"%s\"\n",
              qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );

    // Point Hessian again to the first Hessian
    vars->hess = vars->hess1;

    /* For full-memory Hessian: Restore fallback Hessian if convex combinations
     * were used during the loop */
    if (!param->hessLimMem && maxQP > 2 && matricesChanged )
      {
        double mu = 1.0 / ((double) l);
        double mu1 = 1.0 - mu;
        int nBlocks = (param->whichSecondDerv == 1) ? vars->nBlocks-1 : vars->nBlocks;
        for (int iBlock=0; iBlock<nBlocks; iBlock++ )
          for (int i=0; i<vars->hess[iBlock].M(); i++ )
            for (int j=i; j<vars->hess[iBlock].N(); j++ )
              {
                vars->hess2[iBlock]( i,j ) *= mu;
                vars->hess2[iBlock]( i,j ) += mu1 * vars->hess1[iBlock]( i,j );
              }
      }

    /* Return code depending on qpOASES returnvalue
     * 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if (ret == qpOASES::SUCCESSFUL_RETURN )
      return 0;
    else if (ret == qpOASES::RET_MAX_NWSR_REACHED )
      return 1;
    else if (ret == qpOASES::RET_HESSIAN_NOT_SPD ||
             ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS ||
             ret == qpOASES::RET_QP_UNBOUNDED ||
             ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS )
      return 2;
    else if (ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY ||
             ret == qpOASES::RET_QP_INFEASIBLE ||
             ret == qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY )
      return 3;
    else
      return 4;
  }


  /**
   * Set bounds on the step (in the QP), either according
   * to variable bounds in the NLP or according to
   * trust region box radius
   */
  void Blocksqp::updateStepBounds( bool soc ) {
    int i;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // Bounds on step
    for (i=0; i<nVar; i++ ) {
      if (prob->bl(i) != param->inf )
        vars->deltaBl( i ) = prob->bl( i ) - vars->xi( i );
      else
        vars->deltaBl( i ) = param->inf;

      if (prob->bu(i) != param->inf )
        vars->deltaBu( i ) = prob->bu( i ) - vars->xi( i );
      else
        vars->deltaBu( i ) = param->inf;
    }

    // Bounds on linearized constraints
    for (i=0; i<nCon; i++ ) {
      if (prob->bl( nVar+i ) != param->inf ) {
        vars->deltaBl( nVar+i ) = prob->bl( nVar+i ) - vars->constr( i );
        if (soc ) vars->deltaBl( nVar+i ) += vars->AdeltaXi( i );
      } else {
        vars->deltaBl( nVar+i ) = param->inf;
      }

      if (prob->bu( nVar+i ) != param->inf ) {
        vars->deltaBu( nVar+i ) = prob->bu( nVar+i ) - vars->constr( i );
        if (soc ) vars->deltaBu( nVar+i ) += vars->AdeltaXi( i );
      } else {
        vars->deltaBu( nVar+i ) = param->inf;
      }
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

  SQPstats::SQPstats( PATHSTR myOutpath ) {
    strcpy( outpath, myOutpath );

    itCount = 0;
    qpItTotal = 0;
    qpIterations = 0;
    qpIterations2 = 0;
    qpResolve = 0;
    rejectedSR1 = 0;
    hessSkipped = 0;
    hessDamped = 0;
    averageSizingFactor = 0.0;
    nFunCalls = 0;
    nDerCalls = 0;
    nRestHeurCalls = 0;
    nRestPhaseCalls = 0;

    nTotalUpdates = 0;
    nTotalSkippedUpdates = 0;
  }


  void SQPstats::printProgress( Problemspec *prob, SQPiterate *vars, SQPoptions *param, bool hasConverged ) {
    /*
     * vars->steptype:
     *-1: full step was accepted because it reduces the KKT error although line search failed
     * 0: standard line search step
     * 1: Hessian has been reset to identity
     * 2: feasibility restoration heuristic has been called
     * 3: feasibility restoration phase has been called
     */

    if (itCount == 0 ) {
      if (param->printLevel > 0 ) {
        prob->printInfo();

        // Headline
        printf("%-8s", "   it" );
        printf("%-21s", " qpIt" );
        printf("%-9s","obj" );
        printf("%-11s","feas" );
        printf("%-7s","opt" );
        if (param->printLevel > 1 ) {
          printf("%-11s","|lgrd|" );
          printf("%-9s","|stp|" );
          printf("%-10s","|lstp|" );
        }
        printf("%-8s","alpha" );
        if (param->printLevel > 1 ) {
          printf("%-6s","nSOCS" );
          printf("%-18s","sk, da, sca" );
          printf("%-6s","QPr,mu" );
        }
        printf("\n");

        // Values for first iteration
        printf("%5i  ", itCount );
        printf("%11i ", 0 );
        printf("% 10e  ", vars->obj );
        printf("%-10.2e", vars->cNormS );
        printf("%-10.2e", vars->tol );
        printf("\n");
      }

      if (param->debugLevel > 0 ) {
        // Print everything in a CSV file as well
        fprintf( progressFile, "%23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %i, %i, %23.16e, %i, %23.16e\n",
                 vars->obj, vars->cNormS, vars->tol, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0, 0.0 );
      }
    } else {
      // Every twenty iterations print headline
      if (itCount % 20 == 0 && param->printLevel > 0 ) {
        printf("%-8s", "   it" );
        printf("%-21s", " qpIt" );
        printf("%-9s","obj" );
        printf("%-11s","feas" );
        printf("%-7s","opt" );
        if (param->printLevel > 1 )
          {
            printf("%-11s","|lgrd|" );
            printf("%-9s","|stp|" );
            printf("%-10s","|lstp|" );
          }
        printf("%-8s","alpha" );
        if (param->printLevel > 1 )
          {
            printf("%-6s","nSOCS" );
            printf("%-18s","sk, da, sca" );
            printf("%-6s","QPr,mu" );
          }
        printf("\n");
      }

      // All values
      if (param->printLevel > 0 ) {
        printf("%5i  ", itCount );
        printf("%5i+%5i ", qpIterations, qpIterations2 );
        printf("% 10e  ", vars->obj );
        printf("%-10.2e", vars->cNormS );
        printf("%-10.2e", vars->tol );
        if (param->printLevel > 1 )
          {
            printf("%-10.2e", vars->gradNorm );
            printf("%-10.2e", lInfVectorNorm( vars->deltaXi ) );
            printf("%-10.2e", vars->lambdaStepNorm );
          }

        if ((vars->alpha == 1.0 && vars->steptype != -1) || !param->printColor ) {
          printf("%-9.1e", vars->alpha );
        } else {
          printf("\033[0;36m%-9.1e\033[0m", vars->alpha );
        }

        if (param->printLevel > 1 ) {
          if (vars->nSOCS == 0 || !param->printColor ) {
            printf("%5i", vars->nSOCS );
          } else {
            printf("\033[0;36m%5i\033[0m", vars->nSOCS );
          }
          printf("%3i, %3i, %-9.1e", hessSkipped, hessDamped, averageSizingFactor );
          printf("%i, %-9.1e", qpResolve, l1VectorNorm( vars->deltaH )/vars->nBlocks );
        }
        printf("\n");
      }

      if (param->debugLevel > 0 ) {
        // Print everything in a CSV file as well
        fprintf( progressFile, "%23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %i, %i, %i, %23.16e, %i, %23.16e\n",
                 vars->obj, vars->cNormS, vars->tol, vars->gradNorm, lInfVectorNorm( vars->deltaXi ),
                 vars->lambdaStepNorm, vars->alpha, vars->nSOCS, hessSkipped, hessDamped, averageSizingFactor,
                 qpResolve, l1VectorNorm( vars->deltaH )/vars->nBlocks );

        // Print update sequence
        fprintf( updateFile, "%i\t", qpResolve );
      }
    }

    // Print Debug information
    printDebug( vars, param );

    // Do not accidentally print hessSkipped in the next iteration
    hessSkipped = 0;
    hessDamped = 0;

    // qpIterations = number of iterations for the QP that determines the step, can be a resolve (+SOC)
    // qpIterations2 = number of iterations for a QP which solution was discarded
    qpItTotal += qpIterations;
    qpItTotal += qpIterations2;
    qpIterations = 0;
    qpIterations2 = 0;
    qpResolve = 0;

    if (param->printLevel > 0 ) {
      if (hasConverged && vars->steptype < 2 ) {
        if (param->printColor ) {
          printf("\n\033[1;32m***CONVERGENCE ACHIEVED!***\n\033[0m");
        } else {
          printf("\n***CONVERGENCE ACHIEVED!***\n");
        }
      }
    }
  }


  void SQPstats::initStats( SQPoptions *param ) {
    PATHSTR filename;

    // Open files

    if (param->debugLevel > 0 ) {
      // SQP progress
      strcpy( filename, outpath );
      strcat( filename, "sqpits.csv" );
      progressFile = fopen( filename, "w");

      // Update sequence
      strcpy( filename, outpath );
      strcat( filename, "updatesequence.txt" );
      updateFile = fopen( filename, "w" );
    }

    if (param->debugLevel > 1 ) {
      // Primal variables
      strcpy( filename, outpath );
      strcat( filename, "pv.csv" );
      primalVarsFile = fopen( filename, "w");

      // Dual variables
      strcpy( filename, outpath );
      strcat( filename, "dv.csv" );
      dualVarsFile = fopen( filename, "w");
    }

    itCount = 0;
    qpItTotal = 0;
    qpIterations = 0;
    hessSkipped = 0;
    hessDamped = 0;
    averageSizingFactor = 0.0;
  }


  void SQPstats::printPrimalVars( const Matrix &xi ) {
    for (int i=0; i<xi.M()-1; i++ )
      fprintf( primalVarsFile, "%23.16e ", xi( i ) );
    fprintf( primalVarsFile, "%23.16e\n", xi( xi.M()-1 ) );
  }


  void SQPstats::printDualVars( const Matrix &lambda ) {
    for (int i=0; i<lambda.M()-1; i++ )
      fprintf( dualVarsFile, "%23.16e ", lambda( i ) );
    fprintf( dualVarsFile, "%23.16e\n", lambda( lambda.M()-1 ) );
  }


  void SQPstats::printHessian( int nBlocks, SymMatrix *&hess ) {
    PATHSTR filename;
    int offset, i, j, iBlock, nVar;

    nVar = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++ )
      nVar += hess[iBlock].M();

    SymMatrix fullHessian;
    fullHessian.Dimension( nVar ).Initialize( 0.0 );

    strcpy( filename, outpath );
    strcat( filename, "hes.m" );
    hessFile = fopen( filename, "w");

    offset = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++ )
      {
        for (i=0; i<hess[iBlock].N(); i++ )
          for (j=i; j<hess[iBlock].N(); j++ )
            fullHessian( offset + i, offset + j ) = hess[iBlock]( i,j );

        offset += hess[iBlock].N();
      }

    fprintf( hessFile, "H=" );
    fullHessian.Print( hessFile, 23, 1 );
    fprintf( hessFile, "\n" );
    fclose( hessFile );
  }


  void SQPstats::printHessian( int nVar, double *hesNz, int *hesIndRow, int *hesIndCol ) {
    PATHSTR filename;

    strcpy( filename, outpath );
    strcat( filename, "hes.dat" );
    hessFile = fopen( filename, "w");

    printSparseMatlab( hessFile, nVar, nVar, hesNz, hesIndRow, hesIndCol );

    fprintf( hessFile, "\n" );
    fclose( hessFile );
  }


  void SQPstats::printJacobian( const Matrix &constrJac ) {
    PATHSTR filename;

    strcpy( filename, outpath );
    strcat( filename, "jac.m" );
    jacFile = fopen( filename, "w");

    fprintf( jacFile, "A=" );
    constrJac.Print( jacFile, 23, 1 );
    fprintf( jacFile, "\n" );

    fclose( jacFile );
  }


  void SQPstats::printJacobian( int nCon, int nVar, double *jacNz, int *jacIndRow, int *jacIndCol ) {
    PATHSTR filename;

    strcpy( filename, outpath );
    strcat( filename, "jac.dat" );
    jacFile = fopen( filename, "w");

    printSparseMatlab( jacFile, nCon, nVar, jacNz, jacIndRow, jacIndCol );

    fprintf( jacFile, "\n" );
    fclose( jacFile );
  }


  void SQPstats::printSparseMatlab( FILE *file, int nRow, int nCol, double *nz, int *indRow, int *indCol ) {
    int i, j, count;

    count = 0;
    fprintf( file, "%i %i 0\n", nRow, nCol );
    for (i=0; i<nCol; i++ )
      for (j=indCol[i]; j<indCol[i+1]; j++ )
        {
          // +1 for MATLAB indices!
          fprintf( file, "%i %i %23.16e\n", indRow[count]+1, i+1, nz[count] );
          count++;
        }
  }


  void SQPstats::printDebug( SQPiterate *vars, SQPoptions *param ) {
    if (param->debugLevel > 1 )
      {
        printPrimalVars( vars->xi );
        printDualVars( vars->lambda );
      }
  }


  void SQPstats::finish( SQPoptions *param ) {
    if (param->debugLevel > 0 ) {
      fprintf( progressFile, "\n" );
      fclose( progressFile );
      fprintf( updateFile, "\n" );
      fclose( updateFile );
    }

    if (param->debugLevel > 1 ) {
      fclose( primalVarsFile );
      fclose( dualVarsFile );
    }
  }


  void SQPstats::printCppNull( FILE *outfile, char* varname ) {
    fprintf( outfile, "    double *%s = NULL;\n", varname );
  }


  void SQPstats::printVectorCpp( FILE *outfile, double *vec, int len, char* varname ) {
    int i;

    fprintf( outfile, "    double %s[%i] = { ", varname, len );
    for (i=0; i<len; i++ ) {
      fprintf( outfile, "%23.16e", vec[i] );
      if (i != len-1 )
        fprintf( outfile, ", " );
      if ((i+1) % 10 == 0 )
        fprintf( outfile, "\n          " );
    }
    fprintf( outfile, " };\n\n" );
  }


  void SQPstats::printVectorCpp( FILE *outfile, int *vec, int len, char* varname ) {
    int i;

    fprintf( outfile, "    int %s[%i] = { ", varname, len );
    for (i=0; i<len; i++ ) {
      fprintf( outfile, "%i", vec[i] );
      if (i != len-1 )
        fprintf( outfile, ", " );
      if ((i+1) % 15 == 0 )
        fprintf( outfile, "\n          " );
    }
    fprintf( outfile, " };\n\n" );
  }


  void SQPstats::dumpQPCpp( Problemspec *prob, SQPiterate *vars, qpOASES::SQProblem *qp, int sparseQP ) {
    int i, j;
    PATHSTR filename;
    FILE *outfile;
    int n = prob->nVar;
    int m = prob->nCon;

    // Print dimensions
    strcpy( filename, outpath );
    strcat( filename, "qpoases_dim.dat" );
    outfile = fopen( filename, "w" );
    fprintf( outfile, "%i %i\n", n, m );
    fclose( outfile );

    // Print Hessian
    if (sparseQP ) {
      strcpy( filename, outpath );
      strcat( filename, "qpoases_H_sparse.dat" );
      outfile = fopen( filename, "w" );
      for (i=0; i<prob->nVar+1; i++ )
        fprintf( outfile, "%i ", vars->hessIndCol[i] );
      fprintf( outfile, "\n" );

      for (i=0; i<vars->hessIndCol[prob->nVar]; i++ )
        fprintf( outfile, "%i ", vars->hessIndRow[i] );
      fprintf( outfile, "\n" );

      for (i=0; i<vars->hessIndCol[prob->nVar]; i++ )
        fprintf( outfile, "%23.16e ", vars->hessNz[i] );
      fprintf( outfile, "\n" );
      fclose( outfile );
    }
    strcpy( filename, outpath );
    strcat( filename, "qpoases_H.dat" );
    outfile = fopen( filename, "w" );
    int blockCnt = 0;
    for (i=0; i<n; i++ ) {
      for (j=0; j<n; j++ ) {
        if (i == vars->blockIdx[blockCnt+1] ) blockCnt++;
        if (j >= vars->blockIdx[blockCnt] && j < vars->blockIdx[blockCnt+1] ) {
          fprintf( outfile, "%23.16e ", vars->hess[blockCnt]( i - vars->blockIdx[blockCnt], j - vars->blockIdx[blockCnt] ) );
        } else {
          fprintf( outfile, "0.0 " );
        }
      }
      fprintf( outfile, "\n" );
    }
    fclose( outfile );

    // Print gradient
    strcpy( filename, outpath );
    strcat( filename, "qpoases_g.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<n; i++ )
      fprintf( outfile, "%23.16e ", vars->gradObj( i ) );
    fprintf( outfile, "\n" );
    fclose( outfile );

    // Print Jacobian
    strcpy( filename, outpath );
    strcat( filename, "qpoases_A.dat" );
    outfile = fopen( filename, "w" );
    if (sparseQP ) {
      // Always print dense Jacobian
      Matrix constrJacTemp;
      constrJacTemp.Dimension( prob->nCon, prob->nVar ).Initialize( 0.0 );
      for (i=0; i<prob->nVar; i++ )
        for (j=vars->jacIndCol[i]; j<vars->jacIndCol[i+1]; j++ )
          constrJacTemp( vars->jacIndRow[j], i ) = vars->jacNz[j];
      for (i=0; i<m; i++ ) {
        for (j=0; j<n; j++ )
          fprintf( outfile, "%23.16e ", constrJacTemp( i, j ) );
        fprintf( outfile, "\n" );
      }
      fclose( outfile );
    } else {
      for (i=0; i<m; i++ ) {
        for (j=0; j<n; j++ )
          fprintf( outfile, "%23.16e ", vars->constrJac( i, j ) );
        fprintf( outfile, "\n" );
      }
      fclose( outfile );
    }

    if (sparseQP ) {
      strcpy( filename, outpath );
      strcat( filename, "qpoases_A_sparse.dat" );
      outfile = fopen( filename, "w" );
      for (i=0; i<prob->nVar+1; i++ )
        fprintf( outfile, "%i ", vars->jacIndCol[i] );
      fprintf( outfile, "\n" );

      for (i=0; i<vars->jacIndCol[prob->nVar]; i++ )
        fprintf( outfile, "%i ", vars->jacIndRow[i] );
      fprintf( outfile, "\n" );

      for (i=0; i<vars->jacIndCol[prob->nVar]; i++ )
        fprintf( outfile, "%23.16e ", vars->jacNz[i] );
      fprintf( outfile, "\n" );
      fclose( outfile );
    }

    // Print variable lower bounds
    strcpy( filename, outpath );
    strcat( filename, "qpoases_lb.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<n; i++ )
      fprintf( outfile, "%23.16e ", vars->deltaBl( i ) );
    fprintf( outfile, "\n" );
    fclose( outfile );

    // Print variable upper bounds
    strcpy( filename, outpath );
    strcat( filename, "qpoases_ub.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<n; i++ )
      fprintf( outfile, "%23.16e ", vars->deltaBu( i ) );
    fprintf( outfile, "\n" );
    fclose( outfile );

    // Print constraint lower bounds
    strcpy( filename, outpath );
    strcat( filename, "qpoases_lbA.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<m; i++ )
      fprintf( outfile, "%23.16e ", vars->deltaBl( i+n ) );
    fprintf( outfile, "\n" );
    fclose( outfile );

    // Print constraint upper bounds
    strcpy( filename, outpath );
    strcat( filename, "qpoases_ubA.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<m; i++ )
      fprintf( outfile, "%23.16e ", vars->deltaBu( i+n ) );
    fprintf( outfile, "\n" );
    fclose( outfile );

    // Print active set
    qpOASES::Bounds b;
    qpOASES::Constraints c;
    qp->getBounds( b );
    qp->getConstraints( c );

    strcpy( filename, outpath );
    strcat( filename, "qpoases_as.dat" );
    outfile = fopen( filename, "w" );
    for (i=0; i<n; i++ )
      fprintf( outfile, "%i ", b.getStatus( i ) );
    fprintf( outfile, "\n" );
    for (i=0; i<m; i++ )
      fprintf( outfile, "%i ", c.getStatus( i ) );
    fprintf( outfile, "\n" );
    fclose( outfile );
  }

  void SQPstats::dumpQPMatlab( Problemspec *prob, SQPiterate *vars, int sparseQP ) {
    Matrix temp;
    PATHSTR filename;
    FILE *qpFile;
    FILE *vecFile;

    // Print vectors g, lb, lu, lbA, luA
    strcpy( filename, outpath );
    strcat( filename, "vec.m" );
    vecFile = fopen( filename, "w");

    fprintf( vecFile, "g=" );
    vars->gradObj.Print( vecFile, 23, 1 );
    fprintf( vecFile, "\n\n" );

    temp.Submatrix( vars->deltaBl, prob->nVar, 1, 0, 0 );
    fprintf( vecFile, "lb=" );
    temp.Print( vecFile, 23, 1 );
    fprintf( vecFile, "\n\n" );

    temp.Submatrix( vars->deltaBu, prob->nVar, 1, 0, 0 );
    fprintf( vecFile, "lu=" );
    temp.Print( vecFile, 23, 1 );
    fprintf( vecFile, "\n\n" );

    temp.Submatrix( vars->deltaBl, prob->nCon, 1, prob->nVar, 0 );
    fprintf( vecFile, "lbA=" );
    temp.Print( vecFile, 23, 1 );
    fprintf( vecFile, "\n\n" );

    temp.Submatrix( vars->deltaBu, prob->nCon, 1, prob->nVar, 0 );
    fprintf( vecFile, "luA=" );
    temp.Print( vecFile, 23, 1 );
    fprintf( vecFile, "\n" );

    fclose( vecFile );

    // Print sparse Jacobian and Hessian
    if (sparseQP )
      {
        printJacobian( prob->nCon, prob->nVar, vars->jacNz, vars->jacIndRow, vars->jacIndCol );
        printHessian( prob->nVar, vars->hessNz, vars->hessIndRow, vars->hessIndCol );
      }

    // Print a script that correctly reads everything
    strcpy( filename, outpath );
    strcat( filename, "getqp.m" );
    qpFile = fopen( filename, "w");

    fprintf( qpFile, "%% Read vectors g, lb, lu, lbA, luA\n" );
    fprintf( qpFile, "vec;\n" );
    fprintf( qpFile, "%% Read sparse Jacobian\n" );
    fprintf( qpFile, "load jac.dat\n" );
    fprintf( qpFile, "if jac(1) == 0\n" );
    fprintf( qpFile, "    A = [];\n" );
    fprintf( qpFile, "else\n" );
    fprintf( qpFile, "    A = spconvert( jac );\n" );
    fprintf( qpFile, "end\n" );
    fprintf( qpFile, "%% Read sparse Hessian\n" );
    fprintf( qpFile, "load hes.dat\n" );
    fprintf( qpFile, "H = spconvert( hes );\n" );

    fclose( qpFile );
  }

  RestorationProblem::RestorationProblem( Problemspec *parentProblem, const Matrix &xiReference ) {
    int i, iVar, iCon;

    parent = parentProblem;
    xiRef.Dimension( parent->nVar );
    for (i=0; i<parent->nVar; i++)
      xiRef( i ) = xiReference( i );

    /* nCon slack variables */
    nVar = parent->nVar + parent->nCon;
    nCon = parent->nCon;

    /* Block structure: One additional block for every slack variable */
    nBlocks = parent->nBlocks+nCon;
    blockIdx = new int[nBlocks+1];
    for (i=0; i<parent->nBlocks+1; i++ )
      blockIdx[i] = parent->blockIdx[i];
    for (i=parent->nBlocks+1; i<nBlocks+1; i++ )
      blockIdx[i] = blockIdx[i-1]+1;

    /* Set bounds */
    objLo = 0.0;
    objUp = 1.0e20;

    bl.Dimension( nVar + nCon ).Initialize( -1.0e20 );
    bu.Dimension( nVar + nCon ).Initialize( 1.0e20 );
    for (iVar=0; iVar<parent->nVar; iVar++ )
      {
        bl( iVar ) = parent->bl( iVar );
        bu( iVar ) = parent->bu( iVar );
      }

    for (iCon=0; iCon<parent->nCon; iCon++ )
      {
        bl( nVar+iCon ) = parent->bl( parent->nVar+iCon );
        bu( nVar+iCon ) = parent->bu( parent->nVar+iCon );
      }
  }


  void RestorationProblem::evaluate( const Matrix &xi, const Matrix &lambda,
                                     double *objval, Matrix &constr,
                                     Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                                     SymMatrix *&hess, int dmode, int *info ) {
    int iCon, i;
    double diff, regTerm;
    Matrix xiOrig, slack;

    // The first nVar elements of the variable vector correspond to the variables of the original problem
    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );

    // Evaluate constraints of the original problem
    parent->evaluate( xiOrig, lambda, objval, constr,
                      gradObj, jacNz, jacIndRow, jacIndCol, hess, dmode, info );

    // Subtract slacks
    for (iCon=0; iCon<nCon; iCon++ )
      constr( iCon ) -= slack( iCon );


    /* Evaluate objective: minimize slacks plus deviation from reference point */
    if (dmode < 0 )
      return;

    *objval = 0.0;

    // First part: sum of slack variables
    for (i=0; i<nCon; i++ )
      *objval += slack( i ) * slack( i );
    *objval = 0.5 * rho * (*objval);

    // Second part: regularization term
    regTerm = 0.0;
    for (i=0; i<parent->nVar; i++ )
      {
        diff = xiOrig( i ) - xiRef( i );
        regTerm += diagScale( i ) * diff * diff;
      }
    regTerm = 0.5 * zeta * regTerm;
    *objval += regTerm;

    if (dmode > 0 )
      {// compute objective gradient

        // gradient w.r.t. xi (regularization term)
        for (i=0; i<parent->nVar; i++ )
          gradObj( i ) = zeta * diagScale( i ) * diagScale( i ) * (xiOrig( i ) - xiRef( i ));

        // gradient w.r.t. slack variables
        for (i=parent->nVar; i<nVar; i++ )
          gradObj( i ) = rho * xi( i );
      }

    *info = 0;
  }

  void RestorationProblem::evaluate( const Matrix &xi, const Matrix &lambda,
                                     double *objval, Matrix &constr,
                                     Matrix &gradObj, Matrix &constrJac,
                                     SymMatrix *&hess, int dmode, int *info ) {
    int iCon, i;
    double diff, regTerm;
    Matrix xiOrig, constrJacOrig;
    Matrix slack;

    // The first nVar elements of the variable vector correspond to the variables of the original problem
    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );
    if (dmode != 0 )
      constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Evaluate constraints of the original problem
    parent->evaluate( xiOrig, lambda, objval, constr,
                      gradObj, constrJacOrig, hess, dmode, info );

    // Subtract slacks
    for (iCon=0; iCon<nCon; iCon++ )
      constr( iCon ) -= slack( iCon );


    /* Evaluate objective: minimize slacks plus deviation from reference point */
    if (dmode < 0 )
      return;

    *objval = 0.0;

    // First part: sum of slack variables
    for (i=0; i<nCon; i++ )
      *objval += slack( i ) * slack( i );
    *objval = 0.5 * rho * (*objval);

    // Second part: regularization term
    regTerm = 0.0;
    for (i=0; i<parent->nVar; i++ ) {
      diff = xiOrig( i ) - xiRef( i );
      regTerm += diagScale( i ) * diff * diff;
    }
    regTerm = 0.5 * zeta * regTerm;
    *objval += regTerm;

    if (dmode > 0 ) {
      // compute objective gradient

      // gradient w.r.t. xi (regularization term)
      for (i=0; i<parent->nVar; i++ )
        gradObj( i ) = zeta * diagScale( i ) * diagScale( i ) * (xiOrig( i ) - xiRef( i ));

      // gradient w.r.t. slack variables
      for (i=parent->nVar; i<nVar; i++ )
        gradObj( i ) = rho * slack( i );
    }

    *info = 0;
  }


  void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol ) {
    int i, info;
    double objval;
    Matrix xiOrig, slack, constrRef;

    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );

    // Call initialize of the parent problem. There, the sparse Jacobian is allocated
    double *jacNzOrig = NULL;
    int *jacIndRowOrig = NULL, *jacIndColOrig = NULL, nnz, nnzOrig;
    parent->initialize( xiOrig, lambda, jacNzOrig, jacIndRowOrig, jacIndColOrig );
    nnzOrig = jacIndColOrig[parent->nVar];

    // Copy sparse Jacobian from original problem
    nnz = nnzOrig + nCon;
    jacNz = new double[nnz];
    jacIndRow = new int[nnz + (nVar+1)];
    jacIndCol = jacIndRow + nnz;
    for (i=0; i<nnzOrig; i++ )
      {
        jacNz[i] = jacNzOrig[i];
        jacIndRow[i] = jacIndRowOrig[i];
      }
    for (i=0; i<parent->nVar; i++ ) jacIndCol[i] = jacIndColOrig[i];

    // Jacobian entries for slacks (one nonzero entry per column)
    for (i=nnzOrig; i<nnz; i++ ) {
      jacNz[i] = -1.0;
      jacIndRow[i] = i-nnzOrig;
    }
    for (i=parent->nVar; i<nVar+1; i++ )
      jacIndCol[i] = nnzOrig + i - parent->nVar;

    // The reference point is the starting value for the restoration phase
    for (i=0; i<parent->nVar; i++ )
      xiOrig( i ) = xiRef( i );

    // Initialize slack variables such that the constraints are feasible
    constrRef.Dimension( nCon );
    parent->evaluate( xiOrig, &objval, constrRef, &info );

    for (i=0; i<nCon; i++ ) {
      if (constrRef( i ) <= parent->bl( parent->nVar + i ) )// if lower bound is violated
        slack( i ) = constrRef( i ) - parent->bl( parent->nVar + i );
      else if (constrRef( i ) > parent->bu( parent->nVar + i ) )// if upper bound is violated
        slack( i ) = constrRef( i ) - parent->bu( parent->nVar + i );
    }

    // Set diagonal scaling matrix
    diagScale.Dimension( parent->nVar ).Initialize( 1.0 );
    for (i=0; i<parent->nVar; i++ )
      if (fabs( xiRef( i ) ) > 1.0 )
        diagScale( i ) = 1.0 / fabs( xiRef( i ) );

    // Regularization factor zeta and rho \todo wie setzen?
    zeta = 1.0e-3;
    rho = 1.0e3;

    lambda.Initialize( 0.0 );
  }


  void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac ) {
    int i, info;
    double objval;
    Matrix xiOrig, slack, constrJacOrig, constrRef;

    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );
    constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Call initialize of the parent problem to set up linear constraint matrix correctly
    parent->initialize( xiOrig, lambda, constrJacOrig );

    // Jacobian entries for slacks
    for (i=0; i<parent->nCon; i++ )
      constrJac( i, parent->nVar+i ) = -1.0;

    // The reference point is the starting value for the restoration phase
    for (i=0; i<parent->nVar; i++ )
      xiOrig( i ) = xiRef( i );

    // Initialize slack variables such that the constraints are feasible
    constrRef.Dimension( nCon );
    parent->evaluate( xiOrig, &objval, constrRef, &info );

    for (i=0; i<nCon; i++ ) {
      if (constrRef( i ) <= parent->bl( parent->nVar + i ) )// if lower bound is violated
        slack( i ) = constrRef( i ) - parent->bl( parent->nVar + i );
      else if (constrRef( i ) > parent->bu( parent->nVar + i ) )// if upper bound is violated
        slack( i ) = constrRef( i ) - parent->bu( parent->nVar + i );
    }

    // Set diagonal scaling matrix
    diagScale.Dimension( parent->nVar ).Initialize( 1.0 );
    for (i=0; i<parent->nVar; i++ )
      if (fabs( xiRef( i ) ) > 1.0 )
        diagScale( i ) = 1.0 / fabs( xiRef( i ) );

    // Regularization factor zeta and rho \todo wie setzen?
    zeta = 1.0e-3;
    rho = 1.0e3;

    lambda.Initialize( 0.0 );
  }


  void RestorationProblem::printVariables( const Matrix &xi, const Matrix &lambda, int verbose ) {
    int k;

    printf("\n<|----- Original Variables -----|>\n");
    for (k=0; k<parent->nVar; k++ )
      //printf("%7i: %-30s   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, parent->varNames[k], bl(k), xi(k), bu(k), lambda(k));
      printf("%7i: x%-5i   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, k, bl(k), xi(k), bu(k), lambda(k));
    printf("\n<|----- Slack Variables -----|>\n");
    for (k=parent->nVar; k<nVar; k++ )
      printf("%7i: slack   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, bl(k), xi(k), bu(k), lambda(k));
  }


  void RestorationProblem::printConstraints( const Matrix &constr, const Matrix &lambda ) {
    printf("\n<|----- Constraints -----|>\n");
    for (int k=0; k<nCon; k++ )
      //printf("%5i: %-30s   %7g <= %10.4g <= %7g   |   mul=%10.3g\n", k+1, parent->conNames[parent->nVar+k], bl(nVar+k), constr(k), bu(nVar+k), lambda(nVar+k));
      printf("%5i: c%-5i   %7g <= %10.4g <= %7g   |   mul=%10.3g\n", k+1, k, bl(nVar+k), constr(k), bu(nVar+k), lambda(nVar+k));
  }


  void RestorationProblem::printInfo() {
    printf("Minimum 2-norm NLP to find a point acceptable to the filter\n");
  }

} // namespace blocksqp
