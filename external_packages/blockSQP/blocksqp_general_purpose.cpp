/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_general_purpose.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of general purpose routines for matrix and vector computations.
 */

#include "blocksqp.hpp"

namespace blockSQP
{

/**
 * Compute the inverse of a matrix
 * using LU decomposition (DGETRF and DGETRI)
 */
int inverse( const Matrix &A, Matrix &Ainv )
{
    int i, j;
    int n, ldim, lwork, info = 0;
    int *ipiv;
    double *work;

    for( i=0; i<A.N(); i++ )
        for( j=0; j<A.M(); j++ )
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
int calcEigenvalues( const SymMatrix &B, Matrix &ev )
{
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
    for( i=0; i<dim; i++ )
    {
        radius = 0.0;
        for( j=0; j<dim; j++ )
            if( j != i )
                radius += fabs( B( i,j ) );

        if( B( i,i ) - radius < lambdaMin )
            lambdaMin = B( i,i ) - radius;
    }

    return lambdaMin;
}


/**
 * Compute scalar product aTb
 */
double adotb( const Matrix &a, const Matrix &b )
{
    double norm = 0.0;

    if( a.N() != 1 || b.N() != 1 )
    {
        printf("a or b is not a vector!\n");
    }
    else if( a.M() != b.M() )
    {
        printf("a and b must have the same dimension!\n");
    }
    else
    {
        for( int k=0; k<a.M(); k++ )
            norm += a(k) * b(k);
    }

    return norm;
}

/**
 * Compute the matrix vector product for a column-compressed sparse matrix A with a vector b and store it in result
 */
void Atimesb( double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result )
{
    int nCol = b.M();
    int nRow = result.M();
    int i, k;

    for( i=0; i<nRow; i++ )
        result( i ) = 0.0;

    for( i=0; i<nCol; i++ )
    {
        // k runs over all elements in one column
        for( k=AIndCol[i]; k<AIndCol[i+1]; k++ )
            result( AIndRow[k] ) += Anz[k] * b( i );
    }

}

/**
 * Compute the matrix vector product A*b and store it in result
 */
void Atimesb( const Matrix &A, const Matrix &b, Matrix &result )
{
    result.Initialize( 0.0 );
    for( int i=0; i<A.M(); i++ )
        for( int k=0; k<A.N(); k++ )
            result( i ) += A( i, k ) * b( k );
}

double l1VectorNorm( const Matrix &v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            norm += fabs(v( k ));
    }

    return norm;
}

double l2VectorNorm( const Matrix &v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            norm += v( k )* v( k );
    }

    return sqrt(norm);
}

double lInfVectorNorm( const Matrix &v )
{
    double norm = 0.0;

    if( v.N() != 1 )
    {
        printf("v is not a vector!\n");
    }
    else
    {
        for( int k=0; k<v.M(); k++ )
            if( fabs(v( k )) > norm )
                norm = fabs(v( k ));
    }

    return norm;
}


/**
 * Calculate weighted l1 norm of constraint violations
 */
double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl, const Matrix &weights )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    if( weights.M() < constr.M() + nVar )
    {
        printf("Weight vector too short!\n");
        return 0.0;
    }

    // Weighted violation of simple bounds
    for( i=0; i<nVar; i++ )
    {
        if( xi( i ) > bu( i ) )
            norm += fabs(weights( i )) * (xi( i ) - bu( i ));
        else if( xi( i ) < bl( i ) )
            norm += fabs(weights( i )) * (bl( i ) - xi( i ));
    }

    // Calculate weighted sum of constraint violations
    for( i=0; i<constr.M(); i++ )
    {
        if( constr( i ) > bu( nVar+i ) )
            norm += fabs(weights( nVar+i )) * (constr( i ) - bu( nVar+i ));
        else if( constr( i ) < bl( nVar+i ) )
            norm += fabs(weights( nVar+i )) * (bl( nVar+i ) - constr( i ));
    }

    return norm;
}


/**
 * Calculate l1 norm of constraint violations
 */
double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for( i=0; i<nVar; i++ )
    {
        if( xi( i ) > bu( i ) )
            norm += xi( i ) - bu( i );
        else if( xi( i ) < bl( i ) )
            norm += bl( i ) - xi( i );
    }

    // Calculate sum of constraint violations
    for( i=0; i<constr.M(); i++ )
    {
        if( constr( i ) > bu( nVar+i ) )
            norm += constr( i ) - bu( nVar+i );
        else if( constr( i ) < bl( nVar+i ) )
            norm += bl( nVar+i ) - constr( i );
    }

    return norm;
}


/**
 * Calculate l2 norm of constraint violations
 */
double l2ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();

    // Violation of simple bounds
    for( i=0; i<nVar; i++ )
        if( xi( i ) > bu( i ) )
            norm += xi( i ) - bu( i );
        if( xi( i ) < bl( i ) )
            norm += bl( i ) - xi( i );

    // Calculate sum of constraint violations
    for( i=0; i<constr.M(); i++ )
        if( constr( i ) > bu( nVar+i ) )
            norm += pow(constr( i ) - bu( nVar+i ), 2);
        else if( constr( i ) < bl( nVar+i ) )
            norm += pow(bl( nVar+i ) - constr( i ), 2);

    return sqrt(norm);
}


/**
 * Calculate l_Infinity norm of constraint violations
 */
double lInfConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl )
{
    double norm = 0.0;
    int i;
    int nVar = xi.M();
    int nCon = constr.M();

    // Violation of simple bounds
    for( i=0; i<nVar; i++ )
    {
        if( xi( i ) - bu( i ) > norm )
            norm = xi( i ) - bu( i );
        else if( bl( i ) - xi( i ) > norm )
            norm = bl( i ) - xi( i );
    }

    // Find out the largest constraint violation
    for( i=0; i<nCon; i++ )
    {
        if( constr( i ) - bu( nVar+i ) > norm )
            norm = constr( i ) - bu( nVar+i );
        if( bl( nVar+i ) - constr( i ) > norm )
            norm = bl( nVar+i ) - constr( i );
    }

    return norm;
}

} // namespace blockSQP
