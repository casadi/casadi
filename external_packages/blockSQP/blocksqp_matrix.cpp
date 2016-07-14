/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_matrix.cpp
 * \author Dennis Janka, based on VPLAN's matrix.h by Stefan Koerkel
 * \date 2012-2015
 *
 *  Implementation of Matrix and SymMatrix classes for easy access of
 *  matrix elements.
 */

#include "blocksqp.hpp"

namespace blockSQP
{

// #define MATRIX_DEBUG

void Error( const char *F )
{
    printf("Error: %s\n", F );
    //exit( 1 );
}

/* ----------------------------------------------------------------------- */

int Ccount = 0;
int Dcount = 0;
int Ecount = 0;

/* ----------------------------------------------------------------------- */

int Matrix::malloc( void )
{
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


int Matrix::free( void )
{
    if ( tflag )
        Error("free cannot be called with Submatrix");

    if ( array != NULL )
        delete[] array;

    return 0;
}


double &Matrix::operator()( int i, int j )
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
        Error("Invalid matrix entry");
    #endif

    return array[i+j*ldim];
}

double &Matrix::operator()( int i, int j ) const
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
        Error("Invalid matrix entry");
    #endif

    return array[i+j*ldim];
}

double &Matrix::operator()( int i )
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m )
        Error("Invalid matrix entry");
    #endif

    return array[i];
}

double &Matrix::operator()( int i ) const
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m )
        Error("Invalid matrix entry");
    #endif

    return array[i];
}

//double Matrix::a( int i, int j ) const
//{
    //#ifdef MATRIX_DEBUG
    //if ( i < 0 || i >= m || j < 0 || j >= n )
        //Error("Invalid matrix entry");
    //#endif

    //return array[i+j*ldim];
//}


//double Matrix::a( int i ) const
//{
    //#ifdef MATRIX_DEBUG
    //if ( i < 0 || i >= m )
        //Error("Invalid matrix entry");
    //#endif

    //return array[i];
//}


Matrix::Matrix( int M, int N, int LDIM )
{
    Ccount++;

    m = M;
    n = N;
    ldim = LDIM;
    tflag = 0;

    malloc();
}


Matrix::Matrix( int M, int N, double *ARRAY, int LDIM )
{
    Ccount++;

    m = M;
    n = N;
    array = ARRAY;
    ldim = LDIM;
    tflag = 0;

    if ( ldim < m )
        ldim = m;
}


Matrix::Matrix( const Matrix &A )
{
    int i, j;
    //printf("copy constructor\n");
    Ccount++;

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

Matrix &Matrix::operator=( const Matrix &A )
{
    int i, j;
    //printf("assignment operator\n");
    Ecount++;

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


Matrix::~Matrix( void )
{
    Dcount++;

    if ( !tflag )
        free();
}

/* ----------------------------------------------------------------------- */

int Matrix::M( void ) const
{   return m;
}


int Matrix::N( void ) const
{   return n;
}


int Matrix::LDIM( void ) const
{   return ldim;
}


double *Matrix::ARRAY( void ) const
{   return array;
}


int Matrix::TFLAG( void ) const
{   return tflag;
}

/* ----------------------------------------------------------------------- */

Matrix &Matrix::Dimension( int M, int N, int LDIM )
{
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


Matrix &Matrix::Initialize( double (*f)( int, int ) )
{
    int i, j;

    for ( i = 0; i < m; i++ )
        for ( j = 0; j < n; j++ )
            (*this)(i,j) = f(i,j);

    return *this;
}


Matrix &Matrix::Initialize( double val )
{
    int i, j;

    for ( i = 0; i < m; i++ )
        for ( j = 0; j < n; j++ )
            (*this)(i,j) = val;

    return *this;
}


/* ----------------------------------------------------------------------- */

Matrix &Matrix::Submatrix( const Matrix &A, int M, int N, int i0, int j0 )
{
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


Matrix &Matrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM )
{
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


const Matrix &Matrix::Print( FILE *f, int DIGITS, int flag ) const
{    int i, j;
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



int SymMatrix::malloc( void )
{
    int len;

    len = m*(m+1)/2.0;

    if ( len == 0 )
       array = NULL;
    else
       if ( ( array = new double[len] ) == NULL )
          Error("'new' failed");

    return 0;
}


int SymMatrix::free( void )
{
    if( array != NULL )
        delete[] array;

    return 0;
}


double &SymMatrix::operator()( int i, int j )
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
    Error("Invalid matrix entry");
    #endif

    int pos;

    if( i < j )//reference to upper triangular part
        pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
        pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
}


double &SymMatrix::operator()( int i, int j ) const
{
    #ifdef MATRIX_DEBUG
    if ( i < 0 || i >= m || j < 0 || j >= n )
    Error("Invalid matrix entry");
    #endif

    int pos;

    if( i < j )//reference to upper triangular part
        pos = (int) (j + i*(m - (i+1.0)/2.0));
    else
        pos = (int) (i + j*(m - (j+1.0)/2.0));

    return array[pos];
}


double &SymMatrix::operator()( int i )
{
    #ifdef MATRIX_DEBUG
    if ( i >= m*(m+1)/2.0 )
    Error("Invalid matrix entry");
    #endif

    return array[i];
}


double &SymMatrix::operator()( int i ) const
{
    #ifdef MATRIX_DEBUG
    if ( i >= m*(m+1)/2.0 )
    Error("Invalid matrix entry");
    #endif

    return array[i];
}


//double SymMatrix::a( int i, int j ) const
//{
    //#ifdef MATRIX_DEBUG
    //if ( i < 0 || i >= m || j < 0 || j >= n )
    //Error("Invalid matrix entry");
    //#endif

    //int pos;

    //if( j > i )//reference to upper triangular part
        //pos = (int) (j + i*(m - (i+1.0)/2.0));
    //else
        //pos = (int) (i + j*(m - (j+1.0)/2.0));

    //return array[pos];
//}


//double SymMatrix::a( int i ) const
//{
    //#ifdef MATRIX_DEBUG
    //if ( i >= m*(m+1)/2.0 )
    //Error("Invalid matrix entry");
    //#endif

    //return array[i];
//}


SymMatrix::SymMatrix( int M )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
}


SymMatrix::SymMatrix( int M, double *ARRAY )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
}


SymMatrix::SymMatrix( int M, int N, int LDIM )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
}


SymMatrix::SymMatrix( int M, int N, double *ARRAY, int LDIM )
{
    m = M;
    n = M;
    ldim = M;
    tflag = 0;

    malloc();
    array = ARRAY;
}


SymMatrix::SymMatrix( const Matrix &A )
{
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


SymMatrix::SymMatrix( const SymMatrix &A )
{
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


SymMatrix::~SymMatrix( void )
{
    Dcount++;

    if( !tflag )
        free();
}



SymMatrix &SymMatrix::Dimension( int M )
{
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
}


SymMatrix &SymMatrix::Dimension( int M, int N, int LDIM )
{
    free();
    m = M;
    n = M;
    ldim = M;

    malloc();

    return *this;
}


SymMatrix &SymMatrix::Initialize( double (*f)( int, int ) )
{
    int i, j;

    for ( j=0; j<m; j++ )
        for ( i=j; i<n ; i++ )
            (*this)(i,j) = f(i,j);

    return *this;
}


SymMatrix &SymMatrix::Initialize( double val )
{
    int i, j;

    for ( j=0; j<m; j++ )
        for ( i=j; i<n ; i++ )
            (*this)(i,j) = val;

    return *this;
}


SymMatrix &SymMatrix::Submatrix( const Matrix &A, int M, int N, int i0, int j0 )
{
    Error("SymMatrix doesn't support Submatrix");
    return *this;
}


SymMatrix &SymMatrix::Arraymatrix( int M, double *ARRAY )
{
    if( !tflag )
        free();

    tflag = 1;
    m = M;
    n = M;
    ldim = M;
    array = ARRAY;

    return *this;
}


SymMatrix &SymMatrix::Arraymatrix( int M, int N, double *ARRAY, int LDIM )
{
    if( !tflag )
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


double delta( int i, int j )
{    return (i == j) ? 1.0 : 0.0;
}


Matrix Transpose( const Matrix &A )
{
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


Matrix &Transpose( const Matrix &A, Matrix &T )
{
    int i, j;

    T.Dimension( A.N(), A.M() );

    for ( i = 0; i < A.N(); i++ )
        for ( j = 0; j < A.M(); j++ )
            T(i,j) = A(j,i);
            //T(i,j) = A.a(j,i);

    return T;
}

} // namespace blockSQP
