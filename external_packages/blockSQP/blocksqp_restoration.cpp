/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_restoration.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of RestorationProblem class that describes a
 *  minimum l_2-norm NLP.
 */

#include "blocksqp.hpp"

namespace blockSQP
{

RestorationProblem::RestorationProblem( Problemspec *parentProblem, const Matrix &xiReference )
{
    int i, iVar, iCon;

    parent = parentProblem;
    xiRef.Dimension( parent->nVar );
    for( i=0; i<parent->nVar; i++)
        xiRef( i ) = xiReference( i );

    /* nCon slack variables */
    nVar = parent->nVar + parent->nCon;
    nCon = parent->nCon;

    /* Block structure: One additional block for every slack variable */
    nBlocks = parent->nBlocks+nCon;
    blockIdx = new int[nBlocks+1];
    for( i=0; i<parent->nBlocks+1; i++ )
        blockIdx[i] = parent->blockIdx[i];
    for( i=parent->nBlocks+1; i<nBlocks+1; i++ )
        blockIdx[i] = blockIdx[i-1]+1;

    /* Set bounds */
    objLo = 0.0;
    objUp = 1.0e20;

    bl.Dimension( nVar + nCon ).Initialize( -1.0e20 );
    bu.Dimension( nVar + nCon ).Initialize( 1.0e20 );
    for( iVar=0; iVar<parent->nVar; iVar++ )
    {
        bl( iVar ) = parent->bl( iVar );
        bu( iVar ) = parent->bu( iVar );
    }

    for( iCon=0; iCon<parent->nCon; iCon++ )
    {
        bl( nVar+iCon ) = parent->bl( parent->nVar+iCon );
        bu( nVar+iCon ) = parent->bu( parent->nVar+iCon );
    }
}


void RestorationProblem::evaluate( const Matrix &xi, const Matrix &lambda,
                                   double *objval, Matrix &constr,
                                   Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                                   SymMatrix *&hess, int dmode, int *info )
{
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
    for( iCon=0; iCon<nCon; iCon++ )
        constr( iCon ) -= slack( iCon );


    /* Evaluate objective: minimize slacks plus deviation from reference point */
    if( dmode < 0 )
        return;

    *objval = 0.0;

    // First part: sum of slack variables
    for( i=0; i<nCon; i++ )
        *objval += slack( i ) * slack( i );
    *objval = 0.5 * rho * (*objval);

    // Second part: regularization term
    regTerm = 0.0;
    for( i=0; i<parent->nVar; i++ )
    {
        diff = xiOrig( i ) - xiRef( i );
        regTerm += diagScale( i ) * diff * diff;
    }
    regTerm = 0.5 * zeta * regTerm;
    *objval += regTerm;

    if( dmode > 0 )
    {// compute objective gradient

        // gradient w.r.t. xi (regularization term)
        for( i=0; i<parent->nVar; i++ )
            gradObj( i ) = zeta * diagScale( i ) * diagScale( i ) * (xiOrig( i ) - xiRef( i ));

        // gradient w.r.t. slack variables
        for( i=parent->nVar; i<nVar; i++ )
            gradObj( i ) = rho * xi( i );
    }

    *info = 0;
}

void RestorationProblem::evaluate( const Matrix &xi, const Matrix &lambda,
                                   double *objval, Matrix &constr,
                                   Matrix &gradObj, Matrix &constrJac,
                                   SymMatrix *&hess, int dmode, int *info )
{
    int iCon, i;
    double diff, regTerm;
    Matrix xiOrig, constrJacOrig;
    Matrix slack;

    // The first nVar elements of the variable vector correspond to the variables of the original problem
    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );
    if( dmode != 0 )
        constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Evaluate constraints of the original problem
    parent->evaluate( xiOrig, lambda, objval, constr,
                      gradObj, constrJacOrig, hess, dmode, info );

    // Subtract slacks
    for( iCon=0; iCon<nCon; iCon++ )
        constr( iCon ) -= slack( iCon );


    /* Evaluate objective: minimize slacks plus deviation from reference point */
    if( dmode < 0 )
        return;

    *objval = 0.0;

    // First part: sum of slack variables
    for( i=0; i<nCon; i++ )
        *objval += slack( i ) * slack( i );
    *objval = 0.5 * rho * (*objval);

    // Second part: regularization term
    regTerm = 0.0;
    for( i=0; i<parent->nVar; i++ )
    {
        diff = xiOrig( i ) - xiRef( i );
        regTerm += diagScale( i ) * diff * diff;
    }
    regTerm = 0.5 * zeta * regTerm;
    *objval += regTerm;

    if( dmode > 0 )
    {// compute objective gradient

        // gradient w.r.t. xi (regularization term)
        for( i=0; i<parent->nVar; i++ )
            gradObj( i ) = zeta * diagScale( i ) * diagScale( i ) * (xiOrig( i ) - xiRef( i ));

        // gradient w.r.t. slack variables
        for( i=parent->nVar; i<nVar; i++ )
            gradObj( i ) = rho * slack( i );
    }

    *info = 0;
}


void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol )
{
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
    for( i=0; i<nnzOrig; i++ )
    {
        jacNz[i] = jacNzOrig[i];
        jacIndRow[i] = jacIndRowOrig[i];
    }
    for( i=0; i<parent->nVar; i++ )
        jacIndCol[i] = jacIndColOrig[i];

    // Jacobian entries for slacks (one nonzero entry per column)
    for( i=nnzOrig; i<nnz; i++ )
    {
        jacNz[i] = -1.0;
        jacIndRow[i] = i-nnzOrig;
    }
    for( i=parent->nVar; i<nVar+1; i++ )
        jacIndCol[i] = nnzOrig + i - parent->nVar;

    // The reference point is the starting value for the restoration phase
    for( i=0; i<parent->nVar; i++ )
        xiOrig( i ) = xiRef( i );

    // Initialize slack variables such that the constraints are feasible
    constrRef.Dimension( nCon );
    parent->evaluate( xiOrig, &objval, constrRef, &info );

    for( i=0; i<nCon; i++ )
    {
        if( constrRef( i ) <= parent->bl( parent->nVar + i ) )// if lower bound is violated
            slack( i ) = constrRef( i ) - parent->bl( parent->nVar + i );
        else if( constrRef( i ) > parent->bu( parent->nVar + i ) )// if upper bound is violated
            slack( i ) = constrRef( i ) - parent->bu( parent->nVar + i );
    }

    // Set diagonal scaling matrix
    diagScale.Dimension( parent->nVar ).Initialize( 1.0 );
    for( i=0; i<parent->nVar; i++ )
        if( fabs( xiRef( i ) ) > 1.0 )
            diagScale( i ) = 1.0 / fabs( xiRef( i ) );

    // Regularization factor zeta and rho \todo wie setzen?
    zeta = 1.0e-3;
    rho = 1.0e3;

    lambda.Initialize( 0.0 );
}


void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac )
{
    int i, info;
    double objval;
    Matrix xiOrig, slack, constrJacOrig, constrRef;

    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );
    constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Call initialize of the parent problem to set up linear constraint matrix correctly
    parent->initialize( xiOrig, lambda, constrJacOrig );

    // Jacobian entries for slacks
    for( i=0; i<parent->nCon; i++ )
        constrJac( i, parent->nVar+i ) = -1.0;

    // The reference point is the starting value for the restoration phase
    for( i=0; i<parent->nVar; i++ )
        xiOrig( i ) = xiRef( i );

    // Initialize slack variables such that the constraints are feasible
    constrRef.Dimension( nCon );
    parent->evaluate( xiOrig, &objval, constrRef, &info );

    for( i=0; i<nCon; i++ )
    {
        if( constrRef( i ) <= parent->bl( parent->nVar + i ) )// if lower bound is violated
            slack( i ) = constrRef( i ) - parent->bl( parent->nVar + i );
        else if( constrRef( i ) > parent->bu( parent->nVar + i ) )// if upper bound is violated
            slack( i ) = constrRef( i ) - parent->bu( parent->nVar + i );
    }

    // Set diagonal scaling matrix
    diagScale.Dimension( parent->nVar ).Initialize( 1.0 );
    for( i=0; i<parent->nVar; i++ )
        if( fabs( xiRef( i ) ) > 1.0 )
            diagScale( i ) = 1.0 / fabs( xiRef( i ) );

    // Regularization factor zeta and rho \todo wie setzen?
    zeta = 1.0e-3;
    rho = 1.0e3;

    lambda.Initialize( 0.0 );
}


void RestorationProblem::printVariables( const Matrix &xi, const Matrix &lambda, int verbose )
{
    int k;

    printf("\n<|----- Original Variables -----|>\n");
    for( k=0; k<parent->nVar; k++ )
        //printf("%7i: %-30s   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, parent->varNames[k], bl(k), xi(k), bu(k), lambda(k));
        printf("%7i: x%-5i   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, k, bl(k), xi(k), bu(k), lambda(k));
    printf("\n<|----- Slack Variables -----|>\n");
    for( k=parent->nVar; k<nVar; k++ )
        printf("%7i: slack   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, bl(k), xi(k), bu(k), lambda(k));
}


void RestorationProblem::printConstraints( const Matrix &constr, const Matrix &lambda )
{
    printf("\n<|----- Constraints -----|>\n");
    for( int k=0; k<nCon; k++ )
        //printf("%5i: %-30s   %7g <= %10.4g <= %7g   |   mul=%10.3g\n", k+1, parent->conNames[parent->nVar+k], bl(nVar+k), constr(k), bu(nVar+k), lambda(nVar+k));
        printf("%5i: c%-5i   %7g <= %10.4g <= %7g   |   mul=%10.3g\n", k+1, k, bl(nVar+k), constr(k), bu(nVar+k), lambda(nVar+k));
}


void RestorationProblem::printInfo()
{
    printf("Minimum 2-norm NLP to find a point acceptable to the filter\n");
}

} // namespace blockSQP
