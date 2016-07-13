/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_iterate.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of SQPiterate class that holds all variables that are
 *  updated during one SQP iteration.
 */

#include "blocksqp_options.hpp"
#include "blocksqp_stats.hpp"
#include "blocksqp_iterate.hpp"

namespace blockSQP
{

SQPiterate::SQPiterate( Problemspec* prob, SQPoptions* param, bool full )
{
    int maxblocksize = 1;

    // Set nBlocks structure according to if we use block updates or not
    if( param->blockHess == 0 || prob->nBlocks == 1 )
    {
        nBlocks = 1;
        blockIdx = new int[2];
        blockIdx[0] = 0;
        blockIdx[1] = prob->nVar;
        maxblocksize = prob->nVar;
        param->whichSecondDerv = 0;
    }
    else if( param->blockHess == 2 && prob->nBlocks > 1 )
    {// hybrid strategy: 1 block for constraints, 1 for objective
        nBlocks = 2;
        blockIdx = new int[3];
        blockIdx[0] = 0;
        blockIdx[1] = prob->blockIdx[prob->nBlocks-1];
        blockIdx[2] = prob->nVar;
    }
    else
    {
        nBlocks = prob->nBlocks;
        blockIdx = new int[nBlocks+1];
        for( int k=0; k<nBlocks+1; k++ )
        {
            blockIdx[k] = prob->blockIdx[k];
            if( k > 0 )
                if( blockIdx[k] - blockIdx[k-1] > maxblocksize )
                    maxblocksize = blockIdx[k] - blockIdx[k-1];
        }
    }

    if( param->hessLimMem && param->hessMemsize == 0 )
        param->hessMemsize = maxblocksize;

    allocMin( prob );

    if( !param->sparseQP )
    {
        constrJac.Dimension( prob->nCon, prob->nVar ).Initialize( 0.0 );
        hessNz = new double[prob->nVar*prob->nVar];
    }
    else
        hessNz = NULL;

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

    if( full )
    {
        allocHess( param );
        allocAlg( prob, param );
    }
}


SQPiterate::SQPiterate( const SQPiterate &iter )
{
    int i;

    nBlocks = iter.nBlocks;
    blockIdx = new int[nBlocks+1];
    for( i=0; i<nBlocks+1; i++ )
        blockIdx[i] = iter.blockIdx[i];

    xi = iter.xi;
    lambda = iter.lambda;
    constr = iter.constr;
    gradObj = iter.gradObj;
    gradLagrange = iter.gradLagrange;

    constrJac = iter.constrJac;
    if( iter.jacNz != NULL )
    {
        int nVar = xi.M();
        int nnz = iter.jacIndCol[nVar];

        jacNz = new double[nnz];
        for( i=0; i<nnz; i++ )
            jacNz[i] = iter.jacNz[i];

        jacIndRow = new int[nnz + (nVar+1) + nVar];
        for( i=0; i<nnz + (nVar+1) + nVar; i++ )
            jacIndRow[i] = iter.jacIndRow[i];
        jacIndCol = jacIndRow + nnz;
    }
    else
    {
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
void SQPiterate::allocMin( Problemspec *prob )
{
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


void SQPiterate::allocHess( SQPoptions *param )
{
    int iBlock, varDim;

    // Create one Matrix for one diagonal block in the Hessian
    hess1 = new SymMatrix[nBlocks];
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        varDim = blockIdx[iBlock+1] - blockIdx[iBlock];
        hess1[iBlock].Dimension( varDim ).Initialize( 0.0 );
    }

    // For SR1 or finite differences, maintain two Hessians
    if( param->hessUpdate == 1 || param->hessUpdate == 4 )
    {
        hess2 = new SymMatrix[nBlocks];
        for( iBlock=0; iBlock<nBlocks; iBlock++ )
        {
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
void SQPiterate::convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_ )
{
    if( hessNz == NULL )
        return;
    int count = 0;
    int blockCnt = 0;
    for( int i=0; i<prob->nVar; i++ )
        for( int j=0; j<prob->nVar; j++ )
        {
            if( i == blockIdx[blockCnt+1] )
                blockCnt++;
            if( j >= blockIdx[blockCnt] && j < blockIdx[blockCnt+1] )
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
                                 double *&hessNz_, int *&hessIndRow_, int *&hessIndCol_, int *&hessIndLo_ )
{
    int iBlock, count, colCountTotal, rowOffset, i, j;
    int nnz, nCols, nRows;

    // 1) count nonzero elements
    nnz = 0;
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
        for( i=0; i<hess_[iBlock].N(); i++ )
            for( j=i; j<hess_[iBlock].N(); j++ )
                if( fabs(hess_[iBlock]( i,j )) > eps )
                {
                    nnz++;
                    if( i != j )// off-diagonal elements count twice
                        nnz++;
                }

    if( hessNz_ != NULL ) delete[] hessNz_;
    if( hessIndRow_ != NULL ) delete[] hessIndRow_;

    hessNz_ = new double[nnz];
    hessIndRow_ = new int[nnz + (prob->nVar+1) + prob->nVar];
    hessIndCol_ = hessIndRow_ + nnz;
    hessIndLo_ = hessIndCol_ + (prob->nVar+1);

    // 2) store matrix entries columnwise in hessNz
    count = 0; // runs over all nonzero elements
    colCountTotal = 0; // keep track of position in large matrix
    rowOffset = 0;
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        nCols = hess_[iBlock].N();
        nRows = hess_[iBlock].M();

        for( i=0; i<nCols; i++ )
        {
            // column 'colCountTotal' starts at element 'count'
            hessIndCol_[colCountTotal] = count;

            for( j=0; j<nRows; j++ )
                if( fabs(hess_[iBlock]( i,j )) > eps )
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
    for( j=0; j<prob->nVar; j++ )
    {
        for( i=hessIndCol_[j]; i<hessIndCol_[j+1] && hessIndRow_[i]<j; i++);
        hessIndLo_[j] = i;
    }

    if( count != nnz )
         printf( "Error in convertHessian: %i elements processed, should be %i elements!\n", count, nnz );
 }


/**
 * Allocate memory for additional variables
 * needed by the algorithm
 */
void SQPiterate::allocAlg( Problemspec *prob, SQPoptions *param )
{
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
    for( iBlock=0; iBlock<nBlocks; iBlock++ )
        noUpdateCounter[iBlock] = -1;

    // For selective sizing: for each block save sTs, sTs_, sTy, sTy_
    deltaNorm.Dimension( nBlocks ).Initialize( 1.0 );
    deltaNormOld.Dimension( nBlocks ).Initialize( 1.0 );
    deltaGamma.Dimension( nBlocks ).Initialize( 0.0 );
    deltaGammaOld.Dimension( nBlocks ).Initialize( 0.0 );
}


void SQPiterate::initIterate( SQPoptions* param )
{
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

SQPiterate::~SQPiterate( void )
{
    if( blockIdx != NULL )
        delete[] blockIdx;
    if( noUpdateCounter != NULL )
        delete[] noUpdateCounter;
    if( jacNz != NULL )
        delete[] jacNz;
    if( jacIndRow != NULL )
        delete[] jacIndRow;
    if( hessNz != NULL )
        delete[] hessNz;
    if( hessIndRow != NULL )
        delete[] hessIndRow;
}

} // namespace blockSQP
