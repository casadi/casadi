/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_hess.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of methods of SQPmethod class associated with
 *  computation of Hessian approximations.
 *
 */

#include "blocksqp.hpp"

namespace blockSQP
{

/**
 * Initial Hessian: Identity matrix
 */
void SQPmethod::calcInitialHessian()
{
    int iBlock;

    for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        //if objective derv is computed exactly, don't set the last block!
        if( !(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks-1) )
            calcInitialHessian( iBlock );
}


/**
 * Initial Hessian for one block: Identity matrix
 */
void SQPmethod::calcInitialHessian( int iBlock )
{
    vars->hess[iBlock].Initialize( 0.0 );

    // Each block is a diagonal matrix
    for( int i=0; i<vars->hess[iBlock].M(); i++ )
        vars->hess[iBlock]( i, i ) = param->iniHessDiag;

    // If we maintain 2 Hessians, also reset the second one
    if( vars->hess2 != NULL )
    {
        vars->hess2[iBlock].Initialize( 0.0 );
        for( int i=0; i<vars->hess2[iBlock].M(); i++ )
            vars->hess2[iBlock]( i, i ) = param->iniHessDiag;
    }
}


void SQPmethod::resetHessian()
{
    for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        //if objective derv is computed exactly, don't set the last block!
        if( !(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks - 1) )
            resetHessian( iBlock );
}


void SQPmethod::resetHessian( int iBlock )
{
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
int SQPmethod::calcFiniteDiffHessian()
{
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
    for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        if( vars->blockIdx[iBlock+1] - vars->blockIdx[iBlock] > maxBlock )
            maxBlock = vars->blockIdx[iBlock+1] - vars->blockIdx[iBlock];

    // Compute original Lagrange gradient
    calcLagrangeGradient( vars->lambda, vars->gradObj, vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->gradLagrange, 0 );

    for( iVar = 0; iVar<maxBlock; iVar++ )
    {
        pert.Initialize( 0.0 );

        // Perturb all blocks simultaneously
        for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        {
            idx = vars->blockIdx[iBlock] + iVar;
            // Skip blocks that have less than iVar variables
            if( idx < vars->blockIdx[iBlock+1] )
            {
                pert( idx ) = myDelta * fabs( vars->xi( idx ) );
                pert( idx ) = fmax( pert( idx ), minDelta );

                // If perturbation violates upper bound, try to perturb with negative
                upperVio = vars->xi( idx ) + pert( idx ) - prob->bu( idx );
                if( upperVio > 0 )
                {
                    lowerVio = prob->bl( idx ) -  ( vars->xi( idx ) - pert( idx ) );
                    // If perturbation violates also lower bound, take the largest perturbation possible
                    if( lowerVio > 0 )
                    {
                        if( lowerVio > upperVio )
                            pert( idx ) = -lowerVio;
                        else
                            pert( idx ) = upperVio;
                    }
                    // If perturbation does not violate lower bound, take -computed perturbation
                    else
                    {
                        pert( idx ) = -pert( idx );
                    }
                }
            }
        }

        // Add perturbation
        for( k=0; k<prob->nVar; k++ )
            vars->xi( k ) += pert( k );

        // Compute perturbed Lagrange gradient
        if( param->sparseQP )
        {
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP.constr, varsP.gradObj,
                            varsP.jacNz, varsP.jacIndRow, varsP.jacIndCol, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP.gradObj, varsP.jacNz, varsP.jacIndRow,
                                  varsP.jacIndCol, varsP.gradLagrange, 0 );
        }
        else
        {
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP.constr, varsP.gradObj, varsP.constrJac, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP.gradObj, varsP.constrJac, varsP.gradLagrange, 0 );
        }

        // Compute finite difference approximations: one column in every block
        for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        {
            idx1 = vars->blockIdx[iBlock] + iVar;
            // Skip blocks that have less than iVar variables
            if( idx1 < vars->blockIdx[iBlock+1] )
            {
                for( jVar=iVar; jVar<vars->blockIdx[iBlock+1]-vars->blockIdx[iBlock]; jVar++ )
                {// Take symmetrized matrices
                    idx2 = vars->blockIdx[iBlock] + jVar;
                    vars->hess[iBlock]( iVar, jVar ) =  ( varsP.gradLagrange( idx1 ) - vars->gradLagrange( idx2 ) );
                    vars->hess[iBlock]( iVar, jVar ) += ( varsP.gradLagrange( idx2 ) - vars->gradLagrange( idx1 ) );
                    vars->hess[iBlock]( iVar, jVar ) *= 0.5 / pert( idx1 );
                }
            }
        }

        // Subtract perturbation
        for( k=0; k<prob->nVar; k++ )
            vars->xi( k ) -= pert( k );
    }

    return info;
}


void SQPmethod::sizeInitialHessian( const Matrix &gamma, const Matrix &delta, int iBlock, int option )
{
    int i, j;
    double scale;
    double myEps = 1.0e3 * param->eps;

    if( option == 1 )
    {// Shanno-Phua
        scale = adotb( gamma, gamma ) / fmax( adotb( delta, gamma ), myEps );
    }
    else if( option == 2 )
    {// Oren-Luenberger
        scale = adotb( delta, gamma ) / fmax( adotb( delta, delta ), myEps );
        scale = fmin( scale, 1.0 );
    }
    else if( option == 3 )
    {// Geometric mean of 1 and 2
        scale = sqrt( adotb( gamma, gamma ) / fmax( adotb( delta, delta ), myEps ) );
    }
    else
    {// Invalid option, ignore
        return;
    }

    if( scale > 0.0 )
    {
        scale = fmax( scale, myEps );
        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= scale;
    }
    else
        scale = 1.0;

    // statistics: average sizing factor
    stats->averageSizingFactor += scale;
}


void SQPmethod::sizeHessianCOL( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j;
    double theta, scale, myEps = 1.0e3 * param->eps;
    double deltaNorm, deltaNormOld, deltaGamma, deltaGammaOld, deltaBdelta;

    // Get sTs, sTs_, sTy, sTy_, sTBs
    deltaNorm = vars->deltaNorm(iBlock);
    deltaGamma = vars->deltaGamma(iBlock);
    deltaNormOld = vars->deltaNormOld(iBlock);
    deltaGammaOld = vars->deltaGammaOld(iBlock);
    deltaBdelta = 0.0;
    for( i=0; i<delta.M(); i++ )
        for( j=0; j<delta.M(); j++ )
            deltaBdelta += delta( i ) * vars->hess[iBlock]( i, j ) * delta( j );

    // Centered Oren-Luenberger factor
    if( vars->noUpdateCounter[iBlock] == -1 ) // in the first iteration, this should equal the OL factor
        theta = 1.0;
    else
        theta = fmin( param->colTau1, param->colTau2 * deltaNorm );
    if( deltaNorm > myEps && deltaNormOld > myEps )
    {
        scale = (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaBdelta / deltaNorm;
        if( scale > param->eps )
            scale = ( (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaGamma / deltaNorm ) / scale;
    }
    else
        scale = 1.0;

    // Size only if factor is between zero and one
    if( scale < 1.0 && scale > 0.0 )
    {
        scale = fmax( param->colEps, scale );
        //printf("Sizing value (COL) block %i = %g\n", iBlock, scale );
        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= scale;

        // statistics: average sizing factor
        stats->averageSizingFactor += scale;
    }
    else
        stats->averageSizingFactor += 1.0;
}

/**
 * Apply BFGS or SR1 update blockwise and size blocks
 */
void SQPmethod::calcHessianUpdate( int updateType, int hessScaling )
{
    int iBlock, nBlocks;
    int nVarLocal;
    Matrix smallGamma, smallDelta;
    bool firstIter;

    //if objective derv is computed exactly, don't set the last block!
    if( param->whichSecondDerv == 1 && param->blockHess )
        nBlocks = vars->nBlocks - 1;
    else
        nBlocks = vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->averageSizingFactor = 0.0;

    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
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
        if( hessScaling < 4 && firstIter )
            sizeInitialHessian( smallGamma, smallDelta, iBlock, hessScaling );
        else if( hessScaling == 4 )
            sizeHessianCOL( smallGamma, smallDelta, iBlock );

        // Compute the new update
        if( updateType == 1 )
        {
            calcSR1( smallGamma, smallDelta, iBlock );

            // Prepare to compute fallback update as well
            vars->hess = vars->hess2;

            // Sizing the fallback update
            if( param->fallbackScaling < 4 && firstIter )
                sizeInitialHessian( smallGamma, smallDelta, iBlock, param->fallbackScaling );
            else if( param->fallbackScaling == 4 )
                sizeHessianCOL( smallGamma, smallDelta, iBlock );

            // Compute fallback update
            if( param->fallbackUpdate == 2 )
                calcBFGS( smallGamma, smallDelta, iBlock );

            // Reset pointer
            vars->hess = vars->hess1;
        }
        else if( updateType == 2 )
            calcBFGS( smallGamma, smallDelta, iBlock );

        // If an update is skipped to often, reset Hessian block
        if( vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates )
            resetHessian( iBlock );
    }

    // statistics: average sizing factor
    stats->averageSizingFactor /= nBlocks;
}


void SQPmethod::calcHessianUpdateLimitedMemory( int updateType, int hessScaling )
{
    int iBlock, nBlocks, nVarLocal;
    Matrix smallGamma, smallDelta;
    Matrix gammai, deltai;
    int i, m, pos, posOldest, posNewest;
    int hessDamped, hessSkipped;
    double averageSizingFactor;

    //if objective derv is computed exactly, don't set the last block!
    if( param->whichSecondDerv == 1 && param->blockHess )
        nBlocks = vars->nBlocks - 1;
    else
        nBlocks = vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->hessSkipped = 0;
    stats->averageSizingFactor = 0.0;

    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        nVarLocal = vars->hess[iBlock].M();

        // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
        // i.e. subvectors of gamma and delta from m prev. iterations
        smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
        smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

        // Memory structure
        if( stats->itCount > smallGamma.N() )
        {
            m = smallGamma.N();
            posOldest = stats->itCount % m;
            posNewest = (stats->itCount-1) % m;
        }
        else
        {
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

        for( i=0; i<m; i++ )
        {
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
            if( hessScaling == 4 )
                sizeHessianCOL( gammai, deltai, iBlock );

            // Compute the new update
            if( updateType == 1 )
                calcSR1( gammai, deltai, iBlock );
            else if( updateType == 2 )
                calcBFGS( gammai, deltai, iBlock );

            stats->nTotalUpdates++;

            // Count damping statistics only for the most recent update
            if( pos != posNewest )
            {
                stats->hessDamped = hessDamped;
                stats->hessSkipped = hessSkipped;
                if( hessScaling == 4 )
                    stats->averageSizingFactor = averageSizingFactor;
            }
        }

        // If an update is skipped to often, reset Hessian block
        if( vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates )
            resetHessian( iBlock );
    }//blocks
    stats->averageSizingFactor /= nBlocks;
}


void SQPmethod::calcBFGS( const Matrix &gamma, const Matrix &delta, int iBlock )
{
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
    for( i=0; i<dim; i++ )
    {
        for( k=0; k<dim; k++ )
            Bdelta( i ) += (*B)( i,k ) * delta( k );

        h1 += delta( i ) * Bdelta( i );
        //h2 += delta( i ) * gamma( i );
    }
    h2 = vars->deltaGamma( iBlock );

    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    damped = 0;
    if( param->hessDamp )
        if( h2 < param->hessDampFac * h1 / vars->alpha && fabs( h1 - h2 ) > 1.0e-12 )
        {// At the first iteration h1 and h2 are equal due to COL scaling

            thetaPowell = (1.0 - param->hessDampFac)*h1 / ( h1 - h2 );

            // Redefine gamma and h2 = delta^T * gamma
            h2 = 0.0;
            for( i=0; i<dim; i++ )
            {
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
    if( fabs( h1 ) < myEps || fabs( h2 ) < myEps )
    {// don't perform update because of bad condition, might introduce negative eigenvalues
        //printf("block = %i, h1 = %g, h2 = %g\n", iBlock, h1, h2);
        vars->noUpdateCounter[iBlock]++;
        stats->hessDamped -= damped;
        stats->hessSkipped++;
        stats->nTotalSkippedUpdates++;
    }
    else
    {
        for( i=0; i<dim; i++ )
            for( j=i; j<dim; j++ )
                (*B)( i,j ) = (*B)( i,j ) - Bdelta( i ) * Bdelta( j ) / h1
                                          + gamma2( i ) * gamma2( j ) / h2;

        vars->noUpdateCounter[iBlock] = 0;
    }
}


void SQPmethod::calcSR1( const Matrix &gamma, const Matrix &delta, int iBlock )
{
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
    for( i=0; i<dim; i++ )
    {
        gmBdelta( i ) = gamma( i );
        for( k=0; k<dim; k++ )
            gmBdelta( i ) -= ( (*B)( i,k ) * delta( k ) );

        h += ( gmBdelta( i ) * delta( i ) );
    }

    // B_k+1 = B_k + gmBdelta * gmBdelta^T / h
    if( fabs( h ) < r * l2VectorNorm( delta ) * l2VectorNorm( gmBdelta ) || fabs( h ) < myEps )
    {// Skip update if denominator is too small
        //printf("block %i, h = %23.16e\n", iBlock, h );
        vars->noUpdateCounter[iBlock]++;
        stats->hessSkipped++;
        stats->nTotalSkippedUpdates++;
    }
    else
    {
        for( i=0; i<dim; i++ )
            for( j=i; j<dim; j++ )
                (*B)( i,j ) = (*B)( i,j ) + gmBdelta( i ) * gmBdelta( j ) / h;
        vars->noUpdateCounter[iBlock] = 0;
    }
}


/**
 * Set deltaXi and gamma as a column in the matrix containing
 * the m most recent delta and gamma
 */
void SQPmethod::updateDeltaGamma()
{
    int nVar = vars->gammaMat.M();
    int m = vars->gammaMat.N();

    if( m == 1 )
        return;

    vars->deltaXi.Submatrix( vars->deltaMat, nVar, 1, 0, stats->itCount % m );
    vars->gamma.Submatrix( vars->gammaMat, nVar, 1, 0, stats->itCount % m );
}

} // namespace blockSQP
