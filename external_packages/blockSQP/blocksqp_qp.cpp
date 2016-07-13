/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_qp.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of methods of SQPmethod class associated with
 *  solution of the quadratic subproblems.
 *
 */

#include "blocksqp_iterate.hpp"
#include "blocksqp_options.hpp"
#include "blocksqp_stats.hpp"
#include "blocksqp_method.hpp"
#include "blocksqp_general_purpose.hpp"

namespace blockSQP
{

void SQPmethod::computeNextHessian( int idx, int maxQP )
{
    // Compute fallback update only once
    if( idx == 1 )
    {
        // Switch storage
        vars->hess = vars->hess2;

        // If last block contains exact Hessian, we need to copy it
        if( param->whichSecondDerv == 1 )
            for( int i=0; i<vars->hess[prob->nBlocks-1].M(); i++ )
                for( int j=i; j<vars->hess[prob->nBlocks-1].N(); j++ )
                    vars->hess2[prob->nBlocks-1]( i,j ) = vars->hess1[prob->nBlocks-1]( i,j );

        // Limited memory: compute fallback update only when needed
        if( param->hessLimMem )
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
    if( maxQP > 2 )
    {
        /* Convexification parameter: mu_l = l / (maxQP-1).
         * Compute it only in the first iteration, afterwards update
         * by recursion: mu_l/mu_(l-1) */
        double idxF = (double) idx;
        double mu = (idx==1) ? 1.0 / (maxQP-1) : idxF / (idxF - 1.0);
        double mu1 = 1.0 - mu;
        for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
            for( int i=0; i<vars->hess[iBlock].M(); i++ )
                for( int j=i; j<vars->hess[iBlock].N(); j++ )
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
int SQPmethod::solveQP( Matrix &deltaXi, Matrix &lambdaQP, bool matricesChanged )
{
    Matrix jacT;
    int maxQP, l;
    if( param->globalization == 1 &&
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
    if( matricesChanged )
    {
        if( param->sparseQP )
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
    if( matricesChanged && maxQP > 1 )
        opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    else
        opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    qp->setOptions( opts );

    if( maxQP > 1 )
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
    for( l=0; l<maxQP; l++ )
    {
        /*
         * Compute a new Hessian
         */
        if( l > 0 )
        {// If the solution of the first QP was rejected, consider second Hessian
            stats->qpResolve++;
            *qp = *qpSave;

            computeNextHessian( l, maxQP );
        }

        if( l == maxQP-1 )
        {// Enable inertia correction for supposedly convex QPs, just in case
            opts.enableInertiaCorrection = qpOASES::BT_TRUE;
            qp->setOptions( opts );
        }

        /*
         * Prepare the current Hessian for qpOASES
         */
        if( matricesChanged )
        {
            if( param->sparseQP )
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
        if( param->debugLevel > 2 ) stats->dumpQPCpp( prob, vars, qp, param->sparseQP );
        if( matricesChanged )
        {
            maxIt = param->maxItQP;
            cpuTime = param->maxTimeQP;
            if( qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                qp->getStatus() == qpOASES::QPS_SOLVED )
            {
                ret = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
            }
            else
            {
                ret = qp->init( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
            }
        }
        else if( !matricesChanged ) // Second order correction: H and A do not change
        {
            maxIt = 0.1*param->maxItQP;
            cpuTime = 0.1*param->maxTimeQP;
            ret = qp->hotstart( g, lb, lu, lbA, luA, maxIt, &cpuTime );
        }

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if( l < maxQP-1 && matricesChanged )
        {
            if( ret == qpOASES::SUCCESSFUL_RETURN )
            {
                if( param->sparseQP == 2 )
                    ret = solAna.checkCurvatureOnStronglyActiveConstraints( dynamic_cast<qpOASES::SQProblemSchur*>(qp) );
                else
                    ret = solAna.checkCurvatureOnStronglyActiveConstraints( qp );
            }

            if( ret == qpOASES::SUCCESSFUL_RETURN )
            {// QP was solved successfully and curvature is positive after removing bounds
                stats->qpIterations = maxIt + 1;
                break; // Success!
            }
            else
            {// QP solution is rejected, save statistics
                if( ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
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
    if( param->sparseQP )
        Atimesb( vars->jacNz, vars->jacIndRow, vars->jacIndCol, deltaXi, vars->AdeltaXi );
    else
        Atimesb( vars->constrJac, deltaXi, vars->AdeltaXi );

    // Print qpOASES error code, if any
    if( ret != qpOASES::SUCCESSFUL_RETURN && matricesChanged )
        printf( "qpOASES error message: \"%s\"\n",
                qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );

    // Point Hessian again to the first Hessian
    vars->hess = vars->hess1;

    /* For full-memory Hessian: Restore fallback Hessian if convex combinations
     * were used during the loop */
    if( !param->hessLimMem && maxQP > 2 && matricesChanged )
    {
        double mu = 1.0 / ((double) l);
        double mu1 = 1.0 - mu;
        int nBlocks = (param->whichSecondDerv == 1) ? vars->nBlocks-1 : vars->nBlocks;
        for( int iBlock=0; iBlock<nBlocks; iBlock++ )
            for( int i=0; i<vars->hess[iBlock].M(); i++ )
                for( int j=i; j<vars->hess[iBlock].N(); j++ )
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
    if( ret == qpOASES::SUCCESSFUL_RETURN )
        return 0;
    else if( ret == qpOASES::RET_MAX_NWSR_REACHED )
        return 1;
    else if( ret == qpOASES::RET_HESSIAN_NOT_SPD ||
             ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS ||
             ret == qpOASES::RET_QP_UNBOUNDED ||
             ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS )
        return 2;
    else if( ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY ||
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
void SQPmethod::updateStepBounds( bool soc )
{
    int i;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // Bounds on step
    for( i=0; i<nVar; i++ )
    {
        if( prob->bl(i) != param->inf )
            vars->deltaBl( i ) = prob->bl( i ) - vars->xi( i );
        else
            vars->deltaBl( i ) = param->inf;

        if( prob->bu(i) != param->inf )
            vars->deltaBu( i ) = prob->bu( i ) - vars->xi( i );
        else
            vars->deltaBu( i ) = param->inf;
    }

    // Bounds on linearized constraints
    for( i=0; i<nCon; i++ )
    {
        if( prob->bl( nVar+i ) != param->inf )
        {
            vars->deltaBl( nVar+i ) = prob->bl( nVar+i ) - vars->constr( i );
            if( soc ) vars->deltaBl( nVar+i ) += vars->AdeltaXi( i );
        }
        else
            vars->deltaBl( nVar+i ) = param->inf;

        if( prob->bu( nVar+i ) != param->inf )
        {
            vars->deltaBu( nVar+i ) = prob->bu( nVar+i ) - vars->constr( i );
            if( soc ) vars->deltaBu( nVar+i ) += vars->AdeltaXi( i );
        }
        else
            vars->deltaBu( nVar+i ) = param->inf;
    }
}

} // namespace blockSQP


