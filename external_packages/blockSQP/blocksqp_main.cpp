/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_main.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of SQPmethod class.
 *
 */


#include "blocksqp.hpp"

namespace blockSQP
{

SQPmethod::SQPmethod( Problemspec *problem, SQPoptions *parameters, SQPstats *statistics )
{
    prob = problem;
    param = parameters;
    stats = statistics;

    // Check if there are options that are infeasible and set defaults accordingly
    param->optionsConsistency();

    vars = new SQPiterate( prob, param, 1 );

    if( param->sparseQP < 2 )
    {
        qp = new qpOASES::SQProblem( prob->nVar, prob->nCon );
        qpSave = new qpOASES::SQProblem( prob->nVar, prob->nCon );
    }
    else
    {
        qp = new qpOASES::SQProblemSchur( prob->nVar, prob->nCon, qpOASES::HST_UNKNOWN, 50 );
        qpSave = new qpOASES::SQProblemSchur( prob->nVar, prob->nCon, qpOASES::HST_UNKNOWN, 50 );
    }

    initCalled = false;
}

SQPmethod::~SQPmethod()
{
    delete qp;
    delete qpSave;
    delete vars;
}


void SQPmethod::init()
{
    // Print header and information about the algorithmic parameters
    printInfo( param->printLevel );

    // Open output files
    stats->initStats( param );
    vars->initIterate( param );

    // Initialize filter with pair ( maxConstrViolation, objLowerBound )
    initializeFilter();

    // Set initial values for all xi and set the Jacobian for linear constraints
    if( param->sparseQP )
        prob->initialize( vars->xi, vars->lambda, vars->jacNz, vars->jacIndRow, vars->jacIndCol );
    else
        prob->initialize( vars->xi, vars->lambda, vars->constrJac );

    initCalled = true;
}


int SQPmethod::run( int maxIt, int warmStart )
{
    int it, infoQP = 0, infoEval = 0;
    bool skipLineSearch = false;
    bool hasConverged = false;
    int whichDerv = param->whichSecondDerv;

    if( !initCalled )
    {
        printf("init() must be called before run(). Aborting.\n");
        return -1;
    }

    if( warmStart == 0 || stats->itCount == 0 )
    {
        // SQP iteration 0

        /// Set initial Hessian approximation
        calcInitialHessian();

        /// Evaluate all functions and gradients for xi_0
        if( param->sparseQP )
            prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                            vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->hess, 1+whichDerv, &infoEval );
        else
            prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                            vars->constrJac, vars->hess, 1+whichDerv, &infoEval );
        stats->nDerCalls++;

        /// Check if converged
        hasConverged = calcOptTol();
        stats->printProgress( prob, vars, param, hasConverged );
        if( hasConverged )
            return 0;

        stats->itCount++;
    }

    /*
     * SQP Loop: during first iteration, stats->itCount = 1
     */
    for( it=0; it<maxIt; it++ )
    {
        /// Solve QP subproblem with qpOASES or QPOPT
        updateStepBounds( 0 );
        infoQP = solveQP( vars->deltaXi, vars->lambdaQP );

        if( infoQP == 1 )
        {// 1.) Maximum number of iterations reached
            printf("***Warning! Maximum number of QP iterations exceeded.***\n");
            ;// just continue ...
        }
        else if( infoQP == 2 || infoQP > 3 )
        {// 2.) QP error (e.g., unbounded), solve again with pos.def. diagonal matrix (identity)
            printf("***QP error. Solve again with identity matrix.***\n");
            resetHessian();
            infoQP = solveQP( vars->deltaXi, vars->lambdaQP );
            if( infoQP )
            {// If there is still an error, terminate.
                printf( "***QP error. Stop.***\n" );
                return -1;
            }
            else
                vars->steptype = 1;
        }
        else if( infoQP == 3 )
        {// 3.) QP infeasible, try to restore feasibility
            bool qpError = true;
            skipLineSearch = true; // don't do line search with restoration step

            // Try to reduce constraint violation by heuristic
            if( vars->steptype < 2 )
            {
                printf("***QP infeasible. Trying to reduce constraint violation...");
                qpError = feasibilityRestorationHeuristic();
                if( !qpError )
                {
                    vars->steptype = 2;
                    printf("Success.***\n");
                }
                else
                    printf("Failed.***\n");
            }

            // Invoke feasibility restoration phase
            //if( qpError && vars->steptype < 3 && param->restoreFeas )
            if( qpError && param->restoreFeas && vars->cNorm > 0.01 * param->nlinfeastol )
            {
                printf("***Start feasibility restoration phase.***\n");
                vars->steptype = 3;
                qpError = feasibilityRestorationPhase();
            }

            // If everything failed, abort.
            if( qpError )
            {
                printf( "***QP error. Stop.***\n" );
                return -1;
            }
        }

        /// Determine steplength alpha
        if( param->globalization == 0 || (param->skipFirstGlobalization && stats->itCount == 1) )
        {// No globalization strategy, but reduce step if function cannot be evaluated
            if( fullstep() )
            {
                printf( "***Constraint or objective could not be evaluated at new point. Stop.***\n" );
                return -1;
            }
            vars->steptype = 0;
        }
        else if( param->globalization == 1 && !skipLineSearch )
        {// Filter line search based on Waechter et al., 2006 (Ipopt paper)
            if( filterLineSearch() || vars->reducedStepCount > param->maxConsecReducedSteps )
            {
                // Filter line search did not produce a step. Now there are a few things we can try ...
                bool lsError = true;

                // Heuristic 1: Check if the full step reduces the KKT error by at least kappaF, if so, accept the step.
                lsError = kktErrorReduction( );
                if( !lsError )
                    vars->steptype = -1;

                // Heuristic 2: Try to reduce constraint violation by closing continuity gaps to produce an admissable iterate
                if( lsError && vars->cNorm > 0.01 * param->nlinfeastol && vars->steptype < 2 )
                {// Don't do this twice in a row!

                    printf("***Warning! Steplength too short. Trying to reduce constraint violation...");

                    // Integration over whole time interval
                    lsError = feasibilityRestorationHeuristic( );
                    if( !lsError )
                    {
                        vars->steptype = 2;
                        printf("Success.***\n");
                    }
                    else
                        printf("Failed.***\n");
                }

                // Heuristic 3: Recompute step with a diagonal Hessian
                if( lsError && vars->steptype != 1 && vars->steptype != 2 )
                {// After closing continuity gaps, we already take a step with initial Hessian. If this step is not accepted then this will cause an infinite loop!

                    printf("***Warning! Steplength too short. Trying to find a new step with identity Hessian.***\n");
                    vars->steptype = 1;

                    resetHessian();
                    continue;
                }

                // If this does not yield a successful step, start restoration phase
                if( lsError && vars->cNorm > 0.01 * param->nlinfeastol && param->restoreFeas )
                {
                    printf("***Warning! Steplength too short. Start feasibility restoration phase.***\n");
                    vars->steptype = 3;

                    // Solve NLP with minimum norm objective
                    lsError = feasibilityRestorationPhase( );
                }

                // If everything failed, abort.
                if( lsError )
                {
                    printf( "***Line search error. Stop.***\n" );
                    return -1;
                }
            }
            else
                vars->steptype = 0;
        }

        /// Calculate "old" Lagrange gradient: gamma = dL(xi_k, lambda_k+1)
        calcLagrangeGradient( vars->gamma, 0 );

        /// Evaluate functions and gradients at the new xi
        if( param->sparseQP )
            prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                            vars->jacNz, vars->jacIndRow, vars->jacIndCol, vars->hess, 1+whichDerv, &infoEval );
        else
            prob->evaluate( vars->xi, vars->lambda, &vars->obj, vars->constr, vars->gradObj,
                            vars->constrJac, vars->hess, 1+whichDerv, &infoEval );
        stats->nDerCalls++;

        /// Check if converged
        hasConverged = calcOptTol();

        /// Print one line of output for the current iteration
        stats->printProgress( prob, vars, param, hasConverged );
        if( hasConverged && vars->steptype < 2 )
        {
            stats->itCount++;
            if( param->debugLevel > 2 )
            {
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
        if( param->hessUpdate < 4 && !param->hessLimMem )
            calcHessianUpdate( param->hessUpdate, param->hessScaling );
        else if( param->hessUpdate < 4 && param->hessLimMem )
            calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling );
        else if( param->hessUpdate == 4 )
            calcFiniteDiffHessian( );

        // If limited memory updates  are used, set pointers deltaXi and gamma to the next column in deltaMat and gammaMat
        updateDeltaGamma();

        stats->itCount++;
        skipLineSearch = false;
    }

    return 1;
}


void SQPmethod::finish()
{
    if( initCalled )
        initCalled = false;
    else
    {
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
void SQPmethod::calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, double *jacNz, int *jacIndRow, int *jacIndCol,
                                      Matrix &gradLagrange, int flag )
{
    int iVar, iCon;

    // Objective gradient
    if( flag == 0 )
        for( iVar=0; iVar<prob->nVar; iVar++ )
            gradLagrange( iVar ) = gradObj( iVar );
    else if( flag == 1 )
        for( iVar=0; iVar<prob->nVar; iVar++ )
            gradLagrange( iVar ) = gradObj( iVar ) - gradLagrange( iVar );
    else
        gradLagrange.Initialize( 0.0 );

    // - lambdaT * constrJac
    for( iVar=0; iVar<prob->nVar; iVar++ )
        for( iCon=jacIndCol[iVar]; iCon<jacIndCol[iVar+1]; iCon++ )
            gradLagrange( iVar ) -= lambda( prob->nVar + jacIndRow[iCon] ) * jacNz[iCon];

    // - lambdaT * simpleBounds
    for( iVar=0; iVar<prob->nVar; iVar++ )
        gradLagrange( iVar ) -= lambda( iVar );
}


/**
 * Compute gradient of Lagrangian or difference of Lagrangian gradients (dense version)
 *
 * flag == 0: output dL(xi, lambda)
 * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
 * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
 */
void SQPmethod::calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, const Matrix &constrJac,
                                      Matrix &gradLagrange, int flag )
{
    int iVar, iCon;

    // Objective gradient
    if( flag == 0 )
        for( iVar=0; iVar<prob->nVar; iVar++ )
            gradLagrange( iVar ) = gradObj( iVar );
    else if( flag == 1 )
        for( iVar=0; iVar<prob->nVar; iVar++ )
            gradLagrange( iVar ) = gradObj( iVar ) - gradLagrange( iVar );
    else
        gradLagrange.Initialize( 0.0 );

    // - lambdaT * constrJac
    for( iVar=0; iVar<prob->nVar; iVar++ )
        for( iCon=0; iCon<prob->nCon; iCon++ )
            gradLagrange( iVar ) -= lambda( prob->nVar + iCon ) * constrJac( iCon, iVar );

    // - lambdaT * simpleBounds
    for( iVar=0; iVar<prob->nVar; iVar++ )
        gradLagrange( iVar ) -= lambda( iVar );
}


/**
 * Wrapper if called with standard arguments
 */
void SQPmethod::calcLagrangeGradient( Matrix &gradLagrange, int flag )
{
    if( param->sparseQP )
        calcLagrangeGradient( vars->lambda, vars->gradObj, vars->jacNz, vars->jacIndRow, vars->jacIndCol, gradLagrange, flag );
    else
        calcLagrangeGradient( vars->lambda, vars->gradObj, vars->constrJac, gradLagrange, flag );
}


/**
 * Compute optimality conditions:
 * ||gradLagrange(xi,lambda)||_infty / (1 + ||lambda||_infty) <= TOL
 * and
 * ||constrViolation||_infty / (1 + ||xi||_infty) <= TOL
 */
bool SQPmethod::calcOptTol()
{
    // scaled norm of Lagrangian gradient
    calcLagrangeGradient( vars->gradLagrange, 0 );
    vars->gradNorm = lInfVectorNorm( vars->gradLagrange );
    vars->tol = vars->gradNorm /( 1.0 + lInfVectorNorm( vars->lambda ) );

    // norm of constraint violation
    vars->cNorm  = lInfConstraintNorm( vars->xi, vars->constr, prob->bu, prob->bl );
    vars->cNormS = vars->cNorm /( 1.0 + lInfVectorNorm( vars->xi ) );

    if( vars->tol <= param->opttol && vars->cNormS <= param->nlinfeastol )
        return true;
    else
        return false;
}

void SQPmethod::printInfo( int printLevel )
{
    char hessString1[100];
    char hessString2[100];
    char globString[100];
    char qpString[100];

    if( printLevel == 0 )
        return;

    /* QP Solver */
    if( param->sparseQP == 0 )
        strcpy( qpString, "dense, reduced Hessian factorization" );
    else if( param->sparseQP == 1 )
        strcpy( qpString, "sparse, reduced Hessian factorization" );
    else if( param->sparseQP == 2 )
        strcpy( qpString, "sparse, Schur complement approach" );

    /* Globalization */
    if( param->globalization == 0 )
        strcpy( globString, "none (full step)" );
    else if( param->globalization == 1 )
        strcpy( globString, "filter line search" );

    /* Hessian approximation */
    if( param->blockHess && (param->hessUpdate == 1 || param->hessUpdate == 2) )
        strcpy( hessString1, "block " );
    else
        strcpy( hessString1, "" );

    if( param->hessLimMem && (param->hessUpdate == 1 || param->hessUpdate == 2) )
        strcat( hessString1, "L-" );

    /* Fallback Hessian */
    if( param->hessUpdate == 1 || param->hessUpdate == 4 || (param->hessUpdate == 2 && !param->hessDamp) )
    {
        strcpy( hessString2, hessString1 );

        /* Fallback Hessian update type */
        if( param->fallbackUpdate == 0 )
            strcat( hessString2, "Id" );
        else if( param->fallbackUpdate == 1 )
            strcat( hessString2, "SR1" );
        else if( param->fallbackUpdate == 2 )
            strcat( hessString2, "BFGS" );
        else if( param->fallbackUpdate == 4 )
            strcat( hessString2, "Finite differences" );

        /* Fallback Hessian scaling */
        if( param->fallbackScaling == 1 )
            strcat( hessString2, ", SP" );
        else if( param->fallbackScaling == 2 )
            strcat( hessString2, ", OL" );
        else if( param->fallbackScaling == 3 )
            strcat( hessString2, ", mean" );
        else if( param->fallbackScaling == 4 )
            strcat( hessString2, ", selective sizing" );
    }
    else
        strcpy( hessString2, "-" );

    /* First Hessian update type */
    if( param->hessUpdate == 0 )
        strcat( hessString1, "Id" );
    else if( param->hessUpdate == 1 )
        strcat( hessString1, "SR1" );
    else if( param->hessUpdate == 2 )
        strcat( hessString1, "BFGS" );
    else if( param->hessUpdate == 4 )
        strcat( hessString1, "Finite differences" );

    /* First Hessian scaling */
    if( param->hessScaling == 1 )
        strcat( hessString1, ", SP" );
    else if( param->hessScaling == 2 )
        strcat( hessString1, ", OL" );
    else if( param->hessScaling == 3 )
        strcat( hessString1, ", mean" );
    else if( param->hessScaling == 4 )
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

} // namespace blockSQP
