/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_glob.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of methods of SQPmethod class associated with the
 *  globalization strategy.
 *
 */

#include "blocksqp_iterate.hpp"
#include "blocksqp_options.hpp"
#include "blocksqp_stats.hpp"
#include "blocksqp_method.hpp"
#include "blocksqp_restoration.hpp"
#include "blocksqp_general_purpose.hpp"

namespace blockSQP
{

void SQPmethod::acceptStep( double alpha )
{
    acceptStep( vars->deltaXi, vars->lambdaQP, alpha, 0 );
}

void SQPmethod::acceptStep( const Matrix &deltaXi, const Matrix &lambdaQP, double alpha, int nSOCS )
{
    int k;
    double lStpNorm;

    // Current alpha
    vars->alpha = alpha;
    vars->nSOCS = nSOCS;

    // Set new xi by accepting the current trial step
    for( k=0; k<vars->xi.M(); k++ )
    {
        vars->xi( k ) = vars->trialXi( k );
        vars->deltaXi( k ) = alpha * deltaXi( k );
    }

    // Store the infinity norm of the multiplier step
    vars->lambdaStepNorm = 0.0;
    for( k=0; k<vars->lambda.M(); k++ )
        if( (lStpNorm = fabs( alpha*lambdaQP( k ) - alpha*vars->lambda( k ) )) > vars->lambdaStepNorm )
            vars->lambdaStepNorm = lStpNorm;

    // Set new multipliers
    for( k=0; k<vars->lambda.M(); k++ )
        vars->lambda( k ) = (1.0 - alpha)*vars->lambda( k ) + alpha*lambdaQP( k );

    // Count consecutive reduced steps
    if( vars->alpha < 1.0 )
        vars->reducedStepCount++;
    else
        vars->reducedStepCount = 0;
}

void SQPmethod::reduceStepsize( double *alpha )
{
    *alpha = (*alpha) * 0.5;
}

void SQPmethod::reduceSOCStepsize( double *alphaSOC )
{
    int i;
    int nVar = prob->nVar;

    // Update bounds on linearized constraints for the next SOC QP:
    // That is different from the update for the first SOC QP!
    for( i=0; i<prob->nCon; i++ )
    {
        if( prob->bl( nVar+i ) != param->inf )
            vars->deltaBl( nVar+i ) = (*alphaSOC)*vars->deltaBl( nVar+i ) - vars->constr( i );
        else
            vars->deltaBl( nVar+i ) = param->inf;

        if( prob->bu( nVar+i ) != param->inf )
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
int SQPmethod::fullstep()
{
    double alpha;
    double objTrial, cNormTrial;
    int i, k, info;
    int nVar = prob->nVar;

    // Backtracking line search
    alpha = 1.0;
    for( k=0; k<10; k++ )
    {
        // Compute new trial point
        for( i=0; i<nVar; i++ )
            vars->trialXi( i ) = vars->xi( i ) + alpha * vars->deltaXi( i );

        // Compute problem functions at trial point
        prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
        stats->nFunCalls++;
        cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
        // Reduce step if evaluation fails, if lower bound is violated or if objective or a constraint is NaN
        if( info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
        {
            printf("info=%i, objTrial=%g\n", info, objTrial );
            // evaluation error, reduce stepsize
            reduceStepsize( &alpha );
            continue;
        }
        else
        {
            acceptStep( alpha );
            return 0;
        }
    }// backtracking steps

    return 1;
}


/**
 *
 * Backtracking line search based on a filter
 * as described in Ipopt paper (Waechter 2006)
 *
 */
int SQPmethod::filterLineSearch()
{
    double alpha = 1.0;
    double cNorm, cNormTrial, objTrial, dfTdeltaXi;

    int i, k, info;
    int nVar = prob->nVar;

    // Compute ||constr(xi)|| at old point
    cNorm = lInfConstraintNorm( vars->xi, vars->constr, prob->bu, prob->bl );

    // Backtracking line search
    for( k=0; k<param->maxLineSearch; k++ )
    {
        // Compute new trial point
        for( i=0; i<nVar; i++ )
            vars->trialXi( i ) = vars->xi( i ) + alpha * vars->deltaXi( i );

        // Compute grad(f)^T * deltaXi
        dfTdeltaXi = 0.0;
        for( i=0; i<nVar; i++ )
            dfTdeltaXi += vars->gradObj( i ) * vars->deltaXi( i );

        // Compute objective and at ||constr(trialXi)||_1 at trial point
        prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
        stats->nFunCalls++;
        cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
        // Reduce step if evaluation fails, if lower bound is violated or if objective is NaN
        if( info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
        {
            // evaluation error, reduce stepsize
            reduceStepsize( &alpha );
            continue;
        }

        // Check acceptability to the filter
        if( pairInFilter( cNormTrial, objTrial ) )
        {
            // Trial point is in the prohibited region defined by the filter, try second order correction
            if( secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 0, k ) )
                break; // SOC yielded suitable alpha, stop
            else
            {
                reduceStepsize( &alpha );
                continue;
            }
        }

        // Check sufficient decrease, case I:
        // If we are (almost) feasible and a "switching condition" is satisfied
        // require sufficient progress in the objective instead of bi-objective condition
        if( cNorm <= param->thetaMin )
        {
            // Switching condition, part 1: grad(f)^T * deltaXi < 0 ?
            if( dfTdeltaXi < 0 )
                // Switching condition, part 2: alpha * ( - grad(f)^T * deltaXi )**sF > delta * cNorm**sTheta ?
                if( alpha * pow( (-dfTdeltaXi), param->sF ) > param->delta * pow( cNorm, param->sTheta ) )
                {
                    // Switching conditions hold: Require satisfaction of Armijo condition for objective
                    if( objTrial > vars->obj + param->eta*alpha*dfTdeltaXi )
                    {
                        // Armijo condition violated, try second order correction
                        if( secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 1, k ) )
                            break; // SOC yielded suitable alpha, stop
                        else
                        {
                            reduceStepsize( &alpha );
                            continue;
                        }
                    }
                    else
                    {
                        // found suitable alpha, stop
                        acceptStep( alpha );
                        break;
                    }
                }
        }

        // Check sufficient decrease, case II:
        // Bi-objective (filter) condition
        if( cNormTrial < (1.0 - param->gammaTheta) * cNorm || objTrial < vars->obj - param->gammaF * cNorm )
        {
            // found suitable alpha, stop
            acceptStep( alpha );
            break;
        }
        else
        {
            // Trial point is dominated by current point, try second order correction
            if( secondOrderCorrection( cNorm, cNormTrial, dfTdeltaXi, 0, k ) )
                break; // SOC yielded suitable alpha, stop
            else
            {
                reduceStepsize( &alpha );
                continue;
            }
        }
    }// backtracking steps

    // No step could be found by the line search
    if( k == param->maxLineSearch )
        return 1;

    // Augment the filter if switching condition or Armijo condition does not hold
    if( dfTdeltaXi >= 0 )
        augmentFilter( cNormTrial, objTrial );
    else if( alpha * pow( (-dfTdeltaXi), param->sF ) > param->delta * pow( cNorm, param->sTheta ) )// careful with neg. exponents!
        augmentFilter( cNormTrial, objTrial );
    else if( objTrial <= vars->obj + param->eta*alpha*dfTdeltaXi )
        augmentFilter( cNormTrial, objTrial );


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
bool SQPmethod::secondOrderCorrection( double cNorm, double cNormTrial, double dfTdeltaXi, bool swCond, int it )
{
    // Perform SOC only on the first iteration of backtracking line search
    if( it > 0 )
        return false;
    // If constraint violation of the trialstep is lower than the current one skip SOC
    if( cNormTrial < cNorm )
        return false;

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
    for( k=0; k<param->maxSOCiter; k++ )
    {
        nSOCS++;

        // Update bounds for SOC QP
        updateStepBounds( 1 );

        // Solve SOC QP to obtain new, corrected deltaXi
        // (store in separate vector to avoid conflict with original deltaXi -> need it in linesearch!)
        info = solveQP( deltaXiSOC, lambdaQPSOC, false );
        if( info != 0 )
            return false; // Could not solve QP, abort SOC

        // Set new SOC trial point
        for( i=0; i<nVar; i++ )
            vars->trialXi( i ) = vars->xi( i ) + deltaXiSOC( i );

        // Compute objective and ||constr(trialXiSOC)||_1 at SOC trial point
        prob->evaluate( vars->trialXi, &objTrialSOC, vars->constr, &info );
        stats->nFunCalls++;
        cNormTrialSOC = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
        if( info != 0 || objTrialSOC < prob->objLo || objTrialSOC > prob->objUp || !(objTrialSOC == objTrialSOC) || !(cNormTrialSOC == cNormTrialSOC) )
            return false; // evaluation error, abort SOC

        // Check acceptability to the filter (in SOC)
        if( pairInFilter( cNormTrialSOC, objTrialSOC ) )
            return false; // Trial point is in the prohibited region defined by the filter, abort SOC

        // Check sufficient decrease, case I (in SOC)
        // (Almost feasible and switching condition holds for line search alpha)
        if( cNorm <= param->thetaMin && swCond )
        {
            if( objTrialSOC > vars->obj + param->eta*dfTdeltaXi )
            {
                // Armijo condition does not hold for SOC step, next SOC step

                // If constraint violation gets worse by SOC, abort
                if( cNormTrialSOC > param->kappaSOC * cNormOld )
                    return false;
                else
                    cNormOld = cNormTrialSOC;
                continue;
            }
            else
            {
                // found suitable alpha during SOC, stop
                acceptStep( deltaXiSOC, lambdaQPSOC, 1.0, nSOCS );
                return true;
            }
        }

        // Check sufficient decrease, case II (in SOC)
        if( cNorm > param->thetaMin || !swCond )
        {
            if( cNormTrialSOC < (1.0 - param->gammaTheta) * cNorm || objTrialSOC < vars->obj - param->gammaF * cNorm )
            {
                // found suitable alpha during SOC, stop
                acceptStep( deltaXiSOC, lambdaQPSOC, 1.0, nSOCS );
                return true;
            }
            else
            {
                // Trial point is dominated by current point, next SOC step

                // If constraint violation gets worse by SOC, abort
                if( cNormTrialSOC > param->kappaSOC * cNormOld )
                    return false;
                else
                    cNormOld = cNormTrialSOC;
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
int SQPmethod::feasibilityRestorationPhase()
{
    // No Feasibility restoration phase
    if( param->restoreFeas == 0 )
        return -1;

    stats->nRestPhaseCalls++;

    int ret, it, i, k, info;
    int maxRestIt = 100;
    int warmStart;
    double cNormTrial, objTrial, lStpNorm;
    RestorationProblem *restProb;
    SQPmethod *restMethod;
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
    restMethod = new SQPmethod( restProb, restOpts, restStats );
    restMethod->init();

    // Iterate until a point acceptable to the filter is found
    warmStart = 0;
    for( it=0; it<maxRestIt; it++ )
    {
        // One iteration for minimum norm NLP
        ret = restMethod->run( 1, warmStart );
        warmStart = 1;

        // If restMethod yields error, stop restoration phase
        if( ret == -1 )
            break;

        // Get new xi from the restoration phase
        for( i=0; i<prob->nVar; i++ )
            vars->trialXi( i ) = restMethod->vars->xi( i );

        // Compute objective at trial point
        prob->evaluate( vars->trialXi, &objTrial, vars->constr, &info );
        stats->nFunCalls++;
        cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
        if( info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
            continue;

        // Is this iterate acceptable for the filter?
        if( !pairInFilter( cNormTrial, objTrial ) )
        {
            // success
            printf("Found a point acceptable for the filter.\n");
            ret = 0;
            break;
        }

        // If minimum norm NLP has converged, declare local infeasibility
        if( restMethod->vars->tol < param->opttol && restMethod->vars->cNormS < param->nlinfeastol )
        {
            ret = 1;
            break;
        }
    }

    // Success or locally infeasible
    if( ret == 0 || ret == 1 )
    {
        // Store the infinity norm of the multiplier step
        vars->lambdaStepNorm = 0.0;
        // Compute restoration step
        for( k=0; k<prob->nVar; k++ )
        {
            vars->deltaXi( k ) = vars->xi( k );

            vars->xi( k ) = vars->trialXi( k );

            // Store lInf norm of dual step
            if( (lStpNorm = fabs( restMethod->vars->lambda( k ) - vars->lambda( k ) )) > vars->lambdaStepNorm )
                vars->lambdaStepNorm = lStpNorm;
            vars->lambda( k ) = restMethod->vars->lambda( k );
            vars->lambdaQP( k ) = restMethod->vars->lambdaQP( k );

            vars->deltaXi( k ) -= vars->xi( k );
        }
        for( k=prob->nVar; k<prob->nVar+prob->nCon; k++ )
        {
            // skip the dual variables for the slack variables in the restoration problem
            if( (lStpNorm = fabs( restMethod->vars->lambda( 2*prob->nCon + k ) - vars->lambda( k ) )) > vars->lambdaStepNorm )
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

    if( ret == 1 )
    {
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
int SQPmethod::feasibilityRestorationHeuristic()
{
    stats->nRestHeurCalls++;

    int info, k;
    double cNormTrial;

    info = 0;

    // Call problem specific heuristic to reduce constraint violation.
    // For shooting methods that means setting consistent values for shooting nodes by one forward integration.
    for( k=0; k<prob->nVar; k++ ) // input: last successful step
        vars->trialXi( k ) = vars->xi( k );
    prob->reduceConstrVio( vars->trialXi, &info );
    if( info )// If an error occured in restoration heuristics, abort
        return -1;

    // Compute objective and constraints at the new (hopefully feasible) point
    prob->evaluate( vars->trialXi, &vars->obj, vars->constr, &info );
    stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm( vars->trialXi, vars->constr, prob->bu, prob->bl );
    if( info != 0 || vars->obj < prob->objLo || vars->obj > prob->objUp || !(vars->obj == vars->obj) || !(cNormTrial == cNormTrial) )
        return -1;

    // Is the new point acceptable for the filter?
    if( pairInFilter( cNormTrial, vars->obj ) )
    {
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
    for( k=0; k<prob->nVar; k++ )
    {
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
int SQPmethod::kktErrorReduction( )
{
    int i, info = 0;
    double objTrial, cNormTrial, trialGradNorm, trialTol;
    Matrix trialConstr, trialGradLagrange;

    // Compute new trial point
    for( i=0; i<prob->nVar; i++ )
        vars->trialXi( i ) = vars->xi( i ) + vars->deltaXi( i );

    // Compute objective and ||constr(trialXi)|| at trial point
    trialConstr.Dimension( prob->nCon ).Initialize( 0.0 );
    prob->evaluate( vars->trialXi, &objTrial, trialConstr, &info );
    stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm( vars->trialXi, trialConstr, prob->bu, prob->bl );
    if( info != 0 || objTrial < prob->objLo || objTrial > prob->objUp || !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
    {
        // evaluation error
        return 1;
    }

    // Compute KKT error of the new point

    // scaled norm of Lagrangian gradient
    trialGradLagrange.Dimension( prob->nVar ).Initialize( 0.0 );
    if( param->sparseQP )
        calcLagrangeGradient( vars->lambdaQP, vars->gradObj, vars->jacNz,
                              vars->jacIndRow, vars->jacIndCol, trialGradLagrange, 0 );
    else
        calcLagrangeGradient( vars->lambdaQP, vars->gradObj, vars->constrJac,
                              trialGradLagrange, 0 );

    trialGradNorm = lInfVectorNorm( trialGradLagrange );
    trialTol = trialGradNorm /( 1.0 + lInfVectorNorm( vars->lambdaQP ) );

    if( fmax( cNormTrial, trialTol ) < param->kappaF * fmax( vars->cNorm, vars->tol ) )
    {
        acceptStep( 1.0 );
        return 0;
    }
    else
        return 1;
}

/**
 * Check if current entry is accepted to the filter:
 * (cNorm, obj) in F_k
 */
bool SQPmethod::pairInFilter( double cNorm, double obj )
{
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

    for( iter=filter->begin(); iter!=filter->end(); iter++ )
        if( (cNorm >= (1.0 - param->gammaTheta) * iter->first ||
            (cNorm < 0.01 * param->nlinfeastol && iter->first < 0.01 * param->nlinfeastol ) ) &&
            obj >= iter->second - param->gammaF * iter->first )
        {
            return 1;
        }

    return 0;
}


void SQPmethod::initializeFilter()
{
    std::set< std::pair<double,double> >::iterator iter;
    std::pair<double,double> initPair ( param->thetaMax, prob->objLo );

    // Remove all elements
    iter=vars->filter->begin();
    while (iter != vars->filter->end())
    {
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
void SQPmethod::augmentFilter( double cNorm, double obj )
{
    std::set< std::pair<double,double> >::iterator iter;
    std::pair<double,double> entry ( (1.0 - param->gammaTheta)*cNorm, obj - param->gammaF*cNorm );

    // Augment filter by current element
    vars->filter->insert( entry );

    // Remove dominated elements
    iter=vars->filter->begin();
    while (iter != vars->filter->end())
    {
        //printf(" iter->first=%g, entry.first=%g, iter->second=%g, entry.second=%g\n",iter->first, entry.first, iter->second, entry.second);
        if( iter->first > entry.first && iter->second > entry.second )
        {
            std::set< std::pair<double,double> >::iterator iterToRemove = iter;
            iter++;
            vars->filter->erase( iterToRemove );
        }
        else
            iter++;
    }
    // Print current filter
    //for( iter=vars->filter->begin(); iter!=vars->filter->end(); iter++ )
        //printf("(%g,%g), ", iter->first, iter->second);
    //printf("\n");
}

} // namespace blockSQP
