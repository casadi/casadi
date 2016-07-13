/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_method.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of blockSQP's central SQPmethod class.
 */

#ifndef BLOCKSQP_METHOD_HPP
#define BLOCKSQP_METHOD_HPP

#include "qpOASES.hpp"
#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"
#include "blocksqp_problemspec.hpp"
#include "blocksqp_options.hpp"
#include "blocksqp_iterate.hpp"
#include "blocksqp_stats.hpp"

namespace blockSQP
{

/**
 * \brief Describes an SQP method for a given problem and set of algorithmic options.
 * \author Dennis Janka
 * \date 2012-2015
 */
class SQPmethod
{
    /*
     * Variables
     */
    public:
        Problemspec*             prob;        ///< Problem structure (has to provide evaluation routines)
        SQPiterate*              vars;        ///< All SQP variables for this method
        SQPoptions*              param;       ///< Set of algorithmic options and parameters for this method
        SQPstats*                stats;       ///< Statistics object for current SQP run
        qpOASES::SQProblem*      qp;          ///< qpOASES qp object
        qpOASES::SQProblem*      qpSave;      ///< qpOASES qp object

    private:
        bool                     initCalled;  ///< indicates if init() has been called (necessary for run())

    /*
     * Methods
     */
    public:
        /// Construct a method for a given problem and set of algorithmic options
        SQPmethod( Problemspec *problem, SQPoptions *parameters, SQPstats *statistics );
        ~SQPmethod();
        /// Initialization, has to be called before run
        void init();
        /// Main Loop of SQP method
        int run( int maxIt, int warmStart = 0 );
        /// Call after the last call of run, to close output files etc.
        void finish();
        /// Print information about the SQP method
        void printInfo( int printLevel );
        /// Compute gradient of Lagrangian function (dense version)
        void calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, const Matrix &constrJacFull, Matrix &gradLagrange, int flag );
        /// Compute gradient of Lagrangian function (sparse version)
        void calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, double *jacNz, int *jacIndRow, int *jacIndCol, Matrix &gradLagrange, int flag );
        /// Overloaded function for convenience, uses current variables of SQPiterate vars
        void calcLagrangeGradient( Matrix &gradLagrange, int flag );
        /// Update optimization tolerance (similar to SNOPT) in current iterate
        bool calcOptTol();

        /*
         * Solve QP subproblem
         */
        /// Update the bounds on the current step, i.e. the QP variables
        void updateStepBounds( bool soc );
        /// Solve a QP with QPOPT or qpOASES to obtain a step deltaXi and estimates for the Lagrange multipliers
        int solveQP( Matrix &deltaXi, Matrix &lambdaQP, bool matricesChanged = true );
        /// Compute the next Hessian in the inner loop of increasingly convexified QPs and store it in vars->hess2
        void computeNextHessian( int idx, int maxQP );

        /*
         * Globalization Strategy
         */
        /// No globalization strategy
        int fullstep();
        /// Set new primal dual iterate
        void acceptStep( const Matrix &deltaXi, const Matrix &lambdaQP, double alpha, int nSOCS );
        /// Overloaded function for convenience, uses current variables of SQPiterate vars
        void acceptStep( double alpha );
        /// Reduce stepsize if a step is rejected
        void reduceStepsize( double *alpha );
        /// Determine steplength alpha by a filter based line search similar to IPOPT
        int filterLineSearch();
        /// Remove all entries from filter
        void initializeFilter();
        /// Is a pair (cNorm, obj) in the current filter?
        bool pairInFilter( double cNorm, double obj );
        /// Augment current filter by pair (cNorm, obj)
        void augmentFilter( double cNorm, double obj );
        /// Perform a second order correction step (solve QP)
        bool secondOrderCorrection( double cNorm, double cNormTrial, double dfTdeltaXi, bool swCond, int it );
        /// Reduce stepsize if a second order correction step is rejected
        void reduceSOCStepsize( double *alphaSOC );
        /// Start feasibility restoration heuristic
        int feasibilityRestorationHeuristic();
        /// Start feasibility restoration phase (solve NLP)
        int feasibilityRestorationPhase();
        /// Check if full step reduces KKT error
        int kktErrorReduction( );

        /*
         * Hessian Approximation
         */
        /// Set initial Hessian: Identity matrix
        void calcInitialHessian();
        /// [blockwise] Set initial Hessian: Identity matrix
        void calcInitialHessian( int iBlock );
        /// Reset Hessian to identity and remove past information on Lagrange gradient and steps
        void resetHessian();
        /// [blockwise] Reset Hessian to identity and remove past information on Lagrange gradient and steps
        void resetHessian( int iBlock );
        /// Compute current Hessian approximation by finite differences
        int calcFiniteDiffHessian();
        /// Compute full memory Hessian approximations based on update formulas
        void calcHessianUpdate( int updateType, int hessScaling );
        /// Compute limited memory Hessian approximations based on update formulas
        void calcHessianUpdateLimitedMemory( int updateType, int hessScaling );
        /// [blockwise] Compute new approximation for Hessian by SR1 update
        void calcSR1( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// [blockwise] Compute new approximation for Hessian by BFGS update with Powell modification
        void calcBFGS( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// Set pointer to correct step and Lagrange gradient difference in a limited memory context
        void updateDeltaGamma();

        /*
         * Scaling of Hessian Approximation
         */
        /// [blockwise] Update scalars for COL sizing of Hessian approximation
        void updateScalars( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// [blockwise] Size Hessian using SP, OL, or mean sizing factor
        void sizeInitialHessian( const Matrix &gamma, const Matrix &delta, int iBlock, int option );
        /// [blockwise] Size Hessian using the COL scaling factor
        void sizeHessianCOL( const Matrix &gamma, const Matrix &delta, int iBlock );
};

} // namespace blockSQP

#endif
