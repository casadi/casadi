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
 *  Declaration of SQPiterate class that holds all variables that are
 *  updated during one SQP iteration.
 */

#ifndef BLOCKSQP_ITERATE_HPP
#define BLOCKSQP_ITERATE_HPP

#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"
#include "blocksqp_problemspec.hpp"
#include "blocksqp_options.hpp"

namespace blockSQP
{

/**
 * \brief Holds all variables that are updated during one SQP iteration
 * \author Dennis Janka
 * \date 2012-2015
 */
class SQPiterate
{
    /*
     * Variables
     */
    public:
        double obj;                                   ///< objective value
        double qpObj;                                 ///< objective value of last QP subproblem
        double cNorm;                                 ///< constraint violation
        double cNormS;                                ///< scaled constraint violation
        double gradNorm;                              ///< norm of Lagrangian gradient
        double lambdaStepNorm;                        ///< norm of step in dual variables
        double tol;                                   ///< current optimality tolerance

        Matrix xi;                                    ///< variable vector
        Matrix lambda;                                ///< dual variables
        Matrix constr;                                ///< constraint vector

        Matrix constrJac;                             ///< full constraint Jacobian (not used in sparse mode)
        double *jacNz;                                ///< nonzero elements of Jacobian (length)
        int *jacIndRow;                               ///< row indices (length)
        int *jacIndCol;                               ///< indices to first entry of columns (nCols+1)

        Matrix deltaMat;                              ///< last m primal steps
        Matrix deltaXi;                               ///< alias for current step
        Matrix gradObj;                               ///< gradient of objective
        Matrix gradLagrange;                          ///< gradient of Lagrangian
        Matrix gammaMat;                              ///< Lagrangian gradient differences for last m steps
        Matrix gamma;                                 ///< alias for current Lagrangian gradient

        int nBlocks;                                  ///< number of diagonal blocks in Hessian
        int *blockIdx;                                ///< indices in the variable vector that correspond to diagonal blocks (nBlocks+1)

        SymMatrix *hess;                              ///< [blockwise] pointer to current Hessian of the Lagrangian
        SymMatrix *hess1;                             ///< [blockwise] first Hessian approximation
        SymMatrix *hess2;                             ///< [blockwise] second Hessian approximation (convexified)
        double *hessNz;                               ///< nonzero elements of Hessian (length)
        int *hessIndRow;                              ///< row indices (length)
        int *hessIndCol;                              ///< indices to first entry of columns (nCols+1)
        int *hessIndLo;                               ///< Indices to first entry of lower triangle (including diagonal) (nCols)

        /*
         * Variables for QP solver
         */
        Matrix deltaBl;                               ///< lower bounds for current step
        Matrix deltaBu;                               ///< upper bounds for current step
        Matrix lambdaQP;                              ///< dual variables of QP
        Matrix AdeltaXi;                              ///< product of constraint Jacobian with deltaXi

        /*
         * For modified BFGS updates
         */
        Matrix deltaNorm;                             ///< sTs
        Matrix deltaNormOld;                          ///< (from previous iteration)
        Matrix deltaGamma;                            ///< sTy
        Matrix deltaGammaOld;                         ///< (from previous iteration)
        int *noUpdateCounter;                         ///< count skipped updates for each block

        /*
         * Variables for globalization strategy
         */
        int steptype;                                 ///< is current step a restoration step (1)?
        double alpha;                                 ///< stepsize for line search
        int nSOCS;                                    ///< number of second-order correction steps
        int reducedStepCount;                         ///< count number of consecutive reduced steps,
        Matrix deltaH;                                ///< scalars for inertia correction (filter line search w indef Hessian)
        Matrix trialXi;                               ///< new trial iterate (for line search)
        std::set< std::pair<double,double> > *filter; ///< Filter contains pairs (constrVio, objective)

    /*
     * Methods
     */
    public:
        /// Call allocation and initializing routines
        SQPiterate( Problemspec* prob, SQPoptions* param, bool full );
        SQPiterate( const SQPiterate &iter );
        /// Allocate variables that any SQP code needs
        void allocMin( Problemspec* prob );
        /// Allocate diagonal block Hessian
        void allocHess( SQPoptions* param );
        /// Convert *hess to column compressed sparse format
        void convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_,
                             double *&hessNz_, int *&hessIndRow_, int *&hessIndCol_, int *&hessIndLo_ );
        /// Convert *hess to double array (dense matrix)
        void convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_ );
        /// Allocate variables specifically needed by vmused SQP method
        void allocAlg( Problemspec* prob, SQPoptions* param );
        /// Set initial filter, objective function, tolerances etc.
        void initIterate( SQPoptions* param );
        ~SQPiterate( void );
};

} // namespace blockSQP

#endif
