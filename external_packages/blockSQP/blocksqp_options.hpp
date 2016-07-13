/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_options.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of SQPoptions class that holds all algorithmic options.
 */

#ifndef BLOCKSQP_OPTIONS_HPP
#define BLOCKSQP_OPTIONS_HPP

#include "blocksqp_defs.hpp"

namespace blockSQP
{

/**
 * \brief Contains algorithmic options and parameters for SQPmethod.
 * \author Dennis Janka
 * \date 2012-2015
 */
class SQPoptions
{
    /*
     * Variables
     */
    public:
        int printLevel;                     ///< information about the current iteration
        int printColor;                     ///< use colored terminal output
        int debugLevel;                     ///< amount of debug information that is printed during every iteration
        double eps;                         ///< values smaller than this are regarded as numerically zero
        double inf;                         ///< values larger than this are regarded as numerically infinity
        double opttol;                      ///< optimality tolerance
        double nlinfeastol;                 ///< nonlinear feasibility tolerance

        /* Algorithmic options */
        int sparseQP;                       ///< which qpOASES variant is used (dense/sparse/Schur)
        int globalization;                  ///< Globalization strategy
        int restoreFeas;                    ///< Use feasibility restoration phase
        int maxLineSearch;                  ///< Maximum number of steps in line search
        int maxConsecReducedSteps;          ///< Maximum number of consecutive reduced steps
        int maxConsecSkippedUpdates;        ///< Maximum number of consecutive skipped updates
        int maxItQP;                        ///< Maximum number of QP iterations per SQP iteration
        int blockHess;                      ///< Blockwise Hessian approximation?
        int hessScaling;                    ///< Scaling strategy for Hessian approximation
        int fallbackScaling;                ///< If indefinite update is used, the type of fallback strategy
        double maxTimeQP;                   ///< Maximum number of time in seconds per QP solve per SQP iteration
        double iniHessDiag;                 ///< Initial Hessian guess: diagonal matrix diag(iniHessDiag)
        double colEps;                      ///< epsilon for COL scaling strategy
        double colTau1;                     ///< tau1 for COL scaling strategy
        double colTau2;                     ///< tau2 for COL scaling strategy
        int hessDamp;                       ///< activate Powell damping for BFGS
        double hessDampFac;                 ///< damping factor for BFGS Powell modification
        int hessUpdate;                     ///< Type of Hessian approximation
        int fallbackUpdate;                 ///< If indefinite update is used, the type of fallback strategy
        int hessLimMem;                     ///< Full or limited memory
        int hessMemsize;                    ///< Memory size for L-BFGS updates
        int whichSecondDerv;                ///< For which block should second derivatives be provided by the user
        bool skipFirstGlobalization;        ///< If set to true, no globalization strategy in first iteration is applied
        int convStrategy;                   ///< Convexification strategy
        int maxConvQP;                      ///< How many additional QPs may be solved for convexification per iteration?

        /* Filter line search parameters */
        int maxSOCiter;                     ///< Maximum number of SOC line search iterations
        double gammaTheta;                  ///< see IPOPT paper
        double gammaF;                      ///< see IPOPT paper
        double kappaSOC;                    ///< see IPOPT paper
        double kappaF;                      ///< see IPOPT paper
        double thetaMax;                    ///< see IPOPT paper
        double thetaMin;                    ///< see IPOPT paper
        double delta;                       ///< see IPOPT paper
        double sTheta;                      ///< see IPOPT paper
        double sF;                          ///< see IPOPT paper
        double kappaMinus;                  ///< see IPOPT paper
        double kappaPlus;                   ///< see IPOPT paper
        double kappaPlusMax;                ///< see IPOPT paper
        double deltaH0;                     ///< see IPOPT paper
        double eta;                         ///< see IPOPT paper

    /*
     * Methods
     */
    public:
        SQPoptions();
        /// Some options cannot be used together. In this case set defaults
        void optionsConsistency();
};

} // namespace blockSQP

#endif
