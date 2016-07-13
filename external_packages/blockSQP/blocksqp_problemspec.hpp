/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_problemspec.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of ProblemSpec class to describes an NLP to be solved by blockSQP.
 */

#ifndef BLOCKSQP_PROBLEMSPEC_HPP
#define BLOCKSQP_PROBLEMSPEC_HPP

#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"

namespace blockSQP
{

/**
 * \brief Base class for problem specification as required by SQPmethod.
 * \author Dennis Janka
 * \date 2012-2015
 */
class Problemspec
{
    /*
     * VARIABLES
     */
    public:
        int         nVar;               ///< number of variables
        int         nCon;               ///< number of constraints
        int         nnCon;              ///< number of nonlinear constraints

        double      objLo;              ///< lower bound for objective
        double      objUp;              ///< upper bound for objective
        Matrix      bl;                 ///< lower bounds of variables and constraints
        Matrix      bu;                 ///< upper bounds of variables and constraints

        int         nBlocks;            ///< number of separable blocks of Lagrangian
        int*        blockIdx;           ///< [blockwise] index in the variable vector where a block starts

    /*
     * METHODS
     */
    public:
        Problemspec( ){};
        virtual ~Problemspec( ){};

        /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (dense version).
        virtual void initialize( Matrix &xi,            ///< optimization variables
                                 Matrix &lambda,        ///< Lagrange multipliers
                                 Matrix &constrJac      ///< constraint Jacobian (dense)
                                 ){};

        /// Set initial values for xi (and possibly lambda) and parts of the Jacobian that correspond to linear constraints (sparse version).
        virtual void initialize( Matrix &xi,            ///< optimization variables
                                 Matrix &lambda,        ///< Lagrange multipliers
                                 double *&jacNz,        ///< nonzero elements of constraint Jacobian
                                 int *&jacIndRow,       ///< row indices of nonzero elements
                                 int *&jacIndCol        ///< starting indices of columns
                                 ){};

        /// Evaluate objective, constraints, and derivatives (dense version).
        virtual void evaluate( const Matrix &xi,        ///< optimization variables
                               const Matrix &lambda,    ///< Lagrange multipliers
                               double *objval,          ///< objective function value
                               Matrix &constr,          ///< constraint function values
                               Matrix &gradObj,         ///< gradient of objective
                               Matrix &constrJac,       ///< constraint Jacobian (dense)
                               SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                               int dmode,               ///< derivative mode
                               int *info                ///< error flag
                               ){};

        /// Evaluate objective, constraints, and derivatives (sparse version).
        virtual void evaluate( const Matrix &xi,        ///< optimization variables
                               const Matrix &lambda,    ///< Lagrange multipliers
                               double *objval,          ///< objective function value
                               Matrix &constr,          ///< constraint function values
                               Matrix &gradObj,         ///< gradient of objective
                               double *&jacNz,          ///< nonzero elements of constraint Jacobian
                               int *&jacIndRow,         ///< row indices of nonzero elements
                               int *&jacIndCol,         ///< starting indices of columns
                               SymMatrix *&hess,        ///< Hessian of the Lagrangian (blockwise)
                               int dmode,               ///< derivative mode
                               int *info                ///< error flag
                               ){};

        /// Short cut if no derivatives are needed
        virtual void evaluate( const Matrix &xi,        ///< optimization variables
                               double *objval,          ///< objective function value
                               Matrix &constr,          ///< constraint function values
                               int *info                ///< error flag
                               );

        /*
         * Optional Methods
         */
        /// Problem specific heuristic to reduce constraint violation
        virtual void reduceConstrVio( Matrix &xi,       ///< optimization variables
                                      int *info         ///< error flag
                                      ){ *info = 1; };

        /// Print information about the current problem
        virtual void printInfo(){};
};

} // namespace blockSQP

#endif

