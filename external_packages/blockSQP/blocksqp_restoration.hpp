/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_restoration.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of RestorationProblem class that describes a
 *  minimum l_2-norm NLP.
 */

#ifndef BLOCKSQP_RESTORATION_HPP
#define BLOCKSQP_RESTORATION_HPP

#include "blocksqp_defs.hpp"
#include "blocksqp_problemspec.hpp"

namespace blockSQP
{

/**
 * \brief Describes a minimum l_2-norm NLP for a given parent problem
 *        that is solved during the feasibility restoration phase.
 * \author Dennis Janka
 * \date 2012-2015
 */
class RestorationProblem : public Problemspec
{
    /*
     * CLASS VARIABLES
     */
    public:
        Problemspec *parent;
        Matrix xiRef;
        Matrix diagScale;
        int neq;
        bool *isEqCon;

        double zeta;
        double rho;

    /*
     * METHODS
     */
    public:
        RestorationProblem( Problemspec *parent, const Matrix &xiReference );

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (dense version)
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac );

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (sparse version)
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol );

        /// Evaluate all problem functions and their derivatives (dense version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, Matrix &constrJac,
                               SymMatrix *&hess, int dmode, int *info );

        /// Evaluate all problem functions and their derivatives (sparse version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info );

        virtual void printInfo();
        virtual void printVariables( const Matrix &xi, const Matrix &lambda, int verbose );
        virtual void printConstraints( const Matrix &constr, const Matrix &lambda );
};

} // namespace blockSQP

#endif

