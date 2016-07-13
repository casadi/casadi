/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_problemspec.cpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Implementation of problem-independent methods of ProblemSpec class.
 */

#include "blocksqp_problemspec.hpp"

namespace blockSQP
{

void Problemspec::evaluate( const Matrix &xi, double *objval, Matrix &constr, int *info )
{
    Matrix lambdaDummy, gradObjDummy;
    SymMatrix *hessDummy;
    int dmode = 0;

    Matrix constrJacDummy;
    double *jacNzDummy;
    int *jacIndRowDummy, *jacIndColDummy;
    *info = 0;

    // Try sparse version first
    evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, jacNzDummy, jacIndRowDummy, jacIndColDummy, hessDummy, dmode, info );

    // If sparse version is not implemented, try dense version
    if( info )
        evaluate( xi, lambdaDummy, objval, constr, gradObjDummy, constrJacDummy, hessDummy, dmode, info );
}

} // namespace blockSQP
