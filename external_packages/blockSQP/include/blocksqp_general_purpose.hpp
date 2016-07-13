/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_general_purpose.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of general purpose routines for matrix and vector computations.
 */

#ifndef GENERAL_PURPOSE_HPP
#define GENERAL_PURPOSE_HPP

#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"
#include "blocksqp_lapack.h"

namespace blockSQP
{

double l1VectorNorm( const Matrix &v );
double l2VectorNorm( const Matrix &v );
double lInfVectorNorm( const Matrix &v );

double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl, const Matrix &weights );
double l1ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
double l2ConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );
double lInfConstraintNorm( const Matrix &xi, const Matrix &constr, const Matrix &bu, const Matrix &bl );

double adotb( const Matrix &a, const Matrix &b );
void Atimesb( const Matrix &A, const Matrix &b, Matrix &result );
void Atimesb( double *Anz, int *AIndRow, int *AIndCol, const Matrix &b, Matrix &result );

int calcEigenvalues( const Matrix &B, Matrix &ev );
double estimateSmallestEigenvalue( const Matrix &B );
int inverse( const Matrix &A, Matrix &Ainv );

} // namespace blockSQP

#endif
