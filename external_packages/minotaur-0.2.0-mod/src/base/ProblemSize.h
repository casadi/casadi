//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file ProblemSize.h
 * \brief Declare structure ProblemSize for storing statistics about the problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURPROBLEMSIZE_H
#define MINOTAURPROBLEMSIZE_H

#include "Types.h"

namespace Minotaur {
  /// ProblemSize stores several statistics on the size of the Problem
  struct ProblemSize {
    /// Number of variables in the Problem.
    UInt vars;

    /// Number of constraints in the Problem.
    UInt cons;

    /// Number of objectives in the Problem.
    UInt objs;

    /// Number of binary variables.
    UInt bins;

    /// Number of integer variables not including binary variables
    UInt fixed;

    /// Number of integer variables not including binary variables
    UInt ints;

    /// Number of continuous variables
    UInt conts;

    /// Number of constraints that have a linear function only.
    UInt linCons;

    /// Number of SOS Type 1 constraints.
    UInt SOS1Cons;

    /// Number of SOS Type 2 constraints.
    UInt SOS2Cons;

    /**
     * Number of constraints that have a bilinear function and optionally, a
     * linear function.
     */
    UInt bilinCons;

    /**
     * Number of constraints that have a multilinear function and optionally, 
     * a linear function.
     */
    UInt multilinCons;

    /**
     * Number of constraints that have a quadratic function and optionally, a
     * linear function.
     */
    UInt quadCons;

    /// Number of constraints that have nonlinear function. 
    UInt nonlinCons;

    /// Number of constraints that have a linear function.
    UInt consWithLin;

    /// Number of constraints that have a bilinear function.
    UInt consWithBilin;

    /// Number of constraints that have a multilinear function.
    UInt consWithMultilin;

    /// Number of constraints that have a quadratic function.
    UInt consWithQuad;

    /// Number of constraints that have a nonlinear function.
    UInt consWithNonlin;

    /// Number of terms in the all linear functions in the constraints.
    UInt linTerms;

    /// Number of terms in the all multilinear functions in the constraints.
    UInt multiLinTerms;

    /// Number of terms in the all quadratic functions in the constraints.
    UInt quadTerms;

    /// Number of terms in linear function in the objectives.
    UInt objLinTerms;

    /// Number of terms in quadratic function in the objectives.
    UInt objQuadTerms;

    /// Type of objective: constant, linear, quadratic ...
    FunctionType objType;

  };
  typedef boost::shared_ptr< ProblemSize > ProblemSizePtr;
  typedef boost::shared_ptr< const ProblemSize > ConstProblemSizePtr;
}

#endif
    
// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
