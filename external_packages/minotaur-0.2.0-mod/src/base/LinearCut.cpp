// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file LinearCut.cpp
 * \brief Implement the methods of LinearCut class. 
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Function.h"
#include "LinearCut.h"
#include "LinearFunction.h"
#include "Problem.h"

using namespace Minotaur;


LinearCut::LinearCut()
  : cons_(ConstraintPtr()),
    f_(FunctionPtr()),
    lb_(-INFINITY),
    lf_(LinearFunctionPtr()),
    ub_(INFINITY)
{
  info_.timesEnabled  = 0;
  info_.timesDisabled = 0;
  info_.lastEnabled   = 0;
  info_.lastDisabled  = 0;
}


LinearCut::LinearCut(LinearFunctionPtr lf, double lb, double ub)
  : cons_(ConstraintPtr()),
    lb_(lb),
    lf_(lf),
    ub_(ub)
{
  f_ = (FunctionPtr) new Function(lf_);
}


LinearCut::~LinearCut()
{

}


void LinearCut::applyToProblem(ProblemPtr problem) 
{
  cons_ = problem->newConstraint(f_, lb_, ub_, "linear_cut");
}


void LinearCut::undoToProblem(ProblemPtr)
{
}


void LinearCut::write(std::ostream &out) const
{
  out << lb_ << " <= ";
  lf_->write(out);
  out << " <= " << ub_;
}


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
