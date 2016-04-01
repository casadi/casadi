// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

// /**
// \file SecantMod.cpp
// \brief Implement the modification class SecantMod that modifies a linear
// constraint and the upper bound on the constraint, when the constraint is a
// secant approximation of a concave function. We also modify the bound on
// the auxiliary variable.
// \author Ashutosh Mahajan, Argonne National Laboratory
// */

#include <cmath>

#include "MinotaurConfig.h"
#include "Operations.h"
#include "SecantMod.h"
#include "Variable.h"

using namespace Minotaur;

SecantMod::SecantMod(ConstraintPtr con, LinearFunctionPtr new_lf,
                     double new_rhs, VariablePtr x, BoundType lu, double new_b,
                     VariablePtr y)
{
  double y_lb, y_ub, b2;
  if (lu==Lower) {
    b2 = x->getUb();
    BoundsOnSquare(new_b, b2, y_lb, y_ub);
  } else {
    b2 = x->getLb();
    BoundsOnSquare(b2, new_b, y_lb, y_ub);
  }
  ymod_ = (VarBoundMod2Ptr) new VarBoundMod2(y, y_lb, y_ub);
  xmod_ = (VarBoundModPtr) new VarBoundMod(x, lu, new_b);
  lmod_ = (LinConModPtr) new LinConMod(con, new_lf, -INFINITY, new_rhs);
}


SecantMod::~SecantMod()
{
  ymod_.reset();
  xmod_.reset();
  lmod_.reset();
}


VariablePtr SecantMod::getY()
{
  return ymod_->getVar();
}


void SecantMod::applyToProblem(ProblemPtr problem) 
{
  ymod_->applyToProblem(problem);
  xmod_->applyToProblem(problem);
  lmod_->applyToProblem(problem);
}


void SecantMod::undoToProblem(ProblemPtr problem) 
{
  ymod_->undoToProblem(problem);
  xmod_->undoToProblem(problem);
  lmod_->undoToProblem(problem);
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
