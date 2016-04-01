// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Types.cpp
 * \author Ashutosh Mahajan, Argonne National Laboratory.
 * \brief Define methods specific to Types.h (like compare functions for
 * sorting etc.)
 */

#include "MinotaurConfig.h"

#include "BrCand.h"
#include "BrVarCand.h"
#include "Variable.h"
#include "Types.h"

using namespace Minotaur;

bool Minotaur::CompareVarBrCand::operator()(ConstBrVarCandPtr b1,
                                            ConstBrVarCandPtr b2) const
{
  return (b1->getPCostIndex() < b2->getPCostIndex());
}


bool Minotaur::CompareVariablePtr::operator()(ConstVariablePtr v1,
                                              ConstVariablePtr v2) const
{
  return (v1->getId() < v2->getId());
}

// Serdar added.
bool Minotaur::CompareValueVariablePair::operator() (VariableValuePair v1, 
                                                     VariableValuePair v2) const
{
  return (v2.second < v1.second);
}

bool Minotaur::CompareValueVarInc::operator() (VariableValuePair v1,
                                               VariableValuePair v2) const
{
  return (v1.second < v2.second);
}

bool Minotaur::CompareIntDouble::operator() (id id1, id id2) const
{
  return (id2.second < id1.second);
}

// Serdar ended.

bool Minotaur::CompareVariablePair::operator()(ConstVariablePair tv1,
                                               ConstVariablePair tv2) const 
{
  UInt min1 = tv1.first->getId();
  UInt min2 = tv2.first->getId();
  
  if (min1 == min2) {
    return (tv1.second->getId() < tv2.second->getId());
  }
  return (min1 < min2);
}


FunctionType Minotaur::funcTypesAdd(FunctionType f1, FunctionType f2)
{
  FunctionType type;
  switch (f1) {
   case Nonlinear:
     type = Nonlinear;
     break;
   case Polynomial:
     type = (f2 == Nonlinear) ?  Nonlinear : Polynomial;
     break;
   case Quadratic:
     type = (f2 == Constant || f2 == Linear || f2 == Quadratic) ? 
       Quadratic : f2;
     break;
   case Linear:
     type = (f2 == Constant || f2 == Linear) ?  Linear : f2;
     break;
   case Constant:
     type = f2;
     break;
   default:
     type = Nonlinear;
     break;
  }
  return type;
}


FunctionType Minotaur::funcTypesMult(FunctionType f1, FunctionType f2)
{
  FunctionType type;
  switch (f1) {
  case Nonlinear:
    type = Nonlinear;
    break;
  case Polynomial:
    type = (f2 == Nonlinear) ?  Nonlinear : Polynomial;
    break;
  case Quadratic:
    if (Constant==f2) {
      type = Quadratic;
    } else if (Linear==f2 || Quadratic==f2 || Polynomial==f2) {
      type = Polynomial;
    } else {
      type = Nonlinear;
    }
    break;
  case Linear:
    if (Constant==f2) {
      type = Linear;
    } else if (Linear==f2) {
      type = Quadratic;
    } else if (Quadratic==f2 || Polynomial==f2) {
      type = Polynomial;
    } else {
      type = Nonlinear;
    }
    break;
  case Constant:
    type = f2;
    break;
  default:
    type = Nonlinear;
    break;
  }
  return type;
}


std::string Minotaur::getFunctionTypeString(FunctionType f)
{
  switch(f) {
  case (Constant):
    return "constant";
    break;
  case (Linear):
    return "linear";
    break;
  case (Quadratic):
    return "quadratic";
    break;
  case (Polynomial):
    return "polynomial";
    break;
  case (Nonlinear):
    return "nonlinear";
    break;
  default:
    return "unknown function type";
  }
  return "unknown function type";
}


std::string Minotaur::getProblemTypeString(ProblemType p)
{
  switch(p) {
  case (LP):
    return "LP";
    break;
  case (MILP):
    return "MILP";
    break;
  case (QP):
    return "QP";
    break;
  case (MIQP):
    return "MIQP";
    break;
  case (QCQP):
    return "QCQP";
    break;
  case (MIQCQP):
    return "MIQCQP";
    break;
  case (POLYP):
    return "POLYP";
    break;
  case (MIPOLYP):
    return "MIPOLYP";
    break;
  case (NLP):
    return "NLP";
    break;
  case (MINLP):
    return "MINLP";
    break;
  default:
    return "UnknownProblem";
    break;
  }
}


std::string Minotaur::getSolveStatusString(SolveStatus s)
{
  switch(s) {
  case (NotStarted):
    return "Not started";
  case (Started):
    return "Started";
  case (Restarted):
    return "Restarted";
  case (SolvedOptimal):
    return "Optimal solution found";
  case (SolvedInfeasible):
    return "Detected infeasibility";
  case (SolvedUnbounded):
    return "Detected unboundedness of relaxation";
  case (SolvedGapLimit):
    return "Reached limit on gap";
  case (SolvedSolsLimit):
    return "Reached limit on number of solutions";
  case (IterationLimitReached):
    return "Reached iteration limit";
  case (Interrupted):
    return "Interrupted";
  case (TimeLimitReached):
    return "Reached time limit";
  case (SolLimitReached):
    return "Reached the limit on number of solutions";
  case (Finished):
    return "Finished for some other reason";
  default:
    break;
  }
  return "Unknown solve status";
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
