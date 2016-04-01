// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Relaxation.cpp
 * \brief Define base class Relaxation.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cassert>
#include <cmath>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Function.h"
#include "LinearFunction.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "QuadraticFunction.h"
#include "Relaxation.h"
#include "SOS.h"
#include "Variable.h"

using namespace Minotaur;

Relaxation::Relaxation()
  : p_(ProblemPtr()) // NULL
{
}


Relaxation::Relaxation(ProblemPtr problem)
: p_(problem)
{
  VariablePtr vcopy, v0, v1;
  VariableGroupConstIterator vit;

  ConstConstraintPtr cconstr;
  ConstraintPtr cCopy;

  LinearFunctionPtr lf, lf2;
  NonlinearFunctionPtr nl;
  QuadraticFunctionPtr qf, qf2;
  FunctionPtr fun;
  ObjectivePtr obj;
  UInt i;
  int err;
  VariableIterator vbeg;
  VarVector vvec;

  initialPt_ = 0; //NULL

  // add variables
  i = 0;
  for (VariableConstIterator it=p_->varsBegin(); it!=p_->varsEnd(); 
      ++it, ++i) {
    vcopy = (*it)->clone(i);
    setIndex_(vcopy, i);
    vars_.push_back(vcopy);
  }
  
  vbeg = vars_.begin();
  // add constraints
  for (ConstraintConstIterator cit=p_->consBegin(); cit!=p_->consEnd();
       ++cit) {
    cconstr = *cit;

    // linear part
    lf = cconstr->getLinearFunction();
    if (lf) {
      lf2 = lf->cloneWithVars(vbeg);
    } else {
      lf2 = LinearFunctionPtr(); // NULL
    }

    // quadratic part
    qf2 = QuadraticFunctionPtr(); // NULL
    qf = cconstr->getQuadraticFunction();
    if (qf) {
      qf2 = (QuadraticFunctionPtr) new QuadraticFunction(); 
      for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
          ++it) {
        v0 = vars_[it->first.first->getIndex()];
        v1 = vars_[it->first.second->getIndex()];
        qf2->incTerm(std::make_pair(v0, v1), it->second);
      }
    }

    err = 0;
    nl = cconstr->getNonlinearFunction();
    if (nl) {
      nl = nl->cloneWithVars(vars_.begin(), &err);
      if (err!=0) {
        // ugly hack.
        // pass pointer to the original nonlinear function.
        // fix this using a user option!
        nl = cconstr->getNonlinearFunction();
      }
    }

    // done.
    fun = (FunctionPtr) new Function(lf2, qf2, nl);
    cCopy = newConstraint(fun, cconstr->getLb(), cconstr->getUb(), 
                          cconstr->getName());
  }

  // add SOS1 constraints
  vvec.clear();
  for (SOSConstIterator it=p_->sos1Begin(); it!=p_->sos1End();
       ++it) {
    for (VariableConstIterator it2 = (*it)->varsBegin();
         it2!=(*it)->varsEnd(); ++it2) {
      vvec.push_back(vars_[(*it2)->getIndex()]);
    }
    newSOS((*it)->getNz(), SOS1, (*it)->getWeights(), vvec,
           (*it)->getPriority());
    vvec.clear();
  }

  // add SOS2 constraints
  vvec.clear();
  for (SOSConstIterator it=p_->sos2Begin(); it!=p_->sos2End();
       ++it) {
    for (VariableConstIterator it2 = (*it)->varsBegin();
         it2!=(*it)->varsEnd(); ++it2) {
      vvec.push_back(vars_[(*it2)->getIndex()]);
    }
    newSOS((*it)->getNz(), SOS2, (*it)->getWeights(), vvec,
           (*it)->getPriority());
    vvec.clear();
  }

  // add objective
  obj = p_->getObjective();
  lf = obj->getLinearFunction();
  if (lf) {
    lf2 = (LinearFunctionPtr) new LinearFunction();
    for (vit = lf->termsBegin(); vit != lf->termsEnd(); vit++) {
      lf2->addTerm(vars_[vit->first->getIndex()], vit->second);
    }
  } else {
    lf2 = LinearFunctionPtr(); // NULL
  }
  // quadratic part
  qf2 = QuadraticFunctionPtr(); // NULL
  qf = obj->getQuadraticFunction();
  if (qf) {
    qf2 = (QuadraticFunctionPtr) new QuadraticFunction(); 
    for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); ++it) {
      v0 = vars_[it->first.first->getIndex()];
      v1 = vars_[it->first.second->getIndex()];
      qf2->incTerm(std::make_pair(v0, v1), it->second);
    }
  }
  // pass pointer to the original nonlinear function.
  err = 0;
  nl = obj->getNonlinearFunction();
  if (nl) {
    nl = nl->cloneWithVars(vars_.begin(), &err);
    if (err!=0) {
      // ugly hack.
      // pass pointer to the original nonlinear function.
      // fix this using a user option!
      nl = obj->getNonlinearFunction();
    }
  }


  fun = (FunctionPtr) new Function(lf2, qf2, nl);
  obj_ = (ObjectivePtr) new Objective(fun, obj->getConstant(), Minimize, 
                                      obj->getName()); 

  nextCId_ = cons_.size();
  nextVId_ = vars_.size();
  nativeDer_ = p_->hasNativeDer();
}


VariablePtr Relaxation::getOriginalVar(VariablePtr r_var)
{
  if (p_->getNumVars() > r_var->getIndex()) {
    return p_->getVariable(r_var->getIndex());
  } 
  return VariablePtr();
}


VariablePtr Relaxation::getRelaxationVar(VariablePtr p_var)
{
  return vars_[p_var->getIndex()];
}


void Relaxation::setProblem(ProblemPtr p)
{
  p_ = p;
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
