// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file QuadHandler.cpp
 * \brief Implement a handler for simple quadratic constraints of the form
 * \f$ y = x_1x_2 \f$,
 * and
 * \f$ y = x_1^2 \f$,
 * \author Ashutosh Mahajan, IIT Bombay
 */


#include <cmath>
#include <iostream>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrVarCand.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "LinBil.h"
#include "LinearFunction.h"
#include "LinMods.h"
#include "Logger.h"
#include "Node.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "ProblemSize.h"
#include "QuadHandler.h"
#include "QuadraticFunction.h"
#include "Relaxation.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "Variable.h"

// #define SPEW 1
// TODO:
// * arrange alphabetically
// * Move other classes to other files
// * remove rel->getRelaxationVar(). By default, these variables must be saved
//   in the data-structures.

using namespace Minotaur;

const std::string QuadHandler::me_ = "QuadHandler: ";

QuadHandler::QuadHandler(EnvPtr env, ProblemPtr problem)
: aTol_(1e-5),
  rTol_(1e-4)
{
  p_ = problem; 
  modProb_ = true;
  modRel_ = true;
  logger_  = (LoggerPtr) new Logger((LogLevel) 
                                    env->getOptions()->
                                    findInt("handler_log_level")->getValue());
  resetStats_();
  timer_ = env->getTimer();
}


QuadHandler::~QuadHandler()
{
  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    delete it->second;
  }
  x2Funs_.clear();

  for (LinBilSetIter it=x0x1Funs_.begin(); it!=x0x1Funs_.end(); ++it) {
    delete *it;
  }
  x0x1Funs_.clear();
}


void QuadHandler::addConstraint(ConstraintPtr newcon)
{
  LinearFunctionPtr lf0, lf1, lf;
  QuadraticFunctionPtr qf; 
  NonlinearFunctionPtr nlf;
  VariablePtr y, x0, x1;
  FunctionPtr f;
  LinSqrPtr lx2;
  LinBil* linbil;
  LinBilSetIter biter;

  cons_.push_back(newcon);
  qf = newcon->getQuadraticFunction();

  if (qf) {
    assert(!"cannot yet handle qf in QuadHandler.");
  } else {
    nlf = newcon->getNonlinearFunction();
    lf  = newcon->getLinearFunction();

    assert(lf && nlf);
    assert(1==lf->getNumTerms());
    assert(nlf->numVars()<3);

    y = lf->termsBegin()->first;
    if (nlf->numVars()==1) {
      x0 = *(nlf->varsBegin());
      lx2        = new LinSqr();
      lx2->y     = y;
      lx2->x     = x0;
      lx2->oeCon = ConstraintPtr(); // NULL
      x2Funs_.insert(std::pair<VariablePtr, LinSqrPtr>(x0, lx2));
    } else {
      VariableSet::iterator it = nlf->varsBegin();
      x0 = *it;
      ++it;
      x1 = *it;
      linbil = new LinBil(x0, x1, y);
      x0x1Funs_.insert(linbil);
    }
  }
}


void QuadHandler::findLinPt_(double xval, double yval, double &xl,
                             double &yl)
{
  // The point (xval, yval) satisfies yval < xval^2.
  // We want to find a new point (xl, yl) on the curve y=x^2 so that 
  // the gradient inequality at (xl, yl) cuts off (xval, yval).
  //
  // We will find a point (xl, yl) that is on the parabola y=x^2 and is
  // nearest to (xval, yval), i.e., xl solves
  // min_x (x - xval)^2  + (x^2 - yval)^2.
  // Solving it analytically is tedious as it involves solving a depressed
  // cubic equation. It may require using complex numbers.
  //
  // We merely find the solution using golden section search.
  
  double a, b, mu, la;
  double mu_val, la_val;
  const double alfa = 0.618;
  const double errlim = 1e-4; // don't want too much accuracy.

  if (xval>0) {
    a = sqrt(yval);
    b = xval;
  } else {
    a = xval;
    b = -sqrt(yval);
  }

  mu = a + alfa*(b-a);
  la = b - alfa*(b-a);
  mu_val = (mu-xval)*(mu-xval) + (mu*mu-yval)*(mu*mu-yval);
  la_val = (la-xval)*(la-xval) + (la*la-yval)*(la*la-yval);
  while ( (b-a) > errlim) {
    if (mu_val < la_val) {
      a = la;
      la = mu;
      la_val = mu_val;
      mu = a + alfa*(b-a);
      mu_val = (mu-xval)*(mu-xval) + (mu*mu-yval)*(mu*mu-yval);
    } else {
      b = mu;
      mu = la;
      mu_val = la_val;
      la = b - alfa*(b-a);
      la_val = (la-xval)*(la-xval) + (la*la-yval)*(la*la-yval);
    }
  }
  xl = la;
  yl = la*la;
}


Branches QuadHandler::getBranches(BrCandPtr cand, DoubleVector &x,
                                  RelaxationPtr rel, SolutionPoolPtr)
{
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  VariablePtr v2;
  double value = x[v->getIndex()];
  BranchPtr branch;
  VarBoundModPtr mod;
  Branches branches = (Branches) new BranchPtrVector();

  // can't branch on something that is at its bounds.
  assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);

  // down branch
  branch = (BranchPtr) new Branch();
  if (modProb_) {
    v2 = rel->getOriginalVar(v);
    mod = (VarBoundModPtr) new VarBoundMod(v2, Upper, value);
    branch->addPMod(mod);
  }
  if (modRel_) {
    mod = (VarBoundModPtr) new VarBoundMod(v, Upper, value);
    branch->addRMod(mod);
  }
  branch->setActivity(0.5);// TODO: set this correctly
  branches->push_back(branch);

  // up branch
  branch = (BranchPtr) new Branch();
  if (modProb_) {
    mod    = (VarBoundModPtr) new VarBoundMod(v2, Lower, value);
    branch->addPMod(mod);
  }
  if (modRel_) {
    mod    = (VarBoundModPtr) new VarBoundMod(v, Lower, value);
    branch->addRMod(mod);
  }
  branch->setActivity(0.5); // TODO: set this correctly
  branches->push_back(branch);

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "branching on " << v->getName()
                                       << " <= " << value << " or " 
                                       << " >= " << value << std::endl;
#endif
  return branches;
}


void QuadHandler::getBranchingCandidates(RelaxationPtr rel,
                                         const DoubleVector &x, ModVector &,
                                         BrVarCandSet &cands,
                                         BrCandVector &, bool &is_inf)
{
  double yval, x0val, x1val, udist, ddist;
  BrVarCandPtr br_can;
  VariablePtr x0, x1;
  std::pair<BrVarCandIter, bool> ret;
  bool check;

  is_inf = false;

  // First check if there is a candidate x0 that violates y <= x0^2.
  for (LinSqrMapIter s_it=x2Funs_.begin(); s_it!=x2Funs_.end(); ++s_it) {
    x0 = rel->getRelaxationVar(s_it->first);
    x0val = x[x0->getIndex()];
    x1 = rel->getRelaxationVar(s_it->second->y);
    yval = x[x1->getIndex()];
    if ((yval - x0val*x0val)/(fabs(yval)+1e-6)>1e-4) {
#if SPEW
      logger_->msgStream(LogDebug2) << std::setprecision(9) << me_ 
        << "branching candidate for x^2: " << s_it->first->getName()
        << " value = " << x0val << " aux var: " 
        << s_it->second->y->getName() 
        << " value = " << yval << std::endl;
#endif
      ddist = (yval - x0val*x0val)/
             sqrt(1.0+(x0->getLb()+x0val)*(x0->getLb()+x0val));
      udist = (yval - x0val*x0val)/
             sqrt(1.0+(x0->getUb()+x0val)*(x0->getUb()+x0val));
      br_can = (BrVarCandPtr) new BrVarCand(x0, x0->getIndex(), ddist,
                                            udist);
      ret = cands.insert(br_can);
      if (false == ret.second) { // already exists.
        br_can = *(ret.first);
        br_can->setDist(ddist+br_can->getDDist(), udist+br_can->getDDist());
      }
    }
  }

  // Now check if there is a violated constraint of the form y = x0.x1.
  // If so, add both x0 and x1 to the candidate set.
  for (LinBilSetIter s_it=x0x1Funs_.begin(); s_it!=x0x1Funs_.end(); ++s_it) {
    LinBil* bil = (*s_it);
    x0 = rel->getRelaxationVar(bil->getX0());
    x1 = rel->getRelaxationVar(bil->getX1());
    x0val = x[x0->getIndex()];
    x1val = x[x1->getIndex()];
    yval  = x[bil->getY()->getIndex()];
    if (bil->isViolated(x0val, x1val, yval)) {
      check = false;
      // If a variable is at bounds, then it is not a candidate.
      if (false==isAtBnds_(x0, x0val)) {
        check = true;
        if (x0val*x1val > yval) {
          ddist = (-yval + x0val*x1val)/sqrt(1.0+x0val*x0val+
                                            x1->getUb()*x1->getUb());
          udist = (-yval + x0val*x1val)/sqrt(1.0+x0val*x0val+
                                            x1->getLb()*x1->getLb());
        } else {
          ddist = (yval - x0val*x1val)/sqrt(1.0+x0val*x0val+
                                            x1->getLb()*x1->getLb());
          udist = (yval - x0val*x1val)/sqrt(1.0+x0val*x0val+
                                            x1->getUb()*x1->getUb());
        }
        br_can = (BrVarCandPtr) new BrVarCand(x0, x0->getIndex(), ddist,
                                              udist); 
        ret = cands.insert(br_can);
        if (false == ret.second) { // already exists.
          br_can = *(ret.first);
          br_can->setDist(ddist+br_can->getDDist(), udist+br_can->getUDist());
        }
      }

      if (false==isAtBnds_(x1, x1val)) {
        check = true;
        if (x0val*x1val > yval) {
          ddist = (-yval + x1val*x0val)/sqrt(1.0+x1val*x1val+
                                            x0->getUb()*x0->getUb());
          udist = (-yval + x1val*x0val)/sqrt(1.0+x1val*x1val+
                                            x0->getLb()*x0->getLb());
        } else {
          ddist = (yval - x1val*x0val)/sqrt(1.0+x1val*x1val+
                                            x0->getLb()*x0->getLb());
          udist = (yval - x1val*x0val)/sqrt(1.0+x1val*x1val+
                                            x0->getUb()*x0->getUb());
        }
        br_can = (BrVarCandPtr) new BrVarCand(x1, x1->getIndex(), ddist,
                                              udist); 
        ret = cands.insert(br_can);
        if (false == ret.second) { // already exists.
          br_can = *(ret.first);
          br_can->setDist(ddist+br_can->getDDist(), udist+br_can->getUDist());
        }
#if SPEW
      logger_->msgStream(LogDebug2) << std::setprecision(9) << me_ 
        << "branching candidate for x0x1: " << x0->getName()
        << " = " << x0val << " " << x1->getName() << " = " << x1val 
        << " " << bil->getY()->getName() << " = " << yval << " vio = " 
        << fabs(x0val*x1val-yval) << std::endl;
#endif
      }
      if (false==check) {
        int err;
        logger_->msgStream(LogError) << std::setprecision(15) << me_ 
         << "both variables are at bounds, but we still want to branch on "
         << "a bilinear constraint. " << x0->getName() 
         << " = " << x0val << " " << x1->getName() 
         << " = " << x1val << " " << bil->getY()->getName() 
         << " = " << yval  << " product = " << x0val*x1val << std::endl;
        logger_->msgStream(LogDebug) << bil->getC0()->getActivity(&(x[0]),
                                                                  &err) << " ";
        bil->getC0()->write(logger_->msgStream(LogDebug)); 

        logger_->msgStream(LogDebug) << bil->getC1()->getActivity(&(x[0]),
                                                                  &err) << " " ;
        bil->getC1()->write(logger_->msgStream(LogDebug));
        logger_->msgStream(LogDebug) << bil->getC2()->getActivity(&(x[0]),
                                                                  &err) << " " ;
        bil->getC2()->write(logger_->msgStream(LogDebug));
        logger_->msgStream(LogDebug) << bil->getC3()->getActivity(&(x[0]),
                                                                  &err) << " " ;
        bil->getC3()->write(logger_->msgStream(LogDebug));
      }
    }
  }
}


ModificationPtr QuadHandler::getBrMod(BrCandPtr cand, DoubleVector &x, 
                                      RelaxationPtr rel, BranchDirection dir)
{
  LinModsPtr lmods = (LinModsPtr) new LinMods();
  LinearFunctionPtr lf;
  LinConModPtr lmod;
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  VariablePtr x0, x1, y;
  double x0val, yval, vio, rhs;
  UInt vind = v->getIndex();
  bool found = false;
 
  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    if (vind==it->first->getIndex()) {
      x0val  = x[it->first->getIndex()];
      yval = x[it->second->y->getIndex()];
      vio = x0val*x0val-yval;
      if (vio > aTol_) {
        x0 = rel->getRelaxationVar(it->first);
        y = rel->getRelaxationVar(it->second->y);
        if (dir==DownBranch) {
          lf = getNewSqLf_(x0, y, x0->getLb(), x0val, rhs);
        } else {
          lf = getNewSqLf_(x0, y, x0val, x0->getUb(), rhs);
        }
        lmod = (LinConModPtr) new LinConMod(it->second->oeCon, lf,
                                            -INFINITY, rhs);
        lmods->insert(lmod);
      }
    }
  }
  for (LinBilSetIter it=x0x1Funs_.begin(); it != x0x1Funs_.end(); ++it) {
    if (vind==(*it)->getX0()->getIndex()) {
      if ((*it)->isViolated(&(x[0]), vio)) {
        x0 = rel->getRelaxationVar((*it)->getX0());
        x1 = rel->getRelaxationVar((*it)->getX1());
        found = true;
      }
    } else if (vind==(*it)->getX1()->getIndex()) {
      if ((*it)->isViolated(&(x[0]), vio)) {
        x0 = rel->getRelaxationVar((*it)->getX1());
        x1 = rel->getRelaxationVar((*it)->getX0());
        found = true;
      }
    }
    if (true==found) {
      y = rel->getRelaxationVar((*it)->getY());
      x0val  = x[x0->getIndex()];
      yval   = x[y->getIndex()];
      if (dir==DownBranch) {
        lf = getNewBilLf_(x0, x0->getLb(), x0val,
                          x1, x1->getLb(), x1->getUb(), y, 1, rhs);
        lmod = (LinConModPtr) new LinConMod((*it)->getC1(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
        lf = getNewBilLf_(x0, x0->getLb(), x0val,
                          x1, x1->getLb(), x1->getUb(), y, 3, rhs);
        lmod = (LinConModPtr) new LinConMod((*it)->getC3(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
      } else {
        lf = getNewBilLf_(x0, x0val, x0->getUb(),
                          x1, x1->getLb(), x1->getUb(), y, 0, rhs);
        lmod = (LinConModPtr) new LinConMod((*it)->getC0(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
        lf = getNewBilLf_(x0, x0val, x0->getUb(),
                          x1, x1->getLb(), x1->getUb(), y, 2, rhs);
        lmod = (LinConModPtr) new LinConMod((*it)->getC2(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
      }
      found = false;
    }
  }
  return lmods;
}


std::string QuadHandler::getName() const
{
   return "QuadHandler (Handling quadratic terms of the form y=x1*x2).";
}


LinearFunctionPtr QuadHandler::getNewBilLf_(VariablePtr x0, double lb0,
                                            double ub0, VariablePtr x1,
                                            double lb1, double ub1,
                                            VariablePtr y, int type,
                                            double &rhs)
{
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();

  switch (type) {
  case(0):
    // y >= l0x1 + l1x0 - l1l0
    lf->addTerm(x1, lb0);
    lf->addTerm(x0, lb1);
    lf->addTerm(y, -1.);
    rhs = lb0*lb1;
    break;
  case(1):
     // y >= u0x1 + u1x0 - u1u0
     lf->addTerm(x1, ub0);
     lf->addTerm(x0, ub1);
     lf->addTerm(y, -1.);
     rhs = ub0*ub1;
     break;
  case(2):
     // y <= u1x0 + l0x1 - l0u1
     lf->addTerm(x0, -1.0*ub1);
     lf->addTerm(x1, -1.0*lb0);
     lf->addTerm(y, 1.);
     rhs = -lb0*ub1;
     break;
  case(3):
     // y <= l1x0 + u0x1 - u0l1
     lf->addTerm(x0, -1.0*lb1);
     lf->addTerm(x1, -1.0*ub0);
     lf->addTerm(y, 1.);
     rhs = -ub0*lb1;
     break;

  default:
     assert(!"getNewBilLf_ called with wrong value of i");
     break;
  }
  return lf;
}


LinearFunctionPtr QuadHandler::getNewSqLf_(VariablePtr x, VariablePtr y,
                                           double lb, double ub, double & r)
{
  LinearFunctionPtr lf = LinearFunctionPtr(); // NULL
  r = -ub*lb;
  assert (lb > -1e21 && ub < 1e21); // Can't approximate when unbounded
  if (fabs(ub+lb) > aTol_) {
    lf = (LinearFunctionPtr) new LinearFunction();
    lf->addTerm(y, 1.);
    lf->addTerm(x, -1.*(ub+lb));
  } else {
    lf = (LinearFunctionPtr) new LinearFunction();
    lf->addTerm(y, 1.);
#if SPEW
    logger_->msgStream(LogDebug) << me_ << 
      "warning: generating a bound as a secant constraint." << std::endl;
#endif
  }
  return lf;
}


void QuadHandler::addCut_(VariablePtr x, VariablePtr y, 
                          double xl, double yl, double xval, double yval,
                          RelaxationPtr rel, bool &ifcuts)
{
  // add the cut 2*xl*x - y - yl <= 0
  ifcuts = false;
  if (2*xl*xval - yval - yl > 1e-5 &&
      2*xl*xval - yval > yl*(1+1e-4)) {
    LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
    FunctionPtr f;
    ConstraintPtr c;

    lf->addTerm(x, 2*xl);
    lf->addTerm(y, -1.0);
    f = (FunctionPtr) new Function(lf);
    c = rel->newConstraint(f, -INFINITY, xl*xl);
    ifcuts = true;
    ++sStats_.cuts;
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "new cut added" << std::endl;
    c->write(logger_->msgStream(LogDebug2));
#endif
  } else {
#if SPEW
    logger_->msgStream(LogDebug2) << me_ << "Not adding cut because of "
                                  << "insufficient violation "
                                  << 2*xl*xval - yval - xl*xl << std::endl;
#endif
  }
}


bool QuadHandler::isAtBnds_(ConstVariablePtr x, double xval)
{
  double lb = x->getLb();
  double ub = x->getUb();
  //double rtol = rTol_/10.0;
  return(fabs(xval-lb) < aTol_ || fabs(xval-ub) < aTol_);
         // fabs(xval-lb) < fabs(lb)*rtol || fabs(ub-xval) < fabs(ub)*rtol);
}


bool QuadHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr , bool &,
                             double &inf_meas)
{
  double yval, xval, vio;
  const double *x = sol->getPrimal();
  bool is_feas = true;

  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    // check if y <= x^2
    xval  = x[it->first->getIndex()];
    yval = x[it->second->y->getIndex()];
    vio = fabs(yval - xval*xval);
    if (vio/(fabs(yval) + 1e-6) > 1e-4 && vio > 1e-5) {
      is_feas = false;
      inf_meas += vio;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "no branching candidates for y=x^2" 
                                << std::endl;
#endif

  for (LinBilSetIter it=x0x1Funs_.begin(); it != x0x1Funs_.end(); ++it) {
    if ((*it)->isViolated(x, vio)) {
      is_feas = false;
      inf_meas += vio;
    }
  }

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "no branching candidates for y=x1x2" 
                                << std::endl;
#endif
  return is_feas;
}


SolveStatus QuadHandler::presolve(PreModQ *, bool *changed)
{

  bool is_inf = false;
  SolveStatus status = Finished;

  *changed = false;

  is_inf = varBndsFromCons_(changed);
  if (is_inf) {
    status = SolvedInfeasible;
  }

  if (Started==status) {
    status = Finished;
  }

  return status;
}


bool QuadHandler::presolveNode(RelaxationPtr rel, NodePtr, SolutionPoolPtr,
                               ModVector &p_mods, ModVector &r_mods)
{
  bool lchanged = true;
  bool changed = true;
  bool is_inf = false;
  double stime = timer_->query();

  // visit each quadratic constraint and see if bounds can be improved.
  while (true==changed) {
    ++pStats_.iters;
    changed = false;
    for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
      is_inf = propSqrBnds_(it, rel, modRel_, &lchanged, p_mods, r_mods);
      if (true==lchanged) {
        changed = true;
        lchanged = false;
      }
      if (is_inf) {
        return true;
      }
    }
    for (LinBilSetIter it=x0x1Funs_.begin(); it != x0x1Funs_.end(); ++it) {
      is_inf = propBilBnds_(*it, rel, modRel_, &lchanged, p_mods, r_mods);
      if (true==lchanged) {
        changed = true;
        lchanged = false;
      }
      if (is_inf) {
        return true;
      }
    }
  }

  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    upSqCon_(it->second->oeCon, rel->getRelaxationVar(it->first),
             rel->getRelaxationVar(it->second->y), rel, r_mods);
  }
  for (LinBilSetIter it=x0x1Funs_.begin(); it != x0x1Funs_.end(); ++it) {
    upBilCon_(*it, rel, r_mods);
  }

  if (modRel_) {
    pStats_.nMods += r_mods.size();
  } else {
    pStats_.nMods += p_mods.size();
  }
  //rel->write(std::cout);

  pStats_.timeN += timer_->query()-stime;
  return false;
}


bool QuadHandler::propBilBnds_(LinBil* lx0x1, RelaxationPtr rel,
                               bool mod_rel, bool *changed, ModVector &p_mods,
                               ModVector &r_mods)
{
  VariablePtr x0 = lx0x1->getX0();
  VariablePtr x1 = lx0x1->getX1();
  VariablePtr  y = lx0x1->getY();
  double lb, ub;
  VarBoundMod2Ptr b2mod;
  VarBoundModPtr bmod;

  BoundsOnProduct(x0, x1, lb, ub);
  //x0->write(std::cout);
  //x1->write(std::cout);
  // y->write(std::cout);
  //std::cout << lb << " " << ub << std::endl;
  if (updatePBounds_(y, lb, ub, rel, mod_rel, changed, p_mods, r_mods)<0) {
    return true; // infeasible
  }

  // other direction
  BoundsOnDiv(y->getLb(), y->getUb(), x0->getLb(), x0->getLb(), lb, ub);
  if (updatePBounds_(x1, lb, ub, rel, mod_rel, changed, p_mods, r_mods)<0) {
    return true; // infeasible
  }

  BoundsOnDiv(y->getLb(), y->getUb(), x0->getLb(), x0->getLb(), lb, ub);
  if (updatePBounds_(x0, lb, ub, rel, mod_rel, changed, p_mods, r_mods)<0) {
    return true; // infeasible
  }

  return false;
}

bool QuadHandler::propBilBnds_(LinBil* lx0x1, bool *changed)
{
  double lb, ub;
  VariablePtr x0   = lx0x1->getX0();     // x0 and x1 are variables in p_
  VariablePtr x1   = lx0x1->getX1();
  VariablePtr y    = lx0x1->getY();

  BoundsOnProduct(x0, x1, lb, ub);
  if (updatePBounds_(y, lb, ub, changed) < 0) {
    return true;
  } 

  // reverse
  BoundsOnDiv(y->getLb(), y->getUb(), x0->getLb(), x0->getLb(), lb, ub);
  if (updatePBounds_(x1, lb, ub, changed) < 0) {
    return true;
  } 

  BoundsOnDiv(y->getLb(), y->getUb(), x1->getLb(), x1->getLb(), lb, ub);
  if (updatePBounds_(x0, lb, ub, changed) < 0) {
    return true;
  } 

  return false;
}


bool QuadHandler::propSqrBnds_(LinSqrMapIter lx2, bool *changed)
{
  VariablePtr x   = lx2->first;      // x0 and y are variables in p_
  VariablePtr y   = lx2->second->y;
  double lb, ub;

  BoundsOnSquare(x, lb, ub);
  if (updatePBounds_(y, lb, ub, changed) < 0) {
    return true;
  } 

  // other direction.
  if (y->getUb()>aTol_) {
    ub = sqrt(y->getUb());
    lb = -ub;
    assert(y->getLb()>=0.0); // square of a number.
    if (x->getLb() > -sqrt(y->getLb())+aTol_) {
      lb = sqrt(y->getLb());
    }
    if (updatePBounds_(x, lb, ub, changed) < 0) {
      return true;
    } 
  } else if (y->getUb()<-aTol_) {
    return true;
  } else {
    if (updatePBounds_(x, 0.0, 0.0, changed) < 0) {
      return true;
    } 
  }
  return false;
}


bool QuadHandler::propSqrBnds_(LinSqrMapIter lx2, RelaxationPtr rel,
                               bool mod_rel, bool *changed, ModVector &p_mods,
                               ModVector &r_mods)
{
  double lb, ub;
  VarBoundMod2Ptr b2mod;
  VarBoundModPtr bmod;
  LinearFunctionPtr lf;

  VariablePtr x = lx2->first;      // x and y are variables in p_
  VariablePtr y = lx2->second->y;
  ConstraintPtr con = lx2->second->oeCon;

  BoundsOnSquare(x, lb, ub);
  if (updatePBounds_(y, lb, ub, rel, mod_rel, changed, p_mods, r_mods)<0) {
    return true; // infeasible
  }

  // other direction.
  if (y->getUb()>aTol_) {
    ub = sqrt(y->getUb());
    lb = -ub;
    assert(y->getLb()>=0.0);
    if (x->getLb() > -sqrt(y->getLb())+aTol_ ) {
      lb = sqrt(y->getLb());
    }
    if (updatePBounds_(x, lb, ub, rel, mod_rel, changed, p_mods, r_mods)<0) {
      return true; // infeasible
    }
  } else if (x->getUb()<-aTol_) {
    return true; // infeasible
  } else {
    if (updatePBounds_(x, 0.0, 0.0, rel, mod_rel, changed, p_mods, r_mods)<0) {
      return true; // infeasible
    }
  }

  return false;
}


void QuadHandler::relax_(RelaxationPtr rel, bool *)
{
  double rhs;
  LinearFunctionPtr lf;
  VariablePtr y, x0, x1;
  FunctionPtr f;
  ConstraintPtr c;
  ConstraintVector cons(4);


  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    x0 = rel->getRelaxationVar(it->second->x);
    y  = rel->getRelaxationVar(it->second->y);
    lf = getNewSqLf_(x0, y, x0->getLb(), x0->getUb(), rhs);

    f = (FunctionPtr) new Function(lf);
    it->second->oeCon = rel->newConstraint(f, -INFINITY, rhs);
  }

  for (LinBilSetIter it=x0x1Funs_.begin(); it!=x0x1Funs_.end(); ++it) {
    x0 = rel->getRelaxationVar((*it)->getX0());
    x1 = rel->getRelaxationVar((*it)->getX1());
    y  = rel->getRelaxationVar((*it)->getY());
    for (int i=0; i<4; ++i) {
      lf = getNewBilLf_(x0, x0->getLb(), x0->getUb(),
                        x1, x1->getLb(), x1->getUb(), y, i, rhs);
      f = (FunctionPtr) new Function(lf);
      cons[i] = rel->newConstraint(f, -INFINITY, rhs);
    }
    (*it)->setCons(cons[0], cons[1], cons[2], cons[3]);
  }

  //std::cout << "Writing the relaxation after quadratic relaxation:\n";
  //rel->write(std::cout);
  assert(0 == rel->checkConVars());
  return;
}


void QuadHandler::relaxInitFull(RelaxationPtr rel, bool *is_inf)
{

  relax_(rel, is_inf);
}


void QuadHandler::relaxInitInc(RelaxationPtr rel, bool *is_inf)
{
  relax_(rel, is_inf);
}


void QuadHandler::relaxNodeFull(NodePtr, RelaxationPtr, bool *)
{
  assert(!"QuadHandler::relaxNodeFull not implemented!");
}


void QuadHandler::relaxNodeInc(NodePtr, RelaxationPtr, bool *)
{
  // do nothing. Presolve will take care of tightening bounds
}


void QuadHandler::resetStats_()
{
  pStats_.iters  = 0;
  pStats_.time   = 0.0;
  pStats_.timeN  = 0.0;
  pStats_.vBnd   = 0;
  pStats_.nMods  = 0;
  
  sStats_.iters  = 0;
  sStats_.cuts   = 0;
  sStats_.time   = 0.0;
}


void QuadHandler::separate(ConstSolutionPtr sol, NodePtr , RelaxationPtr rel,
                           CutManager *, SolutionPoolPtr , bool *,
                           SeparationStatus *status)
{
  double yval, xval;
  double yl, xl;
  const double *x = sol->getPrimal();
  bool ifcuts;

  ++sStats_.iters;
  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    xval = x[it->first->getIndex()];
    yval = x[it->second->y->getIndex()];
    if (xval*xval > (1.0+1e-4)*fabs(yval) &&
       fabs(xval*xval-yval) > 1e-5) {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "xval = " << xval
                                    << " yval = " << yval << " violation = "
                                    << xval*xval-yval << std::endl;
#endif
      findLinPt_(xval, yval, xl, yl);
      addCut_(rel->getRelaxationVar(it->first), 
              rel->getRelaxationVar(it->second->y), xl, yl, xval, yval, 
              rel, ifcuts);
      if (true==ifcuts) {
        *status = SepaResolve;
      }
    }
  }
}



int QuadHandler::updatePBounds_(VariablePtr v, double lb, double ub,
                                bool *changed)
{
  if (ub < v->getLb() - aTol_ || lb > v->getUb() + aTol_) { 
    return -1;
  }
  
  if (ub < v->getUb() - aTol_) {
    p_->changeBound(v, Upper, ub);
    *changed = true;
    ++pStats_.vBnd;
  }
  if (lb > v->getLb() + aTol_) {
    p_->changeBound(v, Lower, lb);
    *changed = true;
    ++pStats_.vBnd;
  }

  return 0;
}


int QuadHandler::updatePBounds_(VariablePtr v, double lb, double ub,
                                RelaxationPtr rel, bool mod_rel,
                                bool *changed, ModVector &p_mods,
                                ModVector &r_mods)
{
  VarBoundMod2Ptr b2mod;
  VarBoundModPtr bmod;

  if (lb > v->getUb()+aTol_ || ub < v->getLb()-aTol_) {
    return -1;
  }

  if (lb > v->getLb()+aTol_ && ub < v->getUb()-aTol_) {
    *changed = true;
    ++pStats_.vBnd;
    b2mod  = (VarBoundMod2Ptr) new VarBoundMod2(v, lb, ub);
    b2mod->applyToProblem(p_);
    p_mods.push_back(b2mod);

    if (true == mod_rel) {
      b2mod = (VarBoundMod2Ptr)
              new VarBoundMod2(rel->getRelaxationVar(v), lb, ub);
      b2mod->applyToProblem(rel);
      r_mods.push_back(b2mod);
    }
  } else if (lb > v->getLb()+aTol_) {
    ++pStats_.vBnd;
    *changed = true;
    bmod  = (VarBoundModPtr) new VarBoundMod(v, Lower, lb);
    bmod->applyToProblem(p_);
    p_mods.push_back(bmod);

    if (true == mod_rel) {
      bmod = (VarBoundModPtr)
             new VarBoundMod(rel->getRelaxationVar(v), Lower, lb);
      bmod->applyToProblem(rel);
      r_mods.push_back(bmod);
    }
  } else if (ub < v->getUb()-aTol_) {
    ++pStats_.vBnd;
    *changed = true;
    bmod  = (VarBoundModPtr) new VarBoundMod(v, Upper, ub);
    bmod->applyToProblem(p_);
    p_mods.push_back(bmod);
    
    if (true == mod_rel) {
      bmod  = (VarBoundModPtr)
               new VarBoundMod(rel->getRelaxationVar(v), Upper, ub);
      bmod->applyToProblem(rel);
      r_mods.push_back(bmod);
    }
  } 
  return 0;
}


void QuadHandler::upBilCon_(LinBil* lx0x1, RelaxationPtr rel, ModVector
                            &r_mods)
{
  LinModsPtr   lmods;
  LinConModPtr lmod;
  LinearFunctionPtr lf;
  VariablePtr y  = rel->getRelaxationVar(lx0x1->getY());
  VariablePtr x0 = rel->getRelaxationVar(lx0x1->getX0());
  VariablePtr x1 = rel->getRelaxationVar(lx0x1->getX1());
  double l0 = x0->getLb();
  double u0 = x0->getUb();
  double l1 = x1->getLb();
  double u1 = x1->getUb();
  double a0, a1;
  double rhs;
  ConstraintPtr con;

  // all constraints in the relaxation are of (<= rhs) type.
  // y >= l1x0 + l0x1 - l1l0: binding at (l0, l1), (l0, u1), (u0,l1)
  con = lx0x1->getC0();
  lf = con->getLinearFunction();
  a0 = lf->getWeight(x0);
  a1 = lf->getWeight(x1);
  if (a0*l0 + a1*l1 - l0*l1 < con->getUb() - aTol_ ||
      a0*l0 + a1*u1 - l0*u1 < con->getUb() - aTol_ || 
      a0*u0 + a1*l1 - u0*l1 < con->getUb() - aTol_) {
    lf = getNewBilLf_(x0, l0, u0, x1, l1, u1, y, 0, rhs);
    lmod = (LinConModPtr) new LinConMod(con, lf, -INFINITY, rhs);
    lmod->applyToProblem(rel);
    r_mods.push_back(lmod);
  }

  // y >= u0x1 + u1x0 - u1u0: binding at (l0, u1), (u0, l1), (u0, u1)
  con = lx0x1->getC1();
  lf = con->getLinearFunction();
  a0 = lf->getWeight(x0);
  a1 = lf->getWeight(x1);
  if (a0*l0 + a1*u1 - l0*u1 < con->getUb() - aTol_ ||
      a0*u0 + a1*l1 - u0*l1 < con->getUb() - aTol_ || 
      a0*u0 + a1*u1 - u0*u1 < con->getUb() - aTol_) {
    lf = getNewBilLf_(x0, l0, u0, x1, l1, u1, y, 1, rhs);
    lmod = (LinConModPtr) new LinConMod(con, lf, -INFINITY, rhs);
    lmod->applyToProblem(rel);
    r_mods.push_back(lmod);
  }


  // y <= u1x0 + l0x1 - l0u1: binding at (l0, 11), (l0, u1), (u0, u1)
  con = lx0x1->getC2();
  lf = con->getLinearFunction();
  a0 = lf->getWeight(x0);
  a1 = lf->getWeight(x1);
  if (a0*l0 + a1*l1 + l0*l1 < con->getUb() - aTol_ ||
      a0*l0 + a1*u1 + l0*u1 < con->getUb() - aTol_ || 
      a0*u0 + a1*u1 + u0*u1 < con->getUb() - aTol_) {
    lf = getNewBilLf_(x0, l0, u0, x1, l1, u1, y, 2, rhs);
    lmod = (LinConModPtr) new LinConMod(con, lf, -INFINITY, rhs);
    lmod->applyToProblem(rel);
    r_mods.push_back(lmod);
  }

  // y <= l1x0 + u0x1 - u0l1: binding at (l0, l1), (u0, l1), (u0, u1)
  con = lx0x1->getC3();
  lf = con->getLinearFunction();
  a0 = lf->getWeight(x0);
  a1 = lf->getWeight(x1);
  if (a0*l0 + a1*l1 + l0*l1 < con->getUb() - aTol_ ||
      a0*u0 + a1*l1 + u0*l1 < con->getUb() - aTol_ || 
      a0*u0 + a1*u1 + u0*u1 < con->getUb() - aTol_) {
    lf = getNewBilLf_(x0, l0, u0, x1, l1, u1, y, 3, rhs);
    lmod = (LinConModPtr) new LinConMod(con, lf, -INFINITY, rhs);
    lmod->applyToProblem(rel);
    r_mods.push_back(lmod);
  }
}


void QuadHandler::upSqCon_(ConstraintPtr con, VariablePtr x, VariablePtr y, 
                           RelaxationPtr rel, ModVector &r_mods)
{
  LinearFunctionPtr lf = con->getLinearFunction();
  double a_x = lf->getWeight(x);
  double lb = x->getLb();
  double ub = x->getUb();
  LinConModPtr lmod;
  double rhs;

  assert(fabs(lf->getWeight(y) - 1.0) <= 1e-8);
  // y - (lb+ub)x <= -ub*lb
  if ((lb*lb + a_x*lb < con->getUb() - aTol_) ||
      (ub*ub + a_x*ub < con->getUb() - aTol_)) {
    lf = getNewSqLf_(x, y, x->getLb(), x->getUb(), rhs);
    lmod = (LinConModPtr) new LinConMod(con, lf, -INFINITY, rhs);
    lmod->applyToProblem(rel);
    r_mods.push_back(lmod);
  }
}


bool QuadHandler::varBndsFromCons_(bool *changed)
{
  bool is_inf = false;
  for (LinSqrMapIter it=x2Funs_.begin(); it != x2Funs_.end(); ++it) {
    is_inf = propSqrBnds_(it, changed);
    if (true==is_inf) {
      return true;
    }
  }
  for (LinBilSetIter it=x0x1Funs_.begin(); it != x0x1Funs_.end(); ++it) {
    is_inf = propBilBnds_(*it, changed);
    if (true==is_inf) {
      return true;
    }
  }
  return false;
}


void QuadHandler::writeStats(std::ostream &out) const
{
  out << me_ << "Statistics for presolve by QuadHandler:"        << std::endl
    << me_ << "Number of iterations           = "<< pStats_.iters  << std::endl
    << me_ << "Time taken in initial presolve = "<< pStats_.time   << std::endl
    << me_ << "Time taken in node presolves   = "<< pStats_.timeN  << std::endl
    << me_ << "Times variables tightened      = "<< pStats_.vBnd   << std::endl
    << me_ << "Changes in nodes               = "<< pStats_.nMods  << std::endl
    ;

  out << me_ << "Statistics for separation by QuadHandler:"        << std::endl
    << me_ << "Number of calls to separate    = "<< sStats_.iters  << std::endl
    << me_ << "Number of cuts added           = "<< sStats_.cuts   << std::endl
    << me_ << "Time taken in separation       = "<< sStats_.time   << std::endl
    ;
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
