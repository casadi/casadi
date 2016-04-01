// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file CxQuadHandler.cpp
 * \brief Implement the handler for functions of the general quadratic form 
 * \f$ \sum_i x^TAx \leq b \f$,
 * where \f$A\f$ may be indefinite.
 * \author Ashutosh Mahajan, IIT Bombay
 */

/// TODO:
/// # remove asserts
/// # check for y - x^2 \leq k types of constraints
/// # check for x^2 \leq k types of constraints
/// # check for x^2 \geq k types of constraints
/// # check for x1x2 \geq k types of constraints
/// # check for x1x2 \leq k types of constraints
/// # many times, a secant inequality is just a bound

#include <cmath>
#include <iostream>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrVarCand.h"
#include "Constraint.h"
#include "CxQuadHandler.h"
#include "Environment.h"
#include "Function.h"
#include "LinMods.h"
#include "Logger.h"
#include "SecantMod.h"
#include "Node.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "QuadraticFunction.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "SolutionPool.h"
#include "Variable.h"

//#define SPEW 1

using namespace Minotaur;

const std::string CxQuadHandler::me_ = "CxQuadHandler: ";

CxQuadHandler::CxQuadHandler(EnvPtr env, ProblemPtr problem)
  : eTol_(1e-6)
{
  problem_ = problem; 
  logger_  = (LoggerPtr) new Logger((LogLevel) 
      env->getOptions()->findInt("handler_log_level")->getValue());
}


CxQuadHandler::~CxQuadHandler()
{
  for (VarSecantMapIter it=cvCons_.begin(); it != cvCons_.end(); ++it) {
    delete it->second;
  }
  cvCons_.clear();
  mcCons_.clear();
  brVars_.clear();
}


void CxQuadHandler::relax_(RelaxationPtr rel, bool *)
{
  ObjectivePtr oPtr;
  ConstraintConstIterator c_iter;
  ConstraintPtr cons, sec_cons;

  LinearFunctionPtr lf0, lf1, lf;
  QuadraticFunctionPtr qf, cx_qf0, cx_qf1; 
  FunctionPtr f;

  std::vector< VariablePtr > psqVars, nsqVars;
  VariablePtr v0, v1, v;

  //problem_->write(std::cout);
  // take care of objective.
  oPtr = problem_->getObjective();
  relaxObj_(oPtr, rel);
  //if (oPtr && oPtr->getQuadraticFunction()) {
  //  assert(!"can't handle quadratic objective function yet!");
  //}

  for (c_iter=problem_->consBegin(); c_iter!=problem_->consEnd(); ++c_iter) {
    cons = (*c_iter);
    if (cons->getNonlinearFunction()) {
      assert(!"can not handle nonlinear functions in createQuadRelax yet.");
    }
    qf = cons->getQuadraticFunction();
    if (qf) {
      assert(cons->getUb() < 1e15);
      // for a constraint of type l <= f <= u
      // cx_qf0 is the quadratic part of the relaxation of f <= u side.
      // cx_qf1 is the quadratic part of the relaxation of l <= f side.
      if (cons->getLb() > -1e15) {
        // constraint is ranged. It means non-convex.
        // visit each term of the quadratic function.
        relaxTwoSided_(qf, cons, rel);
      } else {
        // it is a <= constraint. We don't need cx_qf1, lf1.
        relaxOneSided_(qf, cons, rel);
      }
    }
  }
  //rel->write(std::cout);
  return;
}


void CxQuadHandler::relaxInitFull(RelaxationPtr rel, bool *is_inf)
{
  relax_(rel, is_inf);
}


void CxQuadHandler::relaxInitInc(RelaxationPtr rel, bool *is_inf)
{
  relax_(rel, is_inf);
}


void CxQuadHandler::relaxNodeFull(NodePtr, RelaxationPtr, bool *)
{
  assert(!"CxQuadHandler::relaxNodeFull not implemented!");
}


void CxQuadHandler::relaxNodeInc(NodePtr, RelaxationPtr, bool *)
{
  assert(!"CxQuadHandler::relaxNodeInc not implemented!");
}


void CxQuadHandler::relaxTwoSided_(QuadraticFunctionPtr qf, 
    ConstraintPtr cons, RelaxationPtr rel)
{
  VariablePtr v0, v1, v;
  double vlb, vub;
  FunctionPtr f;

  QuadraticFunctionPtr cx_qf0 = (QuadraticFunctionPtr) new QuadraticFunction();
  QuadraticFunctionPtr cx_qf1 = (QuadraticFunctionPtr) new QuadraticFunction();
  LinearFunctionPtr lf0 = (LinearFunctionPtr) new LinearFunction();
  LinearFunctionPtr lf1 = (LinearFunctionPtr) new LinearFunction();

  for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
      ++it) {
    v0 = rel->getVariable(it->first.first->getIndex());
    v1 = rel->getVariable(it->first.second->getIndex());
    if (v0==v1) {
      BoundsOnSquare(v0, vlb, vub);
      // convex and secant.
      if (it->second > 0.) {
        // convex part
        cx_qf0->incTerm(std::make_pair(v0, v1), it->second);
        // secant part
        v = rel->newVariable(vlb, vub, Continuous);
        lf1->incTerm(v, -1.0*it->second);
        addSecant_(v0, v, rel);
      } else {
        // convex part
        cx_qf1->incTerm(std::make_pair(v0, v1), -1.*it->second);
        // secant part
        v = rel->newVariable(vlb, vub, Continuous);
        lf0->incTerm(v, it->second);
        addSecant_(v0, v, rel);
      }
      brVars_.insert(v0);
    } else {
      // McCormick
      v = addMcCormick_(v0, v1, rel);
      lf0->incTerm(v, it->second);
      lf1->incTerm(v, -it->second);
      brVars_.insert(v0);
      brVars_.insert(v1);
    }
  }

  if (cx_qf0->getNumTerms()==0) {
    cx_qf0.reset();
  }
  if (cx_qf1->getNumTerms()==0) {
    cx_qf1.reset();
  }

  (*lf0) += cons->getLinearFunction();
  (*lf1) -= cons->getLinearFunction();

  if (!cx_qf0 && !cx_qf1) {
    // If there are no quadratic terms left, we need to add only one linear
    // constraint.
    f = (FunctionPtr) new Function(lf0);
    rel->newConstraint(f, cons->getLb(), cons->getUb());
  } else if (!cx_qf0 && cx_qf1) {
    f = (FunctionPtr) new Function(lf0);
    rel->newConstraint(f, cons->getLb(), cons->getUb());

    f = (FunctionPtr) new Function(lf1, cx_qf1);
    rel->newConstraint(f, -INFINITY, -1.*cons->getLb());
  } else {
    if (cx_qf0) {
      f = (FunctionPtr) new Function(lf0, cx_qf0);
      rel->newConstraint(f, -INFINITY, cons->getUb());
    } else {
      f = (FunctionPtr) new Function(lf1);
      rel->newConstraint(f, -1.*cons->getUb(), -1.*cons->getLb());
    }

    if (cx_qf1) {
      f = (FunctionPtr) new Function(lf1, cx_qf1);
      rel->newConstraint(f, -INFINITY, -1.*cons->getLb());
    } else {
      f = (FunctionPtr) new Function(lf1);
      rel->newConstraint(f, -cons->getUb(), -cons->getLb());
    }
  }
}


void CxQuadHandler::relaxOneSided_(QuadraticFunctionPtr qf, 
    ConstraintPtr cons, RelaxationPtr rel)
{
  VariablePtr v0, v1, v;
  double vlb, vub;
  FunctionPtr f;

  QuadraticFunctionPtr cx_qf0 = (QuadraticFunctionPtr) new QuadraticFunction();
  LinearFunctionPtr lf0 = (LinearFunctionPtr) new LinearFunction();


  for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
      ++it) {
    v0 = rel->getVariable(it->first.first->getIndex());
    v1 = rel->getVariable(it->first.second->getIndex());
    if (v0==v1) {
      if (it->second > 0.) {
        // convex only 
        cx_qf0->incTerm(std::make_pair(v0, v0), it->second);
      } else {
        // secant only 
        BoundsOnSquare(v0, vlb, vub);
        v = rel->newVariable(vlb, vub, Continuous);
        lf0->incTerm(v, it->second);
        addSecant_(v0, v, rel);
      }
    } else {
      // McCormick
      if (it->second > 0.) {
        // we only need lower-McCormick for y >= x1x2
        v = addMcCormickLower_(v0, v1, rel);
      } else {
        // we only need upper-McCormick for y <= x1x2
        v = addMcCormickUpper_(v0, v1, rel);
      }
      assert(v); // v should not be NULL
      lf0->incTerm(v, it->second);
    }
  }

  if (cx_qf0->getNumTerms()==0) {
    cx_qf0.reset();
  }
  (*lf0) += cons->getLinearFunction();
  f = (FunctionPtr) new Function(lf0, cx_qf0);
  rel->newConstraint(f, -INFINITY, cons->getUb());
}


void CxQuadHandler::relaxObj_(ObjectivePtr obj, RelaxationPtr rel)
{
  QuadraticFunctionPtr qf, cx_qf0;
  LinearFunctionPtr lf0, lf1;
  FunctionPtr f;
  double vlb, vub;
  VariablePtr v,v0,v1;
  if (!obj) {
   return;
  } 

  qf = obj->getQuadraticFunction();
  if (!qf) {
    return;
  }

  lf0 = (LinearFunctionPtr) new LinearFunction();
  cx_qf0 = (QuadraticFunctionPtr) new QuadraticFunction();
  for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
      ++it) {
    v0 = rel->getVariable(it->first.first->getIndex());
    v1 = rel->getVariable(it->first.second->getIndex());
    if (v0==v1) {
      if (it->second > 0.) {
        // convex only 
        cx_qf0->incTerm(std::make_pair(v0, v0), it->second);
      } else {
        // secant only 
        BoundsOnSquare(v0, vlb, vub);
        v = rel->newVariable(vlb, vub, Continuous);
        lf0->incTerm(v, it->second);
        addSecant_(v0, v, rel);
      }
    } else {
      // McCormick
      if (it->second > 0.) {
        // we only need lower-McCormick for y >= x1x2
        v = addMcCormickLower_(v0, v1, rel);
      } else {
        // we only need upper-McCormick for y <= x1x2
        v = addMcCormickUpper_(v0, v1, rel);
      }
      assert(v); // v should not be NULL
      lf0->incTerm(v, it->second);
    }
  }
  if (cx_qf0->getNumTerms()==0) {
    cx_qf0.reset();
  }
  lf1 = obj->getLinearFunction()->cloneWithVars(rel->varsBegin()); 
  (*lf0) += lf1;
  f = (FunctionPtr) new Function(lf0, cx_qf0);
  rel->newObjective(f, obj->getConstant(), Minimize);
}


void CxQuadHandler::addSecant_(VariablePtr x, VariablePtr y, RelaxationPtr rel)
{
  LinearFunctionPtr lf;
  SecantPtr sec;
  ConstraintPtr cons;
  FunctionPtr f;

  double lb = x->getLb();
  double ub = x->getUb();
  double r;

  lf = getNewSecantLf_(x, y, lb, ub, r);

  f = (FunctionPtr) new Function(lf);
  cons = rel->newConstraint(f, -INFINITY, r);

  sec         = new Secant();
  sec->auxVar = y;
  sec->sqVar  = x;
  sec->cons   = cons;

  cvCons_.insert(std::pair<VariablePtr, SecantPtr>(x, sec));
}



LinearFunctionPtr CxQuadHandler::getNewSecantLf_(VariablePtr x, VariablePtr y,
                                                 double & lb, double & ub,
                                                 double & r)
{
  LinearFunctionPtr lf = LinearFunctionPtr(); // NULL
  r = -ub*lb;
  assert (lb > -1e15 && ub < 1e15); // Can't approximate when unbounded
  if (fabs(ub+lb) > eTol_) {
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


VariablePtr CxQuadHandler::addMcCormick_(VariablePtr x0, VariablePtr x1,  
    RelaxationPtr rel)
{
  // add all four McCormick inequalities
  VariablePtr y = addMcCormickLower_(x0, x1, rel);
  y = addMcCormickUpper_(x0, x1, rel);
  return y;
}


VariablePtr CxQuadHandler::addMcCormickUpper_(VariablePtr x0, VariablePtr x1,
    RelaxationPtr rel)
{
  // add two McCormick inequalities for \f$  y \leq x0.x1 \f$.
  McCormickPtr mcc = (McCormickPtr) new McCormick(x0, x1, McCormick::GT);
  McCormickSetIter biter = mcCons_.find(mcc);
  double lb = 0, ub = 0;
  LinearFunctionPtr lf;
  FunctionPtr f;
  ConstraintPtr cons;
  VariablePtr y = VariablePtr();   // NULL
  bool exists = true;

  if (biter==mcCons_.end()) {
    exists = false;
  }

  if (exists) {
    mcc = (*biter);
    y = mcc->getAux();
  } else {
    // add a new bilinear term in our list
    BoundsOnProduct(x0, x1, lb, ub);
    y = rel->newVariable(lb, ub, Continuous);
  }

  assert(y);
  if (!exists || mcc->getSense() == McCormick::LT) {
    // we don't already have constraints for x0x1 <= y. Add them
    // y <= u1x0 + l0x1 - l0u1
    lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), 
        x1->getUb(), y, ub, 2);
    f = (FunctionPtr) new Function(lf);
    cons = rel->newConstraint(f, -INFINITY, ub);
    mcc->setC2(cons);

    // y <= l1x0 + u0x1 - u0l1
    lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), 
        x1->getUb(), y, ub, 3);
    f = (FunctionPtr) new Function(lf);
    cons = rel->newConstraint(f, -INFINITY, ub);
    mcc->setC3(cons);

    if (exists) {
      // we have both kinds of inequalities, so change the sense.
      mcc->setSense(McCormick::EQ);
    } else {
      mcc->setAux(y);
      mcc->setSense(McCormick::GT);
      // add it to the set.
      mcCons_.insert(mcc);
    }
  } else {
    // no need to add any new variables or new constraints.
    y = mcc->getAux();
  }
  return y;
}


VariablePtr CxQuadHandler::addMcCormickLower_(VariablePtr x0, VariablePtr x1, 
    RelaxationPtr rel)
{
  // add two McCormick inequalities for \f$ y \geq x0.x1 \f$.
  McCormickPtr mcc = (McCormickPtr) new McCormick(x0, x1, McCormick::LT);
  McCormickSetIter biter = mcCons_.find(mcc);
  double lb = 0, ub = 0;
  LinearFunctionPtr lf;
  FunctionPtr f;
  ConstraintPtr cons;
  VariablePtr y = VariablePtr();   // NULL
  bool exists = true;

  if (biter==mcCons_.end()) {
    exists = false;
  }

  if (exists) {
    mcc = *(biter);
    y = mcc->getAux();
  } else {
    // add a new bilinear term in our list
    BoundsOnProduct(x0, x1, lb, ub);
    y = rel->newVariable(lb, ub, Continuous);
  }

  assert(y);
  if (!exists || mcc->getSense() == McCormick::GT) {
    // we don't already have constraints for x0x1 <= y. Add them
    // y >= l0x1 + l1x0 - l1l0
    lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
        y, ub, 0);
    f = (FunctionPtr) new Function(lf);
    cons = rel->newConstraint(f, -INFINITY, ub);
    mcc->setC0(cons);

    // y >= u0x1 + u1x0 - u1u0
    lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
        y, ub, 1);
    f = (FunctionPtr) new Function(lf);
    cons = rel->newConstraint(f, -INFINITY, ub);
    mcc->setC1(cons);

    if (exists) {
      // we have both kinds of inequalities, so change the sense.
      mcc->setSense(McCormick::EQ);
    } else {
      mcc->setAux(y);
      mcc->setSense(McCormick::LT);
      // add it to the set.
      mcCons_.insert(mcc);
    }
  } else {
    // no need to add any new variables or new constraints.
    y = mcc->getAux();
  }
  return y;
}


LinearFunctionPtr CxQuadHandler::getMcLf_(VariablePtr x0, double lb0, double ub0,
                                          VariablePtr x1, double lb1, double ub1,
                                          VariablePtr y, double &rhs, UInt i)
{
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  assert(y->getLb() > -1e15 && y->getUb() < 1e15);
  switch(i) {
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
     assert(!"getMcLf_ called with wrong value of i");
     break;
  }
  return lf;
}


bool CxQuadHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr, bool &,
                               double &)
{
  double yval, xval;
  const double *x = sol->getPrimal();

  for (VarSecantMapIter it=cvCons_.begin(); it != cvCons_.end(); ++it) {
    // check if y <= x^2
    xval  = x[it->first->getIndex()];
    yval = x[it->second->auxVar->getIndex()];
    if (yval > xval*xval + eTol_) {
      return false;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "relaxation is secant feasible." 
    << std::endl;
#endif

  for (McCormickSetIter it=mcCons_.begin(); it != mcCons_.end(); ++it) {
    if ((*it)->isViolated(x, eTol_)) {
      return false;
    }
  }

#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "relaxation is McCormick feasible." 
    << std::endl;
#endif
  return true;
}


void CxQuadHandler::getBranchingCandidates(RelaxationPtr, 
                                           const DoubleVector &x, ModVector &,
                                           BrVarCandSet &cands, BrCandVector &,
                                           bool &is_inf)
{
  double yval, x0val, x1val;
  BrVarCandPtr br_can;
  VariablePtr x0, x1;
  UIntSet cand_inds;
  std::pair<UIntSet::iterator, bool> ret;
#if DEBUG
  bool check;
#endif

  is_inf = false;

  // First check if there is a candidate x0 that violates y <= x0^2.
  for (VarSecantMapIter s_it=cvCons_.begin(); s_it!=cvCons_.end(); ++s_it) {
    x0val = x[s_it->first->getIndex()];
    yval = x[s_it->second->auxVar->getIndex()];
    if (yval > x0val*x0val+eTol_) {
#if SPEW
      logger_->msgStream(LogDebug2) << std::setprecision(9) << me_ 
        << "branching candidate for secant: " << s_it->first->getName()
        << " value = " << x0val << " aux var: " 
        << s_it->second->auxVar->getName() 
        << " value = " << yval << std::endl;
      //s_it->first->write(std::cout);
      //s_it->second->auxVar->write(std::cout);
#endif
      ret = cand_inds.insert(s_it->first->getIndex());
      if (true == ret.second) { // does not already exist.
        br_can = (BrVarCandPtr) new BrVarCand(s_it->first, 
                                              s_it->first->getIndex(), 0.5,
                                              0.5);
        cands.insert(br_can);
      }
    }
  }

  // Now check if there is a violated constraint of the form y <=/=/=> x0.x1.
  // If so, add both x0 and x1 to the candidate set.
  for (McCormickSetIter s_it=mcCons_.begin(); s_it!=mcCons_.end(); ++s_it) {
    McCormickPtr bil = (*s_it);
    x0 = bil->getX0();
    x1 = bil->getX1();
    x0val = x[x0->getIndex()];
    x1val = x[x1->getIndex()];
    yval  = x[bil->getAux()->getIndex()];
    if (bil->isViolated(x0val, x1val, yval, eTol_)) {
#if DEBUG
      check = false;
#endif
      // If a variable is at bounds, then it is not a candidate.
      if (x0val < x0->getUb() - eTol_ && x0val > x0->getLb() + eTol_) {
#if DEBUG
        check = true;
#endif
        ret = cand_inds.insert(x0->getIndex());
        if (true == ret.second) { // does not already exist.
          br_can = (BrVarCandPtr) new BrVarCand(x0, x0->getIndex(), 0.5, 0.5); 
          cands.insert(br_can);
        }
      }

      if (x1val < x1->getUb() - eTol_ && x1val > x1->getLb() + eTol_) {
#if DEBUG
        check = true;
#endif
        ret = cand_inds.insert(x1->getIndex());
        if (true == ret.second) { // does not already exist.
          br_can = (BrVarCandPtr) new BrVarCand(x1, x1->getIndex(), 0.5, 0.5); 
          cands.insert(br_can);
        }
      }
      //if (false==check) {
      //  std::cout << std::setprecision(9) << x0->getName() 
      //   << " = " << x0val << " " << x1->getName() 
      //   << " = " << x1val << " " << bil->getAux()->getName() 
      //   << " = " << yval << std::endl;
      //  rel->write(std::cout);
      //}
#if DEBUG
      assert (check); // If both variables are at bounds, the McCormick 
                      // inequalities can not be violated.
#endif
    }
  }
}


ModificationPtr CxQuadHandler::getBrMod(BrCandPtr cand, DoubleVector &xval, 
    RelaxationPtr , BranchDirection dir)
{
  double            lb, ub, lb1, ub1, b2, rhs=0;
  BoundType         lu;
  ConstraintPtr     cons;
  BrVarCandPtr      vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr       x0, x1, v, y;
  VarSecantMapIter  s_it;
  LinearFunctionPtr lf;
  McCormickPtr      mcc;
  SecantModPtr      smod;
  LinModsPtr        lmods;
  LinConModPtr      lmod;
  VarBoundModPtr    bmod;
  VarBoundMod2Ptr   b2mod;

  x0 = vcand->getVar();

  if (dir == DownBranch ) {
    lb    = x0->getLb();
    ub    = xval[x0->getIndex()];
    b2    = ub;
    lu    = Upper;
  } else {
    lb    = xval[x0->getIndex()];
    ub    = x0->getUb();
    b2    = lb;
    lu    = Lower;
  }

  // first find if we have secants associated with x0
  s_it = cvCons_.find(x0);
  if (s_it!=cvCons_.end()) {
    y = s_it->second->auxVar;
    cons = s_it->second->cons;
    lf    = getNewSecantLf_(x0, y, lb, ub, rhs);
    smod = (SecantModPtr) new SecantMod(cons, lf, rhs, x0, lu, b2, y);
    return smod;
  } 

  // also try to find any McCormick inequalities associated with x0
  lmods = (LinModsPtr) new LinMods();
  for (McCormickSetIter it=mcCons_.begin(); it != mcCons_.end(); ++it) {
    mcc = (*it);
    x1 = mcc->getOtherX(x0);
    if (x1) {
      // This term contains x0 and x1. 
      y = mcc->getAux();
      lmods = (LinModsPtr) new LinMods();
      if (mcc->getSense() == McCormick::LT || 
          mcc->getSense() == McCormick::EQ) {
        // y >= l0x1 + l1x0 - l1l0
        lf = getMcLf_(x0, lb, ub, x1, x1->getLb(), x1->getUb(), y, rhs, 0);
        lmod = (LinConModPtr) new LinConMod(mcc->getC0(), lf, -INFINITY, rhs);
        lmods->insert(lmod);

        // y >= u0x1 + u1x0 - u1u0
        lf = getMcLf_(x0, lb, ub, x1, x1->getLb(), x1->getUb(), y, rhs, 1);
        lmod = (LinConModPtr) new LinConMod(mcc->getC1(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
      }

      if (mcc->getSense() == McCormick::GT || 
          mcc->getSense() == McCormick::EQ) {
        // y <= u1x0 + l0x1 - l0u1
        lf = getMcLf_(x0, lb, ub, x1, x1->getLb(), x1->getUb(), y, rhs, 2);
        lmod = (LinConModPtr) new LinConMod(mcc->getC2(), lf, -INFINITY, rhs);
        lmods->insert(lmod);

        // y <= l1x0 + u0x1 - u0l1
        lf = getMcLf_(x0, lb, ub, x1, x1->getLb(), x1->getUb(), y, rhs, 3);
        lmod = (LinConModPtr) new LinConMod(mcc->getC3(), lf, -INFINITY, rhs);
        lmods->insert(lmod);
      }
      BoundsOnProduct(lb, ub, x1->getLb(), x1->getUb(), lb1, ub1);
      b2mod  = (VarBoundMod2Ptr) new VarBoundMod2(y, lb1, ub1);
      lmods->insert(b2mod);
    }
  }
  bmod = (VarBoundModPtr) new VarBoundMod(x0, lu, b2);
  lmods->insert(bmod);
  return lmods;
}


Branches CxQuadHandler::getBranches(BrCandPtr cand, DoubleVector & x,
                                    RelaxationPtr rel, SolutionPoolPtr)
{
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  VariablePtr v2;
  double value = x[v->getIndex()];
  BranchPtr branch;
  Branches branches = (Branches) new BranchPtrVector();
  VarBoundModPtr mod;

  // can't branch on something that is at its bounds.
  assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);

  // down branch
  branch = (BranchPtr) new Branch();
  if (modProb_) {
    mod = (VarBoundModPtr) new VarBoundMod(v, Upper, value);
    branch->addPMod(mod);
  }
  if (modRel_) {
    v2 = rel->getRelaxationVar(v);
    mod = (VarBoundModPtr) new VarBoundMod(v2, Upper, value);
    branch->addRMod(mod);
  }
  branch->setActivity(0.5);// TODO: set this correctly
  branches->push_back(branch);

  // up branch
  branch = (BranchPtr) new Branch();
  if (modProb_) {
    mod    = (VarBoundModPtr) new VarBoundMod(v, Lower, value);
    branch->addPMod(mod);
  }
  if (modRel_) {
    v2 = rel->getRelaxationVar(v);
    mod = (VarBoundModPtr) new VarBoundMod(v2, Lower, value);
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


void CxQuadHandler::separate(ConstSolutionPtr, NodePtr , RelaxationPtr ,
                             CutManager *, SolutionPoolPtr , bool *,
                             SeparationStatus *)
{

}



bool CxQuadHandler::presolveNode(RelaxationPtr rel, NodePtr, SolutionPoolPtr,
                                 ModVector &, ModVector &r_mods)
{
  // visit each concave square constraint and update Secant if the bounds have
  // changed.
  double lb, ub, lb1, ub1, r;
  ConstraintPtr cons;
  VariablePtr x0, x1, y;
  VarBoundMod2Ptr b2mod;
  VarBoundModPtr bmod;
  LinearFunctionPtr lf;
  LinConModPtr lmod;
  McCormickPtr  mcc;

  for (VarSecantMapIter it=cvCons_.begin(); it != cvCons_.end(); ++it) {
    x0   = it->first;
    y    = it->second->auxVar;
    lb   = x0->getLb();
    ub   = x0->getUb();
    cons = it->second->cons;
    BoundsOnSquare(x0, lb1, ub1);
    if (lb1>y->getLb()+eTol_ && ub1<y->getUb()-eTol_) {
      // bounds on y and also the secant approximation can be updated.
      b2mod  = (VarBoundMod2Ptr) new VarBoundMod2(y, lb1, ub1);
      b2mod->applyToProblem(rel);
      r_mods.push_back(b2mod);
      lf   = getNewSecantLf_(x0, y, lb, ub, r);
      lmod = (LinConModPtr) new LinConMod(cons, lf, -INFINITY, r);
      lmod->applyToProblem(rel);
    } else if (lb1>y->getLb()+eTol_) {
      bmod  = (VarBoundModPtr) new VarBoundMod(y, Lower, lb1);
      bmod->applyToProblem(rel);
      r_mods.push_back(bmod);
      lf   = getNewSecantLf_(x0, y, lb, ub, r);
      lmod = (LinConModPtr) new LinConMod(cons, lf, -INFINITY, r);
      lmod->applyToProblem(rel);
    } else if (ub1<y->getUb()-eTol_) {
      bmod  = (VarBoundModPtr) new VarBoundMod(y, Upper, ub1);
      bmod->applyToProblem(rel);
      r_mods.push_back(bmod);
      lf   = getNewSecantLf_(x0, y, lb, ub, r);
      lmod = (LinConModPtr) new LinConMod(cons, lf, -INFINITY, r);
      lmod->applyToProblem(rel);
    } else {
      lf = cons->getLinearFunction();
      if ((cons->getUb()+lb*ub)>1e-8 || 
           (fabs(lf->getWeight(x0)+(lb+ub)) > 1e-8)) {
        // bounds are up to date but the constraint needs update.
        lf   = getNewSecantLf_(x0, y, lb, ub, r);
        lmod = (LinConModPtr) new LinConMod(cons, lf, -INFINITY, r);
        lmod->applyToProblem(rel);
      }
    }
  }

  for (McCormickSetIter it=mcCons_.begin(); it != mcCons_.end(); ++it) {
    mcc = *it;
    x0 = mcc->getX0();
    x1 = mcc->getX1();
     y = mcc->getAux();
    BoundsOnProduct(x0, x1, lb, ub);
    if (lb>y->getLb()+eTol_ && ub<y->getUb()-eTol_) {
      b2mod  = (VarBoundMod2Ptr) new VarBoundMod2(y, lb, ub);
      r_mods.push_back(b2mod);
      b2mod->applyToProblem(rel);

    } else if (lb>y->getLb()+eTol_) {
      bmod  = (VarBoundModPtr) new VarBoundMod(y, Lower, lb);
      r_mods.push_back(bmod);
      bmod->applyToProblem(rel);

    } else if (ub<y->getUb()-eTol_) {
      bmod  = (VarBoundModPtr) new VarBoundMod(y, Upper, ub);
      r_mods.push_back(bmod);
      bmod->applyToProblem(rel);

    } 

    if (mcc->getSense() == McCormick::LT || mcc->getSense() == McCormick::EQ) {
      lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
          y, r, 0);
      lmod = (LinConModPtr) new LinConMod(mcc->getC0(), lf, -INFINITY, r);
      lmod->applyToProblem(rel);

      lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
          y, r, 1);
      lmod = (LinConModPtr) new LinConMod(mcc->getC1(), lf, -INFINITY, r);
      lmod->applyToProblem(rel);
    }

    if (mcc->getSense() == McCormick::GT || mcc->getSense() == McCormick::EQ) {
      lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
          y, r, 2);
      lmod = (LinConModPtr) new LinConMod(mcc->getC2(), lf, -INFINITY, r);
      lmod->applyToProblem(rel);

      lf = getMcLf_(x0, x0->getLb(), x0->getUb(), x1, x1->getLb(), x1->getUb(), 
          y, r, 3);
      lmod = (LinConModPtr) new LinConMod(mcc->getC3(), lf, -INFINITY, r);
      lmod->applyToProblem(rel);
    }
  }
  return false;
}


SolveStatus CxQuadHandler::presolve(PreModQ *, bool *)
{
  return Finished; // disabled for now.
  removeFixed_();
  binToLin_();
  return Finished;
}


void CxQuadHandler::removeFixed_()
{
  ObjectivePtr oPtr = problem_->getObjective();
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  FunctionPtr f;
  ConstraintConstIterator c_iter;
  double c;

  if (oPtr && oPtr->getFunctionType() == Quadratic) {
    c = 0.;
    f = oPtr->getFunction();
    removeFixedFun_(f, lf, &c);
    if (lf->getNumTerms()>0) {
      lf = (LinearFunctionPtr) new LinearFunction();
    }
    if (fabs(c)>eTol_) {
      problem_->addToObj(c);
    }
  }

  for (c_iter=problem_->consBegin(); c_iter!=problem_->consEnd(); ++c_iter) {
    f = (*c_iter)->getFunction();
    c = 0.;
    if (f->getType()==Quadratic) {
      removeFixedFun_(f, lf, &c);
      if (lf->getNumTerms()>0) {
        lf = (LinearFunctionPtr) new LinearFunction();
      }
      if (fabs(c)>eTol_) {
        problem_->addToCons(*c_iter, c);
      }
    }
  }
}


void CxQuadHandler::removeFixedFun_(FunctionPtr f, LinearFunctionPtr lf2, 
                                    double *c)
{
  QuadraticFunctionPtr qf,qf2;
  LinearFunctionPtr lf;
  bool new_lf = false;
  ConstVariablePtr v1, v2;
  double l1, l2, u1, u2;

  qf = f->getQuadraticFunction();
  lf = f->getLinearFunction();
  if (!lf) {
    new_lf = true;
    lf = lf2;
  }
  qf2 = (QuadraticFunctionPtr) new QuadraticFunction();
  *c  = 0;
  for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
      ++it) {
    v1 = it->first.first;
    v2 = it->first.second;
    l1 = v1->getLb();
    l2 = v2->getLb();
    u1 = v1->getUb();
    u2 = v2->getUb();
    if (u1-l1 < eTol_ && u2-l2 < eTol_) {
      // remove
      *c += l1*l2*it->second;
      qf2->addTerm(v1, v2, -1.*(it->second));
    } else if (u1-l1 < eTol_) {
      // remove and add linear term
      lf->incTerm(v2, l1*(it->second));
      qf2->addTerm(v1, v2, -1.*(it->second));
    } else if (u2-l2 < eTol_) {
      // remove and add linear term
      lf->incTerm(v1, l2*(it->second));
      qf2->addTerm(v1, v2, -1.*(it->second));
    }
  }
  (*qf) += qf2;
  if (new_lf==true && lf->getNumTerms()>0) {
    f->add(lf); // lf is cloned and added.
  }
}


void CxQuadHandler::binToLin_()
{
  if (problem_->getSize()->bins > 0) {
    LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
    FunctionPtr f;
    ConstraintConstIterator c_iter;
    ObjectivePtr oPtr = problem_->getObjective();

    if (oPtr && oPtr->getFunctionType() == Quadratic) {
      f = oPtr->getFunction();
      if (f->getType()==Quadratic) {
        binToLinFun_(f, lf);
        if (lf->getNumTerms()>0) {
          lf = (LinearFunctionPtr) new LinearFunction();
        }
      }
      if (lf->getNumTerms()>0) {
        lf = (LinearFunctionPtr) new LinearFunction();
      }
    }
    for (c_iter=problem_->consBegin(); c_iter!=problem_->consEnd(); ++c_iter) {
      f = (*c_iter)->getFunction();
      if (f->getType()==Quadratic) {
        binToLinFun_(f, lf);
        if (lf->getNumTerms()>0) {
          lf = (LinearFunctionPtr) new LinearFunction();
        }
      }
    }
  }
}


void CxQuadHandler::binToLinFun_(FunctionPtr f, LinearFunctionPtr lf2)
{
  QuadraticFunctionPtr qf;
  LinearFunctionPtr lf;
  bool new_lf = false;

  qf = f->getQuadraticFunction();
  lf = f->getLinearFunction();
  if (!lf) {
    new_lf = true;
    lf = lf2;
  }
  for (VariablePairGroupConstIterator it=qf->begin(); it!=qf->end(); 
      ++it) {
    if (it->first.first == it->first.second && it->first.first->getType() == 
        Binary) {
      lf->incTerm(it->first.first, 1*(it->second));
      qf->incTerm(it->first.first, it->first.first, -1.*(it->second));
      assert(!"bug above!");
    }
  }
  if (new_lf==true && lf->getNumTerms()>0) {
    f->add(lf); // lf is cloned and added.
  }
}


std::string CxQuadHandler::getName() const
{
   return "CxQuadHandler (Handling quadratic terms).";
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

McCormick::McCormick(VariablePtr x0, VariablePtr x1, McCormick::Sense sense)
  : s_(sense)
{
  if (x0->getIndex()>x1->getIndex()) {
    x0_ = x1;
    x1_ = x0;
  } else {
    x0_ = x0;
    x1_ = x1;
  }
  y_ = VariablePtr(); // NULL
  c0_ = c1_ = c2_ = c3_ = ConstraintPtr(); // NULL
}


McCormick::~McCormick() 
{
  x0_.reset();
  x1_.reset();
  y_.reset();
  c0_.reset();
  c1_.reset();
  c2_.reset();
  c3_.reset();
}


bool Minotaur::CompareMcCormick::operator()(McCormickPtr b0, McCormickPtr b1)
  const
{
  UInt b0x0 = b0->getX0()->getId();
  UInt b0x1 = b0->getX1()->getId();

  UInt b1x0 = b1->getX0()->getId();
  UInt b1x1 = b1->getX1()->getId();

  if (b0x0 == b1x0) {
    return (b0x1 < b1x1);
  }
  return (b0x0 < b1x0);
}


bool McCormick::isViolated(const double *x, const double &tol) const
{

  double xval = x[x0_->getIndex()] * x[x1_->getIndex()];
  double yval  = x[y_->getIndex()];
  switch (s_) {
   case (LT): 
     if (xval > yval+tol) {
       return true;
     }
     break;
   case (GT):
     if (xval < yval-tol) {
       return true;
     }
     break;
   case (EQ):
     if (fabs(xval - yval) > tol) {
       return true;
     }
     break;
   default:
     break;
  }

  return false;
}


bool McCormick::isViolated(const double &x0val, const double &x1val, 
                           const double &yval, const double &tol) const
{
  double xval = x1val*x0val;
  switch (s_) {
   case (LT): 
     if (xval > yval+tol) {
       return true;
     }
     break;
   case (GT):
     if (xval < yval-tol) {
       return true;
     }
     break;
   case (EQ):
     if (fabs(xval - yval) > tol) {
       return true;
     }
     break;
   default:
     break;
  }
  return false;
}


VariablePtr McCormick::getOtherX(ConstVariablePtr x) const
{
  if (x0_==x) {
   return x1_;
  } else if (x1_==x) {
    return x0_;
  } else {
    VariablePtr v = VariablePtr(); // NULL
    return v;
  }
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
