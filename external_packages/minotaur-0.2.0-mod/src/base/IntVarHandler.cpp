//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file IntVarHandler.cpp
 * \brief Define the IntVarHandler class for handling integer constrained
 * variables. It checks integrality and provides branching candidates. Does
 * not do any presolving and cut-generation.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>

#include "MinotaurConfig.h"
#include "BrVarCand.h"
#include "Branch.h"
#include "Environment.h"
#include "IntVarHandler.h"
#include "Logger.h"
#include "Option.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "VarBoundMod.h"
#include "Variable.h"

//#define SPEW 1

using namespace Minotaur;
const std::string IntVarHandler::me_ = "IntVarHandler: ";

IntVarHandler::IntVarHandler(EnvPtr env, ProblemPtr problem)
  : env_(env)
{
  logger_   = (LoggerPtr) new Logger((LogLevel) env_->getOptions()->
                                     findInt("handler_log_level")->getValue());
  modProb_  = true;
  modRel_   = true;
  intTol_   = env_->getOptions()->findDouble("int_tol")->getValue();
  gDive_    = env_->getOptions()->findBool("guided_dive")->getValue();
  problem_  = problem;
}


IntVarHandler::~IntVarHandler()
{
  problem_.reset();
  env_.reset();
  logger_.reset();
}


bool IntVarHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                               bool &, double &inf_meas)
{
  VariableConstIterator v_iter, v_iter2;
  VariableType v_type;
  double value;
  const double *x = sol->getPrimal();
  bool is_feas = true;

  inf_meas = 0.0;
  for (v_iter=relaxation->varsBegin(); v_iter!=relaxation->varsEnd(); 
       ++v_iter) {
    v_type = (*v_iter)->getType();
    if (v_type==Binary || v_type==Integer) {
      value = x[(*v_iter)->getIndex()];
      if (fabs(value - floor(value+0.5)) > intTol_) {
        is_feas = false;
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "variable " <<
          (*v_iter)->getName() << " has fractional value = " << value <<
          std::endl;
#endif
        inf_meas += fabs(value-floor(value+0.5));
      }
    }
  }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "is_feas = " << is_feas << std::endl;
#endif
  return is_feas;
}


void IntVarHandler::getBranchingCandidates(RelaxationPtr rel, 
                                           const DoubleVector &x,
                                           ModVector &, BrVarCandSet &cands,
                                           BrCandVector &, bool &is_inf)
{
  VariablePtr v;
  VariableType v_type;
  UInt index;
  BrVarCandPtr br_can;

  for (VariableConstIterator it=rel->varsBegin(); it!=rel->varsEnd(); ++it) {
    v = *it;
    v_type = v->getType();
    index = v->getIndex();
    if ((v_type==Binary || v_type==Integer) && 
        fabs(floor(x[index]+0.5) - x[index]) > intTol_) {
      // yes, it can be branched upon.
      br_can = (BrVarCandPtr) new BrVarCand(v, v->getIndex(), 
                                            x[index]-floor(x[index]),
                                            ceil(x[index])-x[index]);
      cands.insert(br_can);
    } 
  }
  is_inf = false;
}


ModificationPtr IntVarHandler::getBrMod(BrCandPtr cand, DoubleVector & x,
                                        RelaxationPtr , BranchDirection dir) 
{
  // TODO: fix this dynamic cast
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  VarBoundModPtr mod;
  double bnd;

  if (dir==DownBranch) {
    bnd = floor(x[v->getIndex()]);
    mod = (VarBoundModPtr) new VarBoundMod(v, Upper, bnd);
  } else {
    bnd = ceil(x[v->getIndex()]);
    mod = (VarBoundModPtr) new VarBoundMod(v, Lower, bnd);
  }
  return mod;
}


Branches IntVarHandler::getBranches(BrCandPtr cand, DoubleVector & x, 
                                    RelaxationPtr rel, SolutionPoolPtr s_pool)
{
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  VariablePtr v2;
  double value = x[v->getIndex()];
  BranchPtr branch1, branch2;
  Branches branches = (Branches) new BranchPtrVector();
  VarBoundModPtr mod;
  SolutionPtr bestsol = s_pool->getBestSolution();

  branch1 = (BranchPtr) new Branch();
  if (modProb_) {
    v2 = rel->getOriginalVar(v);
    mod = (VarBoundModPtr) new VarBoundMod(v2, Upper, floor(value));
    branch1->addPMod(mod);
  }
  if (modRel_) {
    mod = (VarBoundModPtr) new VarBoundMod(v, Upper, floor(value));
    branch1->addRMod(mod);
  }
  branch1->setActivity(value);

  branch2 = (BranchPtr) new Branch();
  if (modProb_) {
    v2 = rel->getOriginalVar(v);
    if (v2) {
      mod = (VarBoundModPtr) new VarBoundMod(v2, Lower, ceil(value));
      branch2->addPMod(mod);
    }
  }
  if (modRel_) {
    mod = (VarBoundModPtr) new VarBoundMod(v, Lower, ceil(value));
    branch2->addRMod(mod);
  }
  branch2->setActivity(value);

  if (true==gDive_ && bestsol) {
    if (bestsol->getPrimal()[v->getIndex()] < x[v->getIndex()]) {
      branches->push_back(branch1);
      branches->push_back(branch2);
    } else {
      branches->push_back(branch2);
      branches->push_back(branch1);
    }
  } else {
    if (cand->getDir() == DownBranch) {
      branches->push_back(branch1);
      branches->push_back(branch2);
    } else {
      branches->push_back(branch2);
      branches->push_back(branch1);
    }
  }
  return branches;
}


double IntVarHandler::getTol() const
{
  return intTol_;
}


bool IntVarHandler::isNeeded()
{
  if (problem_) {
    problem_->calculateSize();
    if (problem_->getSize()->bins > 0 || problem_->getSize()->ints > 0) {
      return true;
    }
  }
  return false;
}


void IntVarHandler::relaxInitFull(RelaxationPtr, bool *is_inf)
{
  *is_inf = false;
}


void IntVarHandler::relaxInitInc(RelaxationPtr , bool *is_inf)
{
  *is_inf = false;
}


void IntVarHandler::relaxNodeFull(NodePtr , RelaxationPtr, bool *is_inf)
{
  *is_inf = false;
}


void IntVarHandler::relaxNodeInc(NodePtr , RelaxationPtr , bool *is_inf)
{
  *is_inf = false;
}


void IntVarHandler::setTol(double tol)
{
  intTol_ = tol;
}


std::string IntVarHandler::getName() const
{
  return "IntVarHandler (Handling integrality of variables).";
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
