//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MaxVioBrancher.cpp
 * \brief Define methods for maximum violation branching candidate.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iomanip>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrCand.h"
#include "BrVarCand.h"
#include "Environment.h"
#include "Handler.h"
#include "Logger.h"
#include "MaxVioBrancher.h"
#include "Modification.h"
#include "Node.h"
#include "Option.h"
#include "Relaxation.h"
#include "Solution.h"
#include "Timer.h"

// #define SPEW 1

using namespace Minotaur;

const std::string MaxVioBrancher::me_ = "MaxVio brancher: "; 

MaxVioBrancher::MaxVioBrancher(EnvPtr env, HandlerVector & handlers) 
: handlers_(handlers),    // Create a copy, the vector is not too big
  rel_(RelaxationPtr())   // NULL
  
{
  timer_ = env->getNewTimer();
  logger_ = (LoggerPtr) new Logger((LogLevel) env->getOptions()->
                                   findInt("br_log_level")->getValue());
  stats_ = new MaxVioBrStats();
  stats_->calls = 0;
  stats_->time = 0.0;
}


MaxVioBrancher::~MaxVioBrancher()
{
  delete stats_;
  delete timer_;
  logger_.reset();
}


Branches MaxVioBrancher::findBranches(RelaxationPtr rel, NodePtr , 
                                      ConstSolutionPtr sol,
                                      SolutionPoolPtr s_pool, 
                                      BrancherStatus & br_status,
                                      ModVector &mods) 
{
  Branches branches;
  BrCandPtr br_can = BrCandPtr(); //NULL
  const double *x = sol->getPrimal();
  bool should_prune = false;

  timer_->start();
  x_.resize(rel->getNumVars());
  std::copy(x, x+rel->getNumVars(), x_.begin());

  //sol->write(std::cout);
  
  ++(stats_->calls);
  rel_ = rel;
  br_status = NotModifiedByBrancher;

  findCandidates_(mods, should_prune);
  if (true==should_prune) {
    br_status = PrunedByBrancher;
  } else {
    br_can = findBestCandidate_();
  }

  if (mods.size()>0) {
    br_status = ModifiedByBrancher;
  } else if (br_can) {
    branches = br_can->getHandler()->getBranches(br_can, x_, rel_, s_pool); 
    for (BranchConstIterator br_iter=branches->begin(); 
         br_iter!=branches->end(); ++br_iter) {
      (*br_iter)->setBrCand(br_can);
    }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "best candidate = "
                               << br_can->getName() 
                               << " Handler: "
                               << br_can->getHandler()->getName() << std::endl;
#endif
  } else if (br_status != PrunedByBrancher) {
    assert(!"problem finding candidate in MaxVioBrancher");
  }
  stats_->time += timer_->query();
  timer_->stop();
  return branches;
}


void MaxVioBrancher::findCandidates_(ModVector &mods, bool &should_prune)
{
  BrVarCandSet cands2;    // Temporary set.
  BrCandVector gencands2; // Temporary set.
  BrCandPtr br_can;
  std::pair<BrVarCandIter, bool> ret;

  should_prune = false;
  cands_.clear();
  gencands_.clear();
  for (HandlerIterator h = handlers_.begin(); h != handlers_.end(); ++h) {
    // ask each handler to give some candidates
    (*h)->getBranchingCandidates(rel_, x_, mods, cands2, gencands2, should_prune);
    if (should_prune) {
      break;
    }
    for (BrVarCandIter it = cands2.begin(); it != cands2.end(); ++it) {
      (*it)->setHandler(*h);
      ret = cands_.insert(*it);
      if (false == ret.second) { // already exists.
        br_can = *(ret.first);
        br_can->setDist((*it)->getDDist()+br_can->getDDist(),
                        (*it)->getDDist()+br_can->getUDist());
      }
    }
    for (BrCandVIter it = gencands2.begin(); it != gencands2.end(); ++it) {
      (*it)->setHandler(*h);
    }
    gencands_.insert(gencands_.end(), gencands2.begin(), gencands2.end());
    cands2.clear();
    gencands2.clear();
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "candidates: " << std::endl;
  for (BrVarCandIter it=cands_.begin(); it!=cands_.end(); ++it) {
    logger_->msgStream(LogDebug)
        << std::setprecision(6) << (*it)->getName() << "\t" 
        << (*it)->getDDist() << "\t" << (*it)->getUDist()
        << std::endl;
  }
#endif
  return;
}


BrCandPtr MaxVioBrancher::findBestCandidate_()
{
  BrCandPtr best_cand = BrCandPtr(); // NULL
  double best_score = -INFINITY;
  double cand_score;
  double lmin = 0.8;
  double lmax = 0.2;

  for (BrVarCandIter it = cands_.begin(); it != cands_.end(); ++it) {
    cand_score = lmin*std::min((*it)->getDDist(),(*it)->getUDist()) + 
                 lmax*std::max((*it)->getDDist(),(*it)->getUDist());
    if (VarOrig==(*it)->getVar()->getSrcType()) {
      cand_score = 0.1*cand_score;
    } 
    if (cand_score > best_score) {
      best_score = cand_score;
      best_cand = (*it);
    }
  }

  for (BrCandVIter it = gencands_.begin(); it != gencands_.end(); ++it) {
    cand_score = lmin*std::min((*it)->getDDist(),(*it)->getUDist()) + 
                 lmax*std::max((*it)->getDDist(),(*it)->getUDist());
    if (cand_score > best_score) {
      best_score = cand_score;
      best_cand = (*it);
    }
  }

  if (best_cand) {
    if (best_cand->getDDist() > best_cand->getUDist()) {
      best_cand->setDir(UpBranch);
    } else {
      best_cand->setDir(DownBranch);
    }
  }

  return best_cand;
}


void MaxVioBrancher::updateAfterLP(NodePtr , ConstSolutionPtr )
{
}


void MaxVioBrancher::writeStats(std::ostream &out) const
{
  if (stats_) {
    out << me_ << "times called = " << stats_->calls << std::endl
      << me_ << "time in brancher = " << stats_->time << std::endl;
  }
}


std::string MaxVioBrancher::getName() const
{
  return "MaxVioBrancher";
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
