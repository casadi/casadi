//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file RandomBrancher.cpp
 * \brief Define methods for Random branching.
 * \author Suresh B, IIT Bombay
 */

#include <iostream>

#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrCand.h"
#include "BrVarCand.h"
#include "Environment.h"
#include "Handler.h"
#include "RandomBrancher.h"
#include "Logger.h"
#include "Option.h"
#include "Problem.h"
#include "Relaxation.h"
#include "Solution.h"
#include "Timer.h"

//#define SPEW 1

using namespace Minotaur;

const std::string RandomBrancher::me_ = "RandomBrancher: "; 

RandomBrancher::RandomBrancher()
  : handlers_(0),
    timer_(0)
{
  logger_ = new Logger(LogInfo);
  stats_ = new RandomBrStats();
  stats_->calls = 0;
  stats_->time = 0.0;
}


RandomBrancher::RandomBrancher(EnvPtr env, HandlerVector handlers)
{
  logger_ = new Logger((LogLevel) 
                       env->getOptions()->findInt("br_log_level")->getValue());
  timer_ = env->getNewTimer();
  stats_ = new RandomBrStats();
  stats_->calls = 0;
  stats_->time = 0.0;
  handlers_ = handlers;
  seed_ = env->getOptions()->findInt("rand_seed")->getValue();
  if (seed_ == 0) {
    srand(time(NULL));
  } else {
    srand(seed_);
  }
}


RandomBrancher::~RandomBrancher()
{
  if (logger_) {
    delete logger_;
  }
  if (timer_) {
    delete timer_;
  }
  if (stats_) {
    delete stats_;
  }
  handlers_.clear();
}


Branches RandomBrancher::findBranches(RelaxationPtr rel, NodePtr ,
                                      ConstSolutionPtr sol,
                                      SolutionPoolPtr s_pool,
                                      BrancherStatus & br_status,
                                      ModVector &mods) 
{
  Branches branches;
  DoubleVector x(rel->getNumVars());
  BrVarCandSet cands;      // candidates from which to choose one.
  BrVarCandSet cands2;     // temporary set.
  BrCandVector gencands;
  BrCandVector gencands2;  // temporary set.
  BrCandPtr best_can = BrCandPtr(); // NULL
  bool is_inf = false;

  timer_->start();
  std::copy(sol->getPrimal(), sol->getPrimal()+rel->getNumVars(), x.begin());
  
  ++(stats_->calls);
  br_status = NotModifiedByBrancher;

  for (HandlerIterator h = handlers_.begin(); h != handlers_.end(); ++h) {
    // ask each handler to give some candidates
    (*h)->getBranchingCandidates(rel, x, mods, cands2, gencands2, is_inf);
    for (BrVarCandIter it = cands2.begin(); it != cands2.end(); ++it) {
      (*it)->setHandler(*h);
    }
    for (BrCandVIter it = gencands2.begin(); it != gencands2.end(); ++it) {
      (*it)->setHandler(*h);
    }
    cands.insert(cands2.begin(), cands2.end());
    gencands.insert(gencands.end(), gencands2.begin(), gencands2.end());
    cands2.clear();
    gencands2.clear();
    if (is_inf || mods.size()>0) {
      cands.clear();
      gencands.clear();
      if (is_inf) {
        br_status = PrunedByBrancher;
        stats_->time += timer_->query();
        timer_->stop();
        return branches;
      }
    }
  }

  if (cands.size() > 0) {
    BrVarCandIter it = cands.begin();

    std::advance(it,rand()%cands.size());
    best_can = *(it);
  } else if (gencands.size() > 0) {
    BrCandVIter it = gencands.begin();

    std::advance(it,rand()%gencands.size());
    best_can = *(it);
  }

  if (best_can) {
    best_can->setDir(DownBranch);

    branches = best_can->getHandler()->getBranches(best_can, x, rel, s_pool); 
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "best candidate = "
                                 << best_can->getName() << std::endl;
#endif
  } else {
    assert(!"problem finding candidate in RandomBrancher");
  }

  stats_->time += timer_->query();
  timer_->stop();
  return branches;
}


void RandomBrancher::writeStats(std::ostream &out) const 
{
  if (stats_) {
    out << me_ << "times called = " << stats_->calls << std::endl
      << me_ << "time in brancher = " << stats_->time << std::endl;
  }
}


std::string RandomBrancher::getName() const
{
  return "RandomBrancher";
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
