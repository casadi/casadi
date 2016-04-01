//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MaxFreqBrancher.cpp
 * \brief Define methods for maximum-frequency branching.
 * \author Suresh B, IIT Bombay
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
#include "MaxFreqBrancher.h"
#include "Modification.h"
#include "Node.h"
#include "Option.h"
#include "Relaxation.h"
#include "Solution.h"
#include "Timer.h"
#include "Variable.h"

//#define SPEW 1

using namespace Minotaur;

const std::string MaxFreqBrancher::me_ = "MaxFreq Brancher: "; 

MaxFreqBrancher::MaxFreqBrancher(EnvPtr env, HandlerVector & handlers) 
: status_(NotModifiedByBrancher),
  rel_(RelaxationPtr()),  // NULL
  handlers_(handlers),     // Create a copy, the vector is not too big
  init_(false)
{
  timer_ = env->getNewTimer();
  logger_ = (LoggerPtr) new Logger((LogLevel) env->getOptions()->
                                   findInt("br_log_level")->getValue());
  zTol_ = 1e-6;
  stats_ = new MaxFreqBrStats();
  stats_->calls = 0;
  stats_->time = 0.0;
}


MaxFreqBrancher::~MaxFreqBrancher()
{
  delete stats_;
  delete timer_;
  logger_.reset();
}


Branches MaxFreqBrancher::findBranches(RelaxationPtr rel, NodePtr , 
                                       ConstSolutionPtr sol,
                                       SolutionPoolPtr s_pool, 
                                       BrancherStatus & br_status,
                                       ModVector &) 
{
  Branches branches;
  BrCandPtr br_can = BrCandPtr(); //NULL
  const double *x = sol->getPrimal();

  timer_->start();
  x_.resize(rel->getNumVars());
  std::copy(x, x+rel->getNumVars(), x_.begin());

  ++(stats_->calls);
  if (!init_) {
    init_ = true;
    initialize(rel);
  }
  rel_ = rel;
  br_status = NotModifiedByBrancher;

  findCandidates_();

  if (status_ == PrunedByBrancher) {
    br_status = status_;
    return branches;
  }

  updateFracCount_();
  updateUnfixedCount_();

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "candidates: " << std::endl;
  for (BrVarCandIter it=cands_.begin(); it!=cands_.end(); ++it) {
    logger_->msgStream(LogDebug)
        << std::setprecision(6) << (*it)->getName() << "\t" 
        << x_[(*it)->getPCostIndex()] << "\t"
        << fracCount_[(*it)->getPCostIndex()] << "\t"
        << unfixedCount_[(*it)->getPCostIndex()] << "\t"
        << std::endl;
  }
#endif

  br_can = findBestCandidate_();
 
  if (br_can) {
    branches = br_can->getHandler()->getBranches(br_can, x_, rel_, s_pool); 
    for (BranchConstIterator br_iter=branches->begin(); 
         br_iter!=branches->end(); ++br_iter) {
      (*br_iter)->setBrCand(br_can);
    }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "best candidate = "
                               << br_can->getName() 
                               << " Handler: "
                               << br_can->getHandler()->getName()
                               << std::endl;
#endif
  } else {
    assert(!"problem finding candidate in MaxFreqBrancher");
  }
  stats_->time += timer_->query();
  timer_->stop();
  return branches;
}


void MaxFreqBrancher::initialize(RelaxationPtr rel)
{
  int n = rel->getNumVars();
  // initialize following vectors to zero.
  fracCount_ = std::vector<UInt>(n,0); 
  unfixedCount_ = std::vector<UInt>(n,0); 
}


void MaxFreqBrancher::findCandidates_()
{
  BrVarCandSet cands2;      // Temporary set.
  BrCandVector gencands;
  ModVector mods;        // handlers may ask to modify the problem.
  bool is_inf = false;
  std::string handlerName;

  cands_.clear();
  for (HandlerIterator h = handlers_.begin(); h != handlers_.end(); ++h) {
    // ask each handler to give some candidates
    (*h)->getBranchingCandidates(rel_, x_, mods, cands2, gencands, is_inf);
    for (BrVarCandIter it = cands2.begin(); it != cands2.end(); ++it) {
      (*it)->setHandler(*h);
    }
    cands_.insert(cands2.begin(), cands2.end());
    if (is_inf) {
      cands2.clear();
      cands_.clear();
      status_ = PrunedByBrancher;
      return;
    }
    cands2.clear();
  }

  return;
}


BrCandPtr MaxFreqBrancher::findBestCandidate_()
{
  BrCandPtr best_cand = BrCandPtr(); // NULL
  double best_score = -1;
  double cand_score;
  UInt index;

#if SPEW
 logger_->msgStream(LogDebug) << me_ << "candidate score from BestCand Fn: "
                                     << std::endl;
#endif
 for (BrVarCandIter it = cands_.begin(); it != cands_.end(); ++it) {
    index = (*it)->getPCostIndex();
    cand_score = (double) fracCount_[index]/((double) unfixedCount_[index]);
#if SPEW
   logger_->msgStream(LogDebug)
       << std::setprecision(6) << (*it)->getName() << "\t"
       << cand_score << "\t"
       << std::endl;
#endif

   if (cand_score > best_score) {
      best_score = cand_score;
      best_cand = (*it);
    }
  }

  if (best_cand) {
      best_cand->setDir(DownBranch);
  }
  return best_cand;
}


void MaxFreqBrancher::updateFracCount_()
{
  UInt index;
  for (BrVarCandIter it = cands_.begin(); it != cands_.end(); ++it) {
    index = (*it)->getPCostIndex();
    fracCount_[index] += 1;
  }
}


void MaxFreqBrancher::updateUnfixedCount_()
{
  VariablePtr v;
  VariableType v_type;
  UInt index;

  for (VariableConstIterator it = rel_->varsBegin(); it != rel_->varsEnd();
       ++it) {
    v = (*it);
    v_type = v->getType();
    index = v->getIndex();
    if ((v_type == Binary || v_type == Integer) &&
        (fabs(v->getLb() - v->getUb()) > zTol_)) {
      unfixedCount_[index] += 1;
    }
  }
}

void MaxFreqBrancher::updateAfterLP(NodePtr , ConstSolutionPtr )
{
}


void MaxFreqBrancher::writeStats(std::ostream &out) const
{
  if (stats_) {
    out << me_ << "times called = " << stats_->calls << std::endl
      << me_ << "time in brancher = " << stats_->time << std::endl;
  }
}


std::string MaxFreqBrancher::getName() const
{
  return "MaxFreqBrancher";
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
