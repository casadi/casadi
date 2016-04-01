// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NLPMultiStart.cpp
 * \brief Multi start heuristic for continuous NLPs
 * \author Jayash Koshal, Argonne National Laboratory
 */

#include <cmath> // for INFINITY

#include "MinotaurConfig.h"
#include "Engine.h"
#include "Variable.h"
#include "Environment.h"
#include "NLPMultiStart.h"
#include "Logger.h"
#include "Node.h"
#include "Operations.h"
#include "Option.h"
#include "SolutionPool.h"
#include "Timer.h"
#include <iomanip>
#include <algorithm>
#include <cstdlib>
using namespace Minotaur;

//#define SPEW 1

const std::string NLPMultiStart::me_ = "NLP Multi-Start Heuristic: ";

NLPMultiStart::NLPMultiStart(EnvPtr env, ProblemPtr p, EnginePtr e)
: e_(e),
  env_(env),
  p_(p)
{
  VariablePtr variable;
  UInt n      = p_->getNumVars();
  distBound_  = 0.0;
  for (UInt i=0; i<n; ++i) {
    variable = p_->getVariable(i);
    distBound_ += pow(variable->getUb()-variable->getLb(),2);
  }
  distBound_ = (distBound_ >= INFINITY) ? 10.0*sqrt(n) : sqrt(distBound_);
  logger_ = (LoggerPtr) new Logger((LogLevel) env_->getOptions()->
                                   findInt("heur_log_level")->getValue());
  random_                   = new double[n];  

  // statistics
  stats_.numNLPs           = 0;
  stats_.numInfeas         = 0;
  stats_.numImprove        = 0;
  stats_.numBadstatus      = 0;
  stats_.time              = 0;
  stats_.iterations        = 0;
  stats_.bestObjValue      = INFINITY;

}


NLPMultiStart::~NLPMultiStart(){
  delete [] random_;
}


void NLPMultiStart::constructInitial_(double* a, const double* b, double rho,
                                      UInt n)
{
  double dist;
  VariablePtr variable;
  double norm;
  VariableConstIterator v_iter;

  for (UInt i=0; i<n; ++i) {
    random_[i] = rand()/double(RAND_MAX) - 0.5;
  }
  norm = sqrt(InnerProduct(random_, random_, n)); 
  for (UInt i=0; i<n; ++i) {
    random_[i] /= norm;
  }

#if SPEW
    logger_->msgStream(LogDebug) << me_ << "rho = " << rho << std::endl;
#endif

  // find how far we need to go.
  dist  = getDistance(a, b, n);
  if (dist <= distBound_) {
    dist = distBound_;
  }
  dist *= rho;

  for (v_iter=p_->varsBegin(); v_iter!=p_->varsEnd(); ++v_iter, 
      ++a, ++b, ++random_) {
    variable = *v_iter;
    // find a point in a random direction outside the ball
    // centered around x* with radius = ||x*-x||
    *a = std::max(std::min(*b + (*random_) * dist, 
          variable->getUb()), variable->getLb());
  }
  random_ -= n;

#if SPEW
  logger_->msgStream(LogDebug2)
    << me_ << "distance to new point = " 
    <<  getDistance(a-n, b-n, n)
    << std::endl;
#endif 
}


void NLPMultiStart::solve(NodePtr, RelaxationPtr, SolutionPoolPtr s_pool)
{
  ConstSolutionPtr sol; 
  EngineStatus status;
  UInt heur_bound                = 10; // no. of times the heuristic should run
  UInt unchanged_obj_count_limit =  3;
  double obj_tol                 = 1e-6; 
  double rho_initial             = 1.1;// amplification factor
  double rho                     = rho_initial;
  Timer *timer                   = env_->getNewTimer();
  UInt n                         = p_->getNumVars();
  double* prev_feasible          = new double[n];
  double* initial_point          = new double[n];

  // start at a random point.
  for (UInt i=0; i<n; ++i){
    initial_point[i] = rand()/double (RAND_MAX);
  }

  srand(1);
  timer->start();
  e_->clear();
  e_->load(p_);
  for (UInt i=0, unchanged_obj_count=0; i < heur_bound &&
       unchanged_obj_count < unchanged_obj_count_limit; ++i) {
    // XXX: ashu to bring this out of the loop.
    p_->setInitialPoint(initial_point);
    status = e_->solve();
    ++(stats_.numNLPs);
    sol = e_->getSolution();
    if (sol->getObjValue() < stats_.bestObjValue - obj_tol) {
      stats_.bestObjValue = sol->getObjValue();
      rho = rho_initial;
      unchanged_obj_count = 0;
      s_pool->addSolution(sol);
      constructInitial_(initial_point, sol->getPrimal(), rho, n);
      std::copy(sol->getPrimal(), sol->getPrimal() + n, prev_feasible); 
      ++(stats_.numImprove);
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "Better solution " 
        << stats_.bestObjValue << std::endl;
#endif
    } else if ((ProvenInfeasible==status || ProvenLocalInfeasible==status ||
        ProvenObjectiveCutOff==status || ProvenFailedCQInfeas==status || 
        FailedInfeas==status) || ((FailedFeas==status || 
        ProvenFailedCQFeas==status || ProvenLocalOptimal==status) && 
        stats_.bestObjValue <= sol->getObjValue() + obj_tol)) {
      rho *= 1.07;
      ++unchanged_obj_count;
      // use previously found feasible solution
      constructInitial_(initial_point, prev_feasible, rho, n);
      ++(stats_.numInfeas);
#if SPEW
      logger_->msgStream(LogDebug) << me_ 
        << "Engine status = " << status << std::endl
        << me_
        << "Optimal solution no better than best known " << sol->getObjValue()  
        << std::endl;
#endif
    } else if (status == ProvenUnbounded) { 
      rho *= 0.9;
      ++unchanged_obj_count;
      ++(stats_.numBadstatus);
      // use previously found feasible solution
      constructInitial_(initial_point, prev_feasible, rho, n);       
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "Unbounded." << std::endl;
#endif
    } else { 
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "Solution found is not optimal" 
        << " solution value = " <<  sol->getObjValue() << std::endl; 
#endif
    }
    stats_.time = timer->query();
  }
  if (timer) {
    delete timer;
  }
  if (initial_point) {
    delete [] initial_point;
  }
  if (prev_feasible) {
    delete [] prev_feasible;
  }
}


void NLPMultiStart::writeStats(std::ostream &out) const
{
  out << me_ << " number of nlps solved                 = " 
    << stats_.numNLPs << std::endl
    << me_ << " number of Infeasible solves           = " 
    << stats_.numInfeas << std::endl
    << me_ << " number of Improvements in objective   = " 
    << stats_.numImprove << std::endl
    << me_ << " number of Bad status(unbounded etc)   = " 
    << stats_.numNLPs << std::endl
    << me_ << " total time taken                      = " 
    << stats_.time << std::endl
    << me_ << " number of iterations                  = " 
    << stats_.iterations << std::endl
    << me_ << " Best Objective Value                  = " 
    << stats_.bestObjValue << std::endl;
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
