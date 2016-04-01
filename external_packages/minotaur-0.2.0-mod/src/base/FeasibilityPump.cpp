//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
 
/**
 * \file FeasibilityPump.cpp
 * \brief Define the Feasibility Pump class for generating a feasible solution 
 * using FP heuristic for MINLPs.
 * \author Jayash Koshal, Argonne National Laboratory
 *
 * Implements the class Feasibility Pump.
 */

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "FeasibilityPump.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Node.h"
#include "Operations.h"
#include "Option.h"
#include "SolutionPool.h"
#include "Objective.h"
#include "ProblemSize.h"
#include "Timer.h"
#include "Types.h"
#include "Variable.h"
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace Minotaur;

//#define SPEW 1

const std::string FeasibilityPump::me_ = "Feasibility Pump: ";

FeasibilityPump::FeasibilityPump(EnvPtr env, ProblemPtr p, EnginePtr e)
: e_(e),
  env_(env),
  intTol_(1e-6),
  nToFlip_(2),
  p_(p),
  stats_(NULL)
{
  // initialize the random vector for hashing
  srand(1);
  VariablePtr variable; 
  
  for (VariableConstIterator v_iter=p_->varsBegin();
        v_iter!=p_->varsEnd(); ++v_iter) {
    variable = *v_iter;
    if (variable->getType() == Binary) {
      bins_.push_back(variable);
      random_.push_back(((double) rand()/(RAND_MAX)));
    } else {
      random_.push_back(0.0);
    }
  }

  roundedSol_.reserve(p_->getNumVars());
  std::fill(roundedSol_.begin(), roundedSol_.end(), 0);
    
  // initialize the logger pointer
  logger_ = (LoggerPtr) new Logger((LogLevel) env_->getOptions()->
     findInt("heur_log_level")->getValue());
  
  // statistics for Feasibilty Pump heuristic
  stats_                    = new FeasPumpStats();
  stats_->numNLPs           = 0;
  stats_->errors            = 0;
  stats_->numCycles         = 0;
  stats_->time              = 0;
  stats_->bestObjValue      = INFINITY;

  timer_                    = env_->getNewTimer();
}


FeasibilityPump::FeasibilityPump(EnvPtr env, ProblemPtr p, EnginePtr nlpe,
    EnginePtr)
: e_(nlpe),
  env_(env),
  intTol_(1e-6),
  nToFlip_(2),
  p_(p),
  stats_(NULL)
{
  // initialize the random vector for hashing
  srand(1);
  VariablePtr variable; 
  
  for (VariableConstIterator v_iter=p_->varsBegin();
        v_iter!=p_->varsEnd(); ++v_iter) {
    variable = *v_iter;
    if (variable->getType() == Binary) {
      bins_.push_back(variable);
      random_.push_back(((double) rand()/(RAND_MAX)));
    } else {
      random_.push_back(0.0);
    }
  }

  roundedSol_.reserve(p_->getNumVars());
  std::fill(roundedSol_.begin(), roundedSol_.end(), 0);
 
  // statistics for Feasibilty Pump heuristic
  stats_                    = new FeasPumpStats();
  stats_->numNLPs           = 0;
  stats_->errors            = 0;
  stats_->numCycles         = 0;
  stats_->time              = 0;
  stats_->bestObjValue      = INFINITY;

  timer_                    = env_->getNewTimer();
 }


FeasibilityPump::~FeasibilityPump()
{
  delete stats_;

  if (timer_) {
    delete timer_;
  }
}


void FeasibilityPump::constructObj_(ProblemPtr prob, ConstSolutionPtr)
{
  double value, lb, ub;
  VariablePtr variable;

  UInt i                  = 0;
  double constant         = 0;
  LinearFunctionPtr lf    = (LinearFunctionPtr) new LinearFunction();
  FunctionPtr f           = (FunctionPtr) new Function(lf);

  for (VariableConstIterator v_iter=prob->varsBegin();
      v_iter!=prob->varsEnd(); ++v_iter, ++i) {
    variable = *v_iter;
    if (variable->getType() == Binary) {
      value  = roundedSol_[i];
      lb     = variable->getLb();
      ub     = variable->getUb();
      if (fabs(value - lb) > intTol_) { 
        lf->addTerm(variable, -1.0);
        constant += ub;
#if SPEW
        logger_->msgStream(LogDebug2) << me_ << "Including variable for UB"
          << std::endl;
        variable->write(logger_->msgStream(LogDebug2));
#endif
      } else if (fabs(value - ub) > intTol_) {
        lf->addTerm(variable, 1.0);
        constant -= lb;
#if SPEW
        logger_->msgStream(LogDebug2) << me_ << "Including variable for LB"
          << std::endl;
        variable->write(logger_->msgStream(LogDebug2));
#endif
      } else {
        // add a new variable with coeff 1
        // add two constraints for absolute value
        lf->addTerm(variable, 0);
#if SPEW
        logger_->msgStream(LogDebug2) << me_ << "Including absolute value"
          << std::endl;
        variable->write(logger_->msgStream(LogDebug2));
#endif
      }      
    }
  }
  prob->changeObj(f, constant);
}


void FeasibilityPump::convertSol_(SolutionPoolPtr s_pool, ConstSolutionPtr sol)
{
  UInt i            = 0;
  UInt numvars      = p_->getNumVars();
  const double* x   = sol->getPrimal();
  double* LB_copy   = new double[numvars];
  double* UB_copy   = new double[numvars];
  int err           = 0;

  saveBounds_(LB_copy, UB_copy, numvars);
  // fix bounds for binary variables
  for (VariableConstIterator v_iter=p_->varsBegin();
      v_iter!=p_->varsEnd(); ++v_iter, ++i) {
    if ((*v_iter)->getType() == Binary) {
      p_->changeBound(i, x[i], x[i]);
    }
  }
  //solve the original problem with modified bounds
  e_->clear();
  e_->load(p_);
  e_->solve();
  ++(stats_->numNLPs);
  x                         = e_->getSolution()->getPrimal();
  stats_->bestObjValue      = p_->getObjValue(x, &err);
  // construct a new solution to the original problem
  ConstSolutionPtr original_sol  = 
    (ConstSolutionPtr) new Solution(stats_->bestObjValue, x, p_);
  s_pool->addSolution(original_sol);
  logger_->msgStream(LogInfo) << me_ << "Adding solution to original sol pool"
    << std::endl; 
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Solution value is = " 
      << stats_->bestObjValue << std::endl
      << me_ << "Feasible solution is" << std::endl;
  original_sol->write(logger_->msgStream(LogDebug2));
#endif

  restoreBounds_(LB_copy, UB_copy, numvars);

  if (LB_copy) {
    delete [] LB_copy;
  }

  if (UB_copy) {
    delete [] UB_copy;
  }
}



bool FeasibilityPump::cycle_(double find_value)
{
  // search for the hash value in the previously visited solution's hash value
  assert(!hashVal_.empty());
  for (UInt i=0; i<hashVal_.size()-1; ++i) {
    if (fabs(find_value - hashVal_[i]) < intTol_) {
      ++(stats_->numCycles); 
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "Cycling detected" << std::endl;
#endif
      return true;
    }
  }
  return false;
}


double FeasibilityPump::hash_()
{
  double hash_value = 0;
  DoubleVector::iterator it_rand, it_sol;
  for (it_rand=random_.begin(), it_sol=roundedSol_.begin(); 
          it_rand!=random_.end(); ++it_rand, ++it_sol) {
    hash_value += (*it_rand) * (*it_sol);
  }
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "Hash value for rounded solution  = "
    << hash_value << std::endl;
#endif
  hashVal_.push_back(hash_value);
  return hash_value;
} 


void FeasibilityPump::implementFP_(const double* x, SolutionPoolPtr s_pool)
{
  ConstSolutionPtr sol; 
  double hash_val;
  UInt n_to_flip;
  UInt k;
  
  bool cont_FP         = true;
  UInt max_iter        = 100;
  UInt max_cycle       = 300;
  UInt min_flip        = 3;
  ProblemPtr prob      = p_->clone();

  e_->load(prob);
  while (cont_FP && stats_->numNLPs < max_iter 
      && stats_->numCycles < max_cycle) {
    constructObj_(prob, sol);
    e_->solve();
    ++(stats_->numNLPs);
    sol = e_->getSolution();
#if SPEW
    prob->write(logger_->msgStream(LogDebug2));
    sol->write(logger_->msgStream(LogDebug2));
#endif
    x          =  sol->getPrimal();
    cont_FP    = isFrac_(x);
    hash_val   = hash_();
    k          = std::max(min_flip, (UInt) ceil(sol->getObjValue()));
    n_to_flip  = std::min(k, p_->getSize()->bins);
    if (cycle_(hash_val)) {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "Cycling detected" << std::endl;
#endif
      perturb_(hash_val, n_to_flip);
    }
  }

  if (!cont_FP) {
    // make solution to the original problem and then add to sol pool
    convertSol_(s_pool, sol);
  }

}


bool FeasibilityPump::isFrac_(const double* x)
{
  VariablePtr variable;
  double value;
  double fractional;
  bool is_frac = false; 
  UInt i = 0;
  UInt num_frac = 0;
  // remove the violated variables from previous solution
  for (VariableConstIterator v_iter=p_->varsBegin(); 
      v_iter!=p_->varsEnd(); ++v_iter, ++i) {
    variable = *v_iter;
    value = x[i];
    if (variable->getType() == Binary || variable->getType() == Integer) {
      fractional = fabs(floor(value + 0.5)-value);
#if SPEW
        variable->write(logger_->msgStream(LogDebug2));
        logger_->msgStream(LogDebug2) << me_  << "value of variable " 
          << i << " is "<< value << std::endl;
#endif
      if (fractional > intTol_) {
        roundedSol_[i] = floor(value + 0.5);
        is_frac = true;
        ++num_frac;
      } else {
      roundedSol_[i] = value;
      } 
    } else {
      roundedSol_[i] = value;
    }
  }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Number of fractionals = " 
    << num_frac << std::endl;
#endif
  return is_frac;
}


void FeasibilityPump::perturb_(double hash_val, UInt n_to_flip)
{
  VariablePtr variable;

  do {
    UInt i = 0;
    VarVector to_flip;
    VarVector::iterator it_flip;
    
    to_flip = selectToFlip_(n_to_flip);
    it_flip = to_flip.begin();
    for (VariableConstIterator v_iter=p_->varsBegin();
        v_iter!=p_->varsEnd(); ++v_iter, ++i) {
      variable = *v_iter;
      if (variable == *it_flip) {
        roundedSol_[i] = (roundedSol_[i] < intTol_) ? 1 : 0;
        ++it_flip;
      }
    }
    hash_val = hash_();
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "Number of variables flipped "
      << n_to_flip << std::endl;
#endif
  } while(cycle_(hash_val));
}


void FeasibilityPump::restoreBounds_(double* LB_copy, double* UB_copy, UInt vars)
{
  for (UInt i=0; i<vars; ++i, ++LB_copy, ++UB_copy) {
    p_->changeBound(i, Lower, *LB_copy);
    p_->changeBound(i, Upper, *UB_copy);
  }
}


void FeasibilityPump::saveBounds_(double* LB_copy, double* UB_copy, UInt vars)
{
  VariablePtr variable;
  for (UInt i=0; i<vars; ++i, ++LB_copy, ++UB_copy) {
    variable = p_->getVariable(i);
    *LB_copy = variable->getLb();
    *UB_copy = variable->getUb();
  }
}


VarVector FeasibilityPump::selectToFlip_(UInt n_to_flip)
{
  double U;
  VarVector bin_to_flip;
  UInt t = 0, m = 0;
  UInt num_bins = bins_.size();
  
  while (m < n_to_flip) {
    U = (double) rand()/RAND_MAX;
    if (U * (num_bins - t) >= (n_to_flip - m)) {
      ++t;
    } else {
#if SPEW
      logger_->msgStream(LogDebug2) << me_ << "Will flip variable" << std::endl;
      bins_[t]->write(logger_->msgStream(LogDebug2));
#endif
      bin_to_flip.push_back(bins_[t]);
      ++t, ++m;
    }
  }
  return bin_to_flip;
}


bool FeasibilityPump::shouldFP_()
{
  ConstProblemSizePtr p_size = p_->getSize();
  FunctionType f_type = p_size->objType;
  if (p_size->nonlinCons == 0 && p_size->ints == 0 
      && (f_type == Linear || f_type == Quadratic || f_type == Constant)) {
    return true;
  } else {
    return false;
  }
}


void FeasibilityPump::solve(NodePtr, RelaxationPtr, SolutionPoolPtr s_pool)
{
  EngineStatus status;
  ConstSolutionPtr sol;
  const double* x;
  if(!shouldFP_()) {
    logger_->msgStream(LogInfo) << me_ << "Skipping" << std::endl;
    return;
  }
  e_->load(p_);
  timer_->start();
  status = e_->solve();
  
  if (status != ProvenOptimal && status != ProvenLocalOptimal
      && status != ProvenFailedCQFeas && status != FailedFeas) {
    return;
  }

  sol = e_->getSolution();
  x   = sol->getPrimal();
#if SPEW
  p_->write(logger_->msgStream(LogDebug2));
  sol->write(logger_->msgStream(LogDebug2));
#endif
  e_->clear();
  logger_->msgStream(LogInfo) << me_ << "Starting" << std::endl;
  // now implement the FP heuristic
  if (isFrac_(x)) {
    implementFP_(x, s_pool);
  } else {
  logger_->msgStream(LogInfo) << me_ << "Adding solution to original sol pool"
    << std::endl; 
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "Feasible Solution found" 
      << std::endl << me_ << "Solution value is "
      << sol->getObjValue() << std::endl;
#endif
    s_pool->addSolution(sol);
  }

  logger_->msgStream(LogInfo) << me_ << "Over" << std::endl;
  stats_->time = timer_->query();
}


void FeasibilityPump::writeStats(std::ostream &out) const
{
  out << me_ << "number of nlps solved         = " << stats_->numNLPs
    << std::endl
    << me_ << "number of cycles              = " << stats_->numCycles
    << std::endl
    << me_ << "numer of errors               = " << stats_->errors
    << std::endl
    << me_ << "total time taken              = " << stats_->time
    << std::endl;
  if (stats_->bestObjValue < INFINITY) {
    out << me_ << "best objective value          = " 
      << stats_->bestObjValue << std::endl;
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
