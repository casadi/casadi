//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
 
/**
 * \file LinFeasPump.cpp
 * \brief Implements the class Feasibility Pump.
 * 
 * \author Jayash Koshal, Argonne National Laboratory
 * 
 * Define the Linear Feasibility Pump class for generating a 
 * feasible solution using FP heuristic for MINLPs.
 */

#include "MinotaurConfig.h"
#include "Engine.h"
#include "Variable.h"
#include "Environment.h"
#include "Logger.h"
#include "Node.h"
#include "Operations.h"
#include "Option.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "Constraint.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Objective.h"
#include "LinFeasPump.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "LinearHandler.h"
#include "QGHandler.h"
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace Minotaur;

//#define SPEW 1

const std::string LinFeasPump::me_ = "Linear Feas Pump: ";

LinFeasPump::LinFeasPump(EnvPtr env, ProblemPtr p, EnginePtr e1,
    EnginePtr e2)
  : FeasibilityPump(env, p, e1, e2),
    gradientObj_(NULL),
    lh_(LinearHandlerPtr()),
    lpE_(e2),
    objConstraint_(ConstraintPtr()),
    objVar_(VariablePtr()),
    olfClone_(LinearFunctionPtr()),
    qh_(QGHandlerPtr()),
    r_(RelaxationPtr())
{
  // allocate space for gradient of objective function
  gradientObj_               = new double[p_->getNumVars()];

  logger_ = (LoggerPtr) new Logger((LogLevel) env_->getOptions()->
     findInt("heur_log_level")->getValue());

  // initialize the statistics
  statsLFP_                  = new LinFeasStats();
  statsLFP_->bestObjValue    = INFINITY;
  statsLFP_->numLPs          = 0;

}


LinFeasPump::~LinFeasPump()
{
  if (gradientObj_) {
    delete [] gradientObj_;
  }
  delete statsLFP_;
}


void LinFeasPump::constructObj_(ProblemPtr, ConstSolutionPtr)
{
  double value, lb, ub;
  VariablePtr variable;
  UInt i;
  //double obj_relaxation;
  double obj_weight;
  double constant;
  FunctionPtr f;
  LinearFunctionPtr olf_mod;

  //obj_relaxation = sol->getPrimal()[r_->getNumVars()-1];
  obj_weight     = 0; //std::min(1/(fabs(obj_relaxation)+1e-7), 0.1);
  obj_weight     = 0.0;
  i              = 0;
  constant       = 0;
  olf_mod        = olfClone_->clone();
  (*olf_mod)    *= obj_weight;

  for (VariableConstIterator v_iter=r_->varsBegin();
      v_iter!=r_->varsEnd(); ++v_iter, ++i) {
    variable = *v_iter;
    if (variable->getType() == Binary) {
      value  = roundedSol_[i];
      lb     = variable->getLb();
      ub     = variable->getUb();
      if (fabs(value - lb) > intTol_) { 
        olf_mod->addTerm(variable,-1.0);
        constant += ub;
      } else if (fabs(value - ub) > intTol_) {
        olf_mod->addTerm(variable,1.0);
        constant -= lb;
      } else {
        // add a new variable with coeff 1
        // add two constraints for absolute value
      }      
    }
  }
  f = (FunctionPtr) new Function(olf_mod);
  r_->changeObj(f, constant);
}


void LinFeasPump::implementFP_(const double*, SolutionPoolPtr s_pool)
{
  ConstSolutionPtr sol;
  const double* x_lp;
  double hash_val, f_nlp;
  UInt n_to_flip, k;
  SeparationStatus sep_status = SepaContinue;
  EngineStatus lp_status      = EngineUnknownStatus;
  NodePtr node                = NodePtr();
  bool to_continue            = true;
  bool is_feasible            = false;
  bool sol_found              = false;
  bool is_prob_infeasible     = prepareLP_();
  bool should_separate        = true;
  UInt max_NLP                = 10;
  UInt max_LP                 = 2000;
  UInt max_cycle              = 1000;
  UInt min_flip               = 2;
  UInt max_non_zero_obj       = 500;
  double inf_meas             = 0.0;
  int err;
  
  e_->setOptionsForSingleSolve();
  lpE_->solve();
  f_nlp = lpE_->getSolutionValue();
  sol   = lpE_->getSolution();
  should_separate = (p_->getObjective()->getFunction()->
                     getNumVars() > max_non_zero_obj) ? false : true;

  while(!is_feasible && stats_->numNLPs < max_NLP && statsLFP_->numLPs < max_LP
        && stats_->numCycles < max_cycle) {
    while(to_continue && statsLFP_->numLPs < max_LP 
        && stats_->numCycles < max_cycle) { 
      sol_found = false;
      constructObj_(r_, sol);
      lp_status = lpE_->solve();
      ++(statsLFP_->numLPs);
      if (!(lp_status == ProvenOptimal || lp_status == ProvenLocalOptimal)) {
        logger_->msgStream(LogDebug) << me_ << "LP relaxation is infeasible." 
          << std::endl;
        return;
      }
      sol         = lpE_->getSolution();
      x_lp        = sol->getPrimal();
      to_continue = isFrac_(x_lp);
      k           = std::max(min_flip, (UInt) ceil(sol->getObjValue()));
      n_to_flip   = std::min(k, p_->getSize()->bins);

      if (!to_continue) {
        is_feasible = qh_->isFeasible(sol, r_, is_prob_infeasible, inf_meas);
        if (is_feasible) {
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "LP soln is feasible "
            << "to NLP" << std::endl;
#endif
          err = 0;
          statsLFP_->bestObjValue = p_->getObjective()->eval(x_lp, &err);
          convertSol_(s_pool, sol);
          sol_found = true;
        } else {
#if SPEW
          logger_->msgStream(LogDebug) << me_ << "LP soln is not feasible "
            << "to NLP." << std::endl;
#endif
        }
        break;
      }
      hash_val = hash_();
      if (cycle_(hash_val)) {
        perturb_(hash_val, n_to_flip);
      }
    } 
    if (!to_continue) {
      if (false==sol_found) {
        ++(stats_->numNLPs);
        qh_->separate(sol, node, r_, 0, s_pool, &sol_found, &sep_status);
        to_continue = true; //reset to continue after separating
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "separation status = " 
          << sep_status << std::endl;
        //r_->write(logger_->msgStream(LogDebug));
        logger_->msgStream(LogDebug) << me_ << 
          "sol_found = " << sol_found << std::endl;
#endif
        if (SepaError==sep_status) {
          break; //exit pump
        }
        if (sol_found) {
          is_feasible = true;
        }
      }
    }
    if (is_feasible) {
#if SPEW
      logger_->msgStream(LogInfo) << me_ << "best solution value = " 
        << s_pool->getBestSolutionValue() << std::endl;
#endif
      if (getSolGap_(f_nlp, s_pool->getBestSolutionValue()) > 10 && 
          true==should_separate) {
        separatingCut_(f_nlp, s_pool);
        is_feasible = false;
      } else {
        break;
      }
    }
  }
  e_->setOptionsForRepeatedSolve();
}


double LinFeasPump::getSolGap_(double f_nlp, double f_feas)
{
  return (f_feas - f_nlp)/(fabs(f_feas) + 1e-6)*100;
}


bool LinFeasPump::prepareLP_()
{
  bool is_inf = false;
 //  std::vector<bool> c_list(p_->getNumCons(), false);
  
  r_  = (RelaxationPtr) new Relaxation();
  lh_ = (LinearHandlerPtr) new LinearHandler(env_, p_);
  qh_ = (QGHandlerPtr)  new QGHandler(env_, p_, e_);
  
  lh_->relaxInitInc(r_, &is_inf);
  if (is_inf == true) {
    return is_inf;
  }

  qh_->relaxInitInc(r_, &is_inf);
  if (is_inf == true) {
    return is_inf;
  }

  if (Constant==p_->getObjective()->getFunctionType() 
      || Linear==p_->getObjective()->getFunctionType()) {
    // all coeffs are zero, but it is not NULL
    olfClone_ = (LinearFunctionPtr) new LinearFunction(); 
  } else {
    olfClone_ = r_->getObjective()->getLinearFunction();
    objVar_   = r_->getVariable(r_->getNumVars()-1);
  }
  
  lpE_->load(r_);

  return is_inf;
}


void LinFeasPump::separatingCut_(double f_nlp, SolutionPoolPtr s_pool)
{
  double new_bnd;
  LinearFunctionPtr lf;
  FunctionPtr obj_f;
  double obj_weight          = 0.6; // increase(<1) to get better improvement
  double f_feas              = s_pool->getBestSolutionValue();

  new_bnd = obj_weight*f_nlp + (1-obj_weight)*f_feas; 
  if (objVar_) {
    r_->changeBound(objVar_, Upper, new_bnd);
    logger_->msgStream(LogDebug) << me_ << "upper bound on objective value is "
                                 << new_bnd << std::endl;
  } else {
    if (!objConstraint_) {
      lf = p_->getObjective()->getFunction()->getLinearFunction()->clone();
      obj_f  = (FunctionPtr) new Function(lf);
      objConstraint_ = r_->newConstraint(obj_f, -INFINITY, new_bnd, 
                                         "obj_sep_cut");
    } else {
      r_->changeBound(objConstraint_, Upper, new_bnd);
    }
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "upper bound on objective value is "
   << new_bnd << std::endl;
#endif
  }

}


bool LinFeasPump::shouldFP_()
{
  ConstConstraintPtr c;

  if (env_->getOptions()->findInt("LinFPump")->getValue()<0) {
    return false;
  }
  if (p_->getSize()->ints > 0) {
    logger_->msgStream(LogInfo) << me_ << "skipping because of integer "
      "variables" << std::endl;
    return false;
  }
  if (p_->getSize()->bins < 31 && 
      0==env_->getOptions()->findInt("LinFPump")->getValue()) {
    logger_->msgStream(LogInfo) << me_ << "skipping because of too few binary "
      "variables" << std::endl;
    return false;
  }

  for (ConstraintConstIterator c_it=p_->consBegin(); 
      c_it!=p_->consEnd(); ++c_it) {
    c = *c_it;
    if (c->getFunctionType() != Linear) {
      if (c->getLb() > -INFINITY && c->getUb() < INFINITY) {
        logger_->msgStream(LogInfo) << me_ << "skipping because of nonlinear "
          << "equality or range constraint" << std::endl;
        return false;
      } 
    }
  }
  
  return true;
}


void LinFeasPump::solve(NodePtr, RelaxationPtr, SolutionPoolPtr s_pool)
{
  const double* x = 0;
  if (!shouldFP_()) {
    return;
  }

  logger_->msgStream(LogInfo) << me_ << "starting" << std::endl;
  timer_->start();
  implementFP_(x, s_pool); 
  logger_->msgStream(LogInfo) << me_ << "over" << std::endl;
  stats_->time = timer_->query();
  timer_->stop();
}


void LinFeasPump::writeStats(std::ostream &out) const
{
  FeasibilityPump::writeStats(out);
  logger_->msgStream(LogInfo) 
    << me_ << "number of LPs solved          = " << statsLFP_->numLPs
    << std::endl;
  if (statsLFP_->bestObjValue < INFINITY) {
    logger_->msgStream(LogInfo) 
      << me_ << "Best objective value          = " << statsLFP_->bestObjValue
      << std::endl;
  }
}
