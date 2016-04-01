//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
 
/**
 * \file MINLPDiving.cpp
 * \brief Define the MINLPDiving class for diving heuristics for MINLPs.
 * \author Jayash Koshal, Argonne National Laboratory
 *
 * Implements the class MINLPDiving.
 */

#include <cmath> // for INFINITY
#include <iomanip>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "Function.h"
#include "LinearHandler.h"
#include "Logger.h"
#include "MINLPDiving.h"
#include "Node.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "ProblemSize.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "VarBoundMod.h"
#include "Variable.h"

using namespace Minotaur;

//#define SPEW 0

const std::string MINLPDiving::me_ = "MINLP Diving Heuristic: "; 

MINLPDiving::MINLPDiving(EnvPtr env, ProblemPtr p, EnginePtr e)
: e_(e), 
  env_(env), 
  gradientObj_(NULL),
  intTol_(1e-5),
  lh_(0),
  maxNLP_(100),
  maxSol_(2), 
  nSelector_(4),
  p_(p), 
  stats_(NULL), 
  timer_(env_->getNewTimer())
{
  for (UInt i=0; i<p_->getNumVars(); ++i) {
    avgDual_.push_back(0.);
  }
    
  // allocate space for gradient of objective function
  gradientObj_              = new double[p_->getNumVars()];

  logger_ = (LoggerPtr) new Logger((LogLevel) env_->getOptions()->
      findInt("heur_log_level")->getValue());

  //DivingheurStats stats_
  stats_                    = new DivingheurStats();
  for (UInt i=0; i<nSelector_; ++i) {
    stats_->numNLPs[i]      = 0;
    stats_->time[i]         = 0;
    stats_->iterations[i]   = 0;
    stats_->numInfeas[i]    = 0;
    stats_->errors[i]       = 0;
    stats_->numSol[i]       = 0;
  }
  stats_->totalNLPs         = 0;
  stats_->totalSol          = 0;
  stats_->numLocal          = 1;
  stats_->best_obj_value    = INFINITY;
  stats_->totalTime         = 0; 
  
}


MINLPDiving::~MINLPDiving()
{
  if (gradientObj_) {
    delete [] gradientObj_;
  }
  if (timer_) {
    delete timer_;
  }
  delete stats_;
  if (lh_) {
    delete lh_;
  }
  lastNodeMods_.clear();
}


void MINLPDiving::backtrack_(UInt n_moded)
{
  std::stack<VarBoundModPtr> temp;
  VarBoundModPtr varmod;
  VarBoundModPtr temp_varmod;
  VariablePtr variable;
  double old_bound_val;

  for (ModVector::reverse_iterator m_iter=lastNodeMods_.rbegin(); 
       m_iter!=lastNodeMods_.rend(); ++m_iter) {
    (*m_iter)->undoToProblem(p_);
  }
  for (UInt i=0; i<n_moded; ++i) {
    varmod = mods_.top();
    variable = varmod->getVar();
    old_bound_val = varmod->getNewVal();
    if (varmod->getLU() == Lower) {
      //previously changed bound was lower
      // decrement the old bound value by 1
      temp_varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, 
                                                     old_bound_val - 1);
    } else {
      // previously changed bound was upper
      // increment the old bound value by 1
      temp_varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, 
                                                     old_bound_val + 1);
    }
    temp.push(temp_varmod); // hold the modification in a temp stack 
    mods_.pop();
  }

  // put the new modification back in the original stack
  while (!temp.empty()) {
    temp_varmod = temp.top();
    temp_varmod->applyToProblem(p_);
    mods_.push(temp_varmod);
    temp.pop();
  }
}


UInt MINLPDiving::FracBounds_(UInt numfrac, const double* x, 
                              Direction d, Order o)
{
  VariablePtr variable;
  VarBoundModPtr varmod;
  UInt id;
  double value;
  double new_bound;
  UInt begin      = 0;
  UInt end        = 0;
  UInt change_wan = (UInt) ceil( (double) numfrac/4);
  UInt changes    = 0;
  ModVector dummy;
  NodePtr node = NodePtr();
  SolutionPoolPtr s_pool = SolutionPoolPtr(); // NULL
  SolveStatus status;

  // get the score of violated variables according to their fractional part
  getScore_(x, Fractional);
  // sort the vector of index according to the fractional part of the variables
  sort_(0,violated_.size()-1);
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Changing bound of " << change_wan 
    << " variables" << std::endl;
#endif 

  switch (o) {
   case (Least) : // least fractional to be fixed => forward counting
     begin = 0; end = change_wan;
     break;

   case (Most) : // most fractional to be fixed => reverse counting
     begin = numfrac-change_wan; end = numfrac;
     break;

   default:
     assert(!"Order unknown");
  }

  for (UInt i=begin; i<end; ++i) {
    id = violated_[i];
    variable = p_->getVariable(id);
    value = x[id];
    new_bound = rounding_(value, d); // get the new bound with direction
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "value of variable " 
      << variable->getName() << " " << value << std::endl;
    logger_->msgStream(LogDebug) << me_ << "value of new bound " 
      << d << "\t" << new_bound << std::endl;
#endif
    // fix the variable by changing the lower or upper bound to new_bound
    if (d == Floor) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    } else if (d == Ceil) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value < new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value > new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    }

    varmod->applyToProblem(p_);
    mods_.push(varmod); 
    lastNodeMods_.push_back(varmod);
    status = Started;
    lh_->simplePresolve(p_, SolutionPoolPtr(), lastNodeMods_, status);

#if SPEW
    logger_->msgStream(LogDebug) << me_ << "changed bound for variable "; 
    variable->write(logger_->msgStream(LogDebug));
    logger_->msgStream(LogDebug) << me_ << "changed bound for index " 
      << id << std::endl << std::endl;
#endif
    ++changes;
  }
  return changes;
}


void MINLPDiving::getScore_(const double* x, Scoretype s)
{
  VariablePtr variable;
  double value;
  UInt i = 0;
  double fractional;
  ObjectivePtr obj;
  UInt numcons;
  double vl_score;
  int error;
  score_.clear(); // remove the score of violated variables from prev sol
  switch (s) {
   case (LexBound)   :
     for (VariableConstIterator v_iter=p_->varsBegin();
         v_iter!=p_->varsEnd(); ++v_iter, ++i) {
       variable = *v_iter;
       if (variable->getType() == Binary || variable->getType() == Integer) {
         value = x[i];
         fractional = fabs(floor(value+0.5)-value);
         if (fractional > MINLPDiving::intTol_) {
           score_.push_back(i);
         }
       }
     }
#if SPEW
     logger_->msgStream(LogDebug) << " Lexicographic Scoring " << std::endl;
#endif
     break;

   case (Fractional) : 
     for (VariableConstIterator v_iter=p_->varsBegin(); 
         v_iter!=p_->varsEnd(); ++v_iter, ++i) {
       variable = *v_iter;
       if (variable->getType() == Binary || variable->getType() == Integer) {
         value = x[i];
         fractional = fabs(floor(value+0.5)-value);
         if (fractional > intTol_) {
           score_.push_back(fractional);
         }
       }
     }
#if SPEW
     logger_->msgStream(LogDebug)  << " Fractional Scoring " << std::endl;
#endif
     break;

   case (VectorLength) : 
     obj = p_->getObjective();
     std::fill(gradientObj_, gradientObj_+p_->getNumVars(), 0);
     obj->evalGradient(x, gradientObj_, &error);
     for (VariableConstIterator v_iter=p_->varsBegin(); 
         v_iter!=p_->varsEnd(); ++v_iter, ++i) {
       variable = *v_iter;
       if (variable->getType() == Binary || variable->getType() == Integer) {
         value = x[i];
         fractional = fabs(floor(value+0.5)-value);
         if (fractional > intTol_) {
           numcons = variable->getNumCons();
           numcons = (numcons > 0) ? numcons : numcons + 1;
           vl_score = gradientObj_[i] * fractional/numcons; 
           score_.push_back(vl_score);
         }
       }
     }
#if SPEW
     logger_->msgStream(LogDebug) << " Vector Length Scoring " << std::endl;
#endif
     break;

   case (ReducedCost) :
     score_.resize(p_->getNumVars());
     std::copy(avgDual_.begin(), avgDual_.end(), score_.begin());
#if SPEW
     logger_->msgStream(LogDebug) << " Reduced Cost Scoring " << std::endl;
#endif
     break;

   default: break;
  }
}


void MINLPDiving::implementDive_(int i, const double* x, SolutionPoolPtr s_pool)
{
  ConstSolutionPtr sol;
  Direction d;
  Order o;
  EngineStatus status;
  UInt backtrack     = 0;
  UInt min_vlength   = 5;
  UInt numfrac       = isFrac_(x);
  FuncPtr f          = selectHeur_(i, d, o);
  UInt n_moded;

  if (f == &MINLPDiving::VectorLength_) {
    if (false==vectorFlag_(min_vlength)) { //vector length diving?
      return; // skip vector length diving
    }
  }

  lastNodeMods_.clear();
  n_moded  = (this->*f)(numfrac, x, d, o);
  while (stats_->totalNLPs < maxNLP_) {
    status = e_->solve();
    ++(stats_->numNLPs[i/8]);
    ++(stats_->totalNLPs);
    if (EngineError == status) {
      e_->clear();  // reset the starting point
      e_->load(p_);
      ++(stats_->errors[i/8]);
    }
    if (status == ProvenLocalOptimal || status == ProvenOptimal 
        || status == ProvenFailedCQFeas || status == FailedFeas) {
      sol = e_->getSolution();
      ++(stats_->numLocal);
      if (stats_->best_obj_value - 1e-6 < sol->getObjValue()) {
#if SPEW
        logger_->msgStream(LogDebug) << me_ 
          << "current solution worse than ub. Returning." << std::endl; 
#endif
        return;
      }
      updateAvgDual_(sol);
      x = sol->getPrimal();
      numfrac = isFrac_(x); // number of fractional vars in current solution
      if (0==numfrac) {  
#if SPEW
        logger_->msgStream(LogDebug) << me_ << "Feasible Solution" << std::endl;
        sol->write(logger_->msgStream(LogDebug2));
#endif
        logger_->msgStream(LogInfo) << me_ << "Updating the solution value to " 
          << sol->getObjValue() << std::endl;
        stats_->best_obj_value = sol->getObjValue();
        s_pool->addSolution(sol);
        ++(stats_->numSol[i/8]);
        ++(stats_->totalSol);
        return;
      } else {
        // dive further down by rounding variable in current solution "x"
        lastNodeMods_.clear();
        n_moded = (this->*f)(numfrac, x, d, o);
        backtrack = 0;
      }
    } else if (0 < backtrack) {
      ++(stats_->numInfeas[i/8]);
      return;
    } else {
      // backtrack to the other child of the parent
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "Backtracking " << std::endl;
#endif
      ++backtrack;
      backtrack_(n_moded);
      ++(stats_->numInfeas[i/8]);
    }
  }
}


UInt MINLPDiving::isFrac_(const double* x)
{
  VariablePtr variable;
  double value;
  UInt i = 0;
  double fractional;
  // remove the violated variables from previous solution
  violated_.clear();
  for (VariableConstIterator v_iter=p_->varsBegin(); 
      v_iter!=p_->varsEnd(); ++v_iter, ++i) {
    variable = *v_iter;
    if (variable->getType() == Binary || variable->getType() == Integer) {
      value = x[i];
      fractional = fabs(floor(value+0.5)-value);
      if (fractional > intTol_) {
        violated_.push_back(i);    // variable i is fractional
        //score_.push_back(fractional);
      }
    }
  }
  return violated_.size();
}


UInt MINLPDiving::LexBounds_(UInt numfrac, const double* x,
    Direction d, Order o)
{
  VariablePtr variable;
  VarBoundModPtr varmod;
  UInt id;
  double new_bound;
  double value;
  UInt begin      = 0;
  UInt end        = 0;
  UInt change_wan = (UInt) ceil( (double) numfrac/4);
  UInt changes    = 0;
  ModVector dummy;
  NodePtr node = NodePtr();
  SolutionPoolPtr s_pool = SolutionPoolPtr(); // NULL
  SolveStatus status;
  
  getScore_(x, LexBound);
  sort_(0, violated_.size()-1);
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Changing bound of " << change_wan 
    << " variables" << std::endl;
#endif 
  // set up begin and end for the selection of variable
  switch (o) {
   case (Least) : // least fractional to be fixed => forward counting
     begin = 0; end = change_wan; 
     break;
   
   case (Most) : // most fractional to be fixed => reverse counting
     begin = numfrac-change_wan; end = numfrac;
     break;

   default:
     assert(!"Order unknown");
  }

  // based on begin and end from switch change the bound of variables
  for (UInt i=begin; i<end; ++i) {
    id = violated_[i];
    variable = p_->getVariable(id);
    value = x[id];
    if (fabs(variable->getUb()-variable->getLb())<1e-7) {
      continue;
    }
    new_bound = rounding_(value, d); // get the new bound with direction
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "value of variable " << 
      variable->getName() << " " << value << std::endl;
    logger_->msgStream(LogDebug) << me_ << "value of new bound " 
      << d << "\t" << new_bound << std::endl;
#endif
    // fix the variable by changing the lower or upper bound to new_bound
    if (d == Floor) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    } else if (d == Ceil) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value < new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value > new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    }

    varmod->applyToProblem(p_);
    mods_.push(varmod); 
    lastNodeMods_.push_back(varmod);
    status = Started;
    lh_->simplePresolve(p_, SolutionPoolPtr(), lastNodeMods_, status);

#if SPEW
    logger_->msgStream(LogDebug) << me_ << "changed bound for variable "; 
    variable->write(logger_->msgStream(LogDebug));
    logger_->msgStream(LogDebug) << std::endl;
    logger_->msgStream(LogDebug) << me_ << "changed bound for index " 
      << id << std::endl << std::endl;
#endif
    ++changes;
  }
  return changes;
}




UInt MINLPDiving::ReducedCost_(UInt numfrac, const double* x, 
    Direction d, Order o)
{
  VariablePtr variable;
  VarBoundModPtr varmod;
  UInt id;
  double new_bound;
  double value;
  UInt begin      = 0;
  UInt end        = 0;
  UInt change_wan = (UInt) ceil( (double) numfrac/4);
  UInt changes    = 0;
  ModVector dummy;
  NodePtr node = NodePtr();
  SolutionPoolPtr s_pool = SolutionPoolPtr(); // NULL
  SolveStatus status;

  getScore_(x, ReducedCost);
  std::setprecision(8);
  sort_(0,violated_.size()-1);
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Changing bound of " << change_wan 
    << " variables" << std::endl;
#endif 
  // set up begin and end for the selection of variable
  switch (o) {
   case (Least) : // least fractional to be fixed => forward counting
     begin = 0; end = change_wan; 
     break;
   
   case (Most) : // most fractional to be fixed => reverse counting
     begin = numfrac-change_wan; end = numfrac;
     break;

   default:
     assert(!"Order unknown");
  }

  // based on begin and end from switch change the bound of variables
  for (UInt i=begin; i<end; ++i) {
    id         = violated_[i];
    variable   = p_->getVariable(id);
    value      = x[id];
    new_bound = rounding_(value, d); // get the new bound with direction
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "value of variable " 
      << variable->getName() << " " << value << std::endl;
    logger_->msgStream(LogDebug) << me_ << "value of new bound " 
      << d << "\t" << new_bound << std::endl;
#endif
    // fix the variable by changing the lower or upper bound to new_bound
    if (d == Floor) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    } else if (d == Ceil) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value < new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value > new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    }

    varmod->applyToProblem(p_);
    mods_.push(varmod); 
    lastNodeMods_.push_back(varmod);
    status = Started;
    lh_->simplePresolve(p_, SolutionPoolPtr(), lastNodeMods_, status);

#if SPEW
    logger_->msgStream(LogDebug) << me_ << "changed bound for variable "; 
    variable->write(logger_->msgStream(LogDebug));
    logger_->msgStream(LogDebug) << std::endl;
    logger_->msgStream(LogDebug) << me_ << "changed bound for index " 
      << id << std::endl << std::endl;
#endif
    ++changes;
  }

  return changes;
}


void MINLPDiving::restoreBounds_(double* LB_copy, double* UB_copy, UInt vars)
{
  for (UInt i=0; i<vars; ++i, ++LB_copy, ++UB_copy) {
    p_->changeBound(i, Lower, *LB_copy);
    p_->changeBound(i, Upper, *UB_copy);
  }
}


double MINLPDiving::rounding_(double value, Direction d)
{
  switch (d) {
       case (Floor)    : return floor(value); 
                         break;
       case (Ceil)     : return ceil(value); 
                         break;
       case (Nearest)  : if (fabs(value-floor(value)) < 0.5) {
                          // floor is the nearest
                          return floor(value);
                         } else {
                           return ceil(value);
                         }
                         break;
       case (Farthest) : if (fabs(value-floor(value)) < 0.5) { 
                           // ceil is the farthest
                           return ceil(value);
                         } else {
                           return floor(value);
                         }
                         break;
      }
  assert(!"unknown case in rounding_");
  return 0;
}


void MINLPDiving::saveBounds_(double* LB_copy, double* UB_copy, UInt vars)
{
  VariablePtr variable;
  for (UInt i=0; i<vars; ++i, ++LB_copy, ++UB_copy) {
    variable = p_->getVariable(i);
    *LB_copy = variable->getLb();
    *UB_copy = variable->getUb();
  }
}


MINLPDiving::FuncPtr MINLPDiving::selectHeur_(int i, Direction &d, Order &o)
{
  switch (i%4) {
   case 0 : d = Floor;
            break;

   case 1 : d = Ceil;
            break;

   case 2 : d = Nearest;
            break;

   case 3 : d = Farthest;
            break;

   default: d = Ceil;
            break;
  }

  if (i%8 < 4) {
    o = Least;
  } else {
    o = Most;
  }

  if (i < 8) {
    return &MINLPDiving::FracBounds_; 
  } else if (i < 16) {
    return &MINLPDiving::VectorLength_;
  } else if (i < 24) {
    return &MINLPDiving::LexBounds_;
  } else {
    return &MINLPDiving::ReducedCost_;
  } 
}


bool MINLPDiving::shouldDive_()
{
  VariablePtr variable;
  
  UInt min_num_bint             = 5;
  ConstProblemSizePtr prob_size = p_->getSize();
  UInt num_bin_int              = prob_size->bins + prob_size->ints;
  int option = env_->getOptions()->findInt("divheur")->getValue();
  
  switch (option) {

   case 0  : if (num_bin_int <= min_num_bint) {
               return false;
             } else {
               return true;
             }
             break;

   case 1  : return true;
             break;

   default : assert(!"Unknown option for MINLPDiving heuristic");
             return false;
             break;
  }

  return false; // if option == -1
}


void MINLPDiving::solve(NodePtr, RelaxationPtr, SolutionPoolPtr s_pool)
{
  ConstSolutionPtr sol;
  EngineStatus status;
  const double* x;
  
  int num_method         = 32;
  UInt numvars           = p_->getNumVars();
  double* root_x;
  double* root_copy;
  double* LB_copy;
  double* UB_copy;
  logger_->msgStream(LogInfo) << me_ << "Starting" << std::endl;
  timer_->start();
  if (!shouldDive_()) {
    logger_->msgStream(LogInfo) << me_ << "Skipping" << std::endl;
    timer_->stop();
    return;
  }
  root_x             = new double[numvars];
  root_copy          = new double[numvars];
  LB_copy            = new double[numvars];
  UB_copy            = new double[numvars];
  e_->clear();
  e_->load(p_);
  e_->setIterationLimit(7000); // try to run for a loooong time.
  status             = e_->solve();
  sol                = e_->getSolution();
  stats_->totalTime  = timer_->query();
  x                  = sol->getPrimal();
  std::copy(x, x + numvars, root_x);
  updateAvgDual_(sol);
  e_->setIterationLimit(200);
  if (status == ProvenLocalInfeasible) {
    logger_->msgStream(LogInfo) << me_ << "Root Infeasible" << std::endl;
  } else if (isFrac_(root_x) == 0) {
    stats_->best_obj_value = sol->getObjValue();
    logger_->msgStream(LogInfo) << me_ << "solution value is " 
      << stats_->best_obj_value << std::endl;
    s_pool->addSolution(sol);
  } else {
    lh_ = new LinearHandler(env_, p_);
    saveBounds_(LB_copy, UB_copy, numvars);
    // loop over the methods starts here
    for (int i=0; i<num_method && stats_->totalSol < maxSol_; ++i) {
      logger_->msgStream(LogDebug) << me_<< "diving method "
        << i << std::endl;
      std::copy(root_x, root_x + numvars, root_copy); 
      x = root_copy;
      implementDive_(i, x, s_pool);
      restoreBounds_(LB_copy, UB_copy, numvars);
      // clear the stack of modification for this heuristic method
      while (!mods_.empty()) {
          mods_.pop();
      }
      if ((i+1)%8 == 0) {
        stats_->time[i/8]  = timer_->query();
        timer_->stop();
        timer_->start();
      }
    } // loop over methods ends here
  }
  e_->resetIterationLimit();
  logger_->msgStream(LogInfo) << me_ << "Over" << std::endl;
  if (root_x){
    delete [] root_x;
  }
  if (LB_copy){
    delete [] LB_copy;
  }
  if (UB_copy){
    delete [] UB_copy;
  }
  if (root_copy){
    delete [] root_copy;
  }
  timer_->stop();
}


void MINLPDiving::sort_(UInt left, UInt right)
{
  UInt i = left;
  UInt j = right;
  double pivot = score_[(left + right)/2];

  while (i <= j) {
    while (score_[i] < pivot) {
      ++i;
    }
    while (score_[j] > pivot) {
      assert (j!=0); // should never come here if j==0, otherwise overflows.
      --j; 
    }
    if (i <= j) {
      std::swap(score_[i],score_[j]);
      std::swap(violated_[i],violated_[j]);
      ++i;
      if (j>0) {
        --j;// should never come here if j==0, otherwise overflows.
      }
    }
  }

  if (left < j) {
    sort_(left,j);
  }
  if (i < right) {
    sort_(i,right);
  }
}


void MINLPDiving::updateAvgDual_(ConstSolutionPtr sol)
{
  DoubleVector::iterator it;
  UInt iter         = stats_->numLocal;
  const double* d   = sol->getDualOfVars();
  for (it=avgDual_.begin(); it!=avgDual_.end(); ++it, ++d) {
    *it = ((*it)*(iter - 1) + *d)/iter;
  }
}


UInt MINLPDiving::VectorLength_(UInt numfrac, const double* x, 
                                Direction d, Order o)
{
  VariablePtr variable;
  VarBoundModPtr varmod;
  UInt id;
  double new_bound;
  double value;
  UInt begin      = 0;
  UInt end        = 0;
  UInt change_wan = (UInt) ceil( (double) numfrac/4);
  UInt changes    = 0;
  ModVector dummy;
  NodePtr node = NodePtr();
  SolutionPoolPtr s_pool = SolutionPoolPtr(); // NULL
  SolveStatus status;

  getScore_(x, VectorLength);
  sort_(0, score_.size()-1);
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "Changing bound of " << change_wan 
    << " variables" << std::endl;
#endif 
  // set up begin and end for the selection of variable
  switch (o) {
   case (Least) : // least fractional to be fixed => forward counting
     begin = 0; end = change_wan; 
     break;
   
   case (Most) : // most fractional to be fixed => reverse counting
     begin = numfrac-change_wan; end = numfrac;
     break;

   default:
     assert(!"Order unknown");
  }

  // based on begin and end from switch change the bound of variables
  for (UInt i=begin; i<end; ++i) {
    id = violated_[i];
    variable = p_->getVariable(id);
    value = x[id];
    new_bound = rounding_(value, d); // get the new bound with direction
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "value of variable " 
      << variable->getName() << " " << value << std::endl;
    logger_->msgStream(LogDebug) << me_ << "value of new bound " 
      << d << "\t" << new_bound << std::endl;
#endif
    // fix the variable by changing the lower or upper bound to new_bound
    if (d == Floor) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    } else if (d == Ceil) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value < new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Lower, new_bound);
    } else if (value > new_bound) {
      varmod = (VarBoundModPtr) new VarBoundMod(variable, Upper, new_bound);
    }

    varmod->applyToProblem(p_);
    mods_.push(varmod); 
    lastNodeMods_.push_back(varmod);
    status = Started;
    lh_->simplePresolve(p_, SolutionPoolPtr(), lastNodeMods_, status);

#if SPEW
    logger_->msgStream(LogDebug) << me_ << "changed bound for variable "; 
    variable->write(logger_->msgStream(LogDebug));
    logger_->msgStream(LogDebug) << std::endl;
    logger_->msgStream(LogDebug) << me_ << "changed bound for index " 
      << id << std::endl << std::endl;
#endif
    ++changes;
  }

  return changes;
}


bool MINLPDiving::vectorFlag_(UInt min_vlength)
{
  UInt num_obj_int      = 0;
  FunctionPtr obj       = p_->getObjective()->getFunction();
  VariablePtr variable;
  // get number of binary and integer variables in the objective
  // if num_obj_int is below a threshold then vector length diving is skipped

  for (VarSetConstIterator it= obj->varsBegin(); it!=obj->varsEnd(); ++it) {
    variable = *it;
    if (variable->getType() == Binary || variable->getType() == Integer) {
      ++num_obj_int ;
    }
  }
  if ((p_->getObjective()->getFunctionType() == Linear 
        || p_->getObjective()->getFunctionType() == Quadratic) 
      && num_obj_int >= min_vlength) {
    return true; 
  } else { 
    return false;
  }
}


void MINLPDiving::writeStats(std::ostream &out) const
{
// write the statistics for MINLP heuristic
  for (UInt i=0; i<nSelector_; ++i) {
    out << me_ << "number of nlps solved      = " << stats_->numNLPs[i] 
      << std::endl
      << me_ << "time taken                 = " << stats_->time[i] 
      << std::endl
      << me_ << "number of solutions found  = " << stats_->numSol[i] 
      << std::endl
      << me_ << "number of Infeasible NLPs  = " << stats_->numInfeas[i] 
      << std::endl
      << me_ << "number of Errors           = " << stats_->errors[i] 
      << std::endl
      << me_ << "number of iterations       = " << stats_->iterations[i] 
      << std::endl << std::endl;
    stats_->totalTime += stats_->time[i];
  }
  if (stats_->best_obj_value < INFINITY) {
    logger_->msgStream(LogInfo) << me_ << "Best feasible sol value    = "
      << stats_->best_obj_value << std::endl;
  }
  logger_->msgStream(LogInfo) << me_ << "Total time taken           = " 
    << stats_->totalTime << std::endl
    << me_ << "Total NLPs solved          = " << stats_->totalNLPs
    << std::endl;
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
