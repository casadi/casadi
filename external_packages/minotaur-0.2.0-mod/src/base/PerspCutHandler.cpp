//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file KnapCovHandler.h
 * \brief Declare the KnapCovHandler class for handling knapsack cover 
 * constraints. It generates the cuts whenever they are needed. 
 * \author Serdar Yildiz, Argonne National Laboratory
 */

#include <cmath>

#include "PerspCutHandler.h"
#include "Option.h"
#include "CutManager.h"

using namespace Minotaur;

typedef std::vector<ConstraintPtr>::const_iterator CCIter;
// const std::string PerspCutHandler::me_ = "PerspCutHandler: ";

PerspCutHandler::PerspCutHandler()
  : env_(EnvPtr()),
    minlp_(ProblemPtr()),
    stats_(0),
    isFeas_(true),
    solAbsTol_(1e-5)
{
  // Logger is na abstract class, find a way to make this work.
  // looger_ = (LoggerPtr) new Logger();
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
}

PerspCutHandler::PerspCutHandler(EnvPtr env, ProblemPtr minlp)
  : env_(env),
    minlp_(minlp),
    stats_(0),
    isFeas_(true),
    solAbsTol_(1e-5)
{
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
  // Initialize logger.
  // Initialize statistics.
  stats_ = new PCStats();
  stats_->cuts = 0;
  stats_->time = 0.0;
}

PerspCutHandler::~PerspCutHandler()
{
  if (stats_) {
    //writeStats(logger_->MsgStream(LogInfo));
    delete stats_;
  }
  //env_.reset();
  //minlp_.reset();
}

bool PerspCutHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr, bool &, 
                                 double &)
{
  // Get primal solution.
  const double *x = sol->getPrimal();
  // Evaluation of constraint for a given solution.
  double activity = 0.0;
  int error = 0;

  // Now, we check all perspective cut constraints if the current 
  // relaxation solution violates any of them.
  // Iterators for perspective constraints.
  CCIter it;
  CCIter begin = cons_.begin();
  CCIter end   = cons_.end();
  // Temporary constraint holder.
  ConstraintPtr cons;
  for (it=begin; it!=end; ++it) {
    cons = *it;
    activity = cons->getActivity(x, &error);
    if (activity > cons->getUb() + solAbsTol_ ||
        activity < cons->getLb() - solAbsTol_) {
      isFeas_ = false;
      return false;
    }
  }
  
  // None of the perspective cut constraints is violated.
  return true;
}

void PerspCutHandler::separate(ConstSolutionPtr sol, NodePtr,
                               RelaxationPtr rel, CutManager * cmanager,
                               SolutionPoolPtr, bool *,
                               SeparationStatus * status)
{
  // Check integer feasibility of sol, must add cuts if it is not integral.
  numvars_ = minlp_->getNumVars();
  VariableType type;
  const double * x = sol->getPrimal();
  // Is the relaxation solution is integer feasible.
  bool isintfeas = true;
  // Iterators for variables.
  VariableConstIterator it;
  VariableConstIterator begin = rel->varsBegin();
  VariableConstIterator end   = rel->varsEnd();
  // Temporary variable holder.
  ConstVariablePtr var;
  // Value of variable.
  double value;
  // Index of variable.
  UInt varindex = 0;

  // Check if integrality is satisfied for each integer variable.
  for (it=begin; it!=end; ++it) {
    var = *it;
    type = var->getType();
    if (type==Binary || type==Integer) {
      varindex = var->getIndex();
      value = x[varindex];
      if (fabs(value - floor(value+0.5)) > intTol_) {
        isintfeas = false;
        break;
      }
    }
  } // end of for loop.

  if (isintfeas == false) {
    // It is more efficient to do integrality check here.
    // Generate perspective cuts from current relaxation.
    PerspCutGeneratorPtr persp = 
      (PerspCutGeneratorPtr) new PerspCutGenerator(rel, sol, env_);
    // Add cuts to relaxation by using cut manager.
    CutVector violatedcuts = persp->getViolatedCutList();
    CutIterator itc;
    CutIterator beginc = violatedcuts.begin();
    CutIterator endc   = violatedcuts.end();
    
    cmanager->addCuts(beginc, endc);
    
    // Update statistics by using return from cover cut generator.
    ConstPerspGenStatsPtr perspstats = persp->getStats();
    stats_->cuts += perspstats->cuts;
  }

  // If any cut is added, relaxation should be resolved.
  if (stats_->cuts >= 1) {
    *status = SepaResolve;
  }
}


std::string PerspCutHandler::getName() const
{
  return "PerspCutHandler (Perspective cuts)";
}

void PerspCutHandler::writeStats(std::ostream &) const
{
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
