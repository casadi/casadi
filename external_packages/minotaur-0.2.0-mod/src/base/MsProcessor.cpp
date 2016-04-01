// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file MsProcessor.cpp
 * \brief Implement simple multi-start node-processor for branch-and-bound
 * \author Ashutosh Mahajan, IIT Bombay
 */
#include <cmath> // for INFINITY

#include "MinotaurConfig.h"
#include "Brancher.h"
#include "Engine.h"
#include "Environment.h"
#include "Handler.h"
#include "MsProcessor.h"
#include "Logger.h"
#include "Node.h"
#include "Option.h"
#include "Modification.h"
#include "Relaxation.h"
#include "SolutionPool.h"

using namespace Minotaur;

//#define SPEW 1

const std::string MsProcessor::me_ = "MsProcessor: ";

MsProcessor::MsProcessor()
  : contOnErr_(false),
    cutOff_(INFINITY),
    engine_(EnginePtr()),
    engineStatus_(EngineUnknownStatus),
    numSolutions_(0),
    relaxation_(RelaxationPtr()),
    ws_(WarmStartPtr())
{
  handlers_.clear();
  logger_ = (LoggerPtr) new Logger(LogInfo);
  stats_.inf = 0;
  stats_.opt = 0;
  stats_.prob = 0;
  stats_.proc = 0;
  stats_.ub = 0;
}


MsProcessor::MsProcessor (EnvPtr env, EnginePtr engine,
                            HandlerVector handlers)
  : contOnErr_(false),
    engine_(engine),
    engineStatus_(EngineUnknownStatus),
    numSolutions_(0),
    relaxation_(RelaxationPtr()),
    ws_(WarmStartPtr())
{
  cutOff_ = env->getOptions()->findDouble("obj_cut_off")->getValue();
  handlers_ = handlers;
  logger_ = (LoggerPtr) new Logger((LogLevel)env->getOptions()->
                                   findInt("node_processor_log_level")->
                                   getValue());
  stats_.bra = 0;
  stats_.inf = 0;
  stats_.opt = 0;
  stats_.prob = 0;
  stats_.proc = 0;
  stats_.ub = 0;
}


MsProcessor::~MsProcessor()
{
  handlers_.clear();
  logger_.reset();
  engine_.reset();
}


bool MsProcessor::foundNewSolution()
{
  return (numSolutions_ > 0);
}


Branches MsProcessor::getBranches()
{
  ++stats_.bra;
  return branches_;
}


WarmStartPtr MsProcessor::getWarmStart()
{
  return ws_;
}


bool MsProcessor::isFeasible_(NodePtr node, ConstSolutionPtr sol, 
                              SolutionPoolPtr s_pool, bool &should_prune)
{
  should_prune = false;
  bool is_feas = true;
  HandlerIterator h;
  double inf_meas;

  // visit each handler and check feasibility. Stop on the first
  // infeasibility.
  for (h = handlers_.begin(); h != handlers_.end(); ++h) {
    is_feas = (*h)->isFeasible(sol, relaxation_, should_prune, inf_meas);
    if (is_feas == false || should_prune == true) {
      break;
    }
  }

  if (is_feas == true && h==handlers_.end()) {
    s_pool->addSolution(sol);
    ++numSolutions_;
    node->setStatus(NodeOptimal);
    ++stats_.opt;
    should_prune = true;
  }
  return is_feas;
}


void MsProcessor::process(NodePtr node, RelaxationPtr rel,
                          SolutionPoolPtr s_pool)
{
  bool should_prune = true;
  bool should_resolve;
  BrancherStatus br_status;
  ConstSolutionPtr sol;
  ModVector mods;
  int iter = 0;

  ++stats_.proc;
  relaxation_ = rel;

#if 0
  double *svar = new double[20];
  bool xfeas = true;
  svar[1-1]   = 0.000000000  ;
  svar[2-1]   = 0.000000000  ;
  svar[3-1]   = 1.042899924  ;
  svar[4-1]   = 0.000000000  ;
  svar[5-1]   = 0.000000000  ;
  svar[6-1]   = 0.000000000  ;
  svar[7-1]   = 0.000000000  ;
  svar[8-1]   = 0.000000000  ;
  svar[9-1]   = 0.000000000  ;
  svar[10-1]  = 0.000000000  ;
  svar[11-1]  = 1.746743790  ;
  svar[12-1]  = 0.000000000  ;
  svar[13-1]  = 0.431470884  ;
  svar[14-1]  = 0.000000000  ;
  svar[15-1]  = 0.000000000  ;
  svar[16-1]  = 4.433050274  ;
  svar[17-1]  = 0.000000000  ;
  svar[18-1]  = 15.858931758 ;
  svar[19-1]  = 0.000000000  ;
  svar[20-1]  = 16.486903370 ;

  for (UInt i=0; i<20; ++i) {
    if (svar[i] > rel->getVariable(i)->getUb()+1e-6 || 
        svar[i] < rel->getVariable(i)->getLb()-1e-6) {
      xfeas = false;
      break;
    }
  }
  if (true==xfeas) {
    std::cout << "xsol feasible in node " << node->getId() << std::endl;
  } else {
    std::cout << "xsol NOT feasible in node " << node->getId() << std::endl;
  }
#endif 

  // loop for branching and resolving if necessary.

  while (true) {
    ++iter;
    should_resolve = false;

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "iteration " << iter << std::endl;
#endif

    //relaxation_->write(std::cout);
    solveRelaxation_();

    sol = engine_->getSolution();

    // check if the relaxation is infeasible or if the cost is too high.
    // In either case we can prune. Also set lb of node.
    should_prune = shouldPrune_(node, sol->getObjValue(), s_pool);
    if (should_prune) {
      break;
    }

    // update pseudo-cost from last branching.
    if (iter == 1) {
      brancher_->updateAfterLP(node, sol);
    }

    // check feasibility. if it is feasible, we can still prune this node.
    isFeasible_(node, sol, s_pool, should_prune);
    if (should_prune) {
      break;
    }


    //relaxation_->write(std::cout);

    //save warm start information before branching. This step is expensive.
    ws_ = engine_->getWarmStartCopy();
    branches_ = brancher_->findBranches(relaxation_, node, sol, s_pool, 
                                        br_status, mods);
    if (br_status==PrunedByBrancher) {

      should_prune = true;
      node->setStatus(NodeInfeasible);
      stats_.inf++;
      break;
    } else if (br_status==ModifiedByBrancher) {
      for (ModificationConstIterator miter=mods.begin(); miter!=mods.end();
           ++miter) {
        node->addPMod(*miter);
        (*miter)->applyToProblem(relaxation_);
      }
      should_resolve = true;
    } 
    if (should_resolve == false) {
      break;
    }
  }
#if 0
  if ((true==should_prune || node->getLb() >-4150) && true==xfeas) {
    std::cout << "problem here!\n";
    std::cout << node->getStatus() << "\n";
    rel->write(std::cout);
    exit(0);
  }
#endif

  return;
}


bool MsProcessor::shouldPrune_(NodePtr node, double solval, 
                               SolutionPoolPtr s_pool)
{
  bool should_prune = false;
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "solution value = " << solval
                                << std::endl; 
#endif
  switch (engineStatus_) {
   case (FailedInfeas):
     logger_->msgStream(LogInfo) << me_ << "failed to converge "
                                 << "(infeasible) in node " << node->getId()
                                 << std::endl;
     node->setStatus(NodeInfeasible);
     should_prune = true;
     ++stats_.inf;
     ++stats_.prob;
     break;
   case (ProvenFailedCQInfeas):
     logger_->msgStream(LogInfo) << me_ << "constraint qualification "
                                 << "violated in node " << node->getId()
                                 << std::endl;
     ++stats_.prob;
   case (ProvenInfeasible):
   case (ProvenLocalInfeasible):
     node->setStatus(NodeInfeasible);
     should_prune = true;
     ++stats_.inf;
     break;

   case (ProvenObjectiveCutOff):
     node->setStatus(NodeHitUb);
     should_prune = true;
     ++stats_.ub;
     break;

   case (ProvenUnbounded):
     should_prune = false;
     logger_->msgStream(LogDebug2) << me_ << "problem relaxation is "
                                   << "unbounded!" << std::endl;
     break;

   case (FailedFeas):
     logger_->msgStream(LogInfo) << me_ << "Failed to converge " 
                                 << "(feasible) in node " << node->getId()
                                 << std::endl;
     if (node->getParent()) {
       node->setLb(node->getParent()->getLb());
     } else {
       node->setLb(-INFINITY);
     }
     node->setStatus(NodeContinue);
     ++stats_.prob;
     break;
   case (ProvenFailedCQFeas):
     logger_->msgStream(LogInfo) << me_ << "constraint qualification "
                                 << "violated in node " << node->getId()
                                 << std::endl;
     if (node->getParent()) {
       node->setLb(node->getParent()->getLb());
     } else {
       node->setLb(-INFINITY);
     }
     node->setStatus(NodeContinue);
     ++stats_.prob;
     break;
   case (EngineIterationLimit):
     ++stats_.prob;
     logger_->msgStream(LogInfo) << me_ << "engine hit iteration limit, "
                                 << "continuing in node " << node->getId()
                                 << std::endl;
     // continue with this node by following ProvenLocalOptimal case.
   case (ProvenLocalOptimal):
   case (ProvenOptimal):
     node->setLb(solval);
     if (solval >= s_pool->getBestSolutionValue() || solval >= cutOff_) {
       node->setStatus(NodeHitUb);
       should_prune = true;
       ++stats_.ub;
     } else {
       should_prune = false;
       node->setStatus(NodeContinue);
     }
     break;
   case (EngineError):
     if (contOnErr_) {
       logger_->msgStream(LogError) << me_ << "engine reports error, "
                                    << " continuing in node " << node->getId()
                                    << std::endl;
       node->setStatus(NodeContinue);
       if (node->getParent()) {
         node->setLb(node->getParent()->getLb());
       } else {
         node->setLb(-INFINITY);
       }
     } else {
       logger_->msgStream(LogError) << me_ << "engine reports error, "
                                    << "pruning node " << node->getId()
                                    << std::endl;
       should_prune = true;
       node->setStatus(NodeInfeasible);
       ++stats_.inf;
     }
     ++stats_.prob;
     break;
   default:
     break;
  }

  return should_prune;
}


void MsProcessor::solveRelaxation_() 
{
  engineStatus_ = EngineError;
  engine_->solve();
  engineStatus_ = engine_->getStatus();
#if SPEW
  logger_->msgStream(LogDebug2) << me_ << "solving relaxation" << std::endl
                                << me_ << "engine status = " 
                                << engine_->getStatusString() << std::endl;
#endif
}


void MsProcessor::writeStats(std::ostream &out) const
{
  out << me_ << "nodes processed     = " << stats_.proc << std::endl 
      << me_ << "nodes branched      = " << stats_.bra << std::endl 
      << me_ << "nodes infeasible    = " << stats_.inf << std::endl 
      << me_ << "nodes optimal       = " << stats_.opt << std::endl 
      << me_ << "nodes hit ub        = " << stats_.ub << std::endl 
      << me_ << "nodes with problems = " << stats_.prob << std::endl 
      ;
}


void MsProcessor::writeStats() const
{
  writeStats(logger_->msgStream(LogNone));
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
