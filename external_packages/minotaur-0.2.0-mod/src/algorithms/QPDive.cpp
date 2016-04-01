//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file QPDive.cpp
 * \brief The main function for solving instances in ampl format (.nl) by
 * using QP based diving in Branch-and-Bound.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <iomanip>
#include <iostream>

#include "MinotaurConfig.h"
#include "BranchAndBound.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "IntVarHandler.h"
#include "LexicoBrancher.h"
#include "LinearHandler.h"
#include "LPEngine.h"
#include "Logger.h"
#include "MaxVioBrancher.h"
#include "NLPEngine.h"
#include "NlPresHandler.h"
#include "Objective.h"
#include "Option.h"
#include "Presolver.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "QPDProcessor.h"
#include "QPDRelaxer.h"
#include "QPEngine.h"
#include "Relaxation.h"
#include "ReliabilityBrancher.h"
#include "Solution.h"
#include "Timer.h"
#include "TreeManager.h"

#include "AMPLInterface.h"

using namespace Minotaur;

BranchAndBound * createBab (EnvPtr env, ProblemPtr p, EnginePtr e, 
                            EnginePtr qe, HandlerVector &handlers)
{
  BranchAndBound *bab = new BranchAndBound(env, p);
  QPDProcessorPtr nproc;
  IntVarHandlerPtr v_hand = (IntVarHandlerPtr) new IntVarHandler(env, p);
  LinearHandlerPtr l_hand = (LinearHandlerPtr) new LinearHandler(env, p);
  BrancherPtr br;
  QPDRelaxerPtr nr;
  const std::string me("qpd: ");

  v_hand->setModFlags(false, true);
  l_hand->setModFlags(false, true);
  handlers.push_back(v_hand);
  handlers.push_back(l_hand);
  nproc = (QPDProcessorPtr) new QPDProcessor(env, p, e, qe, handlers);
  if (env->getOptions()->findString("brancher")->getValue() == "rel") {
    UInt t = 0;
    ReliabilityBrancherPtr rbr = (ReliabilityBrancherPtr)
      new ReliabilityBrancher(env, handlers);
    rbr->setEngine(qe);
    t = (p->getSize()->ints + p->getSize()->bins)/10;
    t = std::max(t, (UInt) 2);
    t = std::min(t, (UInt) 4);
    rbr->setThresh(t);
    env->getLogger()->msgStream(LogExtraInfo) 
      << me << "setting reliability threshhold to " << t << std::endl;
    t = (UInt) p->getSize()->ints + p->getSize()->bins/20+2;
    t = std::min(t, (UInt) 10);
    rbr->setMaxDepth(t);
    env->getLogger()->msgStream(LogExtraInfo)
      << me << "setting reliability maxdepth to " << t << std::endl;
    if (e->getName()=="Filter-SQP") {
      rbr->setIterLim(5);
    }
    env->getLogger()->msgStream(LogExtraInfo)
      << me << "reliability branching iteration limit = " 
      << rbr->getIterLim() << std::endl;
    rbr->setTrustCutoff(false);
    nproc->setBrancher(rbr);
  } else if (env->getOptions()->findString("brancher")->getValue() == 
             "maxvio") {
    br    = (MaxVioBrancherPtr) new MaxVioBrancher(env, handlers);
    nproc->setBrancher(br);
  } else if (env->getOptions()->findString("brancher")->getValue() == "lex") {
    br = (LexicoBrancherPtr) new LexicoBrancher(env, handlers);
    nproc->setBrancher(br);
  }
  bab->setNodeProcessor(nproc);

  nr = (QPDRelaxerPtr) new QPDRelaxer(env, p, qe, e);
  bab->setNodeRelaxer(nr);
  bab->shouldCreateRoot(true);

  return bab;

}


EnginePtr getEngine(EnvPtr env, ProblemPtr p)
{
  EngineFactory *efac = new EngineFactory(env);
  EnginePtr e = EnginePtr(); // NULL

  p->calculateSize();
  if (p->isLinear()) {
    e = efac->getLPEngine();
    if (e) {
      delete efac;
      return e;
    } 
  }

  e = efac->getNLPEngine();
  assert (e || (!"No engine available for this problem."));
  delete efac;
  return e;
}


PresolverPtr presolve(EnvPtr env, ProblemPtr p, size_t ndefs,
                      HandlerVector &handlers)
{
  // create handlers for presolve
  PresolverPtr pres = PresolverPtr(); // NULL
  p->calculateSize();
  if (env->getOptions()->findBool("presolve")->getValue() == true) {
    LinearHandlerPtr lhandler = (LinearHandlerPtr) new LinearHandler(env, p);
    handlers.push_back(lhandler);
    if (p->isQP() || p->isQuadratic() || p->isLinear() ||
        true==env->getOptions()->findBool("use_native_cgraph")->getValue()) {
      lhandler->setPreOptPurgeVars(true);
      lhandler->setPreOptPurgeCons(true);
      lhandler->setPreOptCoeffImp(true);
    } else {
      lhandler->setPreOptPurgeVars(false);
      lhandler->setPreOptPurgeCons(false);
      lhandler->setPreOptCoeffImp(false);
    }
    if (ndefs>0) {
      lhandler->setPreOptDualFix(false);
    } else {
      lhandler->setPreOptDualFix(true);
    }

    if (!p->isLinear() && 
         true==env->getOptions()->findBool("use_native_cgraph")->getValue() && 
         true==env->getOptions()->findBool("nl_presolve")->getValue() 
         ) {
      NlPresHandlerPtr nlhand = (NlPresHandlerPtr) new NlPresHandler(env, p);
      handlers.push_back(nlhand);
    }

    // write the names.
    env->getLogger()->msgStream(LogExtraInfo) << "handlers used in presolve:" << std::endl;
    for (HandlerIterator h = handlers.begin(); h != handlers.end(); 
        ++h) {
      env->getLogger()->msgStream(LogExtraInfo)<<(*h)->getName()<<std::endl;
    }
  }

  pres = (PresolverPtr) new Presolver(p, env, handlers);
  pres->standardize(); 
  if (env->getOptions()->findBool("presolve")->getValue() == true) {
    pres->solve();
  }
  return pres;
}


void showHelp()
{
  std::cout << "QP-Diving algorithm for convex MINLP" << std::endl
            << "Usage:" << std::endl
            << "To show version: qpd -v (or --show_version yes) " << std::endl
            << "To show all options: qpd -= (or --show_options yes)" 
            << std::endl
            << "To solve an instance: qpd --option1 [value] "
            << "--option2 [value] ... " << " .nl-file" << std::endl;
}


int showInfo(EnvPtr env)
{
  OptionDBPtr options = env->getOptions();
  const std::string me("qpd: ");

  // check if user needs help.
  if (options->findBool("show_options")->getValue() ||
      options->findFlag("=")->getValue()) {
    options->write(std::cout);
    return 1;
  }

  if (options->findBool("show_help")->getValue() ||
      options->findFlag("?")->getValue()) {
    showHelp();
    return 1;
  }

  if (options->findBool("show_version")->getValue() ||
      options->findFlag("v")->getValue()) {
    env->getLogger()->msgStream(LogNone) << me << "Minotaur version "
      << env->getVersion() << std::endl;
#if DEBUG
    env->getLogger()->msgStream(LogNone) << me; 
    env->writeFullVersion(env->getLogger()->msgStream(LogNone));
    env->getLogger()->msgStream(LogNone) << std::endl;
#endif
    return 1;
  }

  if (options->findString("problem_file")->getValue()=="") {
    showHelp();
    return 1;
  }

  env->getLogger()->msgStream(LogInfo) << me << "Minotaur version "
    << env->getVersion() << std::endl
    << me << "QP-Diving algorithm for convex MINLP" << std::endl;

  return 0;
}


void writeBnbStatus(EnvPtr env, BranchAndBound *bab, double obj_sense)
{

  const std::string me("qpd: ");
  int err = 0;

  if (bab) {
    env->getLogger()->msgStream(LogInfo)
      << me << std::fixed << std::setprecision(4) 
      << "best solution value = " << obj_sense*bab->getUb() << std::endl
      << me << std::fixed << std::setprecision(4)
      << "best bound estimate from remaining nodes = "
      << obj_sense*bab->getLb() << std::endl
      << me << "gap = " << std::max(0.0,bab->getUb() - bab->getLb())
      << std::endl
      << me << "gap percentage = " << bab->getPerGap() << std::endl
      << me << "time used (s) = " << std::fixed << std::setprecision(2) 
      << env->getTime(err) << std::endl
      << me << "status of branch-and-bound: " 
      << getSolveStatusString(bab->getStatus()) << std::endl;
    env->stopTimer(err); assert(0==err);
  } else {
    env->getLogger()->msgStream(LogInfo)
      << me << std::fixed << std::setprecision(4)
      << "best solution value = " << INFINITY << std::endl
      << me << std::fixed << std::setprecision(4)
      << "best bound estimate from remaining nodes = " << INFINITY << std::endl
      << me << "gap = " << INFINITY << std::endl
      << me << "gap percentage = " << INFINITY << std::endl
      << me << "time used (s) = " << std::fixed << std::setprecision(2) 
      << env->getTime(err) << std::endl 
      << me << "status of branch-and-bound: " 
      << getSolveStatusString(NotStarted) << std::endl;
    env->stopTimer(err); assert(0==err);
  }
}


void writeSol(EnvPtr env, VarVector *orig_v,
              PresolverPtr pres, SolutionPtr sol, SolveStatus status,
              MINOTAUR_AMPL::AMPLInterface* iface)
{
  if (sol) {
    sol = pres->getPostSol(sol);
  }

  if (env->getOptions()->findFlag("AMPL")->getValue()) {
    iface->writeSolution(sol, status);
  } else if (sol && env->getLogger()->getMaxLevel()>=LogExtraInfo) {
    sol->writePrimal(env->getLogger()->msgStream(LogExtraInfo), orig_v);
  }
}


int main(int argc, char** argv)
{
  EnvPtr env      = (EnvPtr) new Environment();
  OptionDBPtr options;

  // interface to AMPL (NULL)
  MINOTAUR_AMPL::AMPLInterfacePtr iface = MINOTAUR_AMPL::AMPLInterfacePtr();  
  ProblemPtr inst;    // instance that needs to be solved
  EnginePtr engine;    
  EnginePtr qe = EnginePtr();
  SolutionPtr sol, sol2;
  BranchAndBound * bab = 0; // the branch-and-bound
  PresolverPtr pres;
  VarVector *orig_v=0;
  const std::string me("qpd: ");
  EngineFactory *efac;
  HandlerVector handlers;
  int err = 0;
  double obj_sense = 1.0;

  // start timing.
  env->startTimer(err);
  if (err) {
    goto CLEANUP;
  }


  options = env->getOptions();
  options->findString("nlp_engine")->setValue("IPOPT");
  iface = (MINOTAUR_AMPL::AMPLInterfacePtr) 
    new MINOTAUR_AMPL::AMPLInterface(env, "qpd");

  // set default options that can be overriden
  options->findBool("presolve")->setValue(true);
  // parse options
  env->readOptions(argc, argv);

  // override some options
  options->findString("interface_type")->setValue("AMPL");
  options->findBool("use_native_cgraph")->setValue(true);

  if (0!=showInfo(env)) {
    goto CLEANUP;
  }

  // load the problem.
  inst = iface->readInstance(options->findString("problem_file")->getValue());
  env->getLogger()->msgStream(LogInfo) << me << "time used in reading instance = " << std::fixed 
    << std::setprecision(2) << env->getTime(err) << std::endl; assert(0==err);

  // display the problem
  inst->calculateSize();
  if (options->findBool("display_problem")->getValue()==true) {
    inst->write(env->getLogger()->msgStream(LogNone));
  }
  if (options->findBool("display_size")->getValue()==true) {
    inst->writeSize(env->getLogger()->msgStream(LogNone));
  }

  // Get the right engine.
  engine = getEngine(env, inst);
  efac = new EngineFactory(env);
  qe = efac->getQPEngine();
  delete efac;
  assert(qe || (!"QP Engine not available!")); 

  // set initial point
  inst->setInitialPoint(iface->getInitialPoint(), 
      inst->getNumVars()-iface->getNumDefs());
   
  if (inst->getObjective() &&
      inst->getObjective()->getObjectiveType()==Maximize) {
    obj_sense = -1.0;
  }
  // get presolver.
  orig_v = new VarVector(inst->varsBegin(), inst->varsEnd());
  pres = presolve(env, inst, iface->getNumDefs(), handlers);
  handlers.clear();

  if (Finished != pres->getStatus() && NotStarted != pres->getStatus()) {
    env->getLogger()->msgStream(LogInfo) << me 
      << "status of presolve: " 
      << getSolveStatusString(pres->getStatus()) << std::endl;
    writeSol(env, orig_v, pres, SolutionPtr(), pres->getStatus(), iface);
    writeBnbStatus(env, bab, obj_sense);
    goto CLEANUP;
  }

  if (false==env->getOptions()->findBool("solve")->getValue()) {
    goto CLEANUP;
  }

  inst->setNativeDer();

  // get branch-and-bound
  bab = createBab(env, inst, engine, qe, handlers);

  // start solving
  bab->solve();
  bab->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  engine->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  qe->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  for (HandlerVector::iterator it=handlers.begin(); it!=handlers.end(); ++it) {
    (*it)->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  }

  writeSol(env, orig_v, pres, bab->getSolution(), bab->getStatus(), iface);
  writeBnbStatus(env, bab, obj_sense);

CLEANUP:
  if (iface) {
    delete iface;
  }
  if (bab) {
    delete bab;
  }
  if (orig_v) {
    delete orig_v;
  }

  return 0;
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
