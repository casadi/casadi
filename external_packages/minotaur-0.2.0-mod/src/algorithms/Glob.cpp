//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Glob.cpp
 * \brief The main function for solving nonconvex optimization instances with
 * quadratic functions in ampl format (.nl).
 * \author Ashutosh Mahajan, IIT Bombay
 */

#include <iomanip>
#include <iostream>

#include "MinotaurConfig.h"
#include "BranchAndBound.h"
#include "Engine.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "LexicoBrancher.h"
#include "LinearHandler.h"
#include "Logger.h"
#include "LPEngine.h"
#include "LPProcessor.h"
#include "MaxVioBrancher.h"
#include "NodeIncRelaxer.h"
#include "NLPEngine.h"
#include "NLPMultiStart.h"
#include "NlPresHandler.h"
#include "Objective.h"
#include "Option.h"
#include "Presolver.h"
#include "ProblemSize.h"
#include "Problem.h"
#include "Relaxation.h"
#include "ReliabilityBrancher.h"
#include "SimpleTransformer.h"
#include "Solution.h"
#include "Timer.h"
#include "Transformer.h"
#include "TreeManager.h"

#include "AMPLInterface.h"

using namespace Minotaur;

BranchAndBound *createBab (EnvPtr env, ProblemPtr p, EnginePtr e, 
                           HandlerVector &handlers);
PresolverPtr createPres(EnvPtr env, ProblemPtr p, size_t ndefs, 
                        HandlerVector &handlers);
EnginePtr getEngine(EnvPtr env);
NLPEnginePtr getNLPEngine(EnvPtr env);
int transform(EnvPtr env, ProblemPtr p, ProblemPtr &newp,
              HandlerVector &handlers);



EnginePtr getEngine(EnvPtr env)
{
  EngineFactory *efac = new EngineFactory(env);
  EnginePtr e = EnginePtr(); // NULL
  const std::string me("mntr-glob: ");
  e = efac->getLPEngine();
  if (!e) {
    env->getLogger()->errStream() << me 
      << "Cannot find an LP engine. Cannot proceed!" << std::endl;
  } 

  delete efac;
  return e;
}


NLPEnginePtr getNLPEngine(EnvPtr env)
{
  EngineFactory *efac = new EngineFactory(env);
  NLPEnginePtr e = NLPEnginePtr(); // NULL
  const std::string me("mntr-glob: ");
  e = efac->getNLPEngine();
  if (!e) {
    env->getLogger()->errStream() << me 
      << "Cannot find an NLP engine. Cannot proceed!" << std::endl;
  } 

  delete efac;
  return e;
}


void loadProblem(EnvPtr env, MINOTAUR_AMPL::AMPLInterface* iface,
                 ProblemPtr &inst, double *obj_sense)
{
  OptionDBPtr options = env->getOptions();
  Timer *timer     = env->getNewTimer();
  const std::string me("mntr-glob: ");

  if (options->findBool("use_native_cgraph")->getValue()==false) {
    options->findBool("use_native_cgraph")->setValue(true); 
    env->getLogger()->msgStream(LogExtraInfo) << me 
      << "Setting value of 'use_native_cgraph option' to True" << std::endl;
  }

  // load the problem.
  timer->start();
  inst = iface->readInstance(options->findString("problem_file")->getValue());
  env->getLogger()->msgStream(LogInfo) << me 
    << "time used in reading instance = " << std::fixed 
    << std::setprecision(2) << timer->query() << std::endl;

  // display the problem
  inst->calculateSize();
  if (options->findBool("display_problem")->getValue()==true) {
    inst->write(env->getLogger()->msgStream(LogNone), 9);
  }
  if (options->findBool("display_size")->getValue()==true) {
    inst->writeSize(env->getLogger()->msgStream(LogNone));
  }

  if (inst->getObjective() &&
      inst->getObjective()->getObjectiveType()==Maximize) {
    *obj_sense = -1.0;
    env->getLogger()->msgStream(LogInfo) << me 
      << "objective sense: maximize (will be converted to Minimize)"
      << std::endl;
  } else {
    *obj_sense = 1.0;
    env->getLogger()->msgStream(LogInfo) << me 
      << "objective sense: minimize" << std::endl;
  }
  delete timer;
}


int transform(EnvPtr env, ProblemPtr p, ProblemPtr &newp,
              HandlerVector &handlers) 
{
  SimpTranPtr trans = SimpTranPtr();
  int status = 0;
  const std::string me("mntr-glob: ");

  handlers.clear();
  trans = (SimpTranPtr) new SimpleTransformer(env, p);
  trans->reformulate(newp, handlers, status);
  
  env->getLogger()->errStream() << me 
    << "handlers used in transformer: " << std::endl;
  for (HandlerVector::iterator it=handlers.begin(); it!=handlers.end();
       ++it) {
    env->getLogger()->errStream() << "  " << (*it)->getName() << std::endl;
  }
  return status;
}


BranchAndBound * createBab (EnvPtr env, ProblemPtr p, EnginePtr e, 
                            HandlerVector &handlers)
{
  BranchAndBound *bab = new BranchAndBound(env, p);
  LPProcessorPtr nproc;
  NodeIncRelaxerPtr nr;
  RelaxationPtr rel;
  BrancherPtr br;
  const std::string me("mntr-glob: ");

  if (env->getOptions()->findString("brancher")->getValue() == "rel") {
    UInt t;
    ReliabilityBrancherPtr rel_br;
    rel_br = (ReliabilityBrancherPtr) new ReliabilityBrancher(env, handlers);
    rel_br->setEngine(e);
    t = (p->getSize()->ints + p->getSize()->bins)/10;
    t = std::max(t, (UInt) 2);
    t = std::min(t, (UInt) 4);
    rel_br->setThresh(t);
    env->getLogger()->msgStream(LogExtraInfo) << me 
      << "setting reliability threshhold to " << t << std::endl;
    t = (UInt) p->getSize()->ints + p->getSize()->bins/20+2;
    t = std::min(t, (UInt) 10);
    rel_br->setMaxDepth(t);
    env->getLogger()->msgStream(LogExtraInfo) << me 
      << "setting reliability maxdepth to " << t << std::endl;
    env->getLogger()->msgStream(LogExtraInfo) << me
      << "reliability branching iteration limit = " 
      << rel_br->getIterLim() << std::endl;
    br = rel_br;
  } else if (env->getOptions()->findString("brancher")->getValue() 
      == "maxvio") {
    MaxVioBrancherPtr mbr = 
      (MaxVioBrancherPtr) new MaxVioBrancher(env, handlers);
    br = mbr;
  } else if (env->getOptions()->findString("brancher")->getValue() 
      == "lex") {
    LexicoBrancherPtr lbr = 
      (LexicoBrancherPtr) new LexicoBrancher(env, handlers);
    br = lbr;
  }
  env->getLogger()->msgStream(LogExtraInfo) << me 
    << "brancher used = " << br->getName() << std::endl;
  nproc = (LPProcessorPtr) new LPProcessor(env, e, handlers);
  nproc->setBrancher(br);
  bab->setNodeProcessor(nproc);

  nr = (NodeIncRelaxerPtr) new NodeIncRelaxer(env, handlers);
  nr->setProblem(p);
  nr->setRelaxation(rel);
  nr->setEngine(e);
  bab->setNodeRelaxer(nr);
  bab->shouldCreateRoot(true);

  if (env->getOptions()->findBool("msheur")->getValue() == true) {
    EnginePtr nlp_e = getNLPEngine(env);
    p->setNativeDer();
    NLPMSPtr ms_heur = (NLPMSPtr) new NLPMultiStart(env, p, nlp_e);
    bab->addPreRootHeur(ms_heur); 
  }

  return bab;
}


PresolverPtr createPres(EnvPtr env, ProblemPtr p, size_t ndefs, 
                        HandlerVector &handlers)
{
  // create handlers for presolve
  PresolverPtr pres = PresolverPtr(); // NULL
  const std::string me("mntr-glob: ");

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
    env->getLogger()->msgStream(LogExtraInfo) << me 
      << "handlers used in presolve:" << std::endl;
    for (HandlerIterator h = handlers.begin(); h != handlers.end(); 
        ++h) {
      env->getLogger()->msgStream(LogExtraInfo) << "  "
        << (*h)->getName()<<std::endl;
    }
  }

  pres = (PresolverPtr) new Presolver(p, env, handlers);
  return pres;
}


void setInitialOptions(EnvPtr env)
{
  OptionDBPtr options = env->getOptions();
  options->findString("interface_type")->setValue("AMPL");
  options->findBool("presolve")->setValue(true);
  options->findBool("nl_presolve")->setValue(true);
  options->findBool("lin_presolve")->setValue(true);
  options->findBool("msheur")->setValue(true);
  options->findString("brancher")->setValue("maxvio");
}

void showHelp()
{
  std::cout << "global optimization for general QCQP" << std::endl
            << "Usage:" << std::endl
            << "To show version: glob -v (or --show_version yes) " << std::endl
            << "To show all options: glob -= (or --show_options yes)" 
            << std::endl
            << "To solve an instance: glob --option1 [value] "
            << "--option2 [value] ... " << " .nl-file" << std::endl;
}


int showInfo(EnvPtr env)
{
  OptionDBPtr options = env->getOptions();
  const std::string me("mntr-glob: ");

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
    env->getLogger()->msgStream(LogNone) << me <<
      "Minotaur version " << env->getVersion() << std::endl;
#if DEBUG
    env->getLogger()->msgStream(LogInfo) << me;
    env->writeFullVersion(env->getLogger()->msgStream(LogInfo));
    env->getLogger()->msgStream(LogInfo) << std::endl;
#endif
    env->getLogger()->msgStream(LogNone) << me 
      << "global optimization for nonconvex QCQP" << std::endl;
    return 1;
  }

  if (options->findString("problem_file")->getValue()=="") {
    showHelp();
    return 1;
  }

  env->getLogger()->msgStream(LogInfo)
    << me << "Minotaur version " << env->getVersion() << std::endl
    << me << "global optimization for nonconvex QCQP" << std::endl;
  return 0;
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


void writeStatus(EnvPtr env, BranchAndBound *bab, double obj_sense)
{

  const std::string me("mntr-glob: ");
  int err = 0;

  if (bab) {
    env->getLogger()->msgStream(LogInfo)
      << me << std::fixed << std::setprecision(4) 
      << "best solution value = " << obj_sense*bab->getUb() << std::endl
      << me << std::fixed << std::setprecision(4)
      << "best bound estimate from remaining nodes = "
      <<  obj_sense*bab->getLb() << std::endl
      << me << "gap = " << std::max(0.0,bab->getUb() - bab->getLb())
      << std::endl
      << me << "gap percentage = " << bab->getPerGap() << std::endl
      << me << "time used = " << std::fixed << std::setprecision(2) 
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
      << me << "time used = " << std::fixed << std::setprecision(2) 
      << env->getTime(err) << std::endl 
      << me << "status of branch-and-bound: " 
      << getSolveStatusString(NotStarted) << std::endl;
    env->stopTimer(err); assert(0==err);
  }
}


int main(int argc, char** argv)
{
  EnvPtr env      = (EnvPtr) new Environment();
  OptionDBPtr options;
  MINOTAUR_AMPL::AMPLInterfacePtr iface;
  ProblemPtr inst;     // instance that needs to be solved.
  EnginePtr engine;    // engine for solving relaxations. 
  SolutionPtr sol, sol2;
  BranchAndBound *bab = 0;
  PresolverPtr pres, pres2;
  const std::string me("mntr-glob: ");
  VarVector *orig_v=0;
  HandlerVector handlers;
  double obj_sense = 1.0;
  int err = 0;
  ProblemPtr newp;

  // start timing.
  env->startTimer(err);

  setInitialOptions(env);

  // Important to setup AMPL Interface first as it adds several options.
  iface = (MINOTAUR_AMPL::AMPLInterfacePtr) 
    new MINOTAUR_AMPL::AMPLInterface(env, "mntr-glob");

  // read user-specified options
  env->readOptions(argc, argv);

  if (0!=showInfo(env)) {
    goto CLEANUP;
  }

  loadProblem(env, iface, inst, &obj_sense);

  // Get the right engine.
  engine = getEngine(env);
  env->getLogger()->msgStream(LogExtraInfo) << me 
    << "engine used = " << engine->getName() << std::endl;

  // get presolver.
  handlers.clear();
  orig_v = new VarVector(inst->varsBegin(), inst->varsEnd());
  pres = createPres(env, inst, iface->getNumDefs(), handlers);
  if (env->getOptions()->findBool("presolve")->getValue() == true) {
    pres->solve();
  }
  handlers.clear();

  err = transform(env, inst, newp, handlers);
  assert(0==err);

  env->getLogger()->msgStream(LogExtraInfo) << me 
    << "Presolving transformed problem ... " << std::endl;
  pres2 = (PresolverPtr) new Presolver(newp, env, handlers);

  pres2->solve();
  env->getLogger()->msgStream(LogExtraInfo) << me 
    << "Finished presolving transformed problem" << std::endl;

  // get branch-and-bound
  bab = createBab(env, newp, engine, handlers);

  if (false==env->getOptions()->findBool("solve")->getValue()) {
    goto CLEANUP;
  }

  // solve
  bab->solve();
  bab->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  engine->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  for (HandlerVector::iterator it=handlers.begin(); it!=handlers.end(); ++it) {
    (*it)->writeStats(env->getLogger()->msgStream(LogExtraInfo));
  }

  writeSol(env, orig_v, pres, bab->getSolution(), bab->getStatus(), iface);
  writeStatus(env, bab, obj_sense);

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
